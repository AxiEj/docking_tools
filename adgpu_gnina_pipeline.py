#!/usr/bin/env python3
"""Run AutoDock-GPU docking followed by gnina CNN rescoring.

Inputs: receptor PDB, ligand MOL2/SDF/PDBQT-compatible molecule, and a box text file.
Designed for robust single-ligand runs in /home/axie/docking_tools.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import shlex
import subprocess
import sys
import textwrap
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence

ROOT = Path(__file__).resolve().parent


@dataclass
class Box:
    center: tuple[float, float, float]
    npts: tuple[int, int, int]
    spacing: float
    size: tuple[float, float, float] | None = None
    source_notes: list[str] | None = None


class PipelineError(RuntimeError):
    pass


def log(msg: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def which(name: str) -> str | None:
    return shutil.which(name)


def require_tool(name: str) -> str:
    path = which(name)
    if not path:
        raise PipelineError(f"Required tool not found in PATH: {name}")
    return path


def split_extra_args(extra: str, label: str) -> list[str]:
    if not extra:
        return []
    try:
        return shlex.split(extra)
    except ValueError as exc:
        raise PipelineError(f"Could not parse {label}; check shell quoting: {extra!r}") from exc


def run_cmd(
    cmd: Sequence[str],
    *,
    cwd: Path,
    log_path: Path,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    cmd_str = " ".join(subprocess.list2cmdline([c]) for c in cmd)
    with log_path.open("w", encoding="utf-8", errors="replace") as fh:
        fh.write(f"$ {cmd_str}\n")
        fh.write(f"# cwd={cwd}\n\n")
        fh.flush()
        proc = subprocess.run(
            list(cmd),
            cwd=str(cwd),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=env,
        )
        fh.write(proc.stdout)
        fh.write(f"\n# exit_code={proc.returncode}\n")
    if check and proc.returncode != 0:
        raise PipelineError(
            f"Command failed with exit {proc.returncode}: {cmd_str}\nLog: {log_path}"
        )
    return proc


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def parse_floats(pattern: str, text: str, count: int) -> tuple[float, ...] | None:
    m = re.search(pattern, text, flags=re.I | re.M)
    if not m:
        return None
    vals = tuple(float(m.group(i + 1)) for i in range(count))
    return vals


def parse_box(path: Path, default_spacing: float = 0.375) -> Box:
    text = read_text(path).replace("\r", "")
    notes: list[str] = []

    npts_match = re.search(
        r"^\s*npts\s+([0-9]+)\s+([0-9]+)\s+([0-9]+)\b", text, flags=re.I | re.M
    )
    npts = tuple(int(npts_match.group(i)) for i in range(1, 4)) if npts_match else None
    if npts:
        notes.append("used explicit AutoGrid npts")

    spacing_match = re.search(r"^\s*spacing\s+([-+0-9.eE]+)\b", text, flags=re.I | re.M)
    spacing = float(spacing_match.group(1)) if spacing_match else default_spacing
    if spacing_match:
        notes.append("used explicit AutoGrid spacing")
    else:
        notes.append(f"spacing missing; used default {default_spacing}")

    gridcenter = parse_floats(
        r"^\s*gridcenter\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)",
        text,
        3,
    )
    if gridcenter:
        center = gridcenter
        notes.append("used explicit AutoGrid gridcenter")
    else:
        cx = re.search(r"--center_x\s+([-+0-9.eE]+)", text, flags=re.I)
        cy = re.search(r"--center_y\s+([-+0-9.eE]+)", text, flags=re.I)
        cz = re.search(r"--center_z\s+([-+0-9.eE]+)", text, flags=re.I)
        if cx and cy and cz:
            center = (float(cx.group(1)), float(cy.group(1)), float(cz.group(1)))
            notes.append("gridcenter missing; used Vina --center values")
        else:
            # LeDock-style three min/max lines after "Binding pocket"
            ledock = re.search(
                r"Binding\s+pocket\s*\n\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\n\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\n\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)",
                text,
                flags=re.I,
            )
            if not ledock:
                raise PipelineError(
                    f"Could not parse box center from {path}. Need gridcenter or --center_x/y/z."
                )
            vals = [float(ledock.group(i)) for i in range(1, 7)]
            center = ((vals[0] + vals[1]) / 2, (vals[2] + vals[3]) / 2, (vals[4] + vals[5]) / 2)
            notes.append("gridcenter missing; derived center from LeDock min/max")

    sx = re.search(r"--size_x\s+([-+0-9.eE]+)", text, flags=re.I)
    sy = re.search(r"--size_y\s+([-+0-9.eE]+)", text, flags=re.I)
    sz = re.search(r"--size_z\s+([-+0-9.eE]+)", text, flags=re.I)
    size = None
    if sx and sy and sz:
        size = (float(sx.group(1)), float(sy.group(1)), float(sz.group(1)))
        notes.append("parsed Vina box size")
    else:
        ledock = re.search(
            r"Binding\s+pocket\s*\n\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\n\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\n\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)",
            text,
            flags=re.I,
        )
        if ledock:
            vals = [float(ledock.group(i)) for i in range(1, 7)]
            size = (vals[1] - vals[0], vals[3] - vals[2], vals[5] - vals[4])
            notes.append("derived Vina-equivalent size from LeDock min/max")

    if not npts:
        if not size:
            raise PipelineError(
                f"Could not parse AutoGrid npts or Vina/LeDock size from {path}."
            )
        npts = tuple(max(1, int(math.ceil(v / spacing))) for v in size)
        notes.append("npts missing; derived npts=ceil(size/spacing)")

    if any(v <= 0 for v in npts) or spacing <= 0:
        raise PipelineError(f"Invalid grid box: npts={npts}, spacing={spacing}")
    if any(v > 1025 for v in npts):
        raise PipelineError(f"npts too large for AutoGrid defaults: {npts}")

    return Box(center=tuple(center), npts=tuple(npts), spacing=spacing, size=size, source_notes=notes)


def normalize_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    data = src.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n").replace("\r", "\n")
    dst.write_text(data, encoding="utf-8")


def pdbqt_atom_types(path: Path) -> list[str]:
    types: list[str] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            tail = line[77:].strip() if len(line) >= 78 else ""
            atom_type = tail.split()[0] if tail else (line.split()[-1] if line.split() else "")
            atom_type = atom_type.strip()
            if atom_type and atom_type not in types:
                types.append(atom_type)
    return types


def atom_count(path: Path) -> int:
    return sum(1 for line in path.read_text(encoding="utf-8", errors="replace").splitlines() if line.startswith(("ATOM", "HETATM")))


def prepare_ligand(args: argparse.Namespace, prep: Path, logs: Path, warnings: list[str]) -> Path:
    ligand_in = Path(args.ligand).resolve()
    ligand_copy = prep / f"ligand_input{ligand_in.suffix.lower()}"
    shutil.copy2(ligand_in, ligand_copy)
    out = prep / "ligand.pdbqt"

    log("Preparing ligand with Meeko")
    proc = run_cmd(
        ["mk_prepare_ligand.py", "-i", str(ligand_copy), "-o", str(out), "--rename_atoms"],
        cwd=prep,
        log_path=logs / "mk_prepare_ligand.log",
        check=False,
    )
    if proc.returncode != 0 or not out.exists() or atom_count(out) == 0:
        warnings.append("mk_prepare_ligand.py failed on ligand input; trying Open Babel -> SDF fallback")
        require_tool("obabel")
        sdf = prep / "ligand_obabel.sdf"
        fmt = ligand_in.suffix.lower().lstrip(".") or "mol2"
        run_cmd(
            ["obabel", f"-i{fmt}", str(ligand_copy), "-osdf", "-O", str(sdf), "-h"],
            cwd=prep,
            log_path=logs / "obabel_ligand_to_sdf.log",
        )
        run_cmd(
            ["mk_prepare_ligand.py", "-i", str(sdf), "-o", str(out), "--rename_atoms"],
            cwd=prep,
            log_path=logs / "mk_prepare_ligand_from_sdf.log",
        )

    if not out.exists() or atom_count(out) == 0:
        raise PipelineError("Ligand PDBQT was not created or contains no atoms")
    types = pdbqt_atom_types(out)
    if not types:
        raise PipelineError("Ligand PDBQT contains no AutoDock atom types")
    log(f"Ligand ready: {out.name}; atoms={atom_count(out)}; types={','.join(types)}")
    return out


def try_prepare_receptor_meeko(
    receptor_pdb: Path, prep: Path, logs: Path, box: Box, mode: str
) -> Path | None:
    outbase = prep / f"receptor_{mode}"
    box_size = box.size if box.size else tuple(n * box.spacing for n in box.npts)
    if mode == "prody":
        cmd = [
            "mk_prepare_receptor.py",
            "-i",
            str(receptor_pdb),
            "-o",
            str(outbase),
            "-p",
            "--box_center",
            *(str(x) for x in box.center),
            "--box_size",
            *(str(x) for x in box_size),
            "-a",
        ]
    else:
        cmd = [
            "mk_prepare_receptor.py",
            "--read_pdb",
            str(receptor_pdb),
            "-o",
            str(outbase),
            "-p",
            "--box_center",
            *(str(x) for x in box.center),
            "--box_size",
            *(str(x) for x in box_size),
            "-a",
        ]
    proc = run_cmd(cmd, cwd=prep, log_path=logs / f"mk_prepare_receptor_{mode}.log", check=False)
    candidate = Path(str(outbase) + ".pdbqt")
    if proc.returncode == 0 and candidate.exists() and atom_count(candidate) > 0:
        return candidate
    return None


def prepare_receptor(
    args: argparse.Namespace, prep: Path, logs: Path, box: Box, warnings: list[str]
) -> tuple[Path, str]:
    receptor_in = Path(args.receptor).resolve()
    receptor_copy = prep / "receptor_input.pdb"
    normalize_copy(receptor_in, receptor_copy)

    final = prep / "receptor.pdbqt"

    # Prefer a prepared same-stem PDBQT when available. This is usually more reliable than
    # re-preparing an already-protonated receptor PDB.
    same_stem = receptor_in.with_suffix(".pdbqt")
    if args.prefer_existing_pdbqt and same_stem.exists():
        normalize_copy(same_stem, final)
        if atom_count(final) > 0 and pdbqt_atom_types(final):
            warnings.append(f"used existing same-stem receptor PDBQT: {same_stem}")
            log(f"Receptor ready from existing PDBQT: {same_stem.name}; atoms={atom_count(final)}")
            return final, "existing_same_stem_pdbqt"
        warnings.append(f"same-stem PDBQT existed but was invalid: {same_stem}")

    log("Preparing receptor with Meeko/ProDy path")
    rec = try_prepare_receptor_meeko(receptor_copy, prep, logs, box, "prody")
    if not rec:
        warnings.append("Meeko ProDy receptor preparation failed; see logs/mk_prepare_receptor_prody.log")
        log("Preparing receptor with Meeko --read_pdb path")
        rec = try_prepare_receptor_meeko(receptor_copy, prep, logs, box, "readpdb")
    if rec:
        normalize_copy(rec, final)
        log(f"Receptor ready from Meeko; atoms={atom_count(final)}")
        return final, "meeko"

    warnings.append("Meeko --read_pdb receptor preparation failed; see logs/mk_prepare_receptor_readpdb.log")

    if (not args.prefer_existing_pdbqt) and same_stem.exists():
        normalize_copy(same_stem, final)
        if atom_count(final) > 0 and pdbqt_atom_types(final):
            warnings.append(f"fallback used existing same-stem receptor PDBQT: {same_stem}")
            log(f"Receptor ready from existing PDBQT fallback: {same_stem.name}; atoms={atom_count(final)}")
            return final, "existing_same_stem_pdbqt_fallback"

    if args.allow_obabel_receptor:
        require_tool("obabel")
        warnings.append(
            "fallback used Open Babel receptor PDBQT generation; charges may be less reliable than a curated receptor PDBQT"
        )
        log("Preparing receptor with Open Babel fallback")
        run_cmd(
            ["obabel", "-ipdb", str(receptor_copy), "-opdbqt", "-O", str(final), "-xr"],
            cwd=prep,
            log_path=logs / "obabel_receptor_to_pdbqt.log",
        )
        if final.exists() and atom_count(final) > 0 and pdbqt_atom_types(final):
            log(f"Receptor ready from Open Babel fallback; atoms={atom_count(final)}")
            return final, "obabel_fallback"

    raise PipelineError(
        "Could not prepare receptor PDBQT. Provide a same-stem .pdbqt next to receptor PDB "
        "or inspect receptor preparation logs."
    )


def write_gpf(prep: Path, receptor_pdbqt: Path, ligand_pdbqt: Path, box: Box) -> Path:
    rec_types = pdbqt_atom_types(receptor_pdbqt)
    lig_types = pdbqt_atom_types(ligand_pdbqt)
    if not rec_types:
        raise PipelineError(f"No receptor atom types found in {receptor_pdbqt}")
    if not lig_types:
        raise PipelineError(f"No ligand atom types found in {ligand_pdbqt}")
    gpf = prep / "receptor.gpf"
    lines = [
        f"npts {box.npts[0]} {box.npts[1]} {box.npts[2]}",
        "gridfld receptor.maps.fld",
        f"spacing {box.spacing:g}",
        "receptor_types " + " ".join(rec_types),
        "ligand_types " + " ".join(lig_types),
        f"receptor {receptor_pdbqt.name}",
        f"gridcenter {box.center[0]:.3f} {box.center[1]:.3f} {box.center[2]:.3f}",
        "smooth 0.500",
    ]
    for t in lig_types:
        lines.append(f"map receptor.{t}.map")
    lines.extend(["elecmap receptor.e.map", "dsolvmap receptor.d.map", "dielectric -42.000"])
    gpf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"GPF written: {gpf}; receptor_types={','.join(rec_types)}; ligand_types={','.join(lig_types)}")
    return gpf


def run_autogrid(prep: Path, logs: Path, gpf: Path) -> Path:
    log("Running AutoGrid4")
    run_cmd(["autogrid4", "-p", gpf.name, "-l", "receptor.glg"], cwd=prep, log_path=logs / "autogrid4.log")
    fld = prep / "receptor.maps.fld"
    if not fld.exists() or fld.stat().st_size == 0:
        raise PipelineError("AutoGrid did not create receptor.maps.fld")
    return fld


def run_adgpu(args: argparse.Namespace, prep: Path, adgpu: Path, logs: Path, fld: Path, ligand: Path) -> Path:
    log("Running AutoDock-GPU")
    adgpu.mkdir(parents=True, exist_ok=True)
    cmd = [
        "autodock_gpu_128wi",
        "--ffile",
        str(Path("../prep") / fld.name),
        "--lfile",
        str(Path("../prep") / ligand.name),
        "--nrun",
        str(args.nrun),
        "--gbest",
        "1",
        "--npdb",
        str(args.npdb),
        "--resnam",
        "ligand_adgpu",
    ]
    if args.seed:
        cmd.extend(["--seed", args.seed])
    if args.devnum:
        cmd.extend(["--devnum", str(args.devnum)])
    if args.nev:
        cmd.extend(["--nev", str(args.nev)])
    if args.no_autostop:
        cmd.extend(["--autostop", "0"])
    cmd.extend(split_extra_args(args.adgpu_extra, "--adgpu-extra"))
    run_cmd(cmd, cwd=adgpu, log_path=logs / "autodock_gpu.log")
    dlg = adgpu / "ligand_adgpu.dlg"
    if not dlg.exists():
        raise PipelineError("AutoDock-GPU did not create ligand_adgpu.dlg")
    return dlg


def split_dlg_run_best(dlg: Path, raw_dir: Path, clean_dir: Path) -> Path:
    log("Splitting AD_GPU DLG run-best poses")
    raw_dir.mkdir(parents=True, exist_ok=True)
    clean_dir.mkdir(parents=True, exist_ok=True)
    text = dlg.read_text(encoding="utf-8", errors="replace").splitlines()
    models: list[tuple[float, int, Path]] = []
    inside = False
    current: list[str] = []
    run_no: int | None = None
    model_no: int | None = None
    energy: float | None = None
    for line in text:
        if line.startswith("DOCKED: MODEL"):
            inside = True
            current = [line.replace("DOCKED: ", "", 1)]
            m = re.search(r"MODEL\s+(\d+)", current[-1])
            model_no = int(m.group(1)) if m else None
            run_no = None
            energy = None
            continue
        if not inside:
            continue
        if "Run =" in line:
            m = re.search(r"Run\s*=\s*(\d+)", line)
            if m:
                run_no = int(m.group(1))
        if "Estimated Free Energy of Binding" in line:
            m = re.search(r"=\s*([-+]?\d+(?:\.\d+)?)\s+kcal/mol", line)
            if m:
                energy = float(m.group(1))
        if line.startswith("DOCKED: ENDMDL"):
            current.append(line.replace("DOCKED: ", "", 1))
            if run_no is None:
                run_no = model_no or len(models) + 1
            if energy is None:
                energy = float("nan")
            raw = raw_dir / f"run{run_no:02d}_adgpu_{energy:+.2f}.pdbqt"
            raw.write_text("\n".join(current) + "\n", encoding="utf-8")
            clean_lines: list[str] = []
            for pose_line in current:
                fields = pose_line.split(None, 1)
                tag = fields[0] if fields else ""
                if tag in {"MODEL", "ENDMDL", "TER"}:
                    continue
                if tag == "USER":
                    pose_line = "REMARK ADGPU " + (fields[1] if len(fields) > 1 else "")
                clean_lines.append(pose_line.rstrip())
            clean = clean_dir / raw.name
            clean.write_text("\n".join(clean_lines) + "\n", encoding="utf-8")
            models.append((energy, run_no, clean))
            inside = False
            current = []
        elif line.startswith("DOCKED: "):
            current.append(line.replace("DOCKED: ", "", 1))
    if not models:
        raise PipelineError(f"No DOCKED models parsed from {dlg}")
    models.sort(key=lambda row: (row[0], row[1]))
    table = clean_dir / "adgpu_run_best_scores.tsv"
    with table.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["adgpu_rank", "run", "adgpu_energy_kcal_mol", "pose_file"])
        for rank, (energy, run_no, clean) in enumerate(models, 1):
            writer.writerow([rank, run_no, energy, clean.name])
    log(f"Parsed {len(models)} run-best poses")
    return table


def parse_gnina_scores(output: str) -> dict[str, str]:
    def grab(pattern: str) -> str:
        m = re.search(pattern, output)
        return m.group(1) if m else ""
    return {
        "gnina_affinity_kcal_mol": grab(r"Affinity:\s+([-+]?\d+(?:\.\d+)?)"),
        "cnnscore": grab(r"CNNscore:\s+([-+]?\d+(?:\.\d+)?)"),
        "cnnaffinity": grab(r"CNNaffinity:\s+([-+]?\d+(?:\.\d+)?)"),
        "cnnvariance": grab(r"CNNvariance:\s+([-+]?\d+(?:\.\d+)?)"),
    }


def score_pose_with_gnina(
    receptor_pdb: Path, pose: Path, log_path: Path, cwd: Path, gnina_extra: str = ""
) -> dict[str, str]:
    cmd = ["gnina", "--receptor", str(receptor_pdb), "--ligand", str(pose), "--score_only", "--cnn_scoring", "rescore"]
    cmd.extend(split_extra_args(gnina_extra, "--gnina-extra"))
    proc = run_cmd(
        cmd,
        cwd=cwd,
        log_path=log_path,
        check=False,
    )
    scores = parse_gnina_scores(proc.stdout)
    if proc.returncode != 0 or not scores.get("cnnscore"):
        scores["error"] = f"gnina_failed_exit_{proc.returncode}"
    else:
        scores["error"] = ""
    return scores


def gnina_rerank(
    args: argparse.Namespace,
    receptor_pdb_for_gnina: Path,
    adgpu: Path,
    gnina_dir: Path,
    logs: Path,
    score_table: Path,
) -> tuple[Path, Path | None, list[str]]:
    log("Running gnina scoring/rerank")
    warnings: list[str] = []
    best_pose = adgpu / "ligand_adgpu-best.pdbqt"
    best_log = logs / "gnina_score_adgpu_best.log"
    if best_pose.exists():
        score_pose_with_gnina(receptor_pdb_for_gnina, best_pose, best_log, cwd=gnina_dir.parent, gnina_extra=args.gnina_extra)
    else:
        warnings.append("AD_GPU best pose file ligand_adgpu-best.pdbqt not found; skipping best-pose score")

    rows_in = []
    with score_table.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows_in = list(reader)
    if args.rerank_limit.lower() != "all":
        try:
            limit = int(args.rerank_limit)
        except ValueError as exc:
            raise PipelineError("--rerank-limit must be 'all' or an integer") from exc
        rows_in = rows_in[:limit]

    clean_dir = score_table.parent
    rows_out: list[dict[str, str]] = []
    for i, row in enumerate(rows_in, 1):
        pose = clean_dir / row["pose_file"]
        log_path = clean_dir / f"{pose.stem}.gnina.log"
        scores = score_pose_with_gnina(receptor_pdb_for_gnina, pose, log_path, cwd=gnina_dir.parent, gnina_extra=args.gnina_extra)
        out_row = {
            "adgpu_rank": row["adgpu_rank"],
            "run": row["run"],
            "adgpu_energy_kcal_mol": row["adgpu_energy_kcal_mol"],
            "pose_file": str(pose.relative_to(gnina_dir.parent)),
            **scores,
        }
        rows_out.append(out_row)
        if scores.get("error"):
            warnings.append(f"gnina failed on {pose.name}: {scores['error']}")
        if i % 10 == 0:
            log(f"gnina scored {i}/{len(rows_in)} poses")

    if not rows_out:
        raise PipelineError("No poses selected for gnina rerank")
    fields = list(rows_out[0].keys())
    all_table = clean_dir / f"gnina_rerank_{len(rows_out)}poses.tsv"
    with all_table.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    ok_rows = [r for r in rows_out if not r.get("error") and r.get("cnnscore")]
    if not ok_rows:
        warnings.append("all gnina pose scores failed")
        return all_table, None, warnings

    def fval(row: dict[str, str], key: str, default: float = -999.0) -> float:
        try:
            return float(row.get(key) or default)
        except ValueError:
            return default

    by_cnn = sorted(ok_rows, key=lambda r: fval(r, "cnnscore"), reverse=True)
    by_cnn_table = clean_dir / f"gnina_rerank_{len(rows_out)}poses_by_cnnscore.tsv"
    with by_cnn_table.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(by_cnn)

    by_cnnaff = sorted(ok_rows, key=lambda r: fval(r, "cnnaffinity"), reverse=True)
    by_cnnaff_table = clean_dir / f"gnina_rerank_{len(rows_out)}poses_by_cnnaffinity.tsv"
    with by_cnnaff_table.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(by_cnnaff)

    return all_table, by_cnn_table, warnings


def write_summary(
    outdir: Path,
    args: argparse.Namespace,
    box: Box,
    receptor_method: str,
    warnings: list[str],
    rerank_table: Path | None,
) -> None:
    summary = outdir / "SUMMARY.md"
    best_cnn_line = ""
    if rerank_table and rerank_table.exists():
        with rerank_table.open("r", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            first = next(reader, None)
        if first:
            best_cnn_line = (
                f"adgpu_rank={first.get('adgpu_rank')} run={first.get('run')} "
                f"adgpu_energy={first.get('adgpu_energy_kcal_mol')} "
                f"gnina_affinity={first.get('gnina_affinity_kcal_mol')} "
                f"CNNscore={first.get('cnnscore')} CNNaffinity={first.get('cnnaffinity')} "
                f"CNNvariance={first.get('cnnvariance')} pose={first.get('pose_file')}"
            )
    lines = [
        "# AutoDock-GPU + gnina pipeline summary",
        "",
        f"Run directory: `{outdir}`",
        f"Receptor input: `{Path(args.receptor).resolve()}`",
        f"Ligand input: `{Path(args.ligand).resolve()}`",
        f"Box input: `{Path(args.box).resolve()}`",
        "",
        "## Parsed box",
        "",
        f"- npts: `{box.npts[0]} {box.npts[1]} {box.npts[2]}`",
        f"- spacing: `{box.spacing}`",
        f"- gridcenter: `{box.center[0]:.3f} {box.center[1]:.3f} {box.center[2]:.3f}`",
        f"- size: `{box.size}`",
        f"- parse notes: `{'; '.join(box.source_notes or [])}`",
        "",
        "## Parameters",
        "",
        f"- AutoDock-GPU nrun: `{args.nrun}`",
        f"- AutoDock-GPU npdb: `{args.npdb}`",
        f"- AD_GPU seed: `{args.seed or 'default/time-based'}`",
        f"- AD_GPU extra args: `{args.adgpu_extra or 'none'}`",
        f"- gnina extra args: `{args.gnina_extra or 'none'}`",
        f"- rerank limit: `{args.rerank_limit}`",
        f"- receptor PDBQT method: `{receptor_method}`",
        "",
        "## Outputs",
        "",
        "- AutoGrid maps: `prep/receptor.maps.fld`",
        "- AD_GPU DLG: `adgpu/ligand_adgpu.dlg`",
        "- AD_GPU best pose: `adgpu/ligand_adgpu-best.pdbqt`",
        "- AD_GPU run-best poses: `gnina/dlg_run_best_clean/run*_adgpu_*.pdbqt`",
        "- gnina rerank by CNNscore: `gnina/dlg_run_best_clean/*_by_cnnscore.tsv`",
        "",
        "## Best by gnina CNNscore",
        "",
        f"```text\n{best_cnn_line or 'not available'}\n```",
        "",
        "## Warnings / notes",
        "",
    ]
    if warnings:
        lines.extend(f"- {w}" for w in warnings)
    else:
        lines.append("- none")
    lines.append("")
    summary.write_text("\n".join(lines), encoding="utf-8")


def format_box_template(
    center: tuple[float, float, float],
    npts: tuple[int, int, int],
    spacing: float,
) -> str:
    return textwrap.dedent(
        f"""\
        *********AutoDock Grid Option*********
        npts {npts[0]} {npts[1]} {npts[2]} # num. grid points in xyz
        spacing {spacing:g} # spacing (A)
        gridcenter {center[0]:.3f} {center[1]:.3f} {center[2]:.3f} # xyz-coordinates or auto
        """
    )


def build_box_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="dock box",
        description="Write an editable AutoDock-style box.txt template for the AD_GPU + gnina pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("output", nargs="?", default="box.txt", help="template path to write")
    p.add_argument(
        "--center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        default=(20.813, 10.963, 22.132),
        help="initial AutoDock gridcenter",
    )
    p.add_argument(
        "--npts",
        nargs=3,
        type=int,
        metavar=("X", "Y", "Z"),
        default=(63, 40, 43),
        help="initial AutoDock grid point counts",
    )
    p.add_argument(
        "--size",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        help="optional Vina-style size in Angstrom; if set, npts=ceil(size/spacing)",
    )
    p.add_argument("--spacing", type=float, default=0.375, help="AutoGrid spacing")
    p.add_argument("--force", action="store_true", help="overwrite output if it already exists")
    return p


def box_main(argv: Sequence[str]) -> int:
    parser = build_box_arg_parser()
    args = parser.parse_args(argv)
    center = tuple(args.center)
    spacing = args.spacing
    if spacing <= 0:
        parser.error("--spacing must be positive")
    if args.size:
        if any(v <= 0 for v in args.size):
            parser.error("--size values must be positive")
        npts = tuple(max(1, int(math.ceil(v / spacing))) for v in args.size)
    else:
        npts = tuple(args.npts)
    if any(v <= 0 for v in npts):
        parser.error("--npts values must be positive")

    path = Path(args.output).expanduser()
    if path.exists() and not args.force:
        parser.error(f"output already exists; use --force to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(format_box_template(center, npts, spacing), encoding="utf-8")
    print(f"Wrote AutoDock box template: {path.resolve()}")
    print(f"Edit npts/spacing/gridcenter, then run: dock receptor.pdb ligand.mol2 {path}")
    return 0


EXTRA_HELP = r"""

Box template shortcut
=====================

After installing the bash alias, create an editable box file anywhere with:
  dock box

Examples:
  dock box
  dock box my_box.txt --center 20.813 10.963 22.132 --size 23.8 15.3 16.1

Extra passthrough quick reference
=================================

Anything in --adgpu-extra is appended verbatim to autodock_gpu_128wi.
Anything in --gnina-extra is appended verbatim to each gnina scoring command:
  gnina --receptor ... --ligand ... --score_only --cnn_scoring rescore [--gnina-extra ...]

Quote the whole extra string:
  --adgpu-extra "--lsit 500 --psize 200 --rmstol 1.5"
  --gnina-extra "--cnn_rotation 24 --cpu 4"

AutoDock-GPU common passthrough options
---------------------------------------
Search / sampling:
  --heuristics 0|1       ligand-based automatic search/evaluation heuristic
  --heurmax N            cap for heuristic max evaluations
  --autostop 0|1         automatic convergence stop
  --asfreq N             AutoStop check frequency in generations
  --nrun N               number of LGA runs (also exposed as script --nrun)
  --nev N                max score evaluations per LGA run (also exposed as script --nev)
  --ngen N               max generations per LGA run
  --lsmet ad|sw          local-search method; ad=ADADELTA, sw=Solis-Wets
  --lsit N               local-search iterations
  --psize N              population size
  --mrat PCT             mutation rate
  --crat PCT             crossover rate
  --lsrat PCT            local-search rate
  --trat PCT             tournament selection rate
  --dmov A               max LGA movement delta in Angstrom
  --dang DEG             max LGA angle delta
  --stopstd X            AutoStop energy std tolerance, kcal/mol
  --initswgens N         initial Solis-Wets generations before selected lsmet

Output / analysis:
  --gbest 0|1            output single best pose as PDBQT (script forces 1)
  --npdb N               output pose PDBQT files from populations (script exposes --npdb)
  --xmloutput 0|1        write XML output
  --dlgoutput 0|1        write DLG output
  --dlg2stdout 0|1       write DLG to stdout
  --clustering 0|1       include clustering in DLG/XML
  --rmstol A             RMSD clustering tolerance
  --contact_analysis 0|1 distance-based contact analysis
  --output-cluster-poses N output up to N poses per cluster

Setup / reproducibility:
  --devnum N             CUDA/OpenCL device; AutoDock-GPU counts from 1 (script exposes --devnum)
  --seed S               one to three comma-separated seeds (script exposes --seed)
  --loadxml FILE         load initial population from XML

Scoring / force-field tweaks:
  --derivtype SPEC       derivative atom types, e.g. C1,C2,C3=C/S4=S/H5=HD
  --modpair SPEC         modify vdW pair params
  --ubmod 0|1|2          unbound model
  --smooth A             vdW smoothing parameter
  --elecmindist A        minimum electrostatic distance
  --modqp 0|1            modified QASP vs AD4 original

Inputs not usually needed in --adgpu-extra because script sets them:
  --lfile, --ffile, --resnam, --nrun, --npdb, --devnum, --seed, --nev, --autostop

Examples:
  --adgpu-extra "--lsit 500 --psize 200 --rmstol 1.5"
  --adgpu-extra "--autostop 0 --nev 2500000 --ngen 42000"
  --adgpu-extra "--contact_analysis 1 --clustering 1"


gnina common passthrough options for this rescoring workflow
-----------------------------------------------------------
CNN scoring:
  --cnn MODEL            built-in CNN model name
  --cnn_model FILE       custom torch CNN model file
  --cnn_rotation N       evaluate multiple rotations of pose; slower, more GPU-visible
  --cnn_scoring MODE     none|rescore|refinement|metrorescore|metrorefine
                         script already sets rescore; extra can override only if gnina honors last flag
  --cnn_verbose          verbose CNN debug output
  --cnn_center_x X       explicit CNN center X
  --cnn_center_y Y       explicit CNN center Y
  --cnn_center_z Z       explicit CNN center Z
  --cnn_empirical_weight W
  --cnn_mix_emp_energy
  --cnn_mix_emp_force

Device / performance:
  --device N             gnina GPU device; gnina counts from 0
  --no_gpu               force CPU; do not use if you want GPU
  --cpu N                CPU threads for non-GPU parts

Pose/result behavior mostly useful outside score_only:
  --pose_sort_order MODE CNNscore|CNNaffinity|Energy
  --min_rmsd_filter A
  --num_modes N
  --exhaustiveness N

Inputs not usually needed in --gnina-extra because script sets them:
  --receptor, --ligand, --score_only, --cnn_scoring rescore

Examples:
  --gnina-extra "--cnn_rotation 24"
  --gnina-extra "--device 0 --cpu 4"
  --gnina-extra "--cnn dense --cnn_rotation 12"

For exact upstream help:
  autodock_gpu_128wi --help
  gnina --help

This script also supports:
  --help-full            print script help plus captured upstream help snippets
"""


class PipelineArgumentParser(argparse.ArgumentParser):
    def format_help(self) -> str:
        return super().format_help() + EXTRA_HELP

def build_arg_parser() -> argparse.ArgumentParser:
    p = PipelineArgumentParser(
        description="Run AutoDock-GPU then gnina rerank for one receptor/ligand/box.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--help-full", action="store_true", help="print script help plus captured upstream AutoDock-GPU/gnina help and exit")
    p.add_argument("receptor", nargs="?", help="receptor PDB file for gnina and receptor preparation")
    p.add_argument("ligand", nargs="?", help="ligand MOL2/SDF/etc. file; MOL2 is expected by user workflow")
    p.add_argument("box", nargs="?", help="box.txt containing AutoGrid npts/spacing/gridcenter or Vina center/size")
    p.add_argument("-o", "--outdir", help="output run directory; default is timestamped next to ligand")
    p.add_argument("--force", action="store_true", help="delete output directory if it already exists")
    p.add_argument("--nrun", type=int, default=50, help="AutoDock-GPU number of LGA runs")
    p.add_argument("--npdb", type=int, default=20, help="AutoDock-GPU population pose output count")
    p.add_argument("--rerank-limit", default="all", help="number of AD_GPU run-best poses to score with gnina, or 'all'")
    p.add_argument(
        "--adgpu-extra",
        default="",
        help="extra raw arguments appended to autodock_gpu_128wi, e.g. '--lsit 500 --psize 200'",
    )
    p.add_argument(
        "--gnina-extra",
        default="",
        help="extra raw arguments appended to every gnina scoring command, e.g. '--cnn_rotation 24'",
    )
    p.add_argument("--seed", help="AutoDock-GPU seed string, e.g. 123 or 1,2,3")
    p.add_argument("--devnum", type=int, help="AutoDock-GPU CUDA device number; AD_GPU counts from 1")
    p.add_argument("--nev", type=int, help="AutoDock-GPU max score evaluations per LGA run")
    p.add_argument("--no-autostop", action="store_true", help="disable AD_GPU autostop")
    p.add_argument("--spacing-default", type=float, default=0.375, help="spacing when box file omits spacing")
    p.add_argument(
        "--no-prefer-existing-pdbqt",
        dest="prefer_existing_pdbqt",
        action="store_false",
        help="do not prefer a same-stem receptor .pdbqt when present",
    )
    p.set_defaults(prefer_existing_pdbqt=True)
    p.add_argument(
        "--no-obabel-receptor-fallback",
        dest="allow_obabel_receptor",
        action="store_false",
        help="do not use Open Babel receptor PDBQT fallback when Meeko fails and no same-stem PDBQT is present",
    )
    p.set_defaults(allow_obabel_receptor=True)
    return p


def print_help_full(parser: argparse.ArgumentParser) -> None:
    print(parser.format_help())
    for tool in ["autodock_gpu_128wi", "gnina"]:
        path = which(tool)
        print("\n" + "=" * 100)
        print(f"Upstream help: {tool}")
        print("=" * 100)
        if not path:
            print(f"{tool} not found in PATH")
            continue
        try:
            proc = subprocess.run([tool, "--help"], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=20)
            print(proc.stdout)
        except Exception as exc:
            print(f"Could not capture {tool} --help: {exc}")


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "box":
        return box_main(argv[1:])

    parser = build_arg_parser()
    args = parser.parse_args(argv)
    if getattr(args, "help_full", False):
        print_help_full(parser)
        return 0
    if not args.receptor or not args.ligand or not args.box:
        parser.error("receptor, ligand, and box are required unless using -h/--help or --help-full")

    receptor = Path(args.receptor).expanduser().resolve()
    ligand = Path(args.ligand).expanduser().resolve()
    box_path = Path(args.box).expanduser().resolve()
    for label, path in [("receptor", receptor), ("ligand", ligand), ("box", box_path)]:
        if not path.exists():
            parser.error(f"{label} file does not exist: {path}")
    if args.nrun <= 0 or args.npdb < 0:
        parser.error("--nrun must be >0 and --npdb must be >=0")

    require_tool("mk_prepare_ligand.py")
    require_tool("mk_prepare_receptor.py")
    require_tool("autogrid4")
    require_tool("autodock_gpu_128wi")
    require_tool("gnina")

    box = parse_box(box_path, default_spacing=args.spacing_default)

    if args.outdir:
        outdir = Path(args.outdir).expanduser().resolve()
    else:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        outdir = ligand.parent / f"adgpu_gnina_{ligand.stem}_{stamp}"
    if outdir.exists():
        if args.force:
            shutil.rmtree(outdir)
        else:
            raise PipelineError(f"Output directory already exists; use --force or choose another -o: {outdir}")
    prep = outdir / "prep"
    adgpu = outdir / "adgpu"
    gnina_dir = outdir / "gnina"
    logs = outdir / "logs"
    for d in [prep, adgpu, gnina_dir, logs]:
        d.mkdir(parents=True, exist_ok=True)

    warnings: list[str] = []
    start = time.time()
    log(f"Output directory: {outdir}")
    try:
        ligand_pdbqt = prepare_ligand(args, prep, logs, warnings)
        receptor_pdbqt, receptor_method = prepare_receptor(args, prep, logs, box, warnings)
        gpf = write_gpf(prep, receptor_pdbqt, ligand_pdbqt, box)
        fld = run_autogrid(prep, logs, gpf)
        dlg = run_adgpu(args, prep, adgpu, logs, fld, ligand_pdbqt)
        raw_dir = gnina_dir / "dlg_run_best_raw"
        clean_dir = gnina_dir / "dlg_run_best_clean"
        score_table = split_dlg_run_best(dlg, raw_dir, clean_dir)
        all_table, by_cnn_table, gnina_warnings = gnina_rerank(
            args, prep / "receptor_input.pdb", adgpu, gnina_dir, logs, score_table
        )
        warnings.extend(gnina_warnings)
        write_summary(outdir, args, box, receptor_method, warnings, by_cnn_table)
    except Exception as exc:
        error_file = outdir / "ERROR.txt"
        error_file.write_text(str(exc) + "\n", encoding="utf-8")
        raise

    elapsed = time.time() - start
    log(f"Complete in {elapsed:.1f}s")
    log(f"Summary: {outdir / 'SUMMARY.md'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except PipelineError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
