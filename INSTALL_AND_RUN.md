# Install and run `docking_tools`

This repository contains the reusable workflow code for a one-step AutoDock-GPU + gnina run. It intentionally does **not** commit local binary installs, downloaded release files, logs, or docking outputs. After cloning, rebuild `downloads/` and `bin/` locally with the steps below.

Tested local target stack:

- Linux / WSL2 with an NVIDIA CUDA-capable GPU
- Python 3.11
- gnina `v1.3.2` CUDA 12.8 build
- AutoDock-GPU `v1.6` CUDA 12 128wi build
- AutoGrid4 `4.2.6` from the Bioconda package
- Meeko command-line tools: `mk_prepare_ligand.py`, `mk_prepare_receptor.py`
- Open Babel for fallback conversion

## 1. Clone

```bash
git clone https://github.com/AxiEj/docking_tools.git
cd docking_tools
```

## 2. Create the Python/tool environment

Use Conda/Mamba if available, because Open Babel is easier to install from conda-forge.

```bash
conda create -n docking_tools python=3.11 -y
conda activate docking_tools
conda install -c conda-forge openbabel numpy scipy biopython -y
python -m pip install meeko prody
```

Check the preparation tools:

```bash
mk_prepare_ligand.py --help | head
mk_prepare_receptor.py --help | head
obabel -V
```

## 3. Download the excluded binaries

The local `downloads/` folder is ignored by Git. Recreate it with the verified URLs:

```bash
mkdir -p downloads bin

curl -L \
  -o downloads/gnina.1.3.2.cuda12.8 \
  https://github.com/gnina/gnina/releases/download/v1.3.2/gnina.1.3.2.cuda12.8

curl -L \
  -o downloads/adgpu-v1.6_linux_x64_cuda12_128wi \
  https://github.com/ccsb-scripps/AutoDock-GPU/releases/download/v1.6/adgpu-v1.6_linux_x64_cuda12_128wi

curl -L \
  -o downloads/adgpu_analysis-v1.6_linux_x64 \
  https://github.com/ccsb-scripps/AutoDock-GPU/releases/download/v1.6/adgpu_analysis-v1.6_linux_x64

curl -L \
  -o downloads/autogrid-4.2.6-h9948957_4.tar.bz2 \
  https://conda.anaconda.org/bioconda/linux-64/autogrid-4.2.6-h9948957_4.tar.bz2
```

Optional integrity checks from the original validated install:

```bash
sha256sum downloads/gnina.1.3.2.cuda12.8 \
  downloads/adgpu-v1.6_linux_x64_cuda12_128wi \
  downloads/adgpu_analysis-v1.6_linux_x64 \
  downloads/autogrid-4.2.6-h9948957_4.tar.bz2
```

Expected hashes:

```text
714ef2928a22c7b20680ccfeade3c2a652a4e452a1f64343da541434528ed9cd  downloads/gnina.1.3.2.cuda12.8
04870ff6cf3451f97540748fa2d20f2b261c9b30bc312ecfc6adb09246dbb276  downloads/adgpu-v1.6_linux_x64_cuda12_128wi
86fe70ea33a8b5e37a7f99b961fc413b40be30cbf69693144a8677dff84bd6f4  downloads/adgpu_analysis-v1.6_linux_x64
d4f0ecaa725ff8c1b4a7c772154e97c96f99c54672f346cbd477f1fdc5a326a1  downloads/autogrid-4.2.6-h9948957_4.tar.bz2
```

## 4. Install binaries into `bin/`

```bash
cp downloads/adgpu-v1.6_linux_x64_cuda12_128wi bin/autodock_gpu_128wi
cp downloads/adgpu_analysis-v1.6_linux_x64 bin/adgpu_analysis
cp downloads/gnina.1.3.2.cuda12.8 bin/gnina.1.3.2.cuda12.8
chmod +x bin/autodock_gpu_128wi bin/adgpu_analysis bin/gnina.1.3.2.cuda12.8

tmpdir="$(mktemp -d)"
tar -xjf downloads/autogrid-4.2.6-h9948957_4.tar.bz2 -C "$tmpdir"
cp "$tmpdir/bin/autogrid4" bin/autogrid4.bioconda-4.2.6-h9948957_4
chmod +x bin/autogrid4.bioconda-4.2.6-h9948957_4
rm -rf "$tmpdir"
```

Create the small wrappers used by the pipeline:

```bash
cat > bin/autogrid4 <<'SH'
#!/usr/bin/env bash
set -euo pipefail
REAL_AUTOGRID4="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/autogrid4.bioconda-4.2.6-h9948957_4"
if [[ $# -eq 0 ]]; then
  exec "$REAL_AUTOGRID4" -h
fi
exec "$REAL_AUTOGRID4" "$@"
SH
chmod +x bin/autogrid4

cat > bin/gnina <<'SH'
#!/usr/bin/env bash
set -euo pipefail
INSTALL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# gnina v1.3.2 CUDA build needs cuDNN 9. If libcudnn.so.9 is not already visible,
# set GNINA_CUDNN_LIB to a directory containing libcudnn.so.9 before running gnina/dock.
if [[ -n "${GNINA_CUDNN_LIB:-}" && -d "$GNINA_CUDNN_LIB" ]]; then
  export LD_LIBRARY_PATH="$GNINA_CUDNN_LIB${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi
exec "$INSTALL_DIR/gnina.1.3.2.cuda12.8" "$@"
SH
chmod +x bin/gnina
```

## 5. Add commands to `PATH`

For the current shell:

```bash
export PATH="$PWD/bin:$PATH"
```

For all future shells, add this to `~/.bashrc` from the repository root:

```bash
cat >> ~/.bashrc <<EOF_BASHRC

# docking_tools
export PATH="$PWD/bin:\$PATH"
alias dock='bash $PWD/run_adgpu_gnina.sh'
EOF_BASHRC
```

Reload:

```bash
source ~/.bashrc
```

Check command discovery:

```bash
which dock || type dock
which gnina autodock_gpu_128wi autogrid4 mk_prepare_ligand.py mk_prepare_receptor.py obabel
```

If `gnina --help` fails with `libcudnn.so.9`, install cuDNN 9 or point `GNINA_CUDNN_LIB` at an existing cuDNN 9 library directory:

```bash
export GNINA_CUDNN_LIB=/path/to/directory/containing/libcudnn.so.9
gnina --help | head
```

## 6. Generate and edit a box file

`dock box` writes an AutoDock Grid template. Vina one-line, LeDock min/max, and BoxCode examples are also included, but commented out and ignored by the parser, so users only need to edit one active format by default.

```bash
dock box
```

Default template starts with a required-change marker:

```text
# !!! MUST CHANGE FOR YOUR REAL SYSTEM !!!
# !!! 下面 AutoDock Grid 三行只是模板默认值；真实对接前必须换成你的口袋参数 !!!

*********AutoDock Grid Option*********
npts 49 75 43 # num. grid points in xyz
spacing 0.375 # spacing (A)
gridcenter 54.809 59.493 30.181 # xyz-coordinates or auto

# *********AutoDock Vina Binding Pocket*********
# --center_x 54.8 --center_y 59.5 --center_z 30.2 --size_x 18.7 --size_y 28.5 --size_z 16.3

# *********LeDock Binding Pocket*********
# Binding pocket
# 45.4 64.2
# 45.2 73.7
# 22.0 38.3

# BoxCode(box_8332) = showbox 45.4, 64.2, 45.2, 73.7, 22.0, 38.3
```

Before a real docking run, replace the active AutoDock Grid values: `npts`, `spacing`, and `gridcenter`. Leave the commented Vina/LeDock examples alone unless you deliberately switch formats.

## 7. Optional smoke test with your own files

The repository does not upload `test/`; local test inputs and run outputs are intentionally ignored. If you have local receptor/ligand/box files, use a small run like this to validate the installation:

```bash
dock \
  /path/to/receptor.pdb \
  /path/to/ligand.mol2 \
  /path/to/box.txt \
  -o /path/to/smoke_run \
  --force \
  --nrun 5 \
  --npdb 0 \
  --seed 123 \
  --rerank-limit all
```

Expected important outputs:

```text
/path/to/smoke_run/SUMMARY.md
/path/to/smoke_run/prep/receptor.maps.fld
/path/to/smoke_run/adgpu/ligand_adgpu.dlg
/path/to/smoke_run/adgpu/ligand_adgpu-best.pdbqt
/path/to/smoke_run/gnina/dlg_run_best_clean/*_by_cnnscore.tsv
```

Run outputs are ignored by Git.

## 8. Real run pattern

```bash
dock \
  /path/to/receptor.pdb \
  /path/to/ligand.mol2 \
  /path/to/box.txt \
  -o /path/to/output_run \
  --nrun 50 \
  --npdb 20 \
  --rerank-limit all \
  --seed 123
```

Useful options:

- `--adgpu-extra "..."`: append raw arguments to `autodock_gpu_128wi`.
- `--gnina-extra "..."`: append raw arguments to each gnina score-only call.
- `--no-prefer-existing-pdbqt`: do not reuse same-stem receptor `.pdbqt`.
- `--no-obabel-receptor-fallback`: fail instead of using Open Babel receptor fallback.
- `--help-full`: show script help plus upstream AutoDock-GPU/gnina help snippets.

## 9. Troubleshooting

- `Required tool not found in PATH`: activate the conda env and ensure `bin/` is on `PATH`.
- `libcudnn.so.9` missing: set `GNINA_CUDNN_LIB` or install a CUDA/cuDNN runtime compatible with gnina `v1.3.2`.
- receptor preparation fails: put a curated same-stem receptor PDBQT next to the receptor PDB, for example `receptor.pdb` + `receptor.pdbqt`.
- box parsing fails: keep the Vina, AutoDock Grid, and LeDock values mutually consistent. The parser accepts AutoDock `npts/spacing/gridcenter`, Vina one-line `--center_* --size_*`, and LeDock `Binding pocket` min/max lines.
