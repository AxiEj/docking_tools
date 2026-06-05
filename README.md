# AutoDock-GPU + gnina CNN Rescoring Workflow

[中文版](README.zh-CN.md) | [Install and run guide](INSTALL_AND_RUN.md) | [Pipeline details](README_ADGPU_GNINA_PIPELINE.md)

This repository provides a local one-step molecular docking workflow: **AutoDock-GPU / AutoGrid4** performs AutoDock4-style sampling, then **gnina CNNscore / CNNaffinity** rescored and ranks the AutoDock-GPU poses.

## Tags / Keywords

`AutoDock-GPU` · `gnina` · `gnina CNN rescoring` · `CNNscore` · `AutoGrid4` · `Meeko` · `Open Babel` · `molecular docking` · `protein-ligand docking` · `virtual screening` · `GPU docking` · `CUDA` · `drug discovery` · `structure-based drug design`

## Why this workflow is useful

- **AutoDock-GPU** accelerates AutoDock4-style search/local optimization on GPUs and keeps the familiar grid/GPF/DLG workflow.
- **gnina CNN rescoring** is literature-backed: the gnina 1.0 paper reports that CNN scoring improved Top1 pose-ranking over AutoDock Vina scoring in redocking and cross-docking benchmarks; with a defined binding pocket, Top1 improved from 58% to 73% for redocking and from 27% to 37% for cross-docking.
- **This repository combines both ideas**: AutoDock-GPU generates run-best poses, the script cleans DLG-derived poses, then gnina runs `--score_only --cnn_scoring rescore` to produce CNNscore-ranked tables for downstream inspection and filtering.

> This is not a claim of universally highest accuracy for every target. The literature support is specifically for gnina CNN rescoring in its benchmark context. Real projects should still use controls, repeat docking, check protonation/conformations, and validate with follow-up methods such as MD, MM-GBSA, or experiments.

## Literature backing

- McNutt A. T. et al. **GNINA 1.0: molecular docking with deep learning.** *Journal of Cheminformatics* 13, 43 (2021). DOI: [10.1186/s13321-021-00522-2](https://doi.org/10.1186/s13321-021-00522-2)
- Santos-Martins D. et al. **Accelerating AutoDock4 with GPUs and Gradient-Based Local Search.** *Journal of Chemical Theory and Computation* 17, 1060-1073 (2021). DOI: [10.1021/acs.jctc.0c01006](https://doi.org/10.1021/acs.jctc.0c01006)

## Quick use

First install dependencies and recreate the Git-ignored `downloads/` and `bin/` directories by following [INSTALL_AND_RUN.md](INSTALL_AND_RUN.md).

After configuring the `dock` alias:

```bash
dock box
dock receptor.pdb ligand.mol2 box.txt
```

## Box template

`dock box` keeps only the AutoDock Grid block active by default. Vina / LeDock / BoxCode examples are retained as commented references, and the parser ignores `#` comments:

```text
# !!! MUST CHANGE FOR YOUR REAL SYSTEM !!!
# !!! The following three AutoDock Grid lines are template defaults; replace them before real docking. !!!

*********AutoDock Grid Option*********
npts 49 75 43 # num. grid points in xyz
spacing 0.375 # spacing (A)
gridcenter 54.809 59.493 30.181 # xyz-coordinates or auto
```

Before a real run, edit only the active `npts`, `spacing`, and `gridcenter` values for your target pocket.

## Repository boundary

This repository tracks only the Python/Bash workflow and documentation. Local test inputs are ignored by Git. It intentionally does not commit:

- local installed binaries: `bin/`
- downloaded release/cache files: `downloads/`
- local test inputs: `test/`
- runtime logs and docking outputs
- OMX / local agent state

This avoids pushing 2GB-scale gnina CUDA binaries and machine-specific artifacts to GitHub. See [INSTALL_AND_RUN.md](INSTALL_AND_RUN.md) for the full setup process.
