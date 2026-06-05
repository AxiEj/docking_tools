# AD_GPU + gnina 一键工作流

标签：`AutoDock-GPU` · `gnina` · `CNNscore` · `CNN rescoring` · `AutoGrid4` · `Meeko` · `Open Babel` · `molecular docking` · `virtual screening` · `GPU docking`

定位：**AutoDock-GPU 负责快速采样，gnina 负责 CNN 重打分排序**。gnina 1.0 论文报告 CNN scoring 重打分在 redocking / cross-docking 的 Top1 pose-ranking 指标上优于 AutoDock Vina scoring；AutoDock-GPU 论文报告其 GPU/梯度局部搜索实现可显著加速 AutoDock4 类搜索。这个流程利用两者优点，但不宣称所有体系无条件最高精度。

文献：

- McNutt A. T. et al. *GNINA 1.0: molecular docking with deep learning.* Journal of Cheminformatics 13, 43 (2021). DOI: [10.1186/s13321-021-00522-2](https://doi.org/10.1186/s13321-021-00522-2)
- Santos-Martins D. et al. *Accelerating AutoDock4 with GPUs and Gradient-Based Local Search.* Journal of Chemical Theory and Computation 17, 1060-1073 (2021). DOI: [10.1021/acs.jctc.0c01006](https://doi.org/10.1021/acs.jctc.0c01006)

详细安装/下载 excluded binaries 的步骤见：[INSTALL_AND_RUN.md](INSTALL_AND_RUN.md)。

主脚本：

```bash
/home/axie/docking_tools/run_adgpu_gnina.sh
```

已在 `~/.bashrc` 里配置后，可以在任何目录直接用：

```bash
dock receptor.pdb ligand.mol2 box.txt
```


## 查看所有参数

```bash
dock -h
```

`-h` 会显示脚本参数，以及 `--adgpu-extra` / `--gnina-extra` 常用可透传参数清单。

如果要看底层工具原始帮助：

```bash
dock --help-full
```

也可以单独看：

```bash
autodock_gpu_128wi --help
gnina --help
```

## 最小输入

```text
受体 PDB:   receptor.pdb
配体 MOL2:  ligand.mol2
盒子文件:   box.txt
```

`box.txt` 支持优先解析：

```text
npts 49 75 43
spacing 0.375
gridcenter 54.809 59.493 30.181
```

如果没有 AutoGrid 格式，也会尝试解析 PyMOL `getbox sele` 这种 Vina 一行格式：

```text
--center_x ... --center_y ... --center_z ... --size_x ... --size_y ... --size_z ...
```

也支持 LeDock：`Binding pocket` 后面三行 min/max。

## 生成 box 模板

在任意目录生成可编辑模板：

```bash
dock box
```

默认写出当前目录的 `box.txt`。默认只启用 AutoDock Grid 三行；Vina/LeDock/BoxCode 作为注释参考保留，解析器会忽略 `#` 注释：

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

你也可以指定文件名和初始参数：

```bash
dock box my_box.txt --center 54.809 59.493 30.181 --npts 49 75 43
```

如果你只有 Vina 的 size，也可以让模板先自动换算 `npts=ceil(size/spacing)`：

```bash
dock box my_box.txt --center 54.809 59.493 30.181 --size 18.7 28.5 16.3
```

之后只需要把 AutoDock Grid 的 `npts`、`spacing`、`gridcenter` 换成目标口袋参数，就可以直接用于 `dock receptor.pdb ligand.mol2 box.txt`。

## 推荐命令

```bash
dock \
  /path/to/receptor.pdb \
  /path/to/ligand.mol2 \
  /path/to/box.txt \
  -o /path/to/output_run \
  --nrun 50 \
  --npdb 20 \
  --rerank-limit all
```

如果想可重复：

```bash
--seed 123
```


## 高级透传参数

如果脚本没显式暴露某个 AutoDock-GPU 或 gnina 参数，可以用：

```bash
--adgpu-extra "--lsit 500 --psize 200 --rmstol 1.5"
--gnina-extra "--cnn_rotation 24 --cpu 4"
```

注意：

- `--adgpu-extra` 会追加到 `autodock_gpu_128wi` 命令末尾。
- `--gnina-extra` 会追加到每一次 `gnina --score_only --cnn_scoring rescore` 命令末尾。
- 如果 extra 里面有空格，整串必须用引号包住。
- 透传参数写错时，底层工具会失败，日志在 `logs/` 或 `gnina/dlg_run_best_clean/*.gnina.log`。

## 输出目录结构

```text
output_run/
  SUMMARY.md
  prep/
    receptor.pdbqt
    ligand.pdbqt
    receptor.gpf
    receptor.maps.fld
    receptor.*.map
  adgpu/
    ligand_adgpu.dlg
    ligand_adgpu.xml
    ligand_adgpu-best.pdbqt
  gnina/
    dlg_run_best_clean/
      run*_adgpu_*.pdbqt
      adgpu_run_best_scores.tsv
      gnina_rerank_*poses.tsv
      gnina_rerank_*poses_by_cnnscore.tsv
      gnina_rerank_*poses_by_cnnaffinity.tsv
  logs/
```

最重要结果表：

```bash
gnina/dlg_run_best_clean/*_by_cnnscore.tsv
```

第一行就是 gnina CNNscore 排名最高的 AD_GPU pose。

## 鲁棒性 / fallback

脚本会自动：

1. 检查 `gnina`, `autodock_gpu_128wi`, `autogrid4`, `mk_prepare_ligand.py`, `mk_prepare_receptor.py` 是否可用。
2. 自动解析 box 文件；缺 `npts` 但有 Vina size 时，会用 `ceil(size/spacing)` 推导。
3. ligand 优先用 Meeko 准备；失败则尝试 Open Babel 转 SDF 后再 Meeko。
4. receptor 优先使用同名已准备好的 `.pdbqt`，例如输入 `receptor.pdb` 时如果旁边有 `receptor.pdbqt`，会直接复用它。
5. 如果没有同名 `.pdbqt`，尝试 Meeko/ProDy 和 Meeko `--read_pdb`。
6. 如果 Meeko 失败，默认允许 Open Babel receptor fallback；这种 fallback 会在 `SUMMARY.md` 里警告，因为 receptor charges 可能不如人工/Meeko/MGLTools 准备可靠。
7. AutoGrid GPF 手写生成，不依赖 Meeko 生成 GPF。
8. 自动拆 AD_GPU DLG 的每个 run-best pose，并清理 `MODEL/ENDMDL/TER` 以适配 gnina。
9. gnina scoring 某个 pose 失败时记录 error，不直接吞掉日志。

## 已通过烟测

测试命令：

```bash
cd /home/axie/docking_tools
./run_adgpu_gnina.sh \
  test/fixed_noligand.pdb \
  test/ligand_954_A501.mol2 \
  test/box.txt \
  -o test/script_smoke_run \
  --force \
  --nrun 5 \
  --npdb 0 \
  --seed 123 \
  --rerank-limit all
```

验证：maps、DLG、best pose、5 行 gnina rerank 表均已生成。
