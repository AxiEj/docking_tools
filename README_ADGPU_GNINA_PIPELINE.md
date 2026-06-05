# AD_GPU + gnina 一键工作流

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
npts 63 40 43
spacing 0.375
gridcenter 20.813 10.963 22.132
```

如果没有 AutoGrid 格式，也会尝试解析 Vina 的：

```text
--center_x ... --center_y ... --center_z ... --size_x ... --size_y ... --size_z ...
```

## 生成 box 模板

在任意目录生成可编辑模板：

```bash
dock box
```

默认写出当前目录的 `box.txt`，内容是 AutoDock 风格：

```text
*********AutoDock Grid Option*********
npts 63 40 43 # num. grid points in xyz
spacing 0.375 # spacing (A)
gridcenter 20.813 10.963 22.132 # xyz-coordinates or auto
```

你也可以指定文件名和初始参数：

```bash
dock box my_box.txt --center 20.813 10.963 22.132 --npts 63 40 43
```

如果你只有 Vina 的 size，也可以让模板先自动换算 `npts=ceil(size/spacing)`：

```bash
dock box my_box.txt --center 20.813 10.963 22.132 --size 23.8 15.3 16.1
```

之后只需要改模板里的具体数值即可直接用于 `dock receptor.pdb ligand.mol2 box.txt`。

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
