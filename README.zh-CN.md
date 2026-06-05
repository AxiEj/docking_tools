# AutoDock-GPU + gnina CNN 重打分工作流

[English version](README.en.md) | [安装与运行指南](INSTALL_AND_RUN.md) | [详细流程说明](README_ADGPU_GNINA_PIPELINE.md)

本仓库是一个本地一键分子对接流程：先用 **AutoDock-GPU / AutoGrid4** 进行 AutoDock4 力场采样，再用 **gnina CNNscore / CNNaffinity** 对 AutoDock-GPU 产生的 pose 做深度学习重打分和排序。

## 标签 / 关键词

`AutoDock-GPU` · `gnina` · `gnina CNN rescoring` · `CNNscore` · `AutoGrid4` · `Meeko` · `Open Babel` · `molecular docking` · `protein-ligand docking` · `virtual screening` · `GPU docking` · `CUDA` · `drug discovery` · `structure-based drug design`

## 为什么这个流程有用

- **AutoDock-GPU**：GPU 加速 AutoDock4 搜索/局部优化，适合把 AutoDock4 风格的 grid/GPF/DLG 工作流跑得更快。
- **gnina CNN 重打分**：gnina 1.0 论文报告，CNN scoring function 对输出 pose 重打分时，在 redocking / cross-docking 的 Top1 pose-ranking 指标上优于 AutoDock Vina scoring；定义 binding pocket 时，Top1 从 58%→73%（redocking）和 27%→37%（cross-docking）。
- **本仓库的组合思路**：保留 AutoDock-GPU 的快速采样和 DLG 输出，再把每个 run-best pose 清理后交给 gnina `--score_only --cnn_scoring rescore` 排名，得到更适合后续人工检查/筛选的 CNNscore 排序表。

> 注意：这里不是宣称“所有体系无条件最高精度”。文献背书的是 gnina CNN 重打分在其基准任务中相对 Vina scoring 的 pose-ranking 改进；真实项目仍建议配合阳性对照、重复 docking、构象/质子化检查和后续 MD/MM-GBSA/实验验证。

## 文献背书

- McNutt A. T. et al. **GNINA 1.0: molecular docking with deep learning.** *Journal of Cheminformatics* 13, 43 (2021). DOI: [10.1186/s13321-021-00522-2](https://doi.org/10.1186/s13321-021-00522-2)
- Santos-Martins D. et al. **Accelerating AutoDock4 with GPUs and Gradient-Based Local Search.** *Journal of Chemical Theory and Computation* 17, 1060-1073 (2021). DOI: [10.1021/acs.jctc.0c01006](https://doi.org/10.1021/acs.jctc.0c01006)

## 快速使用

安装依赖和重建被 Git 忽略的 `downloads/`、`bin/`，见：[INSTALL_AND_RUN.md](INSTALL_AND_RUN.md)。

配置好 `dock` alias 后：

```bash
dock box
dock receptor.pdb ligand.mol2 box.txt
```

## box 模板

`dock box` 默认只启用 AutoDock Grid 三行；Vina / LeDock / BoxCode 作为注释参考保留，解析器会忽略 `#` 注释：

```text
# !!! MUST CHANGE FOR YOUR REAL SYSTEM !!!
# !!! 下面 AutoDock Grid 三行只是模板默认值；真实对接前必须换成你的口袋参数 !!!

*********AutoDock Grid Option*********
npts 49 75 43 # num. grid points in xyz
spacing 0.375 # spacing (A)
gridcenter 54.809 59.493 30.181 # xyz-coordinates or auto
```

真实对接前，只需要把 `npts`、`spacing`、`gridcenter` 换成目标口袋参数。

## 仓库内容边界

本仓库只跟踪 Python/Bash 工作流和文档。本地测试输入也会被 Git 忽略。以下内容故意不提交到 GitHub：

- 本机安装的二进制文件：`bin/`
- 下载缓存：`downloads/`
- 本地测试输入：`test/`
- 运行日志和 docking 输出
- OMX / 本地代理状态

这样可以避免把 2GB 级别的 gnina CUDA 二进制和机器相关产物推到 GitHub。完整安装步骤见：[INSTALL_AND_RUN.md](INSTALL_AND_RUN.md)。
