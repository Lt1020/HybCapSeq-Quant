# HybCapSeq Quant

`hybcapseq_quant` 读取 EsViritu 的 `<sample>.detected_virus.info.tsv`，分别按实验方案的 **Y(seq) 比对碱基数法**和**平均深度标准曲线法**计算建库起始体系浓度，再以体积关系和 PMMoV 回收率换算污水浓度。原有 `ercc-audit` 命令继续保留。

工作目录、EsViritu 结果目录、样本输入文件和输出目录均由运行参数指定，程序不依赖固定磁盘、用户名、日期或样本编号。92 条转录本的配方值来自 Thermo Fisher 官方表，实际稀释倍数、加样量、反应体积和样本投入从各样本输入行读取。

附带的 20261017 历史示例为：商业 ERCC Mix 1 的 92 条转录本，原液稀释 1000 倍、加入 2 µL，建库起始总体积 14 µL；按 RNA 100 ng 取样，同时计入共同进入文库的 DNA 质量。该示例可显式选择，不会默认套用于新批次。

数据库 v3.3 含 96 条 NIST 参考，但商业 Mix 1 不含 ERCC-00007、ERCC-00018、ERCC-00023、ERCC-00128，曲线只使用商业 92 条。

## 安装

建议用独立环境，不需要修改 EsViritu 环境或默认的 `ESVIRITU_DB`。

```bash
conda create -n hybcapseq_quant python=3.12 pip -y
conda activate hybcapseq_quant
git clone https://github.com/Lt1020/HybCapSeq-Quant.git
cd HybCapSeq-Quant
python -m pip install .
hybcapseq_quant --version
```

Python 最低版本 3.10，绘图依赖 Matplotlib。若之前已克隆，先 `git pull --ff-only` 再重新安装。

## 准备任意批次输入

```bash
PROCESSED_DIR="/path/to/current_batch"
ESVIRITU_DIR="/path/to/current_esviritu_results"
INPUT_DIR="$PROCESSED_DIR/HybCapSeq_Quant_inputs"

hybcapseq_quant init \
  --esviritu-dir "$ESVIRITU_DIR" \
  --outdir "$INPUT_DIR"
```

上述路径均为占位示例，替换为本批实际位置。EsViritu 目录可以在工作目录内，也可以在另一块磁盘，子目录名称不限。

生成通用 `samples.tsv`、`ercc92_mix1.tsv`、`input_fields.tsv` 和可选的 `ercc_qc.tsv`。程序从当前结果目录递归提取样本编号，为每个实际样本生成一行；实验投入参数留为 `NA`，需按本批记录填写。没有结果文件时可用 `init --sample SAMPLE_A SAMPLE_B --outdir "$INPUT_DIR"` 指定新编号，或仅用 `init --outdir "$INPUT_DIR"` 生成仅含表头的模板。

初次准备新批次时，需要填写样本类型、核酸浓度/质量、建库提取液体积、ERCC 稀释/加入量和反应总体积。NC 应标记为 `negative_control` 并按实际填写 ERCC 投入，不能根据编号自动假定实验操作。

**需补齐**实际 `extract_elution_ul`、`extract_input_concentrate_ul`、`wastewater_volume_ml`，以及 `pmmov_recovery_fraction` 或完整 PMMoV dPCR 参数。缺失时仍输出可计算的建库体系浓度，污水浓度保留为 `NA`，不会默认回收率为 1。

若确实处理之前的 20261017 数据，可显式使用 `hybcapseq_quant init --preset 20261017 --outdir /path/to/inputs_for_20261017`。该模式才生成已填实际 13 行的 `samples_20261017.tsv`。其中浓缩液 460 µL 换算为 0.460 mL，破壁后 370 µL 与理论洗脱 80 µL 仅保留作记录；不替代实际提取/洗脱体积，原水体积也须确认。

## 批量运行

```bash
hybcapseq_quant \
  --samples "$INPUT_DIR/samples.tsv" \
  --ercc-standards "$INPUT_DIR/ercc92_mix1.tsv" \
  --esviritu-dir "$ESVIRITU_DIR" \
  --processed-dir "$PROCESSED_DIR" \
  --method both \
  --wastewater-unit copies/mL
```

结果自动保存在 `$PROCESSED_DIR/HybCapSeq_Quant_results`。`--outdir` 可指定任何其他位置，并优先于此默认位置；未提供这两个参数时才使用当前运行目录下的 `HybCapSeq_Quant_results`。

程序递归查找各样本的原始检测 TSV，支持共用目录或样本子目录。输入位置优先级为：`--input-tsv` > 本次显式传入的 `--esviritu-dir` > 样本表中的 `info_tsv`。这样切换到新结果目录时，旧表内路径不会覆盖本次命令的选择。同一样本发现多个同名结果时会报错，避免混用批次。

样本表内相对文件路径均相对于该样本表所在目录解析，与执行命令时的当前目录无关。初始化自动填写相对 `info_tsv`，因此整批数据按原相对结构搬迁后仍可使用。脚本入口 [`scripts/run_quant.sh`](scripts/run_quant.sh) 只转发命令参数；旧的带日期脚本名也保留兼容入口，两者都不含固定实验路径。

## 单样本运行

```bash
hybcapseq_quant \
  --samples "$INPUT_DIR/samples.tsv" \
  --sample SAMPLE_A \
  --esviritu-dir "$ESVIRITU_DIR" \
  --outdir "$PROCESSED_DIR/HybCapSeq_Quant_results_single" \
  --method both \
  --wastewater-unit copies/L
```

`--sample` 可列出多个样本，省略则分析样本表中的全部样本。单样本也可用 `--input-tsv /absolute/path/sample.detected_virus.info.tsv` 直接指定文件。

## 输出

两套结果分别在 `method1_yield/<sample>/` 和 `method2_depth/<sample>/` 中保存。每个文件均带样本名前缀：

- `detected_virus.info.quantified.<method>.tsv`：保留原始所有行列，在 `RPKMF` 后新增 `浓度concentration_ERCC` 和 `污水浓度concentration_wastewater`。
- `assembly.quantified.<method>.tsv`：每个定量参考组一行，避免分节段结果重复求和。
- `ERCC_calibration_points.tsv`、`ERCC_regression.tsv`：投入、检出状态、拟合集合、排除原因、斜率、截距、R²、方程及 Y(seq)。
- `ERCC_standard_curve.png/pdf`：原方案坐标的正式标准曲线；`ERCC_reaction_concentration_curve.png/pdf`：同一模型在起始体系浓度尺度的展示图。
- `wastewater_conversion.tsv`、`quantification_audit.json`：换算参数、计算依据、输入校验值和质量状态。

根目录另保存 `run_summary.tsv`、逐样本 ERCC92 投入表和技术重复 CV 汇总。

ERCC 行的 `浓度concentration_ERCC` 为已知加样终浓度；病毒行为模型估计的起始体系浓度，均为 copies/µL。ERCC 在提取后加入，其污水浓度为 `NA`。需要英文列名时加 `--ascii-columns`，得到 `concentration_ERCC` 和 `concentration_wastewater`。

## 计算与质量说明

- 方式一严格使用方案式（9）—（17）：比对碱基分数与投入质量分数求单条 Y(seq)，取合格控制的算术均值回算目标质量与拷贝数。回归不替代 Y(seq) 均值算法。
- 方式二严格使用式（18）—（22）：拟合已知提取液等效浓度与按混合核酸输入标准化的平均深度，再执行式（16）—（17）。
- 核对过的 EsViritu v1.1.2 中，全参考 `mean_coverage × Length` 可恢复其累计覆盖碱基数。`covered_bases` 是覆盖广度分子，RPKMF 和 `read_count×150` 均不能替代该值。
- 方式一的总有效碱基数默认从同目录 `.fastp.json` 读取；上传脚本未使用 EsViritu `-q True`，默认选 `before_filtering`，并核对 read 数与 TSV 一致。其他流程需指定正确统计阶段或填写真实 `total_effective_bases`。
- RNA 和 dsDNA 的分子质量系数按方案为 340、660；部分类型按明确分类识别，不明类型需在 `targets.tsv` 中填写，程序不会猜测。深度法不依赖分子质量系数。
- 仅有 TSV 无法证明高质量唯一比对、PCR 去重复、覆盖均一性、阴性背景或 LoQ 已通过验证。默认保留可计算的估计并注明状态；`--strict-qc` 会在验证条件不足时将正式目标浓度留为 `NA`。
- 默认不外推，不自动删异常点提高 R²，不把最稀的已检出 ERCC 当作已验证 LoQ。

**本次真实测序 TSV 与 PMMoV dPCR 结果尚未提供；仓库中的回归演示为模拟数据，不代表实际实验的曲线性能或浓度。**

详细公式、待填参数和可选 BAM 统计见 [定量计算说明](docs/QUANTIFICATION_zh.md)。

## 软件验证演示

```bash
python examples/quant/generate_demo.py
hybcapseq_quant --samples examples/quant/samples.tsv \
  --esviritu-dir examples/quant --outdir quant-demo
python -m unittest discover -s tests -v
```

模拟 RNA 和分节段 RNA 各投入 10,000 copies。两种方法应回算为 714.285714 copies/µL 起始体系、5000 copies/µL 提取液，以及 61333.333333 copies/mL 污水。这些模拟参数仅验证软件算术和文件流程。

## 数据库与来源

[数据库 v3.3 下载及校验附件](https://github.com/Lt1020/HybCapSeq-Quant/releases/tag/db-v3.3)：官方 DB v3.2.4 加入 NIST ERCC96 的自定义扩展，共 31,003 条参考；不是上游 EsViritu 官方 v3.3。说明及构建配置见 [`databases/v3.3`](databases/v3.3)。

- [Thermo Fisher 官方 ERCC92 浓度表](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt)：随包保存原始文本和 SHA256，公开标称配方不替代具体批次实测赋值。
- [NIST SRM2374 certified v1](https://tsapps.nist.gov/srmext/certificates/documents/SRM2374_Sequence_v1.FASTA)：与 v3.3 相同的参考长度。
- [EsViritu v1.1.2 统计源码](https://github.com/cmmr/EsViritu/blob/v1.1.2/src/EsViritu/esv_funcs.py)。
- 用户实施方案式（5）—（22）及实际记录，文件名与 SHA256 见 `src/hybcapseq_quant/data/provenance.json`。
- Li B et al. AEM 2021;87:e00871-21，DOI: [10.1128/AEM.00871-21](https://doi.org/10.1128/AEM.00871-21)。
- Langenfeld K et al. NAR 2025;53:gkaf118，DOI: [10.1093/nar/gkaf118](https://doi.org/10.1093/nar/gkaf118)。

代码采用 MIT 许可证；上游数据库与官方参考资料的来源及原许可分别保留。
