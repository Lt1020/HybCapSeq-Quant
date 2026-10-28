# HybCapSeq-Quant

[![Quantification verification](https://github.com/Lt1020/HybCapSeq-Quant/actions/workflows/tests.yml/badge.svg)](https://github.com/Lt1020/HybCapSeq-Quant/actions/workflows/tests.yml)

HybCapSeq-Quant 是用于杂交捕获测序的 ERCC 校准定量软件。程序读取 EsViritu 的病毒与 ERCC 联合比对结果，分别实施 **sequencing yield（测序产率法）**和 **depth calibration（深度校准法）**，计算提取液及建库起始体系中的核酸浓度，并依据体积关系和 PMMoV 回收率换算污水浓度。

软件代码、命令参数、模板字段、输出文件及图件标注均使用英文。EsViritu 原始列名按上游格式保留。输入路径完全由调用参数或样本表提供；样本模板不预填实验体积、ERCC 加样量或回收率。

## 安装

需要 Linux、Python 3.10–3.12。推荐使用独立环境：

```bash
git clone https://github.com/Lt1020/HybCapSeq-Quant.git
cd HybCapSeq-Quant
conda env create -f environment.yml
conda activate hybcapseq_quant
python -m pip install .
hybcapseq_quant --version
```

如需从 BAM 提取经过显式过滤的比对统计：

```bash
python -m pip install '.[bam]'
```

论文分析应固定 Git 提交、软件依赖及数据库校验和。运行记录包含版本、参数、输入文件 SHA256 和样本级计算依据。

## 升级与数据库注释安装

已有克隆和安装的环境可以直接更新：

```bash
cd HybCapSeq-Quant
git pull --ff-only
python -m pip install --upgrade .
hybcapseq_quant --version
```

本版本为 **2.2.0**，配套数据库修订为 **v3.3.1**。也可以在当前 Python/conda 环境直接安装固定发行版本：

```bash
python -m pip install --upgrade --no-cache-dir "git+https://github.com/Lt1020/HybCapSeq-Quant.git@v2.2.0"
hybcapseq_quant --version
```

将环境变量 `DB_DIR` 设置为实际联合参考数据库目录。已有 v3.3 用户执行一次数据库升级，然后校验：

```bash
hybcapseq_quant upgrade-db --db-dir "$DB_DIR"
hybcapseq_quant validate-db --db-dir "$DB_DIR"
```

升级在原目录完成，不要求修改目录名。它校验原 metadata、FASTA 和注释指纹，备份被替换文件，再校正 Rochester 乙流 8 条参考的组装总长度并安装匹配注释。重复执行不会重复修改；可加 `--check` 先只检查。原有 FASTA、MMI 及 `ESVIRITU_DB` 环境变量保留，三个数据库文件名保持不变，无需重新比对旧数据。升级后的注释文件仍为：

- `virus_pathogen_database.molecular_annotations.tsv`
- `virus_pathogen_database.molecular_annotations.manifest.json`

目录另保存 `virus_pathogen_database.reference_corrections.json`、v3.3.1校验清单和升级报告。仅当 metadata 已经是 v3.3.1、需要补装注释时使用 `install-references --db-dir "$DB_DIR"`。未知或手动改动的数据库不会被自动覆盖。完整说明见 [升级说明](docs/UPGRADE_2.1.0.md)。

## Accession 分子注释

注释覆盖 v3.3 的 **31,003** 条参考，包括 30,907 条病毒参考和 96 条 ERCC 转录本参考。逐条检索 NCBI Nuccore 的精确 `Accession.version` 并核对长度，再匹配 ICTV 官方 `VMR_MSL41.v1.20260729` 的基因组类型。ICTV 匹配层级依次为 accession、species、genus、family、order；分类层级支持的分子类型必须一致。`f__` 等前缀仅在匹配时移除，原始注释保留。

NCBI 中的序列记录可能描述 DNA 克隆、cDNA 或复制形式，因此同时保留 NCBI 记录类型和 ICTV 基因组类型。病毒基因组当量模型优先采用 ICTV 基因组类型，记录差异及证据，不据此断言个别样本中的实际物理核酸形式。不能确认的链型或质量保留 `NA`。

| 基因组类型 | 参考位置质量系数，g/mol | 条目数 |
|---|---:|---:|
| ssRNA | 340 | 17,868 |
| dsRNA | 680 | 7,579 |
| ssDNA | 330 | 3,740 |
| dsDNA | 660 | 1,542 |
| 部分双链 DNA | NA | 263 |
| RNA，链型未确认 | NA | 11 |

以上系数是明确声明的平均质量近似模型。dsRNA 的 680 为两条 RNA 链各 340 的合计；不是另一份实验测得的常数。参考序列长度表示一条链的位置数；双链系数已包括互补链。注释表分别输出理论摩尔质量（g/mol）、单分子理论质量（g 和 ng），并区分单条参考与完整 assembly 质量。30,729 条参考可赋予理论质量；274 条保留 NA。v3.3.1 已将 Rochester 乙流参考组的声明总长度与8条现存序列的长度和统一为13,859 nt。96 条 ERCC 均可描述为 RNA 转录本参考，商业 Mix 1 浓度只赋予其中 92 条。

下载文件、全部字段、来源与离线重建方法见 [分子注释说明](references/molecular_annotations/README_zh.md)。软件运行无需联网检索 NCBI/ICTV。原始来源记录摘要和校验和随软件分发。

## 准备输入

从实际 EsViritu 结果目录发现样本并生成模板：

```bash
hybcapseq_quant init \
  --esviritu-dir esviritu_results \
  --outdir inputs
```

也可以先按样本名建立空白输入表：

```bash
hybcapseq_quant init --sample SAMPLE_A SAMPLE_B --outdir inputs
```

初始化生成 `samples.tsv`、`ercc_standards.tsv`、`input_fields.tsv`、`ercc_qc.tsv` 和 `targets_template.tsv`。在 `samples.tsv` 中填写实际实验参数。完整字段定义见 [输入字段表](templates/input_fields.tsv) 和 [数据字典](docs/DATA_DICTIONARY.md)。未知值填写 `NA`，布尔值填写 `1` 或 `0`。

| 参数组 | 必须明确的信息 |
|---|---|
| 样本与建库 | 样本标识、样本类型、提取液投入体积、DNA+RNA 总质量与总浓度、建库起始总体积 |
| ERCC | Mix 1 原液总稀释倍数、实际加样体积 |
| 测序信号 | 原始检出表；产率法另需实际分析输入的总碱基数或匹配的 fastp 报告 |
| 污水换算 | 原水体积、总浓缩液体积、实际提取 aliquot 和洗脱体积、PMMoV 回收率或完整配对 dPCR 数据 |
| 定量报告 | 比对与检出质控、ERCC 至靶标的校准可转移性验证、独立确定的 LoQ |

`library_input_mass_ng` 必须包含共同进入文库的 DNA 和 RNA，并与 `extract_nucleic_acid_concentration_ng_ul × library_input_volume_ul` 一致。`reaction_volume_ul` 是反应总体积，不能替代提取液投入体积。使用 fastp 报告时，必须明确填写与实际分析输入一致的 `fastp_section`。

针对不同靶标的分子类型、基因组长度、LoQ 或验证状态，可生成并审阅靶标表：

```bash
hybcapseq_quant prepare-targets \
  --esviritu-dir esviritu_results \
  --db-dir "$DB_DIR" \
  --output inputs/targets.tsv
```

`--db-dir` 从全量注释载入精确 accession 的分子类型和完整 assembly 长度；仅需逐靶标提供实验验证信息时再使用 `--targets`。不使用数据库注释时，产率法必须提供有来源的显式靶标质量赋值表。靶标表中非 `NA` 的验证值和 LoQ 优先于样本表。阴性对照应按真实操作填写 `negative_control` 和 ERCC 投入，不根据样本名称推断加样情况。

## 运行定量

批量处理样本表中的全部样本：

```bash
hybcapseq_quant run \
  --samples inputs/samples.tsv \
  --esviritu-dir esviritu_results \
  --db-dir "$DB_DIR" \
  --outdir quantification_results \
  --method both
```

单样本模式：

```bash
hybcapseq_quant run \
  --samples inputs/samples.tsv \
  --sample SAMPLE_A \
  --input-tsv esviritu_results/SAMPLE_A.detected_virus.info.tsv \
  --db-dir "$DB_DIR" \
  --outdir quantification_sample_a \
  --method both
```

`--sample` 可指定一个或多个样本；省略时处理整张样本表。`--method` 可选 `both`、`sequencing_yield` 或 `depth_calibration`。污水单位默认 `copies/mL`，通过 `--wastewater-unit copies/L` 切换。

需要靶标级验证或实际核酸材料的明确质量赋值时，可附加 `--targets inputs/targets.tsv`。该表中非 NA 的值优先；分子类型与质量系数必须一致，质量赋值须提供来源。

输入位置优先级为 `--input-tsv`、`--esviritu-dir`、样本表 `info_tsv`。程序递归发现原始检出表；同一样本匹配到多个文件时终止。样本表中的相对路径以该表所在目录解析，命令中的相对路径以执行目录解析。

`--outdir` 可指定任意输出位置。省略该参数时，`--processed-dir` 指定的工作目录下生成 `hybcapseq_quant_results`；两个参数均省略时，结果写入当前目录下的 `hybcapseq_quant_results`。非空输出目录不会被覆盖。

## 定量模型与报告规则

两种方法分别拟合自己的校准模型，不使用 RPKMF 替代方案中的比对碱基数或标准化深度。产率法使用合格 ERCC 的产率**算术均值**推算目标质量；深度法按标准化深度回归推算提取液等效浓度。完整公式、变量单位、回归方向、分节段处理和假设见 [METHODS](docs/METHODS.md)。

默认 `--reporting-mode validated`：模型估计值始终写入 `estimated_*` 字段；只有校准、校准范围、LoQ 及外部验证条件满足时，才填写报告浓度列。`alignment_qc_verified`、`detection_qc_verified` 和 `calibration_transfer_verified` 是使用者提供的验证声明，软件不会通过读取 TSV 自动建立这些证据。探索性分析可显式使用 `--reporting-mode exploratory`，结果仍保留验证状态与限制标记。

默认要求至少 4 个相邻浓度等级、R² ≥ 0.80，两种方法均要求正斜率。产率法斜率 0.8–1.2 为诊断范围，不作为默认拦截条件；需要重现旧规则时显式设置 `--min-r2 0.95 --enforce-yield-slope-range`。`exploratory` 模式默认输出可计算的区间外模型估计，并保留 `outside_calibration_range`；`validated` 模式默认不允许外推。可用 `--no-extrapolation` 或 `--allow-extrapolation` 显式选择。最低ERCC浓度和校准下界均不自动构成LoD/LoQ；LoQ只使用明确提供的靶标或样本参数。最低读数、覆盖广度和各阈值应在分析前确定，并在方法部分报告。重复间 CV 默认仅作诊断；需要按预先设定规则排除 ERCC 时，使用 `--filter-replicate-cv --max-replicate-cv 0.35`，且每条 ERCC 至少有 3 个可比较重复。

缺少回收率或体积参数时，污水浓度为 `NA`；不会默认回收率为 1。ERCC 在提取后加入，其反应体系浓度为已知投入值，污水浓度始终为 `NA`。无 ERCC 的样本不自动借用其他样本的校准曲线。

## 输出

两套结果分别存放在 `sequencing_yield/<sample_id>/` 和 `depth_calibration/<sample_id>/`。文件名均使用样本标识与固定英文后缀：

| 文件后缀 | 内容 |
|---|---|
| `.reference_quantification.tsv` | 原始检出表全部行列；在 RPKMF 后插入反应体系浓度和污水浓度，再附加估计值与质控字段 |
| `.target_quantification.tsv` | 每个定量参考组一行，汇总分节段信号；用于目标层面的统计 |
| `.calibration_points.tsv` | 全部 92 条 ERCC 的投入量、检出、纳入/排除、模型坐标和残差 |
| `.calibration_model.tsv` | 回归方程、斜率、截距、R²、残差指标及接受状态 |
| `.calibration_curve.png/.pdf` | 主模型散点图与回归线；PNG 为 300 dpi，PDF 为矢量图 |
| `.reaction_concentration_curve.png/.pdf` | 同一模型在 ERCC 反应体系浓度尺度的展示 |
| `.reaction_concentration_model.tsv` | 上述浓度尺度展示所用的代数转换方程 |
| `.wastewater_conversion.tsv` | 体积换算、PMMoV 回收率及换算状态 |
| `.analysis_manifest.json` | 样本参数、模型、阈值与来源校验和 |

主结果列为 `concentration_reaction_copies_ul` 和 `concentration_wastewater`。后者的单位由 `wastewater_concentration_unit` 给出。`estimated_*` 与报告列必须结合 `value_origin`、`reporting_mode` 和 `quantification_status` 解读。分节段参考的每个检出行可重复显示同一目标浓度；不能将这些行相加。

结果根目录另包含输入快照、`run_summary.tsv`、`replicate_diagnostics.tsv`、`reference_corrections.tsv` 和 `run_manifest.json`。`inputs/molecular_annotations.tsv` 保存本批次检出参考的注释快照；manifest 记录全量注释、来源 manifest 和数据库 metadata 的 SHA256。原检出表的 `Asm_length` 保留不变，实际计算使用 `reference_length_nt`；只有已核实的历史长度值能够按公开修正规则兼容。主结果附带 `reference_correction_ids`、`molecule_type`、`molecular_annotation_status`、`mass_model` 和 `quantification_molecular_mass_g_mol`。运行汇总区分可估计、可报告、超出校准范围和缺少质量的靶标数；质量完整不代表所有结果都满足报告条件。

## 软件验证与方法适用性

```bash
python -m unittest discover -s tests -v
python examples/quant/generate_demo.py
hybcapseq_quant run \
  --samples examples/quant/samples.tsv \
  --esviritu-dir examples/quant \
  --targets examples/quant/targets.tsv \
  --outdir analytical_example_results
```

示例为具有解析真值的模拟数据，所有文件均明确标注为模拟数据。验证范围见 [VALIDATION](VALIDATION.md)。软件测试验证公式和数据处理；捕获偏倚、ERCC 与目标的响应可转移性、PMMoV 的代表性、检出限、定量限和实验测量不确定度须通过独立实验评价。R² 合格本身不证明绝对浓度无偏，也不构成 LoQ。

## 联合参考数据库与来源

[自定义 EsViritu 数据库 v3.3](https://github.com/Lt1020/HybCapSeq-Quant/releases/tag/db-v3.3) 包含官方 v3.2.4 病毒参考与 96 条 NIST ERCC 参考。96 条序列用于联合比对，其中商业 Mix 1 中的 92 条用于浓度赋值。数据库说明及校验和见 [数据库说明](databases/v3.3/README_v3.3.md)。该数据库是本项目扩展版本。

- [EsViritu](https://github.com/cmmr/EsViritu)：联合比对结果及上游字段定义。
- [Thermo Fisher ERCC Mix 1/Mix 2 官方浓度表](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt)：92 条转录本的公开标称浓度；原始文件与 SHA256 随包保存。
- [NIST SRM 2374](https://www-s.nist.gov/srmors/view_detail.cfm?srm=2374)：序列来源及本仓库中的 [逐条参考记录](references/ercc96/README_zh.md)。

软件采用 MIT 许可证。第三方参考材料保留其来源及适用条款。

## 2.2.0 分靶标偏差评价与数据库状态

新增 `fit-correction` 和 `apply-correction`，详见 [英文说明](README.md)及[偏差校正方法](docs/TARGET_CORRECTION.md)。30%为独立留出评价目标，不保证每个靶标都能通过。未经足够留出评价的参数仅供探索性分析；原始浓度不覆盖。

数据库 v3.4 尚缺少实际数字PCR靶区注释和已审核参考FASTA，未生成序列资产或索引，未发布正式版本。当前有效数据库仍为 v3.3.1。详见 [v3.4准备状态](databases/v3.4/README.md)。

