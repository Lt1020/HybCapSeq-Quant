# 捕获测序定量计算说明

对应 `hybcapseq_quant` 0.2.1。两套模型分别拟合，输出起始建库体系浓度、提取液中间浓度与污水浓度，不能把这三个层级混用。

## 工作目录和新批次

`PROCESSED_DIR` 是用户每次运行时选择的工作目录，不是代码中的固定常量。使用 `--processed-dir /本批工作目录`，结果默认写入该目录的 `HybCapSeq_Quant_results`；也可用 `--outdir /任意输出目录` 显式覆盖。输入目录单独用 `--esviritu-dir` 指定，不要求固定名称或位于工作目录内。

`init --esviritu-dir /当前比对结果目录 --outdir /当前输入目录` 默认生成通用 `samples.tsv`，只根据本批结果填写样本编号及相对输入路径，不继承历史样本或投入量。按本批记录补填实验参数后即可运行。没有比对文件时可用 `init --sample SAMPLE_A SAMPLE_B --outdir /当前输入目录`，或仅生成空表头。

已有 20261017 的 13 行实测建库记录仅在 `init --preset 20261017` 时复制。以下 1000 倍、2 µL、14 µL 数值为该历史实验示例；实际计算始终读取样本表中 `ercc_dilution`、`ercc_spike_ul`、`reaction_volume_ul` 的值。

输入文件优先级为 `--input-tsv`、本次 `--esviritu-dir`、样本表 `info_tsv`，依次选择。样本表内相对路径以表所在目录为基准，不以 shell 当前目录为基准；搬迁完整批次时保留相对目录结构即可。

## ERCC 投入换算

设官方 Mix 1 第 i 条浓度为 c_i attomol/µL：

```text
原液 copies/µL = c_i × 10^-18 × NA
工作液 copies/µL = 原液浓度 / 1000
实际投入 N_i = 工作液浓度 × 2
起始体系终浓度 = N_i / 14
提取液等效浓度 A_i = N_i / Vlib
ERCC 质量 m_i(ng) = N_i × L_i × 340 / NA × 10^9
NA = 6.02214076×10^23
```

`Vlib` 是进入文库的实际样本提取液体积，不能用反应总体积 14 µL 代替。相同 ERCC 加样条件下，各样本的反应终浓度相同，但提取液等效浓度随 Vlib 改变。

已整理：

- [92 条逐转录本投入表](../templates/quant/ercc92_amounts_1000x_2ul_14ul.tsv)。
- [22 个梯度汇总](../templates/quant/ercc22_levels_1000x_2ul_14ul.tsv)。
- [13 个样本的建库投入核对表](../templates/quant/library_input_audit_20261017.tsv)。

按官方文本数值，终浓度范围约 1.2306778—2,580,917.4686 copies/µL；92 条合计 ERCC 质量约 0.06326964 ng。未检出控制不进入对数拟合，不添加伪计数。数据库多出的四条 ERCC 若检出，会保留原行并注明不属于商业 Mix 1，不虚构浓度。

## 实际混合核酸输入

```text
CNA = RNA浓度 + DNA浓度
Mlib = CNA × Vlib
MTOT = Mlib + Σm_i
```

例如 261008-HN-02-ZKP-a：RNA 36、DNA 24.8 ng/µL，取样 2.77777778 µL，RNA 100 ng、DNA 68.88888889 ng，总样本核酸 168.88888889 ng。不能把 Mlib 全部填成 100 ng。

NC 实际记录为 14 µL 人源白细胞提取物、无 ERCC，不能自行建立标曲，也不会自动借用其他样本的曲线。旧方案个别文字中的 48 条控制和 15 µL 在本次实现中按最新要求与记录更新为 92 条及 14 µL，两套公式定义不变。

## EsViritu 统计字段

| 字段 | 本流程的解释 |
| --- | --- |
| Length | 单条参考长度 |
| Asm_length | 完整定量参考组长度，分节段时用完整组长度 |
| Assembly | 默认定量参考分组；不同 Assembly 不自动合并 |
| read_count | 比对 read 数，不等同于累计比对碱基数 |
| covered_bases | 至少覆盖一次的参考位置数，是覆盖广度的分子 |
| mean_coverage | v1.1.2 全参考位置的平均覆盖深度 |
| filtered_reads_in_sample | 分析输入 read 数，不是总有效碱基数 |
| RPKMF | 原样保留，不进入两套正式计算公式 |
| Pi | 核苷酸多样性，不是 Li 等内标中的 9 bp 标记校正项 |

核对过的 v1.1.2 源码定义 `mean_coverage=sum(coverage)/Length`，因此 `n_x=mean_coverage×Length` 可恢复该工具统计的累计覆盖碱基数。不要把 `covered_bases` 或 `read_count×150` 代入 n_x。该转换不会额外证明比对已满足方案所要求的唯一性、高质量和去重复条件。

nTOT 优先读取 `total_effective_bases`，否则读取同目录 `.fastp.json`。上传脚本已经在外部质控、去宿主，且未用 EsViritu `-q True`，应使用 `fastp_section=before_filtering`；程序核对统计阶段的 read 数。其他流程须填写正确阶段。

只有 info.tsv、没有总碱基数时，方式一输出缺失参数状态，方式二仍可运行。也可统计实际输入 FASTQ：

```bash
hybcapseq_quant measure-input-bases --sample 261008-HN-02-ZKP-a \
  --r1 /path/to/261008-HN-02-ZKP-a_cleaned_R1.fastq.gz \
  --r2 /path/to/261008-HN-02-ZKP-a_cleaned_R2.fastq.gz
```

将输出的真实 `total_effective_bases` 填回输入表，不按标称 PE150 估算。

## 方式一 测序得率法

对应实施方案式（9）—（17）：

```text
x_i = log10(m_i/MTOT)
y_i = log10(n_i/nTOT)
y_i = a + b*x_i
Y_i = (n_i/nTOT)/(m_i/MTOT)
Y = 合格 ERCC 的 Y_i 算术均值
m_v = MTOT*(n_v/nTOT)/Y
N_v = m_v*10^-9*NA/MW_v
Cext = N_v/Vlib
Creaction = N_v/Vreaction
```

同时报告 Y_i 中位数。回归系数用于检验标准曲线；不以 `10^a` 直接替代 Y，不反解回归线替代方案式（13）的 Y 均值算法。

RNA 按 MW=L×340，dsDNA 按 MW=L×660，符合方案式（14）。特殊分子形式须核对实际定量对象和参考长度约定，可在靶标表显式调整质量系数。RNA ERCC 与 DNA 病毒的定量转移仍需 dPCR 验证，不自动假定二者建库、捕获效率相同。

## 方式二 深度标准曲线法

对应方案式（18）—（22）：

```text
A_i = N_i/Vlib
d_i = n_i/L_i
O_i = d_i*CNA/Mlib
log10(A_i) = α + β*log10(O_i)
O_v = d_v*CNA/Mlib
Cext,D = 10^[α + β*log10(O_v)]
Creaction,D = Cext,D*Vlib/Vreaction
```

CNA 与 Mlib 都必须是 RNA＋DNA 输入，且 Mlib=CNA×Vlib。程序检查质量关系。曲线是本 Panel、本次捕获流程的经验曲线；未提供独立 dPCR 训练结果时不生成或套用病毒特异校正系数 Kv。

正式图保留方案坐标。起始体系展示图由同一模型代数换算：

```text
log10(Creaction) = α + β*log10(CNA/Mlib)
                   + log10(Vlib/Vreaction) + β*log10(d)
```

方式一的反应体系展示线按 Y(seq) 推导，斜率固定为 1。两者均不是另拟合一条曲线替代正式模型。

## 污水浓度换算

两种方法取得 Cext 后，均代入式（16）—（17）：

```text
Cconcentrate(copies/mL) = 1000*Cext*Ve1/Vx1
Cww(copies/mL) = Cconcentrate*Vc/(V0*rPM)
Cww(copies/L) = 1000*Cww(copies/mL)
```

Ve1、Vx1 的单位是 µL，Vc、V0 是 mL，rPM 是回收系数。25% 应填 0.25。使用反应体系浓度时，先乘 Vreaction/Vlib 恢复 Cext；不能直接把 14 µL 体系浓度当成提取液浓度。

可直接填 `pmmov_recovery_fraction`；或填 q0、q1、总稀释倍数 D0、D1、原水提取的 Ve0、Vx0，并设 `pmmov_above_loq=1`，由程序执行式（5）—（8）：

```text
C0 = 1000*q0*D0*Ve0/Vx0
C1 = 1000*q1*D1*Ve1/Vx1
rPM = C1*Vc/(C0*V0)
```

回收率不在 0.01—1.50 范围时按方案标为需复核，不输出污水正式浓度。结果是 PMMoV 代理回收校正值，不代表各病毒的真实回收效率均相同。

本次必须确认的字段：

| 输入列 | 填写内容 |
| --- | --- |
| extract_elution_ul | 浓缩液提取物实际洗脱体积 Ve1 |
| extract_input_concentrate_ul | 实际进入提取的原浓缩液体积 Vx1，排除破壁/裂解试剂 |
| wastewater_volume_ml | 实际处理的原始污水体积 V0 |
| pmmov_recovery_fraction | 同源重复的实测回收系数；也可提供完整 dPCR 参数 |

Vc 已按记录填为 0.460 mL；破壁后加样 370 µL 仅作记录；80 µL 是理论洗脱值，实际值需确认。NC 和提取后加入的 ERCC 不换算污水浓度。

## 靶标参考和分节段计数

按 Assembly 分组；缺失时用 Accession。同组累计比对碱基求和，除以完整 Asm_length 得到组平均深度。未检出节段不补造读段。原始 info 表的各节段会显示相同的组基因组浓度，不能再次相加；应使用 `assembly.quantified` 表获取每组一行的结果。

可生成显式靶标输入，检查分子类型和参考长度：

```bash
hybcapseq_quant prepare-targets \
  --esviritu-dir "$ESVIRITU_DIR" \
  --db-metadata /path/to/esviritu_DB/v3.3/virus_pathogen_database.all_metadata.tsv \
  --output "$INPUT_DIR/targets.tsv"
```

核对 quant_group、quant_reference_length、molecule_type 和 mass_per_nt_g_mol。RNA 可填 340，dsDNA 可填 660。没有明确分类和质量系数时，方式一目标浓度为 NA，深度法仍可计算。部分已附官方来源的 dsDNA 病毒科以及明确属于 Orthornavirae 的 RNA 参考可自动识别，其余不猜测。使用自定义表时增加 `--targets /path/to/targets.tsv`。

## 拟合和结果状态

默认至少 4 个相邻浓度梯度进入曲线，R²≥0.95；方式一斜率还需在 0.80—1.20。ERCC 标准表中的 fit_include 控制全局纳入；逐样本排除用 `--ercc-qc ercc_qc.tsv`，列为 sample_ID、ERCC_ID、use_for_fit、reason。必须记录排除原因，不自动挑选最高 R² 窗口或删异常点。

批量模式在同 batch_ID、replicate_group 内评估 ERCC 单条 Y_i CV；至少 3 次有效观测且 CV>35% 时从拟合集合排除。不足 3 次只报告现有 CV，不声称通过三重复验证。另输出文库 Y(seq) 的重复 CV 汇总。

默认不外推，不用最低检出的 ERCC 推断 LoQ。低于外部提供 LoQ 或超出校准范围时，新增正式浓度列为 NA，未筛选估计只保留在审计列。`--allow-extrapolation` 可显式允许外推，仍标注超出范围。

`quant_status` 区分已知 ERCC 投入、标曲失败、分子类型缺失、超范围、低于 LoQ、比对/检出验证不足和缺少污水参数等状态。默认报告可计算的估计；`--strict-qc` 要求提供已验证的 alignment_qc_verified、detection_qc_verified 和 LoQ，条件不足时不输出正式目标浓度。

## 可选 BAM 统计

```bash
python -m pip install '.[bam]'
hybcapseq_quant collect-bam --bam /path/to/sample.sorted.bam \
  --min-mapq 20 --output /path/to/sample.alignment_metrics.tsv
```

要求索引 BAM；排除 secondary、supplementary、QC-fail、已经标记为 duplicate 的比对，默认 MAPQ≥20，NH 标签存在时要求 NH=1。按匹配覆盖块计数，双端各 read 独立统计，不合并配对重叠。它不识别未标记的 PCR 重复，也不单凭 MAPQ 宣称完成唯一性验证。将统计路径填入 alignment_metrics_tsv 后，两套计算使用同一套重算信号。

## 实验验证范围

本次 13 行样本表来自实际记录，SIMULATED_01 的信号为软件模拟。真实回归方程、R²、斜率和浓度应在实际 TSV 上运行得到。公开 ERCC 浓度是商业标称配方；如有批次实测赋值，应替换为可追溯的浓度表。程序保留输入和来源 SHA256，便于复算。
