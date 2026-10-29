# v3.4 数据库与自动化定量流程

本流程说明数据库构建、重新比对、ERCC 校准及数字 PCR 对照的计算步骤。数据库 v3.4 与软件 2.2.0 分别记录版本；新增公开文献靶区不等同于确认了实际商业试剂盒的扩增区域。

## 1. 数据库组成与测量对象

v3.4 由原 v3.3 的 31,003 条参考、v3.3.1 的 8 处乙流组装长度修订，以及 21 条经过核验的 PCRROI 参考片段组成，共 31,024 条。原有病毒与 ERCC 序列保持不变。新冠新增记录固定为 N1、N2、E、ORF1ab 共 4 条。主库新增片段覆盖 11 个靶标；完整审核资料覆盖 12 个靶标、35 套文献体系和 34 条去重片段。其余 13 条边界不足或存在失配的候选仍随包保存，未混入主比对 FASTA。

PCRROI 是“PCR 目标参考区域”的自定义前缀。片段包含两端引物结合位点、中间扩增区及外侧各 100 bp。起止位置是固定 accession.version 上正链、1 起始且包含两端的坐标。`Length` 与 `Asm_length` 均为片段实际长度；不填写原病毒全基因组长度。逐条来源、参考角色、坐标、引物匹配和失配记录见 `references/pcr_roi/`。

原病毒参考与新增片段存在重叠。联合比对中的竞争归属会改变局部参考和原参考之间的信号分配。不同检测体系的片段也可能彼此重叠。因此，不能将它们的 reads、RPKMF 或浓度直接相加，也不能将 PCRROI 名称视为新增病毒物种。物种检出统计应排除 ERCC 和 PCRROI 重复表示，再按实际分类与节段归并规则汇总。

## 2. 安装、下载与核验

```bash
git clone https://github.com/Lt1020/HybCapSeq-Quant.git
cd HybCapSeq-Quant
conda env create -f environment.yml
conda activate hybcapseq_quant
python -m pip install .
hybcapseq_quant --version

curl --fail --location --output esviritu_db_v3.4.tar.gz \
  https://github.com/Lt1020/HybCapSeq-Quant/releases/download/db-v3.4/esviritu_db_v3.4.tar.gz
curl --fail --location --output esviritu_db_v3.4.sha256 \
  https://github.com/Lt1020/HybCapSeq-Quant/releases/download/db-v3.4/esviritu_db_v3.4.sha256
sha256sum -c esviritu_db_v3.4.sha256
tar -xzf esviritu_db_v3.4.tar.gz
cd v3.4
sha256sum -c SHA256SUMS
cd ..
DB_DIR="$PWD/v3.4"
hybcapseq_quant validate-db --db-dir "$DB_DIR"
```

将实际 EsViritu 命令的 `--db` 设置为该目录，保留完整命令、软件版本、输入 FASTQ 指纹和运行日志。v3.4 已附匹配的分子注释，不再执行针对旧库的 `upgrade-db` 或 `install-references`。新增序列后必须重新比对；旧检出表不能通过改名转成 v3.4 结果。

数据库目录保留原三个文件名，另包含匹配的分子注释、来源及坐标表、索引查询审核表、构建清单、校验和与完整文献片段资料。索引采用 minimap2 2.31-r1302，`sr` 预设，k=21、w=11。构建时逐条核对 FASTA、16 列 metadata 和索引的编号及长度，并核查 96 条 ERCC 和 21 条局部参考的查询命中。该核验检验文件与索引的一致性，不代表真实样本检出性能验证。

## 3. 输入表与信号准备

```bash
hybcapseq_quant init --esviritu-dir esviritu_results --outdir inputs
hybcapseq_quant prepare-targets \
  --esviritu-dir esviritu_results --db-dir "$DB_DIR" \
  --output inputs/targets.tsv
```

检查 `samples.tsv` 中实际投入体积、DNA+RNA 总质量、核酸总浓度、建库起始总体积、ERCC 原液总稀释倍数及加样量。污水换算还需要原水体积、总浓缩液体积、实际提取体积、洗脱体积和独立回收率。未知信息保持 `NA`，不以默认实验体积或回收率替代。总投入质量应与核酸浓度乘以投入体积一致。

数据库含 96 条 ERCC，商业 Mix 1 的标称浓度只赋予其中 92 条。未检出或未参与拟合的 ERCC 仍在审核表中保留。建库后加入的 ERCC 不能校正前期浓缩和提取损失。

在确认上游 `mean_coverage` 是含零覆盖位置的全参考平均深度时，比对碱基信号为 `mean_coverage × Length`。`covered_bases` 用于覆盖广度，不能代替累计比对碱基数。得率法需要与实际分析读段集合一致的总输入碱基数。上游软件或过滤规则改变时，应重新核对信号定义。

如需限制到某一物理 PCR 区域，应预先明确基于 BAM 的区间、读段或片段归属以及多重比对处理规则。当前数据库增加参考序列，但没有自动实现“只计算扩增区域且跨参考去重”的新计数算法；不能将主比对中的 PCRROI reads 直接视为该区域全部 reads。

## 4. ERCC 校准与浓度计算

```bash
hybcapseq_quant run \
  --samples inputs/samples.tsv --esviritu-dir esviritu_results \
  --db-dir "$DB_DIR" --targets inputs/targets.tsv \
  --outdir quantification_results --method both \
  --min-r2 0.80 --reporting-mode exploratory \
  --no-extrapolation --wastewater-unit copies/L
```

示例显式限制在校准范围内并输出探索性结果。正式报告采用 `validated` 模式，并要求独立的比对、检出及校准转移验证和靶标 LoQ。两种模式均不改变模型本身的质量要求。

两种方法均要求至少 4 个相邻标称浓度等级、正斜率和 R²≥0.80。得率法斜率 0.8–1.2 仅作诊断，默认不拦截正斜率模型。代码位置为 `src/hybcapseq_quant/cli.py` 中的 `--min-r2` 与 `--enforce-yield-slope-range`，以及 `core.py` 的 `Settings` 和 `calibrate()`。最低标称浓度或拟合区间下界不自动定义为 LoD/LoQ。

得率法以 ERCC 比对碱基占比与其投入质量占比的比值计算单条得率，再以合格 ERCC 得率的算术均值推算目标质量。深度法拟合标准化平均深度与已知提取液等效拷贝浓度的双对数关系。两者分别输出，不取均值。完整公式、单位和回归方向见 [METHODS.md](METHODS.md)。`run` 当前不实现第三种 RPKMF 回归估计，不应在稿件中将其写成已经部署的功能。

原参考的分子质量注释和修订后组装长度继续适用。新增 PCRROI 的默认质量系数保留 `NA`，得率浓度相应留空；深度法可提供带限定标记的局部参考模型估计。只有明确审核了局部参考的物理测量对象、质量赋值和来源，才在 `targets.tsv` 中填写相应参数。局部片段当量不能直接更名为病毒基因组当量，ERCC RNA 校准转移至其他靶标的准确性也不能由高 R² 推断。

分节段目标使用锁定的完整参考组长度；缺失节段不缩短分母。`target_quantification.tsv` 每组只保留一条，适合靶标汇总。源表节段行重复附带的同组浓度不能相加。新增 ROI 各自独立分组，不与原病毒参考合并为一条浓度。

## 5. 回收率、数字 PCR 对照与偏差评估

污水浓度基于提取液浓度、实际体积和 PMMoV 质量回收率换算。回收率是浓缩前后总量的比值，不是简单浓度比；缺失回收率不以 1 替代。PMMoV 作为代理回收指标的适用性应单独讨论，不能假设它对所有病毒完全等效。

将数字 PCR 与测序结果统一到同一样本、同一单位及同一测量层级，记录是否已经完成回收率和稀释校正。识别阴性、未检出、低于定量限及不可计算值，不能强行填成正浓度参与对数回归。重复文库属于技术重复，应按母样本组织训练和留出评价。

```bash
hybcapseq_quant fit-correction \
  --pairs reviewed_dpcr_pairs.tsv --outdir correction_assessment
hybcapseq_quant apply-correction \
  --pairs reviewed_dpcr_pairs.tsv \
  --parameters correction_assessment/correction_parameters.tsv \
  --fit-scope pooled --output corrected_pairs.tsv
```

校正程序按靶标和方法估计乘法系数，保留原始浓度，评价留一母样本及跨批次的适用性。30% APE 是评价目标，不是通过调整系数保证达成的结果。APE、对数偏差、RMSE、R² 与 Lin CCC 应同时报告样本数和有效配对数。v3.4 新比对结果及局部测量对象需要重新评估，不能复用旧的全参考校正参数。字段、纳入条件和评价规则见 [TARGET_CORRECTION.md](TARGET_CORRECTION.md)。

## 6. 投稿与复现记录

方法部分应报告实际软件提交、数据库版本及 SHA256、EsViritu 和 minimap2 版本、参考范围、读段过滤与计数规则、实际体积、ERCC 赋值、拟合标准、外推策略、回收率计算、重复结构及独立验证方式。结果部分区分数学估计、经校正结果和满足验证条件的报告浓度。未完成的真实 FASTQ 端到端验证、LoD/LoQ 或局部校准转移，不写为已经完成。

英文方法文本见 [MANUSCRIPT_METHODS.md](MANUSCRIPT_METHODS.md)。实际实验参数由研究记录填入，不以本仓库示例替代。自动化检验使用明确标注的模拟数据，不作为实验准确性的证据。
