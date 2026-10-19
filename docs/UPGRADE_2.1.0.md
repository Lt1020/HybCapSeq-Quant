# HybCapSeq-Quant 2.1.0 与数据库 v3.3.1

软件版本为2.1.0，数据库修订为v3.3.1。原数据库目录可以继续使用原名称，内部文件名保持不变。既有 `samples.tsv`、`ercc_standards.tsv` 和原始EsViritu检出表可以继续使用。

## Linux更新

在原来运行定量的Python/conda环境中安装：

```bash
python -m pip install --upgrade --no-cache-dir "git+https://github.com/Lt1020/HybCapSeq-Quant.git@v2.1.0"
hybcapseq_quant --version
```

将`DB_DIR`设为实际联合参考数据库目录，然后执行：

```bash
hybcapseq_quant upgrade-db --db-dir "$DB_DIR"
hybcapseq_quant validate-db --db-dir "$DB_DIR"
```

升级命令已包含注释安装，无需另外执行`install-references --replace`。如需先查看检查结果，使用`upgrade-db --db-dir "$DB_DIR" --check`。

## 数据库修订依据

`Influenza B virus (B/Rochester/0093/2013)`的8条现存参考序列合计13,859 nt，原metadata声明15,617 nt。v3.3.1对该组8条记录的`Asm_length`及对应分子注释实施一致修正，组级ssRNA理论摩尔质量为4,712,060 g/mol。参考集合长度不等同于对天然完整基因组长度的新认定。

所有Accession、单条参考长度、序列内容和分类注释保留。原FASTA的SHA256为`16a4a1b616d3cc9c45bd43f446aa561a4e6868b8f0b9861819b99cf402a2c8ed`；升级时会核对该指纹及8条成员长度。FASTA与MMI不写入、不重建；metadata仅替换已确认的8个单元格，原列顺序、其他行和换行形式保留。

原始metadata SHA256为`bb406d8b177deede56ba3642fa2c5e7abfda56911f3fb3a91893a52da0a39ea3`，修订后为`8c9040347d70ed79a10d869ebd577c338488fdf8a2d973d7780b93790c13c34c`。修订分子注释SHA256为`6c838cfbbb14bd1de6e3eb04812ae04022d4147a4a7097288c796c537d97aa1a`。

升级先校验、暂存并验证全部修订文件，再备份替换；普通写入异常会恢复已替换文件。`hybcapseq_quant_backups`保存原文件。相同修订重复运行返回`already_current`；未知修改的metadata、FASTA或注释会终止升级，不强行覆盖。原v3.3发行清单属于历史记录，新目录状态由新增的`esviritu_db_v3.3.1_manifest.json`及其`.sha256`清单描述。

## 兼容已完成的比对

旧检出表仍保留`Asm_length=15617`。软件只对公开修正规则中身份、单条长度、Assembly均匹配的8个Accession接受这一历史值，计算改用修订参考长度13,859。原检出字段不改写，`reference_length_nt`表示实际定量分母；`reference_corrections.tsv`记录原值、有效值、来源表和修正编号。其他长度冲突仍报错。

每个修订组输出`reference_metadata_corrected`及`reference_correction_ids`，样本与批次manifest记录相应证据。该修订不需要重做FASTQ比对。

## 模型范围和NA

探索性模式默认显示区间外估计，同时保留`outside_calibration_range`和`within_calibration_range=false`。正式模式默认禁止外推。两种模式分别可用`--no-extrapolation`和`--allow-extrapolation`显式选择；外推不会自动补齐外部验证声明。

ERCC最低浓度及拟合下界不作为LoD或LoQ。已明确提供LoQ时仍按该值判定。缺少质量因子、没有有效曲线、没有有效信号或缺少污水换算参数的相应浓度继续保留NA并给出原因；部分双链DNA没有被强行赋予单双链质量。

## 运行

样本表中的路径、实验投入和回收参数不需要因升级重新填写。使用新的空输出目录：

```bash
nohup hybcapseq_quant run --samples inputs/samples.tsv --ercc-standards inputs/ercc_standards.tsv --db-dir "$DB_DIR" --outdir quantification_results_2_1 --method both --reporting-mode exploratory --wastewater-unit copies/L > quantification_2_1.log 2>&1 &
```

输出字段、两种公式、ERCC赋值及污水回算单位沿用现有定义。原始数据及旧分析结果保留。

## 验证

测试覆盖数据库指纹、只检查模式、8个目标单元格、换行保留、序列/索引保持、原文件备份、重复升级、未知改动拒绝、写入失败恢复、历史结果兼容、非修正规则冲突拒绝、组装质量算术及区间外输出模式。

`scripts/verify_reference_revision.py`可重算已有结果压缩包。它仅从原结果恢复上游检出字段和manifest中的总有效碱基数，不调整任何实验参数；随后核对ERCC曲线和所有非修订组估计是否保持一致。批次复核摘要见`validation/reference_revision_2.1.0.json`。
