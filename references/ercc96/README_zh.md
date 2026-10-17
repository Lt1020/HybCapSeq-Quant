# NIST ERCC96参考序列与EsViritu格式metadata

整理日期：2026-10-17（UTC）。数据包版本：1.0.1。

本包以NIST SRM 2374认证文件中的96个ERCC编号为记录集合，提供真实序列、EsViritu所需16列格式的忠实注释表、官方GenBank编号对应表、逐条序列SHA256、版本差异表及可重建的Python脚本。按指定规则设置`Assembly=Name`、`Asm_length=Length`；其余缺少所选官方来源支持的10个metadata字段，统一写为字面值`NA`。

本包是依据官方资料整理的衍生文件，不是NIST或EsViritu官方发布物，也不是已合并的EsViritu v3.3数据库。

## 选择与来源

两个子目录分别对应两个明确的官方版本，请按版本选用，保留各自的FASTA及metadata配对。

| 目录 | 序列来源和处理 | ERCC数 | 总参考长度 | 模糊位点数 |
| --- | --- | ---: | ---: | ---: |
| `data/certified_v1/` | 官方`SRM2374_Sequence_v1.FASTA`逐字节复制 | 96 | 86,319 nt | 11 |
| `data/ambiguities_resolved_v2/` | 官方`SRM2374_ambiguities_resolved_v2.FASTA`按认证版96个ERCC编号提取，仅重排FASTA换行 | 96 | 86,319 nt | 0 |

主要忠实存档版本是`certified_v1`。所有参考长度为273–2022 nt。序列保留NIST发布的A/C/G/T字母表达；保留原有polyA相关尾部及方向，没有改成U，没有修剪尾部、补齐启动子、追加侧翼或反向互补。

NIST认证的是ERCC control inserts的DNA序列。证书Table 2明确说明insert长度包含polyA尾部；它不是在所有转录条件下完整RNA产物的统一长度。因此本包`Length`表示实际提供的参考序列长度，不把它描述为额外测定的完整RNA产物长度。

NIST的消歧义v2原始FASTA以及GTF各有97条记录；额外的一条是2732 nt的载体参考，编号为`080418_Consensus_Vector_Sequence_NIST_SEQUENCING_ASSEMBLY_noRestrict_rev`。它保留在`sources/`原始文件中，但不进入两套ERCC96参考或其metadata。排除详情见`data/excluded_records.tsv`。

官方资料与文献：

1. [NIST SRM 2374 Certificate of Analysis](https://tsapps.nist.gov/srmext/certificates/2374.pdf)，2017-12-06发行，Table 2在第7–10页；第11页列出序列文件版本。
2. [认证序列v1](https://tsapps.nist.gov/srmext/certificates/documents/SRM2374_Sequence_v1.FASTA)。
3. [官方消歧义序列v2](https://tsapps.nist.gov/srmext/certificates/documents/SRM2374_ambiguities_resolved_v2.FASTA)。
4. [官方GTF v1](https://tsapps.nist.gov/srmext/certificates/documents/SRM2374_ambiguities_resolved_v1.GTF)，用于ERCC编号与GenBank编号对应及长度复核。
5. Lee H et al. *External RNA Controls Consortium Beta Version Update*. J Genomics. 2016;4:19–22. [doi:10.7150/jgen.16082](https://doi.org/10.7150/jgen.16082)。该文区分NIST 96条集合与商业92条混合物；商业集合不包含ERCC-00007、ERCC-00018、ERCC-00023、ERCC-00128。本包保留这四条。文中beta材料的历史说明不用于自动改写NIST认证序列。
6. [EsViritu自定义数据库格式](https://github.com/cmmr/EsViritu/blob/42aad529d9112451531c81633171bc9a7e1f630b/docs/custom-database.md)，固定到已核对的上游提交`42aad529d9112451531c81633171bc9a7e1f630b`。

## metadata字段规则

每个参考目录均包含`virus_pathogen_database.fna`及`virus_pathogen_database.all_metadata.tsv`。文件名沿用EsViritu习惯；本包内容为ERCC参考。

| 字段 | 本包内容及依据 |
| --- | --- |
| `Accession` | 官方ERCC编号，如`ERCC-00002`；作为自定义记录ID，与FASTA首字段完全一致。它不是NCBI accession。 |
| `description` | NIST FASTA原始标题去掉`>`；当前官方标题仅有ERCC编号。 |
| `Name` | 官方ERCC编号。 |
| `Segment` | `NA`。 |
| `kingdom` | `NA`。 |
| `phylum` | `NA`。 |
| `tclass` | `NA`；保留EsViritu的准确列名。 |
| `order` | `NA`。 |
| `family` | `NA`。 |
| `genus` | `NA`。 |
| `species` | `NA`。 |
| `subspecies` | `NA`。 |
| `Length` | 按当前FASTA计数的实际长度；96条均与NIST Table 2及GTF一致。 |
| `TaxID` | `NA`；所用NIST FASTA、GTF与证书Table 2没有提供该字段，不根据序列来源物种推断。 |
| `Assembly` | 与同一行的`Name`完全一致，如`ERCC-00002`；每条ERCC单独作为一个计算分组，属于指定的数据库组织规则，不是NCBI assembly编号。 |
| `Asm_length` | 与同一行的`Length`完全一致，为参考序列实际长度的整数；一个计算分组仅包含一条ERCC参考。 |

列顺序严格为：

```text
Accession description Name Segment kingdom phylum tclass order family genus species subspecies Length TaxID Assembly Asm_length
```

实际文件使用制表符，不是上面展示用的空格。`NA`不是0，也不是空字符串。GenBank编号按NIST GTF与证书Table 2原样记录在`sequence_provenance.tsv`及`data/NIST_Table2_records.tsv`，不自行补`.1`等版本后缀。后者另含官方箱位、insert长度、polyA长度及证书报告的GC fraction；GC值按原文保留，不替换为另一种分母计算的值。

这些表保持所列16列的字段结构与ID对应关系。`Assembly`采用逐条ERCC的名称，`Asm_length`采用对应参考的数值长度，使每条ERCC形成独立的计算分组。上述两个字段由已有名称和长度按指定规则生成；官方序列与来源注释继续原样保留。分类字段及`TaxID`等其余10个缺失字段仍为`NA`。本次验证覆盖记录、分组对应及长度一致性，未构建`.mmi`或执行EsViritu全流程。

## 官方来源间的差异

所有坐标均为FASTA参考的1-based坐标。以下是逐字符比较的结果，不能据此断言应修改哪一份官方资料。

| 记录 | 官方资料中的差异 | 本包处理 |
| --- | --- | --- |
| ERCC-00138，第997位 | 认证FASTA为`R`，官方消歧义v2为`T`；`R`通常代表A/G，因此这一替换不在该符号的候选集合内 | 两版分别原样保留，差异表中标记`false` |
| ERCC-00018 | 认证FASTA唯一的`R`实际在第540位，证书Table 2脚注(e)写为599 | 记录差异，不移动或替换认证序列碱基 |

两版共11处碱基差异，全部列于`data/official_version_differences.tsv`；其中10处符合认证符号的IUPAC候选集合，1处为上述ERCC-00138差异。该表记录源文件之间的差异，不表示本包自行做了11次碱基替换。认证版的11个模糊位点位于7条记录中。

## 复现与核验

Linux上仅需Python 3.10+标准库即可从随包原始文件重建：

```bash
cd ERCC96_NIST_EsViritu
sha256sum -c SHA256SUMS
python3 build_ercc96.py --out rebuilt
diff -r data rebuilt
```

`--out`必须是新目录或空目录。脚本固定核验官方源文件SHA256、96个ERCC编号、序列字符、两版长度、GTF映射、Table 2长度、FASTA/metadata逐条ID和长度对应关系；官方文件发生变化会报错，需人工核对版本后再更新固定值。证书PDF随下载包提供；GitHub中通过URL和SHA256引用，重建不依赖PDF解析软件。

每个`sequence_provenance.tsv`中的SHA256以去换行后的原始序列字符串计算；`SHA256SUMS`以实际文件字节计算，两者用途不同。`data/manifest.json`记录来源、版本、核验结果及两处源资料差异。

原始来源中没有本次实验每条ERCC的实际投入浓度、批号或拷贝数，本包没有填入任何浓度值。
