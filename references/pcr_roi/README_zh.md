# 文献诊断靶区参考片段合集（2026-10-29扩展版）

本包覆盖12个靶标，整理35套公开检测体系。本次按用户要求，新冠保留N1、N2、E和ORF1ab共4条；其他靶标收录数及原有序列保持上一版不变。去除1组完全相同的序列后，共34条真实NCBI参考片段。其中21条满足双侧各100 bp及所选引物/探针兼容条件，放在根目录主FASTA中。7条边界不足、6条存在所选参考上的引物/探针失配，分别列入异常目录。`all_candidates/`包含全部34条片段。Sapovirus Oka2006与Oka2019在AY237420.2上提取得到完全相同的序列，FASTA仅保存一份，各套体系的文献和比对记录全部保留。

## 文件

| 文件/目录 | 内容 |
|---|---|
| `virus_pathogen_database.fna` | 21条完整100 bp且所选位点兼容的参考片段 |
| `virus_pathogen_database.all_metadata.tsv` | 与主FASTA逐条对应的16列EsViritu格式metadata |
| `all_candidates/` | 全部34条去重候选及对应metadata，包含下列异常 |
| `exceptions/boundary_limited/` | 7条边界不足的真实片段；其中部分同时存在失配 |
| `exceptions/oligo_mismatch/` | 6条双侧100 bp齐全但存在所选位点失配的候选 |
| `literature_assays.tsv` | 35套体系的完整引物/探针、文献、DOI、参考accession及其角色 |
| `reference_regions.tsv` | 各套体系的参考坐标、局部坐标、两端实际延伸长度和QC状态 |
| `oligo_alignment_audit.tsv` | 每条已发表引物/探针在所选参考中的定位、方向、失配数及代表组分 |
| `target_record_mapping.tsv` | 检测体系与FASTA记录的对应关系 |
| `sequence_aliases.tsv` | 完全相同片段的合并记录 |
| `overlapping_reference_windows.tsv` | 同一accession上重叠的不同参考片段 |
| `literature_review_zh.md` | 按靶标组织的检索说明和选取结果 |
| `sources/` | 原始NCBI短片段、固定版本accession、下载URL、分类快照和来源目录 |
| `scripts/` | Python标准库构建、重新下载及独立核验代码 |
| `SHA256SUMS` | 文件校验和 |

## 提取规则

所有坐标均为固定`accession.version`上的**正链、1起始、两端包含**坐标。扩增子外边界定义为正向引物结合位点的最左端a和反向引物结合位点的最右端b；反向引物按反向互补序列定位。提取`[a-100,b+100]`。例如a=203、b=508时，提取103–608，共506 nt；引物已包含在a–b内，不额外叠加一次引物长度。

源记录边界不足时，仅保存`[max(1,a-100),min(source_length,b+100)]`，并从主FASTA中分出。**有accession并不保证该记录覆盖所需的全部延伸范围。** 本包不补N、不拼接其他毒株、不将参考碱基替换为引物序列。全长片段仍逐碱基保持NCBI参考记录原样。

主FASTA的位点兼容条件：同一参考中至少一个已发表正向引物、一个已发表反向引物，以及至少一个内部探针与参考IUPAC兼容且对应参考位点为ACGT。SYBR体系没有探针要求。混合引物中的其他组分可不匹配该代表毒株，逐条结果均在审计表中。已发表脱氧肌苷I仅作为定位通配符，不能据此推断结合效率。失配候选按既有文献靶区附近的最佳比对定位，不标为精确引物匹配记录。

## accession的含义

`design_reference_accessions`只记录能在原始论文中核实的设计模板；未能核实单一设计accession时明确为NA。`published_reference_accessions`记录文献明确给出的坐标、比对或标准物质参考；`reference_role`说明其角色。`source_accession_version`/`extraction_accession_version`是本次实际下载和提取的固定版本。三者不能相互替代。后续应用论文采用的标准物质参考不被写成原始设计模板。

公开文献体系未必等于本次实际数字PCR试剂盒的检测体系；本包没有对商业试剂的引物身份作推断。若用于比较相同物理靶区，应先核对实际数字PCR试剂/引物信息。

## EsViritu字段

metadata严格保持以下16列，使用制表符分隔：

`Accession description Name Segment kingdom phylum tclass order family genus species subspecies Length TaxID Assembly Asm_length`

FASTA第一段标识与`Accession`完全一致。`Length=Asm_length=实际片段长度`，不是原基因组长度。每个唯一参考片段使用独立的`Name`与`Assembly`，二者相同。`Name`是PCRROI记录标签；实际病毒/分类信息见`target_record_mapping.tsv`及分类列。分类使用NCBI Taxonomy快照，未知字段为NA；`tclass`对应NCBI的class。不同文献体系可能取自同一病毒同一基因区域，不能按记录数统计病毒种数。

本包是短靶区**增量参考文件**，未与v3.3全库合并，也未生成MMI索引。与原有全基因组或其他重叠ROI同时比对会出现多重映射；必须定义去重/归属规则，不能直接把各ROI的read_count或浓度相加。RPKMF若针对ROI计算，应使用ROI有效长度，并与原全基因组RPKMF区分。计算定位核验不能代替数据库特异性、定量性能或实验验证。

## 离线重建

在本包解压目录运行：

```bash
python3 scripts/build_reference_package.py --inputs sources --output rebuilt_reference
python3 scripts/validate_reference_package.py --package rebuilt_reference
```

构建程序只读取本包保存的短参考片段和来源文件，不联网、不设计或修改引物。重新向NCBI核对固定版本源片段，可运行`python3 scripts/refresh_source_windows.py --inputs sources --output refreshed_sources`，随后将`--inputs`改为`refreshed_sources`。下载程序检查原版本和已保存的碱基序列，变化时停止，不自动替换为最新版本。

## 各靶标收录情况

| 靶标 | 文献体系数 | 去重片段数 | 主FASTA | 边界不足 | 失配候选 |
|---|---:|---:|---:|---:|---:|
| PMMoV | 1 | 1 | 1 | 0 | 0 |
| DENV-1 | 3 | 3 | 1 | 1 | 1 |
| SARS-CoV-2 | 4 | 4 | 4 | 0 | 0 |
| Influenza_A | 3 | 3 | 1 | 1 | 1 |
| Influenza_B | 3 | 3 | 1 | 1 | 1 |
| Norovirus_GI | 3 | 3 | 2 | 0 | 1 |
| Norovirus_GII | 3 | 3 | 2 | 0 | 1 |
| Human_astrovirus | 3 | 3 | 2 | 1 | 0 |
| Sapovirus | 3 | 2 | 2 | 0 | 0 |
| Rotavirus_A | 3 | 3 | 0 | 3 | 0 |
| MPXV | 3 | 3 | 3 | 0 | 0 |
| CHIKV | 3 | 3 | 2 | 0 | 1 |

## 关键核查结果

- SARS-CoV-2 CDC N1/N2：本次在MN908947.3上的实际匹配坐标较原文表1提前16 nt。提取采用真实序列匹配坐标；原文坐标另存，N1实际扩增子长度为72 nt。
- SARS-CoV-2 China CDC ORF1ab：新增公开ORF1ab体系，保留原有N1、N2和E。MN908947.3上的扩增子为13342–13460（119 nt），双侧各延伸100 bp后为13242–13560（319 nt），三个引物/探针位点均精确匹配。MN908947.3作为本次实际定位参考，未写成原作者唯一设计模板。
- DENV-1 Johnson2005：M87512.1上实际区间8956–9067，较原表坐标提前17 nt。
- Influenza B Chen2011：原文将目标标作NP，但CY040733.1是segment 8 NS1/NEP。本包按实际NCBI记录标注。
- Ward2004：原文毒株表述与所引accession的现有记录名称不同；U52940.1实际为A/Texas/36/1991的340 nt部分M1记录，AF100376.1实际为B/Victoria/2/87。实际记录名称、失配及边界不足均已保留。
- Le Cann2004：L06802.1在2484结束，正好到达反向引物外边界，下游没有100 bp。公开复用表的AV2与探针位置与其长度顺序不符，本次按实际匹配位置纠正。
- Rotavirus A：本次保留的3套参考均有一侧不足100 bp，故均位于边界异常目录，主FASTA中未放入RVA记录。Zeng2008与Jothikumar2009虽均为NSP3体系，但位于不同区域。
- MPXV G2R位于倒置末端重复区，本包只保存一个实际位点；基因拷贝数与基因组拷贝数不可直接混同。

“常用”的筛选依据是公开诊断用途、独立原始研究中的复用及可核查程度，未进行全球引用量排名。PMMoV本次保留1套已核验常用定量体系；这不表示相关文献只有1篇。Zhang2006/Haramoto2013同一靶区的历史沿用未重复计数。详细检索取舍见文献说明文件。
