# EsViritu 自定义数据库 v3.3（v3.2.4 + ERCC96）

本包为用户自定义的v3.3扩展库，基于Michael Tisza发布的官方EsViritu数据库v3.2.4追加96条ERCC参考，并重新生成minimap2索引。

解压后目录为`v3.3/`，保持原数据库的三个文件名和16列metadata结构：

- `virus_pathogen_database.fna`
- `virus_pathogen_database.all_metadata.tsv`
- `virus_pathogen_database.mmi`

原数据库30,907条记录，加上96条ERCC，共31,003条。参考总长235,370,650 nt。官方原FASTA和metadata的所有字节均作为前缀保留，新增内容仅追加在末尾；原索引由合并后的FASTA重新生成。

ERCC采用上一轮v2下载链接对应的NIST认证v1参考集合，即`ERCC96_NIST_certified_v1.fna`，原始序列含11个IUPAC模糊位点。下载修订号v2与NIST的消歧义序列v2是不同概念。本包没有改写或消歧义化ERCC序列。

ERCC的`Assembly=Name`，`Asm_length=Length`；一个ERCC记录对应一个独立计算分组。其他缺少来源支持的10个字段仍为`NA`，不推断或编造TaxID和分类。原有病毒的metadata值全部保留。

## Linux使用

```bash
tar -xzf esviritu_db_v3.3.tar.gz
```

之后将EsViritu的`--db`参数指向解压后的`v3.3`目录。数据库路径取解压后目录的实际绝对路径。

索引由minimap2 `2.31-r1302`生成，采用`sr`预设，索引参数为k=21、w=11。重建索引的标准命令为：

```bash
cd v3.3
minimap2 -x sr -d virus_pathogen_database.mmi virus_pathogen_database.fna
```

## 已完成的验证

- 上游447,285,297字节压缩包MD5等于发布方提供的`24d85c1ec3cbffff12e921d2f39c91b2`。
- 原有30,907条FASTA与metadata编号和长度逐条一致。
- 合并后31,003条FASTA与metadata编号和长度逐条一致，无编号重复或冲突。
- 新索引的全部31,003个参考编号及长度与合并FASTA一致。
- 从每条ERCC提取一段150 nt无模糊位点的测试片段，96条均命中各自的ERCC参考。
- 压缩包重新解压读取后，三个文件的SHA256与构建结果完全一致。

上述为数据库与索引验证，未使用本次实验的真实FASTQ执行完整EsViritu检测或定量验证。

每个数据库文件及压缩包的字节数、SHA256、来源、工具版本及核验结果见同批提供的`esviritu_db_v3.3_manifest.json`。

## 来源及署名

上游病毒数据库：Michael Tisza，*Virus Pathogen Database for EsViritu*，v3.2.4，DOI: 10.5281/zenodo.17716199，CC BY 4.0。上游数据原样保留，本衍生库的变更是追加ERCC记录并重建索引。

- 数据来源：https://zenodo.org/records/17716199
- 许可：https://creativecommons.org/licenses/by/4.0/
- EsViritu：https://github.com/cmmr/EsViritu
- ERCC官方认证序列：https://tsapps.nist.gov/srmext/certificates/documents/SRM2374_Sequence_v1.FASTA
- NIST SRM 2374证书：https://tsapps.nist.gov/srmext/certificates/2374.pdf

## GitHub发布与可复现构建

Release标签为`db-v3.3`。数据库压缩包、本文档、构建清单和SHA256文件作为Release资产提供；源码目录记录构建工具、来源和校验值。

发布流程固定Python 3.12.14、zlib 1.3.2及minimap2 2.31，按保存的打包元数据重建原交付压缩包。只有三个数据库文件及完整压缩包的SHA256均与本次交付一致时，才会发布Release。已有不同校验值的同名资产不会自动覆盖。

- Release：https://github.com/Lt1020/HybCapSeq-Quant/releases/tag/db-v3.3
- 压缩包：https://github.com/Lt1020/HybCapSeq-Quant/releases/download/db-v3.3/esviritu_db_v3.3.tar.gz
