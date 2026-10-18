# EsViritu v3.3 分子类型与理论质量注释

本参考文件与原始 v3.3 数据库通过精确 Accession 键和 metadata SHA256 绑定。全量注释覆盖 31,003 条参考，使用 NCBI Nuccore 记录和 ICTV 官方 VMR，保留证据、来源差异与缺失原因。

## 文件

| 文件 | 用途 |
|---|---|
| [virus_pathogen_database.molecular_annotations.tsv.gz](../../src/hybcapseq_quant/data/reference/virus_pathogen_database.molecular_annotations.tsv.gz) | 全量注释，解压后文件名去掉 `.gz` |
| [virus_pathogen_database.molecular_annotations.manifest.json](../../src/hybcapseq_quant/data/reference/virus_pathogen_database.molecular_annotations.manifest.json) | 注释、来源及配套 metadata 的校验和与模型定义 |
| [ncbi_accession_summary.jsonl.gz](../../src/hybcapseq_quant/data/reference/ncbi_accession_summary.jsonl.gz) | 每条病毒 Accession 的 NCBI 检索摘要，用于追溯和离线重建 |
| [MOLECULAR_ANNOTATIONS.md](../../docs/MOLECULAR_ANNOTATIONS.md) | 英文字段定义及计算说明 |

推荐通过安装命令将注释直接放入实际数据库目录。先在当前 shell 中将 `DB_DIR` 设置为实际目录：

```bash
hybcapseq_quant install-references --db-dir "$DB_DIR"
hybcapseq_quant validate-db --db-dir "$DB_DIR"
```

最终目录包含原有 `virus_pathogen_database.fna`、`virus_pathogen_database.mmi`、`virus_pathogen_database.all_metadata.tsv`，以及新增的 `.molecular_annotations.tsv` 和 `.molecular_annotations.manifest.json`。两份新增文件必须放在一起，文件名保持不变。不会改写原始三个文件。

## 来源和判定原则

1. 对 30,907 条病毒参考逐条查询 NCBI Nuccore ESummary，严格核对 Accession.version 和长度。ERCC 使用已记录来源的 96 条转录本参考；其中只有 92 条获得商业 Mix 1 浓度赋值。
2. 使用 [ICTV VMR MSL41 v1，2026-07-29](https://ictv.global/sites/default/files/VMR/VMR_MSL41.v1.20260729.xlsx)。先匹配 Accession，再依次匹配 species、genus、family 和 order；分类层级支持的分子化学类型必须一致。匹配层级和证据数量保留在表中。
3. 同时保留 NCBI 序列记录的核酸类型与 ICTV 病毒基因组类型。基因组当量模型优先采用 ICTV 类型。两者不一致时显式记录，不将序列记录中的 DNA/cDNA/复制形式误作样本实际化学组成的证明。
4. 263 条部分双链 DNA 参考保留已知类别 `partially_dsDNA`，因互补链完成程度未知，质量为 NA。11 条参考可确认 RNA，但链型未确认，质量为 NA。不能确认的信息不补填默认值。
5. 单条参考质量可计算，不代表整个 assembly 一定可定量。完整组的长度、分子类型及成员质量均需一致。8 条参考所在的一个组存在上游总长度不一致，两种估计器均输出长度质控标记。

## 质量的含义

采用声明的平均质量模型 `average_reference_position_v1`：ssRNA 340、dsRNA 680、ssDNA 330、dsDNA 660 g/mol/参考位置。其中 dsRNA 的 680 由两条 RNA 链各 340 推导。参考长度只计算一条表示链的位置数，双链系数包含互补链。

理论摩尔质量为 `length_nt × molecular_mass_per_nt_g_mol`；除以阿伏伽德罗常数得到单分子质量 g，再乘 10⁹ 得到 ng。这些值是平均模型下的理论值，不是官方逐条实测分子量，也不表示病毒颗粒质量。使用者如有特定材料形态及其质量证据，可通过有来源的 `--targets` 表显式赋值，并保留相应验证状态。

## 来源重建

正常定量直接使用软件附带的注释，无需联网。需要重建时，安装附加依赖，下载上述指定版本的官方 VMR，并将 `ICTV_VMR` 设置为下载文件，`NCBI_SNAPSHOT` 设置为提供的摘要文件。

```bash
python -m pip install '.[references]'
hybcapseq_quant build-molecular-annotations \
  --db-dir "$DB_DIR" \
  --ictv-vmr "$ICTV_VMR" \
  --ncbi-snapshot "$NCBI_SNAPSHOT" \
  --cache-dir reference_cache \
  --output-dir rebuilt_annotations
```

上述方式使用保存的 NCBI 摘要，避免在线记录更新造成不可重现。省略 `--ncbi-snapshot` 时，程序逐条在线检索并缓存响应，采用共享限速与失败重试。指定 VMR 的 SHA256 必须与声明版本相符。输出目录不得已有同名注释文件；生成后应审阅来源和缺失信息，再用于新的分析。

ICTV 数据来源按其公布的 CC BY 4.0 条款署名。NCBI 摘要记录相应公共数据库 Accession 与来源 URL。质量模型来源与推导见 [METHODS](../../docs/METHODS.md)。
