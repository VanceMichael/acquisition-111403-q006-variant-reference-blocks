# 来源资料约定

群体测序资料接收器。收件包含两类交付：GET /variant 的单个已整理变异对象，以及 GET /shards 列出的参考序列标识与 VCF 文本分片。只做资料接入。

## 交互与身份

GET /variant 返回已整理对象。GET /shards 列出reference、intervals和文件URL；intervals 的 start/end 为闭区间。GET /vcf/:id返回VCF 4.3子集文本。`##INFO`/`##FORMAT`中Number=A/R/G分别按替代等位、全部等位和基因型基数解释，GT中的数字索引当前行REF/ALT，`.`是缺失。`<NON_REF>`配合END为闭区间参考块，不表示某个实际替代等位。分片边界允许重叠；reference不同的分片不能归并。分片头部可以用 `##reference=` 和 `##contig=<ID=,length=>` 声明参考序列与坐标长度；未声明的不作推断。

## 外部变化

正常资料在 `fixtures/source.mjs`。同一来源可能处于这些外部阶段：reference-change、conflicting-block、header-sample-order。集成测试通过来源实例的`setPhase`切换外部状态，这不是线上业务接口。阶段内同一URL的响应可以改变，其余URL保留正常资料；来源的请求轨迹在实例的`requests`中可查，控制方法不对网络开放。未知地址返回404，非法请求返回非成功状态。

## 收件边界

`src/reader.mjs` 的 `collect` 同时收取 /variant 与 /shards 两类交付；`collectVariant` 和 `collectShards` 分别只收其中一类。`POST /collect` 执行一次完整收件，`POST /collect/shards` 只收分片交付，`POST /collect/variant` 只收单个整理对象。

分片交付按以下规则接入，任一不符即拒绝该次分片交付并保留来源原文凭据：

- 核对清单的 reference 与每个分片头部的 `##reference`（若声明）。
- 记录位置必须落在所属分片的闭区间内；参考块 END 也必须在分片区间内。若头部声明 contig 长度，记录末端不得超过该长度。
- 头部字段必须有 `##INFO`/`##FORMAT` 声明，且值的数量与 Number 声明一致；GT 索引不得超过当行 ALT 数量；样本列顺序和数量必须与 `#CHROM` 行一致。
- 多等位记录按替代等位拆成独立记录，`sourceAlleleIndex` 保留来源 ALT 索引。Number=A 取对应值，Number=R 取 REF 与被拆等位，Number=G 不做投影，只按来源基因型顺序保留原值；GT 的原始字符串与来源索引一并保留，`projected` 只把被拆出的等位映射为1，其余替代等位与缺失都标为`.`。
- 每个记录带 `provenance`：分片 URL、分片文本摘要、分片区间、原始行号、原始行文本。重叠分片中的同一记录若内容一致，合并为一条记录并保留多份 provenance；不一致时拆成多条记录，并用 `conflict` 记录显式说明矛盾，不合成确定基因型。
- `coverage` 记录按分片列出显式声明的区间（`explicit`）与区间内未出现记录的区间（`noExplicitRecord`）。这两类都只表示"没有显式声明"，`inference` 固定为 `not-reference`。未读取的分片 `unread=true`、`inference=unknown`，不把未读区域算成参考型。
- `delivery` 记录列出清单重叠、清单空隙、完全重复的区间、头部 reference/contig/样本顺序/字段声明核对结果与告警。

断点恢复：`POST /collect/shards?resume=archive` 会用 `GET /sources` 中同一 reference、同一分片区间的已归档文本替代网络读取；`?maxShards=N` 限制本次从来源读取的分片数，未读分片保留为 unread。归档按响应字节摘要独立存储原文，重开归档时逐份校验；分片来源索引按 URL + reference + 区间 + 摘要去重。

## 既有收件边界

归档按资料标识和内容摘要区分版本，原文按响应字节摘要独立存储，打开归档时逐份校验。来源文件、接收记录和对外查询不共享内存对象。未声明的协议能力不作推断；格式子集以本页和实际来源为准。接入端不加入疾病预测、风险评分或育种评分。
