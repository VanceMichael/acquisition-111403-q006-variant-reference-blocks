# 测序变异文件的参考块与等位基因接收

群体测序资料接收器。读取单个已整理变异对象，以及完整分片测序交付（VCF 文本分片及独立参考序列标识）。只做资料接入。

## 运行

使用 Node.js 22，依赖由 package-lock.json 固定。安装后可执行：

```sh
npm ci --ignore-scripts --no-audit --no-fund
npm test
npm run demo
npm start
```

服务启动本地来源与收件接口，默认端口分别为4011和4010。`POST /collect` 一次领取当前资料（单份整理对象与分片交付）；`POST /collect/variant`只收单份整理对象；`POST /collect/shards`只收分片交付，可加 `?maxShards=N` 限制本次读取数、`?resume=archive` 复用已归档分片。`GET /records`查看结构化记录，`GET /receipts`查看原文凭据，`GET /sources`查看分片来源索引，`GET /health`查看服务状态。`PORT`、`SOURCE_PORT`与`DATA_DIR`可分别设置端口和归档目录。示例领取不会自动遍历其他来源资料。

```sh
docker build -t variant-reference-blocks .
docker run --rm --network none variant-reference-blocks npm test
docker run --rm --network none variant-reference-blocks npm run demo
docker run --rm -p 4010:4010 variant-reference-blocks
```

## 资料组织

`provider/`是独立本地来源，`fixtures/`记录合成资料和来源阶段，`src/vcf.mjs`解析并核对该 VCF 子集的头部声明、基数与边界，`src/shards.mjs`拆分等位、记录来源行、并列覆盖与矛盾声明，`src/reader.mjs`提供单份与分片两种收件，`src/source.mjs`限定来源地址并记录响应字节，`src/archive.mjs`保存按摘要校验的原文、分片来源与资料版本。`test/`覆盖单份收件、分片收件、来源边界、原文完整性、断点恢复和服务联通。运行不依赖外网、数据库服务或真实设备。

完整来源约定见 `docs/source-contract.md`。项目只访问本地合成资料，不携带真实身份、密钥或业务数据。接入端不加入疾病预测、风险评分或育种评分。

## 快速查看

`npm run demo`启动实际来源和收件服务、执行一次收件并打印查询结果及凭据，临时归档在退出时清理。正常完成时间小于90秒；需要保留数据时用`npm start`及独立`DATA_DIR`。演示反映当前产品能力，不代表全部来源资料已被接纳。
