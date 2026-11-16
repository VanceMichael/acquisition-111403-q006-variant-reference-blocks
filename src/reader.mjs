import {collectShardDelivery} from './shards.mjs';

// 现有单份已整理资料；保留原有行为并单独导出。
export async function collectVariant(client) {
 const document=await client("/variant");
 return [{id:document.reference+':'+document.chrom+':'+document.position+':'+document.ref+':'+document.alt,kind:'variant',type:'variant',...document}];
}

// 分片测序交付：清单 + VCF 分片。记录、覆盖声明与矛盾声明一起返回。
export async function collectShards(client, options={}) {
 const result=await collectShardDelivery(client,options);
 return result.records;
}

// 默认收件同时接纳单份整理资料与分片交付。
export async function collect(client, options={}) {
 const variant=await collectVariant(client);
 const delivery=await collectShardDelivery(client,options);
 return [...variant,...delivery.records];
}

export {collectShardDelivery,collectShardDelivery as collectDelivery};
