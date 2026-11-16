import test from 'node:test';
import assert from 'node:assert/strict';
import {mkdtemp,rm,writeFile} from 'node:fs/promises';
import {join} from 'node:path';
import {tmpdir} from 'node:os';
import {startSource} from '../provider/server.mjs';
import {createClient} from '../src/source.mjs';
import {collect as collectAll, collectShards} from '../src/reader.mjs';
import {collectShardDelivery} from '../src/shards.mjs';
import {openArchive,digest} from '../src/archive.mjs';
import {startIntake} from '../src/app.mjs';
import {parseVcf,parseGenotype,genotypeCount} from '../src/vcf.mjs';

const HEADER = {
  v43: '##fileformat=VCFv4.3',
  gt: '##FORMAT=<ID=GT,Number=1,Type=String,Description="基因型">',
  ad: '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="等位深度">',
  pl: '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="基因型似然">',
  end: '##INFO=<ID=END,Number=1,Type=Integer,Description="参考块末位">',
  ac: '##INFO=<ID=AC,Number=A,Type=Integer,Description="等位计数">'
};

function vcf({header=['v43','gt'],columns=['plant1'],lines=[]}) {
  return [...header.map(key=>HEADER[key]),'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+columns.join('\t'),...lines].join('\n')+'\n';
}

function stubClient(routes, origin='http://source.local') {
  const exchanges=[];
  const client=async path=>{
    if(!Object.hasOwn(routes,path)) {
      const error=new Error('来源返回 404: '+path);
      error.status=404;
      throw error;
    }
    const value=routes[path];
    const text=typeof value==='string'?value:JSON.stringify(value);
    exchanges.push({
      url:new URL(path,origin).href,
      status:200,
      contentType:typeof value==='string'?'text/plain':'application/json',
      bytes:Buffer.from(text)
    });
    return value;
  };
  client.origin=origin;
  client.exchanges=exchanges;
  return client;
}

const manifestOf=(intervals,reference='ref1')=>({reference,intervals});

test('多等位记录按来源等位拆开并保留来源索引',async()=>{
  const client=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:100,url:'/vcf/a'}]),
    '/vcf/a':vcf({header:['v43','gt','ad','pl'],lines:[
      '1\t20\t.\tG\tA,C\t.\tPASS\t.\tGT:AD:PL\t1|2:0,6,7:90,80,70,60,0,50'
    ]})
  });
  const result=await collectShardDelivery(client);
  assert.equal(result.variants.length,2);
  const [a,c]=result.variants;
  assert.equal(a.alt,'A');
  assert.equal(a.sourceAlleleIndex,1);
  assert.equal(c.alt,'C');
  assert.equal(c.sourceAlleleIndex,2);
  assert.deepEqual(a.sourceAltAlleles,['A','C']);
  const sample=a.samples[0];
  assert.equal(sample.genotype.raw,'1|2');
  assert.deepEqual(sample.genotype.sourceAlleles,[1,2]);
  assert.deepEqual(sample.genotype.alleleClasses,['subject','other-alt']);
  assert.equal(sample.genotype.projected,'1|.');
  assert.equal(sample.fields.AD,'0,6');
  assert.deepEqual(sample.fieldProjection.AD.keptSourceIndexes,[0,1]);
  assert.equal(sample.fields.PL,null);
  assert.equal(sample.fieldProjection.PL.status,'not-projected');
  assert.deepEqual(sample.fieldProjection.PL.sourceValues,['90','80','70','60','0','50']);
  const cSample=c.samples[0];
  assert.equal(cSample.genotype.projected,'.|1');
  assert.equal(cSample.fields.AD,'0,7');
  assert.equal(cSample.fieldProjection.AD.alleleIndex,2);
});

test('参考块按闭区间记录且未声明位置不推断为参考型',async()=>{
  const client=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:20,url:'/vcf/a'}]),
    '/vcf/a':vcf({header:['v43','gt','end'],lines:[
      '1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=10\tGT\t0/0'
    ]})
  });
  const result=await collectShardDelivery(client);
  assert.equal(result.blocks.length,1);
  const block=result.blocks[0];
  assert.equal(block.position,1);
  assert.equal(block.end,10);
  assert.equal(block.isReferenceBlock,true);
  assert.equal(block.alt,'<NON_REF>');
  assert.equal(result.coverage[0].explicit[0].start,1);
  assert.equal(result.coverage[0].explicit[0].end,10);
  assert.deepEqual(result.coverage[0].referenceBlockCovered,[{start:1,end:10,claimIds:['reference-block:ref1:1:1-10']}]);
  assert.deepEqual(result.coverage[0].noExplicitRecord,[{start:11,end:20}]);
  assert.equal(result.coverage[0].inference,'not-reference');
  const text=JSON.stringify(result.records);
  assert.ok(!/referenceGenotype|inferredGenotype|genotypeFromCoverage/.test(text));
});

test('参考块边界或语法不合规时拒绝接收',async()=>{
  const outside=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:20,url:'/vcf/a'}]),
    '/vcf/a':vcf({header:['v43','gt','end'],lines:['1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=30\tGT\t0/0']})
  });
  await assert.rejects(()=>collectShardDelivery(outside),/超出分片边界/);

  const wrongGenotype=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:20,url:'/vcf/a'}]),
    '/vcf/a':vcf({header:['v43','gt','end'],lines:['1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=10\tGT\t0/1']})
  });
  await assert.rejects(()=>collectShardDelivery(wrongGenotype),/参考块的 GT/);

  const missingEnd=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:20,url:'/vcf/a'}]),
    '/vcf/a':vcf({header:['v43','gt'],lines:['1\t1\t.\tA\t<NON_REF>\t.\tPASS\t.\tGT\t0/0']})
  });
  await assert.rejects(()=>collectShardDelivery(missingEnd),/缺少 END/);
});

test('落于分片边界之外的记录不被接受',async()=>{
  const client=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:20,url:'/vcf/a'}]),
    '/vcf/a':vcf({lines:['1\t25\t.\tG\tA\t.\tPASS\t.\tGT\t0/1']})
  });
  await assert.rejects(()=>collectShardDelivery(client),/超出分片边界/);
});

test('重叠分片中的同一记录保留两份来源并可回溯原始行',async()=>{
  const source=await startSource();
  try {
    const result=await collectShardDelivery(createClient(source.origin));
    assert.equal(result.variants.length,2);
    const a=result.variants.find(variant=>variant.alt==='A');
    assert.equal(a.provenance.length,2);
    assert.deepEqual(a.provenance.map(item=>item.shardUrl).sort(),['/vcf/a','/vcf/b']);
    for(const entry of a.provenance){
      const raw=await (await fetch(source.origin+entry.shardUrl)).text();
      const lines=raw.replace(/\r\n/g,'\n').split('\n');
      assert.equal(lines[entry.lineNumber-1],entry.line);
      assert.equal(entry.sourceAlleleIndex,1);
    }
    assert.equal(a.provenance[0].shardDigest.length,64);
    assert.equal(result.conflicts.length,0);
  } finally { await source.close(); }
});

test('互相矛盾的覆盖声明分别保留且不合成确定基因型',async()=>{
  const source=await startSource();
  try {
    source.setPhase('conflicting-block');
    const result=await collectShardDelivery(createClient(source.origin));
    assert.equal(result.blocks.length,2);
    assert.ok(result.conflicts.some(conflict=>conflict.conflictKind==='variant-reference-block-overlap'));
    const conflicting=result.records.filter(record=>record.conflicting===true);
    assert.ok(conflicting.length>=2);
    for(const claim of conflicting){
      assert.ok(claim.conflictIds.length>0);
      assert.equal(Object.hasOwn(claim,'synthesizedGenotype'),false);
    }
    const block=result.blocks.find(item=>item.position===20);
    const variant=result.variants.find(item=>item.position===20&&item.alt==='A');
    assert.ok(block);
    assert.ok(variant);
    assert.ok(block.conflictIds.some(id=>variant.conflictIds.includes(id)));
    assert.equal(result.delivery.complete,true);
  } finally { await source.close(); }
});

test('参考序列版本与头部声明不一致时拒绝接收',async()=>{
  const client=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:20,url:'/vcf/a'}],'ref2'),
    '/vcf/a':vcf({lines:['1\t5\t.\tG\tA\t.\tPASS\t.\tGT\t0/1']}).replace('##fileformat=VCFv4.3','##fileformat=VCFv4.3\n##reference=ref1')
  });
  await assert.rejects(()=>collectShardDelivery(client),/与分片清单 reference=ref2 不一致/);

  const changed=await startSource();
  try {
    changed.setPhase('reference-change');
    await assert.rejects(()=>collectShardDelivery(createClient(changed.origin),{expectedReference:'ref1'}),/与预期 ref1 不一致/);
  } finally { await changed.close(); }
});

test('头部字段声明缺失或基数不符时拒绝接收',async()=>{
  const undeclared=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:20,url:'/vcf/a'}]),
    '/vcf/a':vcf({header:['v43','gt'],lines:['1\t5\t.\tG\tA,C\t.\tPASS\tAC=1,2\tGT\t0/1']})
  });
  await assert.rejects(()=>collectShardDelivery(undeclared),/INFO 字段未在头部声明：AC/);

  const wrongNumber=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:20,url:'/vcf/a'}]),
    '/vcf/a':vcf({header:['v43','gt','ad'],lines:['1\t5\t.\tG\tA,C\t.\tPASS\t.\tGT:AD\t0/1:5,6']})
  });
  await assert.rejects(()=>collectShardDelivery(wrongNumber),/AD 值的数量应为 3/);

  const wrongGenotypeIndex=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:20,url:'/vcf/a'}]),
    '/vcf/a':vcf({lines:['1\t5\t.\tG\tA\t.\tPASS\t.\tGT\t0/2']})
  });
  await assert.rejects(()=>collectShardDelivery(wrongGenotypeIndex),/超出 ALT 数量/);
});

test('分片头部样本顺序不同会被显式标出而不重排',async()=>{
  const client=stubClient({
    '/shards':manifestOf([
      {chrom:'1',start:1,end:10,url:'/vcf/a'},
      {chrom:'1',start:11,end:20,url:'/vcf/b'}
    ]),
    '/vcf/a':vcf({columns:['plant1','plant2'],lines:['1\t5\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\t0/0']}),
    '/vcf/b':vcf({columns:['plant2','plant1'],lines:['1\t15\t.\tG\tA\t.\tPASS\t.\tGT\t0/0\t0/1']})
  });
  const result=await collectShardDelivery(client);
  assert.equal(result.delivery.headerChecks.sampleOrder.status,'same-samples-different-order');
  assert.deepEqual(result.delivery.headerChecks.sampleOrder.byShard[0].samples,['plant1','plant2']);
  assert.deepEqual(result.delivery.headerChecks.sampleOrder.byShard[1].samples,['plant2','plant1']);
  assert.ok(result.conflicts.some(conflict=>conflict.conflictKind==='sample-order-differs'));
  const shardA=result.records.find(record=>record.kind==='shard'&&record.url==='/vcf/a');
  assert.deepEqual(shardA.header.sampleOrder,['plant1','plant2']);
});

test('分片头部同名字段声明不同会被标出',async()=>{
  const client=stubClient({
    '/shards':manifestOf([
      {chrom:'1',start:1,end:10,url:'/vcf/a'},
      {chrom:'1',start:11,end:20,url:'/vcf/b'}
    ]),
    '/vcf/a':vcf({header:['v43','gt','ad'],lines:['1\t5\t.\tG\tA\t.\tPASS\t.\tGT:AD\t0/1:5,6']}),
    '/vcf/b':[
      HEADER.v43,
      HEADER.gt,
      '##FORMAT=<ID=AD,Number=A,Type=Integer,Description="等位深度">',
      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tplant1',
      '1\t15\t.\tG\tA\t.\tPASS\t.\tGT:AD\t0/1:6'
    ].join('\n')+'\n'
  });
  const result=await collectShardDelivery(client);
  assert.equal(result.delivery.headerChecks.declarationConflicts.length,1);
  assert.equal(result.delivery.headerChecks.declarationConflicts[0].scope,'FORMAT');
  assert.equal(result.delivery.headerChecks.declarationConflicts[0].ID,'AD');
  assert.ok(result.conflicts.some(conflict=>conflict.conflictKind==='header-declaration-differs'));
});

test('中断后重读不会把未读区域算成参考型',async()=>{
  const source=await startSource();
  try {
    const client=createClient(source.origin);
    const partial=await collectShardDelivery(client,{maxShards:1});
    assert.equal(partial.delivery.complete,false);
    assert.equal(partial.delivery.unread.length,1);
    const unread=partial.coverage.find(item=>item.unread);
    assert.equal(unread.inference,'unknown');
    assert.deepEqual(unread.noExplicitRecord,[{start:20,end:30}]);
    assert.match(unread.note,/不推断为参考型/);
    assert.ok(partial.blocks.length===1);

    const aText=await (await fetch(source.origin+'/vcf/a')).text();
    const resumed=await collectShardDelivery(client,{resume:new Map([['/vcf/a',aText]])});
    assert.equal(resumed.delivery.complete,true);
    assert.deepEqual(resumed.delivery.resume.reusedUrls,['/vcf/a']);
    assert.ok(resumed.coverage.every(item=>item.unread===false));
    assert.equal(resumed.blocks.length,1);
  } finally { await source.close(); }
});

test('清单的空隙、重复区间与重叠被区分记录',async()=>{
  const gaps=stubClient({
    '/shards':manifestOf([
      {chrom:'1',start:1,end:20,url:'/vcf/a'},
      {chrom:'1',start:25,end:30,url:'/vcf/b'}
    ]),
    '/vcf/a':vcf({lines:['1\t5\t.\tG\tA\t.\tPASS\t.\tGT\t0/1']}),
    '/vcf/b':vcf({lines:['1\t27\t.\tG\tA\t.\tPASS\t.\tGT\t0/1']})
  });
  const result=await collectShardDelivery(gaps);
  assert.deepEqual(result.delivery.intervalGaps,[{chrom:'1',start:21,end:24}]);
  assert.equal(result.delivery.overlaps.length,0);

  const duplicate=stubClient({
    '/shards':manifestOf([
      {chrom:'1',start:1,end:20,url:'/vcf/a'},
      {chrom:'1',start:1,end:20,url:'/vcf/a'},
      {chrom:'1',start:15,end:25,url:'/vcf/a'}
    ]),
    '/vcf/a':vcf({lines:['1\t18\t.\tG\tA\t.\tPASS\t.\tGT\t0/1']})
  });
  const duplicated=await collectShardDelivery(duplicate);
  assert.equal(duplicated.delivery.duplicateIntervals.length,1);
  assert.equal(duplicated.delivery.overlaps.length,3);
  assert.ok(duplicated.delivery.overlaps.every(overlap=>overlap.sameUrl===true));
  assert.ok(duplicated.delivery.overlaps.some(overlap=>overlap.start===15&&overlap.end===20));
});

test('归档保存分片原文与来源标识供断点恢复',async()=>{
  const directory=await mkdtemp(join(tmpdir(),'shard-archive-'));
  const source=await startSource();
  try {
    const client=createClient(source.origin);
    const result=await collectShardDelivery(client);
    const archive=await openArchive(directory);
    const sources=result.manifest.intervals.map(interval=>{
      const exchange=client.exchanges.find(item=>item.url===source.origin+interval.url);
      return {...interval,reference:result.manifest.reference,hash:digest(exchange.bytes)};
    });
    await archive.commit(result.records,client.exchanges,sources);
    assert.equal(archive.sources().length,2);
    const reopened=await openArchive(directory);
    assert.deepEqual(reopened.sources(),archive.sources());
    for(const sourceEntry of reopened.sources()){
      const bytes=await reopened.raw(sourceEntry.hash);
      assert.equal(digest(bytes),sourceEntry.hash);
      assert.ok(bytes.toString('utf8').startsWith('##fileformat=VCFv4.3'));
    }
    await writeFile(join(directory,'objects',archive.sources()[0].hash),'被修改');
    await assert.rejects(()=>openArchive(directory),/校验失败/);
  } finally {
    await source.close();
    await rm(directory,{recursive:true,force:true});
  }
});

test('接入服务可收件分片、查询覆盖与矛盾并与原文凭据对应',async()=>{
  const directory=await mkdtemp(join(tmpdir(),'shard-service-'));
  const source=await startSource();
  try {
    source.setPhase('conflicting-block');
    const intake=await startIntake({origin:source.origin,directory});
    try {
      source.setPhase('normal');
      const response=await fetch(intake.origin+'/collect/shards',{method:'POST'});
      const receipt=await response.json();
      assert.equal(response.status,200,JSON.stringify(receipt));
      assert.equal(receipt.delivery.complete,true);
      const records=await (await fetch(intake.origin+'/records')).json();
      assert.ok(records.some(record=>record.kind==='coverage'));
      assert.ok(records.some(record=>record.kind==='reference-block'));
      const receipts=await (await fetch(intake.origin+'/receipts')).json();
      assert.ok(receipts.some(item=>item.url.endsWith('/vcf/a')));
      assert.ok(receipts.some(item=>item.url.endsWith('/vcf/b')));
      const sources=await (await fetch(intake.origin+'/sources')).json();
      assert.equal(sources.length,2);
      assert.ok(sources.every(item=>/^[a-f0-9]{64}$/.test(item.hash)));

      source.setPhase('conflicting-block');
      const second=await (await fetch(intake.origin+'/collect/shards',{method:'POST'})).json();
      assert.equal(second.delivery.complete,true);
      const conflicts=await (await fetch(intake.origin+'/records')).json();
      assert.ok(conflicts.some(record=>record.kind==='conflict'));
    } finally { await intake.close(); }
  } finally {
    await source.close();
    await rm(directory,{recursive:true,force:true});
  }
});

test('接入服务完成一次默认收件即包含原有单份资料和分片交付',async()=>{
  const directory=await mkdtemp(join(tmpdir(),'shard-default-'));
  const source=await startSource();
  try {
    const intake=await startIntake({origin:source.origin,directory});
    try {
      const response=await fetch(intake.origin+'/collect',{method:'POST'});
      const receipt=await response.json();
      assert.equal(response.status,200,JSON.stringify(receipt));
      const records=await (await fetch(intake.origin+'/records')).json();
      assert.ok(records.some(record=>record.id==='ref1:1:12:A:C'));
      assert.ok(records.some(record=>record.kind==='variant'&&record.alt==='A'));
      assert.ok(records.some(record=>record.kind==='delivery'));
      assert.ok(records.some(record=>record.kind==='reference-block'));
      assert.ok(records.every(record=>typeof record.id==='string'&&record.id.length));
    } finally { await intake.close(); }
  } finally {
    await source.close();
    await rm(directory,{recursive:true,force:true});
  }
});

test('接入服务的参数与摘要校验对错误输入明确拒绝',async()=>{
  const directory=await mkdtemp(join(tmpdir(),'shard-errors-'));
  const source=await startSource();
  try {
    const intake=await startIntake({origin:source.origin,directory});
    try {
      const bad=await fetch(intake.origin+'/collect/shards?maxShards=-1',{method:'POST'});
      assert.equal(bad.status,400);
      const wrongReference=await fetch(intake.origin+'/collect/shards?reference=ref9',{method:'POST'});
      assert.equal(wrongReference.status,502);
      assert.match((await wrongReference.json()).message,/预期 ref9/);
      const missing=await fetch(intake.origin+'/collect/unknown',{method:'POST'});
      assert.equal(missing.status,404);
      assert.equal((await (await fetch(intake.origin+'/sources')).json()).length,0);
    } finally { await intake.close(); }
  } finally {
    await source.close();
    await rm(directory,{recursive:true,force:true});
  }
});

test('接入记录不包含疾病预测或育种评分',async()=>{
  const source=await startSource();
  try {
    const records=await collectAll(createClient(source.origin));
    const text=JSON.stringify(records);
    for(const forbidden of ['disease','riskScore','breedingScore','prediction','致病','育种评分','表型预测']){
      assert.ok(!text.includes(forbidden),forbidden);
    }
    const allowed=['delivery','shard','coverage','reference-block','variant','conflict'];
    assert.ok(records.every(record=>allowed.includes(record.kind)));
  } finally { await source.close(); }
});

test('单等位记录保留全部计数字段且拒绝同一位点不同 REF',async()=>{
  const single=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:20,url:'/vcf/a'}]),
    '/vcf/a':vcf({header:['v43','gt','ad','pl'],lines:['1\t5\t.\tG\tA\t.\tPASS\t.\tGT:AD:PL\t0/1:4,5:30,0,40']})
  });
  const result=await collectShardDelivery(single);
  assert.equal(result.variants.length,1);
  const sample=result.variants[0].samples[0];
  assert.equal(sample.fields.AD,'4,5');
  assert.equal(sample.fields.PL,'30,0,40');
  assert.deepEqual(sample.fieldProjection,{});
  assert.equal(sample.genotype.projected,'0/1');
  assert.equal(sample.genotype.genotypeKnown,true);

  const differentRef=stubClient({
    '/shards':manifestOf([
      {chrom:'1',start:1,end:20,url:'/vcf/a'},
      {chrom:'1',start:10,end:30,url:'/vcf/b'}
    ]),
    '/vcf/a':vcf({lines:['1\t15\t.\tG\tA\t.\tPASS\t.\tGT\t0/1']}),
    '/vcf/b':vcf({lines:['1\t15\t.\tT\tA\t.\tPASS\t.\tGT\t0/1']})
  });
  const mixed=await collectShardDelivery(differentRef);
  assert.equal(mixed.variants.length,2);
  assert.ok(mixed.conflicts.some(conflict=>conflict.conflictKind==='ref-differs'));
  assert.ok(mixed.variants.every(variant=>variant.conflicting));
});

test('contig 声明长度可核对记录末端', async()=>{
  const valid=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:100,url:'/vcf/a'}]),
    '/vcf/a':[
      HEADER.v43,
      '##contig=<ID=1,length=100>',
      HEADER.gt,
      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tplant1',
      '1\t99\t.\tG\tA\t.\tPASS\t.\tGT\t0/1'
    ].join('\n')+'\n'
  });
  const result=await collectShardDelivery(valid);
  assert.equal(result.delivery.headerChecks.contigs[0].ID,'1');
  assert.deepEqual(result.delivery.headerChecks.contigs[0].lengths,['100']);

  const beyond=stubClient({
    '/shards':manifestOf([{chrom:'1',start:1,end:100,url:'/vcf/a'}]),
    '/vcf/a':[
      HEADER.v43,
      '##contig=<ID=1,length=50>',
      HEADER.gt,
      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tplant1',
      '1\t60\t.\tG\tA\t.\tPASS\t.\tGT\t0/1'
    ].join('\n')+'\n'
  });
  await assert.rejects(()=>collectShardDelivery(beyond),/超出 contig 1 长度/);
});

test('解析器拒绝缺失文件格式、列数不符和未声明FORMAT字段',async()=>{
  assert.throws(()=>parseVcf('1\t1\t.\tA\tC\t.\tPASS\t.'),/fileformat/);
  assert.throws(()=>parseVcf(vcf({lines:['1\t1\t.\tA\tC\t.\tPASS\t.']})),/列/);
  assert.throws(()=>parseVcf(vcf({header:['v43'],lines:['1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1']})),/FORMAT 字段未在头部声明：GT/);
  assert.throws(()=>parseVcf(vcf({lines:['1\t1\t.\tA\tC\t.\tPASS\t.\tGT:XX\t0/1:x']})),/FORMAT 字段未在头部声明：XX/);
  assert.throws(()=>parseVcf(vcf({columns:['a','a'],lines:['1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/1']})),/样本名重复/);
});

test('基因型按等位总数与倍数给出VCF顺序的基因型个数',()=>{
  assert.equal(genotypeCount(2,2),3);
  assert.equal(genotypeCount(3,2),6);
  assert.equal(genotypeCount(3,1),3);
  assert.deepEqual(parseGenotype('1|2').sourceAlleles,undefined);
  assert.deepEqual(parseGenotype('1|2').alleles,[1,2]);
  assert.equal(parseGenotype('1|2').phased,true);
  assert.equal(parseGenotype('0/1').phased,false);
});
