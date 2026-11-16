import {digest} from './archive.mjs';
import {genotypeTuples, parseGenotype, parseVcf} from './vcf.mjs';

// 分片交付只做资料接入：保留声明、明确未声明区域、暴露矛盾，
// 不推断参考型、不做疾病预测或育种评分。

const SHORT = value => digest(Buffer.from(String(value))).slice(0, 12);
const READ_NOTE = '分片区间内未出现记录的位置只表示没有显式声明，不推断为参考型';
const UNREAD_NOTE = '该分片尚未读取，区间状态未知，不推断为参考型';

function invalid(message) {
  throw new Error('分片交付格式错误：' + message);
}

function canonicalIntervals(manifest) {
  return manifest.intervals.map(interval => ({
    chrom:interval.chrom,
    start:interval.start,
    end:interval.end,
    url:interval.url
  }));
}

export function validateManifest(manifest) {
  if (!manifest || typeof manifest !== 'object' || Array.isArray(manifest)) invalid('清单必须是对象');
  if (typeof manifest.reference !== 'string' || !manifest.reference) invalid('清单缺少非空 reference');
  if (!Array.isArray(manifest.intervals) || manifest.intervals.length === 0) invalid('清单缺少 intervals');
  const intervals = [];
  for (const [index,item] of manifest.intervals.entries()) {
    const label = 'intervals[' + index + ']';
    if (!item || typeof item !== 'object' || Array.isArray(item)) invalid(label + ' 必须是对象');
    if (typeof item.chrom !== 'string' || !item.chrom) invalid(label + ' 缺少非空 chrom');
    if (!Number.isInteger(item.start) || item.start < 1) invalid(label + ' 的 start 必须是正整数');
    if (!Number.isInteger(item.end) || item.end < item.start) invalid(label + ' 的 end 必须不小于 start');
    if (typeof item.url !== 'string' || !item.url) invalid(label + ' 缺少非空 url');
    intervals.push({chrom:item.chrom,start:item.start,end:item.end,url:item.url});
  }
  return {reference:manifest.reference,intervals};
}

function intersect(left, right) {
  const start = Math.max(left.start, right.start);
  const end = Math.min(left.end, right.end);
  return start <= end ? {start,end} : null;
}

function overlapsOf(intervals) {
  const result = [];
  for (let leftIndex = 0; leftIndex < intervals.length; leftIndex++) {
    for (let rightIndex = leftIndex + 1; rightIndex < intervals.length; rightIndex++) {
      const left = intervals[leftIndex];
      const right = intervals[rightIndex];
      if (left.chrom !== right.chrom) continue;
      const overlap = intersect(left, right);
      if (!overlap) continue;
      result.push({
        chrom:left.chrom,
        start:overlap.start,
        end:overlap.end,
        left:{start:left.start,end:left.end,url:left.url},
        right:{start:right.start,end:right.end,url:right.url},
        sameUrl:left.url === right.url
      });
    }
  }
  return result;
}

function exactDuplicateIntervals(intervals) {
  const seen = new Set();
  const duplicates = [];
  for (const interval of intervals) {
    const key = JSON.stringify(interval);
    if (seen.has(key)) duplicates.push({...interval});
    else seen.add(key);
  }
  return duplicates;
}

function mergeSpans(spans) {
  const sorted = spans.map(span => ({start:span.start,end:span.end}))
    .sort((left,right)=>left.start-right.start || left.end-right.end);
  const merged = [];
  for (const span of sorted) {
    const last = merged[merged.length - 1];
    if (last && span.start <= last.end + 1) last.end = Math.max(last.end, span.end);
    else merged.push({...span});
  }
  return merged;
}

function subtractSpan(interval, covered) {
  const gaps = [];
  let cursor = interval.start;
  for (const span of mergeSpans(covered)) {
    const start = Math.max(span.start, interval.start);
    const end = Math.min(span.end, interval.end);
    if (end < interval.start || start > interval.end) continue;
    if (start > cursor) gaps.push({start:cursor,end:Math.min(start-1,interval.end)});
    cursor = Math.max(cursor, end + 1);
  }
  if (cursor <= interval.end) gaps.push({start:cursor,end:interval.end});
  return gaps.filter(gap => gap.start <= gap.end);
}

function intervalGaps(intervals) {
  const byChrom = new Map();
  for (const interval of intervals) {
    if (!byChrom.has(interval.chrom)) byChrom.set(interval.chrom, []);
    byChrom.get(interval.chrom).push(interval);
  }
  const gaps = [];
  for (const [chrom, list] of byChrom) {
    const merged = mergeSpans(list);
    for (let index = 1; index < merged.length; index++) {
      if (merged[index].start > merged[index-1].end + 1) {
        gaps.push({chrom,start:merged[index-1].end+1,end:merged[index].start-1});
      }
    }
  }
  return gaps;
}

function alleleClass(allele, subjectIndex) {
  if (allele === '.') return 'missing';
  if (allele === 0) return 'ref';
  if (allele === subjectIndex) return 'subject';
  return 'other-alt';
}

function projectGenotype(raw, subjectIndex) {
  if (raw === undefined) {
    return {
      raw:null,
      absent:true,
      phased:false,
      genotypeKnown:false,
      ploidy:null,
      sourceAlleles:[],
      alleleClasses:[],
      subjectAlleleIndexes:[],
      projected:null,
      alleleMap:{'0':'0',[String(subjectIndex)]:'1'},
      note:'该样本在本行没有 GT 字段，未推断基因型'
    };
  }
  const parsed = parseGenotype(raw);
  const separator = parsed.separators[0] || '';
  const mapped = parsed.alleles.map(allele => alleleClass(allele, subjectIndex) === 'subject' ? '1' : allele === 0 ? '0' : '.');
  const subjectIndexes = parsed.alleles
    .map((allele,index)=>allele === subjectIndex ? index : null)
    .filter(index=>index !== null);
  return {
    raw,
    phased:Boolean(parsed.phased),
    genotypeKnown:parsed.alleles.length > 0 && !parsed.alleles.includes('.'),
    absent:false,
    ploidy:parsed.alleles.length,
    sourceAlleles:parsed.alleles,
    alleleClasses:parsed.alleles.map(allele=>alleleClass(allele, subjectIndex)),
    subjectAlleleIndexes:subjectIndexes,
    projected:parsed.alleles.length ? (separator ? mapped.join(separator) : mapped[0]) : '',
    alleleMap:{'0':'0',[String(subjectIndex)]:'1'},
    note:'projected 只把被拆出的等位映射为 1，其余替代等位与缺失均标为 .；raw 保留原始索引'
  };
}

function projectValue(declaration, raw, alleleCount, alleleIndex, ploidy, label) {
  const number = declaration.Number;
  const alleleTotal = alleleCount + 1;
  if (alleleCount <= 1 || (number !== 'A' && number !== 'R' && number !== 'G')) {
    return {value:raw,projection:{number,status:'copied',label}};
  }
  if (raw === '.') {
    return {value:'.',projection:{number,status:'missing',label,alleleIndex}};
  }
  const values = raw.split(',');
  if (number === 'A') {
    return {
      value:values[alleleIndex-1],
      projection:{
        number,status:'sliced',alleleIndex,sourceValues:values,
        keptSourceIndexes:[alleleIndex],label
      }
    };
  }
  if (number === 'R') {
    return {
      value:[values[0],values[alleleIndex]].join(','),
      projection:{
        number,status:'sliced',alleleIndex,sourceValues:values,
        keptSourceIndexes:[0,alleleIndex],
        droppedSourceIndexes:values.map((value,index)=>index).filter(index=>index !== 0 && index !== alleleIndex),
        label
      }
    };
  }
  const tuples = ploidy ? genotypeTuples(alleleTotal, ploidy) : [];
  const entries = tuples.map((tuple,index)=>({
    sourceAlleles:tuple,
    subjectAlleles:tuple.map(allele=>allele === 0 ? 0 : allele === alleleIndex ? 1 : '.'),
    value:values[index] === undefined ? null : values[index]
  }));
  return {
    value:null,
    projection:{
      number:'G',status:'not-projected',alleleIndex,ploidy,sourceValues:values,entries,label,
      reason:'Number=G 的值依赖全部等位基因与倍数；未计算新的标量值，仅按来源基因型顺序保留'
    }
  };
}

function projectInfo(info, declarations, alleleCount, alleleIndex) {
  const values = {};
  const projection = {};
  for (const [key, raw] of Object.entries(info)) {
    const projected = projectValue(declarations[key], raw, alleleCount, alleleIndex, null, 'INFO/' + key);
    values[key] = projected.value;
    if (projected.projection.status !== 'copied') projection[key] = projected.projection;
  }
  return {values,projection};
}

function projectSample(name, order, fields, formatKeys, declarations, alleleCount, alleleIndex) {
  const genotype = projectGenotype(fields.GT, alleleIndex);
  const ploidy = genotype.genotypeKnown ? genotype.ploidy : null;  const values = {};
  const projection = {};
  for (const key of formatKeys) {
    if (key === 'GT') continue;
    const result = projectValue(declarations[key], fields[key], alleleCount, alleleIndex, ploidy, 'FORMAT/' + key);
    values[key] = result.value;
    if (result.projection.status !== 'copied') projection[key] = result.projection;
  }
  return {name,order,genotype,fields:values,fieldProjection:projection};
}

function provenanceOf(shard, record, alleleIndex) {
  const entry = {
    shardUrl:shard.interval.url,
    shardDigest:shard.hash,
    reference:shard.parsed.reference,
    interval:{chrom:shard.interval.chrom,start:shard.interval.start,end:shard.interval.end},
    lineNumber:record.lineNumber,
    line:record.line
  };
  if (alleleIndex !== undefined) {
    entry.sourceAlleleIndex = alleleIndex;
    entry.sourceAltAlleles = [...record.alts];
    entry.sourceRef = record.ref;
    entry.sourcePosition = record.position;
  }
  return entry;
}

function variantClaim({record,shard,alleleIndex,manifest}) {
  const info = projectInfo(record.info, shard.parsed.infoDeclarations, record.alts.length, alleleIndex);
  const samples = record.sampleOrder.map((name,order)=>projectSample(
    name,order,record.samples[name],record.formatKeys,shard.parsed.formatDeclarations,record.alts.length,alleleIndex
  ));
  const spanEnd = record.info.END !== undefined ? Number(record.info.END) : record.position + Math.max(record.ref.length,1) - 1;
  return {
    kind:'variant',
    type:'variant',
    reference:manifest.reference,
    chrom:record.chrom,
    position:record.position,
    start:record.position,
    ref:record.ref,
    alt:record.alts[alleleIndex-1],
    alleleIndex,
    sourceAlleleIndex:alleleIndex,
    sourceAltCount:record.alts.length,
    sourceAltAlleles:[...record.alts],
    end:spanEnd,
    isReferenceBlock:false,
    referenceBlock:false,
    qual:record.qual,
    filter:record.filter,
    info:info.values,
    infoProjection:info.projection,
    formatKeys:[...record.formatKeys],
    sampleOrder:[...record.sampleOrder],
    samples,
    conflicting:false,
    conflictIds:[],
    provenance:[provenanceOf(shard,record,alleleIndex)]
  };
}

function blockClaim({record,shard,manifest}) {
  const samples = record.sampleOrder.map((name,order)=>({
    name,
    order,
    genotype:record.samples[name].GT === undefined
      ? {raw:null,absent:true,note:'该样本在本行没有 GT 字段'}
      : parseGenotype(record.samples[name].GT),
    fields:Object.fromEntries(record.formatKeys.filter(key=>key !== 'GT').map(key=>[key,record.samples[name][key]])),
    fieldProjection:{},
    note:'参考块样本字段按来源原样保留；参考块不表示某个实际替代等位'
  }));
  return {
    kind:'reference-block',
    type:'reference-block',
    reference:manifest.reference,
    chrom:record.chrom,
    position:record.position,
    start:record.position,
    ref:record.ref,
    alt:'<NON_REF>',
    end:Number(record.info.END),
    closed:true,
    isReferenceBlock:true,
    referenceBlock:true,
    qual:record.qual,
    filter:record.filter,
    info:{...record.info},
    infoProjection:{},
    formatKeys:[...record.formatKeys],
    sampleOrder:[...record.sampleOrder],
    samples,
    conflicting:false,
    conflictIds:[],
    provenance:[provenanceOf(shard,record)]
  };
}

function claimsOfShard(shard, manifest) {
  const claims = [];
  for (const record of shard.parsed.records) {
    if (record.alts.includes('<NON_REF>')) {
      if (record.alts.length !== 1) invalid('第 ' + record.lineNumber + ' 行 <NON_REF> 只能作为唯一替代等位');
      if (record.info.END === undefined) invalid('第 ' + record.lineNumber + ' 行参考块缺少 END');
      for (const name of record.sampleOrder) {
        if (record.samples[name].GT === undefined) continue;
        const genotype = parseGenotype(record.samples[name].GT);
        if (genotype.alleles.some(allele => allele !== 0 && allele !== '.')) {
          invalid('第 ' + record.lineNumber + ' 行参考块的 GT 只能使用 0 或 .');
        }
      }
      claims.push(blockClaim({record,shard,manifest}));
      continue;
    }
    for (let alleleIndex = 1; alleleIndex <= record.alts.length; alleleIndex++) {
      claims.push(variantClaim({record,shard,alleleIndex,manifest}));
    }
  }
  return claims;
}

function genotypeSemantic(genotype) {
  return {
    raw:genotype.raw ?? null,
    absent:Boolean(genotype.absent),
    phased:Boolean(genotype.phased),
    ploidy:genotype.ploidy ?? null,
    alleleClasses:genotype.alleleClasses || []
  };
}

function fieldSemantic(sample) {
  const result = {};
  for (const [key,value] of Object.entries(sample.fields)) {
    const projection = sample.fieldProjection[key];
    if (projection && projection.status === 'missing') {
      result[key] = {status:'missing',value};
    } else if (projection && projection.number === 'G') {
      result[key] = {
        status:'not-projected',
        entries:projection.entries
          .map(entry=>({subjectAlleles:entry.subjectAlleles,value:entry.value}))
          .sort((left,right)=>JSON.stringify(left).localeCompare(JSON.stringify(right)))
      };
    } else {
      result[key] = {status:projection ? projection.status : 'copied',value};
    }
  }
  return result;
}

function sampleSemantic(sample) {
  return {name:sample.name,genotype:genotypeSemantic(sample.genotype),fields:fieldSemantic(sample)};
}

function claimSignature(claim) {
  return JSON.stringify({
    kind:claim.kind,
    reference:claim.reference,
    chrom:claim.chrom,
    position:claim.position,
    ref:claim.ref,
    alt:claim.alt,
    end:claim.end,
    qual:claim.qual,
    filter:claim.filter,
    info:claim.info,
    formatKeys:claim.formatKeys,
    sampleOrder:claim.sampleOrder,
    samples:claim.samples.map(sample=>({...sampleSemantic(sample),order:sample.order}))
  });
}

function mergeClaims(rawClaims) {
  const groups = new Map();
  for (const claim of rawClaims) {
    const signature = claimSignature(claim);
    const existing = groups.get(signature);
    if (existing) {
      for (const entry of claim.provenance) {
        const key = JSON.stringify(entry);
        if (!existing.provenanceKeys.has(key)) {
          existing.provenanceKeys.add(key);
          existing.claim.provenance.push(entry);
        }
      }
      if (claim.kind === 'variant') existing.sourceAltSets.add(JSON.stringify(claim.sourceAltAlleles));
      continue;
    }
    groups.set(signature,{
      claim,
      signature,
      sourceAltSets:claim.kind === 'variant' ? new Set([JSON.stringify(claim.sourceAltAlleles)]) : new Set(),
      provenanceKeys:new Set(claim.provenance.map(entry=>JSON.stringify(entry)))
    });
  }
  return [...groups.values()].map(({claim,signature,sourceAltSets})=>({
    ...claim,
    signature,
    sourceAltSets:[...sourceAltSets]
  }));
}

function compareClaims(left, right) {
  return left.chrom.localeCompare(right.chrom) || left.position - right.position ||
    (left.end || 0) - (right.end || 0) || String(left.alt).localeCompare(String(right.alt)) ||
    (left.ref || '').localeCompare(right.ref || '');
}

function assignIds(groups, baseId) {
  const byBase = new Map();
  for (const group of groups) {
    const base = baseId(group);
    if (!byBase.has(base)) byBase.set(base,[]);
    byBase.get(base).push(group);
  }
  for (const [base, list] of byBase) {
    list.sort((left,right)=>left.signature.localeCompare(right.signature));
    for (const [index,group] of list.entries()) {
      group.id = list.length === 1 ? base : base + ':' + SHORT(group.signature) + ':' + index;
    }
  }
}

function sampleKey(claim, {sorted=true}={}) {
  const samples = claim.samples.map(sampleSemantic);
  if (sorted) samples.sort((left,right)=>left.name.localeCompare(right.name));
  return JSON.stringify(samples);
}

function orderKey(claim) {
  return claim.sampleOrder.join('\u0000');
}

function metadataKey(claim) {
  return JSON.stringify({
    chrom:claim.chrom,position:claim.position,ref:claim.ref,alt:claim.alt,end:claim.end,
    qual:claim.qual,filter:claim.filter,info:claim.info
  });
}

function conflictId(conflict) {
  return 'conflict:' + conflict.reference + ':' + conflict.chrom + ':' + conflict.start + '-' + conflict.end +
    ':' + conflict.conflictKind + ':' + SHORT(JSON.stringify({claims:conflict.claims,sources:conflict.sources,message:conflict.message}));
}

function makeConflict({manifest,conflictKind,chrom,start,end,claims=[],sources=[],message}) {
  const conflict = {
    kind:'conflict',
    conflictKind,
    reference:manifest.reference,
    chrom,
    start,
    end,
    claims:claims.map(claim=>claim.id).sort(),
    sources,
    message,
    synthesis:'未合成确定基因型；矛盾声明分别保留'
  };
  conflict.id = conflictId(conflict);
  return conflict;
}

function pairConflictKind(claims) {
  if (claims.length < 2) return null;
  const sampleKeys = new Set(claims.map(claim=>sampleKey(claim)));
  const orderKeys = new Set(claims.map(orderKey));
  const metadataKeys = new Set(claims.map(metadataKey));
  if (sampleKeys.size > 1) return claims[0].kind === 'reference-block' ? 'reference-block-data-differs' : 'sample-data-differs';
  if (orderKeys.size > 1) return 'sample-order-differs';
  if (metadataKeys.size > 1) return 'record-metadata-differs';
  return null;
}

function detectConflicts({variants,blocks,manifest}) {
  const conflicts = [];
  const variantsByPosition = new Map();
  for (const variant of variants) {
    const key = JSON.stringify([variant.chrom,variant.position]);
    if (!variantsByPosition.has(key)) variantsByPosition.set(key,[]);
    variantsByPosition.get(key).push(variant);
  }
  for (const group of variantsByPosition.values()) {
    const refs = [...new Set(group.map(variant=>variant.ref))];
    if (refs.length < 2) continue;
    const first = group[0];
    conflicts.push(makeConflict({
      manifest,conflictKind:'ref-differs',chrom:first.chrom,start:first.position,end:first.position,
      claims:group,sources:refs,
      message:'同一位点在重叠分片中的 REF 不同，可能来自不同参考序列版本，未合并'
    }));
  }

  const variantBySite = new Map();
  for (const variant of variants) {
    const key = JSON.stringify([variant.chrom,variant.position,variant.ref]);
    if (!variantBySite.has(key)) variantBySite.set(key,[]);
    variantBySite.get(key).push(variant);
  }
  for (const group of variantBySite.values()) {
    const first = group[0];
    const altSets = [...new Set(group.flatMap(variant=>variant.sourceAltSets))];
    if (altSets.length > 1) {
      conflicts.push(makeConflict({
        manifest,conflictKind:'alt-set-differs',chrom:first.chrom,start:first.position,end:first.position,
        claims:group,sources:altSets.map(set=>JSON.parse(set)),message:'同一位点的重叠分片声明了不同的 ALT 集合，未合并'
      }));
    }
    const byAlt = new Map();
    for (const variant of group) {
      if (!byAlt.has(variant.alt)) byAlt.set(variant.alt,[]);
      byAlt.get(variant.alt).push(variant);
    }
    for (const [alt,claims] of byAlt) {
      const kind = pairConflictKind(claims);
      if (!kind) continue;
      conflicts.push(makeConflict({
        manifest,conflictKind:kind,chrom:first.chrom,start:first.position,end:first.position,
        claims,sources:[],
        message:kind === 'sample-data-differs'
          ? '同一等位在重叠分片中的样本基因型或计数字段不一致，未择一采信'
          : kind === 'sample-order-differs'
            ? '同一等位的样本顺序在重叠分片间不同，按文件各自顺序保留'
            : '同一等位的注释或记录字段在重叠分片间不同，未合并'
      }));
    }
  }

  const blockByKey = new Map();
  for (const block of blocks) {
    const key = JSON.stringify([block.chrom,block.position,block.end]);
    if (!blockByKey.has(key)) blockByKey.set(key,[]);
    blockByKey.get(key).push(block);
  }
  for (const group of blockByKey.values()) {
    const kind = pairConflictKind(group);
    if (!kind) continue;
    const first = group[0];
    conflicts.push(makeConflict({
      manifest,conflictKind:kind,chrom:first.chrom,start:first.position,end:first.end,
      claims:group,sources:[],
      message:kind === 'sample-order-differs'
        ? '同一参考块的样本顺序在重叠分片间不同，按文件各自顺序保留'
        : '同一参考块在重叠分片中的声明不同，未合并'
    }));
  }

  for (let leftIndex = 0; leftIndex < blocks.length; leftIndex++) {
    for (let rightIndex = leftIndex + 1; rightIndex < blocks.length; rightIndex++) {
      const left = blocks[leftIndex];
      const right = blocks[rightIndex];
      if (left.chrom !== right.chrom) continue;
      if (left.position === right.position && left.end === right.end) continue;
      const overlap = intersect({start:left.position,end:left.end},{start:right.position,end:right.end});
      if (!overlap) continue;
      conflicts.push(makeConflict({
        manifest,conflictKind:'reference-block-overlap',chrom:left.chrom,start:overlap.start,end:overlap.end,
        claims:[left,right],sources:[],message:'参考块边界重叠但范围不同，未合并为一个块'
      }));
    }
  }

  for (const variant of variants) {
    for (const block of blocks) {
      if (variant.chrom !== block.chrom) continue;
      const overlap = intersect({start:variant.position,end:variant.end},{start:block.position,end:block.end});
      if (!overlap) continue;
      conflicts.push(makeConflict({
        manifest,conflictKind:'variant-reference-block-overlap',chrom:variant.chrom,start:overlap.start,end:overlap.end,
        claims:[variant,block],sources:[],message:'同一位置既有变异声明又有参考块声明，二者分别保留'
      }));
    }
  }

  const seen = new Set();
  const unique = [];
  for (const conflict of conflicts) {
    if (seen.has(conflict.id)) continue;
    seen.add(conflict.id);
    unique.push(conflict);
  }
  unique.sort((left,right)=>left.chrom.localeCompare(right.chrom) || left.start-right.start ||
    left.end-right.end || left.conflictKind.localeCompare(right.conflictKind) || left.id.localeCompare(right.id));
  return unique;
}

function markConflicts(variants, blocks, conflicts) {
  const byId = new Map([...variants,...blocks].map(claim=>[claim.id,claim]));
  for (const conflict of conflicts) {
    if (!conflict.claims.length) continue;
    for (const claimId of conflict.claims) {
      const claim = byId.get(claimId);
      if (!claim) continue;
      claim.conflicting = true;
      if (!claim.conflictIds.includes(conflict.id)) claim.conflictIds.push(conflict.id);
    }
  }
}

function coverageOfShard({shard,manifest,claims}) {
  const {chrom,start,end,url} = shard.interval;
  if (shard.unread) {
    return {
      id:'coverage:' + manifest.reference + ':' + chrom + ':' + start + '-' + end + ':unread:' + SHORT(url),
      kind:'coverage',
      reference:manifest.reference,
      chrom,
      interval:{start,end},
      url,
      unread:true,
      error:shard.error,
      explicit:[],
      referenceBlockCovered:[],
      noExplicitRecord:[{start,end}],
      inference:'unknown',
      note:UNREAD_NOTE
    };
  }
  const spans = claims.map(claim=>({
    start:claim.position,
    end:claim.end,
    kind:claim.kind,
    claimIds:[claim.id]
  })).sort((left,right)=>left.start-right.start || left.end-right.end);
  return {
    id:'coverage:' + manifest.reference + ':' + chrom + ':' + start + '-' + end + ':' + shard.hash,
    kind:'coverage',
    reference:manifest.reference,
    chrom,
    interval:{start,end},
    url,
    unread:false,
    error:null,
    explicit:spans,
    referenceBlockCovered:spans
      .filter(span=>span.kind === 'reference-block')
      .map(span=>({start:span.start,end:span.end,claimIds:[...span.claimIds]})),
    noExplicitRecord:subtractSpan({start,end},spans),
    inference:'not-reference',
    note:READ_NOTE
  };
}

function declarationsOf(parsed, scope) {
  const source = scope === 'INFO' ? parsed.infoDeclarations : parsed.formatDeclarations;
  return Object.values(source).map(item=>({ID:item.ID,Number:item.Number,Type:item.Type}));
}

function headerChecksOf(shards, manifest) {
  const read = shards.filter(shard=>!shard.unread);
  const referenceValues = read.map(shard=>({
    url:shard.interval.url,
    interval:{start:shard.interval.start,end:shard.interval.end},
    declared:shard.parsed.reference,
    matches:shard.parsed.reference === null ? null : shard.parsed.reference === manifest.reference
  }));
  const orderKeyOf = shard=>shard.parsed.sampleOrder.join('\u0000');
  const namesKeyOf = shard=>[...shard.parsed.sampleOrder].sort().join('\u0000');
  const orderKeys = new Set(read.map(orderKeyOf));
  const nameKeys = new Set(read.map(namesKeyOf));
  let sampleOrderStatus = 'consistent';
  if (orderKeys.size > 1) sampleOrderStatus = nameKeys.size === 1 ? 'same-samples-different-order' : 'sample-sets-differ';

  const declarationConflicts = [];
  for (const scope of ['INFO','FORMAT']) {
    const seen = new Map();
    for (const shard of read) {
      for (const declaration of declarationsOf(shard.parsed,scope)) {
        if (!seen.has(declaration.ID)) seen.set(declaration.ID,[]);
        const entry = seen.get(declaration.ID);
        const value = declaration.Number + '/' + declaration.Type;
        const existing = entry.find(item=>item.value === value);
        if (existing) {
          if (!existing.urls.includes(shard.interval.url)) existing.urls.push(shard.interval.url);
        } else {
          entry.push({value,urls:[shard.interval.url]});
        }
      }
    }
    for (const [ID,values] of seen) {
      if (values.length > 1) declarationConflicts.push({scope,ID,declarations:values});
    }
  }

  const contigChecks = [];
  const contigs = new Map();
  for (const shard of read) {
    for (const contig of shard.parsed.contigs) {
      if (!contigs.has(contig.ID)) contigs.set(contig.ID,[]);
      contigs.get(contig.ID).push({url:shard.interval.url,contig});
    }
  }
  for (const [ID,entries] of contigs) {
    const definitions = [...new Set(entries.map(entry=>JSON.stringify(entry.contig)))];
    const values = entries.map(entry=>entry.contig);
    const lengths = [...new Set(values.map(contig=>contig.length).filter(length=>length !== undefined))];
    contigChecks.push({ID,consistent:definitions.length === 1,lengths,definitions:definitions.map(value=>JSON.parse(value))});
  }
  contigChecks.sort((left,right)=>left.ID.localeCompare(right.ID));

  return {
    reference:{
      manifest:manifest.reference,
      declaredEverywhere:referenceValues.length > 0 && referenceValues.every(item=>item.declared !== null),
      values:referenceValues
    },
    contigs:contigChecks,
    sampleOrder:{
      status:sampleOrderStatus,
      byShard:read.map(shard=>({
        url:shard.interval.url,
        interval:{start:shard.interval.start,end:shard.interval.end},
        samples:[...shard.parsed.sampleOrder]
      }))
    },
    declarationConflicts
  };
}

function headerConflictRecords({manifest,headerChecks}) {
  const conflicts = [];
  if (headerChecks.sampleOrder.status !== 'consistent') {
    const sourceShards = headerChecks.sampleOrder.byShard;
    const chroms = [...new Set(sourceShards.map(item=>manifest.intervals.find(interval=>interval.url === item.url)?.chrom).filter(Boolean))];
    conflicts.push(makeConflict({
      manifest,
      conflictKind:headerChecks.sampleOrder.status === 'sample-sets-differ' ? 'sample-set-differs' : 'sample-order-differs',
      chrom:chroms.length === 1 ? chroms[0] : '*',
      start:Math.min(...manifest.intervals.map(interval=>interval.start)),
      end:Math.max(...manifest.intervals.map(interval=>interval.end)),
      claims:[],
      sources:sourceShards.map(item=>({url:item.url,interval:item.interval,samples:item.samples})),
      message:headerChecks.sampleOrder.status === 'sample-sets-differ'
        ? '分片头部的样本集合不同，未按位置对齐'
        : '分片头部的样本顺序不同，未重排；各文件按自身顺序保留'
    }));
  }
  for (const conflict of headerChecks.declarationConflicts) {
    conflicts.push(makeConflict({
      manifest,conflictKind:'header-declaration-differs',chrom:'*',
      start:Math.min(...manifest.intervals.map(interval=>interval.start)),
      end:Math.max(...manifest.intervals.map(interval=>interval.end)),
      claims:[],sources:conflict.declarations,
      message:conflict.scope + '/' + conflict.ID + ' 在不同分片头部的 Number/Type 声明不同'
    }));
  }
  for (const check of headerChecks.contigs) {
    if (check.consistent) continue;
    conflicts.push(makeConflict({
      manifest,conflictKind:'contig-declaration-differs',chrom:check.ID,
      start:Math.min(...manifest.intervals.map(interval=>interval.start)),
      end:Math.max(...manifest.intervals.map(interval=>interval.end)),
      claims:[],sources:check.definitions,
      message:'contig ' + check.ID + ' 在不同分片头部的声明不同'
    }));
  }
  return conflicts;
}

function deliveryRecord({manifest,shards,unread,overlaps,gaps,duplicates,headerChecks,resumeUsed,maxShards}) {
  const payload = {
    reference:manifest.reference,
    intervals:canonicalIntervals(manifest),
    shards:shards.map(shard=>({
      url:shard.interval.url,
      interval:{start:shard.interval.start,end:shard.interval.end},
      chrom:shard.interval.chrom,
      digest:shard.hash,
      status:shard.unread ? 'unread' : 'read'
    }))
  };
  const warnings = [];
  if (unread.length) warnings.push('存在未读取分片，相关区域状态未知');
  if (gaps.length) warnings.push('分片清单区间之间存在未覆盖空隙');
  if (duplicates.length) warnings.push('分片清单存在完全重复的区间');
  if (overlaps.length) warnings.push('分片清单存在重叠区间，重叠处可能重复声明同一记录');
  if (!headerChecks.reference.declaredEverywhere) warnings.push('并非每个分片头部都声明了 reference；已核对的值均与清单一致');
  if (headerChecks.sampleOrder.status !== 'consistent') warnings.push('分片头部样本顺序或样本集合不一致');
  if (headerChecks.declarationConflicts.length) warnings.push('分片头部同名字段声明不一致');
  return {
    id:'delivery:' + manifest.reference + ':' + SHORT(JSON.stringify(payload)),
    kind:'delivery',
    reference:manifest.reference,
    manifest:{reference:manifest.reference,intervals:canonicalIntervals(manifest)},
    complete:unread.length === 0,
    shardCount:shards.length,
    readShardCount:shards.filter(shard=>!shard.unread).length,
    unread:unread.map(shard=>({
      chrom:shard.interval.chrom,
      start:shard.interval.start,
      end:shard.interval.end,
      url:shard.interval.url,
      reason:shard.error || '本次未读取，区间状态未知'
    })),
    overlaps,
    intervalGaps:gaps,
    duplicateIntervals:duplicates,
    headerChecks,
    checks:{
      reference:{
        manifest:manifest.reference,
        allDeclaredValuesMatch:headerChecks.reference.values.every(item=>item.matches !== false),
        declaredEverywhere:headerChecks.reference.declaredEverywhere
      },
      boundaries:{
        intervalGaps:gaps.length,
        duplicateIntervals:duplicates.length,
        overlaps:overlaps.length,
        allRecordsWithinIntervals:true
      },
      headers:{
        declarationConflicts:headerChecks.declarationConflicts.length,
        contigConflicts:headerChecks.contigs.filter(check=>!check.consistent).length
      },
      samples:{status:headerChecks.sampleOrder.status},
      reads:{complete:unread.length === 0,unreadShards:unread.length}
    },
    warnings,
    resume:{
      used:Boolean(resumeUsed),
      maxShards:Number.isFinite(maxShards) ? maxShards : null,
      reusedUrls:shards.filter(shard=>shard.source === 'archived').map(shard=>shard.interval.url)
    },
    inference:'not-reference',
    note:unread.length
      ? '交付未完整读取；未读区域不推断为参考型'
      : '所有分片已读取；分片区间内未出现记录的位置不推断为参考型'
  };
}

function publicClaim(claim) {
  const {signature,sourceAltSets,...rest}=claim;
  return rest;
}

export async function collectShardDelivery(client, options={}) {
  const maxShards = options.maxShards === undefined ? Infinity : options.maxShards;
  if (maxShards !== Infinity && (!Number.isInteger(maxShards) || maxShards < 0)) {
    invalid('maxShards 必须是非负整数');
  }
  if (options.resume !== undefined && !(options.resume instanceof Map) && typeof options.resume !== 'function') {
    invalid('resume 必须是 Map 或函数');
  }
  const manifest = validateManifest(await client('/shards'));
  if (options.expectedReference !== undefined && manifest.reference !== options.expectedReference) {
    invalid('来源 reference=' + manifest.reference + ' 与预期 ' + options.expectedReference + ' 不一致');
  }

  for (const interval of manifest.intervals) {
    const target = new URL(interval.url, client.origin);
    if (target.origin !== client.origin) invalid('分片地址超出来源范围：' + interval.url);
  }

  const resumeMap = new Map();
  if (options.resume instanceof Map) {
    for (const [url, value] of options.resume) resumeMap.set(url, value);
  } else if (typeof options.resume === 'function') {
    for (const url of new Set(manifest.intervals.map(interval=>interval.url))) {
      const value = await options.resume(url, manifest);
      if (value !== undefined && value !== null) resumeMap.set(url, value);
    }
  }
  const resumeUsed = resumeMap.size > 0;

  const texts = new Map();
  const fetchErrors = new Map();
  const remaining = new Set();
  const uniqueUrls = [...new Set(manifest.intervals.map(interval=>interval.url))];
  let fetched = 0;
  for (const url of uniqueUrls) {
    if (resumeMap.has(url)) {
      const value = resumeMap.get(url);
      const text = Buffer.isBuffer(value) ? value.toString('utf8') : value;
      if (typeof text !== 'string') invalid('断点恢复内容必须是文本：' + url);
      texts.set(url,text);
      continue;
    }
    if (fetched >= maxShards) {
      remaining.add(url);
      continue;
    }
    try {
      const text = await client(url);
      if (typeof text !== 'string') invalid('分片必须是 VCF 文本：' + url);
      texts.set(url,text);
      fetched++;
    } catch (error) {
      if (error && error.message && error.message.includes('超出来源范围')) throw error;
      fetchErrors.set(url,error.message);
    }
  }

  const shards = [];
  const claims = [];
  for (const interval of manifest.intervals) {
    if (!texts.has(interval.url)) {
      shards.push({
        interval,
        url:interval.url,
        unread:true,
        error:fetchErrors.get(interval.url) || (remaining.has(interval.url) ? '受本次读取上限限制未读取' : '未读取'),
        hash:null,
        source:'unread',
        parsed:null
      });
      continue;
    }
    const text = texts.get(interval.url);
    const hash = digest(Buffer.from(text));
    const parsed = parseVcf(text,{reference:manifest.reference,shard:interval});
    const shard = {
      interval,
      url:interval.url,
      unread:false,
      error:null,
      hash,
      source:resumeMap.has(interval.url) ? 'archived' : 'source',
      parsed
    };
    shards.push(shard);
    claims.push(...claimsOfShard(shard,manifest));
  }

  const uniqueShards = [];
  const shardKeys = new Set();
  for (const shard of shards) {
    const key = shard.interval.chrom + ':' + shard.interval.start + '-' + shard.interval.end + ':' +
      shard.url + ':' + (shard.hash || 'unread');
    if (shardKeys.has(key)) continue;
    shardKeys.add(key);
    uniqueShards.push(shard);
  }

  const variantGroups = mergeClaims(claims.filter(claim=>claim.kind === 'variant'))
    .sort(compareClaims);
  const blockGroups = mergeClaims(claims.filter(claim=>claim.kind === 'reference-block'))
    .sort(compareClaims);
  assignIds(variantGroups, group=>'variant:' + group.reference + ':' + group.chrom + ':' +
    group.position + ':' + group.ref + ':' + group.alt);
  assignIds(blockGroups, group=>'reference-block:' + group.reference + ':' + group.chrom + ':' +
    group.position + '-' + group.end);
  variantGroups.sort(compareClaims);
  blockGroups.sort(compareClaims);

  const headerChecks = headerChecksOf(uniqueShards,manifest);
  const conflicts = detectConflicts({variants:variantGroups,blocks:blockGroups,manifest});
  conflicts.push(...headerConflictRecords({manifest,headerChecks}));
  conflicts.sort((left,right)=>left.chrom.localeCompare(right.chrom) || left.start-right.start ||
    left.end-right.end || left.conflictKind.localeCompare(right.conflictKind) || left.id.localeCompare(right.id));
  markConflicts(variantGroups,blockGroups,conflicts);

  const coverage = [];
  for (const shard of uniqueShards) {
    const shardClaims = [...variantGroups,...blockGroups].filter(claim=>claim.provenance.some(provenance=>
      provenance.shardUrl === shard.url &&
      provenance.interval.chrom === shard.interval.chrom &&
      provenance.interval.start === shard.interval.start &&
      provenance.interval.end === shard.interval.end
    ));
    coverage.push(coverageOfShard({shard,manifest,claims:shardClaims}));
  }
  coverage.sort((left,right)=>left.chrom.localeCompare(right.chrom) || left.interval.start-right.interval.start ||
    left.interval.end-right.interval.end);

  const shardRecords = uniqueShards.map(shard=>({
    id:'shard:' + manifest.reference + ':' + shard.interval.chrom + ':' + shard.interval.start + '-' +
      shard.interval.end + ':' + (shard.hash || 'unread'),
    kind:'shard',
    reference:manifest.reference,
    chrom:shard.interval.chrom,
    interval:{start:shard.interval.start,end:shard.interval.end},
    url:shard.url,
    status:shard.unread ? 'unread' : 'read',
    source:shard.source,
    digest:shard.hash,
    error:shard.error,
    header:shard.parsed ? {
      fileformat:shard.parsed.fileformat,
      reference:shard.parsed.reference,
      contigs:shard.parsed.contigs,
      sampleOrder:[...shard.parsed.sampleOrder],
      info:declarationsOf(shard.parsed,'INFO'),
      format:declarationsOf(shard.parsed,'FORMAT')
    } : null,
    recordCount:shard.parsed ? shard.parsed.records.length : null,
    inference:shard.unread ? 'unknown' : 'not-reference',
    note:shard.unread ? UNREAD_NOTE : READ_NOTE
  }));

  const unread = uniqueShards.filter(shard=>shard.unread);
  const delivery = deliveryRecord({
    manifest,
    shards:uniqueShards,
    unread,
    overlaps:overlapsOf(manifest.intervals),
    gaps:intervalGaps(manifest.intervals),
    duplicates:exactDuplicateIntervals(manifest.intervals),
    headerChecks,
    resumeUsed,
    maxShards
  });

  const records = [
    delivery,
    ...shardRecords,
    ...coverage,
    ...blockGroups.map(publicClaim),
    ...variantGroups.map(publicClaim),
    ...conflicts
  ];
  return {
    records,
    delivery,
    manifest,
    shards:uniqueShards.map(shard=>({...shard,parsed:undefined})),
    variants:variantGroups.map(publicClaim),
    blocks:blockGroups.map(publicClaim),
    coverage,
    conflicts
  };
}
