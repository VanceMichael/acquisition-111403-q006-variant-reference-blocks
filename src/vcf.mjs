// VCF 4.3 子集解析。这里只做资料接入所需的语法、头部声明与数值基数核对，
// 不推断未声明字段，也不做任何临床或育种解释。

const HEADER_COLUMNS = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'];

function fail(message) {
  throw new Error('VCF 格式错误：' + message);
}

function splitStructured(text) {
  const parts = [];
  let current = '';
  let quoted = false;
  for (const character of text) {
    if (character === '"') quoted = !quoted;
    if (character === ',' && !quoted) {
      parts.push(current);
      current = '';
    } else {
      current += character;
    }
  }
  if (quoted) fail('头部声明中的引号未闭合');
  parts.push(current);
  const fields = {};
  for (const part of parts) {
    const separator = part.indexOf('=');
    if (separator <= 0) fail('头部声明字段无效：' + part);
    const key = part.slice(0, separator).trim();
    let value = part.slice(separator + 1).trim();
    if (value.startsWith('"') && value.endsWith('"')) value = value.slice(1, -1);
    if (Object.hasOwn(fields, key)) fail('头部声明字段重复：' + key);
    fields[key] = value;
  }
  return fields;
}

function parseMetaLine(line) {
  const body = line.slice(2);
  const separator = body.indexOf('=');
  if (separator <= 0) fail('头部行无效：' + line);
  const key = body.slice(0, separator);
  const value = body.slice(separator + 1);
  return {key,value};
}

function parseDeclaration(key, value, lineNumber) {
  if (!value.startsWith('<') || !value.endsWith('>')) fail('第 ' + lineNumber + ' 行 ' + key + ' 声明无效');
  const fields = splitStructured(value.slice(1, -1));
  if (!fields.ID) fail('第 ' + lineNumber + ' 行 ' + key + ' 缺少 ID');
  if (!fields.Number) fail('第 ' + lineNumber + ' 行 ' + key + '=' + fields.ID + ' 缺少 Number');
  if (!fields.Type) fail('第 ' + lineNumber + ' 行 ' + key + '=' + fields.ID + ' 缺少 Type');
  if (!/^(?:[0-9]+|[ARG.])$/.test(fields.Number)) fail('第 ' + lineNumber + ' 行 ' + key + '=' + fields.ID + ' Number 无效');
  if (!['Integer','Float','Flag','Character','String'].includes(fields.Type)) fail('第 ' + lineNumber + ' 行 ' + key + '=' + fields.ID + ' Type 无效');
  if (fields.Type === 'Flag' && fields.Number !== '0') fail('第 ' + lineNumber + ' 行 ' + key + '=' + fields.ID + ' Type=Flag 必须 Number=0');
  return {...fields, sourceLine:lineNumber};
}

function parseInfo(text) {
  const values = {};
  if (text === '' || text === '.') return values;
  for (const part of text.split(';')) {
    if (part === '') continue;
    const separator = part.indexOf('=');
    const key = separator === -1 ? part : part.slice(0, separator);
    if (!key) fail('INFO 字段名为空');
    if (Object.hasOwn(values, key)) fail('INFO 字段重复：' + key);
    values[key] = separator === -1 ? '' : part.slice(separator + 1);
  }
  return values;
}

function numberCardinality(declaration, alleleCount) {
  const number = declaration.Number;
  if (number === 'A') return alleleCount;
  if (number === 'R') return alleleCount + 1;
  if (number === 'G') return null;
  if (number === '.') return null;
  return Number(number);
}

function validateTypedValue(value, type, name) {
  if (value === '.') return;
  if (type === 'Integer' && !/^[+-]?[0-9]+$/.test(value)) fail(name + ' 应为整数：' + value);
  if (type === 'Float' && !/^[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$/.test(value)) fail(name + ' 应为浮点数：' + value);
  if (type === 'Character' && [...value].length !== 1) fail(name + ' 应为单个字符：' + value);
}

function validateInfo(info, declarations, alleleCount) {
  for (const [key, raw] of Object.entries(info)) {
    const declaration = declarations[key];
    if (!declaration) fail('INFO 字段未在头部声明：' + key);
    if (declaration.Type === 'Flag') {
      if (raw !== '') fail('INFO 标志字段不应带值：' + key);
      continue;
    }
    if (raw === '') fail('INFO 字段缺少值：' + key);
    if (raw === '.') continue;
    const values = raw.split(',');
    const expected = numberCardinality(declaration, alleleCount);
    if (expected !== null && values.length !== expected) {
      fail('INFO 字段 ' + key + ' 的值的数量应为 ' + expected + '，实际为 ' + values.length);
    }
    for (const value of values) validateTypedValue(value, declaration.Type, 'INFO ' + key);
  }
}

function validateFormat(keys, sampleMap, declarations, alleleCount, sampleOrder) {
  const seen = new Set();
  for (const key of keys) {
    if (seen.has(key)) fail('FORMAT 字段重复：' + key);
    seen.add(key);
    const declaration = declarations[key];
    if (!declaration) fail('FORMAT 字段未在头部声明：' + key);
    if (declaration.Type === 'Flag') fail('FORMAT 字段不能是 Flag：' + key);
    if (key === 'GT') {
      if (declaration.Number !== '1' || declaration.Type !== 'String') fail('GT 必须声明为 Number=1,Type=String');
    }
  }
  for (const sample of sampleOrder) {
    const fields = sampleMap[sample];
    let ploidy = null;
    for (const key of keys) {
      const declaration = declarations[key];
      const raw = fields[key];
      if (raw === undefined) fail('样本 ' + sample + ' 缺少 FORMAT 字段 ' + key);
      if (raw === '' && (declaration.Type !== 'String' || key === 'GT')) fail('样本 ' + sample + ' 的 ' + key + ' 值为空');
      if (key === 'GT') {
        const genotype = validateGenotypeText(raw, alleleCount, sample);
        if (!genotype.alleles.includes('.') && genotype.alleles.length) ploidy = genotype.alleles.length;
      }
    }
    for (const key of keys) {
      if (key === 'GT') continue;
      const declaration = declarations[key];
      const raw = fields[key];
      if (raw === '.') continue;
      const values = raw.split(',');
      if (declaration.Number === '.') continue;
      if (declaration.Number === 'G') {
        if (raw === '.' || ploidy === null) continue;
        const expected = genotypeCount(alleleCount + 1, ploidy);
        if (values.length !== expected) {
          fail('样本 ' + sample + ' 的 ' + key + ' 值的数量应为 ' + expected + '，实际为 ' + values.length);
        }
      } else {
        const expected = numberCardinality(declaration, alleleCount);
        if (expected !== null && values.length !== expected) {
          fail('样本 ' + sample + ' 的 ' + key + ' 值的数量应为 ' + expected + '，实际为 ' + values.length);
        }
      }
      for (const value of values) validateTypedValue(value, declaration.Type, '样本 ' + sample + ' 的 ' + key);
    }
  }
}

export function parseGenotype(text) {
  if (!text) fail('GT 值为空');
  if (text === '.') return {raw:text, alleles:['.'], separators:[]};
  const separators = [...text].filter(character => character === '/' || character === '|');
  if (separators.length && !separators.every(separator => separator === separators[0])) fail('GT 值混用了 / 和 |：' + text);
  const parts = text.split(/[/|]/);
  if (parts.length === 0 || parts.some(part => part === '')) fail('GT 值无效：' + text);
  const alleles = parts.map(part => {
    if (part === '.') return '.';
    if (!/^[0-9]+$/.test(part)) fail('GT 等位基因索引无效：' + text);
    return Number(part);
  });
  return {raw:text, alleles, separators, phased:separators.length > 0 && separators[0] === '|'};
}

export function validateGenotypeText(text, alleleCount, sample='') {
  const genotype = parseGenotype(text);
  for (const allele of genotype.alleles) {
    if (allele !== '.' && allele > alleleCount) fail('样本 ' + sample + ' 的 GT 索引 ' + allele + ' 超出 ALT 数量 ' + alleleCount);
  }
  return genotype;
}

export function genotypeTuples(alleleTotal, ploidy) {
  if (!Number.isInteger(alleleTotal) || alleleTotal < 1) fail('等位基因总数无效：' + alleleTotal);
  if (!Number.isInteger(ploidy) || ploidy < 1) fail('基因型倍数无效：' + ploidy);
  const result = [];
  const combinations = (limit, length, prefix, output) => {
    if (length === 0) {
      output.push(prefix);
      return;
    }
    for (let value = 0; value <= limit; value++) combinations(limit, length - 1, [...prefix, value], output);
  };
  // VCF 的 Number=G 顺序是 0/0, 0/1, 1/1, 0/2, 1/2, 2/2 ……；
  // 对任意倍数，按最大等位基因分组并按非降序枚举。
  for (let maximum = 0; maximum < alleleTotal; maximum++) {
    const prefixes = [];
    combinations(maximum, ploidy - 1, [], prefixes);
    for (const prefix of prefixes) result.push([...prefix, maximum]);
  }
  return result;
}

export function genotypeCount(alleleTotal, ploidy) {
  return genotypeTuples(alleleTotal, ploidy).length;
}

export function parseVcf(text, {reference, shard} = {}) {
  if (typeof text !== 'string') fail('正文不是文本');
  const lines = text.replace(/\r\n/g, '\n').replace(/\r/g, '\n').split('\n');
  while (lines.length && lines[lines.length - 1] === '') lines.pop();
  if (lines.length === 0) fail('文件为空');
  if (!lines[0].startsWith('##fileformat=VCFv4.3')) fail('首行必须声明 ##fileformat=VCFv4.3');

  const infoDeclarations = {};
  const formatDeclarations = {};
  const contigs = new Map();
  let headerReference = null;
  let columns = null;
  let sampleOrder = [];
  const records = [];

  for (let index = 0; index < lines.length; index++) {
    const line = lines[index];
    const lineNumber = index + 1;
    if (line === '') fail('第 ' + lineNumber + ' 行为空行');
    if (line.startsWith('##')) {
      const {key,value} = parseMetaLine(line);
      if (key === 'fileformat') continue;
      if (key === 'INFO' || key === 'FORMAT') {
        const declaration = parseDeclaration(key, value, lineNumber);
        const target = key === 'INFO' ? infoDeclarations : formatDeclarations;
        if (Object.hasOwn(target, declaration.ID)) fail('第 ' + lineNumber + ' 行 ' + key + '=' + declaration.ID + ' 重复声明');
        target[declaration.ID] = declaration;
        continue;
      }
      if (key === 'reference') {
        if (headerReference !== null && headerReference !== value) fail('头部出现多个不同的 reference');
        headerReference = value;
        continue;
      }
      if (key === 'contig') {
        if (!value.startsWith('<') || !value.endsWith('>')) fail('contig 声明无效');
        const fields = splitStructured(value.slice(1, -1));
        if (!fields.ID) fail('contig 缺少 ID');
        if (contigs.has(fields.ID)) fail('contig 重复声明：' + fields.ID);
        contigs.set(fields.ID, fields);
      }
      continue;
    }
    if (line.startsWith('#CHROM')) {
      if (columns !== null) fail('重复的 #CHROM 行');
      columns = line.split('\t');
      if (columns.length < 8) fail('#CHROM 行缺少必需列');
      for (let position = 0; position < HEADER_COLUMNS.length; position++) {
        if (columns[position] !== HEADER_COLUMNS[position]) fail('#CHROM 行第 ' + (position + 1) + ' 列应为 ' + HEADER_COLUMNS[position]);
      }
      if (columns.length === 9) fail('#CHROM 行有 FORMAT 列但没有样本列');
      if (columns.length > 8) {
        if (columns[8] !== 'FORMAT') fail('#CHROM 行第 9 列应为 FORMAT');
        sampleOrder = columns.slice(9);
        const seen = new Set();
        for (const sample of sampleOrder) {
          if (!sample) fail('#CHROM 行存在空样本名');
          if (seen.has(sample)) fail('样本名重复：' + sample);
          seen.add(sample);
        }
      }
      continue;
    }
    if (line.startsWith('#')) fail('第 ' + lineNumber + ' 行头部无效：' + line);
    if (columns === null) fail('第 ' + lineNumber + ' 行数据出现在 #CHROM 行之前');

    const fields = line.split('\t');
    if (fields.length !== columns.length) {
      fail('第 ' + lineNumber + ' 行有 ' + fields.length + ' 列，头部声明了 ' + columns.length + ' 列');
    }
    const [chrom, positionText, id, ref, altText, qual, filter, infoText] = fields;
    if (!chrom) fail('第 ' + lineNumber + ' 行 CHROM 为空');
    if (!/^[1-9][0-9]*$/.test(positionText)) fail('第 ' + lineNumber + ' 行 POS 不是正整数');
    if (!ref || ref === '.') fail('第 ' + lineNumber + ' 行 REF 为空');
    if (!altText || altText === '.') fail('第 ' + lineNumber + ' 行 ALT 为空');
    const alts = altText.split(',');
    const seenAlts = new Set();
    for (const alt of alts) {
      if (!alt) fail('第 ' + lineNumber + ' 行 ALT 存在空值');
      if (seenAlts.has(alt)) fail('第 ' + lineNumber + ' 行 ALT 重复：' + alt);
      seenAlts.add(alt);
    }
    const info = parseInfo(infoText);
    validateInfo(info, infoDeclarations, alts.length);

    let formatKeys = [];
    const sampleMap = {};
    if (fields.length > 8) {
      formatKeys = fields[8].split(':');
      if (!fields[8] || formatKeys.some(key => !key)) fail('第 ' + lineNumber + ' 行 FORMAT 为空或含空字段');
      for (let sampleIndex = 0; sampleIndex < sampleOrder.length; sampleIndex++) {
        const values = fields[9 + sampleIndex].split(':');
        if (values.length !== formatKeys.length) {
          fail('第 ' + lineNumber + ' 行样本 ' + sampleOrder[sampleIndex] + ' 有 ' + values.length + ' 个值，FORMAT 声明了 ' + formatKeys.length + ' 个字段');
        }
        const mapped = {};
        for (let keyIndex = 0; keyIndex < formatKeys.length; keyIndex++) mapped[formatKeys[keyIndex]] = values[keyIndex];
        sampleMap[sampleOrder[sampleIndex]] = mapped;
      }
      validateFormat(formatKeys, sampleMap, formatDeclarations, alts.length, sampleOrder);
    }

    const position = Number(positionText);
    const contig = contigs.get(chrom);
    if (contig && contig.length !== undefined) {
      if (!/^[1-9][0-9]*$/.test(contig.length)) fail('contig ' + chrom + ' 的 length 无效');
      const contigLength = Number(contig.length);
      const recordEnd = info.END !== undefined ? Number(info.END) : position + Math.max(ref.length, 1) - 1;
      if (recordEnd > contigLength) fail('第 ' + lineNumber + ' 行末端 ' + recordEnd + ' 超出 contig ' + chrom + ' 长度 ' + contigLength);
    }
    if (shard) {
      if (chrom !== shard.chrom) fail('第 ' + lineNumber + ' 行染色体 ' + chrom + ' 与分片 ' + shard.chrom + ' 不一致');
      if (position < shard.start || position > shard.end) fail('第 ' + lineNumber + ' 行位置 ' + position + ' 超出分片边界 ' + shard.start + '-' + shard.end);
      if (info.END !== undefined) {
        const end = Number(info.END);
        if (!Number.isInteger(end) || end < position) fail('第 ' + lineNumber + ' 行 END 无效');
        if (end > shard.end) fail('第 ' + lineNumber + ' 行参考块末端 ' + end + ' 超出分片边界 ' + shard.start + '-' + shard.end);
      } else {
        const spanEnd = position + Math.max(ref.length, 1) - 1;
        if (spanEnd > shard.end) fail('第 ' + lineNumber + ' 行记录末端 ' + spanEnd + ' 超出分片边界 ' + shard.start + '-' + shard.end);
      }
    }

    records.push({
      chrom, position, id, ref, alts, qual, filter, info,
      formatKeys, samples:sampleMap, sampleOrder:[...sampleOrder],
      line, lineNumber
    });
  }

  if (columns === null) fail('缺少 #CHROM 行');
  if (reference !== undefined && headerReference !== null && headerReference !== reference) {
    fail('头部 reference=' + headerReference + ' 与分片清单 reference=' + reference + ' 不一致');
  }
  return {
    fileformat:'VCFv4.3',
    reference:headerReference,
    contigs:[...contigs.values()],
    infoDeclarations,
    formatDeclarations,
    sampleOrder,
    columns,
    records
  };
}
