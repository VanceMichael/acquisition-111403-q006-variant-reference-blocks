import {createServer} from 'node:http';
import {openArchive,digest} from './archive.mjs';
import {createClient} from './source.mjs';
import {collect as collectAll, collectVariant, collectShardDelivery} from './reader.mjs';

function sourcesFor(result, client) {
  const entries=[];
  if(!result.manifest) return entries;
  for(const interval of result.manifest.intervals){
    const url=new URL(interval.url,client.origin).href;
    const exchange=client.exchanges.find(item=>item.url===url&&item.status===200);
    if(!exchange) continue;
    entries.push({
      url:interval.url,
      reference:result.manifest.reference,
      chrom:interval.chrom,
      start:interval.start,
      end:interval.end,
      hash:digest(exchange.bytes)
    });
  }
  return entries;
}

function resumeFromArchive(archive) {
  return async (url,manifest)=>{
    const interval=manifest.intervals.find(item=>item.url===url);
    if(!interval) return undefined;
    const match=archive.sources().find(source=>
      source.url===url &&
      source.reference===manifest.reference &&
      source.chrom===interval.chrom &&
      source.start===interval.start &&
      source.end===interval.end
    );
    if(!match) return undefined;
    return archive.raw(match.hash);
  };
}

function shardOptions(target, archive) {
  const options={};
  const maxShards=target.searchParams.get('maxShards');
  if(maxShards!==null){
    if(!/^[0-9]+$/.test(maxShards)){
      const error=new Error('maxShards 必须是非负整数');
      error.status=400;
      throw error;
    }
    options.maxShards=Number(maxShards);
  }
  const expectedReference=target.searchParams.get('reference');
  if(expectedReference!==null) options.expectedReference=expectedReference;
  if(target.searchParams.get('resume')==='archive') options.resume=resumeFromArchive(archive);
  return options;
}

export async function startIntake({origin,directory,port=0}) {
  const archive=await openArchive(directory);
  const server=createServer(async (request,response)=>{
    const send=(value,status=200)=>{response.writeHead(status,{'content-type':'application/json; charset=utf-8'});response.end(JSON.stringify(value));};
    try {
      const target=new URL(request.url,'http://localhost');
      const path=target.pathname;
      if(request.method==='GET'&&path==='/health'){send({status:'运行中'});return;}
      if(request.method==='GET'&&path==='/records'){send(archive.records());return;}
      if(request.method==='GET'&&path==='/receipts'){send(archive.exchanges());return;}
      if(request.method==='GET'&&path==='/sources'){send(archive.sources());return;}
      if(request.method==='POST'&&path==='/collect/variant'){
        const client=createClient(origin);const records=await collectVariant(client);
        send(await archive.commit(records,client.exchanges));return;
      }
      if(request.method==='POST'&&path==='/collect/shards'){
        const client=createClient(origin);
        const options=shardOptions(target,archive);
        const result=await collectShardDelivery(client,options);
        const receipt=await archive.commit(result.records,client.exchanges,sourcesFor(result,client));
        send({...receipt,delivery:result.delivery});return;
      }
      if(request.method==='POST'&&path==='/collect'){
        const client=createClient(origin);
        const options=shardOptions(target,archive);
        const records=await collectAll(client,options);
        const delivery=records.find(record=>record.kind==='delivery');
        const receipt=await archive.commit(records,client.exchanges,sourcesFor({manifest:delivery?.manifest},client));
        send({...receipt,delivery});return;
      }
      send({message:'接口不存在'},404);
    } catch(error){send({message:error.message},error.status===400?400:502);}
  });
  await new Promise(resolve=>server.listen(port,'0.0.0.0',resolve));
  return {origin:`http://127.0.0.1:${server.address().port}`,archive,close:()=>new Promise(resolve=>server.close(resolve))};
}
