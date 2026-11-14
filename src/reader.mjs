export async function collect(client) {
 const document=await client("/variant");
 return [{id:document.reference+':'+document.chrom+':'+document.position+':'+document.ref+':'+document.alt,...document}];
}
