export const seed={
  "reference": "ref1",
  "chrom": "1",
  "position": 12,
  "ref": "A",
  "alt": "C",
  "samples": {
    "plant1": "0/1"
  }
};
export const seedPath="/variant";
export const routes=[
  {
    "path": "/shards",
    "body": {
      "reference": "ref1",
      "intervals": [
        {
          "chrom": "1",
          "start": 1,
          "end": 20,
          "url": "/vcf/a"
        },
        {
          "chrom": "1",
          "start": 20,
          "end": 30,
          "url": "/vcf/b"
        }
      ]
    },
    "method": "GET",
    "status": 200,
    "query": {},
    "headers": {}
  },
  {
    "path": "/vcf/a",
    "body": "##fileformat=VCFv4.3\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"参考块末位\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"基因型\">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"等位深度\">\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"基因型似然\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tplant1\n1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=10\tGT\t0/0\n1\t20\t.\tG\tA,C\t.\tPASS\t.\tGT:AD:PL\t1|2:0,6,7:90,80,70,60,0,50\n",
    "method": "GET",
    "status": 200,
    "query": {},
    "headers": {}
  },
  {
    "path": "/vcf/b",
    "body": "##fileformat=VCFv4.3\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"基因型\">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"等位深度\">\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"基因型似然\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tplant1\n1\t20\t.\tG\tA,C\t.\tPASS\t.\tGT:AD:PL\t1|2:0,6,7:90,80,70,60,0,50\n",
    "method": "GET",
    "status": 200,
    "query": {},
    "headers": {}
  }
];
export const phases={
  "reference-change": [
    {
      "path": "/shards",
      "body": {
        "reference": "ref2",
        "intervals": [
          {
            "chrom": "1",
            "start": 1,
            "end": 20,
            "url": "/vcf/a"
          }
        ]
      },
      "method": "GET",
      "status": 200,
      "query": {},
      "headers": {}
    }
  ],
  "conflicting-block": [
    {
      "path": "/vcf/b",
      "body": "##fileformat=VCFv4.3\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"参考块末位\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"基因型\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tplant1\n1\t20\t.\tG\t<NON_REF>\t.\tPASS\tEND=30\tGT\t0/0\n",
      "method": "GET",
      "status": 200,
      "query": {},
      "headers": {}
    }
  ],
  "header-sample-order": [
    {
      "path": "/vcf/b",
      "body": "##fileformat=VCFv4.3\n##reference=ref1\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"基因型\">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"等位深度\">\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"基因型似然\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tplant2\tplant1\n1\t20\t.\tG\tA,C\t.\tPASS\t.\tGT:AD:PL\t1|2:0,6,7:90,80,70,60,0,50\t0/2:3,0,4:10,20,30,40,50,60\n",
      "method": "GET",
      "status": 200,
      "query": {},
      "headers": {}
    }
  ]
};
export const portalKind="";
