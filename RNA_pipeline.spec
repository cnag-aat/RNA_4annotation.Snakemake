{
  "all":{
    "name": "{rule}.{name}",
    "time": "0:05:00",
    "queue": "genB,main"
  },
  "index_genome":{
    "name": "{rule}.{name}",
    "time": "6:00:00",
    "threads": 12,
    "queue": "genB,main"
  },
  "star":{
    "name": "{rule}.{name}.{wildcards.file}",
    "time": "12:00:00",
    "queue": "genB,main"
  },
  "align_cDNA":{
    "name": "{rule}.{name}.{wildcards.cdnafile}",
    "time": "24:00:00",
    "queue": "genB,main"
  },
  "align_dRNA":{
    "name": "{rule}.{name}.{wildcards.drnafile}",
    "time": "5:00:00",
    "queue": "genB,main"
  },
  "join_RNAseq":{
    "name": "{rule}.{name}",
    "time": "22:00:00",
    "threads": 8,
    "queue": "genB,main"
  }
}
