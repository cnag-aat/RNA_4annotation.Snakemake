from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

rule star_index:
    input:
      genome = "genome.fa"
    output: 
      genomedir = directory("star_index"),
      ok = "index.ok"
    params:
      nbit= 18
    conda:
      "../envs/star2.7.10a.yaml"
    shell:
      "mkdir -p {output.genomedir};"
      "STAR --runMode genomeGenerate --genomeDir {output.genomedir} --genomeFastaFiles {input.genome} --genomeChrBinNbits {params.nbit};"
      "touch {output.ok};"

rule star:
  input:
    genomedir ="star_index",
    reads = [ "reads.1.fastq.gz", "reads.2.fastq.gz" ]
  output:
    bam = "star/readsAligned.sortedByCoord.out.bam",
  params:
    stardir = "star",
    basename = "reads",
  conda:
    "../envs/star2.7.10a.yaml"
  threads: 4
  shell:
    "mkdir -p {params.stardir};"
    "cd {params.stardir};"
    "STAR --genomeDir {input.genomedir} --readFilesIn {input.reads} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.basename} --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical;"
    "echo 'STAR run completed.';"

rule minimap2:
  input:
    genome = "genome.fa",
    reads = "cdna_reads.fastq.gz"
  output:
    bam =  "cDNA/cdna_reads.sorted.bam",
  params:
    basename = "cdna_reads",
    minimap_opts  = ""
  conda:
    "../envs/minimap2.24.yaml"
  threads: 8
  shell:
    "minimap2 -x splice{params.minimap_opts} -t {threads} -a {input.genome} {input.reads} |  samtools sort -@ {threads}  -O BAM -o {output.bam};"
    "samtools index {output.bam};"
    "echo 'MINIMAP2 run completed.';"