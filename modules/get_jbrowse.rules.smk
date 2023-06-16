from datetime import datetime
import re
import os

date = datetime.now().strftime('%Y%m%d.%H%M%S')

rule get_bw:
  input:
    bam = "RNAseq.bam",
    glen = "assembly.genome"
  output:
    bg = "RNAseq.bg",
    bw = "RNAseq.bw"
  params:
    opt = " -split "
  conda:
    "../envs/bedtools2.30.0.yaml"
  shell:
    "genomeCoverageBed {params.opt} -ibam {input.bam} -bg | sortBed -i -  > {output.bg};"
    "bedGraphToBigWig {output.bg} {input.glen} {output.bw};"
