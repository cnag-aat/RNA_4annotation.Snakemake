rule stringtie:
  input:
    bam = "star/readsAligned.sortedByCoord.out.bam"
  output: 
    models = "stringtie/reads.stringtie.gtf"
  params:
    stringtie_opts = "",
    basename = "reads",
  threads: 4
  conda:
    "../envs/stringtie2.2.1.yaml"
  shell:
    "stringtie {input.bam} {params.stringtie_opts} -l {params.basename} -p {threads} -o {output.models};"

rule join_models:
  input:
    genome = "assemvbly.fa",
    #bams =  "star/readsAligned.sortedByCoord.out.bam" ,
    models = "stringtie/reads.stringtie.gtf" 
  output:
    joined_models = "TACO_assembled.gtf",
  params:
    TACO_opts = "",
    outdir = "TACO_output/"
  threads: 4
  conda:
    "../envs/taco0.7.3.yaml"
  shell:
    "cd {params.outdir};"
    "echo {input.models}| sed 's/\s/\\n/g' > gtf_files.txt;"
    "taco_run -o {params.outdir}workdir -p {threads} --filter-splice-juncs {params.TACO_opts} --ref-genome-fasta {input.genome} gtf_files.txt;"
    "cat {params.outdir}/assembly.gtf | sed 's/expr/FPKM/g' > {params.outdir}/assembly.FPKM.gtf;"