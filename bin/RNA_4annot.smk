shell.prefix("source ~jgomez/init_shell.sh;")

from datetime import datetime

date = datetime.now().strftime('%Y%m%d.%H%M%S')

logs_dir = config["Parameters"]["logs_dir"]
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

name = os.path.basename(config["Inputs"]["genome"])

# file wildcards 
illumina_files = config["Wildcards"]["illumina_fastqs"]
cDNA_files = config["Wildcards"]["cDNA_fastqs"]
dRNA_files = config["Wildcards"]["dRNA_fastqs"]

outputs = []
model_outputs = []
if illumina_files != None:
  outputs.append(expand(config["Illumina"]["star_dir"] + "{file}" + "Aligned.sortedByCoord.out.bam", file=illumina_files.split(',')),)
  model_outputs.append(expand(config["Illumina"]["star_dir"] + "{file}" + ".stringtie.gtf", file=illumina_files.split(',')),)

if cDNA_files != None:
  outputs.append(expand(config["cDNA"]["minimap_dir"] + "{cdnafile}.sorted.bam", cdnafile=cDNA_files.split(',')),)
  model_outputs.append(expand(config["cDNA"]["minimap_dir"] + "{cdnafile}" + ".stringtie.gtf", cdnafile=cDNA_files.split(',')),)

if dRNA_files != None:
  outputs.append(expand(config["dRNA"]["minimap_dir"] + "{drnafile}.sorted.bam", drnafile=dRNA_files.split(',')),)
  model_outputs.append(expand(config["dRNA"]["minimap_dir"] + "{drnafile}" + ".stringtie.gtf", drnafile=dRNA_files.split(',')),)

bams_list = outputs

if config["Inputs"]["bams"] != None:
  ls = str(config["Inputs"]["bams"]) 
  a = ls.replace('[\'','')
  b= a.replace('\']', '')
  bams_list.append(b.split(','))

#print (bams_list)

rule all: ##modify rule all when everything is ready!!
  input:
    outputs,
    model_outputs,
    config["Outputs"]["models"],
    config["Outputs"]["junctions"]
  log:
    logs_dir + str(date) + ".rule_all.out",
    logs_dir + str(date) + ".rule_all.err",

stringtie = config["Parameters"]["stringtiePath"]

if illumina_files != None:
  if config["Illumina"]["stringtie_opts"] == None:
    stringtie_opts = "";
  else:
    stringtie_opts = config["Illumina"]["stringtie_opts"]

  rule index_genome:
    input:
      genome = config["Inputs"]["genome"]
    output: 
      genomedir = directory(config["Illumina"]["genome_dir"]),
      ok = "index.ok"
    params:
      nbit= config["Illumina"]["nbit"]
    log:
      logs_dir + str(date) + ".star_index.out",
      logs_dir + str(date) + ".star_index.err",
    shell:
      "mkdir -p {output.genomedir};"
      "module purge;module load gcc/6.3.0;module load STAR/2.7.2a;"
      "STAR --runMode genomeGenerate --genomeDir {output.genomedir} --genomeFastaFiles {input.genome} --genomeChrBinNbits {params.nbit};"
      "touch {output.ok};"

  if config["Illumina"]["PE"] == True:
    rule star:
      input:
        genomedir = rules.index_genome.output.genomedir,
        read1 = config["Inputs"]["illumina_dir"] + "{file}.1.fastq.gz",
        read2 = config["Inputs"]["illumina_dir"] + "{file}.2.fastq.gz" 
      output:
        bam= config["Illumina"]["star_dir"] + "{file}" + "Aligned.sortedByCoord.out.bam",
        models = config["Illumina"]["star_dir"] + "{file}" + ".stringtie.gtf"
      params:
        stardir = config["Illumina"]["star_dir"],
        basename = "{file}",
        stringtie_opts = stringtie_opts
      threads:
        config["Parameters"]["starCores"]
      log:
        logs_dir + str(date) + ".{file}.star.out",
        logs_dir + str(date) + ".{file}.star.err",
      shell:
        "mkdir -p {params.stardir};"
        "module purge;module load gcc/6.3.0;module load STAR/2.7.2a;"
        "cd {params.stardir};"
        "STAR --genomeDir {input.genomedir} --readFilesIn {input.read1} {input.read2} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.basename} --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical;"
        "echo 'STAR run completed, running Stringtie now...';"
        "{stringtie} {output.bam} {params.stringtie_opts} -l {params.basename} -p {threads} -o {output.models};"
  else:
    rule star:
      input:
        genomedir = rules.index_genome.output.genomedir,
        reads = config["Inputs"]["illumina_dir"] + "{file}.fastq.gz"
      output:
        bam= config["Illumina"]["star_dir"] + "{file}" + "Aligned.sortedByCoord.out.bam",
        models = config["Illumina"]["star_dir"] + "{file}" + ".stringtie.gtf"
      params:
        stardir = config["Illumina"]["star_dir"],
        basename = "{file}",
        stringtie_opts = stringtie_opts
      threads:
        config["Parameters"]["starCores"]
      log:
        logs_dir + str(date) + ".{file}.star.out",
        logs_dir + str(date) + ".{file}.star.err",
      shell:
        "mkdir -p {params.stardir};"
        "module purge;module load gcc/6.3.0;module load STAR/2.7.2a;"
        "cd {params.stardir};"
        "STAR --genomeDir {input.genomedir} --readFilesIn {input.reads} --readFilesCommand zcat --runThreadN {threads} --outFileNamePrefix {params.basename} --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outSAMattrIHstart 0 --outFilterIntronMotifs RemoveNoncanonical;"
        "echo 'STAR run completed, running Stringtie now...';"
        "{stringtie} {output.bam} {params.stringtie_opts} -l {params.basename} -p {threads} -o {output.models};"

if cDNA_files != None:
  if config["cDNA"]["stringtie_opts"] == None:
    stringtie_opts = "";
  else:
    stringtie_opts = config["cDNA"]["stringtie_opts"]

  rule align_cDNA:
    input:
      genome = config["Inputs"]["genome"],
      reads = config["Inputs"]["cDNA_dir"] + "{cdnafile}.fastq.gz"
    output:
      bam =  config["cDNA"]["minimap_dir"] + "{cdnafile}.sorted.bam",
      models = config["cDNA"]["minimap_dir"] + "{cdnafile}.stringtie.gtf"
    params:
      basename = "{cdnafile}",
      stringtie_opts = stringtie_opts
    log:
        logs_dir + str(date) + ".{cdnafile}.minimap2_cDNA.out",
        logs_dir + str(date) + ".{cdnafile}.minimap2_cDNA.err",
    threads: config["Parameters"]["minimapCores"]
    shell:
      "module purge;"
      "module load gcc/6.3.0 zlib;module load MINIMAP2/2.14; module load samtools;"
      "minimap2 -x splice -t {threads} -a {input.genome} {input.reads} |  samtools sort -@ {threads} -T $TMPDIR -O BAM -o {output.bam};"
      "samtools index {output.bam};"
      "echo 'MINIMAP2 run completed, running Stringtie now...';"
      "{stringtie} -o {output.models}  -l {params.basename} -L {params.stringtie_opts} -A {params.basename}.abundance.txt {output.bam};"

if dRNA_files != None:
  if config["dRNA"]["stringtie_opts"] == None:
    stringtie_opts = "";
  else:
    stringtie_opts = config["dRNA"]["stringtie_opts"]

  rule align_dRNA:
    input:
      genome = config["Inputs"]["genome"],
      reads = config["Inputs"]["dRNA_dir"] + "{drnafile}.fastq.gz"
    output:
      bam = config["dRNA"]["minimap_dir"] + "{drnafile}.sorted.bam",
      models = config["dRNA"]["minimap_dir"] + "{drnafile}.stringtie.gtf",
    params:
      basename = "{drnafile}",
      stringtie_opts = stringtie_opts
    log:
        logs_dir + str(date) + ".{drnafile}.minimap2_dRNA.out",
        logs_dir + str(date) + ".{drnafile}.minimap2_dRNA.err"
    threads: config["Parameters"]["minimapCores"]
    shell:
      "module purge;"
      "module load gcc/6.3.0 zlib;module load MINIMAP2/2.14; module load samtools;"
      "minimap2 -x splice -uf -k14 -t {threads} -a {input.genome} {input.reads} |  samtools sort -@ {threads} -T $TMPDIR -O BAM -o {output.bam};"
      "samtools index {output.bam};"
      "echo 'MINIMAP2 run completed, running Stringtie now...';"
      "{stringtie} -o {output.models}  -l {params.basename} -L {params.stringtie_opts} -A {params.basename}.abundance.txt {output.bam};"

rule join_RNAseq:
  input:
    genome = config["Inputs"]["genome"],
    bams = bams_list,
    models = model_outputs
  output:
    joined_models = config["Outputs"]["models"],
    junctions = config["Outputs"]["junctions"]
  params:
    illumina_opts = config["Illumina"]["TACO_opts"],
    cDNA_opts = config["cDNA"]["TACO_opts"],
    dRNA_opts = config["dRNA"]["TACO_opts"],
    TACO_opts = config["Parameters"]["TACO_opts"],
    outdir = config["Outputs"]["TACO_dir"]
  log:
    logs_dir + str(date) + ".joinRNAseq.out",
    logs_dir + str(date) + ".joinRNAseq.err",
  threads: config["Parameters"]["TACOCores"]
  run:
    if illumina_files != None and not os.path.exists("TACO_illumina/assembly.FPKM.gtf"):
      models = expand(config["Illumina"]["star_dir"] + "{file}" + ".stringtie.gtf", file=illumina_files.split(','))
      shell(
        "module purge; module load PYTHON/2.7.5;"
        "echo {models}| sed 's/\s/\\n/g' > illumina_gtf_files.txt;"
        "/apps/TACO/SRC/taco-0.6.3/build/scripts-2.7/taco_run.py -o TACO_illumina -p {threads} --filter-splice-juncs {params.illumina_opts} --ref-genome-fasta {input.genome} illumina_gtf_files.txt;"
        "cat TACO_illumina/assembly.gtf | sed 's/expr/FPKM/g' > TACO_illumina/assembly.FPKM.gtf;"
      )
    if cDNA_files != None:
      models = expand(config["cDNA"]["minimap_dir"] + "{cdnafile}" + ".stringtie.gtf", cdnafile=cDNA_files.split(','))
      shell (
        "module purge; module load PYTHON/2.7.5;"
        "echo {models}| sed 's/\s/\\n/g' > cDNA_gtf_files.txt;"
        "/apps/TACO/SRC/taco-0.6.3/build/scripts-2.7/taco_run.py -o TACO_cDNA -p {threads} --filter-splice-juncs {params.cDNA_opts} --ref-genome-fasta {input.genome} cDNA_gtf_files.txt;"   
        "cat TACO_cDNA/assembly.gtf | sed 's/expr/FPKM/g' > TACO_cDNA/assembly.FPKM.gtf;"
        
      )
    if dRNA_files != None:
      models = expand(config["dRNA"]["minimap_dir"] + "{drnafile}" + ".stringtie.gtf", drnafile=dRNA_files.split(','))
      shell (
        "module purge; module load PYTHON/2.7.5;"
        "echo {models}| sed 's/\s/\\n/g' > dRNA_gtf_files.txt;"
        "/apps/TACO/SRC/taco-0.6.3/build/scripts-2.7/taco_run.py -o TACO_dRNA -p {threads} --filter-splice-juncs {params.dRNA_opts} --ref-genome-fasta {input.genome} dRNA_gtf_files.txt;"   
        "cat TACO_dRNA/assembly.gtf | sed 's/expr/FPKM/g' > TACO_dRNA/assembly.FPKM.gtf;"     
      )
    models = []
    r = 0
    if os.path.exists("TACO_illumina"):
      models.append("TACO_illumina/assembly.FPKM.gtf")
      r+=1
    if os.path.exists("TACO_cDNA"):
      models.append("TACO_cDNA/assembly.FPKM.gtf")
      r+=1
    if os.path.exists("TACO_dRNA"):
      models.append("TACO_dRNA/assembly.FPKM.gtf")
      r+=1
    if r > 1:
      shell(
        "module purge; module load PYTHON/2.7.5;"
        "echo {models}| sed 's/\s/\\n/g' > join_gtf_files.txt;"
        "/apps/TACO/SRC/taco-0.6.3/build/scripts-2.7/taco_run.py -o {params.outdir} -p {threads} --filter-splice-juncs {params.TACO_opts} --ref-genome-fasta {input.genome} join_gtf_files.txt;" 
        "ln -s {params.outdir}/assembly.gtf {output.joined_models};"
      )
    else:
      shell(
        "ln -s {models} {output.joined_models};"
      )
    shell(
      "module purge;conda activate /scratch/project/devel/aateam/src/portcullis;"
      "portcullis full --intron_gff -t {threads} {input.genome} {input.bams};"
      "gawk \'$7!=\"?\"\' portcullis_out/3-filt/portcullis_filtered.pass.junctions.intron.gff3 > {output.junctions};"
      "conda deactivate;"
    )
