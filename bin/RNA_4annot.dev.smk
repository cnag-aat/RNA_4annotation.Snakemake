from datetime import datetime

date = datetime.now().strftime('%Y%m%d.%H%M%S')

module align_reads_workflow:
  snakefile: "../modules/align_reads.smk"
module get_models_workflow:
  snakefile: "../modules/get_models.dev.smk"

logs_dir = config["Parameters"]["logs_dir"]
if not os.path.exists(logs_dir):
    os.makedirs(logs_dir)

name = os.path.basename(config["Inputs"]["genome"])

# file wildcards 
illumina_files = config["Wildcards"]["illumina_fastqs"]
cDNA_files = config["Wildcards"]["cDNA_fastqs"]
dRNA_files = config["Wildcards"]["dRNA_fastqs"]
pb_files = config["Wildcards"]["isoseq_fastqs"]
pb_fasta_files = config["Wildcards"]["isoseq_fastas"]

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

if pb_files != None:
  outputs.append(expand(config["Isoseq"]["minimap_dir"] + "{pbfile}.sorted.bam", pbfile=pb_files.split(',')),)
  model_outputs.append(expand(config["Isoseq"]["minimap_dir"] + "{pbfile}" + ".stringtie.gtf", pbfile=pb_files.split(',')),)

if pb_fasta_files != None:
  outputs.append(expand(config["Isoseq"]["minimap_dir"] + "{pbfile}.sorted.bam", pbfile=pb_fasta_files.split(',')),)
  model_outputs.append(expand(config["Isoseq"]["minimap_dir"] + "{pbfile}" + ".stringtie.gtf", pbfile=pb_fasta_files.split(',')),)

bams_list = outputs

if config["Inputs"]["bams"] != None:
  ls = str(config["Inputs"]["bams"]) 
  a = ls.replace('[\'','')
  b= a.replace('\']', '')
  bams_list.append(b.split(','))

#print (bams_list)

rule all: 
  input:
    outputs,
    model_outputs,
    config["Outputs"]["models"],
   # config["Outputs"]["junctions"]
  log:
    logs_dir + str(date) + ".rule_all.out",
    logs_dir + str(date) + ".rule_all.err",

stringtie = ""
stringtie_in = {}
minimap_in = {}
extensions = {}
minimap_opts = {}
taco_bam_in = {}
taco_models_in = {}
taco_opts = {}

taco_models_in[config["Outputs"]["TACO_dir"]] = []
taco_opts[config["Outputs"]["TACO_dir"]] = config["Parameters"]["TACO_opts"]

if illumina_files != None:
  taco_bam_in[config["Illumina"]["star_dir"]] = []
  taco_models_in[config["Illumina"]["star_dir"]] = []
  taco_models_in[config["Outputs"]["TACO_dir"]].append(config["Illumina"]["star_dir"] + "TACO_output/TACO_assembled.gtf")
  for name in illumina_files.split(','):
    stringtie_in[config["Illumina"]["star_dir"] + name + ".stringtie.gtf"] = config["Illumina"]["stringtie_opts"]
    extensions[config["Illumina"]["star_dir"] + name] = "Aligned.sortedByCoord.out.bam"
    taco_bam_in[config["Illumina"]["star_dir"] ].append(config["Illumina"]["star_dir"] + name + "Aligned.sortedByCoord.out.bam")
    taco_models_in[config["Illumina"]["star_dir"] ].append(config["Illumina"]["star_dir"] + name + ".stringtie.gtf")
  taco_opts[config["Illumina"]["star_dir"]] = config["Illumina"]["TACO_opts"]
  if not os.path.exists(config["Illumina"]["star_dir"] + "/logs"):
    os.makedirs(config["Illumina"]["star_dir"] + "/logs")

  use rule star_index from align_reads_workflow with: 
    input:
      genome = config["Inputs"]["genome"]
    output: 
      genomedir =  directory(config["Illumina"]["genome_dir"]),
      ok = "index.ok"
    params:
      nbit= config["Illumina"]["nbit"]
    log:
      logs_dir + str(date) + ".star_index.j%j.out",
      logs_dir + str(date) + ".star_index.j%j.err",
    benchmark:
      logs_dir + str(date) + ".star_index.benchmark.txt"
    conda:
      "../envs/star2.7.10a.yaml"

  use rule star from align_reads_workflow with:
    input:
      genomedir = config["Illumina"]["genome_dir"],
      reads = lambda wildcards: [ config["Inputs"]["illumina_dir"] + "{file}.1.fastq.gz", config["Inputs"]["illumina_dir"] + "{file}.2.fastq.gz"] \
              if config["Illumina"]["PE"] == True else \
              lambda wildcards: [config["Inputs"]["illumina_dir"] + "{file}.fastq.gz"]              
    output:
      bam = config["Illumina"]["star_dir"] + "{file}" + "Aligned.sortedByCoord.out.bam",
    params:
      stardir = config["Illumina"]["star_dir"],
      basename = "{file}",
    log:
      logs_dir + str(date) + ".{file}.star.j%j.out",
      logs_dir + str(date) + ".{file}.star.j%j.err",
    benchmark:
      logs_dir + str(date) + ".{file}.star.benchmark.txt",
    conda:
      "../envs/star2.7.10a.yaml"
    threads:  config["Parameters"]["starCores"]

if cDNA_files != None:
  taco_bam_in[config["cDNA"]["minimap_dir"] ] = []
  taco_models_in[config["cDNA"]["minimap_dir"] ] = []
  taco_models_in[config["Outputs"]["TACO_dir"]].append(config["cDNA"]["minimap_dir"] + "TACO_output/TACO_assembled.gtf")
  for name in cDNA_files.split(','):
    stringtie_in[config["cDNA"]["minimap_dir"] + name + ".stringtie.gtf"] = config["cDNA"]["stringtie_opts"]
    extensions[config["cDNA"]["minimap_dir"] + name] = ".sorted.bam"
    minimap_in[config["cDNA"]["minimap_dir"] + name + ".sorted.bam"] = config["Inputs"]["cDNA_Dir"] + name + ".fastq.gz"
    minimap_opts[config["cDNA"]["minimap_dir"] + name] =  config["cDNA"]["minimap2_opts"]
    taco_bam_in[config["cDNA"]["minimap_dir"]].append(config["cDNA"]["minimap_dir"] + name + ".sorted.bam")
    taco_models_in[config["cDNA"]["minimap_dir"]].append(config["cDNA"]["minimap_dir"] + name + ".stringtie.gtf")
  taco_opts[config["cDNA"]["minimap_dir"] ] = config["cDNA"]["TACO_opts"],
  if not os.path.exists(config["cDNA"]["minimap_dir"] + "/logs"):
    os.makedirs(config["cDNA"]["minimap_dir"] + "/logs")

if dRNA_files != None:
  taco_bam_in[config["dRNA"]["minimap_dir"] ] = []
  taco_models_in[config["dRNA"]["minimap_dir"] ] = []
  taco_models_in[config["Outputs"]["TACO_dir"]].append(config["dRNA"]["minimap_dir"] + "TACO_output/TACO_assembled.gtf")
  for name in dRNA_files.split(','):
    stringtie_in[config["dRNA"]["minimap_dir"] + name + ".stringtie.gtf"] = config["dRNA"]["stringtie_opts"]
    extensions[config["dRNA"]["minimap_dir"] + name] = ".sorted.bam"
    minimap_in[config["dRNA"]["minimap_dir"] + name + ".sorted.bam"] = config["Inputs"]["dRNA_dir"] + name + ".fastq.gz"
    minimap_opts[config["dRNA"]["minimap_dir"] + name] =  config["dRNA"]["minimap2_opts"]
    taco_bam_in[config["dRNA"]["minimap_dir"]].append(config["dRNA"]["minimap_dir"] + name + ".sorted.bam")
    taco_models_in[config["dRNA"]["minimap_dir"]].append(config["dRNA"]["minimap_dir"] + name + ".stringtie.gtf")
  taco_opts[config["dRNA"]["minimap_dir"] ] = config["dRNA"]["TACO_opts"],
  if not os.path.exists(config["dRNA"]["minimap_dir"] + "/logs"):
    os.makedirs(config["dRNA"]["minimap_dir"] + "/logs")

if pb_fasta_files != None or pb_files != None:
  taco_bam_in[config["Isoseq"]["minimap_dir"]] = []
  taco_models_in[config["Isoseq"]["minimap_dir"]] = []
  taco_models_in[config["Outputs"]["TACO_dir"]].append(config["Isoseq"]["minimap_dir"] + "TACO_output/TACO_assembled.gtf")

if pb_files != None:
  for name in pb_files.split(','):
    stringtie_in[config["Isoseq"]["minimap_dir"] + name + ".stringtie.gtf"] = config["Isoseq"]["stringtie_opts"]
    extensions[config["Isoseq"]["minimap_dir"] + name] = ".sorted.bam"
    minimap_in[config["Isoseq"]["minimap_dir"] + name + ".sorted.bam"] = config["Inputs"]["PB_dir"] + name + ".fastq.gz"
    minimap_opts[config["Isoseq"]["minimap_dir"] + name] =  config["Isoseq"]["minimap2_opts"]
    taco_bam_in[config["Isoseq"]["minimap_dir"]].append(config["Isoseq"]["minimap_dir"] + name + ".sorted.bam")
    taco_models_in[config["Isoseq"]["minimap_dir"]].append(config["Isoseq"]["minimap_dir"] + name + ".stringtie.gtf")
  taco_opts[config["Isoseq"]["minimap_dir"]] = config["Isoseq"]["TACO_opts"],
  if not os.path.exists(config["Isoseq"]["minimap_dir"] + "/logs"):
      os.makedirs(config["Isoseq"]["minimap_dir"] + "/logs")

if pb_fasta_files != None:
  for name in pb_fasta_files.split(','):
    stringtie_in[config["Isoseq"]["minimap_dir"] + name + ".stringtie.gtf"] = config["Isoseq"]["stringtie_opts"]
    extensions[config["Isoseq"]["minimap_dir"] + name] = ".sorted.bam"
    minimap_in[config["Isoseq"]["minimap_dir"] + name + ".sorted.bam"] = config["Inputs"]["PB_dir"] + name + ".fasta"
    minimap_opts[config["Isoseq"]["minimap_dir"] + name] =  config["Isoseq"]["minimap2_opts"]
    taco_bam_in[config["Isoseq"]["minimap_dir"]].append(config["Isoseq"]["minimap_dir"] + name + ".sorted.bam")
    taco_models_in[config["Isoseq"]["minimap_dir"]].append(config["Isoseq"]["minimap_dir"] + name + ".stringtie.gtf")
  taco_opts[config["Isoseq"]["minimap_dir"]] = config["Isoseq"]["TACO_opts"],
  if not os.path.exists(config["Isoseq"]["minimap_dir"] + "/logs"):
    os.makedirs(config["Isoseq"]["minimap_dir"] + "/logs")

use rule minimap2 from align_reads_workflow with:
  input:
    genome = config["Inputs"]["genome"],
    reads = lambda wildcards: minimap_in[wildcards.dir + "/" + wildcards.file + ".sorted.bam"],
  output:
    bam = "{dir}/{file}.sorted.bam",
  params:
    basename = "{file}",
    minimap_opts = lambda wildcards:  minimap_opts[wildcards.dir + "/" + wildcards.file]
  log:
    "{dir}/logs" + str(date) + ".{file}.minimap2_pb.j%j.out",
    "{dir}/logs" + str(date) + ".{file}.minimap2_pb.j%j.err"
  benchmark:
    "{dir}/logs" + str(date) + ".{file}.minimap2_pb.benchmark.txt",
  conda:
    "../envs/minimap2.24.yaml"    
  threads: config["Parameters"]["minimapCores"]

use rule stringtie from get_models_workflow with:
  input:
    bam = lambda wildcards: wildcards.dir + "/" + wildcards.file + extensions[wildcards.dir + "/" + wildcards.file],
  output:
    models = "{dir}/{file}" + ".stringtie.gtf"
  params:
    basename = "{file}",
    stringtie_opts = lambda wildcards: stringtie_in[wildcards.dir + "/" +  wildcards.file + ".stringtie.gtf"],
  threads:
    config["Parameters"]["starCores"]
  log:
    "{dir}/logs/" + str(date) + ".{file}.stringtie.j%j.out",
    "{dir}/logs/" + str(date) + ".{file}.stringtie.j%j.err",
  benchmark:
    "{dir}/logs/" + str(date) + ".{file}.stringtie.benchmark.txt",
  conda:
    "../envs/stringtie2.2.1.yaml"

use rule join_models from get_models_workflow with:
  input:
    genome = config["Inputs"]["genome"],
   # bams = lambda wildcards: taco_bam_in[wildcards.dir],
    models = lambda wildcards: taco_models_in[wildcards.dir]
  output:
    joined_models = "{dir}TACO_output/TACO_assembled.gtf",
  params:
    TACO_opts = lambda wildcards: taco_opts[wildcards.dir],
    outdir = "{dir}"
  threads: 
    config["Parameters"]["TACOCores"]
  log:
    "{dir}logs/" + str(date) + ".TACO.j%j.out",
    "{dir}logs/" + str(date) + ".TACO.j%j.err",
  benchmark:
    "{dir}logs/" + str(date) + ".TACO.benchmark.txt",
  conda:
    "../envs/taco0.7.3.yaml"
  