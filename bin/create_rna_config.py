#!/usr/bin/env python3
import os
import json
import argparse
import sys
import re

#Author: Jessica Gomez, CNAG-CRG.
#Contact email: jessica.gomez@cnag.crg.eu
#Date:20191003

class CreateConfigurationFile(object):
    """Class which manages Configuration file Manager"""
      
    def __init__(self):
        """Class constructor"""
        #GENERAL PARAMETERS
        self.configFile = "RNAseq.config"             #Name of the json configuration file with the pipeline parameters to be created
        self.logs_dir = "logs"                        #Directory to keep all the log files        
        self.starCores = 4
        self.minimapCores = 4
        self.tacoCores = 4
        self.espressoCores = 4
        self.TACO_opts = "--isoform-frac 0 --filter-min-expr 0" 

        #INPUT PARAMETERS
        self.genome = None                            #Reference genome
        self.illumina_dir = None                      #Directory with the RNAseq reads
        self.cDNA_dir = None                          #Directory with the cDNA reads
        self.dRNA_dir = None                          #Directory with the dRNA reads
        self.pb_dir = None                            #Directory with the pacbio isoseq reads
        self.bams = None                              #bam files to get the junctions from them all, do not give this option if they are going to be generated by the pipeline
        self.gene_models = None	                      #gtf models that are going to be combined in the last step of this pipeline. Do not give this option if they are going to be generated by the pipeline
        self.samples_tsv = "samples.tsv"              #Samples tsv file to use with ESPRESSO
        self.glen = None                                #Path to the assembly.genome file
        
        #OUTPUT PARAMETERS
        self.models = "TACO_assembled.gtf"             #Final file with the TACO results
        self.TACO_dir = os.getcwd() + "/"                  #Directory to run TACO
        self.junctions = "alljunctions.final.gff3"     #Final file with the junctions

        #ILLUMINA
        self.genome_dir = "genome"                    #Directory for the genome index   
        self.PE = True
        self.nbit = 18                                #genomeChrBinNbits parameter of STAR
        self.star_dir = "star"                        #Directory for the mapping step index
        self.star_additional_options = ""             #Additional options to run star with
        self.stringtie_illumina_opts = ""             #Options to run stringtie in illumina mappings
        self.TACO_illumina_opts = ""                  #Options to run TACO in illumina models

        #TRIMGALORE PARAMETERS
        self.trim_galore_opts = "--gzip -q 20 --paired --retain_unpaired"
        self.Trim_Illumina_cores = 4                                                              #Number of threads to run the trim Illumina step

        #cDNA
        self.cdna_minimap_dir = "cDNA"                 #Directory for the cDNA mappimgs
        self.stringtie_cDNA_opts = "--conservative -R" #Options to run stringtie in cDNA mappings
        self.minimap2_cDNA_opts = ""
        self.TACO_cDNA_opts = "--isoform-frac 0.01"    #Options to run TACO in cDNA models
    
        #dRNA
        self.drna_minimap_dir = "dRNA"                                       #Directory for the dRNA mappimgs 
        self.stringtie_dRNA_opts = ""                                        #Options to run stringtie in dRNA mappings
        self.minimap2_dRNA_opts = " -uf -k14"
        self.TACO_dRNA_opts = "--isoform-frac 0.01 --filter-min-expr 0.2"    #Options to run TACO in dRNA models

        #PB
        self.pb_minimap_dir = "isoseq"                 #Directory for the pacbio isoseq mappimgs
        self.stringtie_pb_opts = "--conservative -R" #Options to run stringtie in pb isoseq mappings
        self.minimap2_pb_opts = ":hq -uf"
        self.TACO_pb_opts = "--isoform-frac 0.01"    #Options to run TACO in pb isoseq models

        #ESPRESSO
        self.espresso_path = "/software/assembly/src/ESPRESSO/espresso_v_1_3_0_beta/src/"
        self.espresso_outdir = "ESPRESSO_out/"
        
        #WILDCARDS
        self.illumina_fastqs = None                   #List with basename of the illumina fastqs
        self.cDNA_fastqs = None                       #List with basename of the cDNA fastqs  
        self.dRNA_fastqs = None                       #List with basename of the dRNA fastqs
        self.pb_fastqs = None                          #List with basename of the pacbio isoseq fastqs
        self.pb_fastas = None                          #List with basename of the pacbio isoseq fasta files                    

###
        #DICTIONARIES
        self.allParameters = {}
        self.generalParameters = {}
        self.inputParameters = {}
        self.outputParameters = {}
        self.trimgaloreParameters = {}
        self.illuminaParameters = {}
        self.cdnaParameters = {}
        self.drnaParameters = {}
        self.pbParameters = {}
        self.espressoParameters = {}
        self.wildcardParameters = {}

####

    def register_parameter(self, parser):
        """Register all parameters with the given
        argparse parser"""
        self.register_general(parser)
        self.register_input(parser)
        self.register_output(parser)
        self.register_trimgalore(parser)
        self.register_illumina(parser)
        self.register_cdna(parser)
        self.register_drna(parser)
        self.register_pb(parser)
        self.register_espresso(parser)
        self.register_wildcards(parser)

###

    def register_general(self, parser):
        """Register all general parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        general_group = parser.add_argument_group('General Parameters')
        general_group.add_argument('--configFile', dest="configFile", metavar="configFile", default=self.configFile, help='Configuration file with the pipeline parameters to be created. Default %s' % self.configFile)
        general_group.add_argument('--logs-dir', dest="logs_dir", metavar="logs_dir", help='Directory to keep all the log files. Default %s' % self.logs_dir)
        general_group.add_argument('--star-cpu', dest="starCores", metavar="starCores", type=int, default = self.starCores, help='Number of threads to run star. Default %s' % self.starCores)
        general_group.add_argument('--minimap-cpu', dest="minimapCores", metavar="minimapCores", type=int, default = self.minimapCores, help='Number of threads to run Minimap2. Default %s' % self.minimapCores)
        general_group.add_argument('--espresso-cpu', dest="espressoCores", metavar="espressoCores", type=int, default = self.espressoCores, help='Number of threads to run espresso. Default %s' % self.espressoCores)
        general_group.add_argument('--taco-cpu', dest="tacoCores", metavar="tacoCores", type=int, default = self.tacoCores, help='Number of threads to run TACO. Default %s' % self.tacoCores)
        general_group.add_argument('--TACO-all-opts', dest="TACO_opts", metavar="TACO_opts", default=self.TACO_opts, help='Options to run TACO when merging all the datasets. Default %s' % self.TACO_opts)

    def register_input(self, parser):
        """Register all input parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        input_group = parser.add_argument_group('Inputs')
        input_group.add_argument('--genome', dest="genome", metavar="genome", help='Path to the fasta genome. Default %s' % self.genome)
        input_group.add_argument('--glen', dest="glen", metavar="glen", help='Path to the assembly.genome file.')
        input_group.add_argument('--illumina-dir', dest="illumina_dir", help='Directory where the illumina fastqs are stored. Default %s' % self.illumina_dir)
        input_group.add_argument('--cdna-dir', dest="cDNA_dir", help='Directory where the cDNA fastqs are stored. Default %s' % self.cDNA_dir)
        input_group.add_argument('--drna-dir', dest="dRNA_dir", help='Directory where the dRNA fastqs are stored. Default %s' % self.dRNA_dir)
        input_group.add_argument('--pb-dir', dest="pb_dir", help='Directory where the pacbio isoseq fastqs are stored. Default %s' % self.pb_dir)
        input_group.add_argument('--bams', dest="bams", metavar="bams", nargs="+", help='bam files to get the junctions from them all, do not give this option if they are going to be generated by the pipeline. Default %s' % self.bams)
        input_group.add_argument('--gene_models', dest="gene_models", metavar="gene_models",  nargs="+", help='gtf models that are going to be combined in the last step of this pipeline. Do not give this option if they are going to be generated by the pipeline. Default %s' % self.gene_models) 
        input_group.add_argument('--samples-tsv', dest="samples_tsv", metavar="samples_tsv",  help='TSV file describing the samples to use with ESPRESSO. Default %s' % self.samples_tsv) 

    def register_output(self, parser):
        """Register all output parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        output_group = parser.add_argument_group('Outputs')
        output_group.add_argument('--gtf-models', dest="models", metavar="models", default=self.models, help='Path to the final stringtie gtf. Default %s' % self.models)
        output_group.add_argument('--TACO-dir', dest="TACO_dir", metavar="TACO_dir", default=self.TACO_dir, help='Directory to tun TACO. Default %s' % self.TACO_dir)
        output_group.add_argument('--junctions', dest="junctions", metavar="junctions", default=self.junctions, help='Path to the final junctions file. Default %s' % self.junctions)

    def register_illumina(self, parser):
        """Register all the illumina  parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        illumina_group = parser.add_argument_group('Illumina')
        illumina_group.add_argument('--genome-dir', dest="genome_dir", metavar="genome_dir", help='Directory for the genome index. Default %s' % self.genome_dir)
        illumina_group.add_argument('--star-dir', dest="star_dir", metavar="star_dir", help='Directory for the mapping step index. Default %s' % self.star_dir)
        illumina_group.add_argument('--genomeChrBinNbits', dest="nbit", metavar="nbit", type=int, default=self.nbit, help='genomeChrBinNbits parameter of STAR. Default %s' % self.nbit)
        illumina_group.add_argument('--star-opts', dest="star_additional_options", metavar="star_additional_options", help='Additional options to run star with. Default %s' % self.star_additional_options)
        illumina_group.add_argument('--no-pe', dest="PE", action="store_false", help='If specified, the input is not paired-end.')
        illumina_group.add_argument('--stringtie-illum-opts', dest="stringtie_illumina_opts", metavar="stringtie_illumina_opts", default=self.stringtie_illumina_opts, help='Options to run stringtie in illumina mappings. Default %s' % self.stringtie_illumina_opts)
        illumina_group.add_argument('--TACO-illum-opts', dest="TACO_illumina_opts", metavar="TACO_illumina_opts", default=self.TACO_illumina_opts, help='Options to run TACO in illumina mappings. Default %s' % self.TACO_illumina_opts)

    def register_trimgalore(self, parser):
        """Register all trimgalore parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        trimgalore_group = parser.add_argument_group('Trim_Galore')
        trimgalore_group.add_argument('--trim-galore-opts', dest="trim_galore_opts", metavar="trim_galore_opts", default=self.trim_galore_opts, help='Optional parameters for the rule trim_galore. Default %s' % self.trim_galore_opts)
        trimgalore_group.add_argument('--trim-Illumina-cores', type = int, dest="Trim_Illumina_cores", metavar="Trim_Illumina_cores", default=self.Trim_Illumina_cores, help='Number of threads to run the Illumina trimming step. Default %s' % self.Trim_Illumina_cores)

    def register_cdna(self, parser):
        """Register all the cdna parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        cdna_group = parser.add_argument_group('cDNA')
        cdna_group.add_argument('--cdna-mappings', dest="cdna_minimap_dir", metavar="cdna_minimap_dir", help='Directory for the cDNA Minimap2 mappings. Default %s' % self.cdna_minimap_dir)
        cdna_group.add_argument('--minimap2-cdna-opts', dest="minimap2_cDNA_opts", metavar="minimap2_cDNA_opts", default=self.minimap2_cDNA_opts, help='Options to run minimap2 in cDNA mappings. Default %s' % self.minimap2_cDNA_opts)
        cdna_group.add_argument('--stringtie-cdna-opts', dest="stringtie_cDNA_opts", metavar="stringtie_cDNA_opts", default=self.stringtie_cDNA_opts, help='Options to run stringtie in cDNA mappings. Default %s' % self.stringtie_cDNA_opts)
        cdna_group.add_argument('--TACO-cdna-opts', dest="TACO_cDNA_opts", metavar="TACO_cDNA_opts", default=self.TACO_cDNA_opts, help='Options to run TACO in cDNA mappings. Default %s' % self.TACO_cDNA_opts)

    def register_drna(self, parser):
        """Register all the dRNA parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        drna_group = parser.add_argument_group('dRNA')
        drna_group.add_argument('--drna-mappings', dest="drna_minimap_dir", metavar="drna_minimap_dir", help='Directory for the dRNA Minimap2 mappings. Default %s' % self.drna_minimap_dir)
        drna_group.add_argument('--minimap2-drna-opts', dest="minimap2_dRNA_opts", metavar="minimap2_dRNA_opts", default=self.minimap2_dRNA_opts, help='Options to run minimap2 in dRNA mappings. Default %s' % self.minimap2_dRNA_opts)
        drna_group.add_argument('--stringtie-drna-opts', dest="stringtie_dRNA_opts", metavar="stringtie_dRNA_opts", default=self.stringtie_dRNA_opts, help='Options to run stringtie in dRNA mappings. Default %s' % self.stringtie_dRNA_opts)
        drna_group.add_argument('--TACO-drna-opts', dest="TACO_dRNA_opts", metavar="TACO_dRNA_opts", default=self.TACO_dRNA_opts, help='Options to run TACO in dRNA mappings. Default %s' % self.TACO_dRNA_opts)

    def register_pb(self, parser):
        """Register all the pacbio parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        pb_group = parser.add_argument_group('Isoseq')
        pb_group.add_argument('--pb-mappings', dest="pb_minimap_dir", metavar="pb_minimap_dir", help='Directory for the Pacbio Isoseq Minimap2 mappings. Default %s' % self.pb_minimap_dir)
        pb_group.add_argument('--minimap2-pb-opts', dest="minimap2_pb_opts", metavar="minimap2_pb_opts", default=self.minimap2_pb_opts, help='Options to run minimap2 in Pacbio Isoseq mappings. Default %s' % self.minimap2_pb_opts)
        pb_group.add_argument('--stringtie-pb-opts', dest="stringtie_pb_opts", metavar="stringtie_pb_opts", default=self.stringtie_pb_opts, help='Options to run stringtie in Pacbio Isoseq mappings. Default %s' % self.stringtie_pb_opts)
        pb_group.add_argument('--TACO-pb-opts', dest="TACO_pb_opts", metavar="TACO_pb_opts", default=self.TACO_pb_opts, help='Options to run TACO in isoseq mappings. Default %s' % self.TACO_pb_opts)


    def register_espresso(self, parser):
        """Register all the espresso parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        espresso_group = parser.add_argument_group('ESPRESSO')
        espresso_group.add_argument('--espresso-path', dest="espresso_path", metavar="espresso_path", help='Path to the ESPRESSO scripts. Default %s' % self.espresso_path)
        espresso_group.add_argument('--espresso-outdir', dest="espresso_outdir", metavar="espresso_outdir", default=self.espresso_outdir, help='Directory for runnning ESPRESSO. Default %s' % self.espresso_outdir)

    def register_wildcards(self, parser):
        """Register all wildcards parameters with the given
        argparse parser

        parser -- the argparse parser
        """
        wildcards_group = parser.add_argument_group('Wildcards')
        wildcards_group.add_argument('--illumina-reads', dest="illumina_fastqs", metavar="illumina_fastqs", help='List with basename of the illumina fastqs. Default %s' % self.illumina_fastqs)
        wildcards_group.add_argument('--cdna-reads', dest="cDNA_fastqs", metavar="cDNA_fastqs", help='List with basename of the cDNA fastqs. Default %s' % self.cDNA_fastqs)
        wildcards_group.add_argument('--drna-reads', dest="dRNA_fastqs", metavar="dRNA_fastqs", help='List with basename of the dRNA fastqs. Default %s' % self.dRNA_fastqs)
        wildcards_group.add_argument('--pb-reads', dest="pb_fastqs", metavar="pb_fastqs", help='List with basename of the Pacbio Isoseq fastqs. Default %s' % self.pb_fastqs)
        wildcards_group.add_argument('--pb-reads-fasta', dest="pb_fastas", metavar="pb_fastas", help='List with basename of the Pacbio Isoseq fastas. Default %s' % self.pb_fastas)

####

    def check_parameters(self,args):
        """Check parameters consistency
            
        args -- set of parsed arguments"""

        working_dir = os.getcwd() + "/"

        if args.genome == None:
            parser.print_help()
            print ("Sorry! No genome fasta file defined")
            sys.exit(-1)
        else:
            args.genome = os.path.abspath(args.genome) 
            if not os.path.exists(args.genome):
                print (args.genome + " not found")

        if args.illumina_dir != None:
            args.illumina_dir = os.path.abspath(args.illumina_dir) + "/" 

        if args.cDNA_dir != None:
            args.cDNA_dir = os.path.abspath(args.cDNA_dir) + "/" 

        if args.dRNA_dir != None:
            args.dRNA_dir = os.path.abspath(args.dRNA_dir) + "/" 

        if args.pb_dir != None:
            args.pb_dir = os.path.abspath(args.pb_dir) + "/" 
           
        if args.logs_dir:
            args.logs_dir = os.path.abspath(args.logs_dir) + "/"
        else:
            args.logs_dir = working_dir + self.logs_dir + "/"

        if args.genome_dir:
            args.genome_dir = os.path.abspath(args.genome_dir) + "/"
        else:
            args.genome_dir = working_dir + self.genome_dir + "/"

        if args.star_dir:
            args.star_dir = os.path.abspath(args.star_dir) + "/"
        else:
            args.star_dir = working_dir + self.star_dir + "/"

        if args.cdna_minimap_dir:
            args.cdna_minimap_dir = os.path.abspath(args.cdna_minimap_dir) + "/"
        else:
            args.cdna_minimap_dir = working_dir + self.cdna_minimap_dir + "/"

        if args.drna_minimap_dir:
            args.drna_minimap_dir = os.path.abspath(args.drna_minimap_dir) + "/"
        else:
            args.drna_minimap_dir = working_dir + self.drna_minimap_dir + "/"

        if args.pb_minimap_dir:
            args.pb_minimap_dir = os.path.abspath(args.pb_minimap_dir) + "/"
        else:
            args.pb_minimap_dir = working_dir + self.pb_minimap_dir + "/"

        if args.models:
            args.models = os.path.abspath(args.models) 
        else:
            args.models = working_dir + self.models 

        if args.junctions:
            args.junctions = os.path.abspath(args.junctions) 
        else:
            args.junctions = working_dir + self.junctions 

        if args.espresso_path == None:
            args.espresso_path = self.espresso_path
        else:
            args.espresso_path = os.path.abspath(args.espresso_path)

        if args.espresso_outdir == None:
            args.espresso_outdir = working_dir + self.espresso_path
        else:
            args.espresso_outdir = os.path.abspath(args.espresso_outdir)

        if args.glen:
          if not os.path.exists(args.glen):
            print (args.glen) + " not found"
          else:
            args.glen = os.path.abspath(args.glen)
        ##Assign wildcards
        samples=[]
        if args.illumina_fastqs == None and args.illumina_dir != None:
            for r, d, f in os.walk(args.illumina_dir):
                for file in f:
                    if args.PE == True and re.search('.1.fastq.gz', file):
                        a = file.replace('.1.fastq.gz','')
                        if args.illumina_fastqs == None:
                            args.illumina_fastqs = a
                        else:
                            args.illumina_fastqs += "," + a
                        samples.append(args.star_dir + a + "Aligned.sortedByCoord.out.sam\t" + a)
                    elif args.PE == False and re.search('fastq.gz', file):
                        a = file.replace('.fastq.gz','')
                        if args.illumina_fastqs == None:
                            args.illumina_fastqs = a
                        else:
                            args.illumina_fastqs += "," + a
                        samples.append(args.star_dir + a + "Aligned.sortedByCoord.out.sam\t" + a)

        if args.cDNA_fastqs == None and args.cDNA_dir != None:
            for r, d, f in os.walk(args.cDNA_dir):
                for file in f:
                    if re.search('fastq.gz', file):
                        a = file.replace('.fastq.gz','')
                        if args.cDNA_fastqs == None:
                            args.cDNA_fastqs = a
                        else:
                            args.cDNA_fastqs += "," + a
                        samples.append(args.cdna_minimap_dir + a + ".sam\t" + a)

        if args.dRNA_fastqs == None and args.dRNA_dir != None:
            for r, d, f in os.walk(args.dRNA_dir):
                for file in f:
                    if re.search('fastq.gz', file):
                        a = file.replace('.fastq.gz','')
                        if args.dRNA_fastqs == None:
                            args.dRNA_fastqs = a
                        else:
                            args.dRNA_fastqs += "," + a
                        samples.append(args.drna_minimap_dir + a + ".sam\t" + a)

        if args.pb_fastqs == None and args.pb_dir != None:
            for r, d, f in os.walk(args.pb_dir):
                for file in f:
                    if re.search('fastq.gz', file):
                        a = file.replace('.fastq.gz','')
                        if args.pb_fastqs == None:
                            args.pb_fastqs = a
                        else:
                            args.pb_fastqs += "," + a
                        samples.append(args.pb_minimap_dir + a + ".sam\t" + a)

        if args.pb_fastas == None and args.pb_dir != None:
            for r, d, f in os.walk(args.pb_dir):
                for file in f:
                    if re.search('fasta', file):
                        a = file.replace('.fasta','')
                        if args.pb_fastas == None:
                            args.pb_fastas = a
                        else:
                            args.pb_fastas += "," + a
                        samples.append(args.pb_minimap_dir + a + ".sam\t" + a)
        
        if args.samples_tsv == None:
            args.samples_tsv = working_dir + self.samples_tsv
        else:
            args.samples_tsv = os.path.abspath(args.samples_tsv)
        if not os.path.exists(args.samples_tsv):
            f = open(args.samples_tsv, "w")
            for line in samples:
                f.write(line + "\n")
            f.close()


###

    def storeGeneralParameters(self,args):
        """Updates general parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.generalParameters["logs_dir"] = args.logs_dir
        self.generalParameters["starCores"] = args.starCores
        self.generalParameters["minimapCores"] = args.minimapCores
        self.generalParameters["TACOCores"] = args.tacoCores
        self.generalParameters["EspressoCores"] = args.espressoCores
        self.generalParameters["TACO_opts"] = args.TACO_opts
        self.allParameters  ["Parameters"] = self.generalParameters

    def storeInputParameters(self,args):
        """Updates input parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.inputParameters["genome"] = args.genome
        self.inputParameters["genome_lengths"] = args.glen
        self.inputParameters["illumina_dir"] = args.illumina_dir
        self.inputParameters["cDNA_dir"] = args.cDNA_dir
        self.inputParameters["dRNA_dir"] = args.dRNA_dir
        self.inputParameters["PB_dir"] = args.pb_dir
        self.inputParameters["bams"] = args.bams
        self.inputParameters["gene_models"] = args.gene_models
        self.inputParameters["Samples_tsv"] = args.samples_tsv
        self.allParameters ["Inputs"] = self.inputParameters

    def storeOutputParameters(self,args):
        """Updates input parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.outputParameters["models"] = args.models
        self.outputParameters["junctions"] = args.junctions
        self.outputParameters["TACO_dir"] = args.TACO_dir
        self.allParameters ["Outputs"] = self.outputParameters

    def storeTrimgaloreParameters(self,args):
        """Updates the Trim_Galore parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.trimgaloreParameters["options"] = args.trim_galore_opts
        self.trimgaloreParameters["Trim_Illumina_cores"] = args.Trim_Illumina_cores
        self.allParameters ["Trim_Galore"] = self.trimgaloreParameters

    def storeIlluminaParameters(self,args):
        """Updates illumina parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.illuminaParameters["genome_dir"] = args.genome_dir
        self.illuminaParameters["PE"] = args.PE
        self.illuminaParameters["nbit"] = args.nbit
        self.illuminaParameters["additional_opts"] = args.star_additional_options
        self.illuminaParameters["star_dir"] = args.star_dir
        self.illuminaParameters["stringtie_opts"] = args.stringtie_illumina_opts
        self.illuminaParameters["TACO_opts"] = args.TACO_illumina_opts
        self.allParameters ["Illumina"] = self.illuminaParameters

    def storeCdnaParameters(self,args):
        """Updates cDNA parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.cdnaParameters["minimap_dir"] = args.cdna_minimap_dir
        self.cdnaParameters["minimap2_opts"] = args.minimap2_cDNA_opts
        self.cdnaParameters["stringtie_opts"] = args.stringtie_cDNA_opts
        self.cdnaParameters["TACO_opts"] = args.TACO_cDNA_opts
        self.allParameters ["cDNA"] = self.cdnaParameters

    def storeDrnaParameters(self,args):
        """Updates dRNA parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.drnaParameters["minimap_dir"] = args.drna_minimap_dir
        self.drnaParameters["minimap2_opts"] = args.minimap2_dRNA_opts
        self.drnaParameters["stringtie_opts"] = args.stringtie_dRNA_opts
        self.drnaParameters["TACO_opts"] = args.TACO_dRNA_opts
        self.allParameters ["dRNA"] = self.drnaParameters

    def storePbParameters(self,args):
        """Updates pacbio isoseq parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.pbParameters["minimap_dir"] = args.pb_minimap_dir
        self.pbParameters["minimap2_opts"] = args.minimap2_pb_opts
        self.pbParameters["stringtie_opts"] = args.stringtie_pb_opts
        self.pbParameters["TACO_opts"] = args.TACO_pb_opts
        self.allParameters ["Isoseq"] = self.pbParameters

    def storeEspressoParameters(self,args):
        """Updates ESPRESSO parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.espressoParameters["ESPRESSO_path"] = args.espresso_path
        self.espressoParameters["ESPRESSO outdir"] = args.espresso_outdir
        self.allParameters ["ESPRESSO"] = self.espressoParameters

    def storeWildcardParameters(self,args):
        """Updates wildcard parameters to the map of parameters to be store in a JSON file

        args -- set of parsed arguments
        """
        self.wildcardParameters["illumina_fastqs"] = args.illumina_fastqs
        self.wildcardParameters["cDNA_fastqs"] = args.cDNA_fastqs
        self.wildcardParameters["dRNA_fastqs"] = args.dRNA_fastqs
        self.wildcardParameters["isoseq_fastqs"] = args.pb_fastqs
        self.wildcardParameters["isoseq_fastas"] = args.pb_fastas
        self.allParameters ["Wildcards"] = self.wildcardParameters

#####

#1.Create object class Configuration File
configManager = CreateConfigurationFile()

#2.Create object for argument parsinng
parser = argparse.ArgumentParser(prog="create_configuration_file",
                description="Create a configuration json file for the repeat annotation pipeline."
                )     

#2.1 Updates arguments and parsing
configManager.register_parameter(parser)

args = parser.parse_args()

#2.2 Check Parameters
configManager.check_parameters(args)

#3. store arguments to super map structure
configManager.storeGeneralParameters(args)
configManager.storeInputParameters(args)
configManager.storeOutputParameters(args)
configManager.storeTrimgaloreParameters(args)
configManager.storeIlluminaParameters(args)
configManager.storeCdnaParameters(args)
configManager.storeDrnaParameters(args)
configManager.storePbParameters(args)
configManager.storeEspressoParameters(args)
configManager.storeWildcardParameters(args)

###

#4. Store JSON file
with open(args.configFile, 'w') as of:
    json.dump(configManager.allParameters, of, indent=2)


