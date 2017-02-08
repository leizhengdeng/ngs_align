#!/usr/bin/env python
import os
import sys
import argparse
import logging
python_lib = os.path.join(os.getenv("HOME"), "python_lib")
if not python_lib in sys.path: sys.path.append(python_lib)
from fs import FS
from align_factory import ALIGN_FACTORY

class NGS_ALIGN:
	"""
	NGS_ALIGN: scan the input directory containing fq.gz or fastq.gz files to create fq1 fq2 mapping (fq2 is "" if SE)
	for each (fq1, fq2) pair, call ALIGN_FACTORY to create instance of DNA_ALIGN, RNA_ALIGN (base class ALIGN) to do sequencing-type specific mapping
	print the command lines to the file: [prj_dir]/alignment/[sample_name]/pipeline_commands.sh
	create pipeline running shell file:  [prj_dir]/pipeline_run.sh, as well as log file [prj_dir]/ngs_align.log 
	"""
	def __init__(self, args):
		"""
		NGS_ALIGN constructor
		"""
		self.args = args
		self.qsub_cmd_file = os.path.join(args.prj_dir, "pipeline_run.sh")
		log_file = os.path.join(args.prj_dir, "ngs_align.log")
		if os.path.exists(log_file):
			os.remove(log_file)
		log_level_mapping = {"d":logging.DEBUG, "i":logging.INFO, "w":logging.WARNING, "e":logging.ERROR}
		logging.basicConfig(filename=log_file, format='%(levelname)s %(asctime)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S', level=log_level_mapping[args.log_level])
		logging.info("[COMMAND_LINE]: %s" % " ".join(sys.argv))
		self.__create_align_instances()
		
		if args.submit_job:
			my_cmd = "sh " + self.qsub_cmd_file
			print my_cmd
			logging.info(my_cmd)
			os.system("sh " + self.qsub_cmd_file)
			print "Jobs have been submitted!"

	def __create_align_instances(self):
		"""
		create one instance for fq1 [or with fq2]
		"""
		gz_file_list = []
		for root, dirs, files in os.walk(self.args.input_dir, topdown=False):
		    for filename in files:
			if filename.endswith("fastq.gz") or filename.endswith("fq.gz"):
				gz_file_list.append(os.path.join(root, filename))

		fq1_fq2_mapping = {}
		for fq1 in sorted(gz_file_list):
			if self.args.R1 in fq1:
				fq2 = fq1.replace(self.args.R1, self.args.R2)
				if os.path.exists(fq2):
					fq1_fq2_mapping[fq1] = fq2
				else:
					fq1_fq2_mapping[fq1] = ""  # i.e. fq2 = ""
			elif self.args.R2 not in fq1:
				fq1_fq2_mapping[fq1] = ""  # i.e. fq2 = ""

		fs   = FS()
		
		for fq1, fq2 in sorted(fq1_fq2_mapping.iteritems()):
			logging.info("FASTQ(S): %s\t%s" %(fq1, fq2))
			sample_name= os.path.basename(fq1).replace(self.args.R1, "").replace(".fastq.gz", "").replace(".fq.gz", "")
			pipeline_cmd_file = os.path.join(self.args.prj_dir, "alignment", sample_name,"pipeline_commands.sh")
			pipeline_log_file = os.path.join(self.args.prj_dir, "alignment", sample_name,"pipeline.log")
			fs.write(self.qsub_cmd_file, "nohup qsub_cmd.py %s >%s 2>&1 & \n" % (pipeline_cmd_file, pipeline_log_file))
			align = ALIGN_FACTORY().create_instance(self.args, fq1, fq2)			
			align.print_cmd()
		fs.close()		

	def check(self):
		"""
		check if the utility tools and references are already installed
		"""
		#if we don't find some tool, we just print out [WARNING], but still output pipeline_cmds.sh, and put WARNING as comments in pipeline_cmds.sh
		pass	
		

	def __del__(self):
		"""
		NGS_ALIGN destructor 
		"""
		pass
		
	
		
def get_arguments():
	main_description = '''\
	NGS align - create NGS read mapping pipeline command files and the qsub command file
	Author: Zhengdeng Lei (zlei2@uic.edu)
	'''
	help_help = '''\
	show this help message and exit\
	'''
	version_help = '''\
	show the version of this program\
	'''
	input_dir_help = '''\
	input directory containing fq.gz or fastq.gz files
	'''

	prj_dir_help = '''\
	output directory, i.e. project directory, default is the current directory
	pipeline_commands file: [prj_dir]/alignment/[sample_name]/pipeline_commands.sh
	qc file:		[prj_dir]/alignment/[sample_name]/data
	bam file:		[prj_dir]/alignment/[sample_name]/mapped/*.bam

	'''

	genome_help = '''\
	genome, e.g. human, human_hg38, mouse, etc
	'''

	seq_type_help = '''\
	sequencing type: dna, rna, rna3p
	'''

	check_help = '''\
	check if all the utilities and references exist
	'''

	R1_help = '''\
	fastq R1 identifier
	default _R1
	'''

	R2_help = '''\
	fastq R2 identifier
	default _R2
	'''
	mapper_help = '''\
	mapper algorithm to use (bwa, bwamem, star, bowtie2, or tophat2)
	dna	defualt bwamem
	rna	default star
	rna3p	defualt bwamem
	'''

	# GATK/3'RNA-seq bwamem
	# DEseq/edgeR star
	# cufflinks tophat2
	num_thread_help = '''\
	number of threads
	default 8
	'''
	dedup_help = '''\
	remove PCR duplicates,
	default n, i.e. not to remove PCR duplicates
	'''
	opts_help = '''\
	extra options to pass to a utility program
	e.g. -opts "bwamem":"-k 15 -v 4"
	'''	
	log_level_help = '''\
	log level:
	d - logging.DEBUG
	i - logging.INFO
	w - logging.WARNING
	e - logging.ERROR
	default i, i.e. logging.INFO
	'''	
	submit_job_help = '''\
	submit the jobs to PBS cluster
	default y, i.e. submit the jobs
	'''		
	
	arg_parser = argparse.ArgumentParser(description=main_description, formatter_class=argparse.RawTextHelpFormatter, add_help=False)
	arg_parser.register('type','bool', lambda s: str(s).lower() in ['true', '1', 't', 'y', 'yes']) # add type keyword to registries

	###############################
	#    1. required arguments    #
	###############################
	required_group = arg_parser.add_argument_group("required arguments")
	required_group.add_argument("-i", dest="input_dir",  action="store", required=True, default=None, help=input_dir_help)
	required_group.add_argument("-g", dest="genome",   action="store", required=True, default=None, help=genome_help)
	required_group.add_argument("-s", dest="seq_type", action="store", required=True, default=None, help=seq_type_help)

	###############################
	#    2. optional arguments    #
	###############################
	optional_group = arg_parser.add_argument_group("optional arguments")
	optional_group.add_argument("-c", dest="check", default=None, help=check_help)
	optional_group.add_argument("-o", dest="prj_dir", default=os.getcwd(), help=prj_dir_help)
	optional_group.add_argument("-R1", dest="R1", default="_R1", help=R1_help)
	optional_group.add_argument("-R2", dest="R2", default="_R2", help=R2_help)
	optional_group.add_argument("-m",  dest="mapper", default=None, help=mapper_help)
	optional_group.add_argument("-t",  dest="num_thread", type=int, default=8, help=num_thread_help)
	optional_group.add_argument("-d", dest="dedup", type='bool', default=False, help=dedup_help)
	optional_group.add_argument("-opts",  dest="opts", default=None, help=opts_help)
	optional_group.add_argument("-l", dest="log_level", default="i", help=log_level_help)	
	optional_group.add_argument("-q", dest="submit_job", type='bool', default=True, help=submit_job_help)
	
	#optional_group.add_argument("-rg", dest="read_group", type='bool', default=False, help=read_group_help)

	optional_group.add_argument("-h", "--help", action="help", help=help_help)
	optional_group.add_argument("-v", "--version", action="version", version="%(prog)s: version 1.0", help=version_help)

	args = arg_parser.parse_args()
	
	args.genome = args.genome.lower()
	if args.mapper:
		args.mapper = args.mapper.lower()		
	args.input_dir = os.path.abspath(args.input_dir)
	args.prj_dir = os.path.abspath(args.prj_dir)

	if not os.path.exists(args.input_dir):
		arg_parser.print_help()
		print "\n\nError: source directory does not exist!\n"
		sys.exit()

	create_dirs = ["prj_dir"]
	for create_dir in create_dirs:
		dir_path = eval("args." + create_dir)
		if not os.path.exists(dir_path):
			print dir_path, "is created!\n"
			os.makedirs(dir_path)

	return args



def main():
	args = get_arguments()
	ngs_align = NGS_ALIGN(args)
	#factory = ALIGN_FACTORY()
	#ngs = factory.get_sequencing(args)
	#ngs.print_cmd()
		
if __name__ == '__main__':
	main()
