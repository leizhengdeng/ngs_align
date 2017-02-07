#!/usr/bin/env python
import os
import sys
import argparse
import logging
from fs import FS
from cmd_list import CMD_LIST
#import subprocess

class ALIGN:
	"""
	ALIGN: base class of DNA_ALIGN, RNA_ALIGN, RNA3P_ALIGN etc.
	print the command lines to the file: [prj_dir]/alignment/[sample_name]/pipeline_commands.sh
	"""
	def print_cmd(self):
		fs = self.fs
		cmd_file = self.pipeline_cmd_file
		
		fs.write(cmd_file, "#Generating fastqc reports\n")
		for my_cmd in self.do_qc():
			fs.write(cmd_file, my_cmd + "\n")
		fs.write(cmd_file, "\n\n")	

		fs.write(cmd_file, "#do alignment\n")
		for my_cmd in self.do_align():
			fs.write(cmd_file, my_cmd + "\n\n")
		fs.write(cmd_file, "\n\n")	

	def __init__(self, args, fq1, fq2):
   		"""
   		ALIGN constructor
   		"""
   		if not CMD_LIST.initialized:
   	  		CMD_LIST.initialize(args)
		self.fs = FS()
		self.fq1 = fq1
		self.fq2 = fq2

		self.args = args
		self.sample_name= os.path.basename(fq1).replace(args.R1, "").replace(".fastq.gz", "").replace(".fq.gz", "")
		self.alignment_dir	= os.path.join(self.args.prj_dir, "alignment")
		self.sample_dir	= os.path.join(self.args.prj_dir, "alignment", self.sample_name)
		self.data_dir	= os.path.join(self.args.prj_dir, "alignment", self.sample_name, "data")
		self.mapped_dir	= os.path.join(self.args.prj_dir, "alignment", self.sample_name, "mapped")

		create_dirs = [self.alignment_dir, self.sample_dir, self.data_dir, self.mapped_dir]
		for dir_path in create_dirs:
			if not os.path.exists(dir_path):
				os.makedirs(dir_path)

		self.pipeline_cmd_file = os.path.join(self.args.prj_dir, "alignment", self.sample_name,"pipeline_commands.sh")
	
	def do_qc(self):
		"""
		fastqc commands
		return a generator object for fastqc1_cmd and fastqc2_cmd
		"""
		fastqc1_cmd = CMD_LIST.fastqc.replace("OUT_DIR", self.data_dir).replace("FASTQ", self.fq1)
		#p = subprocess.Popen(CMD_LIST.count_reads_fastq_gz.replace("FASTQ", self.fq1), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		#[output, errmsg] = p.communicate()
		#(total_reads, total_bases) = output.split(" ")
		#logging.info("Total reads: %s, total bases: %s  in %s" % (total_reads, total_bases, self.fq1))
		yield fastqc1_cmd
		if self.fq2:
			fastqc2_cmd = CMD_LIST.fastqc.replace("OUT_DIR", self.data_dir).replace("FASTQ", self.fq2)
			yield fastqc2_cmd

	def do_align(self):
		"""
		interface only
		do_align() will be implemented in derived classes e.g. DNA_ALIGN or RNA_ALIGN
		"""
		pass

	def __del__(self):
	   	"""
   		ALIGN destructor
   		"""
		self.fs.close()		
	