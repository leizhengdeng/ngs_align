import sys
import re
import logging
from align import ALIGN
from cmd_list import CMD_LIST

class RNA_ALIGN(ALIGN):
	def __init__(self, args, fq1, fq2):
		"""
		RNA_ALIGN constructor
		"""
		ALIGN.__init__(self, args, fq1, fq2)
		self.args = args
		self.fq1 = fq1
		self.fq2 = fq2
		if self.args.mapper is None:
			# default mapper for rna is star
			self.args.mapper = "star"	
			
	def do_align(self):
		"""
		RNA_ALIGN.do_align: for each command in CMD_LIST, replace with sample name and bam output directory, as well as fq1 and fq2"
		"""
		my_cmd = ""
		cmd_order_list = ["star", "cleanstarbam", "sortbam", "mvbam", "indexbam", "qualimapgff"]
		for m in cmd_order_list:
			replace_mapped = "my_cmd = re.sub(\"MAPPED\", \"%s\", CMD_LIST.%s)" % (self.mapped_dir, m)
			replace_sample_name = "my_cmd = re.sub(\"SAMPLE_NAME\", \"%s\", my_cmd)" % (self.sample_name)
			replace_fastq1 = "my_cmd = re.sub(\"FASTQ1\", \"%s\", my_cmd)" % (self.fq1)
			replace_fastq2 = "my_cmd = re.sub(\"FASTQ2\", \"%s\", my_cmd)" % (self.fq2)
			exec(replace_mapped)
			exec(replace_sample_name)
			exec(replace_fastq1)
			exec(replace_fastq2)	
			logging.debug("[RNA_ALIGN] after replacement: " + my_cmd)
			yield my_cmd
			
	def __del__(self):
		ALIGN.__del__(self)




