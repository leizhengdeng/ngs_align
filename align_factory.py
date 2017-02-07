import sys
import logging

class ALIGN_FACTORY:
	"""align factory"""
	def create_instance(self, args, fq1, fq2):
		logging.debug("create align instance in ALIGN_FACTORY")
		seq_type = args.seq_type.lower()
		from dna_align import DNA_ALIGN
		from rna_align import RNA_ALIGN
		from rna3p_align import RNA3P_ALIGN	

		if seq_type == "dna":
			return DNA_ALIGN(args, fq1, fq2)
		elif seq_type == "rna":
			return RNA_ALIGN(args, fq1, fq2)
		elif seq_type == "rna3p":
			return RNA3P_ALIGN(args, fq1, fq2)
		else:
			print "sequencing type is not supported yet"
			sys.exit()

