import sys
import re
import logging
from util import UTIL
from ref import REF

class CMD_LIST:
	"""
	CMD_LIST: list all commands used in the pipeline
	"""

	# fastqc for raw read quality control
	fastqc = '''
	fastqc -q -t NUM_THREAD -o OUT_DIR FASTQ
	'''
	
	# bwa mem aligner
	# read group can be added using -R
	bwamem = '''
	bwa mem -t NUM_THREAD bwa_GENOME_index FASTQ1 FASTQ2 | samtools view -bS - > MAPPED/SAMPLE_NAME.bam
	'''
	
	# clean bam file, i.e. remove unmapped reads from bam file
	cleanbam = '''
	samtools view -h -F4 MAPPED/SAMPLE_NAME.bam |
	samtools view -bS - > MAPPED/SAMPLE_NAME_tmp.bam
	'''

	mvbam = '''
	mv MAPPED/SAMPLE_NAME_tmp.bam MAPPED/SAMPLE_NAME.bam
	'''

	# clean star bam file, i.e. remove unmapped reads from star bam file
	cleanstarbam = '''
	samtools view -h -F4 MAPPED/SAMPLE_NAMEAligned.out.bam |
	samtools view -bS - > MAPPED/SAMPLE_NAME.bam
	'''

	# sort bam file
	sortbam = '''
	java -jar sortsam TMP_DIR=MAPPED VALIDATION_STRINGENCY=LENIENT 
	I=MAPPED/SAMPLE_NAME.bam O=MAPPED/SAMPLE_NAME_tmp.bam SO=coordinate
	'''

	# create index for bam file
	indexbam = '''
	samtools index MAPPED/SAMPLE_NAME.bam
	'''
	
	# removing PCR duplicates
	dedup = '''
	java -jar markduplicates TMP_DIR=MAPPED VALIDATION_STRINGENCY=LENIENT 
	INPUT=MAPPED/SAMPLE_NAME.bam OUTPUT=MAPPED/SAMPLE_NAME_tmp.bam 
	METRICS_FILE=MAPPED/SAMPLE_NAME_duplicates.txt REMOVE_DUPLICATES=true ASSUME_SORTED=true
	'''
	
	# qualimap to get mapping stats
	qualimap = '''
	qualimap --java-mem-size=32G bamqc -bam MAPPED/SAMPLE_NAME.bam -outdir MAPPED/mapping_QC -nt NUM_THREAD
	'''

	# qualimap with annotation file to get mapping stats
	qualimapgff = '''
	qualimap --java-mem-size=32G bamqc -bam MAPPED/SAMPLE_NAME.bam -outdir MAPPED/mapping_QC -nt NUM_THREAD
	-gff GENOME_gtf
	'''
	
	#star aligner
	star = '''
	star
	--readFilesCommand gzip -dc
	--readFilesIn FASTQ1 FASTQ2	
	--genomeDir star_GENOME_index
	--runThreadN NUM_THREAD
	--outFilterIntronMotifs RemoveNoncanonical 	
	--outFileNamePrefix MAPPED/SAMPLE_NAME
	--outSAMtype BAM Unsorted
	--outSAMunmapped Within
	'''

	#Generate genome:
	#STAR --runThreadN 40 
	#--runMode genomeGenerate 
	#--genomeDir [genome dir with star index] 
	#--genomeFastaFiles [*.fa]
	#--sjdbGTFfile [*.gtf]

	#count_reads_fastq_gz = '''
	#gzip -dc FASTQ | echo $((`wc -l`/4))
	#'''
	
	#count_reads_bam = '''
	#samtools view -F4 MAPPED/SAMPLE_NAME.bam | cut -f1 | sort | uniq | wc -l
	#'''
	
	#trimmomatic = '''
	#java -jar trimmomatic -t NUM_THREAD FASTQ1
	#'''
	
	#samtools mpileup -A -t SP -u -f hg19_chr22.fa SRR2968056_chr22_nodups.bam | bcftools call -vc - > SRR2968056_chr22.vcf
	#bgzip -c SRR2968056_chr22.vcf > SRR2968056_chr22.vcf.gz
	#tabix -p vcf SRR2968056_chr22.vcf.gz
	

	initialized = False
	cmd2comment = {\
	"bwamem":"read mapping with bwa mem",\
	"cleanbam":"clean bam file, i.e. i.e. remove unmapped reads from bam file"\
	}
	
	@staticmethod
	def initialize(args):
		members = [attr for attr in dir(CMD_LIST) if not callable(attr) and not attr.startswith("__")]
		members.remove("initialized")
		members.remove("initialize")
		members.remove("cmd2comment")
		if args.genome not in REF.name2path:
			print "[ERROR] %s(case insensitive) is not found in REF.ini\n\n" % args.genome
			logging.error("[CMD_LIST] %s(case insensitive) is not found in REF.ini\n\n" % args.genome)
			sys.exit()			
		for m in members:
			strip = "CMD_LIST.%s = CMD_LIST.%s.strip()" % (m, m)	
			replace_tab = "CMD_LIST.%s = CMD_LIST.%s.replace(\"\\t\", \" \")" % (m, m)
			replace_linesep = "CMD_LIST.%s = CMD_LIST.%s.replace(\"\\n\", \" \")" % (m, m)
			replace_thread = "CMD_LIST.%s = CMD_LIST.%s.replace(\"NUM_THREAD\", \"%d\")" % (m, m, args.num_thread)
			replace_genome = "CMD_LIST.%s = CMD_LIST.%s.replace(\"GENOME\", \"%s\")" % (m, m, args.genome)
			extra_space = "CMD_LIST.%s = CMD_LIST.%s + \" \"" % (m, m) #extra space is needed in the end of each cmd
			exec(strip)
			exec(replace_tab)
			exec(replace_linesep)					
			exec(replace_thread)
			exec(replace_genome)
			exec(extra_space)

			# replace the util name and reference name by using the information in UTIL.ini and REF.ini
			logging.debug("[CMD_LIST] before using UTIL and REF: " + eval("CMD_LIST.%s" % m))
			for name, path in UTIL.name2path.iteritems():
				replace_util = "CMD_LIST.%s = re.sub(\"%s\s+\", \"%s \", CMD_LIST.%s)" % (m, name, path, m)
				exec(replace_util)
			for name, path in REF.name2path.iteritems():
				replace_ref = "CMD_LIST.%s = re.sub(\"%s\s+\", \"%s \", CMD_LIST.%s)" % (m, name, path, m)
				exec(replace_ref)				
			logging.debug("[CMD_LIST] after uing UTIL and REF: " + eval("CMD_LIST.%s" % m))

		CMD_LIST.initialized = True
		
