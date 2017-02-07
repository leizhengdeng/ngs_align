import sys

class MAPPER:
	# command seperator is 	##################

	star = '''\
	MAPPER
	--readFilesCommand gzip -dc
	--readFilesIn FASTQ1 FASTQ2
	--outSAMtype BAM SortedByCoordinate
	--genomeDir GENOME_DIR
	--runThreadN NUM_THREAD
	--outFileNamePrefix BAM_PREFIX
	
	##################
	
	'''

	bwamem = '''\
	MAPPER
	mem dddd

	##################

	rm ...
	'''

	bowtie2 = '''\
	'''
	
