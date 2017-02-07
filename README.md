# ngs_align
A pipeline for NGS read alignment (fastq to bam) for DNA-seq, RNA-seq:

$ngs_align.py -h
usage: ngs_align.py -i INPUT_DIR -g GENOME -s SEQ_TYPE [-c CHECK] [-o PRJ_DIR]
                    [-R1 R1] [-R2 R2] [-m MAPPER] [-t NUM_THREAD] [-d DEDUP]
                    [-opts OPTS] [-l LOG_LEVEL] [-q SUBMIT_JOB] [-h] [-v]

        NGS align - create NGS read mapping pipeline command files and the qsub command file
        Author: Zhengdeng Lei (zlei2@uic.edu)


required arguments:
  -i INPUT_DIR          input directory containing fq.gz or fastq.gz files

  -g GENOME             genome, e.g. human, human_hg38, mouse, etc

  -s SEQ_TYPE           sequencing type: dna, rna, rna3p


optional arguments:
  -c CHECK              check if all the utilities and references exist

  -o PRJ_DIR            output directory, i.e. project directory, default is the current directory
                        pipeline_commands file: [prj_dir]/alignment/[sample_name]/pipeline_commands.sh
                        qc file:                [prj_dir]/alignment/[sample_name]/data
                        bam file:               [prj_dir]/alignment/[sample_name]/mapped/*.bam


  -R1 R1                fastq R1 identifier
                        default _R1

  -R2 R2                fastq R2 identifier
                        default _R2

  -m MAPPER             mapper algorithm to use (bwa, bwamem, star, bowtie2, or tophat2)
                        dna     defualt bwamem
                        rna     default star
                        rna3p   defualt bwamem

  -t NUM_THREAD         number of threads
                        default 8

  -d DEDUP              remove PCR duplicates,
                        default n, i.e. not to remove PCR duplicates

  -opts OPTS            extra options to pass to a utility program
                        e.g. -opts "bwamem":"-k 15 -v 4"

  -l LOG_LEVEL          log level:
                        d - logging.DEBUG
                        i - logging.INFO
                        w - logging.WARNING
                        e - logging.ERROR
                        default i, i.e. logging.INFO

  -q SUBMIT_JOB         submit the jobs to PBS cluster
                        default y, i.e. submit the jobs

  -h, --help            show this help message and exit
  -v, --version         show the version of this program

