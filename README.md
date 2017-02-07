# ngs_align
A pipeline for NGS read alignment (fastq to bam) for DNA-seq, RNA-seq:
<br>
<br>
$ngs_align.py -h<br>
usage: ngs_align.py -i INPUT_DIR -g GENOME -s SEQ_TYPE [-c CHECK] [-o PRJ_DIR]<br>
                    [-R1 R1] [-R2 R2] [-m MAPPER] [-t NUM_THREAD] [-d DEDUP]<br>
                    [-opts OPTS] [-l LOG_LEVEL] [-q SUBMIT_JOB] [-h] [-v]<br>
<br>
        NGS align - create NGS read mapping pipeline command files and the qsub command file<br>
        Author: Zhengdeng Lei (zlei2@uic.edu)<br>
<br>
<br>
required arguments:<br>
  -i INPUT_DIR          input directory containing fq.gz or fastq.gz files<br>
<br>
  -g GENOME             genome, e.g. human, human_hg38, mouse, etc<br>
<br>
  -s SEQ_TYPE           sequencing type: dna, rna, rna3p<br>
<br>
<br>
optional arguments:<br>
  -c CHECK              check if all the utilities and references exist<br>
<br>
  -o PRJ_DIR            output directory, i.e. project directory, default is the current directory<br>
                        pipeline_commands file: [prj_dir]/alignment/[sample_name]/pipeline_commands.sh<br>
                        qc file:                [prj_dir]/alignment/[sample_name]/data<br>
                        bam file:               [prj_dir]/alignment/[sample_name]/mapped/*.bam<br>
<br>
<br>
  -R1 R1                fastq R1 identifier<br>
                        default _R1<br>
<br>
  -R2 R2                fastq R2 identifier<br>
                        default _R2<br>
<br>
  -m MAPPER             mapper algorithm to use (bwa, bwamem, star, bowtie2, or tophat2)<br>
                        dna     defualt bwamem<br>
                        rna     default star<br>
                        rna3p   defualt bwamem<br>
<br>
  -t NUM_THREAD         number of threads<br>
                        default 8<br>
<br>
  -d DEDUP              remove PCR duplicates,<br>
                        default n, i.e. not to remove PCR duplicates<br>
<br>
  -opts OPTS            extra options to pass to a utility program<br>
                        e.g. -opts "bwamem":"-k 15 -v 4"<br>
<br>
  -l LOG_LEVEL          log level:<br>
                        d - logging.DEBUG<br>
                        i - logging.INFO<br>
                        w - logging.WARNING<br>
                        e - logging.ERROR<br>
                        default i, i.e. logging.INFO<br>
<br>
  -q SUBMIT_JOB         submit the jobs to PBS cluster<br>
                        default y, i.e. submit the jobs<br>
<br>
  -h, --help            show this help message and exit<br>
  -v, --version         show the version of this program<br>

