#Generating fastqc reports
/export/home/clustcrilab/tools/FastQC/fastqc -q -t 8 -o /cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/data /cri/home2/zlei2/projects/2016.05-bxayarath/data/BX-DM-1_S4_R1.fastq.gz 
/export/home/clustcrilab/tools/FastQC/fastqc -q -t 8 -o /cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/data /cri/home2/zlei2/projects/2016.05-bxayarath/data/BX-DM-1_S4_R2.fastq.gz 


#do alignment
/export/home/clustcrilab/tools/bwa-0.7.5a/bwa mem -t 8 /export/home/clustcrilab/reference/igenome_ucsc_hg19/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa /cri/home2/zlei2/projects/2016.05-bxayarath/data/BX-DM-1_S4_R1.fastq.gz /cri/home2/zlei2/projects/2016.05-bxayarath/data/BX-DM-1_S4_R2.fastq.gz | /export/home/clustcrilab/tools/samtools-0.1.19/samtools view -bS - > /cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped/BX-DM-1_S4.bam 

/export/home/clustcrilab/tools/samtools-0.1.19/samtools view -h -F4 /cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped/BX-DM-1_S4.bam |  /export/home/clustcrilab/tools/samtools-0.1.19/samtools view -bS - > /cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped/BX-DM-1_S4_tmp.bam 

mv /cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped/BX-DM-1_S4_tmp.bam /cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped/BX-DM-1_S4.bam 

/export/share/compilers/java-1.8/jdk1.8.0_73/bin/java -jar /export/home/clustcrilab/tools/picard-tools-1.107/SortSam.jar TMP_DIR=/cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped VALIDATION_STRINGENCY=LENIENT   I=/cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped/BX-DM-1_S4.bam O=/cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped/BX-DM-1_S4_tmp.bam SO=coordinate 

mv /cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped/BX-DM-1_S4_tmp.bam /cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped/BX-DM-1_S4.bam 

/export/home/clustcrilab/tools/samtools-0.1.19/samtools index /cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped/BX-DM-1_S4.bam 

/export/home/clustcrilab/tools/qualimap_v0.7.1/qualimap --java-mem-size=32G bamqc -bam /cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped/BX-DM-1_S4.bam -outdir /cri/home2/zlei2/projects/2016.05-bxayarath/alignment/BX-DM-1_S4/mapped/mapping_QC -nt 8 



