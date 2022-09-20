
#!/usr/bin/bash
SAMPLE=test01
GENOME=~/ref/bwa/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
INDEX='~/ref/bwa/pig/Sus_scrofa'
GATK=~/software/gatk/gatk-package-4.1.9.0-local.jar
KNOWSITE1=~/ref/bwa/pig/sus_scrofa.vcf.gz
GATKMEM=6g
fq1=test01_1_out.fastq.gz       #cleandata
fq2=test01_2_out.fastq.gz

#比对(850Mb原始双端数据，2cpu消耗55min)
bwa mem -t 3 -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tLB:WGS\tPL:Illumina" $INDEX $fq1 $fq2 | samtools view -@3 -b > bam/$SAMPLE.bam  #(850Mb原始数据产生2.1g bam文件)
#排序(850Mb原始双端数据，2cpu消耗3min)
samtools sort -@3 bam/$SAMPLE.bam -o bam/$SAMPLE.sort.bam  #(850Mb原始数据产生1.9g bam文件)
#标记PCR重复 (850Mb原始双端数据，2cpu消耗10min，3g运行内存)
java -Xmx$GATKMEM -jar $GATK MarkDuplicates  --TMP_DIR /Work/user/huanong/tmpdir  -I bam/$SAMPLE.sort.bam -O bam/markdup/$SAMPLE.markdup.bam -M bam/markdup/$SAMPLE.markdup_metrics.txt #(850Mb原始数据产生2.6g bam文件)
#创建比对索引文件
samtools index bam/markdup/$SAMPLE.markdup.bam
#FixMateInformation(确认成对 reads 之间的信息是否一致；如有必要，进行修复)(850Mb原始双端数据，2cpu消耗12min，3g运行内存)
java -Xmx$GATKMEM -Djava.io.tmpdir=/Work/user/huanong/tmpdir -jar $GATK  FixMateInformation I=bam/markdup/$SAMPLE.markdup.bam O=bam/markdup/SNP_calling/$SAMPLE.fixed.bam SO=coordinate  TMP_DIR=/Work/user/huanong/tmpdir  #(850Mb原始数据产生2.6g bam文件)
#BQSR质量校正对比  (第一步：850Mb原始数据，2cpu消耗20min,2g运行内存；第二步：850Mb原始数据，2cpu消耗6min,1g运行内存。)
java -Xmx$GATKMEM -jar $GATK BaseRecalibrator -R $GENOME -I bam/markdup/SNP_calling/$SAMPLE.fixed.bam --known-sites $KNOWSITE1  -O bam/markdup/SNP_calling/$SAMPLE.table    
java -Xmx$GATKMEM -jar $GATK ApplyBQSR -R $GENOME -I bam/markdup/SNP_calling/$SAMPLE.fixed.bam -bqsr bam/markdup/SNP_calling/$SAMPLE.table -O bam/markdup/SNP_calling/$SAMPLE.BQSR.bam  #(850Mb原始数据产生2.6g bam文件)
#SNP-calling (850Mb原始数据，2cpu消耗330min,3g运行内存)
java -Xmx$GATKMEM -jar $GATK HaplotypeCaller -R $GENOME -I bam/markdup/SNP_calling/$SAMPLE.BQSR.bam -ERC GVCF -O bam/markdup/SNP_calling/$SAMPLE.g.vcf
#压缩文件
bgzip bam/markdup/SNP_calling/$SAMPLE.g.vcf
#建立索引
tabix -p vcf bam/markdup/SNP_calling/$SAMPLE.g.vcf.gz
#删除无用文件(判断最后结果是否存在，如果存在删除无用文件并输出运行成功，否则输出报错)
if [ -f bam/markdup/SNP_calling/$SAMPLE.g.vcf.gz ]; then rm -f bam/$SAMPLE.bam bam/$SAMPLE.sort.bam  bam/markdup/SNP_calling/$SAMPLE.BQSR.bam bam/markdup/SNP_calling/$SAMPLE.fixed.bam bam/markdup/SNP_calling/$SAMPLE.table  bam/markdup/SNP_calling/$SAMPLE.BQSR.bai ;echo "$SAMPLE is done! " >>SNP_Calliing.info ;else echo "$SAMPLE is error! " >>SNP_Calliing.info ;fi

