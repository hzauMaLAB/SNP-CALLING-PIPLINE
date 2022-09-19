**预先安装的软件有`gatk4.1.9`、`bwa`、`samtools`、`bcftools`**
- 传输文件
```
/Work/user/huanong/software/obsutil cp obs://baiyi/HuaNongShuJu/sus_scrofa_incl_consequences.vcf.gz ./
```
- 先创建储存文件的目录
```
cd /Work/user/huanong/project/Mangalica_clean_reads     #进入梅山猪的数据目录
mkdir -p  bam/markdup/SNP_calling  #创建文件夹
```
- 给参考基因组以及dbsnp数据建立索引
```
#下载参考基因组并建立索引
cd /Work/user/huanong/ref/bwa/
wget http://ftp.ensembl.org/pub/release-106/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
gzip -d Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz  #解压
samtools faidx  Sus_scrofa.Sscrofa11.1.dna.toplevel.fa   #建立索引
#下载dbsnp数据文件
cd /Work/user/huanong/ref/vcf/
wget  http://ftp.ensembl.org/pub/release-106/variation/vcf/sus_scrofa/sus_scrofa.vcf.gz
#处理dbsnp数据文件（因为不处理会报错：info那一列不允许有空格出现）
zcat  sus_scrofa.vcf.gz|grep "#" >header.tmp 
zcat  sus_scrofa.vcf.gz|grep -v "#"|tr " " "_" >body.tmp
rm -f sus_scrofa.vcf.gz
cat header.tmp body.tmp >sus_scrofa.vcf
rm -f *.tmp
#dbsnp数据文件建立索引
bgzip sus_scrofa.vcf
tabix -p vcf sus_scrofa.vcf.gz
```
- bwa建立索引
```
cd /Work/user/huanong/ref/bwa
bwa index -p Sus_scrofa Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
```
- 创建单个样本callsnp示例脚本
```
cd /Work/user/huanong/project  #进入存放质控过的数据目录
vim SNP_CALLtest01.sh
```
```
#!/usr/bin/bash
SAMPLE=test01
GENOME=/Work/user/huanong/ref/bwa/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
INDEX='/Work/user/huanong/ref/bwa/Sus_scrofa'
GATK=/Work/user/huanong/software/gatk/gatk-package-4.1.9.0-local.jar
KNOWSITE1=/Work/user/huanong/ref/vcf/sus_scrofa.vcf.gz
GATKMEM=6g
fq1=test01_1_out.fastq.gz
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
```
- 运行
```
for i in `ls *_1_*|tr "_" "\t"|cut -f1`;do sed "s/test01/$i/g" ../SNP_CALLtest01.sh >SNP_CALL$i.sh;done
for i in `ls *_1_*|tr "_" "\t"|cut -f1`;do qsub -V -cwd -q all.q,fat.q -l vf=6G -pe smp 3 -S /bin/bash  SNP_CALL$i.sh;done
```
- gvcf合并(要等所有样本的gvcf都出来之后，可以分群体也可以一块合并，如果分群体最后也是要一块合并的)
**这一步可以分染色进行，建议不要在我们自己服务器上跑，会非常慢**
```
find *.gz > input.list
nohup java -Xmx50g -jar /public/jychu/soft/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar CombineGVCFs -R /public/jychu/refs/Gallus_gallus.GRCg6a.dna.toplevel.fa -V  input.list  -O ~/chicken_body_size/VCF/ALLSAMPLE.vcf.gz &
```
###接下来是在作重服务器上的操作
- gvcf合并后提取猪的不同染色体的vcf（分染色体）  可以用vcftools和bcftools提取，bcftools提取的时候会和原来的文件一模一样，vcftools提取的时候info那一列好像会改变一些信息，这个问题我没有仔细去研究
1.vcftools提取
```
bsub -J chr1 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 1 --recode --stdout | bgzip -c >pig_chr1.vcf.gz"
bsub -J chr2 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 2 --recode --stdout | bgzip -c >pig_chr2.vcf.gz"
bsub -J chr3 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 3 --recode --stdout | bgzip -c >pig_chr3.vcf.gz"
bsub -J chr4 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 4 --recode --stdout | bgzip -c >pig_chr4.vcf.gz"
bsub -J chr5 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 5 --recode --stdout | bgzip -c >pig_chr5.vcf.gz"
bsub -J chr6 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 6 --recode --stdout | bgzip -c >pig_chr6.vcf.gz"
bsub -J chr7 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 7 --recode --stdout | bgzip -c >pig_chr7.vcf.gz"
bsub -J chr8 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 8 --recode --stdout | bgzip -c >pig_chr8.vcf.gz"
bsub -J chr9 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 9 --recode --stdout | bgzip -c >pig_chr9.vcf.gz"
bsub -J chr10 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 10 --recode --stdout | bgzip -c >pig_chr10.vcf.gz"
bsub -J chr11 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 11 --recode --stdout | bgzip -c >pig_chr11.vcf.gz"
bsub -J chr12 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 12 --recode --stdout | bgzip -c >pig_chr12.vcf.gz"
bsub -J chr13 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 13 --recode --stdout | bgzip -c >pig_chr13.vcf.gz"
bsub -J chr14 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 14 --recode --stdout | bgzip -c >pign_chr14.vcf.gz"
bsub -J chr15 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 15 --recode --stdout | bgzip -c >pig_chr15.vcf.gz"
bsub -J chr16 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 16 --recode --stdout | bgzip -c >pig_chr16.vcf.gz"
bsub -J chr17 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 17 --recode --stdout | bgzip -c >pig_chr17.vcf.gz"
bsub -J chr18 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr 18 --recode --stdout | bgzip -c >pig_chr18.vcf.gz"
bsub -J chrX -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr X --recode --stdout | bgzip -c >pig_chrX.vcf.gz"
bsub -J chrY -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr Y --recode --stdout | bgzip -c >pig_chrY.vcf.gz"
bsub -J chrMT -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "vcftools --gzvcf combine.vcf.gz --chr MT --recode --stdout | bgzip -c >pig_chrMT.vcf.gz"
```
2.bcftools提取
```
vi chrMT.region
MT      1       16775      #tab键分割
bsub -J chrMT -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chrMT.region combine.vcf.gz --output bcftools/chrMT.vcf.gz --output-type z"
vi chrW.region
W    1    6813114
bsub -J chrW -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chrW.region combine.vcf.gz --output bcftools/chrW.vcf.gz --output-type z"
vi chrZ.region
Z    1    82529921
bsub -J chrZ -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chrZ.region combine.vcf.gz --output bcftools/chrZ.vcf.gz --output-type z"
vi chr1.region
1    1    197608386
bsub -J chr1 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr1.region combine.vcf.gz --output bcftools/chr1.vcf.gz --output-type z"
vi chr2.region
2    1    149682049
bsub -J chr2 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr2.region combine.vcf.gz --output bcftools/chr2.vcf.gz --output-type z"
vi chr3.region
3    1    110838418
bsub -J chr3 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr3.region combine.vcf.gz --output bcftools/chr3.vcf.gz --output-type z"
vi chr4.region
4    1    91315245
bsub -J chr4 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr4.region combine.vcf.gz --output bcftools/chr4.vcf.gz --output-type z"
vi chr5.region
5    1    59809098
bsub -J chr5 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr5.region combine.vcf.gz --output bcftools/chr5.vcf.gz --output-type z"
vi chr6.region
6    1    36374701
bsub -J chr6 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr6.region combine.vcf.gz --output bcftools/chr6.vcf.gz --output-type z"
vi chr7.region
7    1    36742308
bsub -J chr7 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr7.region combine.vcf.gz --output bcftools/chr7.vcf.gz --output-type z"
vi chr8.region
8    1    30219446
bsub -J chr8 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr8.region combine.vcf.gz --output bcftools/chr8.vcf.gz --output-type z"
vi chr9.region
9    1    24153086
bsub -J chr9 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr9.region combine.vcf.gz --output bcftools/chr9.vcf.gz --output-type z"
vi chr10.region
10   1    21119840
bsub -J chr10 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr10.region combine.vcf.gz --output bcftools/chr10.vcf.gz --output-type z"
vi chr11.region
11    1    20200042
bsub -J chr11 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr11.region combine.vcf.gz --output bcftools/chr11.vcf.gz --output-type z"
vi chr12.region
12    1    20387278
bsub -J chr12 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr12.region combine.vcf.gz --output bcftools/chr12.vcf.gz --output-type z"
vi chr13.region
13    1    19166714
bsub -J chr13 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr13.region combine.vcf.gz --output bcftools/chr13.vcf.gz --output-type z"
vi chr14.region
14    1    16219308
bsub -J chr14 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr14.region combine.vcf.gz --output bcftools/chr14.vcf.gz --output-type z"
vi chr15.region
15    1    13062184
bsub -J chr15 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr15.region combine.vcf.gz --output bcftools/chr15.vcf.gz --output-type z"
vi chr16.region
16    1    2844601
bsub -J chr16 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr16.region combine.vcf.gz --output bcftools/chr16.vcf.gz --output-type z"
vi chr17.region
17    1    10762512
bsub -J chr17 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr17.region combine.vcf.gz --output bcftools/chr17.vcf.gz --output-type z"
vi chr18.region
18    1    11373140
bsub -J chr18 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr18.region combine.vcf.gz --output bcftools/chr18.vcf.gz --output-type z"
vi chr19.region
19    1    10323212
bsub -J chr19 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr19.region combine.vcf.gz --output bcftools/chr19.vcf.gz --output-type z"
vi chr20.region
20    1    13897287
bsub -J chr20 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr20.region combine.vcf.gz --output bcftools/chr20.vcf.gz --output-type z"
vi chr21.region
21    1    6844979
bsub -J chr21 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr21.region combine.vcf.gz --output bcftools/chr21.vcf.gz --output-type z"
vi chr22.region
22    1    5459462
bsub -J chr22 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr22.region combine.vcf.gz --output bcftools/chr22.vcf.gz --output-type z"
vi chr23.region
23    1    6149580
bsub -J chr23 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr23.region combine.vcf.gz --output bcftools/chr23.vcf.gz --output-type z"
vi chr24.region
24    1    6491222
bsub -J chr24 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr24.region combine.vcf.gz --output bcftools/chr24.vcf.gz --output-type z"
vi chr25.region
25    1    3980610
bsub -J chr25 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr25.region combine.vcf.gz --output bcftools/chr25.vcf.gz --output-type z"
vi chr26.region
26   1    6055710
bsub -J chr26 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr26.region combine.vcf.gz --output bcftools/chr26.vcf.gz --output-type z"
vi chr27.region
27    1    8080432
bsub -J chr27 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr27.region combine.vcf.gz --output bcftools/chr27.vcf.gz --output-type z"
vi chr28.region
28    1    5116882
bsub -J chr28 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr28.region combine.vcf.gz --output bcftools/chr28.vcf.gz --output-type z"
vi chr30.region
30     1        1818525
bsub -J chr30 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr30.region combine.vcf.gz --output bcftools/chr30.vcf.gz --output-type z"
vi chr31.region
31    1    6153034
bsub -J chr31 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr31.region combine.vcf.gz --output bcftools/chr31.vcf.gz --output-type z"
vi chr32.region
32    1    725831
bsub -J chr32 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr32.region combine.vcf.gz --output bcftools/chr32.vcf.gz --output-type z"
vi chr33.region
33    1    7821666
bsub -J chr33 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr33.region combine.vcf.gz --output bcftools/chr33.vcf.gz --output-type z"
```
- 对提取出来的vcf文件建立索引
```
tabix -p vcf chr*.region combine.vcf.gz       #这样概括所有的染色体，这个命令只能一次操作一个文件
```
- 添加注释并genotype
```
wget http://ftp.ensembl.org/pub/release-103/variation/vcf/gallus_gallus/gallus_gallus.vcf.gz
bsub -J chr10 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "tabix -p vcf gallus_gallus.vcf.gz "
bsub -J chr10 -n 8 -R span[hosts=1] -o %J.out -e %J.err -q normal "java -Xmx80g -jar /public/home/hsong/jychu/soft/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar  GenotypeGVCFs -R /public/home/hsong/jychu/refs/Gallus_gallus.GRCg6a.dna.toplevel.fa -V chr10.vcf.gz --dbsnp gallus_gallus.vcf.gz -O  chr10_typed.vcf.gz"
#chr10_typed.vcf.gz 如果没有索引，也要建立索引
```
- 纵向合并vcf染色体文件(合并的文件必须是bgzip压缩过的文件)   这一步也可以先不合并，最后过滤完snp位点的时候再合并
```
bsub -J chr33 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools concat  chrMT_typed.vcf.gz   chrW_typed.vcf.gz  chrZ_typed.vcf.gz  chr1_typed.vcf.gz   chr2_typed.vcf.gz  chr3_typed.vcf.gz chr4_typed.vcf.gz  chr5_typed.vcf.gz  chr6_typed.vcf.gz  chr7_typed.vcf.gz  chr8_typed.vcf.gz  chr9_typed.vcf.gz  chr10_typed.vcf.gz  chr11_typed.vcf.gz chr12_typed.vcf.gz   chr13_typed.vcf.gz  chr14_typed.vcf.gz  chr15_typed.vcf.gz  chr16_typed.vcf.gz  chr17_typed.vcf.gz chr18_typed.vcf.gz  chr19_typed.vcf.gz chr20_typed.vcf.gz chr21_typed.vcf.gz  chr22_typed.vcf.gz chr23_typed.vcf.gz chr24_typed.vcf.gz chr25_typed.vcf.gz chr26_typed.vcf.gz chr27_typed.vcf.gz chr28_typed.vcf.gz chr30_typed.vcf.gz  chr31_typed.vcf.gz chr32_typed.vcf.gz  chr33_typed.vcf.gz -o  chicken_typed.vcf.gz  --output-type z "
```
###后面就可以用我们自己的服务器了
- 为数据建立索引
```
tabix -p vcf chicken_typed.vcf.gz
```
- 选出SNP和indel（建议分染色体）
```
nohup gatk SelectVariants -R /public/jychu/refs/Gallus_gallus.GRCg6a.dna.toplevel.fa  -select-type SNP  -V chicken_typed.vcf.gz -O chicken_typed.snp.vcf.gz &
nohup gatk SelectVariants -R /public/jychu/refs/Gallus_gallus.GRCg6a.dna.toplevel.fa -select-type INDEL -V chicken_typed.vcf.gz -O chicken_typed.indel.vcf.gz 
```
- 对snp作硬过滤
```
gatk VariantFiltration  -R /public/jychu/refs/Gallus_gallus.GRCg6a.dna.toplevel.fa  -V chicken_typed.snp.vcf.gz -cluster 3 -window 10 --filter-expression " QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"  --filter-name "my_snp_filters"  -O chicken_typed.snp.filter.vcf.gz
#10bp以内的variants的数量不应当超过3个，如果超过了则会用SnpCluster标示出来
时间太长，改为分染色体对snp作硬过滤
mkdir chr_hardfilter
nohup bcftools filter -R chrMT.region chicken_typed.snp.vcf.gz --output chr_hardfilter/chrMT_typed.snp.vcf.gz --output-type z &
nohup bcftools filter -R chrW.region chicken_typed.snp.vcf.gz --output chr_hardfilter/chrW_typed.snp.vcf.gz --output-type z &
nohup bcftools filter -R chrZ.region chicken_typed.snp.vcf.gz --output chr_hardfilter/chrZ_typed.snp.vcf.gz --output-type z &
.........
nohup bcftools filter -R chr33.region chicken_typed.snp.vcf.gz --output chr_hardfilter/chr33_typed.snp.vcf.gz --output-type z &
tabix -p vcf chr*_typed.snp.vcf.gz              #同理建立索引
nohup gatk VariantFiltration  -R /public/jychu/refs/Gallus_gallus.GRCg6a.dna.toplevel.fa  -V chrMT_typed.snp.vcf.gz -cluster 3 -window 10 --filter-expression " QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"  --filter-name "my_snp_filters"  -O  chrMT_typed.snp.filter.vcf.gz &
......
nohup gatk VariantFiltration  -R /public/jychu/refs/Gallus_gallus.GRCg6a.dna.toplevel.fa  -V chr33_typed.snp.vcf.gz -cluster 3 -window 10 --filter-expression " QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"  --filter-name "my_snp_filters"  -O chr33_typed.snp.filter.vcf.gz &
```
- 取每个染色体PASS的位点
```
awk '/^#/||$7=="PASS"' chr10_typed.snp.filter.vcf> chr10_typed.snp.filter.PASS.vcf
##去掉没用的注释信息
sed -i '/AADN/d' chrMT_typed.snp.filter.PASS.vcf
```
- **合并常染色体**（bgzip压缩文件）
```
bgzip chr1_typed.snp.filter.PASS.vcf
tabix -p vcf chr1_typed.snp.filter.PASS.vcf.gz
bcftools concat chr1_typed.snp.filter.PASS.vcf.gz chr2_typed.snp.filter.PASS.vcf.gz chr3_typed.snp.filter.PASS.vcf.gz chr4_typed.snp.filter.PASS.vcf.gz chr5_typed.snp.filter.PASS.vcf.gz chr6_typed.snp.filter.PASS.vcf.gz chr7_typed.snp.filter.PASS.vcf.gz chr8_typed.snp.filter.PASS.vcf.gz chr9_typed.snp.filter.PASS.vcf.gz chr10_typed.snp.filter.PASS.vcf.gz chr11_typed.snp.filter.PASS.vcf.gz chr12_typed.snp.filter.PASS.vcf.gz chr13_typed.snp.filter.PASS.vcf.gz chr14_typed.snp.filter.PASS.vcf.gz chr15_typed.snp.filter.PASS.vcf.gz chr16_typed.snp.filter.PASS.vcf.gz chr17_typed.snp.filter.PASS.vcf.gz chr18_typed.snp.filter.PASS.vcf.gz chr19_typed.snp.filter.PASS.vcf.gz chr20_typed.snp.filter.PASS.vcf.gz chr21_typed.snp.filter.PASS.vcf.gz chr22_typed.snp.filter.PASS.vcf.gz chr23_typed.snp.filter.PASS.vcf.gz chr24_typed.snp.filter.PASS.vcf.gz chr25_typed.snp.filter.PASS.vcf.gz chr26_typed.snp.filter.PASS.vcf.gz chr27_typed.snp.filter.PASS.vcf.gz chr28_typed.snp.filter.PASS.vcf.gz chr30_typed.snp.filter.PASS.vcf.gz chr31_typed.snp.filter.PASS.vcf.gz chr32_typed.snp.filter.PASS.vcf.gz chr33_typed.snp.filter.PASS.vcf.gz  -o chicken_PASS.vcf.gz --output-type z
```
- 查看有多少变异
```
bcftools view -H chicken_PASS.vcf.gz | wc -l
```
- vcftools 质控（ 最大等位基因数小于等于2，最小测序深度大于等于5）
```
vcftools --gzvcf  chicken_PASS.vcf.gz --min-alleles 2 --max-alleles 2 --min-meanDP 5  --recode --stdout | bgzip -c >chicken.487.vcf.gz
```
- 对indel做硬过滤
```
tabix -p vcf chicken_typed.indel.vcf.gz 
gatk VariantFiltration -R /public/jychu/refs/Gallus_gallus.GRCg6a.dna.toplevel.fa -V chicken_typed.indel.vcf.gz  --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "my_indel_filters"  -O chicken477.indel.typed.filter.vcf
```
- 选出通过的indel
```
awk '/^#/||$7=="PASS"' chicken477.indel.typed.filter.vcf> chicken477.indel.typed.pass.vcf
```
- 对indel进行等位基因和覆盖深度的质控
```
bgzip chicken477.indel.typed.pass.vcf
mkdir indelvcf
vcftools --gzvcf  chicken477.indel.typed.pass.vcf.gz --max-alleles 2 --min-meanDP 5  --recode --stdout | bgzip -c >indelvcf/chicken.477.indel.vcf.gz
```
