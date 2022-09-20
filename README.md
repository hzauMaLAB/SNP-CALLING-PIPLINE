# SNP-CALLING-PIPLINE （gatk4.1.9版本）

**预先安装的软件有`gatk4.1.9`、`bwa`、`samtools`、`bcftools`**
- 先创建储存文件的目录
```
cd ~
mkdir project    #创建一个项目
cd project
mkdir breed1   #创建要进行snpcalling的品种的文件夹
cd breed1     #进入breed1的数据目录
mkdir -p  bam/markdup/SNP_calling  #创建文件夹
```
- 给参考基因组以及dbsnp数据建立索引
```
#下载参考基因组并建立索引
mkdir -p ~/ref/bwa/pig
cd ~/ref/bwa/pig
wget http://ftp.ensembl.org/pub/release-106/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
gzip -d Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz  #解压
samtools faidx  Sus_scrofa.Sscrofa11.1.dna.toplevel.fa   #建立索引
wget  http://ftp.ensembl.org/pub/release-106/variation/vcf/sus_scrofa/sus_scrofa.vcf.gz  #下载dbsnp数据文件
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
cd cd ~/ref/bwa/pig
bwa index -p Sus_scrofa Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
```
- 运行单个样本callsnp
```
cd ~/project/breed1/cleandata  #进入存放质控过的数据目录\
vim ~/project/SNP_CALLtest01.sh    #创建一个最底层callsnp脚本
```
- 运行
```
for i in `ls *_1_*|tr "_" "\t"|cut -f1`;do sed "s/test01/$i/g" ../../SNP_CALLtest01.sh >SNP_CALL$i.sh;done    #批量替换脚本名
for i in `ls *_1_*|tr "_" "\t"|cut -f1`;do qsub -V -cwd -q all.q,fat.q -l vf=6G -pe smp 3 -S /bin/bash  SNP_CALL$i.sh;done  #qsub提交系统批量提交任务
```
- gvcf合并(要等所有样本的gvcf都出来之后，可以分群体也可以一块合并，如果分群体最后也是要一块合并的)
```
mkdir -p ~/project/allsample
find *.gz > input.list
nohup java -Xmx50g -jar ~/softwore/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar CombineGVCFs -R ~/ref/bwa/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -V  input.list  -O ~~/project/allsample/ALLSAMPLE.vcf.gz &
```
- gvcf合并后提取猪的不同染色体的vcf（分染色体）  可以用vcftools和bcftools提取，bcftools提取的时候会和原来的文件一模一样，vcftools提取的时候info那一列好像会改变一些信息，这个问题我没有仔细去研究
- 1.vcftools提取
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
- 2.bcftools提取
```
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 X Y MT;do echo "$i 1 274330532"|tr " " "\t" > chr$i.region;done      
bsub -J chrMT -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chrMT.region combine.vcf.gz --output bcftools/chrMT.vcf.gz --output-type z"
bsub -J chrW -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chrX.region combine.vcf.gz --output bcftools/chrX.vcf.gz --output-type z"
bsub -J chrZ -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chrY.region combine.vcf.gz --output bcftools/chrY.vcf.gz --output-type z"
bsub -J chr1 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr1.region combine.vcf.gz --output bcftools/chr1.vcf.gz --output-type z"
bsub -J chr2 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr2.region combine.vcf.gz --output bcftools/chr2.vcf.gz --output-type z"
bsub -J chr3 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr3.region combine.vcf.gz --output bcftools/chr3.vcf.gz --output-type z"
bsub -J chr4 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr4.region combine.vcf.gz --output bcftools/chr4.vcf.gz --output-type z"
bsub -J chr5 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr5.region combine.vcf.gz --output bcftools/chr5.vcf.gz --output-type z"
bsub -J chr6 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr6.region combine.vcf.gz --output bcftools/chr6.vcf.gz --output-type z"
bsub -J chr7 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr7.region combine.vcf.gz --output bcftools/chr7.vcf.gz --output-type z"
bsub -J chr8 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr8.region combine.vcf.gz --output bcftools/chr8.vcf.gz --output-type z"
bsub -J chr9 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr9.region combine.vcf.gz --output bcftools/chr9.vcf.gz --output-type z"
bsub -J chr10 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr10.region combine.vcf.gz --output bcftools/chr10.vcf.gz --output-type z"
bsub -J chr11 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr11.region combine.vcf.gz --output bcftools/chr11.vcf.gz --output-type z"
bsub -J chr12 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr12.region combine.vcf.gz --output bcftools/chr12.vcf.gz --output-type z"
bsub -J chr13 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr13.region combine.vcf.gz --output bcftools/chr13.vcf.gz --output-type z"
bsub -J chr14 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr14.region combine.vcf.gz --output bcftools/chr14.vcf.gz --output-type z"
bsub -J chr15 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr15.region combine.vcf.gz --output bcftools/chr15.vcf.gz --output-type z"
bsub -J chr16 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr16.region combine.vcf.gz --output bcftools/chr16.vcf.gz --output-type z"
bsub -J chr17 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr17.region combine.vcf.gz --output bcftools/chr17.vcf.gz --output-type z"
bsub -J chr18 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools filter -R chr18.region combine.vcf.gz --output bcftools/chr18.vcf.gz --output-type z"
```
- 对提取出来的vcf文件建立索引
```
tabix -p vcf chr*.region combine.vcf.gz       #这样概括所有的染色体，这个命令只能一次操作一个文件
```
- 添加注释并genotype
```
wget http://ftp.ensembl.org/pub/release-106/variation/vcf/sus_scrofa/sus_scrofa.vcf.gz
bsub -J chr10 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "tabix -p vcf sus_scrofa.vcf.gz"
bsub -J chr10 -n 8 -R span[hosts=1] -o %J.out -e %J.err -q normal "java -Xmx80g -jar ~/soft/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar  GenotypeGVCFs -R ~/ref/bwa/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -V chr10.vcf.gz --dbsnp sus_scrofa.vcf.gz -O  chr10_typed.vcf.gz"
#chr10_typed.vcf.gz 如果没有索引，也要建立索引
```
- 纵向合并vcf染色体文件(合并的文件必须是bgzip压缩过的文件)   这一步也可以先不合并，最后过滤完snp位点的时候再合并
```
bsub -J chr33 -n 1 -R span[hosts=1] -o %J.out -e %J.err -q normal "/public/home/hsong/jychu/soft/bcftools-1.8/bcftools concat chrMT_typed.vcf.gz  chrX_typed.vcf.gz   chrY_typed.vcf.gz   chr1_typed.vcf.gz   chr2_typed.vcf.gz  chr3_typed.vcf.gz chr4_typed.vcf.gz  chr5_typed.vcf.gz  chr6_typed.vcf.gz  chr7_typed.vcf.gz  chr8_typed.vcf.gz  chr9_typed.vcf.gz  chr10_typed.vcf.gz  chr11_typed.vcf.gz chr12_typed.vcf.gz   chr13_typed.vcf.gz  chr14_typed.vcf.gz  chr15_typed.vcf.gz  chr16_typed.vcf.gz  chr17_typed.vcf.gz chr18_typed.vcf.gz  -o  pig_typed.vcf.gz  -Oz "
```
- 为数据建立索引
```
tabix -p vcf pig_typed.vcf.gz
```
- 选出SNP和indel（建议分染色体）
```
nohup gatk SelectVariants -R ~/ref/bwa/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa  -select-type SNP  -V pig_typed.vcf.gz -O pig_typed.snp.vcf.gz &
nohup gatk SelectVariants -R ~/ref/bwa/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -select-type INDEL -V pig_typed.vcf.gz -O pig_typed.indel.vcf.gz 
```
- 对snp作硬过滤
```
gatk VariantFiltration  -R ~/ref/bwa/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa  -V pig_typed.snp.vcf.gz -cluster 3 -window 10 --filter-expression " QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"  --filter-name "my_snp_filters"  -O pig_typed.snp.filter.vcf.gz
#10bp以内的variants的数量不应当超过3个，如果超过了则会用SnpCluster标示出来
时间太长，改为分染色体对snp作硬过滤
mkdir chr_hardfilter
nohup bcftools filter -R chrMT.region pig_typed.snp.vcf.gz --output chr_hardfilter/chrMT_typed.snp.vcf.gz --output-type z &
.........
nohup bcftools filter -R chr18.region pig_typed.snp.vcf.gz --output chr_hardfilter/chr18_typed.snp.vcf.gz --output-type z &
tabix -p vcf chr*_typed.snp.vcf.gz              #同理建立索引
nohup gatk VariantFiltration  -R ~/ref/bwa/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa  -V chrMT_typed.snp.vcf.gz -cluster 3 -window 10 --filter-expression " QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"  --filter-name "my_snp_filters"  -O  chrMT_typed.snp.filter.vcf.gz &
......
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
..........
bcftools concat chr1_typed.snp.filter.PASS.vcf.gz chr2_typed.snp.filter.PASS.vcf.gz chr3_typed.snp.filter.PASS.vcf.gz chr4_typed.snp.filter.PASS.vcf.gz chr5_typed.snp.filter.PASS.vcf.gz chr6_typed.snp.filter.PASS.vcf.gz chr7_typed.snp.filter.PASS.vcf.gz chr8_typed.snp.filter.PASS.vcf.gz chr9_typed.snp.filter.PASS.vcf.gz chr10_typed.snp.filter.PASS.vcf.gz chr11_typed.snp.filter.PASS.vcf.gz chr12_typed.snp.filter.PASS.vcf.gz chr13_typed.snp.filter.PASS.vcf.gz chr14_typed.snp.filter.PASS.vcf.gz chr15_typed.snp.filter.PASS.vcf.gz chr16_typed.snp.filter.PASS.vcf.gz chr17_typed.snp.filter.PASS.vcf.gz chr18_typed.snp.filter.PASS.vcf.gz  -o pig_PASS.vcf.gz --output-type z
```
- 查看有多少变异
```
bcftools view -H pig_PASS.vcf.gz | wc -l
```
- vcftools 质控（ 最大等位基因数小于等于2，最小测序深度大于等于5）
```
vcftools --gzvcf  pig_PASS.vcf.gz --min-alleles 2 --max-alleles 2 --min-meanDP 5  --recode --stdout | bgzip -c >pig.vcf.gz
```
- 对indel做硬过滤
```
tabix -p vcf pig_typed.indel.vcf.gz 
gatk VariantFiltration -R ~/ref/bwa/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa -V pig_typed.indel.vcf.gz  --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "my_indel_filters"  -O pig.indel.typed.filter.vcf
```
- 选出通过的indel
```
awk '/^#/||$7=="PASS"' pig.indel.typed.filter.vcf> pig.indel.typed.pass.vcf
```
- 对indel进行等位基因和覆盖深度的质控
```
bgzip pig.indel.typed.pass.vcf
mkdir indelvcf
vcftools --gzvcf  pig.indel.typed.pass.vcf.gz --max-alleles 2 --min-meanDP 5  --recode --stdout | bgzip -c >indelvcf/pig.indel.vcf.gz
```
