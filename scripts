https://www.jianshu.com/p/a84cd44bac67
# 进行质控
nohup fastqc -t 4 ./*.gz -o ./fastqc_results > ./fastqc.log 2>&1 &
# 去除低质量的部分
nohup trim_galore --illumina -o clean raw/fastq/*.gz &
nohup trim_galore -q 25 --phred33 --fastqc --length 36 --stringency 3 -o clean raw/fastq/*.gz &
平台：--nextera
# 序列比对
index=/home/duansq/reference/index/hisat/mm10
ls ./*.gz | while read id;do
  hisat2 -p 5 -x $index -U $id -S $(basename $id).bam.hisat.sam
done
# sam2bam
ls *.sam|while read id;do (samtools sort -O bam -@ 5  -o $(basename ${id} ".sam").bam   ${id});done
# gtf
gtf=/home/duansq/reference/gtf/gencode/mm10/gencode.vM24.annotation.gtf.gz
featureCounts -T 5 -p -t exon -g gene_id -a $gtf -o all.id.txt *.bam 1>counts.id.log 2>&1 &






# mm10 gtf gencode
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.annotation.gtf.gz
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.chr_patch_hapl_scaff.annotation.gtf.gz
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.primary_assembly.annotation.gtf.gz
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.basic.annotation.gtf.gz
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.chr_patch_hapl_scaff.basic.annotation.gtf.gz
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.long_noncoding_RNAs.gtf.gz
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.polyAs.gtf.gz
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.2wayconspseudos.gtf.gz
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.tRNAs.gtf.gz
