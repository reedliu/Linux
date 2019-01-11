# One-liners for Bioinformatics

> Reed Liu created on 2019.1.11
>
> This is my collection of useful linux one-liners in daily work

### About Fastq/fasta

##### fastq sequences length distribution => 得到fq文件中序列长度的分布

```shell
$ zcat file.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'  
```

##### reverse complement  => 反向互补

```shell
echo 'ATTGCTATGCTNNNT' | rev | tr 'ACTG' 'TGAC'
```

##### fastq2fasta 

```shell
zcat file.fastq.gz | paste - - - - | perl -ane 'print ">$F[0]\n$F[1]\n";' | gzip -c > file.fasta.gz
```

##### split a multifasta file into single ones with csplit => fasta按>拆分

```shell
# * refers to the number of files 可以选择拆分的文件数量
$ csplit -z -q -n 4 -f test test.fa /\>/ {*}
# OR use awk 一次性全部按>拆分
awk '/^>/{s=++d".fa"} {print > s}' multi.fa
```

##### single line fasta to multi-line of 50 characters in each line => 单行fa变多行

```shell
awk -v FS= '/^>/{print;next}{for (i=0;i<=NF/50;i++) {for (j=1;j<=50;j++) printf "%s", $(i*50 +j); print ""}}' file

# fold -w 50 file
```

##### multi-line fasta to one-line => 一个多行fa文件变单行

```shell
# 方法一：
awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' file.fa
# 方法二：
cat file.fasta | awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}'
```



##### Number of reads in a fastq file => 统计fq中序列数（4行一个序列）

```shell
cat file.fq | echo $((`wc -l`/4))
```



### About Sam/Bam

##### bam2bed

```shell
samtools view file.bam | perl -F'\t' -ane '$strand=($F[1]&16)?"-":"+";$length=1;$tmp=$F[5];$tmp =~ s/(\d+)[MD]/$length+=$1/eg;print "$F[2]\t$F[3]\t".($F[3]+$length)."\t$F[0]\t0\t$strand\n";' > file.bed
```

##### bam2wig

```shell
samtools mpileup -BQ0 file.sorted.bam | perl -pe '($c, $start, undef, $depth) = split;if ($c ne $lastC || $start != $lastStart+1) {print "fixedStep chrom=$c start=$start step=1 span=1\n";}$_ = $depth."\n";($lastC, $lastStart) = ($c, $start);' | gzip -c > file.wig.gz
```

##### 