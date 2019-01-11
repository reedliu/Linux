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
$ echo 'ATTGCTATGCTNNNT' | rev | tr 'ACTG' 'TGAC'
```

##### fastq2fasta 

```shell
$ zcat file.fastq.gz | paste - - - - | perl -ane 'print ">$F[0]\n$F[1]\n";' | gzip -c > file.fasta.gz
```

##### split a multifasta file into single ones with csplit => fasta按>拆分

```shell
# * refers to the number of files 可以选择拆分的文件数量
$ csplit -z -q -n 4 -f test test.fa /\>/ {*}
# OR use awk 一次性全部按>拆分
$ awk '/^>/{s=++d".fa"} {print > s}' multi.fa
```

##### single line fasta to multi-line of 50 characters in each line => 单行fa变多行

```shell
$ awk -v FS= '/^>/{print;next}{for (i=0;i<=NF/50;i++) {for (j=1;j<=50;j++) printf "%s", $(i*50 +j); print ""}}' file

# fold -w 50 file
```

##### multi-line fasta to one-line => 一个多行fa文件变单行

```shell
# 方法一：
$ awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' file.fa
# 方法二：
$ cat file.fasta | awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}'
```

##### Number of reads in a fastq file => 统计fq中序列数（4行一个序列）

```shell
$ cat file.fq | echo $((`wc -l`/4))
```

##### print length of each entry in a multifasta file => 打印fa中每个序列的长度

```shell
$ awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' file.fa
```

##### subsample fastq => 取fq文件的子集（其中0.01是指取出来百分之1的reads）

```shell
$ cat file.fq | paste - - - - | awk 'BEGIN{srand(1234)}{if(rand() < 0.01) print $0}' | tr '\t' '\n' > out.fq
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





### Basis Linux

##### get all folders' size in the current folder => 当前目录下的所有目录大小

```shell
$ du -h --max-depth=1
```

##### exit a dead ssh session => 退出卡死的ssh界面

```shell
$ ~.
```

##### copy large folders fast => 快速拷贝大文件夹

```shell
# copy every file in folder 拷贝目录下的所有文件
rsync -av from_dir/ to_dir
# skip transferred files 跳过已拷贝的文件
rsync -avhP from_dir/ to_dir
```

##### delete with sed

```shell
# delete blank lines
sed /^$/d
# delete the last line
sed $d
```

##### awk join two files with common columns => awk连接有共同列的文件（类似于R的merge函数）

```shell
# http://stackoverflow.com/questions/13258604/join-two-files-using-awk
# file_a.bed： 
chr1	123	aa	b	c	d
chr1	234	a	b	c	d
chr1	345	aa	b	c	d
chr1	456	a	b	c	d
# file_b.bed
xxxx	abcd	chr1	123	aa	c	d	e
yyyy	defg	chr1	345	aa	e	f	g
# 现在想在a的基础上根据a、b共有列来增加b中的新内容
$ awk 'NR==FNR{a[$3,$4,$5]=$1OFS$2;next}{$6=a[$1,$2,$3];print}' OFS='\t' \
file_b.bed file_a.bed
```

