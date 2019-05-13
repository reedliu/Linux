# One-liners for Bioinformatics

> Reed Liu created on 2019.1.11
>
> This is my collection of useful linux one-liners in daily work
> <https://github.com/stephenturner/oneliners>

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

##### find bam in the current folder recursively and copy them to a new dir with 5 CPUs => 拷贝大文件（如bam）到其他文件夹，并用5个线程

```shell
find . -name "*bam" | xargs -P5 -I{} rsync -av {} dest_dir
```

##### group files by extensions => 按照后缀的顺序排序文件

```shell
ll -X
```

##### loop through all the names => 循环语句

```shell
for i in {1..22} X Y 
do
  echo $i
done
# 对于{01..22} 的结果是 01 02 ...
```



### GREP

##### grep fastq reads containing a pattern but maintain the fastq format => 匹配fq中序列并打印

```shell
# 例如要在SP1.fq中找到这段序列的fq格式
# 如果匹配到多个，那么每条序列中间会用--分隔，因此需要用sed去除
$ grep -A 2 -B 1 'TGAGACAACATCT' SP1.fq | sed '/^--$/d' > out.fq
```

### SED 

##### delete with sed => 删除行

```shell
# delete blank lines
sed /^$/d
# delete the last line
sed $d
```

### AWK 

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

# 结果
chr1	123	aa	b	c	xxxx	abcd
chr1	234	a	b	c	
chr1	345	aa	b	c	yyyy	defg
chr1	456	a	b	c	
```

> Explanation:
>
> - **NR==FNR** NR is the current input line number and FNR the current file's line number. The two will be equal only while the 1st file is being read.
> - **OFS**  awk set the output field seperator; while set the input seperator is **-F** 
> - **next** means to proceed for the next line, rather than execute the following { } code block

##### awk to compare two different files and print if matches=> 比较两个文件的指定列，然后打印比对上的行

```shell
# https://unix.stackexchange.com/questions/134829/compare-two-columns-of-different-files-and-print-if-it-matches
# 例如 file1
abc|123|BNY|apple|
cab|234|cyx|orange|
def|kumar|pki|bird|
# file2
abc|123|
kumar|pki|
cab|234
# expected
abc|123|BNY|apple|
cab|234|cyx|orange|

$  awk -F'|' 'NR==FNR{a[$1$2]++;next};a[$1$2] > 0' file2 file1
```

##### conditional operator => 条件判断

基本格式：`var=condition?condition_if_true:condition_if_false`

例如：

```shell
# 现在有这个文件 test
a1	ACTGTCTGTCACTGTGTTGTGATGTTG
a2	ACTTTATATAT
a3	ACTTATATATATATA
a4	ACTTATATATATATA
a5	ACTTTATATATT	
# 我想看看每行序列部分是不是大于14个碱基
$ awk '{print (length($2)>14)?$0">14":$0"<=14"}' test
```

##### get new line => 在原来的内容基础上增加新内容

```shell
$ awk 'BEGIN{while((getline k <"test")>0) print "NEW:"k}{print}' test
```

##### merge multi-fasta into one single fasta => 合并多个fasta文件到一个文件中

```shell
# give a awk script called linearize.awk
$cat >linearize.awk 
# then copy and paste below
/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}

# run the awk script
$ paste <(awk -f linearize.awk file1.fa ) <(awk -f linearize.awk file2.fa  )| tr "\t" "\n" > multi.fa
```

##### 根据id输出序列

```shell
while read -r line; do awk -v pattern=$line -v RS=">" '$0 ~ pattern { printf(">%s", $0); }'  Seq.fasta; done < id.txt > output.fa 
```

需要指定Seq.fasta、id.txt