# FASTQ/A Linux Exercises

> Reed Liu created on 2019.1.11

### 下载数据

```shell
$ wget http://molb7621.github.io/workshop/_downloads/SP1.fq
```

---

### 练习一  统计行数

##### 统计总行数

```shell
$ wc -l SP1.fq 
1000
```

#####  统计fq序列数(每四行一个fq序列)

```shell
$ wc -l Sp1.fq | awk '{print $1/4}'
250
```

##### 错误示范：利用grep找开头字符`@`

```shell
$ grep -c "@" SP1.fq
459
```

> 因为`@`符号除了开头，还会在第四行Phred质量值中出现

---

### 练习二 提取fa序列 

> 因为fq文件中的第二行是fa信息，所以想办法提取出每个fq的第二行

##### 首先我们需要得到行号信息

```shell
$ awk '{print NR}' SP1.fq
#会输出每一行的行号，共1000行
```

##### 然后利用取余/求模符号`%` 

怎么用？例如：取出1-4行，可以这样

```shell
$ awk 'BEGIN {print 1%4}'
$ awk 'BEGIN {print 2%4}'
$ awk 'BEGIN {print 3%4}'
$ awk 'BEGIN {print 4%4}'

#关于BEGIN的作用
#  In the case of BEGIN and END blocks, awk will process the statements only once. You can use the BEGIN block to initialise variables or other routines which only need to be performed once but it can also be used to run Awk commands when there are no files to process. 
```

这样我们就可以对每条fq进行编号，从1-4（这里4用0表示，因为余数是0）

```shell
$ awk '{print NR % 4, $0}' SP1.fq
1 @cluster_807:UMI_GATATG
2 TTTTTCCACACGTAAAATTTATAAACATTTA
3 +
0 83EEA3:620725;DD5CAB:53C3=/472-
1 @cluster_809:UMI_CTTTTA
2 CACAAGGAATATCATTTTATTACTGTAATCA
3 +
0 ?=?A?A>@AC??A@EEEC?EC=?C@C@A?A=
```

##### 这样的话，其中这4行，想要哪行选哪行

```shell
# 如果只是判断条件的话，不需要大括号。大括号的目的是运行命令 action
$ awk 'NR % 4 == 1' SP1.fq | head  # get header line
$ awk 'NR % 4 == 2' SP1.fq | head  # get sequence line
$ awk 'NR % 4 == 3' SP1.fq | head  # get comment line
$ awk 'NR % 4 == 0' SP1.fq | head  # get quality line
```

>  好了，每个fq的第二行选出来就是fa格式了

##### 如果想看看每个fq中的fa序列出现的次数

```shell
$  awk 'NR % 4 ==2' SP1.fq | sort | uniq -c  | wc -l
#意思就是先全部选出来，然后进行A-Z排序，统计排序后两条一致的数目
240
# 有趣的是，只sort不uniq结果是正常的250，因此，有10条序列是重复出现的
```

##### 再看看都是哪些序列出现了重复

```shell
$ awk 'NR % 4 ==2 ' SP1.fq | sort | uniq -c | sort -k1,1nr | head
3 CCCCCCCCCAAATCGGAAAAACACACCCCTA
3 TCCCCCCCCCAAATCGGAAAAACACACCCCC
2 CAGCTTTGCAACCATACTCCCCCCGGAACCC
2 CCCCCCCAGATCGGAAAAGCACACGCCTGAA
2 GCCCCCCCCCAAATCGGAAAAACACACCCCC
2 GGCTTTGCAACCATACTCCCCCCGGAACCCA
2 GGTTGAGCACAGGGTACTTTATTGATGGTAC
2 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
```

##### 如果选1-10bp序列，那么又是哪些重复呢？和上面统计的结果还是同一条吗？

```shell
$ awk 'NR % 4 ==2' SP1.fq | cut -c 1-10 | sort | uniq -c | sort -k1,1nr | head
4 CCCCCCCAGA
4 CCCCCCCCCA
4 GTTTTTTTTT
4 TCCCCCCCCC
4 TTTTTTTTTT
3 CAGCTTTGCA
2 AGCTTTGCAA
2 CCCCCCCCAA
2 CTTTTTTTTT
2 GAGACAGAGT
# 发现取出来10bp，统计的重复结果和对整条fa序列的统计结果不同
```

##### 既然成功提取了第二行的fa序列，那么进阶一下，提取第一行，并提取出其中的cluster信息

```shell
# 原来是这样：@cluster_2:UMI_ATTCCG
# 我们想要这样：@cluster_2
# 也就是提取第一行的1-10个字符
# 第一种：利用substr搭配awk，其中$0表示一整行
# substr($field_index, start, end)
$ awk 'NR % 4 ==1 {print substr($0, 1,10)}' SP1.fq | head
# 第二种：或者通过观察发现，我们要提取的正好是:分隔的第一部分，所有还可以这样：
$ awk -F: 'NR % 4 ==1 {print $1}' SP1.fq | head
```

同样，也可以提取UMI（Unique Molecular Identifyier）信息

```shell
# 第三种：当然，这也是提取特定字符的第三种方法
$ awk 'NR % 4 ==1' SP1.fq | cut -d ':' -f 2 | head -3
UMI_ATTCCG
UMI_CTTTGA
UMI_GGTCAA
```

> 比较以上三种方法，可以发现：substr使用最灵活，可以按字符提取，其他两种需要指定分隔符

**进阶：**利用`substr`，就可以整出来一个小的fq文件（比如trim掉前10bp的序列，注意质量序列也需要trim 前10个字符）

```shell
$ awk 'NR % 4 ==1 {print $0}; NR % 4==2 {print substr($0, 10,length($0))}; NR % 4 ==3 {print $0}; NR % 4 ==0 {print substr($0, 10,length($0))}' SP1.fq | head
@cluster_2:UMI_ATTCCG
CACATAATCTTCAGCCGGGCGC
+
4868>9:67AA<9>65<=>591
@cluster_8:UMI_CTTTGA
AATACTCTCCGAACGGGAGAGC
+
003,-2-22+00-12./.-.4-
@cluster_12:UMI_GGTCAA
GATCATTTTATTGAAGAGCAAG
```

**进阶：**利用awk将fq转为fa

```shell
# 如果只是加一个">"，那么还带着fastq标识符@是不对的
$ awk 'NR %4==1 {print ">"$1}; NR%4==2 {print}' SP1.fq | head
# 因此可以用substr来取
$ awk 'NR %4==1 {print ">"substr($0,2,length($0))}; NR%4==2 {print}' SP1.fq | head -4
>cluster_2:UMI_ATTCCG
TTTCCGGGGCACATAATCTTCAGCCGGGCGC
>cluster_8:UMI_CTTTGA
TATCCTTGCAATACTCTCCGAACGGGAGAGC
```

> 当转成fasta后，就可以用`grep -c ">"` 来统计序列的数目了

**进阶：**从fasta中找到至少四个A的序列，并显示是序列名称

```shell
# 首先要构建一个名称+fa序列的子文件，然后再查找，注意grep的-B表示打印匹配行的前一行；-A表示打印匹配行的后一行
$ awk 'NR % 4 ==1 {print $0}; NR % 4 ==2 {print $0}' SP1.fq | egrep "A{4,}" -B1
```

上面的结果中的四个A都是序列中的吗？会不会在名称中也抓取到了`AAAA`？
因此想要看看是否在名称中出现四个A：

```shell
# 利用awk实现
$ awk 'NR % 4 ==1 {print $0}; NR % 4 ==2 {print $0}' SP1.fq | awk '/UMI/ && /AAAA/'
# 利用grep实现
$ awk 'NR % 4 ==1 {print $0}; NR % 4 ==2 {print $0}' SP1.fq | egrep "UMI.*AAAA"

@cluster_252:UMI_CAAAAG
#还真的有一个！
```





**参考：**

关于BEGIN的作用

https://unix.stackexchange.com/questions/119907/begin-and-end-with-the-awk-command

关于awk

https://likegeeks.com/awk-command/

关于grep打印上下行

https://stackoverflow.com/questions/1072643/how-can-i-make-grep-print-the-lines-below-and-above-each-matching-line