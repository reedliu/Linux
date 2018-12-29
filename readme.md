# 生信Linux常用系统操作

> 整理的都是日常很好用但又不容易记住的
> 许多命令详情可以在[命令大全] (http://man.linuxde.net) 中搜索

#### 关于linux系统 

```shell
lsb_release -a #然后得到详细系统版本
LSB Version:	:core-4.1-amd64:core-4.1-noarch
Distributor ID:	CentOS
Description:	CentOS Linux release 7.4.1708 (Core) 
Release:	7.4.1708
Codename:	Core
head -n 1 /etc/issue   # 查看操作系统版本
```

```shell
uname -a # 得到linux位数
Linux bioinfoplanet 3.10.0-693.2.2.el7.x86_64 #1 SMP Tue Sep 12 22:26:13 UTC 2017 x86_64 x86_64 x86_64 GNU/Linux
free -g # 查看内存
df -h # 查看硬盘总额
du -sh # 查看当前用户用量，或者加特定文件（夹）查看大小
uname -a               # 查看内核/操作系统/CPU信息
cat /proc/cpuinfo| grep process | wc -l  # 查看CPU个数
env                    # 查看环境变量
```

```shell
# linux上的任务管理器
ps -ef
# 如果想查找某个任务（比如找“包含sra”的）
ps -ef | grep "sra"
```

```shell
#批量kill进程
ps -ef | grep "sra" | awk '{print $2}' | while read id;do kill $id; done
```

```shell
#basename命令格式：
basename [pathname] [suffix]
basename [string] [suffix]
 
#suffix为后缀，如果suffix被指定了，basename会将pathname或string中的suffix去掉

$ basename /tmp/test/file.txt
file.txt
$ basename /tmp/test/file.txt .txt
file
```



#### 小命令：

`cd -` 两个路径来回切换 `-`代表前一个目录

`diff` 统计两个文件的不同

```shell
diff log2014.log log2013.log  -y -W 50
#结果返回
2013-02                 2013-02
2014-03               | 2013-03
2013-07                 2013-07
2013-07               | 2013-08
2013-11               <
2013-12               <
# "|"表示前后2个文件内容有不同
#"<"表示后面文件比前面文件少了1行内容
#">"表示后面文件比前面文件多了1行内容
```

`paste` 合并文件的列

```shell
$ cat a.txt     
bioinfo 520 
$ cat b.txt     
doudou huahua 
$ cat c.txt    
welcome hah 
#不加-s
paste a.txt b.txt c.txt 
bioinfo 520 	doudou huahau	welcome haha
#加上-s（serial 　串列进行而非平行处理）
paste -s a.txt b.txt c.txt
bioinfo 520 
doudou huahau
welcome haha
```



`wc` 统计行数、词数、字节数

> 普及一下字节与字符：
> 字节（Bytes）是计量单位，表示数据量大小
> 字符是文字和符号，如A、b、3、@等
> 不同的编码方式，他们对应的关系不同：
> **ASCII码**中，一个 英文字母（不分大小写）占一个字节，中文占两个字节
> **UTF-8编码**中，一个英文字母一个字节，一个中文（包括繁体）占三个字节
> **Unicode编码**中，一个英文两个字节，一个中文两个字节；英文标点占一个字节，中文标识占两个
> **UTF-16编码**中， 英文与中文都要两个字节
> **UTF-32编码**中，任何字符都要4个字节

####  进程管理：

* 查看进程：
  * 自带的命令： `ps aux` 或 `top`
  * 更好看更实用：`htop` （需要安装）
* 进程处理：
  * & 放到后台运行：`COMMAND &` 
  * fg 查看后台进程
  * 终止进程：两种方法
    * `fg %COMMAND` 先将程序放到前台，然后`ctrl + c`
    * htop 中 `F9`直接结束某一进程
  * 假如一开始忘记加&放到后台，后来想放
    可以先`ctrl+z` 暂停，然后 `bg` 放后台

#### 命令传递:

`$()` 先在括号中执行命令，然后返回结果给$前的命令

> 例如：新建当天日期的文件夹：**mkdir $(date +%F)** => 2018-06-06
>
> 这样就省去了先查看日期，再手动建立文件夹的操作了
> %F 显示完整日期年-月-日; 注意date 后面要跟 + 
> 【%D显示月-日-年】
>
> #如果只是mkdir $(date)，结果会将日期拆分成日期、星期、时间等好几个文件夹

或者用``也是一样的

#### 字符替换：

`tr命令`可以对来自标准输入的字符进行替换、压缩和删除

```
# ex.1 改变大小写
echo "HELLO WORLD" | tr 'A-Z' 'a-z'
hello world

# ex.2 -d删除，‘ ’匹配要删除的内容
echo "hello 123 world 456" | tr -d '0-9'
hello  world 

# ex.3 -c ‘0-9 \n’匹配除了数字、空格和换行符以外的字符，然后-d删除
echo aa.,a 1 b#$bb 2 c*/cc 3 ddd 4 | tr -d -c '0-9 \n'
 1  2  3  4
```

#### 获取路径文件名和目录名：

直接命令行输入`basename $(pwd) ` 就得到文件名
`dirname $(pwd)` 得到目录名

> 例如：pwd是/home/bioinfo/ncbi/blast.sh
> basename $(pwd)就得到blast.sh
> dirname $(pwd)就得到 /home/bioinfo/ncbis

#### 重命名：

> 想想吧，Windows电脑上批量处理文件名是一件多么难受的工作

`mv` 可以用，但是只能对单个文件命名，如果有多个，你需要rename

`rename 原字符串 目标字符串 要更改的文件` 就是由这三部分构成

> 支持通配符：原来有bio1, bio2, ..., bio233
> `？` 替代单个字符 rename bio bio0 bio? => 将 bio1-bio9变成bio01-bio09
> `*` 替代多个字符 rename bio bio0 bio* => 将 bio10-bio90变成bio010-bio099
> `[]` 替代[]中的任意单个字符 rename bio bio0 bio[2]* => 将bio200-bio233变成bio0200-bio0233

>  支持正则：
> 批量替换文件名：rename '"s/AA/aa/" * 将所有文件名中的AA替换成aa
> 批量替换后缀： rename "s//.html//.php" * 全部后缀换成.php
> 批量加后缀： rename "s/$//.txt" *
> 批量删除文件后缀： rename "s//.txt//" *

#### 通配符wildcard：

使用通配符的模式叫Globbing(因为最初被应用在一个叫/etc/glob的文件中)

> 1. 包括比如`* ? [] {}`
> 2. 删除目录中r2、d2开头的shell脚本
>    rm -fr {r2, d2}*.sh 
> 3. 拷贝文件到上一个目录【注意末尾的 - 】
>    cp /home/tmp/log[0-9].txt  -

#### 下载操作：

1. wget [option] URL

   ```
   -o 以新的文件名保存【wget默认会以最后一个符合/的后面的字符来命名】
   -c 断点续传
   -b 后台下载 【tail -f wget-log查看下载进程】
   --tries 数字 网络有问题或下载一个大文件也有可能失败，默认重试20次
   -i 同时下载多个文件
   	【首先新建一份下载链接cat > file.txt, 然后输入各个下载地址，再ctr+C退出；然后用wget -i file.txt就可以】
   
   #例如要下载ncbi的nt库
   wget -r -c -nH ftp://ftp.ncbi.nih.gov/blast/db/nt.*
   ```

2. curl :例如我要下载人类基因组22号染色体数据

   ```
   curl http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz | gunzip > chr22.fa
   ```

   关于curl，如果不加-o选项又不使用管道，他会自动输出到屏幕。-o指定输出的文件名

3. prefetch：需要先安装sratoolkit可以下载SRR文件。但是可能有时不稳定
   例如prefetch SRR519926，就自动下载到～/ncbi/public/sra中。

4. axel： 实测速度比wget快两倍，因为他是多线程下载，但是没断点续传选项，下载的数据完整性不如wget

5. ascp：大批量下载SRR/SRA/fq文件

   ```
   ascp -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l 200m 
   anonftp@ftp-private.ncbi.nlm.nih.gov:
   /sra/sra-instant/reads/ByRun/sra/SRR/SRR949/SRR949627/SRR949627.sra 
   # -v: 记录程序在干嘛
   # -i： 提供私钥文件的地址，地址一般是~/.aspera/connect/etc中的asperaweb_id_dsa.openssh文件
   # -k：断点续传，设为1
   # -T：取消加密，否则有时下载出错
   # -l： 最大传输速度，一般200m到500m
   # SRA在ascp的用户名是anonftp，数据的存放地址是ftp-private.ncbi.nlm.nih.gov
   # 数据存放一般都是/sra/sra-instant/reads/ByRun/sra/SRR/SRR***/SRR***###/SRR***###.sra
   ```

   ```
   # 想批量下载SRR173010-SRR173015
   for i in {10..15};do ascp -v -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l 200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR173/SRR1730"$i"/SRR1730"$i".sra ～/data;done
   ```


#### Shell:

**使用多个命令：** 放在一行，分号分隔

**脚本文件第一行：** 指定要用的shell `#!bin/bash`【用！告诉shell用哪个解释器；只用#是注释的意思】

**运行脚本注意：** 直接运行出错可能是没有权限，先改下脚本权限，一般新建完要`chmod u+x` ；shell会在环境变量PATH中查找脚本，可以在命令提示符后面用**绝对路径或者相对路径引用脚本**

**让脚本告诉用户他在干嘛：** 通过`echo` 加一个文本【如果文本中包含字符串，echo后要加双引号】

【**这里为何用双引号不用单引号？**事实上echo两个都支持，但是很大的不同就是单引号是强引用，它会把命令都忽略掉；一般我们想输出字符串就使用双引号即可】

**想让字符串和命令输出在同一行：** `echo -n` ，并且保证字符串末尾有一个空格

**引用变量：** 

- 环境变量：
  `env`显示全局变量， 比如HOME，UID等。要在shell中引用这些变量，在前面加一个美元符号，`$HOME` 、`$UID`

  【**实用：** 有时拿到的服务器中家目录设置可能不尽如人意，`cd ~` 经常跑错地方😱，如何自己修改家录？】
  `HOME="新的路径"`，`echo $HOME`看一下是不是改过来了。当然这也只是临时的。
  只有root用户才能永久修改用户名和家目录，`usermod -l  新用户名 旧用户名` `sudo usermod -d 新家目录 `

- 用户变量：
  在shell脚本中定义的变量，用等号赋值。**注意⚠️：变量、等号、值之间不能有空格** 
  【常常会因为美观问题忽视】

**一个脚本案例：** 
![image-20180714114055274](/Users/reedliu1/Library/Mobile Documents/com~apple~CloudDocs/232 - Past Writing | 历史文章，归档重查/My Writing/生信Linux常用操作/1.png)

**数学运算：** 一般shell脚本中，数字都是被当成字符处理的，所以数学运算需要特定的运算工具

- 方括号--整数运算	例如，`var1=$[2+5] && echo $var1 ` 结果返回7

- bc(bash caculator)浮点运算，支持设置小数点后精度

  ```
  #!/bin/bash
  var=$(echo "scale=4; 2.56/4" | bc)
  echo The result is: $var
  结果返回：The result is: .6400
  ```

  如果进行大量计算，使用<< EOF效果更好【支持多步运算】

  ```
  var=$(bc << EOF
  options【比如设置精度scale=2】
  statements【a1=($var1 * $var2); a2 = ($var2/$var3)】
  expression【运算 a1 * a2】
  EOF)
  ```

**脚本退出码 $?** 

这是个信号，表示脚本运行状态，0是成功，1-255是各种各样的问题

| 退出码 |       意思       |
| :----: | :--------------: |
|   0    |     运行成功     |
|   1    |     未知错误     |
|   2    | shell命令不适合  |
|  126   |   命令不可执行   |
|  127   |    没找到命令    |
|  128   |  无效的退出参数  |
|  130   |  通过ctr+C终止   |
|  255   | 正常范围外的退出 |

#### conda:

**基本配置：** `conda config` 
`--add channels` 添加国内镜像源，提高下载速度
`--show` 查看已有配置

**查看环境:**  `conda info --env` 

**新建环境：** `conda create -n xxx python=2 samtools` 
新建一个名为xxx，python环境为2/3的环境，并安装samtools的软件

**启动环境：** `conda activate xxx` **反激活环境** `conda deactivate xxx`

 **删除环境:**`conda remove -n xxx --all` 或者 `rm -rf /path to conda/xxx`  

**搜索软件：** `conda search xxx` 

**下载某个版本：** `conda install xxx=2.8.0 -y` 其中-y表示确认



#### 从头新建新用户：

1、 新建一个没有家目录的用户bioplanet;
​		--> 结果是该用户的属主、属组仍然为root

> useradd -M bioplanet

2、 复制/etc/skel 为 /home/bioplanet; 
​		—> 添加架构文件为了在/home下表现出来该用户

> cp  -r  /etc/skel   /home/bioplanet
> ls -la /home/bioplanet 看三个架构文件是否生成

3、改变/home/bioplanet 及其内部文件的属主属组均为bioplanet;

> chown -R bioplanet:bioplanet /home/bioplanet
> #ls -ld /home/bioplanet 验证

4、/home/bioplanet及其内部的文件，属组和其他用户没有任何访问权限

> 我们这里只需要改属组和其他，使用 chmod -R go=  /home/bioplanet

5、su bioplanet 新建用户完成, 接下来还要改三样东西
​	设置基本组为bioplanet（5000），附加组为之前存在的mygroup

6、修改passwd:  

> vi /etc/passwd
> 最后一行添加：`bioplanet:x:5000:5000::/home/bioplanet:/bin/bash`

7、修改group：

> vi /etc/group
> 最后一行添加：`bioplanet:x:5000:`, 另外找到mygroup，在后面添加上bioplanet

8、passwd [USERNAME]

#### VIM

1. 非编辑模式下，x删除后一个字符，X删除前一个字符，删除3个字符就是3X/x
2. 不退出vi编辑器查看当前目录的文件：`!pwd` 

#### 关于查看

- 不解压前提下查看`zless -SN *.gz`  【-S：加上分隔符tab；-N：加上行号】

#### linux的control快捷键

- control + 两次i：输入命令输到一半忘了文件夹下需要用到哪个文件，一般操作是返回先去查看文件，然后回头重新输入一遍。这个命令解决了这个问题，可以输入命令数到一半，及时查看需要的文件

#### 关于检查Fastq的Phred值

```shell
cat file.fq | awk 'NR%4==0' | tr -d '\n' | hexdump -v -e'/1 "%u\n"' | sort -nu
```

输出结果 **33-93** (Sanger/Illumina1.8), **64-104**(Illumina1.3 or Illumina1.5) and **59-104** (Solexa)

