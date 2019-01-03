# TCGA的各种信息下载相关

> 刘小泽写于19.1.3

### 关于CNV下载

CNV（Cpoy Number Variant）指的是拷贝数变异，又称拷贝数多态性（Cpoy Number Polymorphism, CNP），大小在1kb~3Mb的DNA片段。覆盖的核苷酸总数远超SNP，临床上将CNV分为：致病性CNV、非致病性CNV和临床意义不明CNV。

肿瘤研究着重看somatic情况，TGCA数据库存放了somatic CNV，这些CNV是经过肿瘤病人的拷贝数变异信息剔除掉对照中的CNV后留下的，可能与癌症相关的CNV。

TGCA中的somatic CNV主要利用Affymetrix SNP6.0 array芯片，得到的数据有4个开放程度，一般都下载level 3的

>  去Broad Institute中下载

```shell
# level3 存放于
https://gdac.broadinstitute.org/runs/stddata__latest/data/

# 想下载BRCA基因，感兴趣的是snp6，并且基于hg19
wget -c -r -np -nH -k -L --cut-dirs 6 -p -A "*snp_6*Level_3*hg19*" https://gdac.broadinstitute.org/runs/stddata__latest/data/BRCA/20160128/
# 下载后的minus的germline_cnv就是癌症相关的somatic CNV
```

