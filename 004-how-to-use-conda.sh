#################### WES Miniconda Prepare ######################
#################### Made by Bioinfoplanet #####################

################First--> installation #####################
wget https://mirrors.ustc.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh
# then bash Miniconda3-latest-Linux-x86_64.sh to install
HOME=/YOUR/HOME_DIR
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3 #-b自动安装 -p安装位置
# after installation, mv conda into your PATH
echo export PATH=$HOME/miniconda3/bin:$PATH>>~/.bashrc
#then reboot terminal OR source ~/.bashrc to activate conda
source ~/.bashrc

################Second--> config ##########################

# Add conda configs
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/free
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes

# 修改为官方源的方法如下：

conda config --remove-key channels

conda config --add channels r 
conda config --add channels conda-forge 
conda config --add channels bioconda


################Third--> DIY your conda ###################
# create conda environment
conda  create -n wes  python=2 bwa
#查看当前conda环境（带*是默认的）
conda info --envs
source activate wes
source deactivate
#conda info --envs 再次查看一下（*是不是移到了wes环境的前面，另外命令提示符首位会出现（wes））
#在wes环境下安装软件
conda install -y sra-tools samtools vcftools  snpeff multiqc qualimap
#【要安装指定版本的软件】=》先搜索版本（比如samtools）
conda search samtools
#然后在软件名后加版本号
conda install samtools=1.8 -y
#列出所有软件列表
conda list
#删除一个环境中某个软件
conda remove -n wes samtools -y
#删除整个环境
conda remove  -n wes --all
#或者直接rm -rf $HOME/miniconda3/envs/wes/

# conda更新
conda update -n base -c defaults conda
################Fourth--> export ###################
#导出当前你的环境，可以分享给别人
onda env export -n wes -f wes.yml
#导入环境
conda env create -f=wes.yml 
#检查是否导入成功
source activate wes
conda list
















