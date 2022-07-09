#!/bin/bash

echo “Valkyrie Weekly Update Start:”
echo `date "+%Y-%m-%d %H:%M:%S"`

#编译脚本
#g++ main.cpp -o Vcf2mut5col.out -std=c++11 -lpthread
#g++ main.cpp -o Mutresult2Sample_mutlistLINUX.out -std=c++11 -lpthread
#g++ main.cpp -o Seq2GenotypeLineageCombine.out -std=c++11 -lpthread
#g++ main.cpp -o GenotypeDistanceCalculator.out -std=c++11 -lpthread

#创建文件夹
#mkdir -p /home/mawentai/Valkyrie/LinuxNew
mkdir -p /home/mawentai/Valkyrie/LinuxNew/GenotypeDetailFile

#下载vcf和metadata
cd /home/mawentai/Valkyrie/LinuxNew
wget -O metadata.tsv.zip --content-disposition https://ngdc.cncb.ac.cn/ncov/genome/export/meta
wget -O 2019-nCoV_total.vcf.gz --content-disposition https://download.cncb.ac.cn/GVM/Coronavirus/vcf/2019-nCoV_total.vcf.gz  
unzip metadata.tsv.zip
gunzip 2019-nCoV_total.vcf.gz

#bigd下载的vcf转mut5col
/home/mawentai/APP/CPlusPlusShell/Vcf2mut5col/Vcf2mut5col.out /home/mawentai/Valkyrie/LinuxNew/2019-nCoV_total.vcf /home/mawentai/Valkyrie/LinuxNew

#mut5col转sample_mutlist
/home/mawentai/APP/CPlusPlusShell/Mutresult2Sample_mutlistLINUX/Mutresult2Sample_mutlistLINUX.out /home/mawentai/Valkyrie/LinuxNew
rm -rf /home/mawentai/Valkyrie/LinuxNew/2019-nCoV_total.vcf &
rm -rf /home/mawentai/Valkyrie/LinuxOld &
cat /home/mawentai/Valkyrie/*.mutlist >> /home/mawentai/Valkyrie/LinuxNew/sample_mutlist.tsv



#去重复，划分基因型
/home/mawentai/APP/CPlusPlusShell/Seq2GenotypeLineageCombine/Seq2GenotypeLineageCombine.out /home/mawentai/Valkyrie/LinuxNew/sample_mutlist.tsv /home/mawentai/Valkyrie/LinuxNew/metadata.tsv /home/mawentai/Valkyrie/LinuxNew

#拷贝之前的突变特征数据
for f in `find /home/mawentai/Valkyrie/Linux -name "Tree_LineageCombined_*"`
do 
\cp -f $f /home/mawentai/Valkyrie/LinuxNew
done

\cp -f /home/mawentai/Valkyrie/Linux/TotalBackmutRate.txt /home/mawentai/Valkyrie/LinuxNew
\cp -f /home/mawentai/Valkyrie/Linux/TotalMutationRate.txt /home/mawentai/Valkyrie/LinuxNew
\cp -f /home/mawentai/Valkyrie/Linux/TreeFileLineageDetail.tsv /home/mawentai/Valkyrie/LinuxNew

#调换文件夹名称
mv /home/mawentai/Valkyrie/Linux /home/mawentai/Valkyrie/LinuxOld
mv /home/mawentai/Valkyrie/LinuxNew /home/mawentai/Valkyrie/Linux

#清空旧数据

