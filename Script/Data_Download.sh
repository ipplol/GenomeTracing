#!/bin/bash

thisfold="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $thisfold

#--------------data update-----------------
mainDir=`dirname $thisfold`
cd "$mainDir/Data"
echo "wget can be slow, if so, please use ftp download manually [filezilla]"
wget -O metadata.tsv.zip --content-disposition https://ngdc.cncb.ac.cn/ncov/genome/export/meta
wget -O 2019-nCoV_total.vcf.gz --content-disposition https://download.cncb.ac.cn/GVM/Coronavirus/vcf/2019-nCoV_total.vcf.gz  
unzip metadata.tsv.zip
gunzip 2019-nCoV_total.vcf.gz

echo "Done"
echo `date "+%Y-%m-%d %H:%M:%S"`