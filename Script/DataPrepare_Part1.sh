#This Sciprt is for Data preparing
#!/bin/bash

echo "Data preparing [Some warnings may show up during compiling which are normal, please relax.]:"
echo `date "+%Y-%m-%d %H:%M:%S"`
thisfold="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $thisfold

#--------------data update-----------------
mainDir=`dirname $thisfold`
cd "$mainDir/Data"
unzip metadata.tsv.zip
gunzip 2019-nCoV_total.vcf.gz

#-------------data preparing--------------
echo "Changing data format. It may take a while. 5~8 hours"

#change vcf to mut5col
$thisfold/Vcf2mut5col.out $mainDir/Data/2019-nCoV_total.vcf $mainDir/Data

#change mut5col to sample_mutlist
$thisfold/Mutresult2Sample_mutlistLINUX.out $mainDir/Data
#rm -rf $mainDir/Data/2019-nCoV_total.vcf &

#Sequence dedup and merge into genotypes
$thisfold/Seq2GenotypeLineageCombine.out $mainDir/Data/sample_mutlist.tsv $mainDir/Data/metadata.tsv $mainDir/Data

#-------------phylogenetic reconstruction------------
#/// This step was designed to run on a compute cluster ///
#/// This step may take 2~3 days ///
#dsub the script of each tree
#if dsub is not available, use bash
for f in `find $mainDir/Data -name "*.sh"`
do 
echo $f
#dsub $f
bash $f
done

echo "Done"
echo `date "+%Y-%m-%d %H:%M:%S"`