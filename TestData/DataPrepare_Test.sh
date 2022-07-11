#This Sciprt is for C++ compiling and Data preparing
#!/bin/bash

echo "Data preparing [Some warnings may show up during compiling which are normal, please relax.]:"
echo `date "+%Y-%m-%d %H:%M:%S"`
thisfold="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mainDir=`dirname $thisfold`
thisfold="$mainDir/Script"
echo $thisfold

#----------------Compiling--------------
#cd $thisfold
#g++ ./Vcf2mut5col/Vcf2mut5col/main.cpp -o Vcf2mut5col.out -std=c++11 -lpthread
#g++ ./Seq2GenotypeLineageCombine/Seq2GenotypeLineageCombine/main.cpp -o Seq2GenotypeLineageCombine.out -std=c++11 -lpthread
#g++ ./Mutresult2Sample_mutlistLINUX/Mutresult2Sample_mutlistLINUX/main.cpp -o Mutresult2Sample_mutlistLINUX.out -std=c++11 -lpthread
#g++ ./MutationRate_Query/MutationRate_Query/main.cpp -o MutationRate_Query.out -std=c++11 -lpthread
#g++ ./MergeHotspotBackmutation/MergeHotspotBackmutation/main.cpp -o MergeHotspotBackmutation.out -std=c++11 -lpthread
#g++ ./LinkageDisequilibriumRMatrix_MapReduce/LinkageDisequilibriumRMatrix_MapReduce/main.cpp -o LinkageDisequilibriumRMatrix_MapReduce.out -std=c++11 -lpthread
#g++ ./GenotypeDistanceCalculator/GenotypeDistanceCalculator/main.cpp -o GenotypeDistanceCalculator.out -std=c++11 -lpthread
#g++ ./Mutation2Sequence/Mutation2Sequence/main.cpp -o Mutation2Sequence.out -std=c++11 -lpthread
#g++ ./FindHotspotBackmutation_IQTreeVersion/FindHotspotBackmutation_IQTreeVersion/main.cpp -o FindHotspotBackmutation_IQTreeVersion.out -std=c++11 -lpthread

#Some steps are skipped since we are running the test data pre-prepared
#--------------data update-----------------
#cd "$mainDir/TestData"
#echo "wget can be slow, if so, please use ftp download manually [filezilla]"
#wget -O metadata.tsv.zip --content-disposition https://ngdc.cncb.ac.cn/ncov/genome/export/meta
#wget -O 2019-nCoV_total.vcf.gz --content-disposition https://download.cncb.ac.cn/GVM/Coronavirus/vcf/2019-nCoV_total.vcf.gz  
#unzip metadata.tsv.zip
#gunzip 2019-nCoV_total.vcf.gz


#-------------data preparing Part 1--------------
#echo "Changing data format. It may take a while. 5~8 hours"

#change vcf to mut5col
#$thisfold/Vcf2mut5col.out $mainDir/TestData/2019-nCoV_total.vcf $mainDir/TestData

#change mut5col to sample_mutlist
#$thisfold/Mutresult2Sample_mutlistLINUX.out $mainDir/TestData
#rm -rf $mainDir/TestData/2019-nCoV_total.vcf &

#Sequence dedup and merge into genotypes
$thisfold/Seq2GenotypeLineageCombine.out $mainDir/TestData/sample_mutlist.tsv $mainDir/TestData/metadata.tsv $mainDir/TestData

#-------------phylogenetic reconstruction------------
#/// This step was designed to run on a compute cluster ///
#dsub the script of each tree
#if dsub is not available, use bash
for f in `find $mainDir/TestData -name "IQTREE*"`
do 
echo $f
#dsub $f
#bash $f
done


#-------------data preparing Part 2--------------
cat $mainDir/TestData/*.mutevent > $mainDir/TestData/Tree_LineageCombined_Total.mutevent

#merge hotspot and backmutation from each tree
$mainDir/Script/MergeHotspotBackmutation.out $mainDir/TestData

#calculate linkage disequilibrium
#The test dataset has too few sequences to estimate the LD
#$mainDir/Script/LinkageDisequilibriumRMatrix_MapReduce.out $mainDir/TestData

#calculate mutation rate
$mainDir/Script/MutationRate_Query.out $mainDir/TestData 3037C/T $mainDir/TestData/tmp.o

#sample query test
$mainDir/Script/GenotypeDistanceCalculator.out $mainDir/TestData EPI_ISL_7734032:8393G/A,29301A/G,28881G/A,2832A/G,28311C/T,28271A/T,27807C/T,27259A/C,28883G/C,26270C/T,25708C/T,25584C/T,5386T/G,24424A/T,24503C/T,241C/T,24130C/A,23604C/A,21762C/T,18163A/G,13195T/C,11537A/G,3037C/T,23075T/C,22578G/A,23403A/G,23525C/T,22673T/C,14408C/T,28882G/A,23599T/G,10449C/A,10029C/T,21C/T,15240C/T,26767T/C,23055A/G,25000C/T,10135T/C,23948G/T,22674C/T,22679T/C,24469T/A,22686C/T,22917T/G,23013A/C,21846C/T,22992G/A,22995C/A,5025C/T,23040A/G,23048G/A,23063A/T,23202C/A $mainDir/TestData/testQueryResult.txt -LDCL

echo "Done"
echo `date "+%Y-%m-%d %H:%M:%S"`