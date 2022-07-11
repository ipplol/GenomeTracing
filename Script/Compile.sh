#该脚本用于编译c++代码，准备数据
#This Sciprt is for C++ compiling and Data preparing
#!/bin/bash

echo "Data preparing [Some warnings may show up during compiling which are normal, please relax.]:"
echo `date "+%Y-%m-%d %H:%M:%S"`
thisfold="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $thisfold

#----------------Compiling--------------
cd $thisfold
g++ ./Vcf2mut5col/Vcf2mut5col/main.cpp -o Vcf2mut5col.out -std=c++11 -lpthread
g++ ./Seq2GenotypeLineageCombine/Seq2GenotypeLineageCombine/main.cpp -o Seq2GenotypeLineageCombine.out -std=c++11 -lpthread
g++ ./Mutresult2Sample_mutlistLINUX/Mutresult2Sample_mutlistLINUX/main.cpp -o Mutresult2Sample_mutlistLINUX.out -std=c++11 -lpthread
g++ ./MutationRate_Query/MutationRate_Query/main.cpp -o MutationRate_Query.out -std=c++11 -lpthread
g++ ./MergeHotspotBackmutation/MergeHotspotBackmutation/main.cpp -o MergeHotspotBackmutation.out -std=c++11 -lpthread
g++ ./LinkageDisequilibriumRMatrix_MapReduce/LinkageDisequilibriumRMatrix_MapReduce/main.cpp -o LinkageDisequilibriumRMatrix_MapReduce.out -std=c++11 -lpthread
g++ ./GenotypeDistanceCalculator/GenotypeDistanceCalculator/main.cpp -o GenotypeDistanceCalculator.out -std=c++11 -lpthread
g++ ./Mutation2Sequence/Mutation2Sequence/main.cpp -o Mutation2Sequence.out -std=c++11 -lpthread
g++ ./FindHotspotBackmutation_IQTreeVersion/FindHotspotBackmutation_IQTreeVersion/main.cpp -o FindHotspotBackmutation_IQTreeVersion.out -std=c++11 -lpthread

echo "Done"
echo `date "+%Y-%m-%d %H:%M:%S"`