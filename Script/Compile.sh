#该脚本用于编译c++代码，准备数据
#This Sciprt is for C++ compiling and Data preparing
#!/bin/bash

echo "Data preparing [Some warnings may show up during compiling which are normal, please relax.]:"
echo `date "+%Y-%m-%d %H:%M:%S"`
thisfold="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $thisfold

#----------------Compiling--------------
cd $thisfold
g++ ./CPPSourceCode/Vcf2mut5col.cpp -o Vcf2mut5col.out -std=c++11 -lpthread
g++ ./CPPSourceCode/Seq2GenotypeLineageCombine.cpp -o Seq2GenotypeLineageCombine.out -std=c++11 -lpthread
g++ ./CPPSourceCode/Mutresult2Sample_mutlistLINUX.cpp -o Mutresult2Sample_mutlistLINUX.out -std=c++11 -lpthread
g++ ./CPPSourceCode/MutationRate_Query.cpp -o MutationRate_Query.out -std=c++11 -lpthread
g++ ./CPPSourceCode/MergeHotspotBackmutation.cpp -o MergeHotspotBackmutation.out -std=c++11 -lpthread
g++ ./CPPSourceCode/LinkageDisequilibriumRMatrix_MapReduce.cpp -o LinkageDisequilibriumRMatrix_MapReduce.out -std=c++11 -lpthread
g++ ./CPPSourceCode/GenotypeDistanceCalculator.cpp -o GenotypeDistanceCalculator.out -std=c++11 -lpthread
g++ ./CPPSourceCode/Mutation2Sequence.cpp -o Mutation2Sequence.out -std=c++11 -lpthread
g++ ./CPPSourceCode/FindHotspotBackmutation_IQTreeVersion.cpp -o FindHotspotBackmutation_IQTreeVersion.out -std=c++11 -lpthread

echo "Done"
echo `date "+%Y-%m-%d %H:%M:%S"`