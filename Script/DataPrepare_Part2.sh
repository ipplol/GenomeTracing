#This Sciprt is for Data preparing and mutation rate calculation
#Please make sure that DataPrepare_Part1.sh is done.
#!/bin/bash

echo "/// Please make sure that DataPrepare_Part1.sh and all phylogenetic scripts were proper finished. ///"
echo `date "+%Y-%m-%d %H:%M:%S"`
thisfold="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mainDir=`dirname $thisfold`
echo $thisfold

#-------------data preparing--------------
cat $mainDir/Data/*.mutevent > $mainDir/Data/Tree_LineageCombined_Total.mutevent

#merge hotspot and backmutation from each tree
$thisfold/MergeHotspotBackmutation.out $mainDir/Data

#calculate linkage disequilibrium
$thisfold/LinkageDisequilibriumRMatrix_MapReduce.out $mainDir/Data

#calculate mutation rate
$thisfold/MutationRate_Query.out $mainDir/Data 3037C/T $mainDir/Data/tmp.o

#sample query test
$thisfold/GenotypeDistanceCalculator.out $mainDir/Data EPI_ISL_7734032:8393G/A,29301A/G,28881G/A,2832A/G,28311C/T,28271A/T,27807C/T,27259A/C,28883G/C,26270C/T,25708C/T,25584C/T,5386T/G,24424A/T,24503C/T,241C/T,24130C/A,23604C/A,21762C/T,18163A/G,13195T/C,11537A/G,3037C/T,23075T/C,22578G/A,23403A/G,23525C/T,22673T/C,14408C/T,28882G/A,23599T/G,10449C/A,10029C/T,21C/T,15240C/T,26767T/C,23055A/G,25000C/T,10135T/C,23948G/T,22674C/T,22679T/C,24469T/A,22686C/T,22917T/G,23013A/C,21846C/T,22992G/A,22995C/A,5025C/T,23040A/G,23048G/A,23063A/T,23202C/A $mainDir/Data/testQueryResult.txt -LDCL

echo "Done"
echo `date "+%Y-%m-%d %H:%M:%S"`