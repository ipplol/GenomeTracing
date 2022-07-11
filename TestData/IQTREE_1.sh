#PBS -q core24
#PBS -l mem=60gb,nodes=1:ppn=12,walltime=900:00:00
#PBS -o Valk1.o
#PBS -e Valk1.e
#HSCHED -s Valk+Phylo+Recon
#PPN limit 12
#!/bin/bash
thisfold="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mainDir=`dirname $thisfold`
cd $thisfold
$mainDir/Script/Mutation2Sequence.out $mainDir/Data/reference.fa $thisfold/inDel.info $thisfold/Tree_LineageCombined_1.tsv $thisfold/Tree_LineageCombined_1.fa
$mainDir/Script/iqtree2 -s Tree_LineageCombined_1.fa -m GTR+F -asr -nt AUTO -o Wuhan-Hu-1
$mainDir/Script/FindHotspotBackmutation_IQTreeVersion.out $thisfold 1
#rm Tree_LineageCombined_1.fa.state
rm Tree_LineageCombined_1.fa.mldist
#rm Tree_LineageCombined_1.fa
