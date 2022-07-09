#bigd下载的vcf转mut5col
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/Vcf2mut5col/Vcf2mut5col/bin/x64/Debug/Vcf2mut5col.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux/2019-nCoV_total.vcf /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux

#mut5col转sample_mutlist
/mnt/c/Users/lilab/OneDrive/BIG/程序/变异位点/Mutresult2Sample_mutlistLINUX/Mutresult2Sample_mutlistLINUX/bin/x64/Debug/Mut5col2Sample_mutlistLINUX.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux

#去重复，划分基因型
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/Seq2GenotypeLineageCombine/Seq2GenotypeLineageCombine/bin/x64/Debug/Seq2GenotypeLineageCombine.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux/sample_mutlist.tsv /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux/metadata.tsv /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux

#dsub
for f in `find /xtdisk/limk_group/mawt/Valkyrie/LinuxNewUpdate -name "*.sh"`
do 
echo $f
dsub $f
done

#dsub
for((i=100;i<=108;i++));
do
echo $i
dsub IQTREE_$i.sh
done

#压缩genotype文件
tar -zcvf GenotypeDetailFile.tar.gz GenotypeDetailFile

#合并mutevent
cat *.mutevent > Tree_LineageCombined_Total.mutevent

#merge不同的树找出来的突变热点，回复突变
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/MergeHotspotBackmutation/MergeHotspotBackmutation/bin/x64/Debug/MergeHotspotBackmutation.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux

#计算连锁不平衡
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/LinkageDisequilibriumRMatrix_MapReduce/LinkageDisequilibriumRMatrix_MapReduce/bin/x64/Debug/LinkageDisequilibriumRMatrix_MapReduce.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux

#计算突变率
/mnt/c/Users/lilab/OneDrive/BIG/程序/Valkyrie/MutationRate_Query/MutationRate_Query/bin/x64/Debug/MutationRate_Query.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux 3037C/T /mnt/g/TMP/out.txt

#解压缩genotype文件
nohup tar -zxvf GenotypeDetailFile.tar.gz &

#ncweb目录
/xtdisk/ncov_group/ncweb/webShils/

#需要的文件
GenotypeDetailTotal.tsv
GenotypeLocationCM.tsv
inDel.info
LinkageDisequilibriumPair_Circox.txt
LinkageDisequilibriumR.matrix

#删除旧的Linux依赖项
nohup rm -rf GenotypeDetailFile &

#测试数据
DeltaI:210G/T,241C/T,410G/T,509GGTCATGTTA/G,3037C/T,5184C/T,5584A/G,9429T/C,9891C/T,11418T/C,11514C/T,11665C/T,13019C/T,14408C/T,15451G/A,16466C/T,21618C/G,21987G/A,22028GAGTTCA/G,22227C/T,22917T/G,22995C/A,23403A/G,23604C/G,24410G/A,25469C/T,26767T/C,27571AGCACTCAATTT/A,27638T/C,27677A/C,27752C/T,28247AGATTTC/A,28270TA/T,28461A/G,28881G/T,29402G/T,29742G/T

#本地测试
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/GenotypeDistanceCalculator/GenotypeDistanceCalculator/bin/x64/Debug/GenotypeDistanceCalculator.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux EPI_ISL_7734032:8393G/A,29301A/G,28881G/A,2832A/G,28311C/T,28271A/T,27807C/T,27259A/C,28883G/C,26270C/T,25708C/T,25584C/T,5386T/G,24424A/T,24503C/T,241C/T,24130C/A,23604C/A,21762C/T,18163A/G,13195T/C,11537A/G,3037C/T,23075T/C,22578G/A,23403A/G,23525C/T,22673T/C,14408C/T,28882G/A,23599T/G,10449C/A,10029C/T,21C/T,15240C/T,26767T/C,23055A/G,25000C/T,10135T/C,23948G/T,22674C/T,22679T/C,24469T/A,22686C/T,22917T/G,23013A/C,21846C/T,22992G/A,22995C/A,5025C/T,23040A/G,23048G/A,23063A/T,23202C/A /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux/result.txt -LDCL

#比较两种画树方法
cat GenotypeDetailTotal.tsv.sampled | while read line
do
#echo $line;
gtname=`echo $line|cut -d ':' -f1`
echo $gtname
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/GenotypeDistanceCalculator/GenotypeDistanceCalculator/bin/x64/Debug/GenotypeDistanceCalculator.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20210901/SepTreeTest/separate $line /mnt/g/VariationMutation/回复突变/iqtreePopulation20210901/SepTreeTest/result/$gtname.sep_result.txt
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/GenotypeDistanceCalculator/GenotypeDistanceCalculator/bin/x64/Debug/GenotypeDistanceCalculator.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20210901/SepTreeTest/single $line /mnt/g/VariationMutation/回复突变/iqtreePopulation20210901/SepTreeTest/result/$gtname.sin_result.txt
done;

cat GenotypeDetailTotal.tsv.sampled |cut -d':' -f1|while read id
do
gtname=${line#:*}
echo $gtname
done;

#拆分合并在一起的fasta
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/FastaUnMerge/FastaUnMerge/bin/x64/Debug/FastaUnMerge.out /mnt/g/tmp/usa/USAJun.tsv

#从fasta序列到突变list
for f in `find /mnt/g/tmp/usa -name "*.fasta"`
do
echo $f
cat /mnt/g/WH/MN908947.fasta $f > $f.fas
mafft --auto $f.fas > $f.mafft
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/MSA2Sample_Mutlist/MSA2Sample_Mutlist/bin/x64/Debug/MSA2Sample_Mutlist.out $f.mafft $f.mutlist
done

#批量找最近序列
for f in `cat /mnt/g/ChinaCoV/CDC2201-2203/fasta/Sample_mutlist.tsv`
do
#echo $f
name=${f%%:*}
name1=`basename $name`
echo $name
echo $name1
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/GenotypeDistanceCalculator/GenotypeDistanceCalculator/bin/x64/Debug/GenotypeDistanceCalculator.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux $f /mnt/g/ChinaCoV/CDC2201-2203/fasta/result/$name1.WM.result
#/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/GenotypeDistanceCalculator/GenotypeDistanceCalculator/bin/x64/Debug/GenotypeDistanceCalculator.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20211101/Linux $f /mnt/g/VariationMutation/Omicron/1130_182/182/$name1.WM.result
#/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/GenotypeDistanceCalculator/GenotypeDistanceCalculator/bin/x64/Debug/GenotypeDistanceCalculator.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20211101/Linux $f /mnt/g/VariationMutation/Omicron/1130_182/182/$name1.LD.result -LD
#/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/GenotypeDistanceCalculator/GenotypeDistanceCalculator/bin/x64/Debug/GenotypeDistanceCalculator.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20210901/Linux $f /mnt/g/VariationMutation/回复突变/iqtreePopulation20210901/Linux/ResultTest/Result/$name1.CL.result -CL
#/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/GenotypeDistanceCalculator/GenotypeDistanceCalculator/bin/x64/Debug/GenotypeDistanceCalculator.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux $f /mnt/g/tmp/$name1.SM.result -SM
#/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/GenotypeDistanceCalculator/GenotypeDistanceCalculator/bin/x64/Debug/GenotypeDistanceCalculator.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux $f /mnt/g/NCBISRA/疑似新冠0215/244/$name1.SM.result -SM
done

#批量找最近序列 考虑覆盖区域
for f in `cat /mnt/g/NCBISRA/疑似新冠0326/mut_result.txt`
do
#echo $f
name=${f%%:*}
name1=`basename $name`
name2=${name1%%.*}
echo $name2
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/GenotypeDistanceCalculator_Coverage/GenotypeDistanceCalculator_Coverage/bin/x64/Node/GenotypeDistanceCalculator_Coverage.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux $f /mnt/g/NCBISRA/疑似新冠0326/result/$name1.MOD.result -MODepth /mnt/g/NCBISRA/疑似新冠0326/readcounts/$name2.readcounts
#/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/GenotypeDistanceCalculator_Coverage/GenotypeDistanceCalculator_Coverage/bin/x64/Node/GenotypeDistanceCalculator_Coverage.out /mnt/g/VariationMutation/回复突变/iqtreePopulation20220630/Linux $f /mnt/g/NCBISRA/疑似新冠0215/244/$name1.SMD.result -SMDepth /mnt/g/NCBISRA/疑似新冠0215/244/readcounts_5x/$name1.readcounts
done

#文件前N行
for f in `find /mnt/g/ChinaCoV/CDC2201-2203/fasta -name "*.result"`
do
sed -n '3,3p' $f > $f.top
done

#iqtree
iqtree -s TotalLength.mafft.cut.fas -m GTR+F -asr -nt AUTO -o Wuhan-Hu-1

#祖先序列突变
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/FindHotspotBackmutation_IQTreeVersion/FindHotspotBackmutation_IQTreeVersion/bin/x64/Debug/FindHotspotBackmutation_IQTreeVersion.out /mnt/g/TMP/Omicron/1130_182 182

#beast
java -Xmx60g -jar /mnt/g/APP/beast/lib/launcher.jar -threads 24 -beagle_CPU -beagle_single -working /mnt/g/VariationMutation/Omicron/fastq/D5_90/DelDel/Tree/TMRCA/CSSSE_HKY_fixclock.xml

for f in `find /mnt/g/VariationMutation/回复突变/iqtreePopulation20210901/Linux/ResultTest/Result -name "*.treefile"`
do 
echo $f
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/PhyloTreeNodeDistanceCalculator/PhyloTreeNodeDistanceCalculator/bin/x64/Debug/PhyloTreeNodeDistanceCalculator.out $f
done

for f in `find /mnt/g/VariationMutation/DVG -name "*.readcounts"`
do 
echo $f
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/CalliSNVfromReadcounts/CalliSNVfromReadcounts/bin/x64/Debug/CalliSNVfromReadcounts.out -i $f
done

for f in `find /mnt/g/tmp -name "*.readcounts"`
do 
echo $f
/mnt/c/Users/lilab/OneDrive/BIG/程序/Linux下的C++脚本/CallConsensusSequence/CallConsensusSequence/bin/x64/Debug/CallConsensusSequence.out -i $f
done

echo 'export PATH='$(pwd)'/octopus/bin:$PATH' >> ~/.bash_profile