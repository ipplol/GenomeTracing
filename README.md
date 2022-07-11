# A fast and accurate method for SARS-CoV-2 genomic tracing

Here store the source codes of the genome tracing platform 

https://ngdc.cncb.ac.cn/ncov/online/tool/genome-tracing/

## Installation

Download the zip file and unzip it. Enter the main fold and,

- **`bash ./Script/Compile.sh`**

This will automatically compile several scripts that have been written in C++.

Then, to check if everything is settled, a simple example is pre-prepared using data from Nextstrain global open.

- **`bash ./TestData/DataPrepare_Test.sh`**

It starts with two files, `sample_mutlist.tsv` and `metadata.tsv` ,in the `TestData` fold.

And it will calculate the mutation rate and other parameters for genomic tracing.

The test query sequence is EPI_ISL_7734032. And the result will be at `./TestData/testQueryResult.txt`

## Database prepare

**One can skip this step and use our pre-calculated database instead.**

**Simply download it from GITHUB, unzip and put them into the `Data` fold.**

The software that finds the closest hit to the query relies on a genotype list and several mutation characteristic files.

To generate those files, two input [metadata.tsv](https://ngdc.cncb.ac.cn/ncov/genome/export/meta) and [2019-nCoV_total.vcf](https://download.cncb.ac.cn/GVM/Coronavirus/vcf/2019-nCoV_total.vcf.gz) ,are needed.

The Data_Download script will try to download them using wget. 

- **`bash ./Script/Data_Download.sh`**

If the download speed is too slow, please download them through ftp manually and put them into the `Data` fold.

Then, the Data_Prepare scripts will perform data format change and mutation characteristic calculation.

- **`bash ./Script/DataPrepare_Part1.sh`**

The Part_1 script generate multiple `IQTREE_*.sh` scripts for phylogenetic reconstruction.

After all of them are finished,

- **`bash ./Script/DataPrepare_Part2.sh`**

	to obtain mutation rates and other files.
    
##  Find the closest hit to your query

- **`./Script/GenotypeDistanceCalculator.out <path to the Data fold> <mutations of your query> <path to the result file> <scoring method>`**

for example:

- **`./Script/GenotypeDistanceCalculator.out ./Data EPI_ISL_7734032:8393G/A,29301A/G,28881G/A,2832A/G,28311C/T,28271A/T,11537A/G,3037C/T ./Data/testQueryResult.txt -LDCL`**

there are 5 scoring methods available, the default one is WM:

**`-WM`** -> Weighted mutation

**`-MO`** -> Unweighted method

**`-LD`** -> Weighted mutation with linkage disequilibrium

**`-CL`** -> Weighted mutation with lineage rate

**`-LDCL`** -> Weighted mutation with both linkage disequilibrium and lineage rate
