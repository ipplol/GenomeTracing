# A fast and accurate method for SARS-CoV-2 genomic tracing

## Installation

Download the zip file and unzip it. Enter the main folder and,

- **`bash ./Script/Compile.sh`**

This will automatically compile several scripts that written in C++.

Then, to check if everything is properly settled, a test file is pre-prepared using data from Nextstrain database.

- **`bash ./TestData/DataPrepare_Test.sh`**

It uses two files, `sample_mutlist.tsv` and `metadata.tsv` ,in the `TestData` fold as input.

And this script will calculate the mutation rate and other metrics for genomic tracing.

The test query sequence is EPI_ISL_7734032. And the result will be found at `./TestData/testQueryResult.txt`

## Database prepare

**One can skip this step and use our pre-calculated database instead.**

**Simply download it from [Google Drive](https://drive.google.com/drive/folders/1feiqGvoKvP9NDxh__GTdzaMiNpTrpoMR?usp=sharing), unzip and put them into the `Data` fold.**

Two input files [metadata.tsv](https://ngdc.cncb.ac.cn/ncov/genome/export/meta) and [2019-nCoV_total.vcf](https://download.cncb.ac.cn/GVM/Coronavirus/vcf/2019-nCoV_total.vcf.gz) ,are needed to build the database.

The Data_Download script can be used to download them using wget.

- **`bash ./Script/Data_Download.sh`**

If the download speed is slow, you can also download them through ftp manually and put them into the `Data` fold.

Then, the Data_Prepare scripts perform data formatting and calculate the mutation metrics for the genome tracing

- **`bash ./Script/DataPrepare_Part1.sh`**

The Part_1 script generate multiple `IQTREE_*.sh` scripts for phylogenetic reconstruction, and implement one by one.

After all trees are successfully constructed, fun the part_2 script to obtain the mutation metrics

- **`bash ./Script/DataPrepare_Part2.sh`**


    
##  Find the closest hit to your query

- **`./Script/GenotypeDistanceCalculator.out <path to the Data folder> <mutations of your query> <path to the result file> <scoring method>`**

for example:

- **`./Script/GenotypeDistanceCalculator.out ./Data EPI_ISL_7734032:8393G/A,29301A/G,28881G/A,2832A/G,28311C/T,28271A/T,11537A/G,3037C/T ./Data/testQueryResult.txt -LDCL`**

there are 5 scoring methods available, the default one is WM:

**`-WM`** -> Weighted mutation

**`-MO`** -> Unweighted method

**`-LD`** -> Weighted mutation with linkage disequilibrium considerred

**`-CL`** -> Lineage-specific weighted mutation

**`-LDCL`** -> Weighted mutation with both linkage disequilibrium and lineage considerred

### *ValidationResut.zip*
It includes files and scripts used in the Validation Section in the manuscript.

`/CollectionLocationTest` 2000 randomly selected sequences

`/Simulation` Simulatied viral genome with known ancestry information

`/TransmissionPair` 563 sequences with highly plausible transmission links inferred from the genomic
