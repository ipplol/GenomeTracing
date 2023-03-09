# A fast and accurate method for SARS-CoV-2 genomic tracing

This website provides access to the codes and data used to deploy our new genomic tracing method, Valkyrie. Note that this method is also available on the webserver of China National Center for Bioinformation at https://ngdc.cncb.ac.cn/ncov/online/tool/genome-tracing/.

Users have the option to perform the analysis either on the online platform or using the standalone version. The virus database is updated monthly on the webserver, while the mutation metrics is updated every three months. However, please be aware that the computational resources used by the webserver are shared among multiple applications, which may cause occasional slow performance or even temporary unresponsiveness.

## Installation

Download the zip file and unzip the file. Enter the main folder and run:

- **`bash ./Script/Compile.sh`**

This will automatically compile several scripts that written in C++.

Then, to examine if everything is properly settled, a test file is pre-prepared using data from Nextstrain database.

- **`bash ./TestData/DataPrepare_Test.sh`**

The script accepts two files, `sample_mutlist.tsv` and `metadata.tsv` in the `TestData` folder as input.

This script calculates the mutation rate and other metrics for genomic tracing.

The test query sequence is EPI_ISL_7734032. And the result will be found at `./TestData/testQueryResult.txt`

## Database prepare

**(This step can be skipped and use our pre-calculated database instead. Simply download it from [Google Drive](https://drive.google.com/drive/folders/1feiqGvoKvP9NDxh__GTdzaMiNpTrpoMR?usp=sharing), unzip and put them into the `Data` folder.)**

Two files [metadata.tsv](https://ngdc.cncb.ac.cn/ncov/genome/export/meta) and [2019-nCoV_total.vcf](https://download.cncb.ac.cn/GVM/Coronavirus/vcf/2019-nCoV_total.vcf.gz) are needed to build the database.

The Data_Download script can be used to download the above two files using wget.

- **`bash ./Script/Data_Download.sh`**

If the download speed is slow, you can also download them manually through ftp and put them into the `Data` folder.

Then, run the following two Data_Prepare scripts to format the data and calculate the mutation metrics used for the genome tracing

- **`bash ./Script/DataPrepare_Part1.sh`**

The Part_1 script generates multiple `IQTREE_*.sh` scripts for phylogenetic reconstruction, please execute one by one.

After successfully constructing all trees, run the part_2 script to obtain the overall mutation metrics

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

Please strickly follow the format of the mutation: seq_ID:POS1Allel1/Allele2,POS2Allel1/Allele2,... or you can use the fasta format sequence as input on the webserver.
