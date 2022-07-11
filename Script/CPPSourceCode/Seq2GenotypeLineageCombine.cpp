#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <thread>
#include <pthread.h>
#include <map>
#include <mutex>
#include <unordered_map>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/syscall.h>
using namespace std;

class Metadata//Metadata类
{
public:
	string strainName;
	string accessionID;
	string lineage;
	string host;
	string country;
	string collectionDate;
	string completeness;
	string quailty;
	string originalMetadata;
};
class GenoType//基因型类
{
public:
	//string rSeqName;//代表性序列的名字
	unordered_map<string, int> mutList;//突变列表 1代表SNP 0代表indel
	vector<string> sampleIDList;//所包含序列的ID
	vector<Metadata> sampleMetadataList;//所包含序列的Metadata
	vector<string> sampleLocation;//基因型样本的国家地区【0】地区【1】数量 上下对齐
	vector<string> sampleCollectionDate;//基因型样本的采样时间
	vector<string> sampleLineage;//基因型样本的Lineage信息

	string representiveSeq;//代表性序列(时间最早)的accession ID
	string highCompleteHuman = "false";//该基因型是否具有全长高质量人类宿主的序列
};
class Lineage//lineage类
{
public:
	string lineageName;//lineage的名字
	//vector<int> genotypeList;改用map机制，使默认按采样时间排序
	//map<"日期+样本名称"，序号>
	unordered_map<string,int> genotypeList;//所包含的基因型(合格的，能用于画树的全长高质量人类宿主
};
static long starttime;//程序开始时间
static vector<GenoType> genoTypeList;//储存所有GenoType的列表
static unordered_map<string, Metadata> metadataMap;//strain名称和对应meta信息的hash表
static vector<Lineage> totalLineageList;//储存所有Lineage的列表
static unordered_map<string, int>lineageMap;//lineage名称和对应在lineage表中的位置hash
static Metadata unknownMetadata;//空metadata数据，全部是unknown
static string metatitle;//metadata的表头
static vector<string> indelList;//被重新编码的indel列表，vector保证有顺序
static unordered_map<string, string> indelMap;//indel和重新编码的indel hash
static int qcpassedCount = 0;

//2021.10.21 为保证indel recoding的替换矩阵和SNP一致，需要从SNP中抽样给recoding的
static unordered_map<string, string> SNPrefA;//从原始SNP到碱基替换的map，ref位置为A 685A/T -> <685A/T,A/T>
static unordered_map<string, string> SNPrefT;
static unordered_map<string, string> SNPrefC;
static unordered_map<string, string> SNPrefG;
static vector<string> randomSubstituteA;//SNPrefA的value值集合，用于随机遍历
static vector<string> randomSubstituteT;
static vector<string> randomSubstituteC;
static vector<string> randomSubstituteG;
//--------------------------------------------

std::mutex mtx;//线程锁

//string split函数
vector<string> splitStr(string str, char delimiter)
{
	vector<string> r;
	string tmpstr;
	int i, j, k;
	vector<int> pointList;
	pointList.emplace_back(-1);
	for (i = 0; i < str.length(); i++)
		if (str[i] == delimiter)pointList.emplace_back(i);
	pointList.emplace_back(i);
	for (i = 1; i < pointList.size(); i++)
		r.emplace_back(str.substr(pointList[i - 1] + 1, (pointList[i] - pointList[i - 1] - 1)));

	return r;
}

//indel重新编码的算法
/*string IndelRecoding(string indel)
{
	int i = 0;
	while (indel[i] >= '0' && indel[i] <= '9') i++;
	vector<string> indel1 = splitStr(indel.substr(i, indel.length() - i), '/');
	if (indel1[0].length() > indel1[1].length())
	{
		if (indel1[0][indel1[0].length() - 1] != indel1[1][indel1[1].length() - 1])
			return "" + to_string(30001 + indelMap.size()) + indel1[0][indel1[0].length() - 1] + "/" + indel1[1][indel1[1].length() - 1];//deletion
		else
		{
			if (indel1[0][indel1[0].length() - 1] != 'A')
				return "" + to_string(30001 + indelMap.size()) + "A/" + indel1[1][indel1[1].length() - 1];//deletion
			else
				return "" + to_string(30001 + indelMap.size()) + "C/" + indel1[1][indel1[1].length() - 1];//deletion
		}
	}
	if (indel1[0].length() < indel1[1].length())
	{
		if (indel1[0][indel1[0].length() - 1] != indel1[1][indel1[1].length() - 1])
			return "" + to_string(30001 + indelMap.size()) + indel1[0][indel1[0].length() - 1] + "/" + indel1[1][indel1[1].length() - 1];//insertion
		else
		{
			if (indel1[0][indel1[0].length() - 1] != 'T')
				return "" + to_string(30001 + indelMap.size()) + "T/" + indel1[1][indel1[1].length() - 1];//deletion
			else
				return "" + to_string(30001 + indelMap.size()) + "G/" + indel1[1][indel1[1].length() - 1];//deletion
		}
	}
	cout << "Indel Recoding ERROR: " << indel <<endl;
	cout << indel.substr(i, indel.length() - i) << ";" << endl;
	cout << indel1[0] << ";" << endl;
	cout << indel1[1] << ";" << endl;
	return "Indel Recoding ERROR";
}*/

//makesure the ATCG substitution matrix is the same
string IndelRecoding(string indel, int GTnumber)//保持和SNP一致的核酸替换矩阵，和基因组一样的ATCG比例，随机抽取SNP来进行替换
{
	int id = syscall(SYS_gettid);//线程编号
	srand(time(0)+id+GTnumber);
	//新冠基因组比例 ATCG约为3：3：2：2
	int index = indelMap.size() % 10;
	if (index >= 0 && index <= 2)
		return "" + to_string(30001 + indelMap.size()) + randomSubstituteA[rand() % randomSubstituteA.size()];
	if (index >= 3 && index <= 5)
		return "" + to_string(30001 + indelMap.size()) + randomSubstituteT[rand() % randomSubstituteT.size()];
	if (index >= 6 && index <= 7)
		return "" + to_string(30001 + indelMap.size()) + randomSubstituteC[rand() % randomSubstituteC.size()];
	if (index >= 8 && index <= 9)
		return "" + to_string(30001 + indelMap.size()) + randomSubstituteG[rand() % randomSubstituteG.size()];
}

//突变到位置 from mutation to position
int Mut2Pos(string mut)
{
	int i = 0;
	while (i < mut.length() && mut[i] <= '9' && mut[i] >= '0')
		i++;
	return(stoi(mut.substr(0,i)));
}

//将突变序列排序后输出 哈希值排序
//sort mutation list by hash value
string MutSort(string mutSeq)
{
	int i, j, k;
	vector<string> mutSeqList = splitStr(mutSeq, ' ');
	map<string, string> mutMap;
	for (i = 0; i < mutSeqList.size(); i++)
		if(mutMap.find(mutSeqList[i])==mutMap.end())
			mutMap.insert(pair<string, string>(mutSeqList[i], mutSeqList[i]));
	string output = "";
	for (auto val = mutMap.begin(); val != mutMap.end(); val++)
		output += val->second + " ";
	return output.substr(0, output.size() - 1);
}

//序列去重复，归纳为GenoType
//sequences dedup to genotypes
int SeqDedup2Genotype(string smlfile)
{
	ifstream read(smlfile);
	if (!read.is_open())
	{
		printf("ERROR /// Datainput: Can Not Open SMLfile ///\n"); return -1;
	}
	string line;
	int i, j, k;
	unordered_map<string, int> mutationMap;//用hash单线程去重复，存储突变string和genoTypeList的位置int
	while (getline(read, line))
	{
		if (line[line.length() - 1] == '\r')
			line = line.substr(0, line.length() - 1);
		transform(line.begin(), line.end(), line.begin(), ::toupper);//全部转大写
		vector<string> line1 = splitStr(line, '\t');
		if (line1.size() > 1)
		{
			line1[1] = MutSort(line1[1]);//排一下序

			if (mutationMap.find(line1[1]) != mutationMap.end())//该突变型的基因型已经存在
			{
				j = mutationMap[line1[1]];
				genoTypeList[j].sampleIDList.emplace_back(line1[0]);
			}
			else
			{
				GenoType a;
				vector<string> line2 = splitStr(line1[1], ' ');
				for (i = 0; i < line2.size(); i++)
					if (line2[i] != "") 
					{
						if (!(line2[i][line2[i].length() - 4] >= '0' && line2[i][line2[i].length() - 4] <= '9'))
							a.mutList.insert(pair<string, int>(line2[i], 0));
						else
							a.mutList.insert(pair<string, int>(line2[i], 1));
					}
				a.sampleIDList.emplace_back(line1[0]);
				genoTypeList.emplace_back(a);
				mutationMap.insert(pair<string, int>(line1[1], genoTypeList.size() - 1));
			}
		}
	}
	read.close();
	return 0;
}


//读入metadata
//read in metadata
int ReadMetadata(string metaDataFile)
{
	ifstream read(metaDataFile);
	if (!read.is_open())
	{
		printf("ERROR /// Datainput: Can Not Open Metadata file ///\n"); return -1;
	}
	unknownMetadata.collectionDate = "MetadataUnknown";
	unknownMetadata.completeness = "MetadataUnknown";
	unknownMetadata.country = "MetadataUnknown";
	unknownMetadata.host = "MetadataUnknown";
	unknownMetadata.lineage = "MetadataUnknown";
	unknownMetadata.quailty = "MetadataUnknown";
	unknownMetadata.strainName = "MetadataUnknown";
	unknownMetadata.originalMetadata = "";
	unknownMetadata.accessionID = "-";
	int i, j, k;
	string line;
	getline(read, metatitle);
	if (metatitle[metatitle.length() - 1] == '\r')
		metatitle = metatitle.substr(0, metatitle.length() - 1);
	while (getline(read, line))
	{
		if (line[line.length() - 1] == '\r')
			line = line.substr(0, line.length() - 1);
		transform(line.begin(), line.end(), line.begin(), ::toupper);//全部转大写
		vector<string> line1 = splitStr(line, '\t');
		Metadata a;
		a.collectionDate = line1[10];
		a.completeness = line1[5];
		vector<string> linez = splitStr(line1[11], '/');
		if (linez[0][linez[0].length() - 1] == ' ')
			linez[0] = linez[0].substr(0, linez[0].length() - 1);
		a.country = linez[0];
		a.host = line1[9];
		a.lineage = line1[4];
		a.quailty = line1[7];
		a.accessionID = line1[1];
		a.strainName = line1[0];
		a.originalMetadata = line;
		metadataMap.insert(pair<string, Metadata>(line1[1], a));//加入hash表 用accession number当key
		metadataMap.insert(pair<string, Metadata>(line1[3], a));//加入hash表 Related ID当备用key
		metadataMap.insert(pair<string, Metadata>(line1[0], a));//加入hash表 VirusStrainName当备用key
	}
}

//多线程遍历genotype找metadata,从start搜到end
//find metadata for each genotype
void ThreadSearchMetadata(int start, int end)
{
	int i, j, k;
	for (i = start; i <= end; i++)
	{
		string QC = "TBD";//是否包含高质量，全长，人类宿主的序列 用于画树
		for (j = 0; j < genoTypeList[i].sampleIDList.size(); j++)
		{
			if (metadataMap.find(genoTypeList[i].sampleIDList[j]) != metadataMap.end())//找到id对应的metadata了
			{
				genoTypeList[i].sampleMetadataList.emplace_back(metadataMap[genoTypeList[i].sampleIDList[j]]);
				if (metadataMap[genoTypeList[i].sampleIDList[j]].completeness == "COMPLETE")
					if (metadataMap[genoTypeList[i].sampleIDList[j]].quailty == "HIGH")
						if (metadataMap[genoTypeList[i].sampleIDList[j]].host == "HOMO SAPIENS")
						{
							QC = "PASS";
							genoTypeList[i].highCompleteHuman = "true";
						}
			}
			else
			{
				genoTypeList[i].sampleMetadataList.emplace_back(unknownMetadata);
			}
		}
		int tmpi; string tmps;
		int unknownCount;
		//------------------------------
		vector<string> locationList; vector<int> locationCount;//整理location，排序
		unknownCount = 0;
		for (j = 0; j < genoTypeList[i].sampleIDList.size(); j++)
		{
			auto val = std::find(locationList.begin(), locationList.end(), genoTypeList[i].sampleMetadataList[j].country);
			if (val != locationList.end())
				locationCount[val - locationList.begin()]++;
			else
			{
				if (genoTypeList[i].sampleMetadataList[j].country != "MetadataUnknown")
				{
					locationList.emplace_back(genoTypeList[i].sampleMetadataList[j].country);
					locationCount.emplace_back(1);
				}
				else
					unknownCount++;
			}
		}
		for(j=0;j<locationCount.size();j++)
			for(k=j+1;k<locationCount.size();k++)//按照含有序列的数目顺序排列
				if (locationCount[j] < locationCount[k])
				{
					tmpi = locationCount[j]; locationCount[j] = locationCount[k]; locationCount[k] = tmpi;
					tmps = locationList[j]; locationList[j] = locationList[k]; locationList[k] = tmps;
				}
		if (unknownCount > 0)
		{
			locationCount.emplace_back(unknownCount);
			locationList.emplace_back("MetadataUnknown");
		}
		string locaA = "" + locationList[0]; string locaB = ""+ to_string(locationCount[0]);
		for (j = 1; j < locationCount.size(); j++)
		{
			locaA += "," + locationList[j]; locaB += "," + to_string(locationCount[j]);
		}
		genoTypeList[i].sampleLocation.emplace_back(locaA);//基因型存在的国家 国家1，国家2，.......
		genoTypeList[i].sampleLocation.emplace_back(locaB);//对应的序列数 7，6，2，1........
		//-------------------------------

		//------------------------------
		vector<string> dateList; vector<int> dateCount;//整理Collection Date，排序
		unknownCount = 0;
		int earlestDate = 99999999;//顺带找一下采样时间最早的序列
		string accessionID = "-";
		for (j = 0; j < genoTypeList[i].sampleIDList.size(); j++)
		{
			vector<string> date1 = splitStr(genoTypeList[i].sampleMetadataList[j].collectionDate, '-');//拆分日期yyyy-mm-dd
			string date = "MetadataUnknown";
			int rSdate = 99999999;
			if (date1[0] != "MetadataUnknown")
				if (date1.size() > 1)
				{
					date = to_string(stoi(date1[0])*100 + stoi(date1[1]));//只保留yyyy 和 mm
					if (date1.size() > 2)
						rSdate = stoi(date1[0]) * 10000 + stoi(date1[1]) * 100 + stoi(date1[2]);
					else
						rSdate = stoi(date1[0]) * 10000 + stoi(date1[1]) * 100;
				}
				else
				{
					date = date1[0] + "00";
					rSdate = stoi(date1[0]) * 10000;
				}
			if (rSdate < earlestDate)
			{
				earlestDate = rSdate;
				accessionID = genoTypeList[i].sampleMetadataList[j].accessionID;
			}
			genoTypeList[i].sampleMetadataList[j].collectionDate = date;
			auto val = std::find(dateList.begin(), dateList.end(), date);
			if (val != dateList.end())
				dateCount[val - dateList.begin()]++;
			else
			{
				if (date != "MetadataUnknown")
				{
					dateList.emplace_back(date);
					dateCount.emplace_back(1);
				}
				else
					unknownCount++;
			}
		}
		for (j = 0; j < dateList.size(); j++)
			for (k = j + 1; k < dateList.size(); k++)//按照日期顺序排列
				if (stoi(dateList[j]) > stoi(dateList[k]))
				{
					tmpi = dateCount[j]; dateCount[j] = dateCount[k]; dateCount[k] = tmpi;
					tmps = dateList[j]; dateList[j] = dateList[k]; dateList[k] = tmps;
				}
		if (unknownCount > 0)
		{
			dateCount.emplace_back(unknownCount);
			dateList.emplace_back("MetadataUnknown");
		}
		locaA = "" + dateList[0]; locaB = "" + to_string(dateCount[0]);
		for (j = 1; j < dateCount.size(); j++)
		{
			locaA += "," + dateList[j]; locaB += "," + to_string(dateCount[j]);
		}
		genoTypeList[i].representiveSeq = accessionID;
		genoTypeList[i].sampleCollectionDate.emplace_back(locaA);//基因型的采样时间 202003，202004，.......
		genoTypeList[i].sampleCollectionDate.emplace_back(locaB);//对应的序列数 7，6，2，1........
		//-------------------------------
		
		//------------------------------
		vector<string> lineageList; vector<int> lineageCount;//整理lineage，排序
		unknownCount = 0;
		for (j = 0; j < genoTypeList[i].sampleIDList.size(); j++)
		{
			auto val = std::find(lineageList.begin(), lineageList.end(), genoTypeList[i].sampleMetadataList[j].lineage);
			if (val != lineageList.end())
				lineageCount[val - lineageList.begin()]++;
			else
			{
				if (genoTypeList[i].sampleMetadataList[j].lineage != "MetadataUnknown")
				{
					lineageList.emplace_back(genoTypeList[i].sampleMetadataList[j].lineage);
					lineageCount.emplace_back(1);
				}
				else
					unknownCount++;
			}
		}
		for (j = 0; j < lineageCount.size(); j++)
			for (k = j + 1; k < lineageCount.size(); k++)//按照含有序列的数目顺序排列
				if (lineageCount[j] < lineageCount[k])
				{
					tmpi = lineageCount[j]; lineageCount[j] = lineageCount[k]; lineageCount[k] = tmpi;
					tmps = lineageList[j]; lineageList[j] = lineageList[k]; lineageList[k] = tmps;
				}
		if (unknownCount > 0)
		{
			lineageCount.emplace_back(unknownCount);
			lineageList.emplace_back("MetadataUnknown");
		}
		locaA = "" + lineageList[0]; locaB = "" + to_string(lineageCount[0]);
		for (j = 1; j < lineageCount.size(); j++)
		{
			locaA += "," + lineageList[j]; locaB += "," + to_string(lineageCount[j]);
		}
		genoTypeList[i].sampleLineage.emplace_back(locaA);//基因型存在的lineage L1，L2，.......
		genoTypeList[i].sampleLineage.emplace_back(locaB);//对应的序列数 7，6，2，1........
		//-------------------------------
		if (QC == "PASS")//如果基因型通过了QC
		{
			mtx.lock();//互斥锁上锁
			qcpassedCount++;
			if (lineageMap.find(lineageList[0]) != lineageMap.end())//该lineage存在于列表中，添加genotype到lineage下面
			{
				totalLineageList[lineageMap[lineageList[0]]].genotypeList.insert(pair<string, int>(genoTypeList[i].sampleCollectionDate[0]+to_string(i), i));
			}
			else
			{
				lineageMap.insert(pair<string, int>(lineageList[0], totalLineageList.size()));//lineage不存在，添加lineage
				Lineage a;
				a.lineageName = lineageList[0];
				a.genotypeList.insert(pair<string,int>(genoTypeList[i].sampleCollectionDate[0] + to_string(i),i));
				totalLineageList.push_back(a);
			}
			//查询突变是否有indel，重编码
			for (auto val = genoTypeList[i].mutList.begin(); val != genoTypeList[i].mutList.end(); val++)
			{
				if (val->second == 0)
					if (indelMap.find(val->first) == indelMap.end())
					{
						indelMap.insert(pair<string, string>(val->first, IndelRecoding(val->first, i)));
						indelList.emplace_back(val->first);
					}
			}
			mtx.unlock();//互斥锁解锁
		}
	}
}


//Genotype找metadata输出
//output genotype with metadata
int FindMetadata(string outputfold)
{
	int i, j, k,m,n;
	ofstream write(outputfold + "/GenotypeDetailTotal.tsv");
	ofstream writeglcm(outputfold + "/GenotypeLocationCM.tsv");
	//mkdir((outputfold + "/GenotypeDetailFile").c_str(), S_IRWXU | S_IROTH);//创建文件夹
	string SSDfold = "/mnt/m";
	mkdir((SSDfold + "/GenotypeDetailFile").c_str(), S_IRWXU | S_IROTH);//创建文件夹
	write << "GenotypeID\tRepresentiveSeq\tGenotypeMutation\tGenotypeSamples\tLocation\tLocationCounts\tLineage\tLineageCounts\tCollectionMonths\tCollectionMonthsCounts\tHighCompleteHuman" << endl;
	writeglcm << "GenotypeID\tLocation\tCollectionMonths\tCollectionMonthsCounts" << endl;

	int threadnum = thread::hardware_concurrency();//查询可用线程数
	printf("Find Available Threads: %d\n", threadnum);
	if (threadnum > 30)threadnum = 30;//设置线程数
	printf("Using Threads: %d\n\n", threadnum);
	thread threadlist[threadnum];
	k = genoTypeList.size() / (threadnum);
	for (i = 0; i < threadnum - 1; i++)
	{
		printf("Thread %d Start: Search from %d to %d\n", i, i * k, (i * k + k - 1));
		threadlist[i] = thread(ThreadSearchMetadata, i * k, (i * k + k - 1));//创建子线程
	}
	printf("Thread %d Start: Search from %d to %d\n", i, i * k, (genoTypeList.size() - 1));
	threadlist[i] = thread(ThreadSearchMetadata, i * k, (genoTypeList.size() - 1));//创建子线程
	for (i = 0; i < threadnum; i++)
	{
		threadlist[i].join();
		printf("Thread %d End, Time: %ld / %ld Seconds\n", i, time(NULL), (time(NULL) - starttime));
	}
	printf("\nStart Writing Genotypes\n");
	//return 0;//!!!
	for (i = 0; i < genoTypeList.size(); i++)
	{

		//-----------输出GenotypeDetailTotal-----------
		string output = "";
		string mutationOutput = "";
		output += to_string(i) + "\t";
		output += genoTypeList[i].representiveSeq + "\t";
		for (auto val = genoTypeList[i].mutList.begin(); val != genoTypeList[i].mutList.end(); val++)
			mutationOutput += val->first + " ";
		mutationOutput = mutationOutput.substr(0, mutationOutput.length() - 1);
		output += mutationOutput + "\t";
		output += to_string(genoTypeList[i].sampleIDList.size()) + "\t";
		output += genoTypeList[i].sampleLocation[0] + "\t" + genoTypeList[i].sampleLocation[1] + "\t";
		output += genoTypeList[i].sampleLineage[0] + "\t" + genoTypeList[i].sampleLineage[1] + "\t";
		output += genoTypeList[i].sampleCollectionDate[0] + "\t" + genoTypeList[i].sampleCollectionDate[1] + "\t";
		output += genoTypeList[i].highCompleteHuman;
		write << output << endl;
	}
	write.close();
	
	for (i = 0; i < genoTypeList.size(); i++)
	{
		//-----------输出GenotypeDetail-----------
		/*string mutationOutput = "";
		ofstream writeDetail(SSDfold + "/GenotypeDetailFile/GenotypeDetail_" + to_string(i) + ".Valk");
		if (!writeDetail.is_open())
		{
			cout << "Can not open " << SSDfold + "/GenotypeDetailFile/GenotypeDetail_" + to_string(i) + ".Valk" << endl;
		}
		writeDetail << "[GenotypeID(ChangesAveTimeUpd):]" << endl;
		writeDetail << to_string(i) << endl;
		writeDetail << "[Genotype Mutation:]" << endl;
		writeDetail << mutationOutput << endl;
		writeDetail << "[Genotype Contained Samples:]" << endl;
		writeDetail << genoTypeList[i].sampleIDList.size() << endl;
		writeDetail << "[Location & Counts:]" << endl;
		writeDetail << genoTypeList[i].sampleLocation[0] << endl;
		writeDetail << genoTypeList[i].sampleLocation[1] << endl;
		writeDetail << "[Lineage & Counts:]" << endl;
		writeDetail << genoTypeList[i].sampleLineage[0] << endl;
		writeDetail << genoTypeList[i].sampleLineage[1] << endl;
		writeDetail << "[Collection Months & Counts:]" << endl;
		writeDetail << genoTypeList[i].sampleCollectionDate[0] << endl;
		writeDetail << genoTypeList[i].sampleCollectionDate[1] << endl;
		writeDetail << "---------------------------------------Sample-Metadata---------------------------------------" << endl;
		writeDetail << metatitle << endl;
		for (j = 0; j < genoTypeList[i].sampleMetadataList.size(); j++)
			writeDetail << genoTypeList[i].sampleMetadataList[j].originalMetadata << endl;
		writeDetail.close();*/
		//-----------输出GenotypeLocationCM-----------
		string output = "";
		vector<string> country = splitStr(genoTypeList[i].sampleLocation[0], ',');
		for (j = 0; j < country.size(); j++)
		{
			if (country[j] == "MetadataUnknown")break;
			output = to_string(i) + "\t";
			output += country[j] + "\t";
			vector<string> collectMonth = splitStr(genoTypeList[i].sampleCollectionDate[0], ',');
			vector<int> monthCount;
			for (k = 0; k < collectMonth.size(); k++)
				monthCount.emplace_back(0);
			for (k = 0; k < genoTypeList[i].sampleMetadataList.size(); k++)
			{
				if (genoTypeList[i].sampleMetadataList[k].country != country[j])continue;
				for (m = 0; m < collectMonth.size(); m++)
					if (genoTypeList[i].sampleMetadataList[k].collectionDate == collectMonth[m])
						monthCount[m]++;
			}
			string monthA = "";
			string monthB = "";
			for (k = 0; k < collectMonth.size(); k++)
			{
				if (collectMonth[k] == "MetadataUnknown")break;
				if (monthCount[k] > 0)
				{
					monthA += collectMonth[k] + ',';
					monthB += to_string(monthCount[k]) + ',';
				}
			}
			if (monthA != "")
			{
				monthA = monthA.substr(0, monthA.length() - 1);
				monthB = monthB.substr(0, monthB.length() - 1);
				writeglcm << output + monthA + "\t" + monthB << endl;
			}
		}

	}
	writeglcm.close();
}

//改掉序列名称中的非法字符
/*
输出突变文件中，序列名称全部用序号替代
*/

//输出Linux脚本
//output linux script
void WriteScript(string outputfold, int treeNum)
{
	ofstream write(outputfold + "/IQTREE_" + to_string(treeNum) + ".sh");
	write << "#PBS -q core24" <<endl;
	write << "#PBS -l mem=60gb,nodes=1:ppn=12,walltime=900:00:00" << endl;
	write << "#PBS -o Valk" << to_string(treeNum) << ".o" << endl;
	write << "#PBS -e Valk" << to_string(treeNum) << ".e" << endl;
	write << "#HSCHED -s Valk+Phylo+Recon" << endl;
	write << "#PPN limit 12" << endl;
	write << "#!/bin/bash" << endl;
	write << "thisfold=\"$( cd \"$( dirname \"${BASH_SOURCE[0]}\" )\" && pwd )\"" << endl;
	write << "mainDir=`dirname $thisfold`" << endl;
	write << "cd $thisfold" << endl;
	write << "$mainDir/Script/Mutation2Sequence.out $mainDir/Data/reference.fa $thisfold/inDel.info $thisfold/Tree_LineageCombined_"<< to_string(treeNum) <<".tsv $thisfold/Tree_LineageCombined_"<<to_string(treeNum)<<".fa" << endl;
	write << "$mainDir/Script/iqtree2 -redo -s Tree_LineageCombined_" << to_string(treeNum) << ".fa -m GTR+F -asr -nt AUTO -o Wuhan-Hu-1" << endl;
	write << "$mainDir/Script/FindHotspotBackmutation_IQTreeVersion.out $thisfold " << to_string(treeNum) << endl;
	write << "#rm Tree_LineageCombined_" << to_string(treeNum) << ".fa.state" << endl;
	write << "rm Tree_LineageCombined_" << to_string(treeNum) << ".fa.mldist" << endl;
	write << "#rm Tree_LineageCombined_" << to_string(treeNum) << ".fa" << endl;
	write.close();
	return;
}

//多线程输出lineage画树文件
//output files for phylogenetic reconstruction
void WriteLineageFile(string outputfold, vector<int> linList, int treeNum)
{
	unordered_map<int, int> valMap;
	ofstream write(outputfold + "/Tree_LineageCombined_" + to_string(treeNum) + ".tsv");
	int i, j, k,m;
	if (linList.size() > 1)//计算采样频率
		k = 1;
	else
	{
		//k = 1;//不抽样
		k = totalLineageList[linList[0]].genotypeList.size() / 10000 + 1;//!!!
	}
	for (i = 0; i < linList.size(); i++)
	{
		m = 0;
		for (auto val = totalLineageList[linList[i]].genotypeList.begin(); 
			val != totalLineageList[linList[i]].genotypeList.end(); val ++,m++)
			if(m%k==0)
			{
				if (valMap.find(val->second) != valMap.end())
					continue;
				string output = "VAL_" + to_string(val->second) + "\t";
				for (auto val1 = genoTypeList[val->second].mutList.begin();
					val1 != genoTypeList[val->second].mutList.end(); val1++)
					if (indelMap.find(val1->first) != indelMap.end())
						output += indelMap[val1->first] + " ";
					else output += val1->first + " ";
				write << output.substr(0, output.length() - 1) << endl;
				valMap.insert(pair<int, int>(val->second, 1));
			}
	}
	WriteScript(outputfold, treeNum);
	write.close();
	return;
}

//合并lineage准备画树
//lineage merge
int TreePrepare(string outputfold)
{
	//先输出inDel.info
	ofstream writeindel(outputfold + "/inDel.info");
	int i, j, k;
	writeindel << "OriginalInDel\tRecodingInDel" << endl;
	for (i = 0; i < indelList.size(); i++)
		writeindel << indelList[i] << "\t" << indelMap[indelList[i]] << endl;
	writeindel.close();

	//输出各个树
	Lineage tmpl;
	for (i = 0; i < totalLineageList.size(); i++)//按序列数降序排列lineage
		for(j=i+1;j<totalLineageList.size();j++)
			if (totalLineageList[i].genotypeList.size() < totalLineageList[j].genotypeList.size())
			{
				tmpl = totalLineageList[i]; totalLineageList[i] = totalLineageList[j]; totalLineageList[j] = tmpl;
			}
	k = 0; int treeCount = 1;
	ofstream writeDetail(outputfold + "/TreeFileLineageDetail.tsv");
	vector<int> lin2Combine;
	vector<thread> outputThreads;
	for (i = 0; i < totalLineageList.size(); i++)//遍历lineage，合并，调用多线程函数输出
	{
		if (totalLineageList[i].genotypeList.size() + k > 10000)
		{
			lin2Combine.push_back(i);
			string output = "Tree_" + to_string(treeCount) + "\t";
			output += to_string(totalLineageList[i].genotypeList.size() + k);
			for (j = 0; j < lin2Combine.size(); j++)
				output += "\t" + totalLineageList[lin2Combine[j]].lineageName;
			writeDetail << output << endl;
			cout << "Writing Tree_" << treeCount << endl;
			outputThreads.push_back(thread(WriteLineageFile,outputfold,lin2Combine,treeCount));
			k = 0; treeCount++;
			lin2Combine.clear();
			if (outputThreads.size() > 19)
			{
				for (j = 0; j < outputThreads.size(); j++)
					outputThreads[j].join();
				outputThreads.clear();
			}
		}
		else
		{
			k += totalLineageList[i].genotypeList.size();
			lin2Combine.push_back(i);
		}
	}
	string output = "Tree_" + to_string(treeCount) + "\t";
	output += to_string(totalLineageList[i].genotypeList.size() + k);
	for (j = 0; j < lin2Combine.size(); j++)
		output += "\t" + totalLineageList[lin2Combine[j]].lineageName;
	writeDetail << output << endl;
	cout << "Writing Tree_" << treeCount << endl;
	outputThreads.push_back(thread(WriteLineageFile, outputfold, lin2Combine, treeCount));
	for (j = 0; j < outputThreads.size(); j++)
		outputThreads[j].join();

	writeDetail.close();
	return 0;
}

//按lineage亲缘关系合并lineage准备画树
//merge lineage by relativeness
int TreePrepareR(string outputfold)
{
	//先输出inDel.info
	ofstream writeindel(outputfold + "/inDel.info");
	int i, j, k;
	writeindel << "OriginalInDel\tRecodingInDel" << endl;
	for (i = 0; i < indelList.size(); i++)
		writeindel << indelList[i] << "\t" << indelMap[indelList[i]] << endl;
	writeindel.close();

	//输出各个树
	Lineage tmpl;
	for (i = 0; i < totalLineageList.size(); i++)//按序列数降序排列lineage
		for (j = i + 1; j < totalLineageList.size(); j++)
			if (totalLineageList[i].genotypeList.size() < totalLineageList[j].genotypeList.size())
			{
				tmpl = totalLineageList[i]; totalLineageList[i] = totalLineageList[j]; totalLineageList[j] = tmpl;
			}
	k = 0; int treeCount = 1;
	ofstream writeDetail(outputfold + "/TreeFileLineageDetail.tsv");
	vector<int> lin2Combine;
	vector<thread> outputThreads;
	for (i = 0; i < totalLineageList.size(); i++)//先输出单个序列数大于【】的lineage，调用多线程函数输出
	{
		if (totalLineageList[i].genotypeList.size() + k > 10000)
		{
			lin2Combine.push_back(i);
			string output = "Tree_" + to_string(treeCount) + "\t";
			output += to_string(totalLineageList[i].genotypeList.size() + k);
			for (j = 0; j < lin2Combine.size(); j++)
				output += "\t" + totalLineageList[lin2Combine[j]].lineageName;
			writeDetail << output << endl;
			cout << "Writing Tree_" << treeCount << endl;
			outputThreads.push_back(thread(WriteLineageFile, outputfold, lin2Combine, treeCount));
			k = 0; treeCount++;
			lin2Combine.clear();
			if (outputThreads.size() > 19)
			{
				for (j = 0; j < outputThreads.size(); j++)
					outputThreads[j].join();
				outputThreads.clear();
			}
		}
		else
		{
			break;
		}
	}
	int smallLineageStart = i;
	for (; i < totalLineageList.size(); i++)//按lineage名称排列需要合并的lineage
		for (j = i + 1; j < totalLineageList.size(); j++)
			if (strcmp(totalLineageList[i].lineageName.c_str(), totalLineageList[j].lineageName.c_str())>0)
			{
				tmpl = totalLineageList[i]; totalLineageList[i] = totalLineageList[j]; totalLineageList[j] = tmpl;
			}
	string linK = "ToBeDecide";
	int count = 0;
	lin2Combine.clear();
	vector<int> lineagesTooSmall;//不能独自画树的进化枝
	for (i = smallLineageStart; i < totalLineageList.size(); i++)
	{
		vector<string> name1 = splitStr(totalLineageList[i].lineageName, '.');
		if (name1[0] != linK)
		{
			if (lin2Combine.size() != 0 && count > 10000)//这一个大进化枝可以独自画树
			{
				vector<int> l2c;
				int countl2c = 0;
				for (j = 0; j < lin2Combine.size(); j++)//根据具体序列的多少拆成多个树
				{
					countl2c += totalLineageList[lin2Combine[j]].genotypeList.size();
					l2c.push_back(lin2Combine[j]);
					count -= totalLineageList[lin2Combine[j]].genotypeList.size();

					if (countl2c > 10000)//现有序列超过【】
					{
						if (count > 10000)//剩下的序列也超过【】，先输出现有的
						{
							string output = "Tree_" + to_string(treeCount) + "\t";
							output += to_string(countl2c);
							for (k = 0; k < l2c.size(); k++)
								output += "\t" + totalLineageList[l2c[k]].lineageName;
							writeDetail << output << endl;
							cout << "Writing Tree_" << treeCount << endl;
							outputThreads.push_back(thread(WriteLineageFile, outputfold, l2c, treeCount));
							countl2c = 0; treeCount++;
							l2c.clear();
						}
						else//剩下的序列没有超过【】，和现有的合并输出
						{
							for (k = j; k < lin2Combine.size(); k++)
								l2c.push_back(lin2Combine[k]);
							string output = "Tree_" + to_string(treeCount) + "\t";
							output += to_string(countl2c + count);
							for (k = 0; k < l2c.size(); k++)
								output += "\t" + totalLineageList[l2c[k]].lineageName;
							writeDetail << output << endl;
							cout << "Writing Tree_" << treeCount << endl;
							outputThreads.push_back(thread(WriteLineageFile, outputfold, l2c, treeCount));
							countl2c = 0; treeCount++;
							l2c.clear();
						}
					}
					
				}
			}
			else
			{
				for (j = 0; j < lin2Combine.size(); j++)
					lineagesTooSmall.push_back(lin2Combine[j]);
			}
			int count = 0;
			lin2Combine.clear();
			linK = name1[0];
		}
		count += totalLineageList[i].genotypeList.size();
		lin2Combine.push_back(i);
	}

	//处理最后的进化枝
	if (lin2Combine.size() != 0 && count > 10000)//这一个大进化枝可以独自画树
	{
		vector<int> l2c;
		int countl2c = 0;
		for (j = 0; j < lin2Combine.size(); j++)//根据具体序列的多少拆成多个树
		{
			countl2c += totalLineageList[lin2Combine[j]].genotypeList.size();
			l2c.push_back(lin2Combine[j]);
			count -= totalLineageList[lin2Combine[j]].genotypeList.size();

			if (countl2c > 10000)//现有序列超过【】
			{
				if (count > 10000)//剩下的序列也超过【】，先输出现有的
				{
					string output = "Tree_" + to_string(treeCount) + "\t";
					output += to_string(countl2c);
					for (k = 0; k < l2c.size(); k++)
						output += "\t" + totalLineageList[l2c[k]].lineageName;
					writeDetail << output << endl;
					cout << "Writing Tree_" << treeCount << endl;
					outputThreads.push_back(thread(WriteLineageFile, outputfold, l2c, treeCount));
					countl2c = 0; treeCount++;
					l2c.clear();
				}
				else//剩下的序列没有超过【】，和现有的合并输出
				{
					for (k = j; k < lin2Combine.size(); k++)
						l2c.push_back(lin2Combine[k]);
					string output = "Tree_" + to_string(treeCount) + "\t";
					output += to_string(countl2c + count);
					for (k = 0; k < l2c.size(); k++)
						output += "\t" + totalLineageList[l2c[k]].lineageName;
					writeDetail << output << endl;
					cout << "Writing Tree_" << treeCount << endl;
					outputThreads.push_back(thread(WriteLineageFile, outputfold, l2c, treeCount));
					countl2c = 0; treeCount++;
					l2c.clear();
				}
			}

		}
	}
	else
	{
		for (j = 0; j < lin2Combine.size(); j++)
			lineagesTooSmall.push_back(lin2Combine[j]);
	}

	//处理之前找出的不能独自画树的进化枝
	k = 0; lin2Combine.clear();
	for (i = 0; i < lineagesTooSmall.size(); i++)//先输出单个序列数大于【】的lineage，调用多线程函数输出
	{
		if (totalLineageList[lineagesTooSmall[i]].genotypeList.size() + k > 10000)
		{
			lin2Combine.push_back(lineagesTooSmall[i]);
			string output = "Tree_" + to_string(treeCount) + "\t";
			output += to_string(totalLineageList[lineagesTooSmall[i]].genotypeList.size() + k);
			for (j = 0; j < lin2Combine.size(); j++)
				output += "\t" + totalLineageList[lin2Combine[j]].lineageName;
			writeDetail << output << endl;
			cout << "Writing Tree_" << treeCount << endl;
			outputThreads.push_back(thread(WriteLineageFile, outputfold, lin2Combine, treeCount));
			k = 0; treeCount++;
			lin2Combine.clear();
			if (outputThreads.size() > 19)
			{
				for (j = 0; j < outputThreads.size(); j++)
					outputThreads[j].join();
				outputThreads.clear();
			}
		}
		else
		{
			lin2Combine.push_back(lineagesTooSmall[i]);
			k += totalLineageList[lineagesTooSmall[i]].genotypeList.size();
		}
	}

	string output = "Tree_" + to_string(treeCount) + "\t";
	output += to_string(totalLineageList[i].genotypeList.size() + k);
	for (j = 0; j < lin2Combine.size(); j++)
		output += "\t" + totalLineageList[lin2Combine[j]].lineageName;
	writeDetail << output << endl;
	cout << "Writing Tree_" << treeCount << endl;
	outputThreads.push_back(thread(WriteLineageFile, outputfold, lin2Combine, treeCount));
	for (j = 0; j < outputThreads.size(); j++)
		outputThreads[j].join();

	writeDetail.close();
	return 0;
}

//Calculate substitution matrix
void SubstitutionMatrixCalculation()//计算替换矩阵用于indel recoding
{
	int i, j, k;
	int AC = 0, AT = 0, AG = 0;
	int TC = 0, TA = 0, TG = 0;
	int CA = 0, CT = 0, CG = 0;
	int GA = 0, GT = 0, GC = 0;
	for (i = 0; i < genoTypeList.size(); i++)
	{
		for (auto val= genoTypeList[i].mutList.begin(); val != genoTypeList[i].mutList.end(); val++)
		{
			if (val->second == 1)//SNP
			{
				if (val->first[val->first.length() - 3] == 'A' && (val->first[val->first.length() - 1] == 'A' || val->first[val->first.length() - 1] == 'T' || val->first[val->first.length() - 1] == 'C' || val->first[val->first.length() - 1] == 'G'))
				{
					if (SNPrefA.find(val->first) == SNPrefA.end())
					{
						string a = "A/" + string(1,val->first[val->first.length() - 1]);
						SNPrefA.insert(pair<string, string>(val->first, a));
						if (val->first[val->first.length() - 1] == 'C')AC++;
						if (val->first[val->first.length() - 1] == 'T')AT++;
						if (val->first[val->first.length() - 1] == 'G')AG++;
					}
				}
				if (val->first[val->first.length() - 3] == 'T' && (val->first[val->first.length() - 1] == 'A' || val->first[val->first.length() - 1] == 'T' || val->first[val->first.length() - 1] == 'C' || val->first[val->first.length() - 1] == 'G'))
				{
					if (SNPrefT.find(val->first) == SNPrefT.end())
					{
						string a = "T/" + string(1, val->first[val->first.length() - 1]);
						SNPrefT.insert(pair<string, string>(val->first, a));
						if (val->first[val->first.length() - 1] == 'C')TC++;
						if (val->first[val->first.length() - 1] == 'A')TA++;
						if (val->first[val->first.length() - 1] == 'G')TG++;
					}
				}
				if (val->first[val->first.length() - 3] == 'C' && (val->first[val->first.length() - 1] == 'A' || val->first[val->first.length() - 1] == 'T' || val->first[val->first.length() - 1] == 'C' || val->first[val->first.length() - 1] == 'G'))
				{
					if (SNPrefC.find(val->first) == SNPrefC.end())
					{
						string a = "C/" + string(1, val->first[val->first.length() - 1]);
						SNPrefC.insert(pair<string, string>(val->first, a));
						if (val->first[val->first.length() - 1] == 'T')CT++;
						if (val->first[val->first.length() - 1] == 'A')CA++;
						if (val->first[val->first.length() - 1] == 'G')CG++;
					}
				}
				if (val->first[val->first.length() - 3] == 'G' && (val->first[val->first.length() - 1] == 'A' || val->first[val->first.length() - 1] == 'T' || val->first[val->first.length() - 1] == 'C' || val->first[val->first.length() - 1] == 'G'))
				{
					if (SNPrefG.find(val->first) == SNPrefG.end())
					{
						string a = "G/" + string(1, val->first[val->first.length() - 1]);
						SNPrefG.insert(pair<string, string>(val->first, a));
						if (val->first[val->first.length() - 1] == 'A')GA++;
						if (val->first[val->first.length() - 1] == 'T')GT++;
						if (val->first[val->first.length() - 1] == 'C')GC++;
					}
				}
			}
		}
	}
	int total = AC+AT+AG+TC+TA+TG+CA+CT+CG+GA+GT+GC;
	cout << "Substitution Matrix:" << endl;
	cout << "A/T: " << ((double)AT / total) << endl;
	cout << "A/C: " << ((double)AC / total) << endl;
	cout << "A/G: " << ((double)AG / total) << endl;
	cout << "T/A: " << ((double)TA / total) << endl;
	cout << "T/C: " << ((double)TC / total) << endl;
	cout << "T/G: " << ((double)TG / total) << endl;
	cout << "C/A: " << ((double)CA / total) << endl;
	cout << "C/T: " << ((double)CT / total) << endl;
	cout << "C/G: " << ((double)CG / total) << endl;
	cout << "G/A: " << ((double)GA / total) << endl;
	cout << "G/T: " << ((double)GT / total) << endl;
	cout << "G/C: " << ((double)GC / total) << endl;

	for (auto valA = SNPrefA.begin(); valA != SNPrefA.end(); valA++)
		randomSubstituteA.push_back(valA->second);
	for (auto valT = SNPrefT.begin(); valT != SNPrefT.end(); valT++)
		randomSubstituteT.push_back(valT->second);
	for (auto valC = SNPrefC.begin(); valC != SNPrefC.end(); valC++)
		randomSubstituteC.push_back(valC->second);
	for (auto valG = SNPrefG.begin(); valG != SNPrefG.end(); valG++)
		randomSubstituteG.push_back(valG->second);
	
	return;
}

int main(int argc, char* argv[])//三个输入 3 input 1.Sample_Mutlist 2.Metadata.tsv 3.Outputfold
{
	if (argc != 4)
	{
		printf("ERROR /// S2GC: 3 Input Requested: 1.Sample_Mutlist 2.Metadata.tsv 3.Outputfold ///\n");
		return -1;
	}
	starttime = time(NULL);//记录程序开始时间
	printf("\nStart Time: %ld\n", starttime);
	time_t now = time(0);
	cout << ctime(&now) << endl;
	string SMLfile = argv[1];
	int i, j, k;

	//序列去重复，归纳为GenoType
	if (SeqDedup2Genotype(SMLfile) == -1)
		return 1;

	printf("\nDedup Finished\nTotal Genotype Number: %ld\nEnd Time: %ld\nTotal Time Used: %ld Seconds\n",genoTypeList.size(), time(NULL), (time(NULL) - starttime));
	
	//计算替换矩阵用于indel recoding
	SubstitutionMatrixCalculation();
	printf("\nSubstitution Matrix Finished\nTime Used: %ld Seconds\n", (time(NULL) - starttime));
	
	//读入metadata
	if (ReadMetadata(argv[2]) == -1)
		return 1;

	printf("\nMetadata Input Finished\nTotal Metadata Number: %ld\nEnd Time: %ld\nTime Used: %ld Seconds\n", metadataMap.size(), time(NULL), (time(NULL) - starttime));

	//Genotype找metadata输出
	if (FindMetadata(argv[3]) == -1)
		return 1;

	cout << "\nGenotypePassedQC: " << qcpassedCount << endl;
	cout << "Writing Trees" << endl;

	//合并lineage准备画树
	if (TreePrepareR(argv[3]) == -1)
		return 1;

	printf("\nDone: \nEnd Time: %ld\nTotal Time Used: %ld Seconds\n", time(NULL), (time(NULL) - starttime));

    return 0;
}