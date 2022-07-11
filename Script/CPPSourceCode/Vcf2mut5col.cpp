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
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;

//string split函数
/*vector<string> splitStr(string str, char delimiter)
{
	vector<string> r;
	string tmpstr;
	while (!str.empty()) {
		int ind = str.find_first_of(delimiter);
		if (ind == -1) {
			r.emplace_back(str);
			str.clear();
		}
		else {
			r.emplace_back(str.substr(0, ind));
			str = str.substr(ind + 1, str.size() - ind - 1);
		}
	}
	return r;
}*/
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

 static string outputfold;
 std::mutex mtx;

class Mut
 {
 public: string POS;
	  string REF;
	  string ALT;
 };
class Sample
{
public: string SampleName;
	  vector<Mut> MutationList;
};
unordered_map<string, int>sampleMap;//sample名称和其在list中位置的hash表
vector<Sample>sampleList;

Mut VariationChop(Mut a)//砍掉重复snp 两端的尾巴 241CG/TG -> 241C/T
{
	int i, j;
	i = a.REF.length()-1;
	j = a.ALT.length() - 1;
	while (i > 0 && j > 0 && a.REF[i] == a.ALT[j])
	{
		i--; j--;
	}
	a.REF = a.REF.substr(0, i + 1);
	a.ALT = a.ALT.substr(0, j + 1);
	return a;
}

 void LineProcess(string line)
 {
	 vector<string> line1 = splitStr(line, '\t');//CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT  SAMPLE
	 vector<string> alt = splitStr(line1[4], ',');
	 int i, j, k;
	 Mut a;
	 a.POS = line1[1];
	 /*if (line1[1] == "28881" || line1[1] == "28882" || line1[1] == "28883")
	 {
		 cout << line1[0] << endl;
		 cout << line1[1] << endl;
		 cout << line1[2] << endl;
		 cout << line1[3] << endl;
		 cout << line1[4] << endl;
		 cout << line1[5] << endl;
		 cout << line1[6] << endl;
		 cout << line1[7] << endl;
		 cout << line1[8] << endl;
	 }*/
	 for (i = 9; i < line1.size(); i++)
	 {
		 if (line1[i] != "0" && line1[i] != ".")
		 {
			 a.REF = line1[3];
			 a.ALT = alt[stoi(line1[i]) - 1];
			 if (a.REF.length() > 1 && a.ALT.length() > 1)
				 a = VariationChop(a);
			 mtx.lock();
			 sampleList[i - 9].MutationList.emplace_back(a);
			 mtx.unlock();
		 }
	 }
	 return;
 }

int main(int argc, char* argv[])//2个输入 vcf文件路径 输出文件路径
{
	cout << "VCF2mut5col Start" << endl;
	int i, j, k;
    if (argc != 3)
    {
        printf("ERROR /// V2M: 2 Input Requested: 1.VCF file 2.Outputfold ///\n");
        return -1;
    }
	ifstream read(argv[1]);
	if (!read.is_open())
	{
		printf("ERROR /// V2M: Can Not Open VCF File\n");
		return -1;
	}
	outputfold = argv[2];
	outputfold += "/mut5col.BIGDver.tsv";
	
	ofstream write(outputfold);
	if (!write.is_open())
	{
		printf("ERROR /// V2M: Can Not Open Output File\n");
		return -1;
	}
	string line;
	while (getline(read, line))
	{
		vector<string> line1 = splitStr(line,'\t');
		if (line1[0] == "##VirusID")
		{
			for (i = 9; i < line1.size(); i++)
			{
				Sample newa;
				newa.SampleName = line1[i];
				sampleList.emplace_back(newa);
			}
			break;
		}
	}
	getline(read, line);
	vector<thread> threadList;
	while (getline(read, line))
	{
		threadList.emplace_back(thread(LineProcess, line));
		if (threadList.size() > 20)
		{
			for (i = 0; i < threadList.size(); i++)
				threadList[i].join();
			threadList.clear();
		}
	} 
	for (i = 0; i < threadList.size(); i++)
		threadList[i].join();
	cout << "Start Writting Results" << endl;
	for (i = 0; i < sampleList.size(); i++)
	{
		
		for (j = 0; j < sampleList[i].MutationList.size(); j++)
		{	
			string output = "BIGD\t";
			output += sampleList[i].MutationList[j].POS + "\t";
			output += to_string(stoi(sampleList[i].MutationList[j].POS)+sampleList[i].MutationList[j].REF.length() -1) + "\t";
			output += sampleList[i].MutationList[j].REF + "\t";
			output += sampleList[i].MutationList[j].ALT + "\t";
			output += sampleList[i].SampleName;
			write << output << endl;
		}
		//cout << i << endl;
	}

	write.close();
    return 0;
}