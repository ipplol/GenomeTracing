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

class Sample
{
public:
    string sampleName;
    vector<int> positionList;
    vector<string> mutation;
};

static vector<Sample> sampleList;

static unordered_map<string, int> mutationMap;//突变和计数，对mut5col去重复

void MutationSort(int ind)//从小到大排序突变
{
        int j, k;
        int i = ind;
        int tmpi; string tmps;
        for (j = 0; j < sampleList[i].positionList.size(); j++)
            for (k = j + 1; k < sampleList[i].positionList.size(); k++)
                if (sampleList[i].positionList[j] > sampleList[i].positionList[k])
                {
                    tmpi = sampleList[i].positionList[j];
                    sampleList[i].positionList[j] = sampleList[i].positionList[k];
                    sampleList[i].positionList[k] = tmpi;

                    tmps = sampleList[i].mutation[j];
                    sampleList[i].mutation[j] = sampleList[i].mutation[k];
                    sampleList[i].mutation[k] = tmps;
                }
        return;
}

int main(int argc, char* argv[])
{
    long starttime = time(NULL);//记录程序开始时间
    printf("\nStart Time: %ld\n", starttime);
    time_t now = time(0);
    cout << ctime(&now) << endl;

    if (argc != 2)
    {
        printf("ERROR /// M2SML: fold needed, e.g. </home/mawentai/APP/Valkyrie/LinuxNewUpdate/20210507> \n");
        return -1;
    }
    int i, j, k;
    string foldff = argv[1];
    ifstream read(foldff + "/mut5col.BIGDver.tsv");
    if (!read.is_open())
    {
        printf("ERROR /// Datainput: Can Not Open mutresult.tsv ///\n");
        cout << foldff + "/mut5col.BIGDver.tsv" << "\n" << endl;
        return -1;
    }
    ofstream write3(foldff + "/sample_mutlist.tsv");//可能缺少没有突变的序列信息！！！
    ofstream write1(foldff + "/mut_countlist.tsv");
    ofstream write2(foldff + "/Sample_list.tsv");
    string line;
    Sample a;
    string name = "START";
    vector<thread> threadList;
    while (getline(read,line))
    {
        vector<string> line1 =  splitStr( line,'\t');

        if (line1[4] == "R" && line1[3] != "-")//改写兼并碱基
            if (line1[3] == "A") line1[4] = "G";
            else line1[4] = "A";
        if (line1[4] == "Y" && line1[3] != "-")
            if (line1[3] == "C") line1[4] = "T";
            else line1[4] = "C";
        if (line1[4] == "M" && line1[3] != "-")
            if (line1[3] == "A") line1[4] = "C";
            else line1[4] = "A";
        if (line1[4] == "K" && line1[3] != "-")
            if (line1[3] == "G") line1[4] = "T";
            else line1[4] = "G";
        if (line1[4] == "S" && line1[3] != "-")
            if (line1[3] == "G") line1[4] = "C";
            else line1[4] = "G";
        if (line1[4] == "W" && line1[3] != "-")
            if (line1[3] == "A") line1[4] = "T";
            else line1[4] = "A";

        string mutation = line1[0] + "\t" + line1[1] + "\t" + line1[2] + "\t" + line1[3] + "\t" + line1[4];
        if (mutationMap.find(mutation) != mutationMap.end())//dedup mutation
        {
            mutationMap[mutation]++;
        }
        else
        {
            mutationMap.insert(pair<string, int>(mutation, 1));
        }

        if (line1[5] != name)
        {
            if (name != "START")
            {
                sampleList.emplace_back(a);
                int index = sampleList.size() - 1;
                threadList.emplace_back(thread(MutationSort, sampleList.size() - 1));
                if (threadList.size() > 20)
                {
                    for (i = 0; i < threadList.size(); i++)
                    {
                        threadList[i].join();
                    }
                    threadList.clear();
                }
            }
            name = line1[5];
            a.positionList.clear();
            a.mutation.clear();
            a.sampleName = line1[5];
            a.positionList.emplace_back(stoi(line1[1]));
            a.mutation.emplace_back(line1[3] + "/" + line1[4]);
            //transform(a.mutation.begin(), a.mutation.end(), a.mutation.begin(), ::toupper);//全部转大写
        }
        else

        {
            a.positionList.emplace_back(stoi(line1[1]));
            a.mutation.emplace_back(line1[3] + "/" + line1[4]);
            //transform(a.mutation.begin(), a.mutation.end(), a.mutation.begin(), ::toupper);//全部转大写
        }
    }
    for (i = 0; i < threadList.size(); i++)
        threadList[i].join();
    for (i = 0; i < sampleList.size(); i++)
    {
        
        string output = sampleList[i].sampleName + "\t";
        for (j = 0; j < sampleList[i].positionList.size(); j++)
            output += to_string(sampleList[i].positionList[j]) + sampleList[i].mutation[j] + " ";
        write3 << output.substr(0, output.length() - 1) << endl;
        write2 << sampleList[i].sampleName << endl;
    }
    for (auto val = mutationMap.begin(); val != mutationMap.end(); val++)
        write1 << val->first + "\t" + to_string(val->second) << endl;


    write3.close();
    write1.close();
    return 0;
}
