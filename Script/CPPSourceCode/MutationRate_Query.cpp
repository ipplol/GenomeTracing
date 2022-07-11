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
#include <mutex>
#include <pthread.h>
#include <map>
#include <unordered_map>
#include <time.h>
#include <math.h>
using namespace std;

//string split函数
vector<string> splitStr(string str, char delimiter)
{
    vector<string> r;
    string tmpstr;
    while (!str.empty()) {
        int ind = str.find_first_of(delimiter);
        if (ind == -1) {
            r.push_back(str);
            str.clear();
        }
        else {
            r.push_back(str.substr(0, ind));
            str = str.substr(ind + 1, str.size() - ind - 1);
        }
    }
    return r;
}

int MutRateCal(string targetfold)//构建依赖项，计算突变率
{
    string line;
    int i, j, k;
    ifstream readtree(targetfold + "/TreeFileLineageDetail.tsv");
    if (!readtree.is_open())cout << "Can Not Open " << targetfold + "/TreeFileLineageDetail.tsv" << endl;
    int totalseq = 0;
    while (getline(readtree,line))
    {
        vector<string> line1 = splitStr(line, '\t');
        totalseq += stoi(line1[1]);
    }
    readtree.close();
    ifstream read(targetfold + "/Tree_LineageCombined_Total.hotspot");
    if (!read.is_open())cout << "Can Not Open " << targetfold + "/Tree_LineageCombined_Total.hotspot" << endl;
    ofstream writehotspot(targetfold + "/TotalMutationRate.txt");
    writehotspot << "Mutation\tWeight\tMutationRate" << endl;
    getline(read, line);
    while (getline(read,line))
    {
        vector<string> line1 = splitStr(line, '\t');
        double Mutrate = stof(line1[1]) / (totalseq - stoi(line1[2]));
        writehotspot << line1[0] << "\t" << line1[1] << "\t" << Mutrate << endl;
    }
    read.close();
    writehotspot.close();

    ifstream readback(targetfold + "/Tree_LineageCombined_Total.backmutation");
    if (!readback.is_open())cout << "Can Not Open " << targetfold + "/Tree_LineageCombined_Total.backmutation" << endl;
    ofstream writeback(targetfold + "/TotalBackmutRate.txt");
    writeback << "Mutation\tBackMutationRate" << endl;
    getline(readback, line);
    while (getline(readback,line))
    {
        vector<string> line1 = splitStr(line, '\t');
        double backmutrate = stof(line1[1]) / stof(line1[2]);
        if (stof(line1[2]) > 50)
            writeback << line1[0] << "\t" << backmutrate << endl;
    }
    readback.close();
    writeback.close();
}

int main(int argc, char* argv[])//三个输入 1.目标文件夹 2.待查询突变 3.输出文件
{
    printf("MutationRate_Query\n");
    if (argc < 4)
    {
        printf("MutationRate_Query Needs 3 input: 1. Target fold, 2. Target Mutation, 3.Output File\n");
        return 1;
    }
    string targetfold = argv[1];
    string mut = argv[2];
    string outputfile = argv[3];
    ifstream readmut1(targetfold + "/TotalMutationRate.txt");
    if (!readmut1.is_open())
    {
        cout << "Can Not Open Mut Rate File, Strat Calculating" << endl;
        MutRateCal(targetfold);
    }
    readmut1.close();
    ifstream readmut(targetfold + "/TotalMutationRate.txt");
    string Weight = "1";
    string MutRate = "NA";
    string BackRate = "NA";
    string line;
    getline(readmut, line);
    while (getline(readmut,line))
    {
        if (line[1] == mut[1])
        {
            vector<string> line1 = splitStr(line, '\t');
            if (line1[0] == mut)
            {
                Weight = line1[1];
                MutRate = line1[2];
                break;
            }
        }
    }
    readmut.close();
    ifstream readback(targetfold + "/TotalBackmutRate.txt");
    getline(readback, line);
    while (getline(readback , line))
    {
        if (line[1] == mut[1])
        {
            vector<string> line1 = splitStr(line, '\t');
            if (line1[0] == mut)
            {
                BackRate = line1[1];
                break;
            }
        }
    }

    ofstream write(outputfile);
    write << "Mutation\tEvents\tMutationRate\tBackMutationRate" << endl;
    write << mut << "\t" << Weight << "\t" << MutRate << "\t" << BackRate << endl;
    write.close();
    readback.close();
    return 0;
}