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

std::mutex MergeCoOccurMTX;//互斥锁，保护同时出现次数的合并
static string workfold;//工作目录
static map<string,int> MutationList;//所有待计算连锁不平衡的突变位点
static unordered_map<string, int> MutationList_hash;//突变和其在列表中的位置hash map
static vector<vector<int>> CoOccurrenceMatrix;//突变同时出现的次数矩阵
class Genotype
{
public:
    vector<string> mutlist;
};
static vector<Genotype> GenotypeList;//所有基因型，被突变查询的序列

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
//读入数据
int ReadData()
{
    string line;
    //ifstream read(workfold + "/Custompoint.tsv");//读入自定义突变
    ifstream read(workfold + "/Tree_LineageCombined_Total.hotspot");//读入突变热点
    if (!read.is_open())
    {
        printf("ERROR /// LDRM: Can Not Open File %s /// \n", (workfold + "/Tree_LineageCombined_Total.hotspot").c_str());
        return 1;
    }
    getline(read, line);
    while (getline(read, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        if (stoi(line1[1]) >= 10 && MutationList.find(line1[0]) == MutationList.end())//防止重复且要求发生次数大于【】
        {
            MutationList.insert(pair<string, int>(line1[0], MutationList.size()));
            //MutationList_hash.insert(pair<string, int>(line1[0], MutationList.size() - 1));
        }
    }
    read.close();
    ifstream read2(workfold + "/Tree_LineageCombined_Total.backmutation");//读入回复突变
    if (!read2.is_open())
    {
        printf("ERROR /// LDRM: Can Not Open File %s /// \n", (workfold + "/Tree_LineageCombined_Total.backmutation").c_str());
        return -1;
    }
    getline(read2, line);
    while (getline(read2, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        if (stoi(line1[1]) >= 8 && MutationList.find(line1[0]) == MutationList.end())//防止重复且要求发生次数大于【】
        {
            MutationList.insert(pair<string, int>(line1[0], MutationList.size()));
            //MutationList_hash.insert(pair<string, int>(line1[0], MutationList.size() - 1));
        }
    }
    for (auto i = MutationList.begin(); i != MutationList.end(); i++)
        MutationList_hash.insert(pair<string,int>(i->first,MutationList_hash.size()));
    read2.close();
    cout << "Total mutation received: " << MutationList_hash.size() << endl;

    /*ifstream read3(workfold + "/GenotypeDetailTotal.tsv");//读入去重复后的基因型
    if (!read3.is_open())
    {
        printf("ERROR /// LDRM: Can Not Open File %s /// \n", (workfold + "/GenotypeDetailTotal.tsv").c_str());
        return -1;
    }
    getline(read3, line);
    while (getline(read3, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        Genotype a;
        a.mutlist = splitStr(line1[2], ' ');
        GenotypeList.push_back(a);
    }
    read3.close();*/

    ifstream read3(workfold + "/Tree_LineageCombined_Total.mutevent");//读入发生的总突变事件
    if (!read3.is_open())
    {
        printf("ERROR /// LDRM: Can Not Open File %s /// \n", (workfold + "/Tree_LineageCombined_Total.mutevent").c_str());
        return -1;
    }
    getline(read3, line);
    while (getline(read3, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        Genotype a;
        a.mutlist = splitStr(line1[1], ' ');
        GenotypeList.push_back(a);
    }
    read3.close();


    return 0;
}
//线程函数计算co-occurence
void MultiThreadCalCoO(int start, int end)
{
    int i, j, k;
    vector<int>a;
    vector<vector<int>> Thread_CoOccurrenceMatrix;//MapReduce 初始化该线程的突变同时出现的次数矩阵
    for (i = 0; i < MutationList_hash.size(); i++)
        a.push_back(0);
    for (i = 0; i < MutationList_hash.size(); i++)
        Thread_CoOccurrenceMatrix.push_back(a);
    for (i = start; i <= end; i++)//遍历GenotypeList
    {
        for (j = 0; j < GenotypeList[i].mutlist.size(); j++)
        {
            if (MutationList_hash.find(GenotypeList[i].mutlist[j]) == MutationList_hash.end())//这个突变不在查询突变列表中
                continue;
            Thread_CoOccurrenceMatrix[MutationList_hash[GenotypeList[i].mutlist[j]]][MutationList_hash[GenotypeList[i].mutlist[j]]]++;
            for (k = j + 1; k < GenotypeList[i].mutlist.size(); k++)
            {
                if (MutationList_hash.find(GenotypeList[i].mutlist[k]) == MutationList_hash.end())//这个突变不在查询突变列表中
                    continue;
                Thread_CoOccurrenceMatrix[MutationList_hash[GenotypeList[i].mutlist[j]]][MutationList_hash[GenotypeList[i].mutlist[k]]]++;
                Thread_CoOccurrenceMatrix[MutationList_hash[GenotypeList[i].mutlist[k]]][MutationList_hash[GenotypeList[i].mutlist[j]]]++;
            }
        }
    }
    
    MergeCoOccurMTX.lock();//互斥锁上锁
    for (i = 0; i < MutationList_hash.size(); i++)
        for (j = 0; j < MutationList_hash.size(); j++)
            CoOccurrenceMatrix[i][j] += Thread_CoOccurrenceMatrix[i][j];
    MergeCoOccurMTX.unlock();//互斥锁解锁
}
//多线程计算co-occurence
void CalCoO()
{
    int i, j, k;
    vector<int>a;//MapReduce 初始化全局突变同时出现的次数矩阵
    for (i = 0; i < MutationList_hash.size(); i++)
        a.push_back(0);
    for (i = 0; i < MutationList_hash.size(); i++)
        CoOccurrenceMatrix.push_back(a);

    int threadnum = thread::hardware_concurrency();//查询可用线程数
    printf("Find Available Threads: %d\n", threadnum);
    if (threadnum > 16)threadnum = 16;//设置线程数
    printf("Using Threads: %d\n\n", threadnum);

    //cout << GenotypeList.size() << endl;
    thread threadlist[threadnum];
    k = GenotypeList.size() / (threadnum);
    for (i = 0; i < threadnum - 1; i++)
    {
        printf("Thread %d Start: Search from %d to %d\n", i, i * k, (i * k + k - 1));
        threadlist[i] = thread(MultiThreadCalCoO, i * k, (i * k + k - 1));//创建子线程
    }
    printf("Thread %d Start: Search from %d to %d\n", i, i * k, (GenotypeList.size() - 1));
    threadlist[i] = thread(MultiThreadCalCoO, i * k, (GenotypeList.size() - 1));//创建子线程
    for (i = 0; i < threadnum; i++)
    {
        threadlist[i].join();
        printf("Thread %d End, Time: %ld \n", i, time(NULL));
    }
}
double ChiSquareTest( int a, int b, int c, int d)//卡方检验
{
    long n = a + b + c + d;
    long adbc = a * d - b * c;
    double k = n / (a + b);
    k = k * (adbc / (c + d));
    k = k * (adbc / (a + c));
    k = k / (b + d);
    return k;
}
double LDPvalue(double PMutMut, int MutAB)//计算连锁不平衡的P值,理论AB数目，理论非AB数目，实际AB，实际非AB
{
    int ABp = PMutMut * GenotypeList.size();
    int nABp = GenotypeList.size() - ABp;

    int nMutAB = GenotypeList.size() - MutAB;
    //cout << ABp << "\t" << nABp << "\t" << MutAB << "\t" << nMutAB << "\t" << endl;
    //cout << ChiSquareTest(ABp, nABp, MutAB, nMutAB) << endl;
    return ChiSquareTest(ABp, nABp, MutAB, nMutAB);
}
void CalLDR()//计算连锁不平衡R
{
    int i, j, k;
    vector<double> a;
    vector<vector<double>> LDRMartix;
    for (i = 0; i < MutationList_hash.size(); i++)
        a.push_back(-1);//-1代表没有共存在过
    for (i = 0; i < MutationList_hash.size(); i++)
        LDRMartix.push_back(a);
    ofstream writeP(workfold + "/LinkageDisequilibriumPair_Circox.txt");
    int linkpair = 1;
    ofstream write(workfold + "/LinkageDisequilibriumR.matrix");
    string output = to_string(MutationList_hash.size());
    if (!write.is_open())
    {
        cout << "ERROR /// Can Not Open OutputFile ///" << endl;
        return;
    }
    
    int totalgenotype = GenotypeList.size();
    map<int,int> LDsites;//记录达到R阈值的位点
    for (auto stli = MutationList.begin(); stli != MutationList.end(); stli++)
    {
        for (auto stlj = stli; stlj != MutationList.end(); stlj++)
        {
            if (stlj == stli)
            {
                LDRMartix[MutationList_hash[stlj->first]][MutationList_hash[stli->first]] = 1.0;
                continue;
            }
            int MutA = CoOccurrenceMatrix[MutationList_hash[stli->first]][MutationList_hash[stli->first]];//A出现的次数
            int MutB = CoOccurrenceMatrix[MutationList_hash[stlj->first]][MutationList_hash[stlj->first]];//B出现的次数
            int MutAB = CoOccurrenceMatrix[MutationList_hash[stli->first]][MutationList_hash[stlj->first]];//AB同时出现的次数
            double PRawRaw = (double)(totalgenotype - MutA - MutB + MutAB) / (totalgenotype);//频率P11
            double PMutMut = (double)(MutAB) / (totalgenotype);//频率P22
            double PRawMut = (double)(MutB - MutAB) / (totalgenotype);//频率P12
            double PMutRaw = (double)(MutA - MutAB) / (totalgenotype);//频率P21
            double D = (double)PRawRaw * PMutMut - PRawMut * PMutRaw;//D = P11*P22 - P12*P21
            double PA = (double)MutA / totalgenotype;
            double PB = (double)MutB / totalgenotype;
            double R = (double)D / sqrt(PA * (1 - PA) * PB * (1 - PB));
            double PAPB = ((double)MutA / (totalgenotype)) * ((double)MutB / (totalgenotype));//理论A·B
            double P1 = (double) MutA / (totalgenotype);//频率P1
            double P2 = 1 - P1;//频率P2
            double Q1 = (double) MutB / (totalgenotype);//频率Q1
            double Q2 = 1 - Q1;//频率Q2
            double Dmax = min(P1 * Q2, P2 * Q1);
            double Dprime = D / Dmax;
            if (R < 0 || (MutA + MutB + MutAB) < 10)R = 0.0;
            if (R >= 0.0 && MutAB >= 1)
            {
                LDRMartix[MutationList_hash[stli->first]][MutationList_hash[stlj->first]] = R;
                LDRMartix[MutationList_hash[stlj->first]][MutationList_hash[stli->first]] = R;
            }
            if (R >= 0.001 && MutAB >= 10)//0.001, 10
            {
                if(LDsites.find(MutationList_hash[stli->first])==LDsites.end())
                    LDsites.insert(pair<int,int>(MutationList_hash[stli->first],1));
                if (LDsites.find(MutationList_hash[stlj->first]) == LDsites.end())
                    LDsites.insert(pair<int, int>(MutationList_hash[stlj->first], 1));

                double Pvalue = LDPvalue(PAPB, MutAB);
                double ABp = ((double)MutA / (totalgenotype)) * ((double)MutB / (totalgenotype)) * GenotypeList.size();
                double nABp = GenotypeList.size() - ABp;

                int nMutAB = GenotypeList.size() - MutAB;
                //cout << R << endl;
                //输出格式是 Link#, 突变1, 突变2, R, Dprime, Pvalue
                writeP << "link" << to_string(linkpair) << "\t" << stli->first << "\t" << stlj->first << "\t" << R << "\t" << Dprime << "\t";
                writeP << ABp << "\t" << nABp << "\t" << MutAB << "\t" << nMutAB << "\t" << endl;
                linkpair++;
            }
        }
    }
    writeP.close();
    //以矩阵方式输出结果
    /*for (auto stli = MutationList.begin(); stli != MutationList.end(); stli++)
    {
        //if (LDsites.find(MutationList_hash[stli->first]) == LDsites.end())continue;
        output += "\t";
        output += stli->first;
    }
    write << output << endl;
    auto stli = MutationList.begin();
    for (i = 0; i < LDRMartix.size(); i++)
    {
        //if (LDsites.find(i) == LDsites.end()) { stli++; continue; }
        output = stli->first;
        for (j = 0; j < LDRMartix.size(); j++)
        {
            //if (LDsites.find(j) == LDsites.end())continue;
            output += "\t";
            if(LDRMartix[i][j]!=0)
                output += to_string(LDRMartix[i][j]);
            else
                output += "NA";
        }
        write << output << endl;
        stli++;
    }*/
    //以单列方式输出结果
    /*cout << "Start Writing Results, Size N*N: " << LDRMartix.size() << "\n" << endl;
    for (i = 0; i < LDRMartix.size(); i++)
        for (j = 0; j < LDRMartix.size(); j++)
            write << LDRMartix[i][j] << endl;*/
    write.close();
}
int main(int argc, char* argv[])//1个输入 Valkyrie/Linux的详细目录
{
    if (argc < 2)
    {
        printf("ERROR /// LDRM: Please Specify Working Fold Dir /// \n");
        return -1;
    }
    workfold = argv[1];
    //workfold = "/mnt/g/VariationMutation/回复突变/iqtreePopulation20210608/Linux";//debug
    time_t starttime = time(0);
    printf("%s Start!\n", "LinkageDisequilibriumRMatrix_MapReduce");
    cout << ctime(&starttime) << endl;
    
    
    if (ReadData() == -1)//读入数据
    {
        printf("ERROR /// LDRM: ReadData Failed /// \n");
        return 1;
    }
    
    CalCoO();//多线程计算co-occurence
    CalLDR();//计算连锁不平衡R

    time_t endtime = time(0);
    printf("%s Done!\n", "LinkageDisequilibriumRMatrix_MapReduce");
    cout << ctime(&endtime) << endl;

    return 0;
}