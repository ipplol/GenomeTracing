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
#include <mutex>
#include <map>
#include <unordered_map>
#include <time.h>
#include <math.h>
#include <algorithm>
using namespace std;

class Tree
{
    public:
        string TreeName;
        unordered_map<string, double> Mutation_Weight;//lineage突变和权重对应表 Lineage mutation and weight correspondence table
};

class Lineage
{
    public: 
         string LineageName;//lineage的名称 lineage name
         int ReferenceWeight;//lineage所具有突变的参考ref权重 Reference weight of mutation of lineage
         int Mutation_Weight_Index;//lineage突变和权重对应的哪一棵树 Which tree does lineage mutation and weight correspond to
          
};

class GenoType
{
    public: 
         int GenoTypeID;//genotype id
         string representiveSeq;//代表性序列 representive Sequence
         vector<string> MutationList;//具有的突变 mutations
         unordered_map<string, int> LineageContained;//所用于计算突变率的进化枝 lineage belonging hash
         vector<string> LineageBelonged;//所处的进化枝 lineage belonging
         double HighestRank;//和Query的最高得分 highest score across lineage
         string HighestRankLineage;//最高得分时所处的Lineage 
         int SharedMutationNum = 0;//和Query共有的突变数量 number of mutations shared with query
};

static long starttime;//程序开始时间 Start time
static int resultTopN = 50;//取各线程前N个结果 Take the first n results of each thread
static vector<Tree> TreeList;//存放树上算出来的突变率权重
static vector<Lineage> LineageList;//存放Lineage Store the mutation rate weight calculated on the tree
static unordered_map<string, int> LineageListMap;//从Lineage名称到其在列表中位置的映射 Mapping from lineage name to its position in the list
static vector<GenoType> GenotypeList;//存放Genotype Store genotype
static vector<string> QueryMutation;//待查询序列的突变 Mutation of the sequence to be queried
static unordered_map<string, int> QueryMutationMap;//把查询序列转成hash Turn query sequence into hash

static vector<GenoType> Resultlist;//各线程的前100个结果 The first 100 results of each thread

static bool UseLinkageDisequilibrium = false;//是否使用连锁不平衡 Whether to use Linkage Disequilibrium
static bool UseLinkageDisequilibrium_RE = false;//是否使用连锁不平衡 Whether to use Linkage Disequilibrium
static bool UseMutationOnly = false;//是否只使用突变，不附加权重 Whether to only use mutation without additional weight
static bool UseLineageRate = false;//是否考虑进化枝，基因组背景值差异 Whether lineage and genome background value differences are considered
static bool UseSharedMutation = false;//只考虑共有突变的个数 Only the number of common mutations is considered

static vector<vector<double>> LDRMatrix;//存放突变连锁不平衡相关性R的矩阵 Matrix storing mutation linkage disequilibrium correlation R
static unordered_map<string, int> LDMutation_hash;//hash表，突变和其在LDR矩阵中的位置 Hash table, mutation and its position in LD_R matrix
static unordered_map<string, double> LDR_hash;//hash表，突变1突变2 和其的LDR值 Hash table, mutation 1, mutation 2 and its LD_R value

std::mutex mtx_GenotypeList;//读取时操作GenotypeList的线程锁 Thread lock of reading in GenotypeList
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

//读取GenotypeList的协程 
//Coroutine of read in
void GenotypeInputIE(vector<string> lineList)
{
    int i,j,k;
    vector<GenoType>thread_genotypeList;
    for (j = 0; j < lineList.size(); j++)
    {
        vector<string> line1 = splitStr(lineList[j], '\t');
        GenoType a;
        a.GenoTypeID = stoi(line1[0]);
        a.representiveSeq = line1[1];
        a.MutationList = splitStr(line1[2], ' ');
        vector<string> linelineage = splitStr(line1[6], ',');
        for (i = 0; i < linelineage.size(); i++)//这里用map把lineage name和其在list中的位置对应 Here, map the lineage name to its position in the list
        {
            a.LineageBelonged.emplace_back(linelineage[i]);
            if (linelineage[i] != "MetadataUnknown" && LineageListMap.find(linelineage[i]) != LineageListMap.end())
                a.LineageContained.insert(pair<string, int>(linelineage[i], LineageListMap[linelineage[i]]));
        }
        a.LineageContained.insert(pair<string, int>("Total", LineageListMap["Total"]));
        a.HighestRank = 0; a.HighestRankLineage = "Total";
        if (UseMutationOnly)a.HighestRank = -999;
        thread_genotypeList.emplace_back(a);
    }

    mtx_GenotypeList.lock();
    for (i = 0; i < thread_genotypeList.size(); i++)
        GenotypeList.emplace_back(thread_genotypeList[i]);
    mtx_GenotypeList.unlock();
    return;
}

//读入输入数据 返回值不为0代表数据读入出错
//Read the input data. If the return value is not 0, it means that there is an error
int Datainput(string fold)
{
    ifstream readGDT(fold + "/GenotypeDetailTotal.tsv");
    ifstream readTLC(fold + "/TreeFileLineageDetail.tsv");
    if(!readGDT.is_open())
    {printf("ERROR /// Datainput: Can Not Open %s/GenotypeDetailTotal.tsv ///\n", fold.c_str()); return -1;}
    if (!readTLC.is_open())
    {printf("ERROR /// Datainput: Can Not Open %s/TreeFileLineageDetail.tsv ///\n", fold.c_str()); return -1;}
    int i, j, k = 0;
    string line;
    Tree treea;
    TreeList.emplace_back(treea);//给tree0占个位置 Take a place for Tree 0

    //读入merge了的总total lineage信息 Read the total lineage information after merge
    ifstream readhotspot(fold + "/Tree_LineageCombined_Total.hotspot");//读入汇总的突变热点,要求降序排列 Read the summarized mutation hotspots, which are required to be arranged in descending order
    if (!readhotspot.is_open())
    {printf("ERROR /// Datainput: Can Not Open %s/Tree_LineageCombined_Total.hotspot ///\n", fold.c_str()); return -1;}
    string treeline;
    getline(readhotspot, treeline);
    getline(readhotspot, treeline);
    unordered_map<string, double> MutationAndWeightTotal;//暂时储存每棵树算出来的权重 Temporarily store the calculated weight of each tree
    vector<string> line2 = splitStr(treeline, '\t');
    MutationAndWeightTotal.insert(pair<string, double>(line2[0], 1));//先把第一行读进来存储最大突变权重 First read in the first line and store the maximum mutation weight
    int refweight = stoi(line2[1]) + 1;//the highest value plus 1 is taken as a reference value
    while (getline(readhotspot, treeline))
    {
        line2 = splitStr(treeline, '\t');
        if(UseMutationOnly == false)
            MutationAndWeightTotal.insert(pair<string, double>(line2[0], refweight - stoi(line2[1])));//转换成unordered_map存储后，原有排序被打乱 the original sorting is disrupted
        else
            MutationAndWeightTotal.insert(pair<string, double>(line2[0], 1));//不附加权重，所有突变的权重都是1 Without additional weight, the weight of all mutations is 1
            //the weight is obtained by subtracting the number of occurrences from the reference value
    }
    readhotspot.close();
    ifstream readbackmut(fold + "/Tree_LineageCombined_Total.backmutation");//读入回复突变,要求降序排列 Read in back mutation, requiring descending order
    if (!readbackmut.is_open())
    {printf("ERROR /// Datainput: Can Not Open %s/Tree_LineageCombined_Total.backmutation ///\n", fold.c_str()); return -1;}
    getline(readbackmut, treeline);
    while (getline(readbackmut, treeline))//读入回复突变数据，重新处理权重 最大扣分1% Read in the back mutation data and reprocess the weight with a maximum deduction of 1%
    {
        line2 = splitStr(treeline, '\t');
        if (MutationAndWeightTotal.find(line2[0]) != MutationAndWeightTotal.end())
        {
            double backmutrate = stof(line2[1]) / stof(line2[2]);
            if (stof(line2[2]) > 50)
            {
                if (backmutrate > 0.01)backmutrate = 0.01;
                if (UseMutationOnly) backmutrate = 0;//Mutation only
                MutationAndWeightTotal[line2[0]] = MutationAndWeightTotal[line2[0]] * (1 - backmutrate);//map[key]=value;
            }
        }
    }
    readbackmut.close();
    Tree treeb;
    treeb.TreeName = "Total";
    treeb.Mutation_Weight = MutationAndWeightTotal;
    TreeList.emplace_back(treeb);//添加新树 add new tree
    Lineage a; 
    a.LineageName = "Total";
    a.Mutation_Weight_Index = TreeList.size() - 1; 
    if (UseMutationOnly == false)
        a.ReferenceWeight = refweight;
    else
        a.ReferenceWeight = 2;
    LineageList.emplace_back(a); LineageListMap.insert(pair<string, int>("Total", LineageList.size() - 1));
    printf("Loading TotalLineage Finished, Time: %ld / %ld Seconds\n", time(NULL), (time(NULL) - starttime));


    //(根据打分方法选择是否考虑进化枝差异，不附加权重就不读入)
	//According to the scoring method, choose whether to consider the difference of lineages
    if (UseMutationOnly == false && UseLineageRate == true)
    {
        int refweighttotal = refweight;
        while (getline(readTLC, line))//读入各个tree 的 lineage信息 read in lineage information of each tree
        {
            vector<string> line1 = splitStr(line, '\t');
            if (line1.size() > 3)continue;//20220704只读入单独画树的lineage only single-lineage-tree will be read
            vector<string> line2 = splitStr(line1[0], '_');
            string treename = line2[1];//Tree_[1]
            ifstream readhotspot(fold + "/Tree_LineageCombined_" + treename + ".hotspot");//读入突变热点,要求降序排列 Read in mutation hotspots, descending order
            if (!readhotspot.is_open())
            {
                printf("ERROR /// Datainput: Can Not Open %s/Tree_LineageCombined_%s.hotspot ///\n", fold.c_str(), treename.c_str()); return -1;
            }
            string treeline;
            getline(readhotspot, treeline);
            getline(readhotspot, treeline);
            unordered_map<string, double> MutationAndWeight;//暂时储存每棵树算出来的权重
            line2 = splitStr(treeline, '\t');
            MutationAndWeight.insert(pair<string, double>(line2[0], 1));//先把第一行读进来存储最大突变权重
            int refweight = stoi(line2[1]) + 1;//the highest value plus 1 is taken as a reference value
            while (getline(readhotspot, treeline))
            {
                line2 = splitStr(treeline, '\t');
                MutationAndWeight.insert(pair<string, double>(line2[0], refweight - stoi(line2[1])));//转换成unordered_map存储后，原有排序被打乱
                //the weight is obtained by subtracting the number of occurrences from the reference value
            }
            readhotspot.close();
            ifstream readbackmut(fold + "/Tree_LineageCombined_" + treename + ".backmutation");//读入回复突变,要求降序排列
            if (!readbackmut.is_open())
            {
                printf("ERROR /// Datainput: Can Not Open %s/Tree_LineageCombined_%s.backmutation ///\n", fold.c_str(), treename.c_str()); return -1;
            }

            getline(readbackmut, treeline);
            while (getline(readbackmut, treeline))//读入回复突变数据，重新处理权重 最大扣分1%
            {
                line2 = splitStr(treeline, '\t');
                if (MutationAndWeight.find(line2[0]) != MutationAndWeight.end())
                {
                    double backmutrate = stof(line2[1]) / stof(line2[2]);
                    if (stof(line2[2]) > 50)
                    {
                        if (backmutrate > 0.01)backmutrate = 0.01;
                        MutationAndWeight[line2[0]] = MutationAndWeight[line2[0]] * (1 - backmutrate);//map[key]=value;
                    }
                }
            }
            readbackmut.close();
            Tree treeb;
            treeb.TreeName = treename;
            if (line1.size() == 3 || line1[3]=="")
                treeb.Mutation_Weight = MutationAndWeight;//添加新树
            else
            {
                treeb.Mutation_Weight = MutationAndWeightTotal;//每个树只能有一个lineage，合并的树按total算 Each tree can only have one lineage, and the merged trees are calculated by total
                refweight = refweighttotal;
            }
            TreeList.push_back(treeb);
            for (i = 2; i < line1.size(); i++)//依次添加新的lineage到list Add new lineages to the list in turn
            {
                Lineage a;
                a.LineageName = line1[i];
                a.Mutation_Weight_Index = TreeList.size() - 1;
                a.ReferenceWeight = refweight;
                LineageList.push_back(a);//加入列表
                LineageListMap.insert(pair<string, int>(line1[i], LineageList.size() - 1));//添加映射
            }
            //printf("TreeLineage %s Input Finish, Time: %ld / %ldSeconds\n", treename.c_str(), time(NULL), (time(NULL) - starttime));
        }
        printf("Loading TreeLineageFile Finished, Time: %ld / %ld Seconds\n", time(NULL), (time(NULL) - starttime));
    }
    //-----------------------------输出权重表------------------
    /*for (i = 0; i < TreeList.size(); i++)
    {
        ofstream writeweight(fold + "/WeightList/Tree_" + TreeList[i].TreeName + ".weightlist");
        for (auto ai = TreeList[i].Mutation_Weight.begin(); ai != TreeList[i].Mutation_Weight.end(); ai++)
            writeweight << ai->first << "\t" << ai->second << endl;
        writeweight.close();
    }*/

    //printf("Lineage Input Finish, Time: %ld / %ld Seconds\n", time(NULL), (time(NULL) - starttime));
    //---------------------------------------------------------

    vector<thread> threadList_GenotypeInputIE;
    vector<string> lineList;
    getline(readGDT, line);
    while (getline(readGDT, line))//读入Genotype信息 read in Genotype information
    {
        lineList.emplace_back(line);
        if (lineList.size() > 10000)
        {
            threadList_GenotypeInputIE.push_back(thread(GenotypeInputIE, lineList));
            lineList.clear();
        }
    }
    threadList_GenotypeInputIE.push_back(thread(GenotypeInputIE, lineList));
    for (i = 0; i < threadList_GenotypeInputIE.size(); i++)
        threadList_GenotypeInputIE[i].join();
    readGDT.close();
    printf("Loading GenotypeFile Finished, Time: %ld / %ld Seconds\n", time(NULL), (time(NULL) - starttime));

    //----------------------------------------------------------------------------
    if (UseLinkageDisequilibrium)//当使用连锁不平衡修正时才读入相关数据，数据较大 read in only when the linkage imbalance correction is used
    {
        printf("UseLinkageDisequilibrium True, Start Loading LDR Matrix\n");
        ifstream readLD(fold + "/LinkageDisequilibriumPair_Circox.txt");
        if (!readLD.is_open())
        {
            printf("ERROR /// Datainput: Can Not Open %s/LinkageDisequilibriumPair_Circox.txt ///\n", fold.c_str()); return -1;
        }
        while (getline(readLD, line))
        {
            vector<string> line1 = splitStr(line, '\t');
            LDR_hash.insert(pair<string, double>(line1[1] + line1[2], stof(line1[3])));
            LDR_hash.insert(pair<string, double>(line1[2] + line1[1], stof(line1[3])));
        }
        readLD.close();
        printf("Loading LD Finished, Time: %ld / %ld Seconds\n", time(NULL), (time(NULL) - starttime));
    }

    return 0;
}

//构建向量组，正交化，长度求和return
//Construct vector group, orthogonalization, return length summation
unordered_map<string, double> LDGramSchmidtProcess(vector<string> MutList, vector<double> WeightList)
{
    int i, j, k;
    double tmpd;
    string tmps;
    unordered_map<string, double> ResultHash;//突变和修改过的权重 Mutation and modified weight

    for (i = 0; i < WeightList.size(); i++)//降序排列 Descending order
        for (j = i + 1; j < WeightList.size(); j++)
            if (WeightList[i] < WeightList[j])
            {
                tmpd = WeightList[i]; WeightList[i] = WeightList[j]; WeightList[j] = tmpd;
                tmps = MutList[i]; MutList[i] = MutList[j]; MutList[j] = tmps;
            }

    vector<vector<long double>> VectorGroup_V;//向量组_V vector group V
    vector<vector<long double>> VectorGroup_X;//向量组_X vector group X
    vector<long double> a;//基底向量 base vector
    
    for (i = 0; i < MutList.size(); i++)//给基底向量赋值，有多少突变，向量就有多少维度 Assign a value to the base vector, and the vector has as many dimensions as there are mutations
        a.emplace_back(0.0);
    for (i = 0; i < MutList.size(); i++)
    {
        long double Vsquare = WeightList[i] * WeightList[i];//垂线Vi的长度的平方 Square of the length of the vertical line Vi
        vector<long double> VectorSum = a;//直角三角形的底边，垂线是Vi，斜边是Xi At the bottom of a right triangle, the vertical line is Vi and the hypotenuse is Xi
        for (j = 0; j < VectorGroup_X.size(); j++)
        {
            long double R;
            if (LDR_hash.find(MutList[i] + MutList[j]) != LDR_hash.end())
            {
                R = LDR_hash[MutList[i] + MutList[j]];
                //cout << R << endl;
            }
            else
                R = 0;
            long double factor = WeightList[i] * R / WeightList[j];
            for (k = 0; k < VectorSum.size(); k++)//每一维度乘上系数（内积，相关系数R）Multiply each dimension by the coefficient (inner product, correlation coefficient R)
            { 
                VectorSum[k] += VectorGroup_X[j][k] * factor;
            }
        }
        for (k = 0; k < VectorSum.size(); k++)
        {
            Vsquare -= VectorSum[k] * VectorSum[k];
            if (Vsquare < 0) { Vsquare = 1; break; }//给连锁不平衡矫正设个最低值 Set a minimum value for LD correction
        }
        vector<long double>VectorVi = a;
        VectorVi[i] = sqrt(Vsquare);//V2 = <0,|V2|,0,....>
        vector<long double>VectorXi = VectorVi;
        for (k = 0; k < VectorSum.size(); k++)
            VectorXi[k] += VectorSum[k];

        VectorGroup_V.emplace_back(VectorVi);
        VectorGroup_X.emplace_back(VectorXi);

        ResultHash.insert(pair<string, double>(MutList[i], sqrt(Vsquare)));
    }

    return ResultHash;
}

//向量组正交化修正连锁不平衡 【1】
//之前的方法，分别计算两个序列，所有突变的LD修改
double LDRCorrection(unordered_map<string,int>GQSharedMutation,vector<double> GenotypeMutWeight,vector<double> QueryMutWeight,int GenotypeID)
{
    //1.Genotype
    unordered_map<string, double> GenotypeHash = LDGramSchmidtProcess(GenotypeList[GenotypeID].MutationList, GenotypeMutWeight);
    //2.Query
    unordered_map<string, double> QueryHash = LDGramSchmidtProcess(QueryMutation, QueryMutWeight);

    double GQW1 = 0, GQW2 = 0;
    double GenotypeWeight = 0, QueryWeight = 0;

    for (auto i = GenotypeHash.begin(); i != GenotypeHash.end(); i++)
    {
        GenotypeWeight += i->second;
        if (GQSharedMutation.find(i->first) != GQSharedMutation.end())GQW1 += i->second;
    }

    for (auto j = QueryHash.begin(); j != QueryHash.end(); j++)
    {
        QueryWeight += j->second;
        if (GQSharedMutation.find(j->first) != GQSharedMutation.end())GQW2 += j->second;
    }

    return 0.5 * (GQW1 / GenotypeWeight) + 0.5 * (GQW2 / QueryWeight);
}
//向量组正交化修正连锁不平衡 【2】
//20210902新方法，分别只计算两个序列差异的突变位点的LD修改，共有的突变不改变
//Correction of linkage disequilibrium by orthogonalization of vector group
//Only the LD modification of the mutation sites with two sequence differences was calculated, and the common mutation did not change
double LDRCorrection_MutOnly(unordered_map<string, int>GQSharedMutation, vector<double> GenotypeMutWeight, vector<double> QueryMutWeight, int GenotypeID)
{
    int i,k; double GQW = 0;
    double GenotypeWeight = 0, QueryWeight = 0;
    unordered_map<string, double> mutationAndWeight;
    for (i = 0; i < GenotypeList[GenotypeID].MutationList.size(); i++)
        mutationAndWeight.insert(pair<string, double>(GenotypeList[GenotypeID].MutationList[i], GenotypeMutWeight[i]));
    k = 0;
    for (auto vali = QueryMutationMap.begin(); vali != QueryMutationMap.end(); vali++,k++)
        if (mutationAndWeight.find(vali->first) == mutationAndWeight.end())
        {
            mutationAndWeight.insert(pair<string, double>(vali->first, QueryMutWeight[k]));
        }
    
    for (auto j = GQSharedMutation.begin(); j != GQSharedMutation.end(); j++)
    {
        GQW += mutationAndWeight[j->first];
    }

    //1.Genotype
    vector<string> GenotypeDifferenceMut;
    vector<double> GenotypeDifferenceMutWeight;
    for (i = 0; i < GenotypeList[GenotypeID].MutationList.size(); i++)
    {
        if (GQSharedMutation.find(GenotypeList[GenotypeID].MutationList[i]) == GQSharedMutation.end())
        {
            GenotypeDifferenceMut.push_back(GenotypeList[GenotypeID].MutationList[i]);
            GenotypeDifferenceMutWeight.push_back(GenotypeMutWeight[i]);
        }
    }
    
    //2.Query
    vector<string> QueryDifferenceMut;
    vector<double> QueryDifferenceMutWeight;
    k = 0;
    for (auto vali = QueryMutationMap.begin(); vali != QueryMutationMap.end(); vali++,k++)
    {
        if (GQSharedMutation.find(vali->first) == GQSharedMutation.end())
        {
            QueryDifferenceMut.push_back(vali->first);
            QueryDifferenceMutWeight.push_back(QueryMutWeight[k]);
        }
    }

    unordered_map<string, double> GenotypeHash = LDGramSchmidtProcess(GenotypeDifferenceMut, GenotypeDifferenceMutWeight);
    unordered_map<string, double> QueryHash = LDGramSchmidtProcess(QueryDifferenceMut, QueryDifferenceMutWeight);

    for (auto j = GenotypeHash.begin(); j != GenotypeHash.end(); j++)
    {
        GenotypeWeight += j->second;
    }
    for (auto j = QueryHash.begin(); j != QueryHash.end(); j++)
    {
        QueryWeight += j->second;
    }
    return 0.5 * ((GQW / (GenotypeWeight + GQW)) + (GQW / (QueryWeight + GQW)));
}


//多线程搜索,每个线程分别搜索genotype list从start到end的位置
//Multi thread search, each thread searches the location of the genotype list from start to end
void ThreadSearch(int start, int end)
{
    int i, j, k;
    for (i = start; i <= end; i++)//GenotypeList
    {
        GenotypeList[i].HighestRank = 0;
        for (auto iterj = GenotypeList[i].LineageContained.begin(); iterj!=GenotypeList[i].LineageContained.end(); iterj++)//Lineage
        {
            k = iterj->second;//lineage index
            double GenotypeWeight = 0;//Genotype所有突变的权重 Weight of all genotype' mutations
            double QueryWeight = 0;//Query所有突变的权重 Weight of all query' mutations
            double GQWeight = 0;//两者共有的突变的权重 Shared mutation

            unordered_map<string,int>GQSharedMutation;//两者共有的突变 LD
            vector<double> GenotypeMutWeight;//Genotype突变的权重 LD
            vector<double> QueryMutWeight;//Query突变的权重 LD
            vector<double> GQMutWeight;//两者共有的突变的权重 LD

            for (auto iterq = QueryMutationMap.begin(); iterq != QueryMutationMap.end(); iterq++)//算Query所有突变的权重
            {
                auto iterl = TreeList[LineageList[k].Mutation_Weight_Index].Mutation_Weight.find(iterq->first);
                if (iterl != TreeList[LineageList[k].Mutation_Weight_Index].Mutation_Weight.end())
                {
                    QueryWeight += iterl->second;
                    QueryMutWeight.emplace_back(iterl->second);
                }
                else
                {
                    QueryWeight += LineageList[k].ReferenceWeight - 1;
                    QueryMutWeight.emplace_back(LineageList[k].ReferenceWeight - 1);
                }
            }

            for (j = 0; j < GenotypeList[i].MutationList.size(); j++)
            {
                double weight1 = 0;
                auto iterl = TreeList[LineageList[k].Mutation_Weight_Index].Mutation_Weight.find(GenotypeList[i].MutationList[j]);//算Genotype和GQ所有突变的权重
                if (iterl != TreeList[LineageList[k].Mutation_Weight_Index].Mutation_Weight.end())
                {
                    weight1 += iterl->second;
                    GenotypeMutWeight.emplace_back(iterl->second);
                }
                else
                {
                    weight1 += LineageList[k].ReferenceWeight - 1;
                    GenotypeMutWeight.emplace_back(LineageList[k].ReferenceWeight - 1);
                }
                GenotypeWeight += weight1;
                if (QueryMutationMap.find(GenotypeList[i].MutationList[j]) != QueryMutationMap.end())//该突变genotype query共有
                {
                    GQWeight += weight1;
                    GenotypeList[i].SharedMutationNum++;
                    GQMutWeight.emplace_back(weight1);
                    GQSharedMutation.insert(pair<string,int>(GenotypeList[i].MutationList[j],1));
                }
            }
            
            double RankLineage = 0.5 * (GQWeight / GenotypeWeight) + 0.5 * (GQWeight / QueryWeight);//当前lineage下的得分
            
            if (UseSharedMutation)//只考虑共有的突变数 Only the common mutation number is considered
            {
                RankLineage = GenotypeList[i].SharedMutationNum;
            }
            if (UseMutationOnly)//只考虑突变个数差异 Only consider the difference in the number of mutations
            {
                //20210513，MO方法改成全部使用相同的突变权重（1）
                // 因此注释了下面两行，以突变差异数代替打分的语句
                //RankLineage = 0.0 + (double)GQMutWeight.size() + (double)GQMutWeight.size();
                //RankLineage = RankLineage -(double)GenotypeMutWeight.size() - (double)QueryMutWeight.size();
            }
            if (UseLinkageDisequilibrium_RE)//考虑连锁不平衡的修正 Consider the correction of LD
            {
                //调用子线程修正
                // 20220514 LDCL面对日益增长的序列太慢了，决定只对WM或者CL排名靠前的结果进行LD重算
                // 分两次进行，先算CL，得分最高的排序到前50
                // 再对这前50重算一次,重算的时候再把UseLinkageDisequilibrium_RE调整为true
                //double LDCorrectRank = LDRCorrection(GQSharedMutation, GenotypeMutWeight, QueryMutWeight, i);//修正所有位点
                double LDCorrectRank = LDRCorrection_MutOnly(GQSharedMutation, GenotypeMutWeight, QueryMutWeight, i);//只修正差异位点
                RankLineage = LDCorrectRank;
            }
            

            if (RankLineage > GenotypeList[i].HighestRank)
            {
                GenotypeList[i].HighestRank = RankLineage;
                GenotypeList[i].HighestRankLineage = iterj->first;
            }
        }
        GenotypeList[i].SharedMutationNum /= GenotypeList[i].LineageContained.size();
    }
    //分线程Rank排序，只冒泡前50个 Rank sorting by thread
    GenoType tmp;
    for (i = start; i <= end && i < start + resultTopN; i++)
        for (j = i + 1; j <= end; j++)
        {
            if (GenotypeList[i].HighestRank < GenotypeList[j].HighestRank)
            {
                tmp = GenotypeList[i]; GenotypeList[i] = GenotypeList[j]; GenotypeList[j] = tmp;
            }
            else//当得分相同的情况下，按突变距离再次排序 When the scores are the same, sort again according to the mutation distance
            {
                if (GenotypeList[i].HighestRank == GenotypeList[j].HighestRank)
                {
                    int mutdi = abs((int)(GenotypeList[i].MutationList.size() + QueryMutation.size() - GenotypeList[i].SharedMutationNum * 2));
                    int mutdj = abs((int)(GenotypeList[j].MutationList.size() + QueryMutation.size() - GenotypeList[j].SharedMutationNum * 2));
                    if (mutdi > mutdj)
                    {
                        tmp = GenotypeList[i]; GenotypeList[i] = GenotypeList[j]; GenotypeList[j] = tmp;
                    }
                }
            }
            
        }
}

//开始搜索
//Start search
void GenotypeSearch()
{
    int i, j, k;
    for (i = 0; i < QueryMutation.size(); i++)
        QueryMutationMap.insert(pair<string, int>(QueryMutation[i], i+1));
    int threadnum = thread::hardware_concurrency();//查询可用线程数 Query the number of available threads
    printf("Find Available Threads: %d\n", threadnum);
    if(threadnum > 30)threadnum = 30;//设置线程数
    printf("Using Threads: %d\n\n", threadnum);
    if (UseLinkageDisequilibrium)
        printf("/// UseLinkageDisequilibrium True ///\n");

    thread threadlist[threadnum];
    k = GenotypeList.size() / (threadnum);
    for (i = 0; i < threadnum - 1; i++)
    {
        printf("Thread %d Start: Search from %d to %d\n", i, i*k, (i*k+k-1));
        threadlist[i] = thread(ThreadSearch, i * k, (i * k + k - 1));//创建子线程
    }
    printf("Thread %d Start: Search from %d to %d\n", i, i * k, (GenotypeList.size() - 1));
    threadlist[i] = thread(ThreadSearch, i * k, (GenotypeList.size() - 1));//创建子线程
    for (i = 0; i < threadnum; i++)
    {
        threadlist[i].join();
        printf("Thread %d End, Time: %ld / %ld Seconds\n", i, time(NULL), (time(NULL) - starttime));
    }
    if (UseLinkageDisequilibrium == true)//需要连锁不平衡调整 重跑各线程前50 If LD correction is required, rerun the top 50 result of each thread
    {
        thread threadlist2[threadnum];
        int TopN = resultTopN + 30;
        if (k - 1 < TopN)
            TopN = resultTopN + 30;
        printf("\nStart LD Adjustment _\n");
        UseLinkageDisequilibrium_RE = true;
        for (i = 0; i < threadnum - 1; i++)
        {
            printf("Thread %d Start: Search from %d to %d\n", i, i * k, (i * k + TopN));
            threadlist2[i] = thread(ThreadSearch, i * k, (i * k + TopN));//创建子线程
        }
        printf("Thread %d Start: Search from %d to %d\n", i, i * k, (i * k + TopN));
        if((i * k + TopN) < (GenotypeList.size() - 1))
            threadlist2[i] = thread(ThreadSearch, i * k, (i * k + TopN));//创建子线程
        else
            threadlist2[i] = thread(ThreadSearch, i * k, (GenotypeList.size() - 1));//创建子线程
        for (i = 0; i < threadnum; i++)
        {
            threadlist2[i].join();
            printf("Thread %d End, Time: %ld / %ld Seconds\n", i, time(NULL), (time(NULL) - starttime));
        }
    }
    printf("\nSearch Finished\nStart Sorting Results\n");
    //根据结果Rank排序 sort the results by Rank
    GenoType tmp;
    
    if (k > resultTopN + 30)
    {
        for (i = 0; i < threadnum; i++)//各线程结果取前50个 Take the top 50 results of each thread
        {
            for (j = i * k; j <= i * k + resultTopN + 30 && j <= (GenotypeList.size() - 1); j++)
                Resultlist.emplace_back(GenotypeList[j]);
        }
    }
    else
        for (i = 0; i < GenotypeList.size(); i++)//各线程结果取前50个 Take the top 50 results of each thread
        {
            Resultlist.emplace_back(GenotypeList[i]);
        }

    cout << "Resultlist Size" << endl;
    cout << Resultlist.size() << endl;
    for(i=0;i< resultTopN && i < Resultlist.size();i++)//冒泡排序 sort 
        for (j = i + 1; j < Resultlist.size(); j++)
        {
            if (Resultlist[i].HighestRank < Resultlist[j].HighestRank)
            {
                tmp = Resultlist[i]; Resultlist[i] = Resultlist[j]; Resultlist[j] = tmp;
            }
            else//当得分相同的情况下，按突变距离再次排序 When the scores are the same, sort again according to the mutation distance
            {
                if (Resultlist[i].HighestRank == Resultlist[j].HighestRank)
                {
                    int mutdi = abs((int)(Resultlist[i].MutationList.size() + QueryMutation.size() - Resultlist[i].SharedMutationNum * 2));
                    int mutdj = abs((int)(Resultlist[j].MutationList.size() + QueryMutation.size() - Resultlist[j].SharedMutationNum * 2));
                    if (mutdi > mutdj)
                    {
                        tmp = Resultlist[i]; Resultlist[i] = Resultlist[j]; Resultlist[j] = tmp;
                    }
                }
            }
        }

    printf("Sorting Finished, Time: %ld\n", time(NULL));
}

int main(int argc, char *argv[])//4个输入，目标文件夹，默认路径Valkyrie/Linux; 待查询序列; 输出文件路径；是否使用连锁不平衡（默认不用）
{
    //string fold = "/xtdisk/limk_group/mawt/Valkyrie/LinuxTest";
    string fold = "/home/mawentai/APP/Valkyrie/Linux";
    string outputfold;//输出文件,不是路径
    ifstream readtitle("/home/mawentai/APP/Valkyrie/README.txt");
    //ifstream readtitle("/xtdisk/limk_group/mawt/Valkyrie/README.txt");
    printf("\n\n\n\n\n");
    string title;
    while (getline(readtitle, title))
        printf("%s\n", title.c_str());
    readtitle.close();
    printf("\nValkyrie Core Version: 0.4.5\n");
    starttime = time(NULL);//记录程序开始时间
    printf("\nStart Time: %ld\n", starttime);
    time_t now = time(0);
    cout << ctime(&now) << endl;
    printf("\n\nWelcome to Vallhalla!\n\nA Valkyrie Has Been Assigned To Help You Find The Best Match\n\n");
    printf("\n%s Hello, Nice to meet you!\n\n", "Message from your Valkyrie:");
    if (argc == 1)
    {
        printf("ERROR /// GDC: Query Sequence Needed ///\n");
        return -1;
    }
    
    string QuerySeq;
    if (argc > 3) { fold = argv[1]; QuerySeq = argv[2]; outputfold = argv[3]; }
    else if (argc == 3) { QuerySeq = argv[1]; outputfold = argv[2]; }
    else { printf("Please Check Your Input: Mutation And Outputdir!\n"); return -1; }
    if (argc == 5)
    {
        string method = argv[4];
        if (method == "-LD")
        {
            UseLinkageDisequilibrium = true;//根据输入判断是否使用连锁不平衡修正
            printf("-LD: Use Linkage Disequilibrium\n");
        }
        else
            if (method == "-MO")
            {
                UseMutationOnly = true;//根据输入判断是否使用不附加权重
                printf("-MO: Use Mutation Number Distance Only\n");
            }
            else
                if (method == "-CL")
                {
                    UseLineageRate = true;//根据输入判断是否使用进化枝特异的突变率
                    printf("-CL: Use Lineage Mutation Rate\n");
                }
                else
                    if(method == "-SM")
                {
                    UseSharedMutation = true;//根据输入判断是否使用共有突变
                    printf("-SM: Use Shared Mutation Number\n");
                }
                else
                    if (method == "-LDCL")
                    {
                        UseLineageRate = true;//根据输入判断是否使用进化枝特异的突变率
                        UseLinkageDisequilibrium = true;//根据输入判断是否使用连锁不平衡修正
                        printf("-LDCL: Use Linkage Disequilibrium With Lineage Mutation Rate\n");
                    }
                    else
                    {
                        printf("ERROR /// GDC: Unknown Methods ERROR ///\n");
                        printf("-LD: Use Linkage Disequilibrium\n");
                        printf("-MO: Use Mutation Number Distance Only\n");
                        printf("-CL: Use Lineage Mutation Rate\n");
                        printf("-SM: Use Shared Mutation Number\n");
                        printf("-LDCL: Use Linkage Disequilibrium With Lineage Mutation Rate\n");
                        return 1;
                    }

    }
    //string QuerySeq = "Dalian722:241C/T,3037C/T,14408C/T,23403A/G,28881G/A,28882G/A,28883G/C,2091C/T,5128A/G,8360A/G,13860C/T,19839T/C,19999G/T,28905C/T";//Debug

    vector<string> QuerySeq1 = splitStr(QuerySeq, ':');
    if (QuerySeq1.size() < 2)return 1;//wrong input;
    printf("Valkyrie Recevied your Query Seq: %s\n", QuerySeq1[0].c_str());
    if (QuerySeq1[1].length() == 0)
    {
        printf("Your Query Seq Do Not Contain Any Mutations\n");
        printf("\n%s Please stop teasing me!\n", "Message from your Valkyrie:");
        printf("Query End\n");
    }
    
    int i,j;
    vector<string> QuerySeq2 = splitStr(QuerySeq1[1], ',');
    for (i = 0; i < QuerySeq2.size(); i++)
        if (QuerySeq2[i] != "" && QuerySeq2[i] != " ")QueryMutation.push_back(QuerySeq2[i]);//过滤输入文件

    printf("Your Query Seq Contains %d Mutation\n", QueryMutation.size());

    printf("\n\nStart Loading Relevant Files\n");//开始读入依赖文件
        if (Datainput(fold) != 0)
        {
            printf("Loading Relevant Files Failed\n");
            printf("\n%s Please check Relevant Files.\n", "Message from your Valkyrie:");
            return -1;
        }
    
    printf("Loading Relevant Files Success, Time: %ld\n\n",time(NULL));

    /*vector<string> Muttest; Muttest.emplace_back("28881G/A"); Muttest.emplace_back("28882G/A"); Muttest.emplace_back("28883G/C");
    vector<double> Weighttest; Weighttest.emplace_back(100); Weighttest.emplace_back(100); Weighttest.emplace_back(100);
    double LDD = LDGramSchmidtProcess(Muttest, Weighttest);*/
    GenotypeSearch();//遍历genotype计算距离
    
    printf("\n%s Search ended，results in printing, please wait.\n", "Message from your Valkyrie:");
    
    GenoType tmp;
    ofstream write(outputfold);

    ofstream writeseqmut((string)(outputfold + ".mutlist"));//输出突变准备画树
    string outputseqmut = QuerySeq1[0] + "\t";
    for (i = 0; i < QueryMutation.size(); i++)
        outputseqmut += QueryMutation[i] + " ";
    outputseqmut += "\n";

    if (!write.is_open())
    {
        printf("ERROR /// Output: Can Not Open %s ///\n", outputfold.c_str()); return -1;
    }
    write << "QuerySeqName\tSeqID\tAccession\tScore\tLineage\tMutation\tSubject private\tQuery private" << endl;
    double highestrank = Resultlist[49].HighestRank;
    for (i = 0; i < resultTopN || (Resultlist[i].HighestRank >= highestrank && highestrank!=0); i++)
    {
        if (i >= Resultlist.size())break;
        map<string, int> muthash;
        int MutDis = 0;
        string MutGenotypeOnly = "";
        string MutQueryOnly = "";
        //---------------

        for (j = 0; j < Resultlist[i].MutationList.size(); j++)//Genotype特有的突变
        {
            muthash.insert(pair<string, int>(Resultlist[i].MutationList[j], j));
            if (QueryMutationMap.find(Resultlist[i].MutationList[j]) == QueryMutationMap.end())
            {
                MutGenotypeOnly += Resultlist[i].MutationList[j] + " ";
                MutDis++;
            }
        }
        //---------------

        for (j = 0; j < QueryMutation.size(); j++)//Query特有的突变
            if (muthash.find(QueryMutation[j]) == muthash.end())
            {
                MutQueryOnly += QueryMutation[j] + " ";
                MutDis++;
            }
        if (MutGenotypeOnly != "")MutGenotypeOnly = MutGenotypeOnly.substr(0, MutGenotypeOnly.length() - 1);
        if (MutQueryOnly != "")MutQueryOnly = MutQueryOnly.substr(0, MutQueryOnly.length() - 1);
        //---------------

        string lineageout = "";
        for (int valj = 0; valj < Resultlist[i].LineageBelonged.size(); valj++)//所属Lineage
            if (Resultlist[i].LineageBelonged[valj] != "MetadataUnknown")
            {
                lineageout += Resultlist[i].LineageBelonged[valj] + " ";
                break;
            }
        if (lineageout == "")lineageout = "Unknown";

        

        write << QuerySeq1[0] << "\tV" << Resultlist[i].GenoTypeID << "\t" << Resultlist[i].representiveSeq << "\t" << Resultlist[i].HighestRank << "\t" << lineageout << "\t" << MutDis << "\t" << MutGenotypeOnly << "\t" << MutQueryOnly << endl;
    
        outputseqmut += "V" + to_string(Resultlist[i].GenoTypeID) + "\t";
        for (j = 0; j < Resultlist[i].MutationList.size(); j++)
            outputseqmut += Resultlist[i].MutationList[j] + " ";
        outputseqmut += "\n";
    }
    write.close();
    writeseqmut << outputseqmut << endl;
    writeseqmut.close();

    printf("Printing Ended\n");
    printf("\nBest Match Find For Yor Query: Genotype %d\n", Resultlist[0].GenoTypeID);
    printf("\nQuery Ended\nEnd Time: %ld\nTotal Time Used: %ld Seconds\n",time(NULL),(time(NULL)-starttime));
    now = time(0);
    cout << ctime(&now) << endl;

    vector<string> testMutlist;
    vector<double> testWeightList;
    return 0;
}