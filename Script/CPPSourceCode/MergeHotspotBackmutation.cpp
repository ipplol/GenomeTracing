#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
//#include <thread>
//#include <pthread.h>
#include <map>
#include <unordered_map>
#include <time.h>
using namespace std;

class Mut
{
public:
    string Mutation = "";
    int HappenTimes = 0;
    int ShowPassTimes = 0;
    int InheritanceTimes = 0;
    int TreeNodeNumber = 0;
};
//vector indexof
int IndexOf(vector<string> v, string sub)
{
    int i;
    for (i = 0; i < v.size(); i++)
        if (v[i] == sub)return i;
    return -1;
}
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

int main(int argc, char* argv[])//一个输入，目标文件所在的目录
{
    printf("Merge Hotspot And Backmutation From Tree_LineageCombined\n");
    
    vector<string> HotspotFiles;
    vector<string> BackmutationFiles;
    if (argc < 2) {
        printf("ERROR /// MHB: Target Fold Needed /// \n");
        return 0;
    }
    string treefile = argv[1];
    string fold = argv[1];
    //string treefile = "/home/mawentai/Valkyrie/LinuxTest";//debug
    //string fold = "/home/mawentai/Valkyrie/LinuxTest"; //debug
    treefile += "/TreeFileLineageDetail.tsv";
    ifstream readtree(treefile);
    int i, j, k;
    string line;
    while (getline(readtree, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        vector<string> line2 = splitStr(line1[0], '_');
        HotspotFiles.push_back(fold + "/Tree_LineageCombined_" + line2[1] + ".hotspot"); //Tree_LineageCombined_9.hotspot
        BackmutationFiles.push_back(fold + "/Tree_LineageCombined_" + line2[1] + ".backmutation");//Tree_LineageCombined_9.backmutation
    }
    readtree.close();
    vector<Mut> hotspotlist;
    vector<string> hotspots;
    vector<Mut> backmutlist;
    vector<string> backmuts;
    for (i = 0; i < HotspotFiles.size(); i++)
    {
        ifstream readh(HotspotFiles[i]);//hotspot_i
        ifstream readb(BackmutationFiles[i]);//backmutation_i
        getline(readh, line); getline(readb, line);

        if (readh.is_open())
        {
            while (getline(readh, line))//first read hotspot_i
            {
                vector<string> line1 = splitStr(line, '\t');
                if (IndexOf(hotspots, line1[0]) == -1)//do not contain
                {
                    hotspots.push_back(line1[0]);
                    Mut a;
                    a.Mutation = line1[0];
                    a.HappenTimes = stoi(line1[1]);
                    a.ShowPassTimes = stoi(line1[2]);
                    a.InheritanceTimes = stoi(line1[3]);
                    a.TreeNodeNumber = stoi(line1[4]);
                    hotspotlist.push_back(a);
                }
                else
                {
                    k = IndexOf(hotspots, line1[0]);
                    hotspotlist[k].HappenTimes += stoi(line1[1]);
                    hotspotlist[k].ShowPassTimes += stoi(line1[2]);
                    hotspotlist[k].InheritanceTimes += stoi(line1[3]);
                    hotspotlist[k].TreeNodeNumber += stoi(line1[4]);
                }
            }
        }
        else
            cout << "Can not find file: " << HotspotFiles[i] << endl;

        if (readb.is_open())
        {
            while (getline(readb, line))//second read backmutation_i
            {
                vector<string> line1 = splitStr(line, '\t');
                if (IndexOf(backmuts, line1[0]) == -1)//do not contain
                {
                    backmuts.push_back(line1[0]);
                    Mut a;
                    a.Mutation = line1[0];
                    a.HappenTimes = stoi(line1[1]);
                    a.ShowPassTimes = stoi(line1[2]);
                    backmutlist.push_back(a);
                }
                else
                {
                    k = IndexOf(backmuts, line1[0]);
                    backmutlist[k].HappenTimes += stoi(line1[1]);
                    backmutlist[k].ShowPassTimes += stoi(line1[2]);
                }
            }
        }
        else
            cout << "Can not find file: " << BackmutationFiles[i] << endl;

        readb.close();
        readh.close();
    }
    Mut tmp;
    for (i = 0; i < hotspotlist.size(); i++)//sort by happentimes
        for(j=i+1;j<hotspotlist.size();j++)
            if (hotspotlist[i].HappenTimes < hotspotlist[j].HappenTimes)
            {
                tmp = hotspotlist[i]; hotspotlist[i] = hotspotlist[j]; hotspotlist[j] = tmp;
            }
    for (i = 0; i < backmutlist.size(); i++)//sort by happentimes
        for (j = i + 1; j < backmutlist.size(); j++)
            if (backmutlist[i].HappenTimes < backmutlist[j].HappenTimes)
            {
                tmp = backmutlist[i]; backmutlist[i] = backmutlist[j]; backmutlist[j] = tmp;
            }

    ofstream writeh(fold + "/Tree_LineageCombined_Total.hotspot");
    ofstream writeb(fold + "/Tree_LineageCombined_Total.backmutation");
    writeh << "Hotspot\tHappenedTimes\tNumofSeqContains\tInheritanceTimes\tTreeNodeNumber" << endl;
    writeb << "Backmutation\tHappenedTimes\tInheritanceTimes" << endl;
    for (i = 0; i < hotspotlist.size(); i++)
        writeh << hotspotlist[i].Mutation << "\t"<< hotspotlist[i].HappenTimes << "\t" << hotspotlist[i].ShowPassTimes << "\t" << hotspotlist[i].InheritanceTimes << "\t" << hotspotlist[i].TreeNodeNumber << endl;
    for (j = 0; j < backmutlist.size(); j++)
        writeb << backmutlist[j].Mutation << "\t" << backmutlist[j].HappenTimes << "\t" << backmutlist[j].ShowPassTimes << endl;
    writeb.close();
    writeh.close();
    printf("MHB Done\n");
    return 1;
}