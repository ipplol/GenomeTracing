#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

vector<string> splitStr(string str, char delimiter) //string split函数
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

int main(int argc, char *argv[])//需要四个输入  1.reference 2.indel 3.mutation 4.output path
{
    if (argc < 5)
    {
        cout << "ERROR: /// Mut2Seq v8 require 4 input file by order: 1.reference 2.indel 3.mutation 4.output path ///" << endl;
        return 0;
    }
     
    ifstream readref(argv[1]);//读入参考基因组
    ifstream readindel(argv[2]);//读入indel.info
    ifstream readmut(argv[3]);//读入mutation
    ofstream writefa(argv[4]);//输出序列
    //ifstream readref("/xtdisk/limk_group/mawt/Valkyrie/Linux/reference.fa");//读入参考基因组
    //ifstream readindel("/xtdisk/limk_group/mawt/Valkyrie/Linux/inDel.info");//读入indel.info
    //ifstream readmut("/xtdisk/limk_group/mawt/Valkyrie/Linux/Tree_LineageCombined_1.tsv");//读入mutation
    //ofstream writefa("/xtdisk/limk_group/mawt/Valkyrie/Linux/Tree_LineageCombined_1.fa");//输出序列

    if (!readref.is_open()) //没找到文件1
    {
        cout << "ERROR: /// Mutation2Sequence: Can Not Open Reference Genome. ///" << endl; return -1;
    }
    if (!readindel.is_open()) //没找到文件2
    {
        cout << "ERROR: /// Mutation2Sequence: Can Not Open indel.info. ///" << endl; return -1;
    }
    if (!readmut.is_open()) //没找到文件3
    {
        cout << "ERROR: /// Mutation2Sequence: Can Not Open Mutation File. ///" << endl; return -1;
    }
    if (!writefa.is_open()) //没找到文件4
    {
        cout << "ERROR: /// Mutation2Sequence: Can Not Create Output File. ///" << endl; return -1;
    }

    //读入参考基因组
    string reference = ">";
    string line;
    getline(readref, line);
    getline(readref, line);
    reference += line;//29903 original

    //读入indel，延长基因组
    int i, j, k=0;
    string genomeDLC = "";
    getline(readindel, line);
    while (getline(readindel, line))
    {
        k++;
        genomeDLC += line[line.length() - 3];
    }

    //读入mutation，转换，输出
    if (k != 0) 
    { 
        for (i = 29904; i <= 30000; i++)reference += "N";
        reference += genomeDLC; 
    }
    writefa << ">Wuhan-Hu-1\n" << reference.substr(1, reference.length() - 1) << endl;
    while (getline(readmut,line))
    {
        if (line != "")
        {
            vector<string> line1;
            vector<string> mut1;
            if (line.find(':') == line.npos && line.find('\t') != line.npos)
            {
                line1 = splitStr(line, '\t');
                if(line1.size()>1)
                mut1 = splitStr(line1[1], ' ');
            }
            else
                if (line.find('\t') == line.npos && line.find(':') != line.npos)
                {
                    line1 = splitStr(line, ':');
                    if (line1.size() > 1)
                    mut1 = splitStr(line1[1], ',');
                }
                else
                    continue;
            string outfa = reference;
            for (i = 0; i < mut1.size(); i++)
            {
                string pos = mut1[i].substr(0, mut1[i].length() - 3);
                if(pos[pos.size() - 1]<='9' && pos[pos.size() - 1] >= '0')
                    if(stoi(pos)<outfa.size())
                        outfa[stoi(pos)] = mut1[i][mut1[i].length() - 1];
            }
            writefa << ">" + line1[0] + "\n" << outfa.substr(1, outfa.length() - 1) << endl;
        }
    }

    readindel.close();
    readmut.close();
    readref.close();
    writefa.close();
    return 0;
}