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
#include <algorithm>
using namespace std;

class TreeNode
{
public:
    string Id;//树节点的名称 Name of the node
    int Num;//节点的序号 Number of the node
    vector<string> Mutlist;//节点序列上的突变 Mutation of the node
    int Father;//节点父节点的序号 Number of its father node
    vector<int> ChildList;//节点的孩子的序号 Number of its childern node
    double branchLength = 0;//枝长 branch length
};

static vector<TreeNode> NodeList;//存储所有节点 save all nodes
static vector<string> BackWardMut;//回复突变细节list detail list of back mutations
static vector<string> BackMutList;//回复突变list list of back mutations
static vector<int> HotMutPassCount;//突变在树上发生的传递次数 happened times of mutation
static vector<int> BackMutPassCount;//回复突变在树上发生的传递次数 inheritage times of back mutation
static vector<int> BackMutCount;//回复突变的数量list number of back mutation
static vector<string> HotMutList;//热点突变list list of recurrent mutation
static vector<string> MutEventList;//突变事件list list of mutation event
static vector<int> HotMutCount;//热点突变的数量list number of recurrent mutation
static vector<string> leafname;//叶子节点的名字 name of leaf node
static vector<string> leafmut;//叶子节点的突变 mutation of leaf node
static vector<double> brachLengthList;//所有的枝长 all branch length
static double branchLengthThreshold;//99%枝长范围，用于过滤 a 99% threshold to filter branch length

static string workfold;//工作目录 working dir
static string Lineagename;//搜索对象的id Tree_LineageCombined_【1】 1 the id of tree for searching

static unordered_map<string, string> inDelRecodeList;//<重新编码的indel，原始indel> <recoded indel, original indel>
static ofstream writemut;//CSSMut.tsv
static ofstream writehot_mutlist;//Hot.mutlist


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

//读入重建的祖先序列数据
//Read in reconstructed ancestor sequences
int Readin()
{
    int i, j, k;
    string indelgenome = "";
    //读入indel.info 还原重新编码了的indel
	//Read in indel.info, recover recoded indel
    ifstream readindel(workfold + "inDel.info");
    if (!readindel.is_open()) { printf("ERROR /// FHB: Can Not Open %s ///\n", (workfold + "inDel.info").c_str()); return -1; }
    string line;
    getline(readindel, line);
    while (getline(readindel, line))
    {
        vector<string> line1 = splitStr(line, '\t');
        inDelRecodeList.insert(pair<string, string>(line1[1], line1[0]));
        indelgenome += line1[1][line1[1].length() - 3];
    }
    readindel.close();

    //先读入叶子及其突变 【用来制造序列的突变文件】
	//Read in leaf node and its mutations
    ifstream readleaf(workfold + "Tree_LineageCombined_" + Lineagename + ".tsv");
    if (!readleaf.is_open()) { printf("ERROR /// FHB: Can Not Open %s", (workfold + "Tree_LineageCombined_" + Lineagename + ".tsv").c_str()); return -1; }
    while (getline(readleaf, line))
    {
        vector<string> linex = splitStr(line, '\t');
        vector<string> liney; 
        if(linex.size()>1)
        liney = splitStr(linex[1], ' ');
        string mutrecode = "";
        for (i = 0; i < liney.size(); i++)//查找是否包含了重编码的indel，如果有则替换回去 find if contains any recoded indel, if so replace back
            if (inDelRecodeList.find(liney[i]) != inDelRecodeList.end())
                mutrecode += inDelRecodeList[liney[i]] + ' ';
            else
                mutrecode += liney[i] + ' ';
        leafname.push_back(linex[0]); leafmut.push_back(mutrecode.substr(0,mutrecode.length() - 1));
    }
    readleaf.close();

    //读入IQtree 中间节点 【iqtree推测的祖先序列】
	//Read in mid node on the tree
    ifstream read(workfold + "Tree_LineageCombined_" + Lineagename + ".fa.state");
    if(!read.is_open()) { printf("ERROR /// FHB: Can Not Open %s /// \n", ("Tree_LineageCombined_" + Lineagename + ".fa.state").c_str()); return -1; }
    ifstream readref(workfold + "reference.fa");
    if (!readref.is_open()) { printf("ERROR /// FHB: Can Not Open %s ///\n", (workfold + "reference.fa").c_str()); return -1; } 
    
    getline(readref,line); getline(readref, line);
    string reference = "+" + line; readref.close();//读入参考基因组 read in reference genome
    for (i = 29904; i <= 30000; i++)
        reference += '-';//基因组延长 extend the reference genome
    reference += indelgenome;
    
    getline(read,line);
    while (line.find("#")!=-1)
        getline(read, line);
    getline(read, line);
    string nodename = "start";
    TreeNode b;
    NodeList.push_back(b);
    while (getline(read, line))//中间节点突变 mid node mutation
    {
        TreeNode c;
        vector<string> a1 = splitStr(line, '\t');
        nodename = a1[0]; c.Id = nodename;
        while (a1[0] == nodename)//节点 node
        {
            if (stoi(a1[1]) <= 29903)
                if (a1[2][0] != reference[stoi(a1[1])])
                    c.Mutlist.push_back(a1[1] + reference[stoi(a1[1])] + "/" + a1[2]);
            if (stoi(a1[1]) >= 30001)//重编码位置 the position of recoding
            {
                if (a1[2][0] != reference[stoi(a1[1])])
                    c.Mutlist.push_back(inDelRecodeList[a1[1] + reference[stoi(a1[1])] + "/" + a1[2]]);
            }
            if (!getline(read, line)) break;
            a1 = splitStr(line, '\t');
        }
        NodeList.push_back(c);
    }

    read.close();
    return 1;
}

//连接节点关系，创建树结构
//connect all nodes, restore the tree topological structure
void BuildConnect(string nwk, int father)
{
    int i, j, k;
    if (nwk.find("Node")==-1)//切到叶子节点了 cut till leaf node
    {
        TreeNode d;
        d.Father = father;
        vector<string> a3 =splitStr(nwk,':');
        if (a3.size() > 1 && a3[1][a3[1].length() - 1] <= '9' && a3[1][a3[1].length() - 1] >= '0')
        {
            d.branchLength = stof(a3[1]);
            brachLengthList.push_back(d.branchLength);
        }
        int n = IndexOf(leafname,a3[0]);
        d.Id = a3[0];
        if (n != -1)
        {
            vector<string> a =splitStr( leafmut[n],' ');
            for (i = 0; i < a.size(); i++)
                if (a[i] != "") d.Mutlist.push_back(a[i]);
        }

        NodeList[father].ChildList.push_back(NodeList.size());
        NodeList.push_back(d);
        return;
    }
    else
    {
        TreeNode e;
        k = nwk.length() - 1;
        while (nwk[k] != ')') k--;
        string thisnode = nwk.substr(k + 1, nwk.length() - k - 1);
        vector<string> a2 = splitStr(thisnode, ':');
        if (a2.size() > 1 && a2[1][a2[1].length() - 1] <= '9' && a2[1][a2[1].length() - 1] >= '0')
        {
            e.branchLength = stof(a2[1]);
            brachLengthList.push_back(e.branchLength);
        }
        if(a2[0].find('/')==-1)//iqtree的bootstrap用“/”号加在node后面 如： )Node85/32:0.0000010000)
            e.Id = a2[0];
        else
        {
            vector<string> a3 = splitStr(a2[0], '/');
            e.Id = a3[0];
        }
        for (j = 0; j < NodeList.size(); j++)
            if (NodeList[j].Id == e.Id) { break; }
        if (j == NodeList.size())
            printf("ERROR /// FHB: Can Not Find Node_%s In The Node List, Please Check .state File ///\n", e.Id.c_str());
        NodeList[father].ChildList.push_back(j);
        //NodeList.Add(e);//node在初始化阶段已经加进去了
        NodeList[j].Father = father;
        vector<int> douhao;//逗号 ','
        int nkuohao = 0;//处于括号的层数 number of layers in parentheses
        for (i = 1; i < k; i++)
        {
            if (nwk[i] == '(') nkuohao++; if (nwk[i] == ')') nkuohao--;
            if (nwk[i] == ',' && nkuohao == 0) douhao.push_back(i);//找到所有亲儿子括号的位置find the position of brackets for all sons
        }
        if (douhao.size() == 1)//该节点两个孩子 This node has two children
        {
            BuildConnect(nwk.substr(1, douhao[0] - 1), j);
            BuildConnect(nwk.substr(douhao[douhao.size() - 1] + 1, k - douhao[douhao.size() - 1] - 1), j);
        }
        else
            if (douhao.size() > 1)//该节点好多孩子？？？More than two ???
            {
                BuildConnect(nwk.substr(1, douhao[0] - 1), j);
                for (i = 0; i < douhao.size() - 1; i++)
                    BuildConnect(nwk.substr(douhao[i] + 1, douhao[i + 1] - douhao[i] - 1), j);
                BuildConnect(nwk.substr(douhao[douhao.size() - 1] + 1, k - douhao[douhao.size() - 1] - 1), j);
            }
    }
    return;
}

//递归查找回复突变 Recursive search for backward mutation
void BackMutSearch(int NodeNow, int father)
{
    if (NodeList[NodeNow].ChildList.size() == 0) return;
    int i, j, k, m;
    string mutout = NodeList[NodeNow].Id + "\t";
    for (i = 0; i < NodeList[NodeNow].Mutlist.size(); i++)
        mutout += NodeList[NodeNow].Mutlist[i] + " ";
    //writemut << mutout << endl;
    for (i = 0; i < NodeList[NodeNow].ChildList.size(); i++)
    {
        if (NodeList[NodeList[NodeNow].ChildList[i]].branchLength > branchLengthThreshold)
            continue;//99%枝长过滤 99% branch length filtration

        int p = NodeList[NodeNow].ChildList[i];
        k = 0;
        for (j = 0; j < NodeList[NodeNow].Mutlist.size(); j++)
        {
            if (IndexOf( NodeList[p].Mutlist,NodeList[NodeNow].Mutlist[j])==-1)
                k++;
        }
        //string outputa = "" + NodeList[NodeNow].Id + "\t";
        string outputa = "";
        for (j = 0; j < NodeList[p].Mutlist.size(); j++)
            if (IndexOf( NodeList[NodeNow].Mutlist,(NodeList[p].Mutlist[j]))==-1)
            {
                if (IndexOf(HotMutList,NodeList[p].Mutlist[j])!=-1)
                    HotMutCount[IndexOf(HotMutList, NodeList[p].Mutlist[j])]++;
                else
                    {
                        HotMutList.push_back(NodeList[p].Mutlist[j]);
                        HotMutCount.push_back(1);
                    }

                outputa += NodeList[p].Mutlist[j] + " ";

            }
        if (outputa!="")
            MutEventList.push_back("Node"+ NodeList[NodeNow].Id +"->Node"+NodeList[p].Id+"\t"+outputa.substr(0, outputa.length() - 1));
        if (k == 1)//在传递过程中最多允许1个回复突变/1代 A maximum of 1 reverse mutation per generation is allowed during transmission
        {
            for (j = 0; j < NodeList[NodeNow].Mutlist.size(); j++)
                if (IndexOf( NodeList[p].Mutlist,NodeList[NodeNow].Mutlist[j])==-1)
                {
                    string output = NodeList[NodeNow].Id + "\t";
                    for (k = 0; k < NodeList[NodeNow].Mutlist.size(); k++)
                        output += NodeList[NodeNow].Mutlist[k] + " ";
                    output += "\t->\t" + NodeList[p].Id + "\t";
                    for (k = 0; k < NodeList[p].Mutlist.size(); k++)
                        output += NodeList[p].Mutlist[k] + " ";
                    output += "\t" + NodeList[NodeNow].Mutlist[j];
                    BackWardMut.push_back(output);
                    if (IndexOf( BackMutList,NodeList[NodeNow].Mutlist[j])!=-1)
                        BackMutCount[IndexOf(BackMutList, NodeList[NodeNow].Mutlist[j])]++;
                    else
                        {
                            BackMutList.push_back(NodeList[NodeNow].Mutlist[j]);
                            BackMutCount.push_back(1);
                        }
                }
        }
    }
    for (i = 0; i < NodeList[NodeNow].ChildList.size(); i++)
        BackMutSearch(NodeList[NodeNow].ChildList[i], NodeNow);
    return;
}

//递归查找突变的传递次数
//Pass times of recursive search mutation
void Passcount()
{
    int i, j, k, m;
    for (j = 1; j < NodeList.size(); j++)
    {
        for (k = 0; k < BackMutList.size(); k++)
        {   
            //加上子节点的数量表示传递次数，叶子节点没有孩子 Adding the number of child nodes indicates the number of passes, and the leaf node has no children
            if (IndexOf(NodeList[j].Mutlist, BackMutList[k]) != -1)
                BackMutPassCount[k] += NodeList[j].ChildList.size();
        }
        for (k = 0; k < HotMutList.size(); k++)
        {
            if (IndexOf(NodeList[j].Mutlist, HotMutList[k]) != -1)
                HotMutPassCount[k] += NodeList[j].ChildList.size();
        }
    }
    return;
}

int main(int argc, char* argv[])//两个输入 工作目录 树的编号 Two input: working directory and index of the tree
{
    if (argc < 3) { printf("ERROR /// FHB: Input Need Workfold And Treename ///\n"); return -1; }
    workfold = argv[1];
    Lineagename = argv[2];
    //workfold = "/home/mawentai/Valkyrie/LinuxTest"; Lineagename = "1";//debug

    printf("Find Hotspot Backmutation Start Tree_%s\n",Lineagename.c_str());
    workfold += "/";
    /*ofstream writemut(workfold + "LC_" + Lineagename + "CSSmut.tsv");
    ofstream writehot_mutlist(workfold + "LC_" + Lineagename + "_Hot.mutlist");
    if (!writemut.is_open())
    {
        printf("ERROR /// FHB: Can Not Open %s ///\n", (workfold + "LC_" + Lineagename + "CSSmut.tsv").c_str()); return -1;
    }
    if (!writehot_mutlist.is_open())
    {
        printf("ERROR /// FHB: Can Not Open %s ///\n", (workfold + "LC_" + Lineagename + "_Hot.mutlist").c_str()); return -1;
    }*/
    //读入iqtree产生的树文件 Read in the tree file generated by iqtree
    ifstream read(workfold + "Tree_LineageCombined_" + Lineagename + ".fa.treefile");
    if (!read.is_open())
    {
        printf("ERROR /// FHB: Can Not Open Treefile ///\n"); return -1;
    }
    //输出 output
    //ofstream write(workfold + "Tree_LineageCombined_" + Lineagename + "_回复突变.txt");
    ofstream writebacklist(workfold + "Tree_LineageCombined_" + Lineagename + ".backmutation");
    ofstream writehot(workfold + "Tree_LineageCombined_" + Lineagename + ".hotspot");
    ofstream writeMutEvent(workfold + "Tree_LineageCombined_" + Lineagename + ".mutevent");
    ofstream writeNodeMut(workfold + "Tree_LineageCombined_" + Lineagename + ".nodemutation");
    writebacklist<< "Backmutation\tHappenedTimes\tInheritanceTimes" <<endl;//回复突变\t发生次数\t突变在序列中出现的次数
    writehot<< "Hotspot\tHappenedTimes\tNumofSeqContains\tInheritanceTimes\tTreeNodeNumber" <<endl;
    //write<<"Father\tMut\t\tChild\tMut\t回复突变"<<endl;
    string nwk;
    getline(read, nwk);
    if (Readin() == -1)//读入重建的祖先序列数据 Read in the reconstructed ancestor sequence data
    {
        printf("Data Input Failed\n");
        return -1;
    }
    int root = NodeList.size();
    BuildConnect(nwk.substr(0, nwk.length() - 1), 0);//重建树的结构 Rebuild the structure of the tree

    int i, j, k;
    double tmpd;
    //排序枝长列表，确定99%阈值 Sort the branch length list to determine the 99% threshold
    sort(brachLengthList.begin(), brachLengthList.end());
    branchLengthThreshold = brachLengthList[int( double(brachLengthList.size()) * 0.99)];
    cout << "Branch Length Threshold: " << branchLengthThreshold << endl;
    for (k = 0; k < NodeList.size(); k++)
        if (NodeList[k].Id == "Node1") break;
    if (k == NodeList.size())
    {
        printf("Where is your Node1 ?\n");
        return -1;
    }
    BackMutSearch(k, -1);//递归查找回复突变 Recursive search for back mutations
    //for (i = 0; i < BackWardMut.size(); i++)
    //    write<<BackWardMut[i]<<endl;
    //write.close();
    for (i = 0; i < MutEventList.size(); i++)//输出突变事件 Output mutation event
        writeMutEvent << MutEventList[i] << endl;
    writeMutEvent.close();

    string tmpstring; int tmpint;
    for (i = 0; i < BackMutList.size(); i++)
        for (j = i + 1; j < BackMutList.size(); j++)
            if (BackMutCount[i] < BackMutCount[j])
            {
                tmpstring = BackMutList[i]; BackMutList[i] = BackMutList[j]; BackMutList[j] = tmpstring;
                tmpint = BackMutCount[i]; BackMutCount[i] = BackMutCount[j]; BackMutCount[j] = tmpint;
            }
    for (i = 0; i < HotMutList.size(); i++)
        for (j = i + 1; j < HotMutList.size(); j++)
            if (HotMutCount[i] < HotMutCount[j])
            {
                tmpstring = HotMutList[i]; HotMutList[i] = HotMutList[j]; HotMutList[j] = tmpstring;
                tmpint = HotMutCount[i]; HotMutCount[i] = HotMutCount[j]; HotMutCount[j] = tmpint;
            }
    for (i = 0; i < BackMutList.size(); i++)
        BackMutPassCount.push_back(0);
    for (i = 0; i < HotMutList.size(); i++)
        HotMutPassCount.push_back(0);
    Passcount();
    for (i = 0; i < BackMutList.size(); i++)//输出回复突变 output back mutations
    {
        writebacklist<<BackMutList[i] + "\t" <<BackMutCount[i] << "\t" << BackMutPassCount[i]<<endl;
    }
    for (i = 0; i < HotMutCount.size(); i++)//输出突变热点【改成输出所有突变】output all mutations
    {
        if (HotMutCount[i] > 0)
        {
            int n = 0;
            for (j = 0; j < leafmut.size(); j++)
                if (leafmut[j].find(HotMutList[i])!=-1) n++;
            writehot<<HotMutList[i] << "\t" << HotMutCount[i] << "\t" << n<< "\t" << HotMutPassCount[i] << "\t" << NodeList.size() << endl;
        }
    }
    for (i = 0; i < NodeList.size(); i++)//输出节点所具有的突变 output mutations of each node
    {
        string output = NodeList[i].Id + "\t";
        for (j=0;j<NodeList[i].Mutlist.size();j++)
        {
            output += NodeList[i].Mutlist[j] + " ";
        }
        writeNodeMut << output << endl;
    }

    writehot.close();
    writebacklist.close();
    writemut.close();
    writehot_mutlist.close();
    writeNodeMut.close();
    return 0;
}