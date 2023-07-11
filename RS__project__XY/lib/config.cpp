//
// Created by 11379 on 2020-10-30.
//

#include <fstream>
#include <sstream>
#include "config.h"
#include "pugixml.hpp"
#include "ChrOper.h"
#include "tools.h"

void DeleteFirstLineInFile(string fileName) {
    vector<string> VecContent;
    string StrLine;
    ifstream iFile(fileName);
    if (!iFile) {
        cout << fileName << "filelist open failed!\n";
        exit(0);
    }
    //Read all content in the document into "VecContent"
    while (iFile) {
        getline(iFile, StrLine);
        VecContent.push_back(StrLine);
    }
    iFile.close();
    VecContent.erase(VecContent.begin()); // delete the first line
    ofstream oFile(fileName);
    if (!oFile) {
        cout << fileName << "filelist open failed!\n";
        exit(0);
    }
    vector<string>::const_iterator iter = VecContent.begin();
    //{Rewrite the contents of "vecContent" into the file}.
    for (; VecContent.end() != iter; ++iter) {
        oFile.write((*iter).c_str(), (*iter).size());
        oFile << '\n';
    }
    oFile.close();
}

void ClearALL() {
    Tasks.clear();
    TaskListInLevel.clear();
    LevelIdOfTask.clear();
    VMTypes.clear();
    ParChildTranFileSizeSum.clear();
    AdjMatrix.clear();
    Descendants.clear();
    Ancestors.clear();
//    population.clear();
}

void MDMOGA_ConfigInputParameter() {
    inputParameter.NumOfChromPerPop = ceil(30+0.3*comConst.NumOfTasks);
    inputParameter.MutationProb = 0.2;
    inputParameter.NumOfIterations = 50;
//    inputParameter.theta_eff = 0.5;         //theta_eff为0，时间越短越好；theta_eff为1，成本越短越好--初始化
}

void MOELS_ConfigInputParameter() {
    inputParameter.NumOfChromPerPop = 50;
    inputParameter.NumOfIterations = 50;
    inputParameter.MutationProb = 0.7;
    inputParameter.CrossProb = 0.7;
    inputParameter.Sbx = 10;
    inputParameter.PolyMu = 20;
}

void EMS_ConfigInputParameter() {
    inputParameter.NumOfChromPerPop = 50;
    inputParameter.MutationProb = 1.0 / comConst.NumOfTasks;    //变异率为 1/n
    inputParameter.CrossProb = 1;
    inputParameter.NumOfIterations = 100;
}

int ReadID(string id) {
    int ret = 0;
    for (int i = 0; i < id.length(); ++i) {
        if (id[i] <= '9' && id[i] >= '0')
            ret = ret * 10 + id[i] - '0';
    }
    return ret;
}

//｛从文件里面读取模型数据存入数据结构中}
void ReadFile(string XMLFile, string ResAlcFile) {
    string FilePath = "../data";
    string xmlpath = FilePath + "/" + XMLFile;
    pugi::xml_document doc;
    doc.load_file((xmlpath).c_str());  //int w = doc.load_file((xmlpath).c_str());
    pugi::xml_node root = doc.child("adag");
    for (pugi::xml_node job = root.child("job"); job; job = job.next_sibling("job")) {
        Task task = Task();
        task.length = atof(job.attribute("runtime").value());
        for (pugi::xml_node uses = job.child("uses"); uses; uses = uses.next_sibling("uses")) {
            vfile file = vfile();
            file.source = -1;
            file.FileName = uses.attribute("file").value();
            file.size = atof(uses.attribute("size").value());
            if (!strcmp(uses.attribute("link").value(), "input")) { //输入文件
                task.IFile.push_back(file);
            } else {  //输出文件
                task.OFile.push_back(file);
            }
        }
        Tasks.push_back(task);
    }
    comConst.NumOfTasks = Tasks.size();

    // read the info of relation between tasks
    for (pugi::xml_node child = root.child("child"); child; child = child.next_sibling("child")) {
        int child_index = ReadID(child.attribute("ref").value());
        for (pugi::xml_node parent = child.child("parent"); parent; parent = parent.next_sibling("parent")) {
            int parent_index = ReadID(parent.attribute("ref").value());
            Tasks[child_index].parents.insert(parent_index);
            Tasks[parent_index].children.insert(child_index);
        }
    }

    ParChildTranFileSizeSum.resize(comConst.NumOfTasks, vector<double>(comConst.NumOfTasks, 0));
    for (int i = 0; i < Tasks.size(); ++i) {
        if (Tasks[i].parents.size() == 0) continue;
        for (int p = 0; p < Tasks[i].IFile.size(); ++p) {
            string iName = Tasks[i].IFile[p].FileName;
            int flag = 0;
            for (int prnTsk: Tasks[i].parents) {
                for (int q = 0; q < Tasks[prnTsk].OFile.size(); ++q) {
                    string oName = Tasks[prnTsk].OFile[q].FileName;
                    //判断文件名相同
                    if (iName.compare(oName) == 0) {    // judge whether two file names are the same; 0: same; -1: not same
                        ParChildTranFileSizeSum[prnTsk][i] += Tasks[i].IFile[p].size;
                        Tasks[i].IFile[p].source = prnTsk;
                        flag = 1;
                        break;
                    }
                }
                if (flag == 1) break;  //if this statement exists, the "source" records the first parent task which has the same outfile name as the inputfile name, otherwise it records the last parent task
            }
        }
    }

    char line_1[1024]={0};
    ifstream iFile("../data/VMType.txt",ios::in); //建议VMType.txt中的编号也从“0”开始编号-xy，这样可以与VMTypes[]中的下标一致-xy，
    if(!iFile){
        cout<<"VMType file cannot open.";
        exit(0);
    }

    vector<int> typeIdLst; //按虚拟机类型文件中的顺序存放虚拟机类型
    while (iFile.getline(line_1, sizeof(line_1))){
        stringstream word(line_1);
        string ans;
        word>>ans;
        if(ans == "#"){
            continue;
        }

        VMType vmType;
        int type = stoi(ans);
        vmType.type = type;
        word >> vmType.ps;
        word >> vmType.bw;
        vmType.bw = vmType.bw * 1000000 / 8;  //把带宽由MBit每秒转化为byte每秒
        word >> vmType.UnitCost;
        word >> vmType.FixedCost;
        word >> vmType.ChargeUnit;
        word >> vmType.MinLeaseTime;
        word >> vmType.BootTime;
        VMTypes[type] = vmType;       //VMTypes是map数据结构
        typeIdLst.push_back(type);
        VMTypeSet.insert(type);
        for (auto tsk: Tasks) {      //初始化NoElgRsc为包括所有虚拟机类型
            tsk.NoElgRsc.insert(type);
        }
    }
    for (map<int, VMType>::iterator it = VMTypes.begin(); it != VMTypes.end(); it++) {  //测试用-xy
        if ((*it).first != (*it).second.type ) { //测试下标与存储的类型是否一致
            cout << "VMTypes is wrong!" << endl; exit(2);
        }
    }
    comConst.NumOfVMTypes = VMTypes.size();
    //read vm file to task data structure
    char line[4096] = {0};
    int taskIndex = -1;  //建议把0改成-1；-xy，已改-xy!
    int cnt = 0;
    string ResAlcFilePath = FilePath + "/" + ResAlcFile;
    ifstream fin(ResAlcFilePath, ios::in);
    if (!fin) {
        cout << "Error at open ResAlcFile" << endl; // ResAlcFile与 VMTypes文件中的资源是一行对应一行的！即ResAlcFile中第一行数据为VMTypes中第一行描述的资源可以处理的任务编号，
        exit(0);
    } else {
        while (fin.getline(line, sizeof(line))) {
            if (cnt == VMTypes.size()) {
                cout << endl << "the number of VM Types in ResAlcFile is more than that saved in VMTypes";
                exit(0);
            }
            stringstream word(line);
            while (1) {
                word >> taskIndex;
                if (word.fail()) break;
                Tasks[taskIndex].ElgRsc.insert(typeIdLst[cnt]);    //ElgRsc中存放的是虚拟机类型；
                Tasks[taskIndex].NoElgRsc.erase(typeIdLst[cnt]);   //在NoElgRsc中删除这个可以处理taskIndex的该虚拟机类型
                VMTypes[typeIdLst[cnt]].ElgTsk.insert(taskIndex);
            }
            ++cnt;
        }
        if (cnt != VMTypes.size()) {
            cout << endl << "the number of VM Types in ResAlcFile is less than that saved in VMtypes";
            exit(0);
        }
    }
    fin.close();
}

//任务合并
void TaskCombine() {
    vector<vector<int>> MT(comConst.NumOfTasks, vector<int>(0));
    vector<int> flag(comConst.NumOfTasks, 0);
    vector<int> oid(comConst.NumOfTasks, 0);
    vector<int> nid(comConst.NumOfTasks, 0);
    vector<int> TaskOrderList(comConst.NumOfTasks);
    vector<Task> TaskC;

    TaskOrderList = GenerateTopoOrder_Ran();    //根据拓扑排序生成调度顺序
    for (int i = 0; i < comConst.NumOfTasks; i++) {
        MT[i].push_back(i);                    //每一层的MT中存放的是任务i
    }

    //任务合并   从前向后找
    for (int i : TaskOrderList) {
        if (flag[i] == 0) {
            while (Tasks[i].children.size() == 1 && Tasks[*(Tasks[i].children.begin())].parents.size() == 1 && Tasks[i].ElgRsc == Tasks[*(Tasks[i].children.begin())].ElgRsc) {
                int icl = *(Tasks[i].children.begin());    //记录被合并到任务i的子任务
                MT[i].push_back(icl);                      //向MT集合中添加可与任务i合并的子任务
                MT[icl].clear();                           //将被合并的任务从所属MT中删除
                flag[icl] = 1;
                Tasks[i].length += Tasks[icl].length;      //length = length i+
                Tasks[i].children = Tasks[icl].children;   //SCi = SCi+ 将子任务中集合的子任务全部插入到父任务集合中去
                Tasks[i].OFile = Tasks[icl].OFile;         //OFLi = OFLi+ 将子任务中集合的输出文件全部插入到父任务集合中去，此处不应该是简单的替换，还需要考虑工作流的输出文件！-xy
                for (int tsk: Tasks[icl].children) {  //PRi++ = PRi++ - ti+ + ti//遍历每一个子任务的子任务
//                    replace(Tasks[tsk].parents.begin(),Tasks[tsk].parents.end(), icl, i);
                    Tasks[tsk].parents.erase(icl);
                    Tasks[tsk].parents.insert(i);
                }
            }
        }
    }

    //重新编号
    int N = 0;
    for (int i = 0; i < comConst.NumOfTasks; i++) {
        if (!MT[i].empty()) {
            oid[N] = i;             //记录新的任务编号对应的第一个老任务编号
            for (int j : MT[i]) {
                nid[j] = N;         //记录每个老任务编号对应的新任务编号
            }
            N = N + 1;
        }
    }

    //根据新编号更新任务信息
    for (int n = 0; n < N; n++) {
        Task task = Task();                    //生成一个新的任务对象
        for (int i : Tasks[oid[n]].parents) {  //i是父任务的老编号
            task.parents.insert(nid[i]);       //要把老编号转换成新编号再添加到父任务集合中去
        }

        for (int i : Tasks[oid[n]].children) {
            task.children.insert(nid[i]);
        }
        //把新的参数赋值给新生成的对象
        task.IFile = Tasks[oid[n]].IFile;
        task.OFile = Tasks[oid[n]].OFile;
        task.length = Tasks[oid[n]].length;
        task.ElgRsc = Tasks[oid[n]].ElgRsc;
        TaskC.push_back(task);
    }
    //更新任务数量
    comConst.NumOfTasks = TaskC.size();

    //计算父子任务之间传输文件大小
    ParChildTranFileSizeSum.resize(comConst.NumOfTasks, vector<double>(comConst.NumOfTasks, 0)); //建议要保存合并前父子任务之间传输文件大小-xy-2022.01.09
    for (int i = 0; i < TaskC.size(); ++i) {
        if (TaskC[i].parents.size() == 0) continue;
        for (int p = 0; p < TaskC[i].IFile.size(); ++p) {
            string iName = TaskC[i].IFile[p].FileName;
            int flag = 0;
            for (int prnTsk: TaskC[i].parents) {
                for (int q = 0; q < TaskC[prnTsk].OFile.size(); ++q) {
                    string oName = TaskC[prnTsk].OFile[q].FileName;
                    //判断文件名相同
                    if (iName.compare(oName) == 0) {
                        ParChildTranFileSizeSum[prnTsk][i] += TaskC[i].IFile[p].size;
                        TaskC[i].IFile[p].source = prnTsk;
                        flag = 1;
                        break;
                    }
                }
                if (flag == 1) break;
            }
        }
    }

//    for(int j = 0; j < Rscs.size(); ++j){
//        set<int> TemElgTsk;
//        for(int i = 0; i < Rscs[j].ElgTsk.size(); ++i){
//            TemElgTsk.insert(nid[Rscs[j].ElgTsk[i]]);
//        }
//        Rscs[j].ElgTsk.clear();
//        Rscs[j].ElgTsk.insert(Rscs[j].ElgTsk.begin(), TemElgTsk.begin(), TemElgTsk.end());
//    }
    //合并前的在TaskO中存储，合并后的还用Tasks
    Oid = oid;          //保存合并时老任务编号对应的新任务编号
    Nid = nid;          //保存合并时新任务编号对应的老任务编号
    TaskO = Tasks;      //保存合并前的任务关系
    Tasks = TaskC;      //更新合并后的任务关系
}
