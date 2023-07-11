//
// Created by 11379 on 2020-10-30.
//

#ifndef FRAME_CLASSANDVARDEFINE_H
#include "math.h"
#include <vector>
#include <string>
#include <set>
#include <list>
#include <map>
#define FRAME_CLASSANDVARDEFINE_H
#endif //FRAME_CLASSANDVARDEFINE_H

using namespace std;
//{file}
class vfile {
public:
    string FileName;     //file name
    int source;          //the source of file, -1:from the shared server; i: from task i
    double size;         //the size of file
};

//{task}
class Task {
public:
    double length;          //the length of task
    set<int> ElgRsc;        //the set of resources which are eligible to perform this task
    set<int> NoElgRsc;      //the set of resources which are not eligible to  perform this task
    set<int> parents;       //the set of parent tasks
    set<int> children;      //the set of child tasks
    vector<vfile> IFile;    //the set of input files
    vector<vfile> OFile;    //the set of output files
};

//虚拟机类型
class VMType {
public:
    int type;
    set<int> ElgTsk;               //该类型虚拟机可以处理的任务集合
    double ps, bw;                    //该类型的虚拟机的处理能力和带宽
    double UnitCost, FixedCost;       //该类型的虚拟机的单位时间成本，起租成本-xy2
    double MinLeaseTime;              //该类型的虚拟机的最小租用时间
    double ChargeUnit;                //该类型的虚拟机的最小计费单位
    double BootTime;                  //该类型的虚拟机的启动时间-xy2，BootTime
};
//虚拟机实例
class Instance {
public:
    int id;                          //实例编号
    int type;                        //该实例的类型
    list<int> AssignTskLst;          //分配给该实例处理的任务列表，存放的是对应任务的编号
    double LST, LET;                 //最早开始使用时间(即开始租用时间)，最后结束使用的时间（即归还时间）
    set<double> ITL;                 //空闲时间段
    int state;                       //0:没有被使用；1:被使用；2:被合并
};

class chromosome {
public:
    vector<int> TskSchLst;     //任务调度顺序
    vector<int> InsAlcLst;     //虚拟机分配实例列表：任务到实例的映射-xy2
    vector<int> TypeLst;       //虚拟机类型列表
    list<list<int>> Code;      //set集合中每一个list对应一个虚拟机实例，list中的第一位存放该虚拟机实例的类型，其后表示分配给该实例的任务及任务在该实例上的调度顺序-xy2
    vector<double> TskSchPart;
    vector<double> VTskSchPart;
    vector<double> TypeAlcPart;
    vector<double> VTypeAlcPart;
    vector<double> StartTime;  //任务开始时间；计算同一实例任务最小开始时间（开始租用时间）时需要。-xy2
    vector<double> EndTime;    //任务完成时间
    double FitnessValue;       //Fitness
    double MakeSpan;
    double Cost;
    double Distance;
    double Weight;            //MakeSpan的权重
    bool IsFrw;   //if IsFrw = true, TskSchLst is Forward Topological Sort, otherwise, it is Backward Topological Sort;
    bool IsIFBDI;
    bool IsLBCAI;

    // 当返回｛false,false｝时拒绝插入，当返回｛false,true｝或｛true,false｝时则接受插入，当返回｛true,true｝时没定义如何处理，要杜绝出现此类情况-xy
    bool operator<(const chromosome &otherChromosome)const {
        if (fabs(this->MakeSpan - otherChromosome.MakeSpan)<1e-6 && fabs(this->Cost - otherChromosome.Cost)<1e-6) { //根据MakeSpan和Cost去重
            return false;
        } else {
            if (fabs(this->MakeSpan - otherChromosome.MakeSpan)>=1e-6) { //先按MakeSpan从小到大排序
                return this->MakeSpan + 1e-6 < otherChromosome.MakeSpan;
            } else {                                                    //如果MakeSpan相同则按Cost从小到大排序
                return this->Cost + 1e-6 < otherChromosome.Cost;
            }
        }
    }
    bool operator==(const chromosome &otherChromosome)const {
        return (fabs(this->MakeSpan-otherChromosome.MakeSpan) < 1e-6) && (fabs(this->Cost - otherChromosome.Cost) < 1e-6);
    }

    //{桂的实现-好像有点问题}
//    bool operator<(const chromosome &otherChromosome)const {
//        return (this->MakeSpan < otherChromosome.MakeSpan) || (this->Cost < otherChromosome.Cost);
//    }
    //为什么去掉重载==会报错？-xy 应为在find中会调用==操作符
//    bool operator==(const chromosome &otherChromosome)const {
//        return (this->MakeSpan == otherChromosome.MakeSpan) && (this->Cost == otherChromosome.Cost);
//    }

    //{sort chromosome and remove the same chromosome according to fitness value}
//    bool operator<(const chromosome &otherChromosome)const {
//        return this->FitnessValue > otherChromosome.FitnessValue + 1e-6;
//    }

//    //{sort chromosome according to fitness value and remove the same chromosome according to code }
//    bool operator<(const chromosome &otherChromosome)const {
//        int flag = -1;
//        for (int i = 0; i < this->TskSchLst.size(); i++) {
//            if (this->RscAlcLst[i] != otherChromosome.RscAlcLst[i]) {
//                flag = i;
//                break;
//            }
//            if (this->TskSchLst[i] != otherChromosome.TskSchLst[i]) {
//                flag = i;
//                break;
//            }
//        }
//        if (flag == -1) {
//            return false;
//        } else {
//            if(fabs(this->FitnessValue-otherChromosome.FitnessValue)<1e-6) {
//                if(this->RscAlcLst[flag] != otherChromosome.RscAlcLst[flag]){
//                    return this->RscAlcLst[flag]<otherChromosome.RscAlcLst[flag];
//                }else {
//                    return this->TskSchLst[flag]<otherChromosome.TskSchLst[flag];
//                }
//            }else {
//                return this->FitnessValue < otherChromosome.FitnessValue;
//            }
//        }
//    }
};

class InputParamet {
public:
    double InertiaWeight;    //惯性权重
    double c1;
    double c2;
    int NumOfChromPerPop;    //每个种群的个体体数量
    int NumOfIterations;     //迭代次数
    float CrossProb;         //交叉率
    float MutationProb;      //变异率
    float theta1;            //资源实例类型概率模型的更新速率
    float theta2;            //任务调度顺序概率的更新速率
    float theta3;            //任务分组概率模型或任务分配给资源实例概率的更新速率
    float prp;               //伪随机比例规则系数
    float eta;               //启发式信息权重系数
    float RunTimeRatioOfStg1;//第一阶段运行时间比例
    float Sbx;
    float PolyMu;
};

class ComConst {
public:
    int NumOfTasks;             //任务数量
    int NumOfVMTypes;             //虚拟机数量
};

extern vector<Task> Tasks;
extern vector<Task> TaskO;
extern vector<int> Oid;
extern vector<int> Nid;
//extern vector<double> TaskDeadline;  //即任务的最晚完成时间LFT?-xy
extern vector<vector<int>> TaskListInLevel; //任务层次列表
extern vector<vector<double>> ParChildTranFileSizeSum;
extern vector<int> LevelIdOfTask;      //the level of task
extern InputParamet inputParameter;
extern ComConst comConst;
extern map<int, VMType> VMTypes;
extern set<int> VMTypeSet;
//extern double Budget;
//extern double Deadline;
//extern double Max_Cost;
//extern double Max_MakeSpan;
extern map<string, double> MaxBdg; //最大预算
extern map<string, double> MinBdg; //最小预算
extern map<string, double> MaxDT;  //最大截止时间
extern map<string, double> MinDT;  //最小截止时间

extern vector<set<int> > Descendants;
extern vector<set<int> > Ancestors;
extern vector<vector<int>> AdjMatrix;