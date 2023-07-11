//
// Created by 11379 on 2020-10-30.
//

#include "classAndVarDefine.h"
vector<Task> Tasks;
vector<Task> TaskO;
vector<int> Oid;
vector<int> Nid;
//vector<double> TaskDeadline;
vector<vector<int> > TaskListInLevel;
vector<vector<double> > ParChildTranFileSizeSum;
vector<int> LevelIdOfTask;
vector<chromosome> population;
//vector<vector<chromosome>> populations;
InputParamet inputParameter;
ComConst comConst;
//vector<VMType> VMTypes;
map<int, VMType> VMTypes;
set<int> VMTypeSet;
//double Budget;                      //在MOELS中代表最小成本
//double Deadline;                    //在MOELS中代表最小时间
//double Max_Cost;
//double Max_MakeSpan;
map<string, double> MaxBdg;
map<string, double> MinBdg;
map<string, double> MaxDT;
map<string, double> MinDT;
vector<set<int> > Descendants;
vector<set<int> > Ancestors;
vector<vector<int>> AdjMatrix;