
#include <ChrOper.h>
#include <fstream>
#include <sstream>
#include "tools.h"

void WriteOrgResultToTxt (string AlgName, string MdlName, int iterations, pair<double, vector<chromosome>>& result) {
    ofstream  outFile("../result_" + AlgName + ".txt",ios::app);
    if(!outFile){
        cout << "result_" + AlgName + ".txt" << "can not be opened" << endl; exit(1);
    } else {
        outFile << "# " << AlgName << ", "<< MdlName <<", runtime: " << result.first << ", iterations: " << iterations << endl;
        outFile << "makespan: ";
        for(int i = 0; i < result.second.size(); ++i)  outFile << result.second[i].MakeSpan << " ";
        outFile << endl << "cost: ";
        for(int i = 0; i < result.second.size(); ++i)  outFile << result.second[i].Cost << " ";
        outFile << endl;
    }
    outFile.close();
}

void WriteOrgResultToTxt (string AlgName, string MdlName, pair<double, vector<chromosome>>& result) {
    ofstream  outFile("../result_" + AlgName + ".txt",ios::app);
    if(!outFile){
        cout << "result_" + AlgName + ".txt" << "can not be opened" << endl; exit(1);
    } else {
        outFile << "# " << AlgName << ", "<< MdlName <<", runtime: " << result.first << endl;
        outFile << "makespan: ";
        for(int i = 0; i < result.second.size(); ++i)  outFile << result.second[i].MakeSpan << " ";
        outFile << endl << "cost: ";
        for(int i = 0; i < result.second.size(); ++i)  outFile << result.second[i].Cost << " ";
        outFile << endl;
    }
    outFile.close();
}


map<string, double> CalucateHV(const map<string, pair<double,vector<chromosome>>>& result){
    double MaxMakespan = 0, MaxCost = 0;
    map<string, double> HVs;
    for (const auto & alg: result) {  //遍历模型中每一个算法, 获得参考点的坐标，即完工时间和成本的最大值
        //由于帕累托解集，是按完工时间升序（成本降序）排序的，因此只需考虑第一个和最后一个值即可
        if (alg.second.second[alg.second.second.size()-1].MakeSpan - PrecisionValue > MaxMakespan) MaxMakespan = alg.second.second[alg.second.second.size()-1].MakeSpan;
        if (alg.second.second[0].Cost > MaxCost) MaxCost = alg.second.second[0].Cost;
    }
    MaxCost = 1.1 * MaxCost; MaxMakespan = 1.1 * MaxMakespan; //使第一个和最后一个点也会对HV的值产生影响
    for (const auto & alg: result) {
        double HV = 0;
        double CrnCost = MaxCost;
        for (const auto & chr : alg.second.second) { //遍历帕累托解集中的每个解
            if (CrnCost + PrecisionValue < chr.Cost) { //测试用，如果不是按成本降序排列的，则报错！-xy
                cout << endl <<"the ParetoFront is wrong in CalucateHV !"; exit(2);
            }
            HV += (CrnCost - chr.Cost) * (MaxMakespan - chr.MakeSpan);
            CrnCost = chr.Cost;
        }
        HVs[alg.first] = HV;
    }
    return HVs;
}

map<string, double> CalucateHV(const map<string, vector<pair<double,double>>>& result){
    double MaxMakespan = 0, MaxCost = 0;
    map<string, double> HVs;
    for (const auto & alg: result) {  //遍历模型中每一个算法, 获得参考点的坐标，即完工时间和成本的最大值
        //由于帕累托解集，是按完工时间升序（成本降序）排序的，因此只需考虑第一个和最后一个值即可
        if (alg.second[alg.second.size()-1].first - PrecisionValue > MaxMakespan) MaxMakespan = alg.second[alg.second.size()-1].first;
        if (alg.second[0].second > MaxCost) MaxCost = alg.second[0].second;
    }
    for (const auto & alg: result) {
        double HV = 0;
        double CrnCost = MaxCost;
        for (const auto & data: alg.second) { //遍历帕累托解集中的每个解
            if (CrnCost + PrecisionValue < data.second) { //测试用，如果不是按成本降序排列的，则报错！-xy
                cout << endl <<"the ParetoFront is wrong!"; exit(2);
            }
            HV += (CrnCost - data.second) * (MaxMakespan - data.first);
            CrnCost = data.second;
        }
        HVs[alg.first] = HV;
    }
    return HVs;
}

map<string, double> CalucateIGD(const map<string, pair<double,vector<chromosome>>>& result){
    double MAXMakespan = 0,MAXcost = 0;
    for (const auto & alg: result) {   //遍历模型中每一个算法,合并所有帕累托解，并去除相同解
        if(alg.second.second[alg.second.second.size()-1].MakeSpan - PrecisionValue > MAXMakespan) MAXMakespan = alg.second.second[alg.second.second.size()-1].MakeSpan;
        if(alg.second.second[0].Cost > MAXcost) MAXcost = alg.second.second[0].Cost;
    }
    map<string, double> IGDs;
    if (result.empty())  return IGDs;  //如果result为空，则直接返回一个空的IGDs；
    set<chromosome> mrgPop;
    for (const auto & alg: result) {   //遍历模型中每一个算法,合并所有帕累托解，并去除相同解
        mrgPop.insert(alg.second.second.begin(), alg.second.second.end());
    }
    vector<chromosome> refPF = GnrParetoFront(mrgPop); //生成参考帕累托解集refPF；
//    cout<<endl;  //测试用-xy
//    for (chromosome chr: refPF) {
//        cout<<chr.MakeSpan<<", "<<chr.Cost<<";  ";
//    }
    for (const auto & alg: result) {
        double IGD = 0;
        for (chromosome refChr: refPF) {
            double minDst = InfInt;
            for (chromosome chr: alg.second.second) {
                double dst = sqrt(pow(refChr.MakeSpan/MAXMakespan - chr.MakeSpan/MAXMakespan, 2) + pow(refChr.Cost/MAXcost - chr.Cost/MAXcost, 2));
                if ( dst + PrecisionValue < minDst ) minDst = dst;
            }
            IGD += minDst;
        }
        IGDs[alg.first] = IGD / refPF.size();
    }
    return IGDs;
}

map<string, double> CalucateIGD(const map<string, vector<pair<double,double>>>& result){
    map<string, double> IGDs;
    set<pair<double,double>> mrgSlt;
    for (const auto & alg: result) {   //遍历模型中每一个算法,合并所有帕累托解
        mrgSlt.insert(alg.second.begin(), alg.second.end());
    }
    vector<pair<double,double>> refPF = GnrParetoFront(mrgSlt); //生成参考帕累托解集refPF；

    for (const auto & alg: result) {
        double IGD = 0;
        for (pair<double,double> refSlt: refPF) {
            double minDst = InfInt;
            for (pair<double,double> crnSlt: alg.second) {
                double dst = sqrt(pow(refSlt.first - crnSlt.first, 2) + pow(refSlt.second - crnSlt.second, 2));
                if ( dst + PrecisionValue < minDst ) minDst = dst;
            }
            IGD += minDst;
        }
        IGDs[alg.first] = IGD / refPF.size();
    }
    return IGDs;
}

map<string, map<string, double>> CalucateCM(const map<string, pair<double,vector<chromosome>>>& result) {
    map<string, map<string, double>> CM;
    for (const auto & alg1: result) {
        for (const auto & alg2: result) {
            CM[alg1.first][alg2.first] = 0;
            if (alg1 == alg2) continue;
            for (chromosome chr2: alg2.second.second) {
                for (chromosome chr1: alg1.second.second) {
                    if ( (chr1.MakeSpan + PrecisionValue < chr2.MakeSpan && chr1.Cost - PrecisionValue < chr2.Cost) ||
                         (chr1.MakeSpan - PrecisionValue < chr2.MakeSpan && chr1.Cost + PrecisionValue < chr2.Cost) ) {
                        CM[alg1.first][alg2.first] += 1;
                        break;
                    }
                }
            }
            CM[alg1.first][alg2.first] =  CM[alg1.first][alg2.first] / alg2.second.second.size();
            if (CM[alg1.first][alg2.first] - PrecisionValue > 1) {
                cout << "CM[alg1.first][alg2.first] is large than 1 !"<<endl;
            }
            CM[alg2.first][alg1.first] = 0;
            for (chromosome chr1: alg1.second.second) {
                for (chromosome chr2: alg2.second.second) {
                    if ( (chr2.MakeSpan + PrecisionValue < chr1.MakeSpan && chr2.Cost - PrecisionValue < chr1.Cost) ||
                         (chr2.MakeSpan - PrecisionValue < chr1.MakeSpan && chr2.Cost + PrecisionValue < chr1.Cost) ) {
                        CM[alg2.first][alg1.first] += 1;
                        break;
                    }
                }
            }
            CM[alg2.first][alg1.first] =  CM[alg2.first][alg1.first] / alg1.second.second.size();
            if (CM[alg2.first][alg1.first] - PrecisionValue > 1) {
                cout << "CM[alg2.first][alg1.first] is larger than 1 !"<<endl;
            }

        }
    }
    return CM;
}

map<string, map<string, double>> CalucateCM(const map<string, vector<pair<double,double>>>& result) {
    map<string, map<string, double>> CM;
    for (const auto & alg1: result) {
        for (const auto & alg2: result) {
            CM[alg1.first][alg2.first] = 0;
            if (alg1 == alg2) continue;
            for (pair<double,double> Slt2: alg2.second) {
                for (pair<double,double> Slt1: alg1.second) {
                    if ( (Slt1.first + PrecisionValue < Slt2.first && Slt1.second - PrecisionValue < Slt2.second) ||
                         (Slt1.first - PrecisionValue < Slt2.first && Slt1.second + PrecisionValue < Slt2.second) ) {
                        CM[alg1.first][alg2.first] += 1;
                        break;
                    }
                }
            }
            CM[alg1.first][alg2.first] =  CM[alg1.first][alg2.first] / alg2.second.size();
            for (pair<double,double> Slt1: alg1.second) {
                for (pair<double,double> Slt2: alg2.second) {
                    if ( (Slt2.first + PrecisionValue < Slt1.first && Slt2.second - PrecisionValue < Slt1.second) ||
                         (Slt2.first - PrecisionValue < Slt1.first && Slt2.second + PrecisionValue < Slt1.second) ) {
                        CM[alg2.first][alg1.first] += 1;
                        break;
                    }
                }
            }
            CM[alg2.first][alg1.first] =  CM[alg2.first][alg1.first] / alg1.second.size();
        }
    }
    return CM;
}

//计算任务层次列表
void CalculateLevelList() {
    LevelIdOfTask.resize(comConst.NumOfTasks);
    vector<int> inDegree;      //variables for recording the number of parent tasks whose level have not been calculated;
    vector<int> ReadyTasks;    //a set for recording the index of tasks whose inDegree is equal to 0;
    inDegree.assign(comConst.NumOfTasks, 0);
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        inDegree[i] = Tasks[i].parents.size();
        if (inDegree[i] == 0) ReadyTasks.push_back(i);
    }

    int maxLevel = 0;
    while (!ReadyTasks.empty()) {
        int v = ReadyTasks[0];
        LevelIdOfTask[v] = 0;
        for (int prnTsk: Tasks[v].parents) {
            if (LevelIdOfTask[prnTsk] >= LevelIdOfTask[v]) {
                LevelIdOfTask[v] = LevelIdOfTask[prnTsk] + 1;
            }
        }
        if(LevelIdOfTask[v] + 1>maxLevel) {
            maxLevel = LevelIdOfTask[v] + 1;
            TaskListInLevel.resize(maxLevel);
        }
        TaskListInLevel[LevelIdOfTask[v]].push_back(v);
        ReadyTasks.erase(ReadyTasks.begin());
        for (int chlTsk: Tasks[v].children) {
            inDegree[chlTsk]--;
            if (inDegree[chlTsk] == 0) {
                ReadyTasks.push_back(chlTsk);
            }
        }
    }
}

bool sortPopOnDistanceByDesc(chromosome& a, chromosome& b) {
    return a.Distance > b.Distance;
}

bool SortValueByDescend(pair<int,double>& a, pair<int,double>& b) {
    return a.second > b.second;
}

bool SortValueByAscend(pair<double,double>& a, pair<double,double>& b) {
    return a.second < b.second;
}

bool sortPopByFitAsc(chromosome& a, chromosome& b) {
    return a.FitnessValue < b.FitnessValue;
}

bool sortPopByCostDesc(chromosome& a, chromosome& b){
    return a.Cost > b.Cost;
}

bool sortPopByCostAsc(chromosome& a, chromosome& b) {
    return a.Cost < b.Cost;
}

bool sortPopByMakespanDesc(chromosome& a, chromosome& b){
    return a.MakeSpan > b.MakeSpan;
}

bool sortPopByMakespanAsc(chromosome& a, chromosome& b) {
    return a.MakeSpan < b.MakeSpan;
}

bool sortPopByMsCstAsc_fstMs(chromosome& a, chromosome& b) {
    //对种群中的个体按MakeSpan从小到大排序，如果MakeSpan相同则Cost小的排在前；
    if (fabs(a.MakeSpan - b.MakeSpan) > PrecisionValue) { //如果MakeSpan不同，则MakeSpan小的排在前面
        return a.MakeSpan < b.MakeSpan;
    } else {                                              //如果MakeSpan相同，则Cost小的排在前面
        return a.Cost < b.Cost;
    }
}

//对fitness中的元素进行排序，并把位置/下标index存放在ind中，如：ind[0]存放的是fitness中最小元素在fitness中的位置/下标
void IndexSortByValueAsc(vector<int>& ind, vector<double>& fitness) {
    for (int i = 0; i < ind.size(); ++i) {
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&fitness](int v1, int v2) { return fitness[v1] < fitness[v2]; });
}

//对tskLst中的元素按其在value中的值降序排序，并把结果（tskLst中的元素）存放在ind中，如：tskLst={3,0,1}，value[0]={0.1,1.2,2.1,3.6,4.8},则排序后的ind= {3,1,0}
void IndexSortByValueDesc(vector<int>& ind, vector<int>& tskLst, vector<double>& value) {
    for (int i = 0; i < ind.size(); ++i) {
        ind[i] = tskLst[i];
    }
    sort(ind.begin(), ind.end(), [&value](int v1, int v2) { return value[v1] > value[v2]; });
}


void IndexSortByValueOnAsc(vector<int>& ind, vector<int>& value) {
    for (int i = 0; i < ind.size(); ++i) {
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&value](int v1, int v2) { return value[v1] < value[v2]; });
}


void IndexSortByValueDesc(vector<int>& ind, vector<double>& fitness) {
    for (int i = 0; i < ind.size(); i++) {
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&fitness](int v1, int v2) { return fitness[v1] > fitness[v2]; });
}

bool sortPopOnFitValueByDesc(chromosome& a, chromosome& b) {
    return a.FitnessValue > b.FitnessValue;
}

/*计算把所有任务分配给处理速度最慢的同一个虚拟机实例时的makespan*/
double estMaxMakespan_HPSO(){
    /*找处理速度最慢的虚拟机*/
    double minPs = InfInt;
    VMType vmType;
    for(auto ty : VMTypes){
        if (ty.second.ps + PrecisionValue < minPs) {
            minPs = ty.second.ps;
            vmType = ty.second;
        }
    }
    /*计算所有任务的长度*/
    double totalLength = 0;
    for(auto tsk: Tasks){
        totalLength += tsk.length;
    }
    return totalLength/vmType.ps;
}

/*计算分别为每个任务分配一个处理速度最快的虚拟机实例，同时忽略传输时间时的makespan*/
double estMinMakespan_HPSO(){
    /*找处理速度最快的虚拟机*/
    double maxPs = 0;
    VMType vmType;
    for(auto ty : VMTypes){
        if (ty.second.ps - PrecisionValue > maxPs) {
            maxPs = ty.second.ps;
            vmType = ty.second;
        }
    }
    vector<vector<double>> transferTime(comConst.NumOfTasks,vector<double>(comConst.NumOfTasks,0));
    vector<double> tskExeTime(comConst.NumOfTasks, 0), Rank_b(comConst.NumOfTasks, 0);
    /*计算每个任务的执行时间*/
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        tskExeTime[i] = Tasks[i].length / vmType.ps;
    }
    /*Rank_b中的最大值即为makespan*/
    Calculate_Rank_b(Rank_b, tskExeTime, transferTime); //transferTime中的所有元素为0，即忽略了传输时间
    double makespan = 0;
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        if (makespan + PrecisionValue < Rank_b[i])
            makespan = Rank_b[i];
    }
    return makespan;
}

//double CheapSchedule(){
//    //找处理能力最慢的虚拟机
//    double ps = InfInt;
//    VMType VM;
//    for(int j = 0; j < VMTypes.size(); ++j){
//        if(VMTypes[j].ps < ps){
//            ps = VMTypes[j].ps;
//            VM = VMTypes[j];
//        }
//    }
//    //根据拓扑排序生成一个调度顺序
//    vector<int> Order(comConst.NumOfTasks, -1);
//    Order = GenerateTopoOrder_Ran();
//    //把任务全部在此虚拟机上运行的完成时间
//    double FinalEndTime = 0;
//    for(int taskIndex : Order){
//        double Execute = Tasks[taskIndex].length / VM.ps;
//        FinalEndTime += Execute;
//    }
//    Budget = ceil(FinalEndTime/VM.ChargeUnit);   //有问题-xy
//    return FinalEndTime;
//}

//函数重载
//double CheapSchedule(chromosome &ch){
//    //找处理能力最慢的虚拟机
//    double ps = InfInt;
//    VMType VM;
//    int type ;
//    for(int j = 0; j < VMTypes.size(); ++j){
//        if(VMTypes[j].ps < ps){
//            ps = VMTypes[j].ps;
//            VM = VMTypes[j];
//            type = j;
//        }
//    }
//    ch.TskSchLst = GenerateTopoOrder_Ran();
//    //把任务全部在此虚拟机上运行的完成时间
//    double FinalEndTime = 0;
//    for(int taskIndex: ch.TskSchLst){
////        ch.StartTime[taskIndex] = FinalEndTime;
//        double Execute = Tasks[taskIndex].length / VM.ps;
//        FinalEndTime += Execute;
////        ch.EndTime[taskIndex] = FinalEndTime;
//        ch.InsAlcLst[taskIndex] = 0;
//        ch.TypeLst[taskIndex] = type;
//    }
////    Max_MakeSpan = FinalEndTime + 1;   //边界的设置，防止精度问题适应度变得很大
//    int Sch_Cost = ceil(FinalEndTime/VM.ChargeUnit);
//    Budget = Sch_Cost * VM.UnitCost;
//    ch.Cost = Sch_Cost * VM.UnitCost; ch.MakeSpan = FinalEndTime;
//    return Sch_Cost * VM.UnitCost;
//}

double CheapSchedule_xy(){
    //把所有任务分配给同一个虚拟机时处理费用最小的成本
    double minEstBdg = InfInt;
    double TotalLength = 0;
    for(auto tsk: Tasks){
        TotalLength += tsk.length;
    }
    for(auto ty : VMTypes){
        double lt = XY_MAX(ty.second.BootTime + TotalLength/ty.second.ps, ty.second.MinLeaseTime);
        double bdg = ty.second.FixedCost + ceil((lt-PrecisionValue)/ty.second.ChargeUnit) * ty.second.UnitCost;
        if (bdg +  PrecisionValue < minEstBdg ) {
            minEstBdg = bdg;
        }
    }
    return minEstBdg;
}

//函数重载
double CheapSchedule_xy(chromosome &ch){
    //把所有任务分配给处理费用最小的同一个虚拟机时的调度方案及其成本,并初始化个体
    double minEstBdg = InfInt;
    int minTy = -1;
    for(auto ty : VMTypes){
        double TotalLength = 0;
        for(auto tsk: Tasks){
            TotalLength += tsk.length;
        }
        double lt = XY_MAX(ty.second.BootTime + TotalLength/ty.second.ps, ty.second.MinLeaseTime);
        double bdg = ty.second.FixedCost + ceil((lt-PrecisionValue)/ty.second.ChargeUnit) * ty.second.UnitCost;
        if (bdg +  PrecisionValue < minEstBdg ) {
            minEstBdg = bdg;
            minTy = ty.first;
        }
    }
    ch.TskSchLst = GenerateTopoOrder_Ran();
    double CurTime = VMTypes[minTy].BootTime;
    for(int taskIndex: ch.TskSchLst){
        ch.StartTime[taskIndex] = CurTime;
        CurTime += Tasks[taskIndex].length / VMTypes[minTy].ps;
        ch.EndTime[taskIndex] = CurTime;
        ch.InsAlcLst[taskIndex] = 0;
    }
    ch.TypeLst[0] = minTy;
    ch.MakeSpan = CurTime;
    ch.Cost = minEstBdg;
    return minEstBdg;
}

double CheapSchedule_MOELS(){
    //把所有任务分配给最便宜虚拟机的一个实例时的成本
    double minCost = InfInt;
    VMType minCstTy;
    double TotalLength = 0;
    for(auto tsk: Tasks){
        TotalLength += tsk.length;
    }
    for(auto ty : VMTypes){
        if (ty.second.UnitCost + PrecisionValue < minCost) {
            minCost = ty.second.UnitCost;
            minCstTy = ty.second;
        }
    }
    double lt = XY_MAX(minCstTy.BootTime + TotalLength/minCstTy.ps, minCstTy.MinLeaseTime);
    double bdg = minCstTy.FixedCost + ceil((lt-PrecisionValue)/minCstTy.ChargeUnit) * minCstTy.UnitCost;
    return bdg;
}

//double FastSchedule(){ //HPSO原文中明确提到忽略传输时间，可以优化，需要修改！-xy
//    //找处理能力最快的虚拟机
//    double ps = 0;
//    VMType VM;
//    for(int j = 0; j < VMTypes.size(); ++j){
//        if(ps < VMTypes[j].ps){
//            ps = VMTypes[j].ps;
//            VM = VMTypes[j];
//        }
//    }
//    //租用任务个数个此类型的虚拟机
//    vector<VMType> Ins(comConst.NumOfTasks, VM);
//    //生成一个调度顺序
//    vector<int> Order(comConst.NumOfTasks, -1);
//    vector<double> ww(comConst.NumOfTasks, 0);
//    vector<vector<double>> cc(comConst.NumOfTasks, vector<double>(comConst.NumOfTasks,0));
//    vector<double> rank_b(comConst.NumOfTasks, 0);
//    W_Cal_average(ww);
//    G_Cal_average(cc);
//    Calculate_Rank_b(rank_b,ww,cc);
//    IndexSortByValueDesc(Order,rank_b);
//    //计算最小完成时间
//    vector<int> flag(Ins.size(), 0);
//    vector<set<double> > vatL;//虚拟机初始化
//    double max_EndTi = -100000;
//    vector<int> InsAlcLst(comConst.NumOfTasks, -1);
//    vector<double> EndTime(comConst.NumOfTasks, -1);
//    vector<double> StartTime(comConst.NumOfTasks, -1);
//
//    for (int j = 0; j < Ins.size(); ++j) {
//        set<double> a;
//        a.insert(Ins[j].BootTime);
//        a.insert(InfInt * 1.0);
//        vatL.push_back(a);
//    }
//    for(int taskIndex : Order) {
//        int vmIndex = -1;
//        double finalStartTime = 0;
//        double finalEndTime = 100000000000;
//        for (int j = 0; j < (*max_element(InsAlcLst.begin(),InsAlcLst.end())+2); ++j) {//后边是相同类型新的虚拟机，只需遍历一个新的
//            double rt = 0;
//            for (int i = 0; i < Tasks[taskIndex].parents.size(); ++i) {
//                int parentIndex = Tasks[taskIndex].parents[i];
//                int parentVMIndex = InsAlcLst[parentIndex];
//                double Max = EndTime[parentIndex];
//                if (j != parentVMIndex) {
//                    double sum = ParChildTranFileSizeSum[parentIndex][taskIndex];
//                    Max += sum / XY_MIN(Ins[j].bw, Ins[parentVMIndex].bw);
//                }
//                if (rt < Max) {
//                    rt = Max;
//                }
//            }
//            //执行时间
//            double ExeTime = Tasks[taskIndex].length / Ins[j].ps;
//            double StartTime = 0;
//            double EndTime = 0;
//            //搜索可用时间片段
//            set<double>::iterator pre = vatL[j].begin();
//            set<double>::iterator post = vatL[j].begin();
//            ++post;
//            //搜索可得时间段
//            while (post != vatL[j].end()) {
//                if ((*post - *pre) >= ExeTime && rt <= (*post) - ExeTime) {
//                    if (flag[j] == 0) {
//                        StartTime = XY_MAX(Ins[j].BootTime, rt);
//                        break;
//                    } else {
//                        StartTime = XY_MAX(*pre, rt);
//                        break;
//                    }
//                } else {
//                    ++pre;
//                    ++pre;
//                    ++post;
//                    ++post;
//                }
//            }
//            EndTime = StartTime + ExeTime;
//            //找最快完成
//            if (EndTime < finalEndTime) {
//                finalStartTime = StartTime;
//                finalEndTime = EndTime;
//                vmIndex = j;
//            }
//        }
//
//
//        StartTime[taskIndex] = finalStartTime;
//        EndTime[taskIndex] = finalEndTime;
//        InsAlcLst[taskIndex] = vmIndex;
//        flag[vmIndex] = 1;
//
//        //在虚拟机VMIndex的时间列表段中查找taskIndex的开始时间
//        //如果有，则不插入，并在vatL[VMIndex]删除该元素；如果没有，则插入
//        if (vatL[vmIndex].find(finalStartTime) != vatL[vmIndex].end()) {
//            vatL[vmIndex].erase(finalStartTime);
//        } else {
//            vatL[vmIndex].insert(finalStartTime);
//        }
//
//        if (vatL[vmIndex].find(EndTime[taskIndex]) != vatL[vmIndex].end()) {
//            vatL[vmIndex].erase(EndTime[taskIndex]);
//        } else {
//            vatL[vmIndex].insert(EndTime[taskIndex]);
//        }
//        max_EndTi = XY_MAX(max_EndTi, finalEndTime);
//    }
//    Deadline = max_EndTi;
//    return max_EndTi;
//}

//函数重载
//double FastSchedule(chromosome &ch){
//    //找处理能力最快的虚拟机
//    double ps = 0;
//    int type;
//    VMType VM;
//    for(int j = 0; j < VMTypes.size(); ++j){
//        if(ps < VMTypes[j].ps){
//            ps = VMTypes[j].ps;
//            VM = VMTypes[j];
//            type = j;
//        }
//    }
//    //租用任务个数个此类型的虚拟机
//    vector<VMType> Ins(comConst.NumOfTasks, VM);
//    //任务调度顺序生成
//    vector<double> ww(comConst.NumOfTasks, 0);
//    vector<vector<double>> cc(comConst.NumOfTasks, vector<double>(comConst.NumOfTasks,0));
//    vector<double> rank_b(comConst.NumOfTasks, 0);
//    W_Cal_average(ww);                      //任务平均执行时间
//    G_Cal_average(cc);                      //父子任务平均传输时间
//    Calculate_Rank_b(rank_b,ww,cc);         //rank_b 存放的是优先级
//    IndexSortByValueDesc(ch.TskSchLst,rank_b);
//    //计算最小完成时间
//    vector<set<double> > vatL;//虚拟机初始化
//    double max_EndTi = -100000;
//    vector<double> EndTime(comConst.NumOfTasks, -1);
//
//    for (int j = 0; j < Ins.size(); ++j) {
//        set<double> a;
//        a.insert(Ins[j].BootTime);
//        a.insert(InfInt * 1.0);
//        vatL.push_back(a);
//    }
//    for(int taskIndex: ch.TskSchLst) {
//        int vmIndex = -1;
//        double finalStartTime = 0;
//        double finalEndTime = 100000000000;
//        for (int j = 0; j < Ins.size(); ++j) {      //这里可以优化，虚拟机是同一种类型，不用每次都遍历全部虚拟机
//            double rt = 0;
//            for (int i = 0; i < Tasks[taskIndex].parents.size(); ++i) {
//                int parentIndex = Tasks[taskIndex].parents[i];
//                int parentVMIndex = ch.InsAlcLst[parentIndex];
//                double Max = EndTime[parentIndex];
//                if (j != parentVMIndex) {
//                    double sum = ParChildTranFileSizeSum[parentIndex][taskIndex];
//                    Max += sum / XY_MIN(Ins[j].bw, Ins[parentVMIndex].bw);
//                }
//                if (rt < Max) {
//                    rt = Max;
//                }
//            }
//            //执行时间
//            double ExeTime = Tasks[taskIndex].length / Ins[j].ps;
//            double StartTime = 0;
//            double EndTime = 0;
//            //搜索可用时间片段
//            set<double>::iterator pre = vatL[j].begin();
//            set<double>::iterator post = vatL[j].begin();
//            ++post;
//            //搜索可得时间段
//            while (post != vatL[j].end()) {
//                if((*post - *pre) >= ExeTime && rt <= (*post)-ExeTime) {
//                    StartTime = XY_MAX(*pre, rt);
//                    break;
//                } else {
//                    ++pre;
//                    ++pre;
//                    ++post;
//                    ++post;
//                }
//            }
//            EndTime = StartTime + ExeTime;
//            //找最快完成
//            if (EndTime < finalEndTime) {
//                finalStartTime = StartTime;
//                finalEndTime = EndTime;
//                vmIndex = j;
//            }
//        }
//
//        ch.InsAlcLst[taskIndex] = vmIndex;
//        ch.TypeLst[vmIndex] = type;
//        EndTime[taskIndex] = finalEndTime;
//
//        //在虚拟机VMIndex的时间列表段中查找taskIndex的开始时间
//        //如果有，则不插入，并在vatL[VMIndex]删除该元素；如果没有，则插入
//        if (vatL[vmIndex].find(finalStartTime) != vatL[vmIndex].end()) {
//            vatL[vmIndex].erase(finalStartTime);
//        } else {
//            vatL[vmIndex].insert(finalStartTime);
//        }
//
//        if (vatL[vmIndex].find(EndTime[taskIndex]) != vatL[vmIndex].end()) {
//            vatL[vmIndex].erase(EndTime[taskIndex]);
//        } else {
//            vatL[vmIndex].insert(EndTime[taskIndex]);
//        }
//        max_EndTi = XY_MAX(max_EndTi, finalEndTime);
//    }
//    //把虚拟机类型变为不为-1的随机一个数
//    for(int j = *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end())+1; j < comConst.NumOfTasks; ++j){
//        ch.TypeLst[j] = rand() % comConst.NumOfVMTypes;
//    }
//
//    Deadline = max_EndTi;
//    ch.MakeSpan = max_EndTi;
//    return max_EndTi;
//}

/*为每个任务分配一个处理速度最快的虚拟机实例，同时考虑了传输时间时的makespan*/
double FastSchedule_xy(vector<double>& rank_b){
    //找处理速度最快的虚拟机
    double ps_max = 0;
    int ty_max = -1;
    for(auto ty: VMTypes){
        if(ps_max < ty.second.ps){
            ps_max = ty.second.ps;
            ty_max = ty.first;
        }
    }
    chromosome chrom;
    initChrom(chrom);
    chrom.MakeSpan = -1;
    vector<set<double>> ITL;  //虚拟机初始化
    int rscNumUsed = 0;       //记录已使用的实例数量
    for (int j = 0; j < comConst.NumOfTasks; ++j) {
        set<double> a;
        a.insert(VMTypes[ty_max].BootTime);
        a.insert(InfInt * 1.0);
        ITL.push_back(a);
    }
    for (int j = 0; j < comConst.NumOfTasks; ++j ) {
        chrom.TypeLst[j] = ty_max;
    }
    IndexSortByValueDesc(chrom.TskSchLst,rank_b);
    for(int i = 0; i < comConst.NumOfTasks; ++i) {
        int tskId = chrom.TskSchLst[i];
        int vmId = -1;
        double fnlStartTime = 0;
        double fnlEndTime = InfInt;
        //虚拟机是同一种类型，不用每次都遍历全部虚拟机, 且使用的虚拟机的实例编号是连续的，实例j使用了那么编号小于j的实例也必被使用-xy
        for (int j = 0; j < rscNumUsed + 1; ++j) {
            double rt = 0;
            for (int prnId: Tasks[tskId].parents) {
                int prnVMId = chrom.InsAlcLst[prnId];
                double fft = chrom.EndTime[prnId];
                if (j != prnVMId) {
                    fft += ParChildTranFileSizeSum[prnId][tskId]/ VMTypes[ty_max].bw;
                }
                if (rt + PrecisionValue < fft) {
                    rt = fft;
                }
            }
            double ExeTime = Tasks[tskId].length / VMTypes[ty_max].ps;
            double StartTime = FindIdleTimeSlot(ITL[j],ExeTime,rt);
            double EndTime = StartTime + ExeTime;
            if (EndTime + PrecisionValue < fnlEndTime) {
                fnlStartTime = StartTime;
                fnlEndTime = EndTime;
                vmId = j;
            }
        }
        chrom.StartTime[tskId] = fnlStartTime;
        chrom.EndTime[tskId] = fnlEndTime;
        chrom.MakeSpan = XY_MAX(chrom.MakeSpan,fnlEndTime);
        chrom.InsAlcLst[tskId] = vmId;
        rscNumUsed = XY_MAX(rscNumUsed,vmId+1);
        UpdateITL(ITL[vmId], fnlStartTime, fnlEndTime);
    }
    chrom.Cost = 0;
    for (int k = 0; k < rscNumUsed; ++k) {
        if (ITL[k].size() > 2 || fabs(*(ITL[k].begin())- VMTypes[ty_max].BootTime) > PrecisionValue) {
            set<double>::iterator pFrt = ITL[k].begin(), pLst = ITL[k].end();
            pLst--;  pLst--;
            double frtTime = VMTypes[ty_max].BootTime, lstTime = *pLst;
            if (fabs(VMTypes[ty_max].BootTime - *pFrt) < PrecisionValue) {//第1个点等于BootTime，那么第二个点才是真正开始执行任务的起点
                pFrt++;
                frtTime = *pFrt;
            }
            double LeaseTime = XY_MAX(lstTime - frtTime + VMTypes[ty_max].BootTime, VMTypes[ty_max].MinLeaseTime);
            chrom.Cost = chrom.Cost + VMTypes[ty_max].FixedCost + VMTypes[ty_max].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[ty_max].ChargeUnit);
        }
    }
//    double ct = chrom.Cost, ms = chrom.MakeSpan;
//    DecodeChrom(chrom, true);
//    if (fabs(chrom.MakeSpan - ms) > PrecisionValue || fabs(chrom.Cost - ct) > PrecisionValue )
//        cout << "it is wrong in FastSchedule_xy";
    return  chrom.MakeSpan;
}

//函数重载
double FastSchedule_xy(vector<double>& rank_b, chromosome& chrom){
    //找处理速度最快的虚拟机
    double ps_max = 0;
    int ty_max = -1;
    for(auto ty: VMTypes){
        if(ps_max < ty.second.ps){
            ps_max = ty.second.ps;
            ty_max = ty.first;
        }
    }
//    chromosome chrom;
    initChrom(chrom);
    chrom.MakeSpan = -1;
    vector<set<double>> ITL;  //虚拟机初始化
    int rscNumUsed = 0;       //记录已使用的实例数量
    for (int j = 0; j < comConst.NumOfTasks; ++j) {
        set<double> a;
        a.insert(VMTypes[ty_max].BootTime);
        a.insert(InfInt * 1.0);
        ITL.push_back(a);
    }
    for (int j = 0; j < comConst.NumOfTasks; ++j ) {
        chrom.TypeLst[j] = ty_max;
    }
    IndexSortByValueDesc(chrom.TskSchLst,rank_b);
    for(int i = 0; i < comConst.NumOfTasks; ++i) {
        int tskId = chrom.TskSchLst[i];
        int vmId = -1;
        double fnlStartTime = 0;
        double fnlEndTime = InfInt;
        //虚拟机是同一种类型，不用每次都遍历全部虚拟机, 且使用的虚拟机的实例编号是连续的，实例j使用了那么编号小于j的实例也必被使用-xy
        for (int j = 0; j < rscNumUsed + 1; ++j) {
            double rt = 0;
            for (int prnId: Tasks[tskId].parents) {
                int prnVMId = chrom.InsAlcLst[prnId];
                double fft = chrom.EndTime[prnId];
                if (j != prnVMId) {
                    fft += ParChildTranFileSizeSum[prnId][tskId]/ VMTypes[ty_max].bw;
                }
                if (rt + PrecisionValue < fft) {
                    rt = fft;
                }
            }
            double ExeTime = Tasks[tskId].length / VMTypes[ty_max].ps;
            double StartTime = FindIdleTimeSlot(ITL[j],ExeTime,rt);
            double EndTime = StartTime + ExeTime;
            if (EndTime + PrecisionValue < fnlEndTime) {
                fnlStartTime = StartTime;
                fnlEndTime = EndTime;
                vmId = j;
            }
        }
        chrom.StartTime[tskId] = fnlStartTime;
        chrom.EndTime[tskId] = fnlEndTime;
        chrom.MakeSpan = XY_MAX(chrom.MakeSpan,fnlEndTime);
        chrom.InsAlcLst[tskId] = vmId;
        rscNumUsed = XY_MAX(rscNumUsed,vmId+1);
        UpdateITL(ITL[vmId], fnlStartTime, fnlEndTime);
    }
    chrom.Cost = 0;
    for (int k = 0; k < rscNumUsed; ++k) {
        if (ITL[k].size() > 2 || fabs(*(ITL[k].begin())- VMTypes[ty_max].BootTime) > PrecisionValue) {
            set<double>::iterator pFrt = ITL[k].begin(), pLst = ITL[k].end();
            pLst--;  pLst--;
            double frtTime = VMTypes[ty_max].BootTime, lstTime = *pLst;
            if (fabs(VMTypes[ty_max].BootTime - *pFrt) < PrecisionValue) {//第1个点等于BootTime，那么第二个点才是真正开始执行任务的起点
                pFrt++;
                frtTime = *pFrt;
            }
            double LeaseTime = XY_MAX(lstTime - frtTime + VMTypes[ty_max].BootTime, VMTypes[ty_max].MinLeaseTime);
            chrom.Cost = chrom.Cost + VMTypes[ty_max].FixedCost + VMTypes[ty_max].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[ty_max].ChargeUnit);
        }
    }
    return  chrom.MakeSpan;
}

double FastSchedule_MOELS(vector<double>& rank_b){
    //找处理速度最快的虚拟机
    double ps_max = 0;
    int ty_max = -1;
    for(auto ty: VMTypes){
        if(ps_max < ty.second.ps){
            ps_max = ty.second.ps;
            ty_max = ty.first;
        }
    }
    chromosome chrom;
    initChrom(chrom);
    //为每个任务分配一个处理速度最快的虚拟机实例
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        chrom.InsAlcLst[i] = i;
        chrom.TypeLst[i] = ty_max;
    }
    IndexSortByValueDesc(chrom.TskSchLst,rank_b);
    DecodeChrom(chrom);
    return chrom.MakeSpan;
}

//计算每个任务的子孙任务
void CalculateDescendants() {
    Descendants.resize(comConst.NumOfTasks);
    for(int i = TaskListInLevel.size()-2; i >= 0; --i) {
        for(int taskId: TaskListInLevel[i]) {
            for(int childId: Tasks[taskId].children) {
                Descendants[taskId].insert(childId);
                Descendants[taskId].insert(Descendants[childId].begin(),Descendants[childId].end());
            }
        }
    }
}

//计算每个任务的祖先任务
void CalculateAncestors() {
    Ancestors.resize(comConst.NumOfTasks);
    for(int i = 1; i < TaskListInLevel.size(); ++i) {//可以从i=1开始，-xy3  //已改 -G
        for(int taskId: TaskListInLevel[i]) {
            for(int parentId: Tasks[taskId].parents) {
                Ancestors[taskId].insert(parentId);
                Ancestors[taskId].insert(Ancestors[parentId].begin(),Ancestors[parentId].end());
            }
        }
    }
}

vector<int> GenerateTopoOrder_Ran() {
    vector<int> result;
    vector<int> inDegree(comConst.NumOfTasks, 0);  //保存节点入度
    list<int> stk;       //保存当前入度为0的节点编号,即就绪任务的编号,改用链表实现了，效率应更高-xy
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        inDegree[i] = Tasks[i].parents.size();
        if (inDegree[i] == 0) stk.push_back(i);
    }

    while (!stk.empty()) {
        int rnd = rand() % stk.size();
        list<int>::iterator iter =  stk.begin();
        advance(iter,rnd);
        int tskId = *iter;
        stk.erase(iter);
        //遍历tskId的子任务（后继节点）
        for (int childId: Tasks[tskId].children) {
            --inDegree[childId];
            if (inDegree[childId] == 0) stk.push_back(childId);
        }
        result.push_back(tskId);
    }
    if (result.size() != comConst.NumOfTasks) {
        cout << "图中存在环路，不能进行拓扑排序！" << endl;
    }
    return result;
}

vector<int> GenerateTopoOrder_lev() {
    vector<int> tskLst;
    for (int i = 0; i < TaskListInLevel.size(); i++) {
        vector<int> tem = TaskListInLevel[i];
        random_shuffle(tem.begin(), tem.end());              //随机重排i层次的任务顺序
        tskLst.insert(tskLst.end(), tem.begin(), tem.end()); //i层次的任务添加到a的尾部， 改成了整体添加，效率应更高-xy
    }
    return tskLst;
}

vector<int> GenerateTopoOrder_lev_bt(vector<double>& rank_bt) {
    vector<int> tskLst;
    for (int i = 0; i < TaskListInLevel.size(); i++) {
        vector<int> TskSortedLst(TaskListInLevel[i].size());
        IndexSortByValueDesc(TskSortedLst, TaskListInLevel[i], rank_bt); //对i层中的任务按其rank_bt值降序排列；
        tskLst.insert(tskLst.end(), TskSortedLst.begin(), TskSortedLst.end()); //i层次的任务添加到a的尾部
    }
    return tskLst;
}

//计算临接矩阵，如果AdjMatrix[i][j]=1表示任务j是任务i的子任务 MU LGA HGA
void CalucateAdjMatrix(){
    AdjMatrix.resize(comConst.NumOfTasks);
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        AdjMatrix[i].resize(comConst.NumOfTasks);
        for (int j: Tasks[i].children) {
            AdjMatrix[i][j] = 1;
        }
    }
}

//把UT中其父任务都已被调度的任务移到RT中
void UTtoRT(set<int>& UT, vector<int>& RT, vector<vector<int>>& linkMatrix) {
    for (set<int>::iterator it = UT.begin(); it != UT.end();) {
        int sum = 0;
        //查看是否有父任务
        for (int i = 0; i < comConst.NumOfTasks; ++i) {
            sum += linkMatrix[i][*it];
        }
        if (!sum) {
            RT.push_back(*it);
            it = UT.erase(it);//删除当前指针指向的元素，自动返回下一元素的指针。
        } else
            ++it;
    }
}

vector<double> CalculateDeadline(double Deadline){
    //平均带宽
    double sum = 0;
    for(int j = 0; j < VMTypes.size(); ++j){
        sum += VMTypes[j].bw;
    }
    double AvrBw = sum / VMTypes.size();
    //最快虚拟机处理能力
    double FastVM = 0;
    for(int j = 0; j < VMTypes.size(); ++j){
        if(FastVM < VMTypes[j].ps)
            FastVM = VMTypes[j].ps;
    }
    vector<double> TaskDeadline(comConst.NumOfTasks, -1);
    //计算每个任务的截止日期 --- 把没有子任务的任务截止时间设置为最大的截止时间
    for (int tsk: TaskListInLevel[TaskListInLevel.size()-1]) {
        TaskDeadline[tsk] = Deadline;
    }
    for(int i = TaskListInLevel.size()-2; i >=0; --i){
        for (int tsk: TaskListInLevel[i]) {
            TaskDeadline[tsk] = Deadline;
            for (int child: Tasks[tsk].children) {
                double CompareTime = TaskDeadline[child] - Tasks[child].length/FastVM - ParChildTranFileSizeSum[tsk][child]/AvrBw;
                if(TaskDeadline[tsk] - PrecisionValue > CompareTime){
                    TaskDeadline[tsk] = CompareTime;
                }
            }
        }
    }
    return TaskDeadline;
}

vector<double> calculateTskLFT(double Deadline){
    double sumBw = 0, maxBw = 0, maxPs = 0;
    for(auto vmTy: VMTypes){
        sumBw += vmTy.second.bw;
//        if (maxBw + PrecisionValue < vmTy.second.bw) {
//            maxBw = vmTy.second.bw;
//        }
        if (maxPs + PrecisionValue < vmTy.second.ps) {
            maxPs = vmTy.second.ps;
        }
    }
    double AvrBw = sumBw / VMTypes.size();  //平均带宽, HPSO原文中所有虚拟机的带宽是相同的，因此此处也可考虑使用最大带宽
    vector<double> tskLFT(comConst.NumOfTasks, -1);
    //计算每个任务的截止日期 --- 把没有子任务的任务截止时间设置为最大的截止时间
    for (int tsk: TaskListInLevel[TaskListInLevel.size()-1]) {
        tskLFT[tsk] = Deadline;
    }
    for(int i = TaskListInLevel.size()-2; i >=0; --i){
        for (int tsk: TaskListInLevel[i]) {
            tskLFT[tsk] = Deadline;
            for (int child: Tasks[tsk].children) {
                double CompareTime = tskLFT[child] - Tasks[child].length/maxPs - ParChildTranFileSizeSum[tsk][child]/AvrBw;
                if(tskLFT[tsk] - PrecisionValue > CompareTime){
                    tskLFT[tsk] = CompareTime;
                }
            }
        }
    }
    return tskLFT;
}

bool BetterParticle(chromosome& ch1, chromosome& ch2, double& Deadline) {
    if(ch1.MakeSpan - PrecisionValue < Deadline && ch2.MakeSpan - PrecisionValue > Deadline){        //ch1可行，ch2不可行
        return false;
    } else if (ch2.MakeSpan - PrecisionValue < Deadline && ch1.MakeSpan - PrecisionValue > Deadline){ //ch2可行，ch1不可行
        ch1 = ch2;
        return true;
    } else if (ch1.MakeSpan - PrecisionValue < Deadline && ch2.MakeSpan - PrecisionValue < Deadline){ //ch1和ch2都可行
        if(ch1.Cost < ch2.Cost) {
            return false;
        } else {
            ch1 = ch2;
            return true;
        }
    } else {   //ch1和ch2不可行
        if(ch1.MakeSpan < ch2.MakeSpan) {
            return false;
        } else {
            ch1 = ch2;
            return true;
        }
    }
}
chromosome BetterPartcle(chromosome& p1, chromosome& p2, double& TotalDeadline) {
    if(p1.MakeSpan <= TotalDeadline && p2.MakeSpan > TotalDeadline){
        return p1;
    } else if(p1.MakeSpan > TotalDeadline && p2.MakeSpan <= TotalDeadline){
        return p2;
    } else {//此处还要区分两个都是可行个体，还是都是不可行个体-xy
        if(p1.Cost < p2.Cost) return p1;
        else return p2;
    }
}

//{generate a random number in [Start,End) }
double RandomDouble(int start, int end) {
    double ret = rand() % ((end - start) * 1000) / 1000.0 + start;
    return ret;
}

//锦标赛法选择2个数字,即两个不同的个体(此处选适应值大的个体,此方法是适应值越大越好) --- 四选二
void selection_Tournament(vector<chromosome>& population, int& parent_1, int& parent_2) {
    int p1 = rand() % inputParameter.NumOfChromPerPop;
    int p2 = rand() % inputParameter.NumOfChromPerPop;
    while (p1 == p2)
        p2 = rand() % inputParameter.NumOfChromPerPop;
    if (population[p1].FitnessValue > population[p2].FitnessValue)
        parent_1 = p1;
    else
        parent_1 = p2;
    parent_2 = parent_1;
    //再用锦标赛法选出一个与parent1不同的parent2
    while (parent_2 == parent_1) {
        p1 = rand() % inputParameter.NumOfChromPerPop;
        p2 = rand() % inputParameter.NumOfChromPerPop;
        while (p1 == p2)
            p2 = rand() % inputParameter.NumOfChromPerPop;
        if (population[p1].FitnessValue > population[p2].FitnessValue)
            parent_2 = p1;
        else
            parent_2 = p2;
    }
}

int selectChrom_Tournament(int& NumOfChromPerPop) {
    /* 锦标赛法选择1个数字(个体下标)，要求种群中的个体从优到劣排序，即个体下标越小个体越优*/
    int p1 = rand() % NumOfChromPerPop;
    int p2 = rand() % NumOfChromPerPop;
    while (p1 == p2)
        p2 = rand() % NumOfChromPerPop;
    return XY_MIN(p1, p2);
}

//从合成的种群中选取pareto非支配解集集合
vector<chromosome> ParetoFront(vector<chromosome>& population){
    vector<chromosome> NewPopulation;
    vector<chromosome> NDVec;
    //从population中删除NDSet中的元素
    //1.判断NewPopulation是否满足种群容量，满足，则退出不再进行计算，否则继续添加非支配解
    while(NewPopulation.size() < inputParameter.NumOfChromPerPop){
        //2.从population中删除已添加到NDSet中的元素
        for(chromosome chrom : NDVec){
            population.erase(find(population.begin(), population.end(), chrom));
        }
        NDVec.clear();
        //每次计算一层的非支配解以及拥挤距离排序
        SelectParetoFront(population, NDVec);
        NewPopulation.insert(NewPopulation.end(), NDVec.begin(), NDVec.end());
    }
    return NewPopulation;
}

void SelectParetoFront(vector<chromosome>& population, vector<chromosome>& NDVec){
    for(int i = 0; i < population.size(); ++i) {
        bool dominted_sign = true;
        if (NDVec.empty()) {
            NDVec.push_back(population[i]);
        } else {
            //如果population[i]两个目标都大于或等于NDVec[j]则不添加到NDVec，如果两个目标都比NDVec[j]小则添加到NDVec，并删除NDVec[j]
            //如果population[i]两个目标其中一个比NDVec[j]大，另一个比NDVec[j]小，则添加到NDVec，不删除NDVec[j]
            for (int j = 0; j < NDVec.size(); ++j) {
                if ((NDVec[j].MakeSpan <= population[i].MakeSpan && NDVec[j].Cost <= population[i].Cost)) {
                    dominted_sign = false;
                    break;
                } else if (population[i].MakeSpan <= NDVec[j].MakeSpan && population[i].Cost <= NDVec[j].Cost) { //相等条件上一个if已处理
                    vector<chromosome>::iterator iter = NDVec.begin();
                    advance(iter, j);
                    NDVec.erase(iter);
                    j--;                //删除后指向了下一个位置
                }
            }
            if (dominted_sign) {
                NDVec.push_back(population[i]);
            }
        }
    }
    //根据拥挤度对NDSet容器中的个体进行排序
    //1.根据成本进行升序排序
    sort(NDVec.begin(), NDVec.end(), sortPopByCostAsc);
    //计算每个个体的拥挤度
    //1.1 取出最小的和最大的适应度值(成本)作为最小拥挤度和最大拥挤度
    double MinCostDis = NDVec[0].Cost;
    double MaxCostDis = NDVec[NDVec.size()-1].Cost;
    //1.2 取出最小的和最大的适应度值(工期)作为最小拥挤度和最大拥挤度
    double MinMakeSpanDis = NDVec[NDVec.size()-1].MakeSpan;
    double MaxMakeSpanDis = NDVec[0].MakeSpan;
    //1.3 把目标值最小的和最大拥挤距离设为无穷大，其他初始为0
    for(int i = 0 ; i < NDVec.size(); ++i){
        NDVec[i].Distance = 0;
    }
    NDVec[0].Distance = InfInt;
    NDVec[NDVec.size()-1].Distance = InfInt;
    //1.4 计算中间部分的拥挤距离
    for(int i = 1; i < (NDVec.size()-1); ++i){
        NDVec[i].Distance += (NDVec[i+1].Cost - NDVec[i-1].Cost) / (MaxCostDis - MinCostDis) + (NDVec[i-1].MakeSpan - NDVec[i+1].MakeSpan) / (MaxMakeSpanDis - MinMakeSpanDis);
    }
    //2.按拥挤距离从大到小排序
    sort(NDVec.begin(), NDVec.end(), sortPopOnDistanceByDesc);
}

vector<chromosome> GnrParetoFront(set<chromosome>& pop){
    //要求pop是set结构，其中的个体是按makespan从小到大排序的（如果makespan相等则cost小的排在前），这样任何个体不可能支配排在其前面的个体！
    //因此当前要考察的当前个体的ms必大于等于已在PF中的任何个体的ms，同时其也不可能支配已在PF中的任何个体；
    vector<chromosome> PF;
    PF.push_back(*(pop.begin()));
    pop.erase(pop.begin());
    for ( set<chromosome>::iterator it = pop.begin(); it != pop.end(); ) {
        if ( ((*it).MakeSpan + PrecisionValue < PF[PF.size()-1].MakeSpan) ||
             (fabs((*it).MakeSpan - PF[PF.size()-1].MakeSpan) < PrecisionValue && (*it).Cost + PrecisionValue <PF[PF.size()-1].Cost ) ) {
            cout << endl << "the pop does not meet the sorting requirements!";
        } //测试用-xy
        if ((*it).Cost + PrecisionValue < PF[PF.size()-1].Cost) { //找到第一个成本小于PF中的最后一个个体，即找到一个非支配解
            PF.push_back(*it);
            pop.erase(it++); //等效于 set<chromosome>::iterator temIt = it; it++; pop.erase(temIt); //要先把指针移到下一个元素，然后删除当前要删除的元素！！
        } else {
            it++;
        }
    }
    return PF;
}

vector<chromosome> GnrParetoFront(vector<chromosome>& pop){
    //要求pop是vector结构，其中的个体是按makespan从小到大排序的（如果makespan相等则cost小的排在前），这样任何个体不可能支配排在其前面的个体！
    //因此当前要考察的当前个体pop[n]的ms必大于等于已在PF中的任何个体的ms，同时其也不可能支配已在PF中的任何个体；
    vector<chromosome> PF;
    PF.push_back(pop[0]);
    pop.erase(pop.begin());
    for ( int n = 0; n < pop.size(); ++n) {
        if ( (pop[n].MakeSpan + PrecisionValue < PF[PF.size()-1].MakeSpan) ||
             (fabs(pop[n].MakeSpan - PF[PF.size()-1].MakeSpan) < PrecisionValue && pop[n].Cost + PrecisionValue <PF[PF.size()-1].Cost ) ) {
            cout << endl << "the pop does not meet the sorting requirements!";
        }
        if (pop[n].Cost + PrecisionValue < PF[PF.size()-1].Cost) { //找到第一个成本小于PF中最后一个个体，即找到一个非支配解
            PF.push_back(pop[n]);
            pop.erase(pop.begin() + n);
            --n;
        }
    }
    return PF;
}

vector<pair<double, double>> GnrParetoFront(set<pair<double, double>>& solutions){
    //要求solutions是set结构，其中的解（数据对）是按makespan（first）从小到大排序的（如果makespan相等则cost（second）小的排在前），这样任何解不可能支配排在其前面的解！
    //因此当前要考察的当前解的ms必大于等于已在PF中的任何解的ms，同时其也不可能支配已在PF中的任何解；
    vector<pair<double, double>> PF;
    PF.push_back(*(solutions.begin()));
    solutions.erase(solutions.begin());
    for ( set<pair<double, double>>::iterator it = solutions.begin(); it != solutions.end(); ) {
        if ( ((*it).first + PrecisionValue < PF[PF.size()-1].first) ||
             (fabs((*it).first - PF[PF.size()-1].first) < PrecisionValue && (*it).second + PrecisionValue <PF[PF.size()-1].second ) ) { //测试用-xy
            cout << endl << "the pop does not meet the sorting requirements!";
        }
        if ((*it).second + PrecisionValue < PF[PF.size()-1].second) { //找到第一个成本小于PF中的最后一个个体，即找到一个非支配解
            PF.push_back(*it);
            solutions.erase(it++);
//            set<chromosome>::iterator temIt = it; it++; pop.erase(temIt);
        } else {
            it++;
        }
    }
    return PF;
}

vector<chromosome> GnrParetoFront(list<chromosome>& pop){
    //由于pop是vector结构，因此为提高效率个体在放到PFVec中的同时从pop中删除该个体
    int num = pop.size();
    list<chromosome> orgPop=pop;
    vector<chromosome> PFVec;
    PFVec.push_back(*(pop.begin()));
    pop.erase(pop.begin());
    //由于pop是list结构，其中的个体是按makespan从小到大排序的（如果makespan相等则cost小的在前）
    //因此当前要考察的当前个体（*it）的ms必大于等于已在NDVec中的任何个体的ms，同时其也不可能支配已在NDVec中的任何个体；
    for(list<chromosome>::iterator it = pop.begin(); it != pop.end(); ++it) {
        bool dominated_sign = true;
        for (int j = PFVec.size() - 1; j >= 0; --j) {
            if ((*it).MakeSpan + PrecisionValue < PFVec[j].MakeSpan ) {
                cout.precision(2);
                cout<<fixed<<endl<<"pop: "<<(*it).MakeSpan<< ", "<<(*it).Cost<<";  ND: "<<PFVec[j].MakeSpan <<", "<<PFVec[j].Cost<<endl<<endl;
                for (int k=0; k<PFVec.size(); ++k) {
                    cout<<fixed<<PFVec[k].MakeSpan<<", "<<PFVec[k].Cost<<";  ";
                }
                cout<<endl<<endl;
                for(list<chromosome>::iterator it1 = pop.begin(); it1 != pop.end();++it1 ) {
                    cout<<fixed<<(*it1).MakeSpan<<", "<<(*it1).Cost<<";  ";
                }
                cout<<endl<<endl;
                for(list<chromosome>::iterator it1 = orgPop.begin(); it1 != orgPop.end();++it1 ) {
                    cout<<fixed<<(*it1).MakeSpan<<", "<<(*it1).Cost<<";  ";
                }
                exit(5);
            }
            if ( PFVec[j].Cost - PrecisionValue < (*it).Cost) {
                dominated_sign = false;
                break;
            } else if ((*it).MakeSpan + PrecisionValue < PFVec[j].MakeSpan && (*it).Cost + PrecisionValue < PFVec[j].Cost) { //相等条件上一个if已处理
                cout<<endl<<"it is wrong in GnrParetoFront, the sort of original Pop is wrong!";
                pop.push_back(PFVec[j]);
                PFVec.erase(PFVec.begin()+j);
            }
        }
        if (dominated_sign) {
            PFVec.push_back(*it);
            pop.erase(it);
            --it;
        }
    }
    if (num != pop.size() + PFVec.size()) {
        cout<<endl<<"number is not equal in GnrParetoFront";
        exit(3);
    }
    return PFVec;
}

vector<chromosome> GnrParetoFront_setOld(set<chromosome>& pop){
    set<chromosome> orgPop = pop;
    vector<chromosome> PFVec;
    set<chromosome>::iterator it = pop.begin();
    PFVec.push_back(*it);
//    set<chromosome>::iterator it1 = pop.begin();
    ++it;
//    pop.erase(it1);
    //由于pop是set结构，其中的个体是按makespan从小到大排序的（如果makespan相等则cost小的在前）
    //因此当前要考察的当前个体（*it）的ms必大于等于已在NDVec中的任何个体的ms，同时其也不可能支配已在NDVec中的任何个体；
    for(; it != pop.end(); ) {
        bool dominated_sign = true;
        for (int j = PFVec.size() - 1; j >= 0; --j) {
            if ((*it).MakeSpan + PrecisionValue < PFVec[j].MakeSpan ) { //测试用-xy
                cout.precision(2);
                cout<<fixed<<endl<<"pop: "<<(*it).MakeSpan<< ", "<<(*it).Cost<<";  ND: "<<PFVec[j].MakeSpan <<", "<<PFVec[j].Cost<<endl<<endl;
                for (int k=0; k<PFVec.size(); ++k) {
                    cout<<fixed<<PFVec[k].MakeSpan<<", "<<PFVec[k].Cost<<";  ";
                }
                cout<<endl<<endl;
                for(set<chromosome>::iterator it1 = pop.begin(); it1 != pop.end();++it1 ) {
                    cout<<fixed<<(*it1).MakeSpan<<", "<<(*it1).Cost<<";  ";
                }
                cout<<endl<<endl;
                for(set<chromosome>::iterator it1 = orgPop.begin(); it1 != orgPop.end();++it1 ) {
                    cout<<fixed<<(*it1).MakeSpan<<", "<<(*it1).Cost<<";  ";
                }
                cout<<endl<<"it is wrong in GnrParetoFront";
                exit(5);
            }
            if (PFVec[j].Cost - PrecisionValue < (*it).Cost) { //由于在PFVec中的任何个体的ms必小于等于（*it）的ms，因此makespan可以不用判断-xy，只需判断cost是否小于等于即可
                dominated_sign = false;
                break;
            } else if ((*it).MakeSpan + PrecisionValue < PFVec[j].MakeSpan && (*it).Cost + PrecisionValue < PFVec[j].Cost) { //测试用-xy
                cout<<endl<<"it is wrong in GnrParetoFront";
                exit(6);
            }
        }
        if (dominated_sign) {
            PFVec.push_back(*it);
            it++;
//            set<chromosome>::iterator temIt = it; it++; pop.erase(temIt);
//            pop.erase(it++);
        } else {
            it++;
        }
    }
//    for (int i=0; i<PFVec.size()-1; ++i) {  //测试用-xy
//        if (PFVec[i].MakeSpan - PrecisionValue > PFVec[i+1].MakeSpan || PFVec[i].Cost + PrecisionValue < PFVec[i+1].Cost) {
//            cout<<endl<<"PFVec is wrong in sort";
//            exit(6);
//        }
//    }
    return PFVec;
}

vector<chromosome> GnrParetoFront_vct_old(vector<chromosome>& pop){
    //由于pop是vector结构，因此为提高效率个体在放到PFVec中的同时从pop中删除该个体
    vector<chromosome> orgPop = pop;
    vector<chromosome> PFVec;
    PFVec.push_back(pop[0]);
    pop.erase(pop.begin());
    //由于pop是vector结构，其中的个体是按makespan从小到大排序的（如果makespan相等则cost小的在前）*/
    //因此当前要考察的当前个体pop[n]的ms必大于等于已在PFVec中的任何个体的ms，同时其也不可能支配已在PFVec中的任何个体；
    for(int n = 0; n < pop.size(); ++n) {
        bool dominated_sign = true;
        for (int j = PFVec.size() - 1; j >= 0; --j) {
            if (pop[n].MakeSpan + PrecisionValue < PFVec[j].MakeSpan ) { //测试用-xy
                cout.precision(2);
                cout<<fixed <<endl<<"pop: "<<pop[n].MakeSpan<< ", "<<pop[n].Cost<<";  PF: "<<PFVec[j].MakeSpan <<", "<<PFVec[n].Cost<<endl<<endl;
                for (int k = 0; k < PFVec.size(); ++k) {
                    cout << fixed << PFVec[k].MakeSpan<<", " << PFVec[k].Cost<<";  ";
                }
                cout<<endl<<endl;
                for(chromosome ch: pop ) {
                    cout << fixed << ch.MakeSpan << ", " << ch.Cost<<";  ";
                }
                cout<<endl<<endl;
                for(chromosome ch: orgPop ) {
                    cout<< fixed << ch.MakeSpan << ", " << ch.Cost<<";  ";
                }
                exit(5);
            }
            if ( PFVec[j].Cost - PrecisionValue < pop[n].Cost) {
                dominated_sign = false;
                break;
            }
        }
        if (dominated_sign) {
            PFVec.push_back(pop[n]);
            pop.erase(pop.begin() + n);
            --n;
        }
    }
    return PFVec;
}

void sortParetoFrontByCrowdingDistance (vector<chromosome>& PFVec){
    if (PFVec.size()> 2){
        /* 计算每个个体的拥挤度, 要求PFVec中的个体是按makespan从小到大排序的 */
        double MaxCst = PFVec[0].Cost, MinCst = PFVec[PFVec.size()-1].Cost;//1.1 取出最小的和最大的适应度值(成本)作为最小拥挤度和最大拥挤度
        double MaxMS = PFVec[PFVec.size()-1].MakeSpan, MinMS = PFVec[0].MakeSpan;//1.2 取出最小的和最大的适应度值(工期)作为最小拥挤度和最大拥挤度
        /* 1.3 把目标值最小的和最大拥挤距离设为无穷大，其他初始为0 */
        for(int i = 1; i < PFVec.size()-1; ++i)   PFVec[i].Distance = 0;
        PFVec[0].Distance = InfInt;
        PFVec[PFVec.size()-1].Distance = InfInt;
        /* 1.4 计算中间部分的拥挤距离 */
        for(int i = 1; i < (PFVec.size()-1); ++i){
            PFVec[i].Distance += (PFVec[i-1].Cost - PFVec[i+1].Cost) / (MaxCst - MinCst) + (PFVec[i+1].MakeSpan - PFVec[i-1].MakeSpan) / (MaxMS - MinMS);
        }
        sort(PFVec.begin(), PFVec.end(), sortPopOnDistanceByDesc);  //2.按拥挤距离从大到小排序
    }
}

double clcCstOfInsLeased (double strTimeOfUse, double endTimeOfUse, int VMType){
    double LeaseTime = XY_MAX(endTimeOfUse - strTimeOfUse + VMTypes[VMType].BootTime, VMTypes[VMType].MinLeaseTime);
    double cost =  VMTypes[VMType].FixedCost + VMTypes[VMType].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[VMType].ChargeUnit);
    return cost;
}
