//
// Created by xieyi on 2023/3/29.
//
#include "GenOper.h"
#include "ChrOper.h"
#include "tools.h"
//#include <list>
//using namespace std;

void initChrom(chromosome &chrom) {
    chrom.Code.clear();
    chrom.TskSchLst.resize(comConst.NumOfTasks, -1);
    chrom.InsAlcLst.resize(comConst.NumOfTasks, -1);
    chrom.TypeLst.resize(comConst.NumOfTasks, -1);
    chrom.TskSchPart.resize(comConst.NumOfTasks, -1);
    chrom.TypeAlcPart.resize(comConst.NumOfTasks, -1);
    chrom.VTskSchPart.resize(comConst.NumOfTasks,0);
    chrom.VTypeAlcPart.resize(comConst.NumOfTasks,0);
    chrom.StartTime.resize(comConst.NumOfTasks, -1);
    chrom.EndTime.resize(comConst.NumOfTasks, -1);
    chrom.FitnessValue= 0;
    chrom.MakeSpan = 0;
    chrom.Cost = 0;
    chrom.Distance = 0;
    chrom.Weight = 0;
    chrom.IsFrw = true;
    chrom.IsIFBDI = false;
    chrom.IsLBCAI = false;
}


void GnrRscAlc_Random(chromosome& chrom){ //InsAlcLst中生成的虚拟机实例编号是乱序的
    int b;
    do {
        set<int> avlTypes;
        for (int i = 0; i < comConst.NumOfTasks; ++i) {  //为每个实例随机生成实例的类型
            map<int, VMType>::iterator it = VMTypes.begin();
            int rnd = rand()%(comConst.NumOfVMTypes);
            advance(it, rnd);
            chrom.TypeLst[i] = (*it).first;
            avlTypes.insert((*it).first);
        }
        if (avlTypes.size() == comConst.NumOfVMTypes)  break;
        for(int i = 0; i < comConst.NumOfTasks; ++i){   //检查每个任务是否都存在可以处理它的虚拟机
            b = 0;
            for (int ty: Tasks[i].ElgRsc){
                if ( avlTypes.find(ty) != avlTypes.end()){
                    b = 1; break;
                }
            }
            if (b == 0) break;
        }
    } while (b == 0);

    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        vector<int> AvailInstList;
        for (int ty: Tasks[i].ElgRsc) {
            for (int k = 0; k < comConst.NumOfTasks; ++k) {
                if (ty == chrom.TypeLst[k]){
                    AvailInstList.push_back(k);
                }
            }
        }
        int rnd = rand()%(AvailInstList.size());
        chrom.InsAlcLst[i] = AvailInstList[rnd];
    }
}

void GnrRscAlc_BenefitRatio_MOELS(chromosome &chrom, double refCst, double refMs){
    //基于效益比生成个体,生成的个体中InsAlcLst[i]的实例编号是0开始的连续编号，但不一定满足InsAlcLst[i]<=max{InsAlcLst[i']}+1(i'<i)，另外租用时间只包含任务执行时间
    vector<set<double>> ITL;        //已经租用的实例的可得时间段列表
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        chrom.TypeLst[i] = -1;      //类型列表
        chrom.InsAlcLst[i] = -1;    //实例列表
    }
    chrom.Cost = 0; chrom.MakeSpan = 0; //初始化完工时间和成本
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        double FnlStartTime, FnlEndTime; //记录任务最终的开始和完成时间，即分配给具有最小效益比的虚拟机实例时的任务开始时间和完成时间；
        double MinEffRatio = InfInt;
        double TskStartTime, TskEndTime;
        double costAdded;
        int FnlVMid = -1, FnlTyId = -1;
        int TskId = chrom.TskSchLst[i];
        vector<int> AvailVMs;     //存放可以处理任务TskId的虚拟机实例编号
        for (int j = 0; j < ITL.size(); ++j) {
            if (Tasks[TskId].ElgRsc.find(chrom.TypeLst[j]) != Tasks[TskId].ElgRsc.end()) {
                AvailVMs.push_back(j);
            }
        }
        for(int vmId: AvailVMs){  //遍历已存在的虚拟机实例
            double ExeTime = Tasks[TskId].length / VMTypes[chrom.TypeLst[vmId]].ps;
            //计算任务TskId的就绪时间
            double rt = 0;
            for (int prnId: Tasks[TskId].parents) {
                int prnVMId = chrom.InsAlcLst[prnId];
                double fft = chrom.EndTime[prnId];
                if (vmId != prnVMId) {
                    fft += ParChildTranFileSizeSum[prnId][TskId] / XY_MIN(VMTypes[chrom.TypeLst[vmId]].bw, VMTypes[chrom.TypeLst[prnVMId]].bw);
                }
                if (rt + PrecisionValue < fft) {
                    rt = fft;
                }
            }
            TskStartTime = FindIdleTimeSlot(ITL[vmId], ExeTime, rt); //搜索可用时间片段获得任务开始时间
            TskEndTime = TskStartTime + ExeTime;
            double frtTime1, lstTime1, LeaseTime1;      //没有处理任务前虚拟机的开始使用时间，结束使用时间，租用时间归还时间
            double frtTime2, lstTime2, LeaseTime2;      //处理任务后虚拟机的开始使用时间，结束使用时间，租用时间归还时间
            if (ITL[vmId].size() > 2 || fabs( *(ITL[vmId].begin()) - VMTypes[chrom.TypeLst[vmId]].BootTime ) >= PrecisionValue) {
                set<double>::iterator pFrt = ITL[vmId].begin(), pLst = ITL[vmId].end();
                pLst--;  pLst--;
                frtTime1 = VMTypes[chrom.TypeLst[vmId]].BootTime, lstTime1 = *pLst;
                if (fabs(VMTypes[chrom.TypeLst[vmId]].BootTime - *pFrt) < PrecisionValue) {
                    pFrt++;
                    frtTime1 = *pFrt;
                }
            }else {
                cout << endl << "a VM being not used is added! this is wrong!"; exit(0);
            }
            //MOELS中是把任务分配给使部分工作流（所有已分配任务）的归一化时间成本最小的实例；
            // 部分工作流的归一化时间成本最小的实例 = 先前所有已分配任务的归一化时间成本 + 当前任务分配后增加的归一化时间成本；
            // 因此把任务分配给使部分工作流的归一化时间成本最小的实例，即等价于把任务分配给其分配后增加的归一化时间成本最小的实例；
            LeaseTime1 = XY_MAX(lstTime1 - frtTime1 + VMTypes[chrom.TypeLst[vmId]].BootTime, VMTypes[chrom.TypeLst[vmId]].MinLeaseTime);
            frtTime2 = XY_MIN(TskStartTime, frtTime1);
            lstTime2 = XY_MAX(TskEndTime, lstTime1);
            LeaseTime2 = XY_MAX(lstTime2 - frtTime2 + VMTypes[chrom.TypeLst[vmId]].BootTime, VMTypes[chrom.TypeLst[vmId]].MinLeaseTime);
            double temCstAdded = VMTypes[chrom.TypeLst[vmId]].UnitCost *
                                 (ceil((LeaseTime2-PrecisionValue)/VMTypes[chrom.TypeLst[vmId]].ChargeUnit)-ceil((LeaseTime1-PrecisionValue)/VMTypes[chrom.TypeLst[vmId]].ChargeUnit));
            double CrnEffRatio = (1-chrom.Weight) * temCstAdded / refCst + chrom.Weight * (XY_MAX(chrom.MakeSpan,TskEndTime) - chrom.MakeSpan) / refMs;
            if (CrnEffRatio + PrecisionValue < MinEffRatio) {
                FnlStartTime = TskStartTime;
                FnlEndTime = TskEndTime;
                FnlVMid = vmId;
                costAdded = temCstAdded;
                MinEffRatio = CrnEffRatio;
            }
        }
        for (int ty: Tasks[TskId].ElgRsc) {   //遍历可以使用的新增虚拟机实例
            double ExeTime = Tasks[TskId].length / VMTypes[ty].ps;    //执行时间
            double rt = 0;
            for (int prnId: Tasks[TskId].parents) {
                int prnVMId = chrom.InsAlcLst[prnId];
                double fft = chrom.EndTime[prnId] + ParChildTranFileSizeSum[prnId][TskId]/XY_MIN(VMTypes[ty].bw, VMTypes[chrom.TypeLst[prnVMId]].bw);
                if (rt + PrecisionValue < fft) {
                    rt = fft;
                }
            }
            TskStartTime = XY_MAX(rt,VMTypes[ty].BootTime);
            TskEndTime = TskStartTime + ExeTime; //开始时间(启动时间和就绪时间取最大)加上执行时间
            double LT = XY_MAX(ExeTime+VMTypes[ty].BootTime, VMTypes[ty].MinLeaseTime);
            double temCstAdded = VMTypes[ty].FixedCost + VMTypes[ty].UnitCost * ceil((LT-PrecisionValue)/VMTypes[ty].ChargeUnit);
            double CrnEffRatio = (1-chrom.Weight) * temCstAdded / refCst + chrom.Weight * (XY_MAX(chrom.MakeSpan, TskEndTime) - chrom.MakeSpan) / refMs;
            if (CrnEffRatio + PrecisionValue < MinEffRatio) {
                FnlStartTime = TskStartTime;
                FnlEndTime = TskEndTime;
                FnlTyId = ty;
                costAdded = temCstAdded;
                MinEffRatio = CrnEffRatio;
            }
        }
        if (FnlTyId != -1) {   //表明新增一个类型为FnlTyId的虚拟机实例时综合效益比是最好的
            FnlVMid = ITL.size();
            chrom.TypeLst[FnlVMid] = FnlTyId;
            set<double> a;
            a.insert(VMTypes[FnlTyId].BootTime);  //虚拟机实例初始化
            a.insert(InfInt * 1.0);
            ITL.push_back(a);
        }
        chrom.StartTime[TskId] = FnlStartTime;
        chrom.EndTime[TskId] = FnlEndTime;
        chrom.InsAlcLst[TskId] = FnlVMid;
        chrom.MakeSpan = XY_MAX(chrom.MakeSpan, FnlEndTime);
        UpdateITL(ITL[FnlVMid], FnlStartTime, FnlEndTime);
        chrom.Cost += costAdded;  //根据增量计算成本
    }
}

//void GnrRscAlc_BenefitRatio(chromosome &chrom){       //基于效益比 桂奉献实现的方法
//    vector<set<double> > vatL;
//    double max_EndTi = -1000000;
//    int INS = 0;    //因为INS里面存储的值与下标相同，把vector改为int类型来计算
//
//    for (int m = 0; m < comConst.NumOfTasks; ++m) {
//        chrom.TypeLst[m] = -1;      //类型列表
//        chrom.InsAlcLst[m] = -1;    //实例列表
//    }
//
//    for(int i = 0; i < comConst.NumOfTasks; ++i){
//        int vmIndex = -1;
//        int type = -1;
//        vector<double> FinalendTime;
//        vector<double>  FinalstartTime;
//        double TaskIndexStarTime;
//        double TaskIndexEndTime;
//        int TaskIndex = chrom.TskSchLst[i];
//        vector<int> AvailVMIndex;                           //在INS中 && 可获得
//        for (int i = 0; i < INS; ++i) {
//            for (int i1 = 0; i1 < Tasks[TaskIndex].ElgRsc.size(); ++i1) {
//                if (chrom.TypeLst[i] == Tasks[TaskIndex].ElgRsc[i1]) {
//                    AvailVMIndex.push_back(i);
//                    break;
//                }
//            }
//        }
//        //找已执行过的任务的最大完成时间
//        double FiMax = 0; //记录完成时间(已执行过的任务)
//        for (int k = 0; k < i; k++) {
//            double EndTime = chrom.EndTime[chrom.TskSchLst[k]];
//            if (FiMax < EndTime) {
//                FiMax = EndTime;
//            }
//        }
//        //为了计算虚拟机归还时间，找任务的最大输出，传输大小只算此任务所在虚拟机即可
//        double MaxTranChildren = 0;
//        for(int cl = 0; cl < Tasks[TaskIndex].children.size(); ++cl){
//            double Tran = ParChildTranFileSizeSum[TaskIndex][Tasks[TaskIndex].children[cl]];
//            if(MaxTranChildren < Tran)
//                MaxTranChildren = Tran;
//        }
//        vector<double > ComprehensiveRatio;
//        //合并到以前的虚拟机里
//        for(int j = 0; j < AvailVMIndex.size(); ++j){       //任务taskIndex可获得且在INS集合中的虚拟机实例的集合
//            int v = AvailVMIndex[j];
//            double DurTime = Tasks[TaskIndex].length / VMTypes[chrom.TypeLst[v]].ps;
//
//            //计算任务taskIndex就绪时间
//            double rt = 0;
//            for (int i = 0; i < Tasks[TaskIndex].parents.size(); i++) {
//                int parentIndex = Tasks[TaskIndex].parents[i];
//                int parentVMIndex = chrom.InsAlcLst[parentIndex];
//                double ParEndTime = chrom.EndTime[parentIndex];
//                if (v != parentVMIndex) {
//                    double sum = ParChildTranFileSizeSum[parentIndex][TaskIndex];
//                    ParEndTime += sum / (XY_MIN(VMTypes[chrom.TypeLst[v]].bw, VMTypes[chrom.TypeLst[parentVMIndex]].bw));
//                }
//                if (rt < ParEndTime) {
//                    rt = ParEndTime;
//                }
//            }
//
//            //搜索可用时间片段
//            set<double>::iterator pre  = vatL[v].begin();
//            set<double>::iterator post = vatL[v].begin();
//            ++post;
//            while(post != vatL[v].end()) {
//                if((*post - *pre) >= DurTime && rt <= (*post)-DurTime) {
//                    TaskIndexStarTime = XY_MAX(*pre, rt);
//                    break;
//                } else {
//                    ++pre;
//                    ++pre;
//                    ++post;
//                    ++post;
//                }
//            }
//            TaskIndexEndTime = TaskIndexStarTime + DurTime;
//            FinalstartTime.push_back(TaskIndexStarTime);
//            FinalendTime.push_back(TaskIndexEndTime);
//
//            double Rnt = 0;                                 //虚拟机的归还时间
//            double Hrt = 100000000000;                      //虚拟机的开始使用时间-还没有分配编号为k的资源实例的开始租用时间
//
//            vector<int> handleTask;
//            for (int n = 0; n < i; ++n) {                   //找相同实例的任务，放在同一个容器handleTask
//                if (chrom.InsAlcLst[chrom.TskSchLst[n]] == v)
//                    handleTask.push_back(chrom.TskSchLst[n]);
//            }
//            for (int k = 0; k < handleTask.size(); ++k) {   //Rnt的计算 归还时间
//                double time =  chrom.EndTime[handleTask[k]];
//                if (Rnt < time) Rnt = time;
//            }
//
//            for (int l = 0; l < handleTask.size(); ++l) {   //Hrt的计算 开始租用时间
//                double hrt = chrom.StartTime[handleTask[l]] - VMTypes[chrom.TypeLst[v]].BootTime;
//                if (hrt < Hrt) Hrt = hrt;
//            }
//
//            double Lt = Rnt + MaxTranChildren/VMTypes[v].bw - Hrt;  //表示该任务还没有分配给虚拟机的租用时间
//            double Lt1 = XY_MAX(TaskIndexEndTime, Rnt) - XY_MIN(TaskIndexStarTime - VMTypes[chrom.TypeLst[v]].BootTime, Hrt);
//
//            double comprehensiveratio = (inputParameter.theta_eff * VMTypes[chrom.TypeLst[v]].UnitCost *
//                                         ((ceil(XY_MAX(Lt1, VMTypes[chrom.TypeLst[v]].MinLeaseTime) /
//                                                VMTypes[chrom.TypeLst[v]].ChargeUnit)) -
//                                          (ceil(XY_MAX(Lt, VMTypes[chrom.TypeLst[v]].MinLeaseTime) /
//                                                VMTypes[chrom.TypeLst[v]].ChargeUnit)))) / Budget
//                                        + (1 - inputParameter.theta_eff) *
//                                          (XY_MAX(FiMax, TaskIndexEndTime) - FiMax) / Deadline;
//            ComprehensiveRatio.push_back(comprehensiveratio);
//        }
//
//        for (int j = 0; j < Tasks[TaskIndex].ElgRsc.size(); ++j) {  //任务taskIndex可得虚拟机类型的循环
//            int beta= INS;
//            int RecType = Tasks[TaskIndex].ElgRsc[j];
//            double DurTime = Tasks[TaskIndex].length / VMTypes[RecType].ps;    //执行时间
//            double rt = 0;
//            for (int i = 0;  i < Tasks[TaskIndex].parents.size();  ++i) {
//                int parentIndex = Tasks[TaskIndex].parents[i];
//                int parentVMIndex = chrom.InsAlcLst[parentIndex];
//                double Max = chrom.EndTime[parentIndex];
//                Max += ParChildTranFileSizeSum[parentIndex][TaskIndex] / (XY_MIN(VMTypes[RecType].bw, VMTypes[chrom.TypeLst[parentVMIndex]].bw));
//                if (rt < Max) {
//                    rt = Max;
//                }
//            }
//            TaskIndexStarTime = rt;
//            TaskIndexEndTime = XY_MAX(rt,VMTypes[RecType].BootTime) + DurTime; //开始时间(启动时间和就绪时间取最大)加上执行时间
//            double Lt = TaskIndexEndTime + MaxTranChildren/VMTypes[RecType].bw;
//
//            double comprehensiveratio = (inputParameter.theta_eff * (VMTypes[RecType].FixedCost +VMTypes[RecType].UnitCost
//                                                                                                 * (ceil(XY_MAX(VMTypes[RecType].MinLeaseTime,Lt)/
//                                                                                                         VMTypes[RecType].ChargeUnit)))) / Budget  + (1 - inputParameter.theta_eff) *
//                                                                                                                                                     (XY_MAX(FiMax, TaskIndexEndTime) - FiMax) / Deadline;
//
//            FinalstartTime.push_back(TaskIndexStarTime);
//            FinalendTime.push_back(TaskIndexEndTime);
//            ComprehensiveRatio.push_back(comprehensiveratio);
//        }
//        vector<int> ind(ComprehensiveRatio.size());
//        IndexSortAsc(ind, ComprehensiveRatio);               //对综合效益比排序，ind[0]存放的是综合效益比最小的下标
//
//        /**
//         * 任务先从INS&&可获得 里先找，这时下标是小于AvailVMIndex.size()，如果此下标是效益比最小的，那么不用添加新的实例
//         * 如果此下标大于AvailVMIndex.size()且是效益比最小的，那么添加的是新的实例
//         */
//        if(ind[0] >= AvailVMIndex.size()){                              //新添加一个实例   //ind[0]就是要添加的那个实例(综合效益比最小)
//            vmIndex = INS;                                       //ind[0]如果大于AvailVMIndex.size()即表示新增的vmIndex
//            type = Tasks[TaskIndex].ElgRsc[ind[0]-AvailVMIndex.size()]; //记录添加的类型
//            //在此初始化一个虚拟机实例
//            /**
//             * 每当INS新增一个虚拟机实例，vatL就会加一维，所以INS和vatL的下标是同步变化的，一一对应
//             */
//            INS++;                                                      //新增了一个实例，INS值加1
//            set<double> a;
//            a.insert(VMTypes[type].BootTime);                           //虚拟机初始化
//            a.insert(InfInt * 1.0);
//            vatL.push_back(a);
//            chrom.TypeLst[vmIndex]=type;                                //当新添加实例时，才会给此实例赋于类型值
//        }
//        if(ind[0]<AvailVMIndex.size()){                                 //没有添加新实例
//            vmIndex = AvailVMIndex[ind[0]];                             //找出vmIndex下面要用到
//        }
//        chrom.StartTime[TaskIndex] = FinalstartTime[ind[0]];
//        chrom.EndTime[TaskIndex] = FinalendTime[ind[0]];
//        chrom.InsAlcLst[TaskIndex] = vmIndex;
//
//        max_EndTi = XY_MAX(max_EndTi, chrom.EndTime[TaskIndex]);
//        if(vatL[vmIndex].find(chrom.StartTime[TaskIndex]) != vatL[vmIndex].end()) {
//            vatL[vmIndex].erase(chrom.StartTime[TaskIndex]);
//        } else {
//            vatL[vmIndex].insert(chrom.StartTime[TaskIndex]);
//        }
//
//        if(vatL[vmIndex].find(chrom.EndTime[TaskIndex]) != vatL[vmIndex].end()) {
//            vatL[vmIndex].erase(chrom.EndTime[TaskIndex]);
//        } else {
//            vatL[vmIndex].insert(chrom.EndTime[TaskIndex]);
//        }
//    }
//
//    //一维变二维
//    GnrTwoDimCoding(chrom);
//    //工期
//    chrom.MakeSpan = max_EndTi;
//    //计算成本
//    double cost =0;
//    list<int> List;
//    set<list<int> >::iterator st;
//    list<int>::iterator lst;
//    for(int k = 0; k < chrom.Code.size(); ++k){
//        st = chrom.Code.begin();
//        advance(st, k);             //使指针st指向Code中第k个list
//        List = *st;                    //把Code中第k个list赋值给此新的集合list
//        //虚拟机开始租用时间的计算（Hrt）
//        lst = List.begin();   lst++;                //该虚拟机处理的第一个任务
//
//        double Hrt = chrom.StartTime[*lst] - VMTypes[*List.begin()].BootTime;
//
//        //虚拟机归还时间计算（Rnt）
//        lst = List.end(); lst--;        //该虚拟机处理的最后一个任务
//        double Rnt = chrom.EndTime[*lst];
//
//        double Lt = Rnt - Hrt;          //最晚归还时间减去最早开始时间为虚拟机的租用时间
//        double spantime = XY_MAX(Lt, VMTypes[*List.begin()].MinLeaseTime);
//        cost += VMTypes[*List.begin()].FixedCost + (ceil((spantime) / VMTypes[*List.begin()].ChargeUnit)) * VMTypes[*List.begin()].UnitCost;
//    }
//    chrom.Cost = cost;
//}

//void GnrRscAlc_BenefitRatioNew(chromosome &chrom){
//    vector<set<double>> vatL;
//    int NumOfUsedIns = 0;       //已使用的实例数量，同时也用来表示对应的实例编号
//    for (int m = 0; m < comConst.NumOfTasks; ++m) {
//        chrom.TypeLst[m] = -1;      //类型列表
//        chrom.InsAlcLst[m] = -1;    //实例列表
//    }
//    for (int i = 0; i < comConst.NumOfTasks; ++i) {
//        int InsId = -1;
//        vector<double> AllPsbEndTime;
//        vector<double>  AllPsbStartTime;
//        double TaskStartTime;
//        double TaskEndTime;
//        int TaskId = chrom.TskSchLst[i];
//        vector<int> ElgInsId;
//        for (int i = 0; i < NumOfUsedIns; ++i) {
//            for (int ty: Tasks[TaskId].ElgRsc) {
//                if (chrom.TypeLst[i] == ty) {
//                    ElgInsId.push_back(i);
//                    break;
//                }
//            }
//        }
//
//        vector<pair<double,double>>InsLTAndAddCost; //可能分配的实例租用时间集合和可能增加的成本集合
//        vector<double > AllComprehensiveRatio;
//        vector<double>CurLT(comConst.NumOfTasks,0); // 记录各个实例的当前租用时间，初始化为0
//        for(int j = 0; j < ElgInsId.size(); ++j){
//            int v = ElgInsId[j];
//            double DurTime = Tasks[TaskId].length / VMTypes[chrom.TypeLst[v]].ps;
//            double rt = 0;
//            for (int i = 0; i < Tasks[TaskId].parents.size(); i++) {
//                int ParentId = Tasks[TaskId].parents[i];
//                int InsIdOfParent = chrom.InsAlcLst[ParentId];
//                double ParEndTime = chrom.EndTime[ParentId];
//                if (v != InsIdOfParent) {
//                    ParEndTime += ParChildTranFileSizeSum[InsIdOfParent][TaskId] / (XY_MIN(VMTypes[chrom.TypeLst[v]].bw, VMTypes[chrom.TypeLst[InsIdOfParent]].bw));
//                }
//                if (rt < ParEndTime) {
//                    rt = ParEndTime;
//                }
//            }
//            //搜索可用时间片段
//            set<double>::iterator pre  = vatL[v].begin();
//            set<double>::iterator post = vatL[v].begin();
//            ++post;
//            while(post != vatL[v].end()) {
//                if((*post - *pre) >= DurTime && rt <= (*post)-DurTime) {
//                    TaskStartTime = XY_MAX(*pre, rt);
//                    break;
//                } else {
//                    ++pre;
//                    ++pre;
//                    ++post;
//                    ++post;
//                }
//            }
//            TaskEndTime = TaskStartTime + DurTime;
//            AllPsbStartTime.push_back(TaskStartTime);
//            AllPsbEndTime.push_back(TaskEndTime);
//            double TemMaxMakeSpan = XY_MAX(chrom.MakeSpan, TaskEndTime);
//            double Rnt = 0;                                 //虚拟机的归还时间
//            double Hrt = 100000000000;                      //虚拟机的开始使用时间-还没有分配编号为k的资源实例的开始租用时间
//            vector<int> handleTask;
//            double TemMaxCost;
//            for (int n = 0; n < i; ++n) {
//                if (chrom.InsAlcLst[chrom.TskSchLst[n]] == v)
//                    handleTask.push_back(chrom.TskSchLst[n]);
//            }
//            for (int k = 0; k < handleTask.size(); ++k) {   //Rnt的计算 归还时间
//                double time =  chrom.EndTime[handleTask[k]];
//                if (Rnt < time) Rnt = time;
//            }
//            Rnt = XY_MAX(Rnt, TaskEndTime);                  //此正在计算还没有加入染色体的任务也要算上
//            for (int l = 0; l < handleTask.size(); ++l) {   //Hrt的计算 开始租用时间
//                double hrt = chrom.StartTime[handleTask[l]] - VMTypes[chrom.TypeLst[v]].BootTime;
//                if (hrt < Hrt) Hrt = hrt;
//            }
//            Hrt = XY_MIN(Hrt, TaskStartTime - VMTypes[chrom.TypeLst[v]].BootTime);
//
//            double Lt = Rnt - Hrt;         //当前任务执行后虚拟机租用时间
//            double AddMaxCost =  (ceil((XY_MAX(Lt,VMTypes[chrom.TypeLst[v]].MinLeaseTime))/VMTypes[chrom.TypeLst[v]].ChargeUnit) - (ceil(XY_MAX(CurLT[NumOfUsedIns],VMTypes[chrom.TypeLst[v]].MinLeaseTime))/VMTypes[chrom.TypeLst[v]].ChargeUnit)) *VMTypes[chrom.TypeLst[v]].UnitCost;
//            InsLTAndAddCost.emplace_back(Lt,AddMaxCost);
//            TemMaxCost = chrom.Cost + AddMaxCost;
//            double Comprehensiveratio =  (chrom.Weight * TemMaxMakeSpan) / Deadline + (1 - chrom.Weight) * TemMaxCost / Budget;
//            AllComprehensiveRatio.push_back(Comprehensiveratio);
//        }
//
//        for (int RecType: Tasks[TaskId].ElgRsc){
//            double DurTime = Tasks[TaskId].length / VMTypes[RecType].ps;    //执行时间
//            double rt = 0;
//            for (int i = 0; i < Tasks[TaskId].parents.size(); i++) {
//                int ParentId = Tasks[TaskId].parents[i];
//                int InsIdOfParent = chrom.InsAlcLst[ParentId];
//                double Max = chrom.EndTime[ParentId];
//                double TemInputTime = ParChildTranFileSizeSum[ParentId][TaskId] / (XY_MIN(VMTypes[RecType].bw, VMTypes[chrom.TypeLst[InsIdOfParent]].bw));
//                if (rt < Max) {
//                    rt = Max;
//                }
//            }
//            TaskStartTime = XY_MAX(rt,VMTypes[RecType].BootTime);
//            TaskEndTime = TaskStartTime + DurTime;
//
//            double TemMaxMakeSpan = XY_MAX(chrom.MakeSpan, TaskEndTime);
//            double Hrt = TaskStartTime - VMTypes[RecType].BootTime;
//            double Rnt = TaskEndTime;
//            double Lt = Rnt - Hrt;
//            double AddMaxCost = ceil((XY_MAX(Lt,VMTypes[RecType].MinLeaseTime))/VMTypes[RecType].ChargeUnit) *VMTypes[RecType].UnitCost + VMTypes[RecType].FixedCost;
//            InsLTAndAddCost.emplace_back(Lt,AddMaxCost);
//            double TemMaxCost = chrom.Cost + AddMaxCost;
//            double Comprehensiveratio =  (chrom.Weight * TemMaxMakeSpan) / Deadline + (1 - chrom.Weight) * TemMaxCost / Budget;
//            AllComprehensiveRatio.push_back(Comprehensiveratio);
//            AllPsbStartTime.push_back(TaskStartTime);
//            AllPsbEndTime.push_back(TaskEndTime);
//        }
//        vector<int> ind(AllComprehensiveRatio.size());
//        IndexSortByValueAsc(ind, AllComprehensiveRatio);               //对综合效益比排序，ind[0]存放的是综合效益比最小的下标
//        int type = -1;
//        if(ind[0] >= ElgInsId.size()){                              //新添加一个实例   //ind[0]就是要添加的那个实例(综合效益比最小)
//            InsId = NumOfUsedIns;                                       //ind[0]如果大于AvailVMIndex.size()即表示新增的vmIndex
//            set<int>::iterator it = Tasks[TaskId].ElgRsc.begin();
//            advance(it,ind[0]-ElgInsId.size());
//            type = *it; //记录添加的类型
//            //在此初始化一个虚拟机实例
//            /**
//             * 每当INS新增一个虚拟机实例，vatL就会加一维，所以INS和vatL的下标是同步变化的，一一对应
//             */
//            NumOfUsedIns++;                                                      //新增了一个实例，INS值加1
//            set<double> a;
//            a.insert(VMTypes[type].BootTime);                           //虚拟机初始化
//            a.insert(InfInt * 1.0);
//            vatL.push_back(a);
//            chrom.TypeLst[InsId]=type;                                //当新添加实例时，才会给此实例赋于类型值
//        }
//        if(ind[0]<ElgInsId.size()){                                 //没有添加新实例
//            InsId = ElgInsId[ind[0]];                             //找出vmIndex下面要用到
//        }
//        CurLT[InsId] +=  InsLTAndAddCost[ind[0]].first; // 当前实例的租用时间
//        chrom.StartTime[TaskId] = AllPsbStartTime[ind[0]];
//        chrom.EndTime[TaskId] = AllPsbEndTime[ind[0]];
//        chrom.InsAlcLst[TaskId] = InsId;
//
//        chrom.MakeSpan = XY_MAX(chrom.MakeSpan, chrom.EndTime[TaskId]);
//        chrom.Cost += InsLTAndAddCost[ind[0]].second;  // 当前已租用的成本
//
//        if(vatL[InsId].find(chrom.StartTime[TaskId]) != vatL[InsId].end()) {
//            vatL[InsId].erase(chrom.StartTime[TaskId]);
//        } else {
//            vatL[InsId].insert(chrom.StartTime[TaskId]);
//        }
//
//        if(vatL[InsId].find(chrom.EndTime[TaskId]) != vatL[InsId].end()) {
//            vatL[InsId].erase(chrom.EndTime[TaskId]);
//        } else {
//            vatL[InsId].insert(chrom.EndTime[TaskId]);
//        }
//    }
////    //测试用-xy
////    chromosome temChr = chrom;
////    DecodeChrom(temChr, true);
////    if ( fabs(temChr.MakeSpan - chrom.MakeSpan) > PrecisionValue|| fabs(temChr.Cost - chrom.Cost) > PrecisionValue ) {
////        cout << endl << "Original: "<< chrom.MakeSpan << ","<<chrom.Cost<<"; Decoded:"<< temChr.MakeSpan << ","<<temChr.Cost;
////        exit(5);
////    }
//    chrom.FitnessValue = (1 - (chrom.MakeSpan - Deadline) / (Max_MakeSpan - Deadline)) * (1 - (chrom.Cost - Budget) / (Max_Cost - Budget));
//}

/*计算平均任务执行时间*/
void ClcAvrTskExcTime(vector<double>& w) {
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        double sum = 0;
        for (int ty: Tasks[i].ElgRsc )
            sum += 1.0 / VMTypes[ty].ps;
        w[i] = Tasks[i].length * sum / Tasks[i].ElgRsc.size();
    }
}

/*使用平均带宽计算两个任务之间文件的平均传输时间*/
void ClcAvrCmnTime_AB(vector<vector<double>>& c) {
    //使用平均带宽计算平均数据传输时间-xy
    //MOELS、HPSO算法中虚拟机之间的带宽认为是相同的
    double aveBw = 0;
    for(auto elm : VMTypes){
        aveBw += elm.second.bw;
    }
    aveBw =  aveBw / VMTypes.size();
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        for (int k = 0; k < comConst.NumOfTasks; ++k) {
            c[i][k] = ParChildTranFileSizeSum[i][k] / aveBw;
        }
    }
}

/*考虑无需传输文件及有资格处理任务的虚拟机类型可能不同的情况下计算两个任务之间文件的平均传输时间*/
void ClcAvrCmnTime(vector<vector<double>>& c) { //重写-xy
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        for (int parent: Tasks[i].parents) {
            double sum1 = 0;
            for (int ty2: Tasks[i].ElgRsc) {
                for (int ty1: Tasks[parent].ElgRsc) {
                    if (ty1 == ty2) {
                        continue;
                    } else {
                        sum1 += ParChildTranFileSizeSum[parent][i] / XY_MIN(VMTypes[ty1].bw, VMTypes[ty2].bw);
                    }
                }
            }
            c[parent][i] = sum1 / (double) (Tasks[i].ElgRsc.size() * Tasks[parent].ElgRsc.size());
        }
    }
}

void Calculate_Rank_t(vector<double>& RankList, vector<double>& w, vector<vector<double>>& c) {
    for(int i = 1 ;i < TaskListInLevel.size(); ++i){
        for (int TskId: TaskListInLevel[i]) {
            for (int prnTskId: Tasks[TskId].parents) {
                double rnk = w[prnTskId] + c[prnTskId][TskId] + RankList[prnTskId];
                if (RankList[TskId] < rnk) {
                    RankList[TskId] = rnk;
                }
            }
        }
    }
}

//Rank_b，独立I/O流系统；
void Calculate_Rank_b(vector<double>& RankList, vector<double>& w, vector<vector<double>>& c){
    for (int taskId: TaskListInLevel[TaskListInLevel.size()-1]) {
        RankList[taskId] = w[taskId];
    }
    for(int i = TaskListInLevel.size()-2; i >=0 ;--i){
        for (int taskId: TaskListInLevel[i]) {
            double ChildMaxRank = 0;
            for (int childId: Tasks[taskId].children) {
                if(ChildMaxRank < RankList[childId] + c[taskId][childId] ){
                    ChildMaxRank = RankList[childId] + c[taskId][childId];
                }
            }
            RankList[taskId] = w[taskId] + ChildMaxRank;
        }
    }
}

void Calculate_Rank_b_GRPHEFT(vector<double>& RankList, vector<double>& w){
    //GRPHEFT的原始程序中把带宽设置为固定值39321600
    for (int i = 0; i < TaskListInLevel[TaskListInLevel.size()-1].size(); ++i) {
        int taskId=TaskListInLevel[TaskListInLevel.size()-1][i];
        RankList[taskId] = w[taskId];
    }
    for(int i =TaskListInLevel.size()-2 ;i >=0 ;--i){
        for (int j = 0; j < TaskListInLevel[i].size(); ++j) {
            int tsk = TaskListInLevel[i][j];
            double childMaxRankc = 0;
            for (int chlTsk: Tasks[tsk].children) {
                double CompareObject = RankList[chlTsk] + ParChildTranFileSizeSum[tsk][chlTsk] / 39321600;
                if(childMaxRankc < CompareObject ){
                    childMaxRankc = CompareObject;
                }
            }
            RankList[tsk] = w[tsk] + childMaxRankc;    //任务平均执行时间加上所有子任务中最大那个rank值
        }
    }
}

void Calculate_Rank_b(vector<double>& RankList, vector<double>& w){
    double aveBw = 0;
    for(auto elm : VMTypes) aveBw += elm.second.bw;
    aveBw =  aveBw / VMTypes.size();
    for (int i = 0; i < TaskListInLevel[TaskListInLevel.size()-1].size(); ++i) {
        int taskId=TaskListInLevel[TaskListInLevel.size()-1][i];
        RankList[taskId] = w[taskId];
    }
    for(int i =TaskListInLevel.size()-2 ;i >=0 ;--i){
        for (int j = 0; j < TaskListInLevel[i].size(); ++j) {
            int tsk = TaskListInLevel[i][j];
            double childMaxRankc = 0;
            for (int chlTsk: Tasks[tsk].children) {
                double CompareObject = RankList[chlTsk] + ParChildTranFileSizeSum[tsk][chlTsk] / aveBw;
                if(childMaxRankc < CompareObject ){
                    childMaxRankc = CompareObject;
                }
            }
            RankList[tsk] = w[tsk] + childMaxRankc;    //任务平均执行时间加上所有子任务中最大那个rank值
        }
    }
}

vector<chromosome> GnrChrsByGRPHEFT(double estMinBdg, double estMaxBdg, vector<int>& TskSchLstByRnk_b){
    //在选取生成的虚拟机实例集合Y时没有考虑实例对任务的支持情况，但在进行实例分配时考虑了实例对任务的支持情况，
    //因此本算法适用于资源可得率小于1的情况，但可能会出现很多情况下找不到解（return 0）
    //此外在选取生成的虚拟机实例集合Y时也没有考虑固定起租成本和最小起租时间，算法虽然能用于考虑固定起租成本和最小起租时间，但效果可能会变差！
    vector<VMType> TypeListBySrt; //所有相互非支配的虚拟机类型集合按效益（一般来说ps越大效益越好）从高到低排序
    for (auto vmty: VMTypes) {
        TypeListBySrt.push_back(vmty.second);
    }
    //按处理能力对虚拟机类型从大到小进行排序, 如果处理能力一样大则成本小的排在前面！
    for (int i = 0; i < TypeListBySrt.size(); ++i) {
        for (int j = 1; j < TypeListBySrt.size()-i; ++j) {
            if ( (TypeListBySrt[j-1].ps + PrecisionValue < TypeListBySrt[j].ps) ||
                 (fabs(TypeListBySrt[j-1].ps - TypeListBySrt[j].ps) < PrecisionValue) && TypeListBySrt[j-1].UnitCost - PrecisionValue > TypeListBySrt[j].UnitCost ) {
                VMType tem = TypeListBySrt[j-1];
                TypeListBySrt[j-1] = TypeListBySrt[j]; TypeListBySrt[j] = tem;
            }
        }
    }
    //去除所有支配的虚拟机类型形成相互间非支配的虚拟机类型列表
    int crnPt = 1;   //当前要更新数据的位置；
    int cmpPt = 1;   //比较的开始位置
    for ( ; crnPt < TypeListBySrt.size(); ++crnPt) {
        int k = cmpPt;
        for (; k < TypeListBySrt.size(); ++k ) {
            if (TypeListBySrt[k].UnitCost + PrecisionValue < TypeListBySrt[crnPt-1].UnitCost) { //找到第一个成本小于当前考察的（即i-1）
                TypeListBySrt[crnPt] = TypeListBySrt[k];
                cmpPt = k + 1;  //记录下次比较的起始位置
                break;
            }
        }
        if (k==TypeListBySrt.size()) break;
    }
    TypeListBySrt.erase(TypeListBySrt.begin() + crnPt,TypeListBySrt.end());
    //对相互间非支配的虚拟机类型按效益比从大到小进行排序
    for (int i = 0; i < TypeListBySrt.size(); ++i) {
        for (int j = 1; j < TypeListBySrt.size()-i; ++j) {
            if (TypeListBySrt[j-1].ps/TypeListBySrt[j-1].UnitCost + PrecisionValue < TypeListBySrt[j].ps/TypeListBySrt[j].UnitCost) {
                VMType tem = TypeListBySrt[j-1];
                TypeListBySrt[j-1] = TypeListBySrt[j]; TypeListBySrt[j] = tem;
            }
        }
    }
    double MinUintCost = TypeListBySrt[0].UnitCost;
    for (int j = 1; j < TypeListBySrt.size(); ++j) {
        if (TypeListBySrt[j].UnitCost + PrecisionValue < MinUintCost)
            MinUintCost = TypeListBySrt[j].UnitCost;
    }
    if (estMinBdg + PrecisionValue < MinUintCost ) {
        cout << endl << "The budget is less than the minimum unit cost of VM in GrnChromsByGRPHEFT, and No feasible solution!";
        exit(2);
    }
    set<chromosome> result;
    double stepValue = 1.0 * (estMaxBdg - estMinBdg) / comConst.NumOfTasks;
    double currValue = estMinBdg;
    while (currValue - PrecisionValue < estMaxBdg){
        chromosome BestChrom;
        initChrom(BestChrom);
        BestChrom.MakeSpan = InfInt; BestChrom.Cost = InfInt;
        vector<Instance> LstY;
        //通过while循环确定最终方案
        for (int m = 0; m < TypeListBySrt.size(); ++m) {
            double Remaining_Budget = currValue;
            chromosome chrom;
            initChrom(chrom);
            chrom.TskSchLst = TskSchLstByRnk_b; //生成任务调度顺序
            vector<Instance> Y;          //选取的虚拟机实例集合
            map<int,int> NumOfInsUnused; //记录每种类型的实例的未使用数量
            for (int j = 0; j < TypeListBySrt.size(); ++j) {
                NumOfInsUnused[TypeListBySrt[j].type] = - 9999999;
            }
            int InsId = 0, j = m;
            double temMinUintCost = TypeListBySrt[j].UnitCost;
            for (int k = j+1; k < TypeListBySrt.size(); ++k) {
                if (TypeListBySrt[k].UnitCost + PrecisionValue < temMinUintCost)
                    temMinUintCost = TypeListBySrt[k].UnitCost;
            }
            //{根据预算选取虚拟机建立Y集合}，在确保每个任务都有虚拟处理的情况下，考虑按效率比生成Y集合-xy
            while (Remaining_Budget + PrecisionValue > temMinUintCost) {
                //存储的顺序是ps从小到大，所以从后向前选，-应该改成按效率比而不是ps！ 原程序也是这样吗？-xy; 存储的顺序已改为根据效率比从大到小-xy；
                NumOfInsUnused[TypeListBySrt[j].type] = floor((Remaining_Budget + PrecisionValue) / TypeListBySrt[j].UnitCost); //不考虑最小起租时间，或UnitCost大于等于最小起租时间
                //按此生成实例集合Y的方法必定满足以下条件：任一实例的最小租用单位成本必定大于所有其后不同于该类型实例的最小租用单位成本成本之和；
                //这也意味着：当前实例要增加一个租用时间单元而又不存在同类型的还未使用的补偿实例时，且其前面的所有实例都已使用时，该方案必为不可行方案！
                for (int k = 0; k < NumOfInsUnused[TypeListBySrt[j].type]; k++){
                    Instance TemInst;
                    TemInst.id = InsId; ++InsId;
                    TemInst.type = TypeListBySrt[j].type;
                    TemInst.ITL.insert(TypeListBySrt[j].BootTime); TemInst.ITL.insert(InfInt * 1.0);
                    TemInst.LET = 0, TemInst.LST = InfInt;
                    TemInst.state = 0;
                    Y.push_back(TemInst);
                }
                Remaining_Budget -= TypeListBySrt[j].UnitCost * NumOfInsUnused[TypeListBySrt[j].type];
                ++j;
            }
            //方法2：按完成时间从小到大分配给不增加成本或存在补偿实例的实例，否则返回无可行调度方案
            if (Modified_HEFT_XY_Imp(chrom, Y, NumOfInsUnused) == 0) {
                continue;
            }

            if(chrom.Cost - PrecisionValue < currValue && chrom.MakeSpan - PrecisionValue < BestChrom.MakeSpan){
                BestChrom = chrom; LstY = Y;
            }
        }
        if (BestChrom.MakeSpan + PrecisionValue < InfInt)
            result.insert(BestChrom);
        currValue += stepValue;
    }
    vector<chromosome> ParetoFront = GnrParetoFront(result);
    return ParetoFront;
}


//void GnrTwoDimCoding_gfx(chromosome &chrom){ //实例编号可以不连续
//    set<int> vmSet;
//    for(int i = 0; i < chrom.InsAlcLst.size(); ++i){
//        vmSet.insert(chrom.InsAlcLst[i]);
//    }
//    for( set<int>::iterator st = vmSet.begin(); st != vmSet.end(); ++st){
//        list<int> GL;
//        GL.insert(GL.begin(), chrom.TypeLst[*st]);       //将此实例类型添加到GT第一个位置上
//        for (int i = 0; i < comConst.NumOfTasks; ++i) {
//            if(chrom.InsAlcLst[chrom.TskSchLst[i]] == *st){
//                GL.insert(GL.end(), chrom.TskSchLst[i]);
//            }
//        }
//        chrom.Code.insert(GL);  //Code是set结构，其中元素GL首先按第一个元素（实例类型）排序（小的排前面），如果实例类型相同，则再按第二个元素（第一个任务）排序（小的在前）
//    }
//}

void GnrTwoDimCoding(chromosome &chrom){ //一维编码中的实例编号可以不连续
    set<int> insSet;
    map<int,int> InsToIndex;
    chrom.Code.clear();  //使用了push_back方法，因此为确保安全，在生成新的二维编码前要清空旧的二维编码！！-xy
//    vector<pair<int,int>> InsAndTypeVct;
    for(int i = 0; i < chrom.InsAlcLst.size(); ++i){
        insSet.insert(chrom.InsAlcLst[i]); //使用了set结构，故编号小的实例在前面
    }
    int index = 0;
    for (int ins: insSet ) {
        list<int> GL;
        GL.push_back(chrom.TypeLst[ins]);
        InsToIndex[ins] = index;  //InsToIndex[ins]为原一维编码中编号为ins的实例重新编号后的编号
        ++index;
        chrom.Code.push_back(GL); //Code是list结构，根据实例编号排序的
    }

    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        int ind = InsToIndex[chrom.InsAlcLst[chrom.TskSchLst[i]]]; //获得分配给第i个调度任务的实例的新编号
        list<list<int>>::iterator it = chrom.Code.begin();
        advance(it, ind);
        (*it).push_back(chrom.TskSchLst[i]);
    }
}

bool GnrOneDimCoding(chromosome &chrom){
    //生成的一维编码中，实例是从0连续编号的，但不一定满足RSGA1中实例分配编码要求的InsAlcLst[i]<=max{InsAlcLst[0],...,InsAlcLst[i-1]}+1-xy
    vector<list<int>> vecList;
    list<int>::iterator lst;
    vecList.assign(chrom.Code.begin(), chrom.Code.end());
    vector<int> InDg(comConst.NumOfTasks);
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        InDg[i] = Tasks[i].parents.size();
    }
//    chrom.TskSchLst.resize(comConst.NumOfTasks, -1);  chrom.InsAlcLst.resize(comConst.NumOfTasks, -1); chrom.TypeLst.resize(comConst.NumOfTasks, -1);
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        chrom.TskSchLst[i] = -1; chrom.InsAlcLst[i] = -1; chrom.TypeLst[i] = -1;
    }
    for(int j = 0; j < vecList.size(); ++j){
        chrom.TypeLst[j] = *(vecList[j].begin()); //对虚拟机实例的类型进行赋值，虚拟机实例的编号根据Code中的顺序依次从0开始编号
    }

    for(int i = 0; i < comConst.NumOfTasks; ++i){
        int j = 0;//xy2
        for(; j < vecList.size(); j++){
            if(vecList[j].size() > 1){
                lst = vecList[j].begin(); lst++;
                if (InDg[*lst] == 0) {
                    chrom.TskSchLst[i] = *lst;
                    chrom.InsAlcLst[*lst] = j;
                    for(int children: Tasks[*lst].children){
                        InDg[children] -= 1;    //更新其子任务的还没有被调度的父任务数量
                    }
                    vecList[j].erase(lst);      //删除虚拟机j上的已被调度的任务
                    break;
                }
            }
        }
        if (j==vecList.size()) {   //如果出现死锁情况，函数返回false
//            cout<<endl<<"TwoDimCoding is illegal, GnrOneDimCoding is fail";
            return false;
        }
    }
    return true;                    //如果转换成功，函数值返回true
}

void RecodeInsAndType(chromosome &chrom) {
    //对虚拟机实例及类型进行重新编码，满足RSGA1的编码要求，即：InsAlcLst[0]=0; InsAlcLst[i]<=max{InsAlcLst[0],...,InsAlcLst[i-1]}+1
    //同时对没有使用的实例的类型初始化为-1；
    chromosome oldChrom = chrom; //测试用
    map<int,list<int>> tskLstOfIns;
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        tskLstOfIns[chrom.InsAlcLst[i]].push_back(i);
    }
    set<list<int>> SetTskTy;
    for (auto elm: tskLstOfIns) {  //利用set进行排序, 即按list中的首元素（任务编号）从小到大进行排序-xy
        list<int> TsksAndTy = elm.second;
        TsksAndTy.push_back(chrom.TypeLst[elm.first]); //把处理这些任务的虚拟机类型添加到队列尾部
        SetTskTy.insert(TsksAndTy);
    }
    int nmb = 0;
    for (set<list<int>>::iterator it = SetTskTy.begin(); it != SetTskTy.end(); ++it) {
        list<int> lst = *it;
        list<int>::iterator lstIt = lst.end();  --lstIt;
        chrom.TypeLst[nmb] = *lstIt;
        for (list<int>::iterator pt = lst.begin(); pt != lstIt; ++pt) {
            chrom.InsAlcLst[*pt] = nmb;
        }
        ++nmb;
    }
    for (int j = tskLstOfIns.size(); j < comConst.NumOfTasks; ++j) {  //对没有使用的实例的类型初始化为-1;
        chrom.TypeLst[j] = -1;
    }
    // /*测试重新编码是否正确， 即：InsAlcLst[0]=0; InsAlcLst[i] <= max{InsAlcLst[0],...,InsAlcLst[i-1]}+1*/
    int maxId = 0;
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        if (chrom.InsAlcLst[i] > maxId) {
            cout<< endl<< "recode is wrong for insId"; exit(2);
        }
        maxId = XY_MAX(maxId,chrom.InsAlcLst[i]+1);
    }
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        if (chrom.TypeLst[chrom.InsAlcLst[i]] != oldChrom.TypeLst[oldChrom.InsAlcLst[i]]) {
            cout<<endl<<"recode is wrong for type"; exit(2);
        }
    }
}

//chromosome GnrChrByGRP_HEFT(){
//    //1、根据预算选择虚拟机实例
//    vector<VMType> InstTypeList;            //虚拟机类型集合
//    InstTypeList = VMTypes;                 //获得所有虚拟机类型
//    vector<VMType>::iterator vePo;          //InstTypeList 的一个迭代器
//    chromosome BestChrom;
//    BestChrom.MakeSpan = InfInt;
//    //生成任务调度顺序
//    vector<int> TskSchLst(comConst.NumOfTasks);
//    vector<double> ww(comConst.NumOfTasks, 0);
//    vector<double> rank_b(comConst.NumOfTasks, 0);
//    W_Cal_average(ww);
//    Calculate_Rank_b(rank_b,ww);
//    IndexSortDesc(TskSchLst,rank_b);
//
//    while(!InstTypeList.empty()){
//        double Remaining_Budget = Budget;
//        chromosome chrom;
//        initChrom(chrom);
//        chrom.TskSchLst = TskSchLst;
//        vector<VMType> Y;                   //最终选取的虚拟机实例集合
//        int n = 0;
//
//        //根据预算选取虚拟机建立Y集合
//        for(int j = InstTypeList.size() - 1; j < InstTypeList.size(); --j){     //存储的顺序是ps从小到大，所以从后向前选
//            if(Remaining_Budget >= InstTypeList[0].UnitCost){
//                n = floor(Remaining_Budget / InstTypeList[j].UnitCost);
//                for(int k = 0; k < n; k++){
//                    if(Y.size() >= comConst.NumOfTasks)
//                        break;
//                    Y.push_back(InstTypeList[j]);
//                }
//                Remaining_Budget -= InstTypeList[j].UnitCost * n;
//            } else
//                break;
//        }
//
//        //3、根据所选的虚拟机实例生成一种调度方案
//            Modified_HEFT(chrom, Y);
//
//        if(chrom.Cost <= Budget && chrom.MakeSpan <= BestChrom.MakeSpan){
//            BestChrom = chrom;
//        }
//        InstTypeList.erase(InstTypeList.end()-1);
//    }
//    //计算个体适应度值
//    if (BestChrom.MakeSpan > Deadline || BestChrom.Cost > Budget){
//        BestChrom.FitnessValue = BestChrom.MakeSpan / Deadline + BestChrom.Cost / Budget;
//    } else {
//        BestChrom.FitnessValue = (1 - inputParameter.theta) * BestChrom.MakeSpan / Deadline + inputParameter.theta * BestChrom.Cost / Budget;
//    }
//    return BestChrom;
//}

double FindIdleTimeSlot(set<double>& ITLofRscId, double& ExeTime, double& ReadyTime){
    set<double>::iterator pre  = ITLofRscId.begin();
    set<double>::iterator post = ITLofRscId.begin();
    ++post;
    while(post != ITLofRscId.end()) {
        if((*post - *pre) > ExeTime - PrecisionValue && ReadyTime - PrecisionValue < (*post)-ExeTime) {
            return  XY_MAX(*pre, ReadyTime);
        } else {
            ++pre; ++pre; ++post; ++post;
        }
    }
}

void UpdateITL(set<double>& ITLofRscId, double& StartTime, double& EndTime){
    int flag1 = 0, flag2 = 0;
    for (double x: ITLofRscId) {
        if ( fabs( x - StartTime ) < PrecisionValue ) {
            ITLofRscId.erase(x);
            flag1 = 1;
            break;
        }
    }
    if (flag1 == 0)
        ITLofRscId.insert(StartTime);

    for (double x: ITLofRscId) {
        if ( fabs( x - EndTime ) < PrecisionValue ) {
            ITLofRscId.erase(x);
            flag2 = 1;
            break;
        }
    }
    if (flag2 == 0)
        ITLofRscId.insert(EndTime);

//    if(ITLofRscId.find(StartTime) != ITLofRscId.end()) {
//        ITLofRscId.erase(StartTime);
//    } else {
//        ITLofRscId.insert(StartTime);
//    }
//    if(ITLofRscId.find(EndTime) != ITLofRscId.end()) {
//        ITLofRscId.erase(EndTime);
//    } else {
//        ITLofRscId.insert(EndTime);
//    }

}

void Modified_HEFT_XY(chromosome &chrom, vector<Instance> &Y, map<int,int> &NumOfInsUnused){
    double ms = -1;
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        chrom.InsAlcLst[i] = -1; chrom.TypeLst[i] = -1;
    }
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        int InsId1 = -1, InsId2 = -1;
        int taskId = chrom.TskSchLst[i];
        //选择虚拟机
        double fnlEndTime1 = 100000000000, fnlEndTime2 = 100000000000, fnlStrTime1 = 0, fnlStrTime2 = 0;
        int NumOfNdAdd = -1;
        for (int j = 0; j < Y.size(); ++j) {
            //此实例的类型可以处理该任务,并且未被合并
            if(Y[j].state != 2 && find(Tasks[taskId].ElgRsc.begin(), Tasks[taskId].ElgRsc.end(), Y[j].type) != Tasks[taskId].ElgRsc.end()){
                double rt = 0;
                for (int parentId: Tasks[taskId].parents) {
                    int parentInsId = chrom.InsAlcLst[parentId];
                    double Max = chrom.EndTime[parentId];
                    if(j != parentInsId){ //Y中的下标表示实例编号，
                        Max += ParChildTranFileSizeSum[parentId][taskId] / XY_MIN(VMTypes[Y[j].type].bw, VMTypes[Y[parentInsId].type].bw); //由于存储类型的VMTypes的下标（索引）与类型相同，故用Y[j].type代替下标（索引）-xy
                    }
                    if (rt + PrecisionValue < Max){
                        rt = Max;
                    }
                }
                //执行时间
                double ExeTime = Tasks[taskId].length / VMTypes[Y[j].type].ps; //表明数据传输不需要占用虚拟机！原始程序也是这样的吗？-xy
                //搜索可用时间片段
                double StrTime = FindIdleTimeSlot(Y[j].ITL, ExeTime, rt);
                double EndTime = StrTime + ExeTime;
                double NewLST = XY_MIN(StrTime - VMTypes[Y[j].type].BootTime, Y[j].LST);
                double NewLET = XY_MAX(EndTime, Y[j].LET);
                int TemNumOfNdAdd = ceil((NewLET - NewLST - PrecisionValue)/VMTypes[Y[j].type].ChargeUnit) - Y[j].state * ceil((Y[j].LET - Y[j].LST - PrecisionValue)/VMTypes[Y[j].type].ChargeUnit);
//                int TemNumOfNdAdd = ceil((NewLET - NewLST)/VMTypes[Y[j].type].ChargeUnit) - XY_MAX(Y[j].state * ceil((Y[j].LET - Y[j].LST)/VMTypes[Y[j].type].ChargeUnit), 1);
                //找最快完成
                if ( EndTime + PrecisionValue < fnlEndTime1 ) {
                    fnlStrTime1 = StrTime; fnlEndTime1 = EndTime;
                    InsId1 = j;
                    NumOfNdAdd = TemNumOfNdAdd;
                }
                //找不改变成本（不需要补偿实例）情况下执行时间最小的那个虚拟机
                if ( EndTime + PrecisionValue < fnlEndTime2 && TemNumOfNdAdd == 0 ){
                    fnlStrTime2 = StrTime; fnlEndTime2 = EndTime;
                    InsId2 = j;
                }
            }
        }

        if(fnlEndTime1 + PrecisionValue < fnlEndTime2){
            if(NumOfNdAdd <= NumOfInsUnused[Y[InsId1].type]){  //在Y中找是否有足够的补偿实例
                int x = NumOfNdAdd;
                if ( Y[InsId1].state == 0 )  x = NumOfNdAdd -1;
                if ( x > 0 ) {
                    for(int j = Y.size() -1; j >= 0; --j){  //在Y中合并所需数量的还未被使用的同类找
                        if( Y[j].type == Y[InsId1].type && Y[j].state == 0 && j != InsId1 ){
                            Y[j].state = 2; --x;
                            if (x == 0 ) break;
                        }
                    }
                }
            } else {
                if(InsId2 != -1){         //I*min != null
                    chrom.StartTime[taskId] = fnlStrTime2;
                    chrom.EndTime[taskId] = fnlEndTime2;
                    Y[InsId2].AssignTskLst.push_back(taskId);
                    chrom.InsAlcLst[taskId] = InsId2;
                    UpdateITL(Y[InsId2].ITL, fnlStrTime2, fnlEndTime2);
                    Y[InsId2].LST = XY_MIN(fnlStrTime2 - VMTypes[Y[InsId2].type].BootTime, Y[InsId2].LST);
                    Y[InsId2].LET = XY_MAX(fnlEndTime2, Y[InsId2].LET);
                    Y[InsId2].state = 1;
                    ms = XY_MAX(ms, fnlEndTime2);
                }
            }
        }
        if(chrom.InsAlcLst[taskId] == -1){    //ti is not assign  未分配此任务实例为-1
            chrom.StartTime[taskId] = fnlStrTime1;
            chrom.EndTime[taskId] = fnlEndTime1;
            Y[InsId1].AssignTskLst.push_back(taskId);
            chrom.InsAlcLst[taskId] = InsId1;
            UpdateITL(Y[InsId1].ITL, fnlStrTime1, fnlEndTime1);
            Y[InsId1].LST = XY_MIN(fnlStrTime1 - VMTypes[Y[InsId1].type].BootTime, Y[InsId1].LST);
            Y[InsId1].LET = XY_MAX(fnlEndTime1, Y[InsId1].LET);
            Y[InsId1].state = 1;
            NumOfInsUnused[Y[InsId1].type] = NumOfInsUnused[Y[InsId1].type] - NumOfNdAdd;
            if ( NumOfInsUnused[Y[InsId1].type] < 0 ) {
                for ( auto Ins: Y) {
                    if (Ins.state == 0 && Ins.type == Y[InsId1].type)
                        Ins.state = 2;
                }
            }
            ms = XY_MAX(ms, fnlEndTime1);
        }
    }
    chromosome oldChrom1 = chrom;

    //计算成本，对实例重新编码，更新染色体编码
    double cost =0;
    int nmb = 0;
    for ( int j = 0; j < Y.size(); ++j ) {
        if ( Y[j].state == 1 ) {
            double LeaseTime = XY_MAX(Y[j].LET - Y[j].LST, VMTypes[Y[j].type].MinLeaseTime);
            cost = cost + VMTypes[Y[j].type].FixedCost + VMTypes[Y[j].type].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[Y[j].type].ChargeUnit);
            for (int TskId: Y[j].AssignTskLst) {
                chrom.InsAlcLst[TskId] = nmb;
            }
            chrom.TypeLst[nmb] = Y[j].type;
            ++nmb;
        }
    }
    chromosome oldChrom2 = chrom;

    DecodeChrom(chrom); //解码验证

    if ( fabs(chrom.MakeSpan - ms) >= PrecisionValue ) {
        cout << endl << "oldChrom1: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << oldChrom1.InsAlcLst[i] << ", ";
        }
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << oldChrom1.TypeLst[i] << ", ";
        }
        cout << endl << "Chrom: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << chrom.InsAlcLst[i] << ", ";
        }
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << chrom.TypeLst[i] << ", ";
        }
        cout << endl << "calucating makespan is wrong!" ;
        cout << endl << chrom.MakeSpan << ",  ms: " << ms;
    }

    if ( fabs(chrom.Cost - cost) >= PrecisionValue ) {
        cout << endl << "oldChrom1: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << oldChrom1.InsAlcLst[i] << ", ";
        }
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << oldChrom1.TypeLst[i] << ", ";
        }
        cout << endl << "Chrom: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << chrom.InsAlcLst[i] << ", ";
        }
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << chrom.TypeLst[i] << ", ";
        }
        cout << endl << "calucating cost is wrong!" ;
        cout << endl << "decode: " << chrom.Cost << ",  cs: " << cost;
        cout << endl << "Y[j]: " ;
        for ( int j = 0; j < Y.size(); ++j ) {
            if (Y[j].state == 1 ) {
                cout << endl << j << ": ";
                for (auto x: Y[j].ITL) {
                    cout << x << "; ";
                }
            }
        }
        DecodeChrom(oldChrom2); //解码验证
    }

    chrom.MakeSpan = ms;
    chrom.Cost = cost;
}

int Modified_HEFT_XY_1(chromosome &chrom, vector<Instance> &Y, map<int,int> &NumOfInsUnused){
    //按完成时间从小到大分配给不增加成本的实例，如果不存在不增加成本的实例，则完成时间从小到大分配给存在补偿实例的实例，否则返回无可行调度方案；
    double ms = -1;
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        chrom.InsAlcLst[i] = -1; chrom.TypeLst[i] = -1;
    }
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        int InsId1 = -1, InsId2 = -1, InsId3 = -1;
        int taskId = chrom.TskSchLst[i];

        //选择虚拟机
        double fnlEndTime1 = InfInt, fnlStrTime1 = 0; //record the earliest finish time and corresponding start time
        double fnlEndTime2 = InfInt, fnlStrTime2 = 0; //record the earliest finish time without changing the cost
        double fnlEndTime3 = InfInt, fnlStrTime3 = 0; //record the earliest finish time with adding the cost but existing enough unused instances
        int NumOfNdAdd = -1, NumOfNdAdd3 = -1;
        for (int j = 0; j < Y.size(); ++j) {
            //此实例的类型可以处理该任务,并且未被合并
            if(Y[j].state != 2 && find(Tasks[taskId].ElgRsc.begin(), Tasks[taskId].ElgRsc.end(), Y[j].type) != Tasks[taskId].ElgRsc.end()){
                double rt = 0;
                for (int parentId: Tasks[taskId].parents) {
                    int parentInsId = chrom.InsAlcLst[parentId];
                    double Max = chrom.EndTime[parentId];
                    if(j != parentInsId){ //Y中的下标表示实例编号，
                        Max += ParChildTranFileSizeSum[parentId][taskId] / XY_MIN(VMTypes[Y[j].type].bw, VMTypes[Y[parentInsId].type].bw); //由于存储类型的VMTypes的下标（索引）与类型相同，故用Y[j].type代替下标（索引）-xy
                    }
                    if (rt + PrecisionValue < Max){
                        rt = Max;
                    }
                }
                //执行时间
                double ExeTime = Tasks[taskId].length / VMTypes[Y[j].type].ps; //表明数据传输不需要占用虚拟机！原始程序也是这样的吗？-xy
                //搜索可用时间片段
                double StrTime = FindIdleTimeSlot(Y[j].ITL, ExeTime, rt);
                double EndTime = StrTime + ExeTime;
                double NewLST = XY_MIN(StrTime - VMTypes[Y[j].type].BootTime, Y[j].LST);
                double NewLET = XY_MAX(EndTime, Y[j].LET);
                int TemNumOfNdAdd = ceil((NewLET - NewLST - PrecisionValue)/VMTypes[Y[j].type].ChargeUnit) - Y[j].state * ceil((Y[j].LET - Y[j].LST - PrecisionValue)/VMTypes[Y[j].type].ChargeUnit);
//                int TemNumOfNdAdd = ceil((NewLET - NewLST)/VMTypes[Y[j].type].ChargeUnit) - XY_MAX(Y[j].state * ceil((Y[j].LET - Y[j].LST)/VMTypes[Y[j].type].ChargeUnit), 1);
                //找最快完成
//                if ( EndTime + PrecisionValue < fnlEndTime1 ) {
//                    fnlStrTime1 = StrTime; fnlEndTime1 = EndTime;
//                    InsId1 = j;
//                    NumOfNdAdd = TemNumOfNdAdd;
//                }
                //找不改变成本（不需要补偿实例）情况下执行时间最小的那个虚拟机
                if ( EndTime + PrecisionValue < fnlEndTime2 && TemNumOfNdAdd == 0 ){
                    fnlStrTime2 = StrTime; fnlEndTime2 = EndTime;
                    InsId2 = j;
                }
                //找存在补偿实例情况下执行时间最小的那个虚拟机
                if ( EndTime + PrecisionValue < fnlEndTime3 && TemNumOfNdAdd <= NumOfInsUnused[Y[j].type] ){
                    fnlStrTime3 = StrTime; fnlEndTime3 = EndTime;
                    InsId3 = j;
                    NumOfNdAdd3 = TemNumOfNdAdd;
                }
            }
        }

        if(InsId2 != -1){         //I*min != null, 存在不存在不增加成本的实例
            chrom.StartTime[taskId] = fnlStrTime2;
            chrom.EndTime[taskId] = fnlEndTime2;
            Y[InsId2].AssignTskLst.push_back(taskId);
            chrom.InsAlcLst[taskId] = InsId2;
            UpdateITL(Y[InsId2].ITL, fnlStrTime2, fnlEndTime2);
            Y[InsId2].LST = XY_MIN(fnlStrTime2 - VMTypes[Y[InsId2].type].BootTime, Y[InsId2].LST);
            Y[InsId2].LET = XY_MAX(fnlEndTime2, Y[InsId2].LET);
            Y[InsId2].state = 1;
            ms = XY_MAX(ms, fnlEndTime2);
        } else {
            if(InsId3 != -1) { // 存在补偿实例
                int x = NumOfNdAdd3;
                if ( Y[InsId3].state == 0 )  x = NumOfNdAdd3 -1;
                if ( x > 0 ) {
                    for(int j = Y.size() -1; j >= 0; --j){  //在Y中合并所需数量的还未被使用的同类
                        if( Y[j].type == Y[InsId3].type && Y[j].state == 0 && j != InsId3 ){
                            Y[j].state = 2; --x;
                            if (x == 0 ) break;
                        }
                    }
                }
                chrom.StartTime[taskId] = fnlStrTime3;
                chrom.EndTime[taskId] = fnlEndTime3;
                Y[InsId3].AssignTskLst.push_back(taskId);
                chrom.InsAlcLst[taskId] = InsId3;
                UpdateITL(Y[InsId3].ITL, fnlStrTime3, fnlEndTime3);
                Y[InsId3].LST = XY_MIN(fnlStrTime3 - VMTypes[Y[InsId3].type].BootTime, Y[InsId3].LST);
                Y[InsId3].LET = XY_MAX(fnlEndTime3, Y[InsId3].LET);
                Y[InsId3].state = 1;
                NumOfInsUnused[Y[InsId3].type] = NumOfInsUnused[Y[InsId3].type] - NumOfNdAdd3;
                if ( NumOfInsUnused[Y[InsId3].type] < 0 ) {
                    cout << endl << "NumOfInsUnused is less then 0; ";
                }
                ms = XY_MAX(ms, fnlEndTime3);
            } else { //不存在补偿实例,分配给最早完成的
                return 0;
//                chrom.StartTime[taskId] = fnlStrTime1;
//                chrom.EndTime[taskId] = fnlEndTime1;
//                Y[InsId1].AssignTskLst.push_back(taskId);
//                chrom.InsAlcLst[taskId] = InsId1;
//                UpdateITL(Y[InsId1].ITL, fnlStrTime1, fnlEndTime1);
//                Y[InsId1].LST = XY_MIN(fnlStrTime1 - VMTypes[Y[InsId1].type].BootTime, Y[InsId1].LST);
//                Y[InsId1].LET = XY_MAX(fnlEndTime1, Y[InsId1].LET);
//                Y[InsId1].state = 1;
//                NumOfInsUnused[Y[InsId1].type] = NumOfInsUnused[Y[InsId1].type] - NumOfNdAdd;
//                for ( auto Ins: Y) {
//                    if (Ins.state == 0 && Ins.type == Y[InsId1].type)
//                        Ins.state = 2;
//                }
//                ms = XY_MAX(ms, fnlEndTime1);
            }
        }
    }
    chromosome oldChrom1 = chrom;

    //计算成本，对实例重新编码，更新染色体编码
    double cost =0;
    int nmb = 0;
    for ( int j = 0; j < Y.size(); ++j ) {
        if ( Y[j].state == 1 ) {
            double LeaseTime = XY_MAX(Y[j].LET - Y[j].LST, VMTypes[Y[j].type].MinLeaseTime);
            cost = cost + VMTypes[Y[j].type].FixedCost + VMTypes[Y[j].type].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[Y[j].type].ChargeUnit);
            for (int TskId: Y[j].AssignTskLst) {
                chrom.InsAlcLst[TskId] = nmb;
            }
            chrom.TypeLst[nmb] = Y[j].type;
            ++nmb;
        }
    }
    chromosome oldChrom2 = chrom;

    DecodeChrom(chrom); //解码验证

    if ( fabs(chrom.MakeSpan - ms) >= PrecisionValue ) {
        cout << endl << "oldChrom1: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << oldChrom1.InsAlcLst[i] << ", ";
        }
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << oldChrom1.TypeLst[i] << ", ";
        }
        cout << endl << "Chrom: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << chrom.InsAlcLst[i] << ", ";
        }
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << chrom.TypeLst[i] << ", ";
        }
        cout << endl << "calucating makespan is wrong!" << endl;
        cout << chrom.MakeSpan << ",  ms: " << ms <<endl;
    }

    if ( fabs(chrom.Cost - cost) >= PrecisionValue ) {
        cout << "oldChrom1: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << oldChrom1.InsAlcLst[i] << ", ";
        }
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << oldChrom1.TypeLst[i] << ", ";
        }
        cout << endl << "Chrom: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << chrom.InsAlcLst[i] << ", ";
        }
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++) {
            cout << chrom.TypeLst[i] << ", ";
        }
        cout << endl << "calucating cost is wrong!" << endl;
        cout << "decode: " << chrom.Cost << ",  cs: " << cost <<endl;
        cout <<  "Y[j]: " ;
        for ( int j = 0; j < Y.size(); ++j ) {
            if (Y[j].state == 1 ) {
                cout << endl << j << ": ";
                for (auto x: Y[j].ITL) {
                    cout << x << "; ";
                }
            }
        }
        cout << endl;
        DecodeChrom(oldChrom2); //解码验证
    }

    chrom.MakeSpan = ms;
    chrom.Cost = cost;
    return 1;
}

int Modified_HEFT_XY_2(chromosome &chrom, vector<Instance> &Y, map<int,int> &NumOfInsUnused){
    //按完成时间从小到大分配给不增加成本或存在补偿实例的实例，否则返回无可行调度方案；
    double ms = -1;
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        chrom.InsAlcLst[i] = -1; chrom.TypeLst[i] = -1;
    }
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        int InsId2 = -1, InsId3 = -1;
        int taskId = chrom.TskSchLst[i];

        //选择虚拟机
        double fnlEndTime2 = InfInt, fnlStrTime2 = 0; //record the earliest finish time without changing the cost
        double fnlEndTime3 = InfInt, fnlStrTime3 = 0; //record the earliest finish time with adding the cost but existing enough unused instances
        int NumOfNdAdd3 = -1;
        for (int j = 0; j < Y.size(); ++j) {
            //此实例的类型可以处理该任务,并且未被合并
            if (Y[j].state != 2 && Tasks[taskId].ElgRsc.find(Y[j].type) != Tasks[taskId].ElgRsc.end()) {
                double rt = 0;
                for (int parentId: Tasks[taskId].parents) {
                    int parentInsId = chrom.InsAlcLst[parentId];
                    double Max = chrom.EndTime[parentId];
                    if (j != parentInsId) { //Y中的下标表示实例编号，
                        Max += ParChildTranFileSizeSum[parentId][taskId] / XY_MIN(VMTypes[Y[j].type].bw, VMTypes[Y[parentInsId].type].bw);
                    }
                    if (rt + PrecisionValue < Max) {
                        rt = Max;
                    }
                }
                //执行时间
                double ExeTime = Tasks[taskId].length / VMTypes[Y[j].type].ps; //表明数据传输不需要占用虚拟机！原始程序也是这样的吗？-xy
                //搜索可用时间片段
                double StrTime = FindIdleTimeSlot(Y[j].ITL, ExeTime, rt);
                double EndTime = StrTime + ExeTime;
                double NewLST = XY_MIN(StrTime - VMTypes[Y[j].type].BootTime, Y[j].LST);
                double NewLET = XY_MAX(EndTime, Y[j].LET);
                int TemNumOfNdAdd = ceil((NewLET - NewLST - PrecisionValue) / VMTypes[Y[j].type].ChargeUnit) -
                                    Y[j].state *
                                    ceil((Y[j].LET - Y[j].LST - PrecisionValue) / VMTypes[Y[j].type].ChargeUnit);
                //找不改变成本（不需要补偿实例）情况下执行时间最小的那个虚拟机
                if ( EndTime + PrecisionValue < fnlEndTime2 && TemNumOfNdAdd == 0 ){
                    fnlStrTime2 = StrTime; fnlEndTime2 = EndTime;
                    InsId2 = j;
                }
                //找不改变成本或存在补偿实例情况下执行时间最小的那个虚拟机
                if (EndTime + PrecisionValue < fnlEndTime3 && TemNumOfNdAdd <= NumOfInsUnused[Y[j].type]) {
                    fnlStrTime3 = StrTime;
                    fnlEndTime3 = EndTime;
                    InsId3 = j;
                    NumOfNdAdd3 = TemNumOfNdAdd;
                }
            }
        }

        if (fabs(fnlEndTime2 - fnlEndTime3) < PrecisionValue && fnlEndTime2 < InfInt-1) { //存在补偿实例与不增加成本的最小时间相同
            chrom.StartTime[taskId] = fnlStrTime2;
            chrom.EndTime[taskId] = fnlEndTime2;
            Y[InsId2].AssignTskLst.push_back(taskId);
            chrom.InsAlcLst[taskId] = InsId2; //虽然后面重新编码时有赋值，但此处也必须要有，因为在计算rt时需要
            UpdateITL(Y[InsId2].ITL, fnlStrTime2, fnlEndTime2);
            Y[InsId2].LST = XY_MIN(fnlStrTime2 - VMTypes[Y[InsId2].type].BootTime, Y[InsId2].LST);
            Y[InsId2].LET = XY_MAX(fnlEndTime2, Y[InsId2].LET);
            if (Y[InsId2].state !=1 ) { //测试用-xy
                cout<<endl<<"InsId2 is not used, it wrong!"; exit(2);
            }
            Y[InsId2].state = 1;  //可以不需要！
            ms = XY_MAX(ms, fnlEndTime2);
        } else {
            if (InsId3 != -1) { // 存在补偿实例或不增加成本
                int x = NumOfNdAdd3;
                if (Y[InsId3].state == 0) x = NumOfNdAdd3 - 1;
                if (x > 0) {
                    for (int j = Y.size() - 1; j >= 0; --j) {  //在Y中合并所需数量的还未被使用的同类
                        if (Y[j].type == Y[InsId3].type && Y[j].state == 0 && j != InsId3) {
                            Y[j].state = 2;
                            --x;
                            if (x == 0) break;
                        }
                    }
                }
                chrom.StartTime[taskId] = fnlStrTime3;
                chrom.EndTime[taskId] = fnlEndTime3;
                Y[InsId3].AssignTskLst.push_back(taskId);
                chrom.InsAlcLst[taskId] = InsId3;
                UpdateITL(Y[InsId3].ITL, fnlStrTime3, fnlEndTime3);
                Y[InsId3].LST = XY_MIN(fnlStrTime3 - VMTypes[Y[InsId3].type].BootTime, Y[InsId3].LST);
                Y[InsId3].LET = XY_MAX(fnlEndTime3, Y[InsId3].LET);
                Y[InsId3].state = 1;
                NumOfInsUnused[Y[InsId3].type] = NumOfInsUnused[Y[InsId3].type] - NumOfNdAdd3;
                if (NumOfInsUnused[Y[InsId3].type] < 0) {
                    cout << endl << "NumOfInsUnused is less then 0; ";
                }
                ms = XY_MAX(ms, fnlEndTime3);
            } else { //不存在补偿实例或不增加成本，即无可行方案
                return 0;
            }
        }
    }
    chromosome oldChrom1 = chrom;
    //计算成本，对实例重新编码，更新染色体编码，使用的实例是0开始连续编码，但不一定满足RSGA1的编码要求。未使用的实例的类型为-1。
    double cost =0;
    int nmb = 0;
    for ( int j = 0; j < Y.size(); ++j ) {
        if ( Y[j].state == 1 ) {
            double LeaseTime = XY_MAX(Y[j].LET - Y[j].LST, VMTypes[Y[j].type].MinLeaseTime);
            cost = cost + VMTypes[Y[j].type].FixedCost + VMTypes[Y[j].type].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[Y[j].type].ChargeUnit);
            for (int TskId: Y[j].AssignTskLst) {
                chrom.InsAlcLst[TskId] = nmb;  //使用的实例是0开始连续编码，但不一定满足RSGA1的编码要求：InsAlcLst[0]=0; InsAlcLst[i]<=max{InsAlcLst[0],...,InsAlcLst[i-1]}+1；
            }
            chrom.TypeLst[nmb] = Y[j].type;
            ++nmb;
        }
    }
    DecodeChrom(chrom); //解码验证，测试用-xy

    if ( fabs(chrom.MakeSpan - ms) >= PrecisionValue ) {
        cout << endl << "oldChrom1: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)   cout << oldChrom1.InsAlcLst[i] << ", ";
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)  cout << oldChrom1.TypeLst[i] << ", ";
        cout << endl << "Chrom: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)  cout << chrom.InsAlcLst[i] << ", ";
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)  cout << chrom.TypeLst[i] << ", ";
        cout << endl << "calucating makespan is wrong!" << endl;
        cout << chrom.MakeSpan << ",  ms: " << ms <<endl;
    }

    if ( fabs(chrom.Cost - cost) >= PrecisionValue ) {
        cout << "oldChrom1: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)  cout << oldChrom1.InsAlcLst[i] << ", ";
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)  cout << oldChrom1.TypeLst[i] << ", ";
        cout << endl << "Chrom: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)   cout << chrom.InsAlcLst[i] << ", ";
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)  cout << chrom.TypeLst[i] << ", ";
        cout << endl << "calucating cost is wrong!" << endl;
        cout << "decode: " << chrom.Cost << ",  cs: " << cost << endl;
        cout <<  "Y[j]: " ;
        for ( int j = 0; j < Y.size(); ++j ) {
            if (Y[j].state == 1 ) {
                cout << endl << j << ": ";
                for (auto x: Y[j].ITL) {
                    cout << x << "; ";
                }
            }
        }
        cout << endl;
    }

    chrom.MakeSpan = ms;
    chrom.Cost = cost;
    return 1;
}

int Modified_HEFT_XY_Imp(chromosome &chrom, vector<Instance> &Y, map<int,int> &NumOfInsUnused){
    //对Modified_HEFT_XY_2的改进，在遍历访问实例时增加了当前遍历的实例与最近访问遍历的实例是否类型相同且都是未被使用的判断，如果是则可以跳过当前的实例，以提高算法效率！
    //按完成时间从小到大分配给不增加成本或存在补偿实例的实例，否则返回无可行调度方案；
    double ms = -1;
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        chrom.InsAlcLst[i] = -1; chrom.TypeLst[i] = -1;
    }
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        int InsId2 = -1, InsId3 = -1;
        int taskId = chrom.TskSchLst[i];

        //选择虚拟机
        double fnlEndTime2 = InfInt, fnlStrTime2 = 0; //record the earliest finish time without changing the cost
        double fnlEndTime3 = InfInt, fnlStrTime3 = 0; //record the earliest finish time with adding the cost but existing enough unused instances
        int NumOfNdAdd3 = -1;
        int lstTy = -1, flg = 1;  //lstTy记录上次最近遍历的实例的类型，flg记录上次最近遍历的实例是否被使用过，如果被使用过则其为1，否则为0
        for (int j = 0; j < Y.size(); ++j) {
            if (Y[j].type == lstTy && Y[j].state == 0 && flg == 0) continue; //如果当前考察的实例与上次最近考察的类型相同且都未被使用，则跳过，以提高效率-xy
            if (Y[j].state != 2 && Tasks[taskId].ElgRsc.find(Y[j].type) != Tasks[taskId].ElgRsc.end()) { //此实例的类型可以处理该任务,并且未被合并
                lstTy = Y[j].type;
                if (Y[j].state == 0)
                    flg = 0;
                else
                    flg = 1;
                double rt = 0;
                for (int parentId: Tasks[taskId].parents) {
                    int parentInsId = chrom.InsAlcLst[parentId];
                    double Max = chrom.EndTime[parentId];
                    if (j != parentInsId) { //Y中的下标表示实例编号，
                        Max += ParChildTranFileSizeSum[parentId][taskId] / XY_MIN(VMTypes[Y[j].type].bw, VMTypes[Y[parentInsId].type].bw);
                    }
                    if (rt + PrecisionValue < Max) {
                        rt = Max;
                    }
                }
                //执行时间
                double ExeTime = Tasks[taskId].length / VMTypes[Y[j].type].ps; //表明数据传输不需要占用虚拟机！原始程序也是这样的吗？-xy
                //搜索可用时间片段
                double StrTime = FindIdleTimeSlot(Y[j].ITL, ExeTime, rt);
                double EndTime = StrTime + ExeTime;
                double NewLST = XY_MIN(StrTime - VMTypes[Y[j].type].BootTime, Y[j].LST);
                double NewLET = XY_MAX(EndTime, Y[j].LET);
                int TemNumOfNdAdd = ceil((NewLET - NewLST - PrecisionValue) / VMTypes[Y[j].type].ChargeUnit) -
                                    Y[j].state *
                                    ceil((Y[j].LET - Y[j].LST - PrecisionValue) / VMTypes[Y[j].type].ChargeUnit);
                //找不改变成本（不需要补偿实例）情况下执行时间最小的那个虚拟机
                if ( EndTime + PrecisionValue < fnlEndTime2 && TemNumOfNdAdd == 0 ){
                    fnlStrTime2 = StrTime; fnlEndTime2 = EndTime;
                    InsId2 = j;
                }
                //找不改变成本或存在补偿实例情况下执行时间最小的那个虚拟机
                if (EndTime + PrecisionValue < fnlEndTime3 && TemNumOfNdAdd <= NumOfInsUnused[Y[j].type]) {
                    fnlStrTime3 = StrTime; fnlEndTime3 = EndTime;
                    InsId3 = j;
                    NumOfNdAdd3 = TemNumOfNdAdd;
                }
            }
        }

        if (fabs(fnlEndTime2 - fnlEndTime3) < PrecisionValue && fnlEndTime2 < InfInt - 1) { //存在补偿实例与不增加成本的最小时间相同
            chrom.StartTime[taskId] = fnlStrTime2;
            chrom.EndTime[taskId] = fnlEndTime2;
            Y[InsId2].AssignTskLst.push_back(taskId);
            chrom.InsAlcLst[taskId] = InsId2; //虽然后面重新编码时有赋值，但此处也必须要有，因为在计算rt时需要
            UpdateITL(Y[InsId2].ITL, fnlStrTime2, fnlEndTime2);
            Y[InsId2].LST = XY_MIN(fnlStrTime2 - VMTypes[Y[InsId2].type].BootTime, Y[InsId2].LST);
            Y[InsId2].LET = XY_MAX(fnlEndTime2, Y[InsId2].LET);
            if (Y[InsId2].state != 1 ) { //测试用-xy
                cout<<endl<<"InsId2 is not used, it wrong!"; exit(2);
            }
            Y[InsId2].state = 1;  //可以不需要！
            ms = XY_MAX(ms, fnlEndTime2);
        } else {
            if (InsId3 != -1) { // 存在补偿实例或不增加成本
                int x = NumOfNdAdd3;
                if (Y[InsId3].state == 0) x = NumOfNdAdd3 - 1;
                if (x > 0) {
                    for (int j = Y.size() - 1; j >= 0; --j) {  //在Y中合并所需数量的还未被使用的同类
                        if (Y[j].type == Y[InsId3].type && Y[j].state == 0 && j != InsId3) {
                            Y[j].state = 2;
                            --x;
                            if (x == 0) break;
                        }
                    }
                }
                chrom.StartTime[taskId] = fnlStrTime3;
                chrom.EndTime[taskId] = fnlEndTime3;
                Y[InsId3].AssignTskLst.push_back(taskId);
                chrom.InsAlcLst[taskId] = InsId3;
                UpdateITL(Y[InsId3].ITL, fnlStrTime3, fnlEndTime3);
                Y[InsId3].LST = XY_MIN(fnlStrTime3 - VMTypes[Y[InsId3].type].BootTime, Y[InsId3].LST);
                Y[InsId3].LET = XY_MAX(fnlEndTime3, Y[InsId3].LET);
                Y[InsId3].state = 1;
                NumOfInsUnused[Y[InsId3].type] = NumOfInsUnused[Y[InsId3].type] - NumOfNdAdd3;
                if (NumOfInsUnused[Y[InsId3].type] < 0) {
                    cout << endl << "NumOfInsUnused is less then 0; ";
                }
                ms = XY_MAX(ms, fnlEndTime3);
            } else { //不存在补偿实例或不增加成本，即无可行方案
                return 0;
            }
        }
    }
    chromosome oldChrom1 = chrom;
    //计算成本，对实例重新编码，更新染色体编码，使用的实例是0开始连续编码，但不一定满足RSGA1的编码要求。未使用的实例的类型为-1。
    double cost =0;
    int nmb = 0;
    for ( int j = 0; j < Y.size(); ++j ) {
        if ( Y[j].state == 1 ) {
            double LeaseTime = XY_MAX(Y[j].LET - Y[j].LST, VMTypes[Y[j].type].MinLeaseTime);
            cost = cost + VMTypes[Y[j].type].FixedCost + VMTypes[Y[j].type].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[Y[j].type].ChargeUnit);
            for (int TskId: Y[j].AssignTskLst) {
                chrom.InsAlcLst[TskId] = nmb;  //使用的实例是0开始连续编码，但不一定满足RSGA1的编码要求：InsAlcLst[0]=0; InsAlcLst[i]<=max{InsAlcLst[0],...,InsAlcLst[i-1]}+1；
            }
            chrom.TypeLst[nmb] = Y[j].type;
            ++nmb;
        }
    }
    DecodeChrom(chrom); //解码验证，测试用-xy

    if ( fabs(chrom.MakeSpan - ms) >= PrecisionValue ) {
        cout << endl << "oldChrom1: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)   cout << oldChrom1.InsAlcLst[i] << ", ";
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)  cout << oldChrom1.TypeLst[i] << ", ";
        cout << endl << "Chrom: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)  cout << chrom.InsAlcLst[i] << ", ";
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)  cout << chrom.TypeLst[i] << ", ";
        cout << endl << "calucating makespan is wrong!" << endl;
        cout << chrom.MakeSpan << ",  ms: " << ms <<endl;
    }

    if ( fabs(chrom.Cost - cost) >= PrecisionValue ) {
        cout << "oldChrom1: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)  cout << oldChrom1.InsAlcLst[i] << ", ";
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)  cout << oldChrom1.TypeLst[i] << ", ";
        cout << endl << "Chrom: " << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)   cout << chrom.InsAlcLst[i] << ", ";
        cout << endl;
        for (int i = 0; i < comConst.NumOfTasks; i++)  cout << chrom.TypeLst[i] << ", ";
        cout << endl << "calucating cost is wrong!" << endl;
        cout << "decode: " << chrom.Cost << ",  cs: " << cost << endl;
        cout <<  "Y[j]: " ;
        for ( int j = 0; j < Y.size(); ++j ) {
            if (Y[j].state == 1 ) {
                cout << endl << j << ": ";
                for (auto x: Y[j].ITL) {
                    cout << x << "; ";
                }
            }
        }
        cout << endl;
    }

    chrom.MakeSpan = ms;
    chrom.Cost = cost;
    return 1;
}

//void Modified_HEFT_f(chromosome &chrom, vector<VMType> &Y){
//    vector<set<double> > vatL;//虚拟机初始化
//    vector<int> flag(Y.size(), 0);
//    double max_EndTi = -100000;
//
//    for(int j = 0; j < comConst.NumOfTasks; ++j){
//        chrom.InsAlcLst[j] = -1;
//        chrom.TypeLst[j] = -1;
//    }
//
//    vector<int> InsTimeLine(Y.size(), 3600);
//
//    for (int j = 0; j < Y.size(); ++j) {
//        set<double> a;
//        a.insert(Y[j].BootTime);
//        a.insert(InfInt * 1.0);
//        vatL.push_back(a);
//    }
//
//    for (int i = 0; i < comConst.NumOfTasks; ++i) {
//        int vmIndex1 = -1;
//        int vmIndex2 = -1;
//        int taskIndex = chrom.TskSchLst[i];
//        //选择虚拟机
//        double finalEndTime1 = InfInt;
//        double finalEndTime2 = InfInt;
//        double finalStartTime1 = 0;
//        double finalStartTime2 = 0;
//        for (int j = 0; j < Y.size(); ++j) {
//            //先判断此类型是否为该任务可得的虚拟机类型
//            if(find(Tasks[taskIndex].ElgRsc.begin(), Tasks[taskIndex].ElgRsc.end(), Y[j].id - 1) != Tasks[taskIndex].ElgRsc.end()){
//                double rt = 0;
//                for (int i = 0; i < Tasks[taskIndex].parents.size(); ++i) {
//                    int parentIndex = Tasks[taskIndex].parents[i];
//                    int parentVMIndex = chrom.InsAlcLst[parentIndex];
//                    double Max = chrom.EndTime[parentIndex];
//                    if(j != parentVMIndex){
//                        double sum = ParChildTranFileSizeSum[parentIndex][taskIndex];
//                        Max += sum / XY_MIN(Y[j].bw, Y[parentVMIndex].bw);
//                    }
//                    if (rt < Max){
//                        rt = Max;
//                    }
//                }
//                //执行时间
//                double ExeTime = Tasks[taskIndex].length / Y[j].ps;
//                double StartTime = 0;
//                double EndTime = 0;
//                //搜索可用时间片段
//                set<double>::iterator pre  = vatL[j].begin();
//                set<double>::iterator post = vatL[j].begin();
//                ++post;
//                //源代码有问题，只判断第一段gap，中间gap跳过去了
//                if((*post - *pre) >= ExeTime && rt <= (*post)-ExeTime) {
//                    StartTime = XY_MAX(*pre, rt);
//                } else {
//                    pre = vatL[j].end();
//                    pre--; pre--;
//                    StartTime = XY_MAX(*pre, rt);
//                }
//                EndTime = StartTime + ExeTime;
//                //找最快完成
//                if (EndTime < finalEndTime1) {
//                    finalStartTime1 = StartTime;
//                    finalEndTime1 = EndTime;
//                    vmIndex1 = j;
//                }
//                //找在不改变成本的情况下，执行时间最小的那个虚拟机
//                if(EndTime < finalEndTime2 && ceil(EndTime) <= InsTimeLine[j]){
//                    finalStartTime2 = StartTime;
//                    finalEndTime2 = EndTime;
//                    vmIndex2 = j;
//                }
//            } else
//                continue;
//        }
//        if(vmIndex1 != vmIndex2){
//            //在Y中找是否有与vmIndex1类型相同且未使用过的的虚拟机
//            int Unused_SameType = -1;
//            for(int k = 0; k < Y.size(); ++k){
//                if(Y[vmIndex1].id == Y[k].id && flag[k] == 0){
//                    Unused_SameType = k;
//                    /**
//                     * 虚拟机是从前向后用的(ps高的排在前面)，不会出现跳过一个未使用的虚拟机而使用后面未使用的虚拟机
//                     * 去除一个虚拟机时，在Y中去除掉的那个虚拟机位置后面的那些虚拟机不会有任务执行
//                     * 所以 Y中少了一个元素不用去更新染色体已分配实例的下标
//                     */
//                    break;
//                }
//            }
//            if(Unused_SameType != -1){
//                Y.erase(Y.begin() + Unused_SameType);           //删除第k个元素
//                InsTimeLine.erase(InsTimeLine.begin() + Unused_SameType);   //删除一个实例
//                InsTimeLine[vmIndex1] += 3600;                  //删除了一个相同类型的虚拟机实例，则给此虚拟机实例增加利用时间
//            } else {
//                if(flag[vmIndex2] == 1){                        //I*min != null   Assign ti to I*min
//                    chrom.StartTime[taskIndex] = finalStartTime2;
//                    chrom.EndTime[taskIndex] = finalEndTime2;
//                    chrom.InsAlcLst[taskIndex] = vmIndex2;
//
//                    //在虚拟机VMIndex的时间列表段中查找taskIndex的开始时间
//                    //如果有，则不插入，并在vatL[VMIndex]删除该元素；如果没有，则插入
//                    if(vatL[vmIndex2].find(finalStartTime2) != vatL[vmIndex2].end()) {
//                        vatL[vmIndex2].erase(finalStartTime2);
//                    } else {
//                        vatL[vmIndex2].insert(finalStartTime2);
//                    }
//
//                    if(vatL[vmIndex2].find(finalEndTime2) != vatL[vmIndex2].end()) {
//                        vatL[vmIndex2].erase(finalEndTime2);
//                    } else {
//                        vatL[vmIndex2].insert(finalEndTime2);
//                    }
//                    flag[vmIndex2] = 1;
//                    max_EndTi = XY_MAX(max_EndTi, finalEndTime2);
//                }
//            }
//        }
//        if(chrom.InsAlcLst[taskIndex] == -1){               //if ti is not assign  未分配此任务实例为-1
//            chrom.StartTime[taskIndex] = finalStartTime1;
//            chrom.EndTime[taskIndex] = finalEndTime1;
//            chrom.InsAlcLst[taskIndex] = vmIndex1;
//            chrom.TypeLst[vmIndex1] = Y[vmIndex1].id - 1;
//
//            //在虚拟机VMIndex的时间列表段中查找taskIndex的开始时间
//            //如果有，则不插入，并在vatL[VMIndex]删除该元素；如果没有，则插入
//            if(vatL[vmIndex1].find(finalStartTime1) != vatL[vmIndex1].end()) {
//                vatL[vmIndex1].erase(finalStartTime1);
//            } else {
//                vatL[vmIndex1].insert(finalStartTime1);
//            }
//
//            if(vatL[vmIndex1].find(finalEndTime1) != vatL[vmIndex1].end()) {
//                vatL[vmIndex1].erase(finalEndTime1);
//            } else {
//                vatL[vmIndex1].insert(finalEndTime1);
//            }
//            flag[vmIndex1] = 1;
//            max_EndTi = XY_MAX(max_EndTi, finalEndTime1);
//        }
//    }
//    //一维变二维
////    GnrTwoDimCoding(chrom);
//    //工期
//    chrom.MakeSpan = max_EndTi;
//
//    //计算成本
//    double cost =0;
//    for (int k = 0; k <= *max_element(chrom.InsAlcLst.begin(),chrom.InsAlcLst.end()) ; ++k){
//        double Hrt = 10000000000;
//        double Rnt = 0;
//        vector<int> HandleTask;
//
//        for (int i = 0; i < comConst.NumOfTasks; ++i) {         //找到相同实例的任务
//            if(chrom.InsAlcLst[i] == k) HandleTask.push_back(i);//里面存放是的任务
//        }
//        if(HandleTask.size() == 0){                             //实例不连续情况
//            continue;
//        } else {
//            for (int i = 0; i < HandleTask.size(); ++i) {       //虚拟机开始租用时间的计算（Hrt）
//                //Hrt为任务的开始时间 减去 启动时间(仅考虑任务的启动时间) - 任务传输时间（传输时间占用虚拟机）
//                double hrt = chrom.StartTime[HandleTask[i]] - VMTypes[chrom.TypeLst[k]].RootTime;
//                if (hrt < Hrt) Hrt = hrt;
//            }
//            for (int i = 0; i < HandleTask.size(); ++i) {       //虚拟机归还时间计算（Rnt）
//                //Rnt时间的计算是该任务完成时间
//                if (Rnt < chrom.EndTime[HandleTask[i]])
//                    Rnt = chrom.EndTime[HandleTask[i]];
//            }
//            double Lt = Rnt - Hrt;                              //最晚归还时间减去最早开始时间为虚拟机的租用时间
//            double SpanTime = XY_MAX(Lt, VMTypes[chrom.TypeLst[k]].MinLeaseTime);
//            cost += VMTypes[chrom.TypeLst[k]].FixedCost + (ceil((SpanTime) / VMTypes[chrom.TypeLst[k]].ChargeUnit)) * VMTypes[chrom.TypeLst[k]].UnitCost;
//        }
//    }
//    chrom.Cost = cost;
//}

void GenericParticleByRank(chromosome& chrom, vector<double>& Rank, vector<double>& tskLFT) {
    //第一段虚拟机类型  ---  在此加条件 --- 资源可得不可得情况，TypePart[]下标为任务
    default_random_engine e(time(0));
    uniform_real_distribution<double> u(0,comConst.NumOfVMTypes - 0.51); //生成[0,4.49]之间的随机小数
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        chrom.TypeAlcPart[i] = u(e);
        //如果更新后的例子不是虚拟机可得的
        int Rsc = round(chrom.TypeAlcPart[i]);
        if(Tasks[i].ElgRsc.find(Rsc) == Tasks[i].ElgRsc.end()){
            //取一个可得类型中最近的类型编号
            if(Tasks[i].ElgRsc.size() == 1){
                chrom.TypeAlcPart[i] = chrom.TypeAlcPart[i] - Rsc + *(Tasks[i].ElgRsc.begin());
            } else {
                //查找最近的类型编号--此处有问题，没有考虑小于第一个或大于最后一个的情况！已需要修改-xy 2023-1-6
                //int TemRsc = SearchRsc(i, Rsc);
                //particle.TypePart[i] = particle.TypePart[i] + TemRsc - Rsc;
                set<int>::iterator itEnd=Tasks[i].ElgRsc.end(); --itEnd;
                if (Rsc <= *(Tasks[i].ElgRsc.begin()) ) {
                    chrom.TypeAlcPart[i] = *(Tasks[i].ElgRsc.begin());
                } else if (Rsc >= *itEnd){
                    chrom.TypeAlcPart[i] = *itEnd;
                } else {
                    for (set<int>::iterator it=Tasks[i].ElgRsc.begin(); it!=itEnd; ++it ){
                        set<int>::iterator itNxt = it; ++itNxt;
                        if (*it < Rsc && Rsc < *itNxt ) {
                            if ( *itNxt - Rsc < Rsc - *it ) {
                                chrom.TypeAlcPart[i] = *itNxt;
                            } else {
                                chrom.TypeAlcPart[i] = *it;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }

    //第二段调度顺序
    chrom.TskSchPart = Rank;
    DecodeParticle(chrom,tskLFT);
}

void GnrTypeAlcPartByRnd(chromosome& chrom) {
    //*为每个任务随机生成一个类型*//
    for(int i = 0; i < comConst.NumOfTasks; ++i){
//        //对应HPSO原文中虚拟机类型编码/取值范围为[1,m]即[0,m-1]，此方式要求VMType.txt中的虚拟机类型从0开始编号且连续-xy
//        chrom.TypePart[i] = rand()%(comConst.NumOfVMTypes-1) + rand()%1001/1000.0;
//        //对应虚拟机类型编码/取值范围为(-0.5,m-1+0.5)，此方式要求VMType.txt中的虚拟机类型从0开始编号且连续-xy
//        chrom.TypePart[i] = rand()%comConst.NumOfVMTypes - 0.49999 + rand()%1000/1000.0;
        int rnd = rand()%Tasks[i].ElgRsc.size();
        set<int>::iterator it = Tasks[i].ElgRsc.begin();
        advance(it,rnd);
        chrom.TypeAlcPart[i] = *it + rand()%1000/1000.0 - 0.49999;
    }
}

//修复粒子调度顺序 -- priority为整数的情况
void RepairPriorityAndGnrSchOrd_gui(chromosome& chrom) {
    vector<int> V, Q;       //里面存放的是任务
    vector<double> N(comConst.NumOfTasks, -1);
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        N[i] = round(chrom.TskSchPart[i]);
    }
    vector<int> upr(comConst.NumOfTasks,0);
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i] == 0) {
            Q.push_back(i);
        }
    }
    while (V.size() != comConst.NumOfTasks){
        //记录V中所有任务最大的Priority
        int MaxV = 0;
        for(int i = 0; i < V.size(); ++i){
            if(MaxV < N[V[i]])
                MaxV = N[V[i]];
        }
        for(int i = 0; i < Q.size(); ++i){   //遍历Q中的每一个任务
            int Task = Q[i];
            //找TaskIndex父任务最大的Priotity
            int MaxP = 0;
            for(int prnTsk: Tasks[Task].parents){
                if(MaxP < N[prnTsk]){
                    MaxP = N[prnTsk];
                }
            }
            if(N[Task] <= MaxP && Tasks[Task].parents.size() != 0){
                N[Task] = MaxP + 1;
            }
            for(int i1 = 0; i1 < V.size(); ++i1){
                if(N[Task] == N[V[i1]]){
                    N[Task] = MaxV + 1;
                    MaxV += 1;
                    break;
                }
            }
            V.push_back(Task);
        }
        //找出所有父任务都已被调度的任务添加到Q中
        vector<int> TemQ;
        for(int i = 0; i < Q.size(); ++i){
            int task = Q[i];
            for (int chlTsk: Tasks[task].children) {
                upr[chlTsk] = upr[chlTsk] - 1;
                if (upr[chlTsk] == 0) {
                    TemQ.push_back(chlTsk);
                }
            }
        }
        Q = TemQ;
    }
    //基于粒子得到一个可行的调度顺序 -- particle.TskSchLst即为调度顺序
    IndexSortByValueAsc(chrom.TskSchLst, N);
}

void RepairPriorityAndGnrSchOrd(chromosome& chrom) {
    vector<int> V, Q;
    vector<int> N(comConst.NumOfTasks, -1);
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        N[i] = round(chrom.TskSchPart[i]);
    }
    vector<int> upr(comConst.NumOfTasks, 0);
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0) {
            Q.push_back(i);
        }
    }
    int MaxV = -1;
    while (V.size() != comConst.NumOfTasks) {
        for (int TaskId: Q) {
            int MaxP = -1;
            for (int parent: Tasks[TaskId].parents) { //获得父任务中的最大值
                if (MaxP < N[parent])  MaxP = N[parent];
            }
            if (N[TaskId] <= MaxP) {
                N[TaskId] = MaxP + 1;
            }
            for (int tsk: V) {
                if (N[TaskId] == N[tsk])  {
                    N[TaskId] = MaxV + 1;
                    MaxV += 1;
                    break;
                }
            }
            MaxV = XY_MAX(N[TaskId], MaxV);
            V.push_back(TaskId);
        }
        vector<int> TemQ;
        for (int taskId: Q) {
            for (int childId: Tasks[taskId].children) {
                upr[childId] = upr[childId] - 1;
                if (upr[childId] == 0) {
                    TemQ.push_back(childId);
                }
            }
        }
        Q = TemQ;
    }
//    chrom.TskSchPart = N;
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        chrom.TskSchPart[i] = N[i];
    }
    IndexSortByValueOnAsc(chrom.TskSchLst, N);
}

void UpdateParticle(chromosome &ch,chromosome &Pbest, chromosome &Gbest, double &runtime, double &SchTime){
//    inputParameter.InertiaWeight = 0.09 * (1-(runtime / SchTime)) + 0.01; //HPSO原文设置为：0.1-0.01, 用gui的Pbest设置更新方法效果好
    inputParameter.InertiaWeight = 0.1 * (1-(runtime / SchTime)) + 0.9;   //按正确的Pbest设置更新方法，设置为1-0.9效果好-xy
    inputParameter.c1 = 2.0 * (1-(runtime / SchTime));                    //HPSO原文设置为：2.0-0
    inputParameter.c2 = 2.0 * (runtime / SchTime);                        //HPSO原文设置为：0-2.0，而不是2.0-0,已修正-xy
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        double r1 = RandomDouble(0,1), r2 = RandomDouble(0,1);
        ch.VTskSchPart[i] = inputParameter.InertiaWeight * ch.VTskSchPart[i] + inputParameter.c1 * r1 * (Pbest.TskSchPart[i] - ch.TskSchPart[i])
                            + inputParameter.c2 * r2 * (Gbest.TskSchPart[i] - ch.TskSchPart[i]);
        ch.TskSchPart[i] += ch.VTskSchPart[i];
        ch.VTypeAlcPart[i] = inputParameter.InertiaWeight * ch.VTypeAlcPart[i] + inputParameter.c1 * r1 * (Pbest.TypeAlcPart[i] - ch.TypeAlcPart[i])
                             + inputParameter.c2 * r2 * (Gbest.TypeAlcPart[i] - ch.TypeAlcPart[i]);
        ch.TypeAlcPart[i] += ch.VTypeAlcPart[i];
    }
    RepairTypeAlcPart(ch); //修复虚拟机类型分配粒子
}

void UpdateParticle_gui(chromosome &chrom, chromosome &Pbest, chromosome &Gbest, int t){
    inputParameter.InertiaWeight = 0.09 * (1-((double)t / inputParameter.NumOfIterations)) + 0.01;
    inputParameter.c1 = 2.0 * (1-((double)t / inputParameter.NumOfIterations));
    inputParameter.c2 = 2.0 * (1-((double)t / inputParameter.NumOfIterations));
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        int task = chrom.TskSchLst[i];
        double r1 = RandomDouble(0,1);
        double r2 = RandomDouble(0,1);

        chrom.VTskSchPart[task] = inputParameter.InertiaWeight * chrom.VTskSchPart[task] + inputParameter.c1 * r1 * (Pbest.TskSchPart[task] - chrom.TskSchPart[task])
                                  + inputParameter.c2 * r2 * (Gbest.TskSchPart[task] - chrom.TskSchPart[task]);
        chrom.TskSchPart[task] += chrom.VTskSchPart[task];

        chrom.VTypeAlcPart[task] = inputParameter.InertiaWeight * chrom.VTypeAlcPart[task] + inputParameter.c1 * r1 * (Pbest.TypeAlcPart[task] - chrom.TypeAlcPart[task])
                                  + inputParameter.c2 * r2 * (Gbest.TypeAlcPart[task] - chrom.TypeAlcPart[task]);
        chrom.TypeAlcPart[task] += chrom.VTypeAlcPart[task];
        //类型的边界最大和最小为[0~comConst.NumOfRsc]  如果超出了边界，则令其为第一个可得虚拟机类型
        int Rsc = round(chrom.TypeAlcPart[task]);
//        if(Rsc < 0) particle.TypePart[task] = Tasks[task].ElgRsc[0];
//        if(Rsc > comConst.NumOfRsc - 1) particle.TypePart[task] = Tasks[task].ElgRsc[Tasks[task].ElgRsc.size()-1];
        if ( Rsc<*(Tasks[task].ElgRsc.begin()) ) chrom.TypeAlcPart[task] = *(Tasks[task].ElgRsc.begin());
        if ( Rsc>*(Tasks[task].ElgRsc.end()) ) chrom.TypeAlcPart[task] = *(Tasks[task].ElgRsc.end());

        Rsc = round(chrom.TypeAlcPart[task]);
        //如果更新后的例子不是虚拟机可得的
        if(find(Tasks[task].ElgRsc.begin(), Tasks[task].ElgRsc.end(), Rsc) == Tasks[task].ElgRsc.end()){
            //取一个可得类型中最近的类型编号
            if(Tasks[task].ElgRsc.size() == 1){
                chrom.TypeAlcPart[task] = chrom.TypeAlcPart[task] + *(Tasks[task].ElgRsc.begin()) - Rsc;
            } else {
                //查找最近的类型编号
//                int TemRsc = SearchRsc(task, Rsc);
//                particle.TypePart[task] = particle.TypePart[task] + TemRsc - Rsc;
                set<int>::iterator itEnd=Tasks[task].ElgRsc.end(); --itEnd;
                if (Rsc <= *(Tasks[task].ElgRsc.begin()) ) {
                    chrom.TypeAlcPart[task] = *(Tasks[task].ElgRsc.begin());
                } else if (Rsc >= *itEnd){
                    chrom.TypeAlcPart[task] = *itEnd;
                } else {
                    for (set<int>::iterator it=Tasks[task].ElgRsc.begin(); it!=itEnd; ++it ){
                        set<int>::iterator itNxt = it; ++itNxt;
                        if (*it < Rsc && Rsc < *itNxt ) {
                            if ( *itNxt - Rsc < Rsc - *it ) {
                                chrom.TypeAlcPart[task] = *itNxt;
                            } else {
                                chrom.TypeAlcPart[task] = *it;
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
}

void RepairTypeAlcPart(chromosome& ch) {
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        int RscId = round(ch.TypeAlcPart[i]);
        if(RscId < *(Tasks[i].ElgRsc.begin())) { //超出下限的处理
            ch.TypeAlcPart[i] = *(Tasks[i].ElgRsc.begin());
            continue;
        }
        set<int>::iterator it = Tasks[i].ElgRsc.end(); --it;
        if(RscId > *it) {                       //超出上限的处理
            ch.TypeAlcPart[i] = *it;
            continue;
        }
        if(Tasks[i].ElgRsc.find(RscId) == Tasks[i].ElgRsc.end()){ //不存在的处理
            if(Tasks[i].ElgRsc.size() == 1) {
                ch.TypeAlcPart[i] = *(Tasks[i].ElgRsc.begin());
            } else {
                int TemRscId = FindNearestRscId(i, ch.TypeAlcPart[i]);
                ch.TypeAlcPart[i] = TemRscId;
            }
            continue;
        }
    }
}

int FindNearestRscId(int TaskId, double value) {
    set<int>::iterator lstIt = Tasks[TaskId].ElgRsc.end(); --lstIt;
    for (set<int>::iterator it = Tasks[TaskId].ElgRsc.begin(); it != lstIt; ++it ){
        set<int>::iterator nxtIt = it; ++nxtIt;
        if (*it < value && value < *nxtIt ) {
            if ( *nxtIt - value < value - *it ) {
                return *nxtIt;
            } else {
                return *it;
            }
        }
    }
}

void DecodeParticle(chromosome& chrom,vector<double>& tskLFT) {
    RepairPriorityAndGnrSchOrd(chrom); //修复并生成任务调度顺序
    int INS = 0;
    chrom.Cost = 0; chrom.MakeSpan = 0;
    vector<set<double>> ITL;
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        set<int> X, Y; //存放父任务执行的虚拟机实例，用set去重，X:serial instances，Y:parallel instances;
        int Task = chrom.TskSchLst[i];                 //按调度顺序取出一个任务Task
        int InsType = round(chrom.TypeAlcPart[Task]);  //获得处理任务Task的虚拟机类型；有问题-xy！ “TemTypeLst[i]”应该为TemTypeLst[Task];
        for(int i1 = 0; i1 < i; ++i1){ //遍历所有已调度（分配虚拟机实例）的任务
            if(chrom.TypeLst[chrom.InsAlcLst[chrom.TskSchLst[i1]]] == InsType){   //如果执行任务的虚拟机类型和Task所要求的类型相同，有问题-xy，“TemTypeLst[i1]”应该为“TemTypeLst[chrom.TskSchLst[i1]]”
                //如果父任务在此实例中，则把该任务的父任务执行的虚拟机添加到集合X中，否则添加到集合Y中
                if (find(Tasks[Task].parents.begin(), Tasks[Task].parents.end(), chrom.TskSchLst[i1]) != Tasks[Task].parents.end()){  //后续判断如果父任务集合为空，是否报错-不会-xy
                    X.insert(chrom.InsAlcLst[chrom.TskSchLst[i1]]);
                } else {
                    Y.insert(chrom.InsAlcLst[chrom.TskSchLst[i1]]);
                }
            }
        }
        map<int,pair<double,double>> H;
        for(int ins: X){ //从集合X中找出可以在截止时间内完成任务Task的实例
            //X中的实例类型与任务Task要求的类型是相同的（在生成X,Y时已有此判断），因此无需再加判断-xy
            //判断是否在截止日期(LFT)内完成，如果是，则添加到H集合中
            double DurTime = Tasks[Task].length / VMTypes[InsType].ps;
            double rt = FindTaskReadyTime(Task, ins, chrom);                  //计算任务就绪时间
            double temTaskStrTime = FindIdleTimeSlot(ITL[ins], DurTime, rt);  //搜索可用时间片段并计算开始时间
            double temTaskEndTime = temTaskStrTime + DurTime;                 //计算任务的完成时间
            if(temTaskEndTime - PrecisionValue < tskLFT[Task]){         //在截止日期(LFT)内完成则添加到H集合中
                pair<double,double> TskStrAndEndTime(temTaskStrTime,temTaskEndTime);
                H[ins] = TskStrAndEndTime;
            }
        }
        if(H.empty()){  //如果集合H为空，继续从集合Y中找
            for(int ins: Y){
                //Y中的实例类型与任务Task要求的类型是相同的（在生成X,Y时已有此判断），因此无需再加判断-xy
                double DurTime = Tasks[Task].length / VMTypes[InsType].ps;
                double rt = FindTaskReadyTime(Task, ins, chrom);                 //计算任务就绪时间
                double temTaskStrTime = FindIdleTimeSlot(ITL[ins], DurTime, rt); //搜索可用时间片段并计算开始时间
                double temTaskEndTime = temTaskStrTime + DurTime;                //计算任务完成时间
                if(temTaskEndTime - PrecisionValue < tskLFT[Task]){        //在截止日期(LFT)内完成则添加到H集合中
                    pair<double,double> TskStrAndEndTime(temTaskStrTime,temTaskEndTime);
                    H[ins] = TskStrAndEndTime;
                }
            }
        }
        int RscIndex = -1;
        double TaskStartTime, TaskEndTime;
        if(!H.empty()){  //如果H不为空，则选择一个执行成本增加最小的虚拟机实例
            double minAddedCost = InfInt;  //保存最小增加成本
            for(auto elm: H){
                int ins = elm.first;
                double temTaskStartTime = elm.second.first, temTaskEndTime = elm.second.second;
                double Rtn = 0, Hrt = InfInt; //虚拟机的归还时间和开始使用时间-还没有分配编号为k的资源实例的开始租用时间
                vector<int> handleTask;
                for (int n = 0; n < i; ++n) {
                    if (chrom.InsAlcLst[chrom.TskSchLst[n]] == ins)
                        handleTask.push_back(chrom.TskSchLst[n]);
                }
                for (int tsk: handleTask) {   //Rnt的计算 归还时间
                    if (Rtn + PrecisionValue < chrom.EndTime[tsk]) Rtn = chrom.EndTime[tsk];
                }
                for (int tsk: handleTask) {   //Hrt的计算 开始租用时间
                    double hrt = chrom.StartTime[tsk] - VMTypes[chrom.TypeLst[ins]].BootTime;
                    if (hrt + PrecisionValue < Hrt) Hrt = hrt;
                }

//                double frtTime1, lstTime1, LeaseTime1;      //没有处理任务前虚拟机的开始使用时间，结束使用时间，租用时间归还时间
//                double frtTime2, lstTime2, LeaseTime2;      //处理任务后虚拟机的开始使用时间，结束使用时间，租用时间归还时间
//                if (ITL[ins].size() > 2 || fabs( *(ITL[ins].begin()) - VMTypes[chrom.TypeLst[ins]].BootTime ) >= PrecisionValue) {
//                    set<double>::iterator pFrt = ITL[ins].begin(), pLst = ITL[ins].end();
//                    pLst--;  pLst--;
//                    frtTime1 = VMTypes[chrom.TypeLst[ins]].BootTime, lstTime1 = *pLst;
//                    if (fabs(VMTypes[chrom.TypeLst[ins]].BootTime - *pFrt) < PrecisionValue) {
//                        pFrt++;
//                        frtTime1 = *pFrt;
//                    }
//                }else {
//                    cout << endl << "a VM being not used is added! this is wrong!"; exit(5);
//                }
//                if (fabs(frtTime1-Hrt)>PrecisionValue || fabs(lstTime1-Rtn)>PrecisionValue) {
//                    cout<<endl<<"hirt or rent time is not equal ";
//                    exit(6);
//                }

                double bLT = Rtn - Hrt - PrecisionValue;  //分配该任务之前的租用时间
                double aLT = XY_MAX(temTaskEndTime, Rtn) - XY_MIN(temTaskStartTime - VMTypes[chrom.TypeLst[ins]].BootTime, Hrt) -PrecisionValue; //分配该任务之后的租用时间
                //计算增加的成本
                double TempAddCost = ceil(aLT / VMTypes[chrom.TypeLst[ins]].ChargeUnit) * VMTypes[chrom.TypeLst[ins]].UnitCost
                                     - ceil(bLT / VMTypes[chrom.TypeLst[ins]].ChargeUnit) * VMTypes[chrom.TypeLst[ins]].UnitCost;
                if(minAddedCost - PrecisionValue > TempAddCost){
                    minAddedCost = TempAddCost;
                    TaskStartTime = temTaskStartTime;  //这里有问题,缺少对应的TaskStartTime和TaskEndTime的保存-xy, 已修改
                    TaskEndTime = temTaskEndTime;
                    RscIndex = ins;
                }
            }
            chrom.Cost += minAddedCost;   //基于增量计算成本
        } else { // /*H为空，添加一个新的实例*/
            RscIndex = INS;  INS++;       //新增了一个实例，INS值加1
            chrom.TypeLst[RscIndex] = InsType;
            set<double> a;
            a.insert(VMTypes[InsType].BootTime);
            a.insert(InfInt * 1.0);
            ITL.push_back(a);
            double rt = 0;
            for (int prn: Tasks[Task].parents) { //计算任务taskIndex就绪时间
                int prnInsId = chrom.InsAlcLst[prn];
                double fft = chrom.EndTime[prn] + ParChildTranFileSizeSum[prn][Task] / (XY_MIN(VMTypes[InsType].bw, VMTypes[chrom.TypeLst[prnInsId]].bw));
                if (rt < fft)  rt = fft;
            }
            double DurTime = Tasks[Task].length / VMTypes[InsType].ps;
            TaskStartTime = XY_MAX(VMTypes[InsType].BootTime, rt);
            TaskEndTime = TaskStartTime + DurTime;
            double lt = XY_MAX(DurTime + VMTypes[InsType].BootTime, VMTypes[InsType].MinLeaseTime) - PrecisionValue;
            chrom.Cost += VMTypes[InsType].FixedCost + VMTypes[InsType].UnitCost * ceil( lt / VMTypes[InsType].ChargeUnit);  //有问题-xy,已修改，
        }
        chrom.StartTime[Task] = TaskStartTime;
        chrom.EndTime[Task] = TaskEndTime;
        chrom.InsAlcLst[Task] = RscIndex;
        chrom.MakeSpan = XY_MAX(chrom.MakeSpan, chrom.EndTime[Task]);
        UpdateITL(ITL[RscIndex], TaskStartTime, TaskEndTime);
    }
    /*测试用*/
    chromosome temch = chrom;
    DecodeChrom(temch);
    if (fabs(temch.Cost - chrom.Cost) > PrecisionValue || fabs(temch.MakeSpan - chrom.MakeSpan) > PrecisionValue) {
        cout <<endl<< "gui: "<< chrom.MakeSpan <<"; "<<chrom.Cost;
        cout <<endl<< "ok: "<< temch.MakeSpan <<"; "<<temch.Cost;
        exit(7);
    }
}

//void DecodeParticle_gui(chromosome& chrom) {
//    RepairPriorityAndGnrSchOrd(chrom); //修复并生成任务调度顺序
//    vector<double> TemTypeLst(comConst.NumOfTasks, 0);
//    for(int i = 0; i < comConst.NumOfTasks; ++i){ //获得虚拟机类型
//        TemTypeLst[i] = round(chrom.TypeAlcPart[i]);
//    }
//    int INS = 0;
//    chrom.Cost = 0;
//    vector<set<double>> vatL;
//    for(int i = 0; i < comConst.NumOfTasks; ++i){
//        int RscIndex = -1;
//        double TaskStartTime;
//        double TaskEndTime;
//        set<int> X, Y; //存放父任务执行的虚拟机实例，用set去重，X:serial instances，Y:parallel instances;
//        vector<int> H;
//        //取出该任务及其类型
//        int Task = chrom.TskSchLst[i];
//        int InsType = TemTypeLst[i];
//        for(int i1 = 0; i1 < i; ++i1){
//            if(TemTypeLst[i1] == InsType){   //如果遍历的任务执行的虚拟机类型和该任务所要求的类型相同
//                //如果父任务在此实例中，则把该任务的父任务执行的虚拟机添加到集合X中，否则添加到集合Y中
//                if (find(Tasks[Task].parents.begin(), Tasks[Task].parents.end(), chrom.TskSchLst[i1]) != Tasks[Task].parents.end()){  //后续判断如果父任务集合为空，是否报错
//                    X.insert(chrom.InsAlcLst[chrom.TskSchLst[i1]]);
//                } else {
//                    Y.insert(chrom.InsAlcLst[chrom.TskSchLst[i1]]);
//                }
//            }
//        }
//        //从集合X中找出任务Task可以使用的实例（类型可用且在截止时间内完成）
//        for(int ins: X){
//            //如果该实例的类型是任务Task可得的类型 && 执行的时间小于deadline
//            if(Tasks[Task].ElgRsc.find(chrom.TypeLst[ins]) != Tasks[Task].ElgRsc.end()){
//                //判断是否在截止日期内完成，如果是，则添加到H集合中
//                double DurTime = Tasks[Task].length / VMTypes[InsType].ps;
//                double rt = FindTaskReadyTime(Task, ins, chrom);    //计算任务taskIndex就绪时间
//                TaskStartTime = FindIdleTimeSlot(vatL[ins], DurTime, rt);    //搜索可用时间片段
//                TaskEndTime = TaskStartTime + DurTime;  //该任务的完成时间
//                if(TaskEndTime - PrecisionValue < TaskDeadline[Task]){
//                    H.push_back(ins);
//                }
//            }
//        }
//        if(H.empty()){  //如果集合H为空，继续从集合Y中找
//            for(int ins: Y){
//                //如果该实例的类型是任务Task可得的类型 && 执行的时间小于deadline
//                if(Tasks[Task].ElgRsc.find(chrom.TypeLst[ins]) != Tasks[Task].ElgRsc.end()){
//                    //判断是否在截止日期内完成，如果是，则添加到H集合中
//                    double DurTime = Tasks[Task].length / VMTypes[InsType].ps;
//                    double rt = FindTaskReadyTime(Task, ins, chrom);   //计算任务taskIndex就绪时间
//                    TaskStartTime = FindIdleTimeSlot(vatL[ins], DurTime, rt);   //搜索可用时间片段
//                    TaskEndTime = TaskStartTime + DurTime;                      //计算任务完成时间
//                    if(TaskEndTime - PrecisionValue < TaskDeadline[Task]){
//                        H.push_back(ins);
//                    }
//                }
//            }
//        }
//        double LeaseCost = InfInt;  //租用增加成本最小成本记录
//        if(!H.empty()){  //如果H不为空，则选择一个执行成本增加最小的虚拟机
//            for(int ins: H){
//                double DurTime = Tasks[Task].length / VMTypes[InsType].ps;
//                double rt = FindTaskReadyTime(Task, ins, chrom);    //计算任务taskIndex就绪时间
//                TaskStartTime = FindIdleTimeSlot(vatL[ins], DurTime, rt);    //搜索可用时间片段
//                TaskEndTime = TaskStartTime + DurTime;
//                double Rnt = 0, Hrt = InfInt; //虚拟机的归还时间和开始使用时间-还没有分配编号为k的资源实例的开始租用时间
//                vector<int> handleTask;
//                for (int n = 0; n < i; ++n) {
//                    if (chrom.InsAlcLst[chrom.TskSchLst[n]] == ins)
//                        handleTask.push_back(chrom.TskSchLst[n]);
//                }
//                for (int k = 0; k < handleTask.size(); ++k) {   //Rnt的计算 归还时间
//                    double time =  chrom.EndTime[handleTask[k]];
//                    if (Rnt < time) Rnt = time;
//                }
//                for (int l = 0; l < handleTask.size(); ++l) {   //Hrt的计算 开始租用时间
//                    double hrt = chrom.StartTime[handleTask[l]] - VMTypes[chrom.TypeLst[ins]].BootTime;
//                    if (hrt < Hrt) Hrt = hrt;
//                }
//
//                double Lt = Rnt - Hrt;  //分配该任务之前的租用时间
//                double Lt1 = XY_MAX(TaskEndTime, Rnt) - XY_MIN(TaskStartTime - VMTypes[chrom.TypeLst[ins]].BootTime, Hrt); //分配该任务之后的租用时间
//                //计算增加的成本
//                double TempAddCost = ceil(Lt1 / VMTypes[chrom.TypeLst[ins]].ChargeUnit) * VMTypes[chrom.TypeLst[ins]].UnitCost
//                                     - ceil(Lt / VMTypes[chrom.TypeLst[ins]].ChargeUnit) * VMTypes[chrom.TypeLst[ins]].UnitCost;
//                if(LeaseCost > TempAddCost){
//                    LeaseCost = TempAddCost;
//                    //这里有问题缺少对应的TaskStartTime和TaskEndTime的保存-xy
//                    RscIndex = ins;
//                }
//            }
//            chrom.Cost += LeaseCost; //着重
//        } else {
//            //添加一个新的实例
//            RscIndex = INS;
//            chrom.TypeLst[RscIndex] = InsType;
//            chrom.Cost += VMTypes[InsType].UnitCost;  //有问题-xy，
//            INS++;                                                      //新增了一个实例，INS值加1
//            set<double> a;
//            a.insert(VMTypes[InsType].BootTime);                           //虚拟机初始化
//            a.insert(InfInt * 1.0);
//            vatL.push_back(a);
//
//            //开始时间和结束时间的计算
//            double DurTime = Tasks[Task].length / VMTypes[InsType].ps;
//
//            //计算任务taskIndex就绪时间
//            double rt = 0;
//            for (int i = 0; i < Tasks[Task].parents.size(); i++) {
//                int parentIndex = Tasks[Task].parents[i];
//                int parentVMIndex = chrom.InsAlcLst[parentIndex];
//                double ParEndTime = chrom.EndTime[parentIndex];
//                double sum = ParChildTranFileSizeSum[parentIndex][Task];
//                ParEndTime += sum / (XY_MIN(VMTypes[InsType].bw, VMTypes[chrom.TypeLst[parentVMIndex]].bw));
//                if (rt < ParEndTime) {
//                    rt = ParEndTime;
//                }
//            }
//            TaskStartTime = FindIdleTimeSlot(vatL[RscIndex], DurTime, rt);    //搜索可用时间片段
//            TaskEndTime = TaskStartTime + DurTime;
//        }
//        chrom.StartTime[Task] = TaskStartTime;
//        chrom.EndTime[Task] = TaskEndTime;
//        chrom.InsAlcLst[Task] = RscIndex;
//        chrom.MakeSpan = XY_MAX(chrom.MakeSpan, chrom.EndTime[Task]);
//
//        UpdateITL(vatL[RscIndex], TaskStartTime, TaskEndTime);
//    }
//}

//calculate ready time
double FindTaskReadyTime(int tsk, int ins, chromosome &chrom){
    //计算任务taskIndex就绪时间
    double rt = 0;
    for (int prn: Tasks[tsk].parents) {
        int prnIns = chrom.InsAlcLst[prn];
        double PrnEndTime = chrom.EndTime[prn];
        if (ins != prnIns) {
            PrnEndTime += ParChildTranFileSizeSum[prn][tsk] / (XY_MIN(VMTypes[chrom.TypeLst[ins]].bw, VMTypes[chrom.TypeLst[prnIns]].bw));
        }
        if (rt < PrnEndTime) {
            rt = PrnEndTime;
        }
    }
    return rt;
}

void DecodeChrom(chromosome& ch){
    vector<set<double>> ITL;
    ch.MakeSpan = 0;  //个体有可能被多次使用，为确保正确，要初始化为0-xy
    int maxInsId = *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end());  //获得最大的实例编号
    if(ch.IsFrw){
        for (int j = 0; j <= maxInsId; ++j) {
            set<double> a;
            a.insert(VMTypes[ch.TypeLst[j]].BootTime);//正向解码时初始化为boottime；反向解码时初始为0，但在计算makespan时每个资源实例的最大完成时间上加boottime；-xy2 //-G 已改
            a.insert(InfInt * 1.0);
            ITL.push_back(a);
        }
    } else {
        for (int j = 0; j <= maxInsId; ++j) {
            set<double> a;
            a.insert(0.0);
            a.insert(InfInt * 1.0);
            ITL.push_back(a);
        }
    }
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        int crnTsk = ch.TskSchLst[i], crnIns = ch.InsAlcLst[crnTsk];
        double ReadyTime = 0;
        if(ch.IsFrw) { //forward-loading
            for (int prnTsk: Tasks[crnTsk].parents) {
                int prnIns = ch.InsAlcLst[prnTsk];
                double fft = ch.EndTime[prnTsk];
                if(crnIns != prnIns) {
                    fft = fft + ParChildTranFileSizeSum[prnTsk][crnTsk] / XY_MIN(VMTypes[ch.TypeLst[crnIns]].bw, VMTypes[ch.TypeLst[prnIns]].bw);
                }
                if (ReadyTime + PrecisionValue < fft) {
                    ReadyTime = fft;
                }
            }
        } else { //backward-loading
            for (int chlTsk: Tasks[crnTsk].children) {
                int chlIns = ch.InsAlcLst[chlTsk];
                double fft = ch.EndTime[chlTsk];
                if(crnIns != chlIns) {
                    fft = fft + ParChildTranFileSizeSum[crnTsk][chlTsk] / XY_MIN(VMTypes[ch.TypeLst[crnIns]].bw, VMTypes[ch.TypeLst[chlIns]].bw);
                }
                if (ReadyTime + PrecisionValue < fft) {
                    ReadyTime = fft;
                }
            }
        }
        double ExecutionTime = Tasks[crnTsk].length / VMTypes[ch.TypeLst[crnIns]].ps;
        ch.StartTime[crnTsk] = FindIdleTimeSlot(ITL[crnIns], ExecutionTime, ReadyTime);
        ch.EndTime[crnTsk] = ch.StartTime[crnTsk] + ExecutionTime;
        if(ch.IsFrw){
            if (ch.MakeSpan + PrecisionValue < ch.EndTime[crnTsk]) {
                ch.MakeSpan = ch.EndTime[crnTsk];
            }
        } else {
            if (ch.MakeSpan + PrecisionValue < ch.EndTime[crnTsk] + VMTypes[ch.TypeLst[crnIns]].BootTime) {
                ch.MakeSpan = ch.EndTime[crnTsk] + VMTypes[ch.TypeLst[crnIns]].BootTime;
            }
        }
        UpdateITL(ITL[crnIns], ch.StartTime[crnTsk], ch.EndTime[crnTsk]);  //{update ITL}
    }
    //计算成本，根据ITL计算
    double cost1 = 0, cost2 = 0;
    for (int j = 0; j < ITL.size(); ++j) {
        double bt = 0;
        if (ch.IsFrw)
            bt = VMTypes[ch.TypeLst[j]].BootTime;
        if (ITL[j].size() > 2 || fabs( *(ITL[j].begin()) - bt ) >= PrecisionValue) { //有2个时间点以上或第1个时间点不等于boottime，那么虚拟机必定被使用了
            set<double>::iterator pFrt = ITL[j].begin(), pLst = ITL[j].end();
            pLst--;  pLst--;
            double frtTime = VMTypes[ch.TypeLst[j]].BootTime, lstTime = *pLst;
            if (fabs(bt - *pFrt) < PrecisionValue) {//第1个点等于BootTime，那么第二个点才是才是真正开始执行任务的起点
                pFrt++;
                frtTime = *pFrt;
            }
            double LeaseTime = XY_MAX(lstTime - frtTime + VMTypes[ch.TypeLst[j]].BootTime, VMTypes[ch.TypeLst[j]].MinLeaseTime);
            cost1 = cost1 + VMTypes[ch.TypeLst[j]].FixedCost + VMTypes[ch.TypeLst[j]].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[ch.TypeLst[j]].ChargeUnit);
        }
    }
    //计算成本，根据任务的开始和完成时间计算
    vector<vector<pair<double, double>>> TskLstOfInstance( maxInsId + 1 );
    for (int i = 0; i < comConst.NumOfTasks; ++i) { //把任务及其开始时间添加到各个实例的任务列表中去
        pair<double, double> tm(ch.StartTime[i], ch.EndTime[i]);
        TskLstOfInstance[ch.InsAlcLst[i]].push_back(tm);
    }
    for ( int k = 0; k <= maxInsId; ++k ) {//实例任务列表中的元素按开始时间从小到大进行排序并计算成本
        if (!TskLstOfInstance[k].empty()) {
            sort(TskLstOfInstance[k].begin(), TskLstOfInstance[k].end(), SortValueByAscend);
            vector<pair<double, double>>::iterator pFrt=TskLstOfInstance[k].begin(), pLst = TskLstOfInstance[k].end();
            pLst--;
            double LeaseTime = XY_MAX((*pLst).second - (*pFrt).first + VMTypes[ch.TypeLst[k]].BootTime, VMTypes[ch.TypeLst[k]].MinLeaseTime);
            cost2 = cost2 + VMTypes[ch.TypeLst[k]].FixedCost + VMTypes[ch.TypeLst[k]].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[ch.TypeLst[k]].ChargeUnit);
        }
    }

    if ( fabs(cost1 - cost2) >= PrecisionValue ) //测试用
        cout << endl << "The results of the two methods to calculate cost are inconsistent. cost1:  " << cost1 << "; cost2" <<cost2 << endl;

    ch.Cost = cost1;
}

void DecodeChrom_NI(chromosome& ch){
    if (ch.IsFrw == false) {
        cout<< endl<< "ch.IsFrw is false in DecodeChrom_NI"; exit(2);
    }
    ch.MakeSpan = 0;  //个体有可能被多次使用，为确保正确，要初始化为0-xy
    int maxInsId = *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end());  //获得最大的实例编号
    vector<double> insAvlTime(maxInsId+1, 0);
    vector<double> insStrTime(maxInsId+1, -1);
    if(ch.IsFrw){
        for (int j = 0; j <= maxInsId; ++j) {
            insAvlTime[j] = VMTypes[ch.TypeLst[j]].BootTime; //正向解码时初始化为boottime；反向解码时初始为0，但在计算makespan时每个资源实例的最大完成时间上加boottime；
        }
    }
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        int crnTsk = ch.TskSchLst[i], crnIns = ch.InsAlcLst[crnTsk];
        double ReadyTime = 0;
        if(ch.IsFrw) { //forward-loading
            for (int prnTsk: Tasks[crnTsk].parents) {
                int prnIns = ch.InsAlcLst[prnTsk];
                double fft = ch.EndTime[prnTsk];
                if(crnIns != prnIns) {
                    fft = fft + ParChildTranFileSizeSum[prnTsk][crnTsk] / XY_MIN(VMTypes[ch.TypeLst[crnIns]].bw, VMTypes[ch.TypeLst[prnIns]].bw);
                }
                if (ReadyTime + PrecisionValue < fft) {
                    ReadyTime = fft;
                }
            }
        } else { //backward-loading
            for (int chlTsk: Tasks[crnTsk].children) {
                int chlIns = ch.InsAlcLst[chlTsk];
                double fft = ch.EndTime[chlTsk];
                if(crnIns != chlIns) {
                    fft = fft + ParChildTranFileSizeSum[crnTsk][chlTsk] / XY_MIN(VMTypes[ch.TypeLst[crnIns]].bw, VMTypes[ch.TypeLst[chlIns]].bw);
                }
                if (ReadyTime + PrecisionValue < fft) {
                    ReadyTime = fft;
                }
            }
        }
        double ExecutionTime = Tasks[crnTsk].length / VMTypes[ch.TypeLst[crnIns]].ps;
        ch.StartTime[crnTsk] = XY_MAX(ReadyTime, insAvlTime[crnIns]);
        if (fabs(insStrTime[crnIns] + 1) < PrecisionValue) { //如果第一次使用crnIns，则记录其开始使用的时间
            insStrTime[crnIns] = ch.StartTime[crnTsk];
        }
        ch.EndTime[crnTsk] = ch.StartTime[crnTsk] + ExecutionTime;
        if(ch.IsFrw){
            if (ch.MakeSpan + PrecisionValue < ch.EndTime[crnTsk]) {
                ch.MakeSpan = ch.EndTime[crnTsk];
            }
        } else {
            if (ch.MakeSpan + PrecisionValue < ch.EndTime[crnTsk] + VMTypes[ch.TypeLst[crnIns]].BootTime) {
                ch.MakeSpan = ch.EndTime[crnTsk] + VMTypes[ch.TypeLst[crnIns]].BootTime;
            }
        }
        insAvlTime[crnIns] = ch.EndTime[crnTsk];  //update insAvlTime
    }
    //计算成本，根据ITL计算
    double cost1 = 0, cost2 = 0;
    for (int j = 0; j <= maxInsId; ++j) {
        if ( fabs(insStrTime[j] + 1) > PrecisionValue ) { //insStrTime[j]的值不为初始值-1，那么该实例必定被使用了
            double LeaseTime = XY_MAX(insAvlTime[j] - insStrTime[j] + VMTypes[ch.TypeLst[j]].BootTime, VMTypes[ch.TypeLst[j]].MinLeaseTime);
            cost1 = cost1 + VMTypes[ch.TypeLst[j]].FixedCost + VMTypes[ch.TypeLst[j]].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[ch.TypeLst[j]].ChargeUnit);
        }
    }
    //计算成本，根据任务的开始和完成时间计算
    vector<vector<pair<double, double>>> TskLstOfInstance( maxInsId + 1 );
    for (int i = 0; i < comConst.NumOfTasks; ++i) { //把任务及其开始时间添加到各个实例的任务列表中去
        pair<double, double> tm(ch.StartTime[i], ch.EndTime[i]);
        TskLstOfInstance[ch.InsAlcLst[i]].push_back(tm);
    }
    for ( int k = 0; k <= maxInsId; ++k ) {//实例任务列表中的元素按开始时间从小到大进行排序并计算成本
        if (!TskLstOfInstance[k].empty()) {
            sort(TskLstOfInstance[k].begin(), TskLstOfInstance[k].end(), SortValueByAscend);
            vector<pair<double, double>>::iterator pFrt=TskLstOfInstance[k].begin(), pLst = TskLstOfInstance[k].end();
            pLst--;
            double LeaseTime = XY_MAX((*pLst).second - (*pFrt).first + VMTypes[ch.TypeLst[k]].BootTime, VMTypes[ch.TypeLst[k]].MinLeaseTime);
            cost2 = cost2 + VMTypes[ch.TypeLst[k]].FixedCost + VMTypes[ch.TypeLst[k]].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[ch.TypeLst[k]].ChargeUnit);
        }
    }

    if ( fabs(cost1 - cost2) >= PrecisionValue ) //测试用
        cout << endl << "The results of the two methods to calculate cost are inconsistent. cost1: " << cost1 << "; cost2: " <<cost2 << endl;

    ch.Cost = cost1;
}

int IFBDI(chromosome& ch) {
    if (ch.IsIFBDI) return 0;
    chromosome NewChrom = ch;
    chromosome OldChrom;
    int nmb = 0;
    do {
        OldChrom = NewChrom;
        vector<int> ind(comConst.NumOfTasks);
        IndexSortByValueDesc(NewChrom.TskSchLst, OldChrom.EndTime);
        NewChrom.IsFrw = !(OldChrom.IsFrw);
        DecodeChrom(NewChrom);
        ++nmb;
    } while ( (NewChrom.Cost + PrecisionValue < OldChrom.Cost && NewChrom.MakeSpan - PrecisionValue < OldChrom.MakeSpan)||
              (NewChrom.Cost - PrecisionValue < OldChrom.Cost && NewChrom.MakeSpan + PrecisionValue < OldChrom.MakeSpan));
//    /* 取优的（可以是反向个体）*/
//    if (fabs(NewChrom.Cost - OldChrom.Cost) < PrecisionValue && fabs(NewChrom.MakeSpan - OldChrom.MakeSpan) < PrecisionValue && NewChrom.IsFrw)  //相等取正向个体
//        ch = NewChrom;
//    else  //不相等取小的
//        ch = OldChrom;
    /*取正向个体*/
//    if ( NewChrom.IsFrw)  //取正向个体
//        ch = NewChrom;
//    else
//        ch = OldChrom;
    if ( OldChrom.IsFrw == false)  {
        IndexSortByValueDesc(OldChrom.TskSchLst, OldChrom.EndTime);
        OldChrom.IsFrw = true;
        DecodeChrom(OldChrom);
    }
    ch = OldChrom;
    ch.IsIFBDI = true;
    return nmb;  //nmb>1表明有改进
}

int LBCAI(chromosome& ch) {
    if (ch.IsLBCAI)  return 0;
    chromosome OldCh = ch;
    int nmbOfInsUsed = *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end()) + 1;
    vector<double> load(nmbOfInsUsed, 0);
    //{calculate the loads of instances,find out the set TSK[j] of tasks allocated to instance j; }
    vector<vector<int>> TSK(nmbOfInsUsed);
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        load[ch.InsAlcLst[i]] += Tasks[i].length / VMTypes[ch.TypeLst[ch.InsAlcLst[i]]].ps;
        TSK[ch.InsAlcLst[i]].push_back(i);
    }
    vector<int> ind(nmbOfInsUsed);
    IndexSortByValueAsc(ind, load);     //sorting according to loads
    int insWithMinLd = ind[0];          //find out the instance with the lowest load;
    set<int> ST;
    if (fabs(load[insWithMinLd]) < PrecisionValue) {//测试用-xy
        cout<<endl<<"load is equal to 0, it is wrong in LBCAI !"; exit(2);
    }
    //traverse the tasks allocated to the instance with lowest load and add their parents and children to set ST
    for (int tsk: TSK[insWithMinLd]) {
        ST.insert(Tasks[tsk].children.begin(),Tasks[tsk].children.end());
        ST.insert(Tasks[tsk].parents.begin(),Tasks[tsk].parents.end());
    }
    //delete the tasks which have been allocated the instance with lowest load
    for (int tsk: TSK[insWithMinLd]) {
        ST.erase(tsk);
    }
    //delete the tasks which can not be performed by the instance with lowest load
    for (auto iter = ST.begin(); iter != ST.end();) {
        if (VMTypes[ch.TypeLst[insWithMinLd]].ElgTsk.find(*iter) ==  VMTypes[ch.TypeLst[insWithMinLd]].ElgTsk.end()){
            //iter = ST.erase(iter); //erase函数删除一个元素后返回的是下一个元素的指针！因此此句的作用和下一句是一样的-xy
            ST.erase(iter++);//等效于 set<int>::iterator temIt = iter; iter++; pop.erase(temIt); //要先把指针移到下一个元素，然后删除当前要删除的元素！！
        } else
            ++iter;
    }
    if(ST.empty()){
        ST.insert(VMTypes[ch.TypeLst[insWithMinLd]].ElgTsk.begin(), VMTypes[ch.TypeLst[insWithMinLd]].ElgTsk.end());
    }
    //Sort the tasks in ST according to the load of the resource to which the task is allocated
    vector<pair<int, double >> t;
    for (auto tsk: ST) {
        t.push_back(pair<int, double>(tsk, load[ch.InsAlcLst[tsk]]));
    }
    sort(t.begin(), t.end(), SortValueByDescend);
    ch.InsAlcLst[t[0].first] = insWithMinLd;
    RecodeInsAndType(ch);  //由于InsAlcLst发生了改变，因此需要重新编码
//    for (int j = comConst.NumOfTasks -1; j >=0 ; --j) {  //对没有使用的实例的类型进行随机赋值
//        if (ch.TypeLst[j] == -1) {
//            set<int>::iterator it = VMTypeSet.begin();
//            int rnd = rand()%VMTypeSet.size();
//            advance(it, rnd);
//            ch.TypeLst[j] = *it;
//        } else
//            break;
//    }
    DecodeChrom(ch);
    ch.IsIFBDI = false;
    IFBDI(ch);
    ch.IsLBCAI = true;
    if ((ch.Cost + PrecisionValue < OldCh.Cost && ch.MakeSpan - PrecisionValue < OldCh.MakeSpan)||
        (ch.Cost - PrecisionValue < OldCh.Cost && ch.MakeSpan + PrecisionValue < OldCh.MakeSpan)) {
        return 2;
    } else {
        ch = OldCh;
        return 1;
    }
}

int RPI_ShfTsk1(chromosome& ch){
    //实例中的首任务尽可能晚的执行！
    int maxInsId = *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end());
    vector<map<double,int>> tmTskLstOfIns( maxInsId + 1 );  //为了利用set结构的排序功能把StartTime[i]作为key，即：map[StartTime[i]] = i；
    for (int i = 0; i < comConst.NumOfTasks; ++i) { //把任务及其开始时间添加到各个实例的时间任务列表中去
        tmTskLstOfIns[ch.InsAlcLst[i]][ch.StartTime[i]] = {i};
    }
    double cst = 0;
    for ( int k = 0; k <= maxInsId; ++k ) {//遍历每个实例
        if (!tmTskLstOfIns[k].empty()) {
            map<double,int>::iterator it = tmTskLstOfIns[k].begin();
            int crnTsk = (*it).second;   ++it;
            if (it != tmTskLstOfIns[k].end() && (*it).first - PrecisionValue > ch.EndTime[crnTsk] ) { //第一个任务存在向后移动的可能性
                double slkTime = InfInt;
                for (int chlTsk: Tasks[crnTsk].children) {
                    double st = ch.StartTime[chlTsk] - ch.EndTime[crnTsk];
                    if (k != ch.InsAlcLst[chlTsk])
                        st -= ParChildTranFileSizeSum[crnTsk][chlTsk]/XY_MIN(VMTypes[ch.TypeLst[k]].bw, VMTypes[ch.TypeLst[ch.InsAlcLst[chlTsk]]].bw);
                    if (st + PrecisionValue < 0) { //测试用-xy
                        cout<<endl<<"slack time is less 0" ; exit(2);
                    }
                    if (st + PrecisionValue < slkTime)
                        slkTime = st;
                }
                if (slkTime > PrecisionValue) {
                    double ShiftTime = XY_MIN(slkTime, (*it).first - ch.EndTime[crnTsk]);
                    ch.StartTime[crnTsk] += ShiftTime;
                    ch.EndTime[crnTsk] += ShiftTime;
//                    cout << endl << "there is a slack time ! ";
                }
            }
            map<double,int>::iterator pLst = tmTskLstOfIns[k].end(); --pLst;
            int lstTsk = (*pLst).second;
            double LeaseTime = XY_MAX(ch.EndTime[lstTsk] - ch.StartTime[crnTsk] + VMTypes[ch.TypeLst[k]].BootTime, VMTypes[ch.TypeLst[k]].MinLeaseTime);
            cst = cst + VMTypes[ch.TypeLst[k]].FixedCost + VMTypes[ch.TypeLst[k]].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[ch.TypeLst[k]].ChargeUnit);
        }
    }
    if (cst + PrecisionValue < ch.Cost) {
        ch.Cost = cst;
//        cout << endl << "RPI_S is Ok ! ";
        return 1;
    }
    return 0;
}

int RPI_ShfTsk2(chromosome& ch){
    //除了实例中的最后一个任务，其余的任务都尽可能晚的执行！
//    if (ch.IsFrw == false) {
//        cout<< endl<< "ch.IsFrw is false in RscPrvImp_A"; exit(2);
//    }
    int maxInsId = *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end());
    vector<map<double,int>> tmTskLstOfIns( maxInsId + 1 ); //为利用set的排序功能把StartTime[i]作为key，即：map[StartTime[i]] = i；
    for (int i = 0; i < comConst.NumOfTasks; ++i) {        //把任务及其开始时间添加到各个实例的时间任务列表中去
        tmTskLstOfIns[ch.InsAlcLst[i]][ch.StartTime[i]] = {i};
    }
    set<int> lstTskOfIns;   //记录所有实例其要处理的最后一个任务；
    for ( int k = 0; k <= maxInsId; ++k ) {
        if (!tmTskLstOfIns[k].empty()) {
            map<double,int>::iterator it = tmTskLstOfIns[k].end(); --it;
            lstTskOfIns.insert((*it).second);
        }
    }
    vector<int> tskSchLst(comConst.NumOfTasks);
    IndexSortByValueDesc(tskSchLst, ch.EndTime); //根据完成时间反向计算
    vector<double> oldStrTime = ch.StartTime;    //保存原始的任务开始时间
    for (int i = 1; i < comConst.NumOfTasks; ++i) {
        int crnTsk = tskSchLst[i];
        if (lstTskOfIns.find(crnTsk) == lstTskOfIns.end()){  //crnTsk不是实例要处理的最后一个任务
            map<double,int>::iterator it = tmTskLstOfIns[ch.InsAlcLst[crnTsk]].find(oldStrTime[crnTsk]);  //查找当前任务的位置，查找时要用原始的任务开始时间，
            if (it == tmTskLstOfIns[crnTsk].end()){ //测试用-xy
                cout<<"crnTsk and its start time do not find in RscPrvImp_A!"; exit(2);
            }
            ++it;  //it指向的是同一个资源实例上当前任务的下一个任务
            if (ch.StartTime[(*it).second] - PrecisionValue > ch.EndTime[crnTsk]) {
                double slkTime = InfInt;
                for (int chlTsk: Tasks[crnTsk].children) {
                    double st = ch.StartTime[chlTsk] - ch.EndTime[crnTsk];
                    if (ch.InsAlcLst[crnTsk] != ch.InsAlcLst[chlTsk])
                        st -= ParChildTranFileSizeSum[crnTsk][chlTsk]/XY_MIN(VMTypes[ch.TypeLst[ch.InsAlcLst[crnTsk]]].bw, VMTypes[ch.TypeLst[ch.InsAlcLst[chlTsk]]].bw);
                    if (st + PrecisionValue < 0) { //测试用-xy
                        cout<<endl<<"slack time is less 0 in  RPI_ShfTsk2, it is wrong !" ; exit(2);
                    }
                    if (st + PrecisionValue < slkTime)
                        slkTime = st;
                }
                if (slkTime > PrecisionValue) {
                    double ShiftTime = XY_MIN(slkTime, ch.StartTime[(*it).second] - ch.EndTime[crnTsk]);
                    ch.StartTime[crnTsk] += ShiftTime;
                    ch.EndTime[crnTsk] += ShiftTime;
//                    cout << endl << "there is a slack time ! ";
                }
            }
        }
    }
    double cst = 0;
    for ( int k = 0; k <= maxInsId; ++k ) { //遍历每个实例
        if (!tmTskLstOfIns[k].empty()) {
            map<double,int>::iterator pFrt = tmTskLstOfIns[k].begin(), pLst = tmTskLstOfIns[k].end(); --pLst;
            int frtTsk = (*pFrt).second, lstTsk = (*pLst).second;
            double LeaseTime = XY_MAX(ch.EndTime[lstTsk] - ch.StartTime[frtTsk] + VMTypes[ch.TypeLst[k]].BootTime, VMTypes[ch.TypeLst[k]].MinLeaseTime);
            cst = cst + VMTypes[ch.TypeLst[k]].FixedCost + VMTypes[ch.TypeLst[k]].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[ch.TypeLst[k]].ChargeUnit);
        }
    }
    if (cst + PrecisionValue < ch.Cost) {
        ch.Cost = cst;
//        cout<<endl<<"RPI_M is Ok ! ";
        return 1;
    }
    return 0;
}

int RPI_SplIns1(chromosome& ch) { //对每个实例只进行一次分裂操作；
//    if (ch.IsFrw == false) {
//        cout<< endl<< "ch.IsFrw is false in RPI_IS"; exit(2);
//    }
    int maxInsId = *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end());
    int crnMaxInsId = maxInsId;
    vector<map<double,int>> tmTskLstOfIns(maxInsId + 1); //为利用set的排序功能把StartTime[i]作为key，即：map[StartTime[i]] = i；
    for (int i = 0; i < comConst.NumOfTasks; ++i) {      //把任务及其开始时间添加到各个实例的时间任务列表中去
        tmTskLstOfIns[ch.InsAlcLst[i]][ch.StartTime[i]] = {i};
    }
    for ( int k = 0; k <= maxInsId; ++k ) {
        if (!tmTskLstOfIns[k].empty()) {
            map<double,int>::iterator crnIt = tmTskLstOfIns[k].begin(), lstIt = tmTskLstOfIns[k].end(); --lstIt;
            for (; crnIt != lstIt ; ++crnIt) {
                double crnCmpUntTime = VMTypes[ch.TypeLst[k]].ChargeUnit *
                        ceil((ch.EndTime[(*crnIt).second] - (*(tmTskLstOfIns[k].begin())).first + VMTypes[ch.TypeLst[k]].BootTime)/ VMTypes[ch.TypeLst[k]].ChargeUnit);
                map<double,int>::iterator nxtIt = crnIt; ++nxtIt;
                //不考虑实例的固定起租成本，即适用于FixedCost=0的情况，下一个任务的开始时间落在其前一个任务所在的租用时间段内，则跳过！
                if ((*nxtIt).first - VMTypes[ch.TypeLst[k]].BootTime + PrecisionValue < (*(tmTskLstOfIns[k].begin())).first - VMTypes[ch.TypeLst[k]].BootTime + crnCmpUntTime) continue;
                int flg = 0;
                for (map<double,int>::iterator it1 = nxtIt; it1 != tmTskLstOfIns[k].end(); ++it1) {
                    for (map<double,int>::iterator it2 = tmTskLstOfIns[k].begin(); it2 != nxtIt; ++it2) {
                        if (Tasks[(*it1).second].parents.find((*it2).second) != Tasks[(*it1).second].parents.end()) {
                            if ((*it1).first - ch.EndTime[(*it2).second] + PrecisionValue < ParChildTranFileSizeSum[(*it2).second][(*it1).second]/VMTypes[ch.TypeLst[k]].bw) {
                                flg = 1;  break;
                            }
                        }
                    }
                    if (flg == 1)  break;
                }
                if (flg == 1) continue;
                // /*进行实例分裂操作*/
                ++crnMaxInsId;
                for (map<double,int>::iterator it = nxtIt; it != tmTskLstOfIns[k].end(); ++it) ch.InsAlcLst[(*it).second] = crnMaxInsId;
                ch.TypeLst[crnMaxInsId] = ch.TypeLst[k];
                break;
            }
        }
    }
    if (crnMaxInsId > maxInsId) {
//        cout<<endl<<"instance split is ok!";
        double oldMs = ch.MakeSpan, oldCst = ch.Cost;
        RecodeInsAndType(ch);
        DecodeChrom_NI(ch);
        if (fabs(oldMs - ch.MakeSpan) > PrecisionValue || oldCst + PrecisionValue < ch.Cost)
            cout<<endl<<oldMs<<", "<<oldCst<<"; New: "<<ch.MakeSpan<<", "<<ch.Cost;
        chromosome temCh = ch;
        DecodeChrom_NI(temCh);
        if (temCh.Cost + PrecisionValue < ch.Cost) {
            ch = temCh;
            cout<<endl<<oldMs<<", "<<oldCst<<"; New: "<<ch.MakeSpan<<", "<<ch.Cost;
        }
        return 1;
    } else
        return 0;
}

int RPI_SplIns2(chromosome& ch) {  //对每个实例可进行多次分裂操作
//    if (ch.IsFrw == false) {
//        cout<< endl<< "ch.IsFrw is false in RPI_IS"; exit(2);
//    }
    int maxInsId = *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end());
    int crnMaxInsId = maxInsId;
    vector<map<double,int>> tmTskLstOfIns(maxInsId + 1); //为利用set的排序功能把StartTime[i]作为key，即：map[StartTime[i]] = i；
    for (int i = 0; i < comConst.NumOfTasks; ++i) {      //把任务及其开始时间添加到各个实例的时间任务列表中去
        tmTskLstOfIns[ch.InsAlcLst[i]][ch.StartTime[i]] = {i};
    }
    for ( int k = 0; k <= maxInsId; ++k ) {
        if (!tmTskLstOfIns[k].empty()) {
            map<double,int>::iterator crnIt = tmTskLstOfIns[k].begin(), lstIt = tmTskLstOfIns[k].end(); --lstIt;
            map<double,int>::iterator newStartIt = tmTskLstOfIns[k].begin(); //newStartIt记录了分裂后的新的开始，即分裂出来的实例的首任务的指针
            for (; crnIt != lstIt ; ++crnIt) {
                double crnCmpUntTime = VMTypes[ch.TypeLst[k]].ChargeUnit *
                                       ceil(XY_MAX(ch.EndTime[(*crnIt).second]-(*newStartIt).first+VMTypes[ch.TypeLst[k]].BootTime, VMTypes[ch.TypeLst[k]].MinLeaseTime)/VMTypes[ch.TypeLst[k]].ChargeUnit);
                map<double,int>::iterator nxtIt = crnIt; ++nxtIt;
                //不考虑实例的固定起租成本，即适用于FixedCost=0的情况，下一个任务的开始时间落在其前一个任务所在的租用时间段内，则跳过！
                if ((*nxtIt).first - VMTypes[ch.TypeLst[k]].BootTime - PrecisionValue < (*newStartIt).first - VMTypes[ch.TypeLst[k]].BootTime + crnCmpUntTime) continue;
                int flg = 0;
                for (map<double,int>::iterator it1 = nxtIt; it1 != tmTskLstOfIns[k].end(); ++it1) {
                    for (map<double,int>::iterator it2 = newStartIt; it2 != nxtIt; ++it2) {
                        if (Tasks[(*it1).second].parents.find((*it2).second) != Tasks[(*it1).second].parents.end()) {
                            if ((*it1).first - ch.EndTime[(*it2).second] + PrecisionValue < ParChildTranFileSizeSum[(*it2).second][(*it1).second]/VMTypes[ch.TypeLst[k]].bw) {
                                flg = 1;  break;
                            }
                        }
                    }
                    if (flg == 1)  break;
                }
                if (flg == 1) continue;
                // /*flg=0表明找到了一个空闲间隙可以进行实例分裂，以下进行实例分裂操作*/
                ++crnMaxInsId;
                for (map<double,int>::iterator it = nxtIt; it != tmTskLstOfIns[k].end(); ++it) ch.InsAlcLst[(*it).second] = crnMaxInsId;
                ch.TypeLst[crnMaxInsId] = ch.TypeLst[k];
                newStartIt = nxtIt;  //如果进行了分裂，那么要更新分裂后实例的首任务的指针；
            }
        }
    }
    if (crnMaxInsId > maxInsId) {
//        cout<<endl<<"instance split is ok!";
        double oldMs = ch.MakeSpan, oldCst = ch.Cost;
        RecodeInsAndType(ch);
        DecodeChrom_NI(ch);
        if (fabs(oldMs - ch.MakeSpan) > PrecisionValue || oldCst + PrecisionValue < ch.Cost)
            cout<<endl<<oldMs<<", "<<oldCst<<"; New: "<<ch.MakeSpan<<", "<<ch.Cost;
        chromosome temCh = ch;
        DecodeChrom_NI(temCh);
        if (temCh.Cost + PrecisionValue < ch.Cost) {
            ch = temCh;
            cout<<endl<<oldMs<<", "<<oldCst<<"; New: "<<ch.MakeSpan<<", "<<ch.Cost;
        }
        return 1;
    } else
        return 0;
}

void InitProModelOfOdr1(vector<vector<double>>& PMO, vector<int>& NumOfAncestors, vector<int>& NumOfNonDescendants, vector<double>& Rank_b) {
    //case1：PMO[i][k] represents the probability that the k-th scheduled task is task i
    vector<int> STS(comConst.NumOfTasks ,0);
    for(int i = 0; i < comConst.NumOfTasks; ++i) {
        int left  = NumOfAncestors[i], right = NumOfNonDescendants[i];
        for(int j = left; j < right; ++j) {
            PMO[i][j] = 1;  ++STS[j];
        }
    }
    for(int j = 0; j < comConst.NumOfTasks; ++j) {
        for(int i = 0; i < comConst.NumOfTasks; ++i) {
            PMO[i][j] = PMO[i][j] / STS[j];
        }
    }
//    for(int j = 0; j < comConst.NumOfTsk; ++j) { //initializing PMO based levels or  Rank_b
////        int sum = 0;
//        double sum = 0.0;
//        for(int i = 0; i < comConst.NumOfTsk; ++i) {
////            sum = sum + (TskLstInLvl.size() - LevelIdOfTask[i]) * PMS[i][j];
//            sum = sum + Rank_b[i] * PMS[i][j];
//        }
//        for(int i = 0; i < comConst.NumOfTsk; ++i) {
////            PMS[i][j] =  (TskLstInLvl.size() - LevelIdOfTask[i]) * PMS[i][j] / sum;
//            PMS[i][j] = Rank_b[i] * PMS[i][j] / sum;
//        }
//    }
//    for (int  i = 0; i < comConst.NumOfTasks; ++i) { //测试用-xy
//        double sum = 0;
//        for (int j = 0; j < comConst.NumOfTasks; ++j) {
//            sum += PMO[j][i];
//        }
//        if (fabs(sum - 1) > PrecisionValue) {
//            cout<<endl<<"Intializing PMO is wrong!"; exit(2);
//        }
//    }
}

void InitProModelOfOdr2(vector<vector<double>>& PMO) {
    //case2：PMO[i][k] represents the probability that task i is scheduled before task k
    for(int i = 0; i < comConst.NumOfTasks; ++i) {
        PMO[i][i] = 1;                            //用了相乘来计算选择概率，故定义为1
        for (int k = i+1; k < comConst.NumOfTasks; ++k) {
            if (Descendants[i].find(k) != Descendants[i].end()) {
                PMO[i][k] = 1; PMO[k][i] = 0;    //k is a descendant of i, namely i is an ancestor of k
                continue;
            }
            if (Ancestors[i].find(k) != Ancestors[i].end()) {
                PMO[i][k] = 0; PMO[k][i] = 1;    //k is an ancestor of i,
                continue;
            }
//            PMS[i][k] = 0.5;  PMS[k][i] = 0.5;
            PMO[i][k] = (1+ 1.0*(LevelIdOfTask[k]-LevelIdOfTask[i])/TaskListInLevel.size()) * 0.5;  PMO[k][i] = 1 - PMO[i][k];  //There dose not exist dependency between i and i1;
        }
    }
}

void InitProModelOfGrp(vector<vector<double>>& PMG) {
    for(int i = 0; i < comConst.NumOfTasks; ++i) {
        for(int k = 0; k < i; ++k) {  //此处要求虚拟机类型编号从0开始连续编号；可以考虑把第二维定义为map结构，这样可以不从0开始连续编号
            PMG[i][k] =  1.0 / comConst.NumOfTasks;
        }
        PMG[i][i] = 1 - i * 1.0 / comConst.NumOfTasks;
    }
//    for (int  i = 0; i < comConst.NumOfTasks; ++i) { //测试用-xy
//        double sum = 0;
//        for (int j = 0; j < comConst.NumOfTasks; ++j) {
//            sum += PMG[i][j];
//        }
//        if (fabs(sum - 1) > PrecisionValue) {
//            cout<<endl<<"initializing PMG is wrong!"; exit(2);
//        }
//    }
}

void InitProModelOfRsc1(vector<map<int, double>>& PMR) {
    //case1：PMR[i][k] represents the probability that the type of resource allocated to task i is k
    for(int i = 0; i < comConst.NumOfTasks; ++i) {
        for (int k: VMTypeSet) {
            PMR[i][k] = 0;
        }
    }
    for(int i = 0; i < comConst.NumOfTasks; ++i) {
        for(int k : Tasks[i].ElgRsc) {
            PMR[i][k] = 1.0 / Tasks[i].ElgRsc.size();
        }
    }
    for (int  i = 0; i < comConst.NumOfTasks; ++i) { //测试用-xy
        double sum = 0;
        for (int k: VMTypeSet) {
            sum += PMR[i][k];
        }
        if (fabs(sum - 1) > PrecisionValue) {
            cout<<endl<<"initializing PMR is wrong!"; exit(2);
        }
    }
}

void InitProModelOfRsc2(vector<map<int, double>>& PMR) {
    //case2：PMR[i][k] represents the probability that the type of instance i is k;
    for(int j = 0; j < comConst.NumOfTasks; ++j) {
        for (int k: VMTypeSet) {
            PMR[j][k] = 1.0/comConst.NumOfVMTypes;
        }
    }
    for (int  j = 0; j < comConst.NumOfTasks; ++j) { //测试用-xy
        double sum = 0;
        for (int k: VMTypeSet) {
            sum += PMR[j][k];
        }
        if (fabs(sum - 1) > PrecisionValue) {
            cout<<endl<<"initializing PMR is wrong in InitProModelOfRsc2 !"; exit(2);
        }
    }
}

chromosome GnrTskLstOfChr1_prp(vector<vector<double>>& PMO, vector<double>& eta_TSO) {
    //case1：PMO[i][k] represents the probability that the k-th scheduled task is task i
    chromosome chrom;
    initChrom(chrom);
    vector<int> upr(comConst.NumOfTasks,0);
    list<int> RTI;  //也可以考虑用set结构，效果差不多！-xy
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i] == 0)  RTI.push_back(i);
    }
    for (int i = 0; i < comConst.NumOfTasks; i++) {
        double sum = 0;
        for(int k : RTI){
            sum += PMO[k][i] * eta_TSO[k];
        }
        vector<double> SltProb(comConst.NumOfTasks);
        for (int k : RTI) {
            SltProb[k] = PMO[k][i] * eta_TSO[k] / sum;
        }
        if (RandomDouble(0,1) < inputParameter.prp) {
            double MaxPrt = -1;
            for (int k : RTI)  {
                if (MaxPrt + PrecisionValue < SltProb[k]) {
                    MaxPrt = SltProb[k];
                    chrom.TskSchLst[i] = k;
                }
            }
        } else {
            double rnd = rand()%1000 / 1000.0;
            double ProbSum = 0;
            for (int k : RTI) {
                ProbSum += SltProb[k];
                if (rnd + PrecisionValue < ProbSum) {
                    chrom.TskSchLst[i] = k;
                    break;
                }
            }
        }
        RTI.erase(find(RTI.begin(), RTI.end(), chrom.TskSchLst[i]));
        for (int chlTsk: Tasks[chrom.TskSchLst[i]].children) {
            upr[chlTsk]--;
            if (upr[chlTsk] == 0){
                RTI.push_back(chlTsk);
            }
        }
    }
    return chrom;
}

chromosome GnrTskLstOfChr2_prp(vector<vector<double>>& PMO, vector<double>& eta_TSO) {
    //case2：PMO[i][k] represents the probability that task i is scheduled before task k
    chromosome chrom;
    initChrom(chrom);
    list<int> RTI;
    vector<int> upr(comConst.NumOfTasks);
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i] == 0)  RTI.push_back(i);
    }
    for (int i = 0; i < comConst.NumOfTasks; i++) {
        double sum = 0;
        vector<double> sample = eta_TSO;
        for (int i1 : RTI){
            for (int i2: RTI) {
                sample[i1] = sample[i1] * PMO[i1][i2];
            }
            sum = sum + sample[i1];
        }
        if (rand()%1000/1000.0 < inputParameter.prp ) { //select the task with the largest priority
            double MaxPrt = -1;
            for (int i1 : RTI)  {
                if (MaxPrt + PrecisionValue < sample[i1]) {
                    MaxPrt = sample[i1];
                    chrom.TskSchLst[i] = i1;
                }
            }
        } else {      //select a task using roulette according to their selection probabilities
            double rnd = rand()%1000 / 1000.0;
            double CumProb = 0;
            for (int i1: RTI) {
                CumProb = CumProb + (sample[i1] / sum);
                if (rnd + PrecisionValue < CumProb ) {
                    chrom.TskSchLst[i] = i1;
                    break;
                }
            }
        }
        RTI.erase(find(RTI.begin(),RTI.end(),chrom.TskSchLst[i]));
        for (int chlTsk: Tasks[chrom.TskSchLst[i]].children) {
            upr[chlTsk]--;
            if (upr[chlTsk] == 0){
                RTI.push_back(chlTsk);
            }
        }
    }
    return chrom;
}

void GnrInsLstOfChr(chromosome& chrom, vector<vector<double>>& PMG) {
    //分组生成实例分配列表，分配的实例编号满足RSGA1要求，即：InsAlcLst[0]=0; InsAlcLst[i] <= max{InsAlcLst[0],...,InsAlcLst[i-1]}+1
    vector<set<int>> NoElgRscOfGrp;  //记录每个组（实例）当前不能选用的资源类型，例如：NoElgRscOfGrp[j]表示实例/组j当前不能选用的资源类型的编号集合
    chrom.InsAlcLst[0] = 0;
    NoElgRscOfGrp.push_back(Tasks[0].NoElgRsc);
    for (int i = 1; i < comConst.NumOfTasks; i++) {
        set<int> tskSlct = {i};     //存放可以与任务i同组的任务
        for (int k = 0; k < i; ++k) {  //找出可以同组（有共同可使用资源实例）的任务
            set<int> allNoElgRsc = NoElgRscOfGrp[chrom.InsAlcLst[k]]; //获得任务k所在组（分配的实例）的不可使用的资源类型集合
            allNoElgRsc.insert(Tasks[i].NoElgRsc.begin(),Tasks[i].NoElgRsc.end());
            if (allNoElgRsc.size() < comConst.NumOfVMTypes) {
                tskSlct.insert(k);
            }
        }
        double sum = 0;
        for (int tsk: tskSlct) {
            sum += PMG[i][tsk];
        }
        double cmlProb = 0, rnd = rand()%1000 / 1000.0;
        for (int tsk: tskSlct) {
            cmlProb += PMG[i][tsk]/sum;
            if(rnd + PrecisionValue < cmlProb) {
                if (tsk != i) {
                    chrom.InsAlcLst[i] = chrom.InsAlcLst[tsk]; //任务i和任务tsk同组，即分配给同一个实例
                    NoElgRscOfGrp[chrom.InsAlcLst[tsk]].insert(Tasks[i].NoElgRsc.begin(),Tasks[i].NoElgRsc.end());
                } else {  //i与其前面（编号小于i）的所有任务都不同组，生成一个新组
                    chrom.InsAlcLst[i] = NoElgRscOfGrp.size(); //给任务i分配一个新的实例
                    NoElgRscOfGrp.push_back(Tasks[i].NoElgRsc); //初始化这个新实例当前不可使用的资源类型
                }
                break;
            }
        }
    }
//    int maxInsId = 0;
//    for (int i = 0; i < comConst.NumOfTasks; ++i) {   //测试用-xy
//        if (chrom.InsAlcLst[i] == maxInsId) {
//            maxInsId++;
//        }
//        if (chrom.InsAlcLst[i] > maxInsId) {
//            cout<<endl<< "the InsAlcLst generated by GnrInsLstOfChr is wrong !"; exit(2);
//        }
//    }
}

void GnrTypeLstOfChr1(chromosome& chrom, vector<map<int, double>>& PMR) {  //生成虚拟机类型列表
//    ///* 从该组/实例可以选取的资源类型集合中随机选择一个，适用于资源可得率小于1的情况 */
//    map<int,set<int>> NoElgRsc;        //记录实例不可使用的资源类型，
//    for (int i = 0; i < comConst.NumOfTasks; ++i) {
//        NoElgRsc[chrom.InsAlcLst[i]].insert(Tasks[i].NoElgRsc.begin(),Tasks[i].NoElgRsc.end()); //获得每个实例不可使用的资源类型；
//    }
//    for (auto elm: NoElgRsc) {    //从可使用的资源类型中随机选择一个
//        set<int> avlElgRsc = VMTypeSet;
//        avlElgRsc.erase(elm.second.begin(),elm.second.end()); //获得每个实例可使用的资源类型；
//        if (!avlElgRsc.empty()){
//            int rnb = rand()%avlElgRsc.size();
//            set<int>::iterator it = avlElgRsc.begin();
//            advance(it, rnb);
//            chrom.TypeLst[elm.first] = *it;
//        } else {
//            cout<<endl<<"the Instance/Group allocated is wrong ! there is no type that can be used for all tasks in the group";
//            exit(2);
//        }
//    }

    //case1：PMR[i][k] represents the probability that the type of resource allocated to task i is k; 给任务i分配的实例的类型为k的概率
/* 根据概率的和计算资源类型的被选择的概率，事先需要获得该组（实例）可以选取的资源类型集合，适用于资源可得率小于1的情况 */
    map<int,vector<int>> TskLstInIns; //记录个体中每个实例的任务列表，即同组的任务集合，且编号小的任务排在编号大的前面
    for (int i = 0; i < comConst.NumOfTasks; ++i) {  //把每个任务放到其对应实例的任务列表中，即对任务进行分组；
        TskLstInIns[chrom.InsAlcLst[i]].push_back(i);
    }
    map<int,set<int>> NoElgRsc;       //记录每个资源实例(组)不可使用的资源类型集合；
    for (int i = 0; i < comConst.NumOfTasks; ++i) { //获得每个资源实例不可使用的资源类型集合；
        NoElgRsc[chrom.InsAlcLst[i]].insert(Tasks[i].NoElgRsc.begin(),Tasks[i].NoElgRsc.end());
    }
    for (auto elm: NoElgRsc) {        //遍历所有的资源实例；
        set<int> avlElgRsc = VMTypeSet;
        avlElgRsc.erase(elm.second.begin(),elm.second.end()); //获得资源实例可使用的资源类型集合；
        if (!avlElgRsc.empty()){      //存在可以共同使用的资源类型
            map<int, double> sum;
            double tSum = 0;
            for (int k: avlElgRsc) {
                sum[k] = 0;
                for (int i: TskLstInIns[elm.first]) { //计算该组中任务选择类型k的概率之和
                    sum[k] += PMR[i][k];
                }
                tSum += sum[k];
            }
            double cmlPrb = 0, rnd = rand() % 1000 / 1000.0;
            for (int k: avlElgRsc) {
                cmlPrb += sum[k]/tSum;
                if (rnd + PrecisionValue < cmlPrb) {
                    chrom.TypeLst[elm.first] = k;
                    break;
                }
            }
        } else {
            cout<<endl<<"the Instance/Group allocated is wrong ! there is no type that can be used for all tasks in the group";
            exit(2);
        }
    }

/* 根据概率的和计算资源类型的被选择的概率，仅适用于任一资源可以处理任一任务（资源可得率为1）的情况 */
//    map<int,vector<int>> TskLstInIns;   //记录个体中每个实例的任务列表，即同组的任务集合，且编号小的任务排在大的前面
//    for (int i = 0; i < comConst.NumOfTasks; ++i) { //把每个任务放到其对应实例的任务列表中，即对任务进行分组；
//        TskLstInIns[chrom.InsAlcLst[i]].push_back(i);
//    }
//    for (auto elm: TskLstInIns) {
//        map<int, double> sum;
//        double tSum = 0;
//        for (int k: VMTypeSet) { //计算选择类型k的概率
//            sum[k] = 0;
//            for (int i: elm.second) {
//                sum[k] +=  PMR[i][k];
//            }
//            tSum += sum[k];
//        }
//        if (fabs(tSum) < PrecisionValue) {
//            cout<<endl<<"the sum of probabilities should not be equal to 0!"; exit(2);
//        }
//        double cmlPrb = 0, rnd = rand() % 1000 / 1000.0;
//        for (int k: VMTypeSet) {
//            cmlPrb += sum[k]/tSum;
//            if (rnd + PrecisionValue < cmlPrb) {
//                chrom.TypeLst[elm.first] = k;
//                break;
//            }
//        }
//    }

///* 根据概率的乘积计算资源类型的被选择的概率，可以不用获得每组（每个资源实例）可以选取的资源类型集合,适用于资源可得率小于1的情况 */
//    map<int,vector<int>> TskLstInIns;        //记录个体中每个实例的任务列表，即同组的任务集合，且编号小的任务排在大的前面
//    for (int i = 0; i < comConst.NumOfTasks; ++i) {
//        TskLstInIns[chrom.InsAlcLst[i]].push_back(i); //把每个任务放到其对应实例的任务列表中，即对任务进行分组；
//    }
//    for (auto elm: TskLstInIns) {
//        map<int,double> prd;
//        double tSum = 0;
//        for (int k: VMTypeSet) { //计算选择类型k的概率
//            prd[k] = 1;
//            for (int i: elm.second) {
//                prd[k] = prd[k] * PMR[i][k] * comConst.NumOfVMTypes; //放大comConst.NumOfVMTypes倍，防止prd[k]的值很小，被计算机认为是0
//            }
//            tSum = tSum + prd[k];
//        }
//        if (fabs(tSum) < PrecisionValue) {
//            cout<<endl<<"the sum of probabilities should not be equal to 0!"; exit(2);
//        }
//        double cmlPrb = 0, rnd = rand() % 1000 / 1000.0;
//        for (int k: VMTypeSet) {
//            cmlPrb += prd[k]/tSum;
//            if (rnd + PrecisionValue < cmlPrb) {
//                chrom.TypeLst[elm.first] = k;
//                break;
//            }
//        }
//    }
}

void GnrTypeLstOfChr2(chromosome& chrom, vector<map<int, double>>& PMR) {  //生成虚拟机类型列表
    //case2：PMR[i][k] represents the probability that the type of instance i is k; 实例i的类型为k的概率；
    /* 事先需要获得该组（实例）可以选取的资源类型集合，适用于资源可得率小于1的情况 */
    map<int,set<int>> NoElgRsc;   //记录每个资源实例(组)不可使用的资源类型集合；
    for (int i = 0; i < comConst.NumOfTasks; ++i) { //获得每个资源实例不可使用的资源类型集合；
        NoElgRsc[chrom.InsAlcLst[i]].insert(Tasks[i].NoElgRsc.begin(),Tasks[i].NoElgRsc.end());
    }
    for (auto elm: NoElgRsc) {    //遍历所有的资源实例；
        set<int> avlElgRsc = VMTypeSet;
        avlElgRsc.erase(elm.second.begin(),elm.second.end()); //获得每个资源实例可使用的资源类型集合；
        if (!avlElgRsc.empty()){  //存在可以共同使用的资源类型
            double Sum = 0;
            for (int k: avlElgRsc) {
                Sum += PMR[elm.first][k];
            }
            double cmlPrb = 0, rnd = rand() % 1000 / 1000.0;
            for (int k: avlElgRsc) {
                cmlPrb += PMR[elm.first][k]/Sum;
                if (rnd + PrecisionValue < cmlPrb) {
                    chrom.TypeLst[elm.first] = k;
                    break;
                }
            }
        } else {
            cout<<endl<<"the Instance/Group allocated is wrong ! there is no type that can be used for all tasks in the group";
            exit(2);
        }
    }
//    for (int j = comConst.NumOfTasks - 1; j >= 0; --j) {  //对没有赋值的实例类型进行随机赋值，确保更新PMR时不会出错！
//        if (chrom.TypeLst[j] == -1) {
//            int rnb = rand() % VMTypeSet.size();
//            set<int>::iterator it = VMTypeSet.begin();
//            advance(it, rnb);
//            chrom.TypeLst[j] = *it;
//        } else {
//            break; //由于在EDA中实例编号是从0开始连续编号的，故从后向前找到第一个已对实例类型进行赋值的，就可以退出循环；
//        }
//    }
//    for (int j = 0; j < comConst.NumOfTasks -1; ++j) {  //测试用-xy
//        if (chrom.TypeLst[j] < 0) {
//            cout<<endl<<"the type allocated is wrong in GnrTypeLstOfChr2!"; exit(2);
//        }
//    }
}

void UpdatePMO1(vector<vector<double>>& PMO, vector<chromosome>& elitPop){
    //case1：PMO[i][k] represents the probability that the k-th scheduled task is task i
    for(int i = 0; i < comConst.NumOfTasks; ++i) {
        for(int j = 0; j < comConst.NumOfTasks; ++j) {
            int count = 0;
            for(int n = 0; n < elitPop.size(); ++n) {
                if(elitPop[n].TskSchLst[i] == j) {
                    ++count;
                }
            }
            PMO[j][i] = (1 - inputParameter.theta1) * PMO[j][i] + inputParameter.theta1 * count / elitPop.size();
        }
    }
    for (int  i = 0; i < comConst.NumOfTasks; ++i) { //测试用-xy
        double sum = 0;
        for (int j = 0; j < comConst.NumOfTasks; ++j) {
            sum += PMO[j][i];
        }
        if (fabs(sum - 1) > PrecisionValue) {
            cout<<endl<<"updating PMO is wrong!"; exit(2);
        }
    }
}

void UpdatePMO2(vector<vector<double>>& PMO, vector<chromosome>& elitPop){
    //case2：PMO[i][k] represents the probability that task i is scheduled before task k
    vector<vector<int>> count(comConst.NumOfTasks, vector<int> (comConst.NumOfTasks, 0));
    for (int n = 0; n < elitPop.size(); ++n) {
        for (int i = 0; i < comConst.NumOfTasks; ++i) {
            for (int k = i + 1; k < comConst.NumOfTasks; ++k) {
                ++count[elitPop[n].TskSchLst[i]][elitPop[n].TskSchLst[k]];
            }
        }
    }
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        for (int k = i + 1; k < comConst.NumOfTasks; ++k) {
            PMO[i][k] = (1 - inputParameter.theta1) * PMO[i][k] + inputParameter.theta1 * (1.0 * count[i][k] / elitPop.size());
            PMO[k][i] = (1 - inputParameter.theta1) * PMO[k][i] + inputParameter.theta1 * (1.0 * count[k][i] / elitPop.size());
            if (abs(PMO[k][i]+PMO[i][k]-1) > PrecisionValue ) {
                cout << "Updating PMO is wrong in UpdatePMO2";
            }
            PMO[k][i] = 1 - PMO[i][k];
        }
    }
}

void UpdatePMG(vector<vector<double>>& PMG, vector<chromosome>& elitPop){
    vector<map<int,vector<int>>> chrInsTsk; //记录elitPop中每个个体的实例-任务列表（任务分组）；
    for (int n = 0; n < elitPop.size(); ++n) {
        map<int,vector<int>> insTsk;        //记录个体的实例-任务列表，即同组的任务集合，且编号小的在前面
        for (int i = 0; i < comConst.NumOfTasks; ++i) {
            insTsk[elitPop[n].InsAlcLst[i]].push_back(i); //把每个任务放到其对应实例的任务列表中，即对任务进行分组；
        }
        chrInsTsk.push_back(insTsk);
    }
    vector<vector<vector<double>>> vctChi;  //记录所有个体的chi值
    for (int n = 0; n < elitPop.size(); ++n) {
        vector<vector<double>> chi(comConst.NumOfTasks, vector<double>(comConst.NumOfTasks, 0));
        for (auto elm : chrInsTsk[n]) {
            chi[elm.second[0]][elm.second[0]] = 1;
            for (int i = elm.second.size()-1; i >0; --i) {
                for (int j = i-1; j >=0 ; --j) {
                    chi[elm.second[i]][elm.second[j]] = 1.0/i;
                }
            }
        }
        vctChi.push_back(chi);
    }
    vector<vector<double>> sumChi(comConst.NumOfTasks, vector<double>(comConst.NumOfTasks, 0));
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        for (int j = 0; j <=i; ++j) {
            for (int n = 0; n < elitPop.size(); ++n) {
                sumChi[i][j] += vctChi[n][i][j];
            }
           PMG[i][j] = (1 - inputParameter.theta2) * PMG[i][j] + inputParameter.theta2 * sumChi[i][j] / elitPop.size();
        }
    }
    for (int  i = 0; i < comConst.NumOfTasks; ++i) { //测试用-xy
        double sum = 0;
        for (int j = 0; j < comConst.NumOfTasks; ++j) {
            sum += PMG[i][j];
        }
        if (fabs(sum - 1) > PrecisionValue) {
            cout<<endl<<"updating PMG is wrong!"; exit(2);
        }
    }
}

void UpdatePMR1(vector<map<int,double>>& PMR, vector<chromosome>& elitPop){
    //case1：PMR[i][k] represents the probability that the type of resource allocated to task i is k; 给任务i分配的实例的类型为k的概率
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        for (int k: VMTypeSet) {
            int count = 0;
            for (int n = 0; n < elitPop.size(); ++n) {
                if (elitPop[n].TypeLst[elitPop[n].InsAlcLst[i]] == k) {
                    ++count;
                }
            }
            PMR[i][k] = (1 - inputParameter.theta3) * PMR[i][k] + inputParameter.theta3 *  count / elitPop.size();
        }
    }
    for (int  i = 0; i < comConst.NumOfTasks; ++i) { //测试用-xy
        double sum = 0;
        for (int k: VMTypeSet) {
            sum += PMR[i][k];
        }
        if (fabs(sum - 1) > PrecisionValue) {
            cout<<endl<<"updating PMR is wrong!"; exit(2);
        }
    }
}

void UpdatePMR2(vector<map<int,double>>& PMR, vector<chromosome>& elitPop){
    //case2：PMR[i][k] represents the probability that the type of instance i is k;
    vector<int> nmb(comConst.NumOfTasks, 0);
    for (chromosome ch: elitPop) {  //未使用的实例的类型不参与概率模型PMR的更新-xy
        int maxInsId = *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end());
        for (int j = maxInsId + 1; j < comConst.NumOfTasks; ++j) {
            ++nmb[j];
//            if (ch.TypeLst[j] != -1)  {  //测试用-xy
//                cout << endl << "it is wrong in UpdatePMR2!"; exit(2);
//            }
        }
    }
    for (int j = 0; j < comConst.NumOfTasks; ++j) {
        if (nmb[j] == elitPop.size()) continue;
        for (int k: VMTypeSet) {
            int count = 0;
            for (int n = 0; n < elitPop.size(); ++n) {
                if (elitPop[n].TypeLst[j] == k) {
                    ++count;
                }
            }
            PMR[j][k] = (1 - inputParameter.theta3) * PMR[j][k] + inputParameter.theta3 * count / (elitPop.size() - nmb[j]);
        }
    }
    for (int  j = 0; j < comConst.NumOfTasks; ++j) { //测试用-xy
        double sum = 0;
        for (int k: VMTypeSet) {
            sum += PMR[j][k];
        }
        if (fabs(sum - 1) > PrecisionValue) {
            cout<<endl<<"updating PMR is wrong in UpdatePMR2!"; exit(2);
        }
    }
}
