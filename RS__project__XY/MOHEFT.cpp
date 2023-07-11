//
// Created by xieyi on 2023/5/13.
//

#include <fstream>
#include "config.h"
#include "tools.h"
#include "GenOper.h"
#include "ChrOper.h"
#include "MOHEFT.h"
vector<chromosome> runMOHEFT(string XMLFile, string ResAlcFile, double& SchTime){
/*  要求虚拟机类型编号连续且从0开始，虚拟机在进行启动所需的时间内已经可以接受父任务传输文件 */
    clock_t start = clock();         //记录程序执行时间
    ReadFile(XMLFile, ResAlcFile);  //读取文件
    CalculateLevelList();           //计算任务层次值
    vector<double> ww(comConst.NumOfTasks, 0);
    vector<double> rank_b(comConst.NumOfTasks, 0);
    ClcAvrTskExcTime(ww);
    Calculate_Rank_b(rank_b,ww);    //计算rank中在计算通信时间时使用了虚拟机类型的平均带宽
    int Numofchrom = 10+ ceil(0.1 * comConst.NumOfTasks);            //每次迭代保留的个体数量，也是最终的帕累托解集大小
    chromosome chrom;
    initChrom(chrom);
    IndexSortByValueDesc(chrom.TskSchLst,rank_b); //生成任务调度顺序
    struct chromIns{
        chromosome chr;         //储存个体
        vector<Instance> insVct;//储存部署的虚拟机实例
        bool operator<(const chromIns &otherChromIns)const { //sort函数排序时需要用到
            if (fabs(this->chr.MakeSpan - otherChromIns.chr.MakeSpan) >= PrecisionValue) { //先按MakeSpan从小到大排序
                return this->chr.MakeSpan < otherChromIns.chr.MakeSpan;
            } else {                                                                       //如果MakeSpan相同则按Cost从小到大排序
                return this->chr.Cost < otherChromIns.chr.Cost;
            }
        }
    };

    vector<chromIns>result;//储存每次迭代的结果
    chromIns chrIns;
    chrIns.chr = chrom;
    result.push_back(chrIns);

    for (int i = 0; i < comConst.NumOfTasks; ++i) {//代码在计算成本时时按照最短租赁时间为单位租用虚拟机，计算工期则按已有任务的最晚完成时间
        int tsk = chrom.TskSchLst[i];
        vector<chromIns> temResult;//用于储存每次迭代后产生的所有可能解与其对应的虚拟机集合，后续要进行筛选并将结果输出至result
        for (auto  rsl:result) {   //遍历上一次循环输出的结果个体
            for (int j: Tasks[tsk].ElgRsc) {  //遍历所有可用的虚拟机类型，任务将部署在新虚拟机实例上 -xy2
                chromosome temChrom = rsl.chr;//将遍历的个体输出至新变量，再修改新变量，因为每个个体需要遍历所有可行方案，不能直接修改！
                temChrom.InsAlcLst[tsk] = rsl.insVct.size();//新虚拟机的编号必然是已使用虚拟机的数量（从0开始）
                temChrom.TypeLst[rsl.insVct.size()] = j;
                double readyTime = 0;        //任务就绪时间
                for (int prn: Tasks[tsk].parents) {
                    int prnIns = temChrom.InsAlcLst[prn];
                    //因为是新创建一个虚拟机，父任务和子任务不可能在同一个虚拟机
                    double fft = temChrom.EndTime[prn] + ParChildTranFileSizeSum[prn][tsk] / XY_MIN(VMTypes[j].bw, VMTypes[temChrom.TypeLst[prnIns]].bw);
                    if (readyTime + PrecisionValue < fft)  readyTime = fft;
                }
                temChrom.StartTime[tsk] = XY_MAX(readyTime,VMTypes[j].BootTime);//因为是虚拟机处理的第一个任务，最早可开始时间就是其开始时间
                temChrom.EndTime[tsk] = temChrom.StartTime[tsk] + Tasks[tsk].length / VMTypes[j].ps;
                temChrom.MakeSpan = XY_MAX(temChrom.EndTime[tsk], temChrom.MakeSpan);
                double LeaseTime = XY_MAX(Tasks[tsk].length / VMTypes[j].ps + VMTypes[j].BootTime, VMTypes[j].MinLeaseTime);//新增加的虚拟机租用时间就是虚拟机的启动时间加上该任务的执行时间
                temChrom.Cost += VMTypes[j].FixedCost + VMTypes[j].UnitCost * ceil((LeaseTime - PrecisionValue) / VMTypes[j].ChargeUnit);//起租成本只能在添加虚拟机的时候加入！
                Instance temIns;
                temIns.id = rsl.insVct.size(); temIns.type = j;
                temIns.ITL.insert(InfInt * 1.0); temIns.ITL.insert(temChrom.EndTime[tsk]);
                if (fabs(VMTypes[j].BootTime - temChrom.StartTime[tsk]) > PrecisionValue){//如果不相同，则在启动时间到任务的开始时间之间也为可得时间
                    temIns.ITL.insert(VMTypes[j].BootTime);  temIns.ITL.insert(temChrom.StartTime[tsk]);
                }//可以直接初始化，这样就可以没有UpdateITL-xy 已改-szj
                temIns.LET = temChrom.EndTime[tsk], temIns.LST = temChrom.StartTime[tsk] - VMTypes[j].BootTime;
                vector<Instance> temInsVct = rsl.insVct;  temInsVct.push_back(temIns);
                temResult.push_back({temChrom,temInsVct});
            }
            for (int j = 0; j < rsl.insVct.size(); ++j) {//遍历所有已使用的虚拟机实例
                chromosome temChrom = rsl.chr;
                int InsType = temChrom.TypeLst[j];       //此虚拟机的类型
                if (Tasks[tsk].ElgRsc.find(InsType) != Tasks[tsk].ElgRsc.end()){//此虚拟机类型可以部署该任务
                    temChrom.InsAlcLst[tsk] = j;
                    double readyTime = 0;               //任务就绪时间
                    for (int prn: Tasks[tsk].parents) {
                        int prnIns = temChrom.InsAlcLst[prn];
                        double fft = temChrom.EndTime[prn];
                        //因为是已有虚拟机，考虑父任务和子任务在同一个虚拟机的情况
                        if (j != prnIns) {
                            fft += ParChildTranFileSizeSum[prn][tsk] / XY_MIN(VMTypes[InsType].bw,VMTypes[temChrom.TypeLst[prnIns]].bw);
                        }
                        if (readyTime + PrecisionValue < fft) readyTime = fft;
                    }
                    //因为是已有的虚拟机，不需要考虑启动时间
                    double ExcTime = Tasks[tsk].length / VMTypes[InsType].ps;
                    double StrTime = FindIdleTimeSlot(rsl.insVct[j].ITL, ExcTime, readyTime);
                    double EndTime = StrTime + ExcTime;
                    temChrom.StartTime[tsk] = StrTime;
                    temChrom.EndTime[tsk] = EndTime;
                    double NewLST = XY_MIN(StrTime,rsl.insVct[j].LST);
                    double NewLET = XY_MAX(EndTime,rsl.insVct[j].LET);
                    double oldLT = XY_MAX(rsl.insVct[j].LET - rsl.insVct[j].LST, VMTypes[InsType].MinLeaseTime);//未插入此任务前的虚拟机实例租用时间
                    double newLT = XY_MAX(NewLET - NewLST, VMTypes[InsType].MinLeaseTime);//插入此任务后的虚拟机实例租用时间
                    double CstAdded = VMTypes[InsType].UnitCost *
                                      (ceil((newLT-PrecisionValue)/VMTypes[InsType].ChargeUnit) - ceil((oldLT-PrecisionValue)/VMTypes[InsType].ChargeUnit));
                    temChrom.Cost += CstAdded; //更正-xy
                    temChrom.MakeSpan = XY_MAX(EndTime,temChrom.MakeSpan);//工期只可能不变或变为此次任务的结束时间  -xy2
                    vector<Instance> temInsVct = rsl.insVct;//通过中间变量将原个体的虚拟机实例转存并修改，成为输出结果
                    UpdateITL(temInsVct[j].ITL, StrTime, EndTime);
                    temInsVct[j].LST = NewLST; temInsVct[j].LET = NewLET;//更新虚拟机的最早开始时间和最晚结束时间
                    temResult.push_back({temChrom, temInsVct});
                }
            }
        }

        if ( temResult.size() > Numofchrom){        //只有数量达到需要筛选的时候才进入拥挤度排序输出，否则直接输出至下次循环
            sort(temResult.begin(),temResult.end());//在生成帕累托前沿前必须对种群中的个体按makespan从小到大排序，相同makespan则按cost小到大排序
            vector<chromIns> newResult;             //储存筛选出的帕累托前沿
            while(newResult.size() < Numofchrom){
                vector<chromIns> PF;                //计算第一层非支配解/帕累托前沿（ParetoFront）
                PF.push_back(temResult[0]);
                temResult.erase(temResult.begin());
                for (int n = 0; n < temResult.size(); ++n) {
                    if ( (temResult[n].chr.MakeSpan + PrecisionValue < PF[PF.size()-1].chr.MakeSpan) ||
                         (fabs(temResult[n].chr.MakeSpan - PF[PF.size()-1].chr.MakeSpan) < PrecisionValue && temResult[n].chr.Cost + PrecisionValue <PF[PF.size()-1].chr.Cost ) ) {
                        cout << endl << "the pop does not meet the sorting requirements!";
                    } //测试用-xy
                    if (temResult[n].chr.Cost + PrecisionValue < PF[PF.size()-1].chr.Cost) { //找到第一个成本小于PF中的最后一个个体（最后一个个体的成本是PF中成本最小的个体），即找到一个非支配解
                        PF.push_back(temResult[n]);
                        temResult.erase(temResult.begin() + n); //删除当前元素，指针移到下一个元素
                        --n;
                    }
                }
                if(PF.size() + newResult.size() < Numofchrom ){//如果插入当前的帕累托前沿，仍然还不够数，在无需排序，直接插入
                    newResult.insert(newResult.end(),PF.begin(),PF.end());
                }else{  //计算拥挤距离并按拥挤距离从大到小进行排序
                    double MaxCst = PF[0].chr.Cost, MinCst = PF[PF.size()-1].chr.Cost;     //1.1 取出最小的和最大的适应度值(成本)作为最小拥挤度和最大拥挤度
                    double MaxMS = PF[PF.size()-1].chr.MakeSpan, MinMS = PF[0].chr.MakeSpan;//1.2 取出最小的和最大的适应度值(工期)作为最小拥挤度和最大拥挤度
                    /* 1.3 把目标值最小的和最大拥挤距离设为无穷大，其他初始为0 */
                    for(int j = 1; j < PF.size()-1; ++j)   PF[j].chr.Distance = 0;
                    PF[0].chr.Distance = InfInt;  PF[PF.size()-1].chr.Distance = InfInt;
                    for(int j = 1; j < (PF.size()-1); ++j){ //1.4 计算中间部分的拥挤距离
                        PF[j].chr.Distance += (PF[j-1].chr.Cost - PF[j+1].chr.Cost) / (MaxCst - MinCst) + (PF[j+1].chr.MakeSpan - PF[j-1].chr.MakeSpan) / (MaxMS - MinMS);
                    }
                    sort(PF.begin(), PF.end(), [](const chromIns& v1,const chromIns& v2){return v1.chr.Distance > v2.chr.Distance;});  //2.按拥挤距离从大到小排序
                    int NeedNum = Numofchrom - newResult.size();                      //需要补充的个体数量
                    newResult.insert(newResult.end(),PF.begin(),PF.begin() + NeedNum);// -xy2
                    if(newResult.size()!= Numofchrom){
                        cout<<" 拥挤距离修改错误"<<endl;
                    }
                }
            }
            temResult = newResult; //取前Numofchrom个个体
        }
        result = temResult;       //本次迭代的输出结果成为下次迭代的输入
    }

    for (auto & rsl:result) {     //测试用-xy
        chromosome temCh = rsl.chr;
        DecodeChrom(temCh);
        if (fabs(temCh.Cost - rsl.chr.Cost) > PrecisionValue || fabs(temCh.MakeSpan - rsl.chr.MakeSpan) > PrecisionValue) {
            cout << endl << rsl.chr.Cost << ", "<<rsl.chr.MakeSpan<<"; Decode: "<<temCh.Cost<<", "<<temCh.MakeSpan;
        }
    }
    set<chromosome>temParetoFront;
    for (auto & rsl:result)  temParetoFront.insert(rsl.chr);
    vector<chromosome> ParetoFront = GnrParetoFront(temParetoFront); //ParetoFront中的个体已经按成本降序排了
    SchTime = double (clock() - start) / CLOCKS_PER_SEC;
    return ParetoFront;
}
