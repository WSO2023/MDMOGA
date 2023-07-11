//
// Created by xieyi on 2023/6/7.
//

#include <fstream>
#include "config.h"
#include "tools.h"
#include "GenOper.h"
#include "ChrOper.h"
#include "HMCMB.h"

vector<chromosome> runHMCMB(string XMLFile, string ResAlcFile, double& SchTime){
    clock_t start =clock();        //记录程序执行时间
    ReadFile(XMLFile, ResAlcFile); //读取文件
    CalculateLevelList();          //计算任务层次值
    CalculateDescendants();        //计算子孙任务
    CalculateAncestors();          //计算祖先任务
    vector<double> ww(comConst.NumOfTasks, 0);
    vector<vector<double>> cc(comConst.NumOfTasks, vector<double>(comConst.NumOfTasks,0));
    vector<double> rank_b(comConst.NumOfTasks, 0), rank_t(comConst.NumOfTasks, 0), rank_bt(comConst.NumOfTasks, 0);
    ClcAvrTskExcTime(ww);
    ClcAvrCmnTime(cc);  //考虑无需传输文件及有资格处理任务的虚拟机类型可能不同的情况下计算两个任务之间文件的平均传输时间
    Calculate_Rank_b(rank_b, ww, cc);
    Calculate_Rank_t(rank_t, ww, cc);
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
       rank_bt[i] = rank_b[i] + rank_t[i];
    }

    // /*估计最小预算和最小完工时间*/
    double refCst = CheapSchedule_xy();     //把所有任务分配给同一个虚拟机实例时处理费用最小的成本作为参考成本
    double refMs = FastSchedule_xy(rank_b); //把所有任务分配给无需等待的速度最快的虚拟机实例时的完工时间作为参考makespan
    vector<int> TskSchLst_b(comConst.NumOfTasks);
    vector<int> TskSchLst_t(comConst.NumOfTasks);
    IndexSortByValueDesc(TskSchLst_b, rank_b);
    IndexSortByValueAsc(TskSchLst_t, rank_t);

    // /*根据rank值产生一个任务调度顺序并根据效益比计算得出一个方案*/
    vector<chromosome> Pop;
    int nmbOfHEFT = 100+ 2 * comConst.NumOfTasks;
    for (int n = 0; n < nmbOfHEFT; ++n) {
        chromosome chrom ;
        initChrom(chrom);
        chrom.Weight = n / double(nmbOfHEFT-1);          //根据效益比生成资源配置方案时需要权重系数
        chrom.TskSchLst = TskSchLst_b;
        GnrRscAlc_BenefitRatio_MOELS(chrom,refCst,refMs);
        Pop.push_back(chrom);
        chrom.TskSchLst = TskSchLst_t;
        GnrRscAlc_BenefitRatio_MOELS(chrom,refCst,refMs);
        Pop.push_back(chrom);
        chrom.TskSchLst = GenerateTopoOrder_lev();
//        chrom.TskSchLst = GenerateTopoOrder_lev_bt(rank_bt);
        GnrRscAlc_BenefitRatio_MOELS(chrom,refCst,refMs);
        Pop.push_back(chrom);
    }
    sort(Pop.begin(), Pop.end(), sortPopByMsCstAsc_fstMs); //在生成ParetoFront前必须对种群中的个体按makespan从小到大排序，相同makespan则按cost小到大排序
    vector<chromosome> paretoFront = GnrParetoFront(Pop);
    SchTime = double (clock() - start) / CLOCKS_PER_SEC;
    return paretoFront;
}
