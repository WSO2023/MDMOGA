//
// Created by 11379 on 2021-08-04.
//

#include "MOELS.h"
#include <fstream>
#include "config.h"
#include "tools.h"
#include "GenOper.h"
#include "ChrOper.h"

vector<chromosome> runMOELS(string XMLFile, string ResAlcFile, double& SchTime, int& iteration){
    clock_t start = clock();       //记录程序开始执行时间
    double runtime;
    ReadFile(XMLFile, ResAlcFile);//读取文件
    MOELS_ConfigInputParameter(); //设置参数配置
    CalculateLevelList();         //计算任务层次值
    vector<double> ww(comConst.NumOfTasks, 0);
    vector<double> rank_b(comConst.NumOfTasks, 0);
    ClcAvrTskExcTime(ww);                      //任务平均执行时间
    Calculate_Rank_b(rank_b,ww);               //计算rank中在计算通信时间时使用了虚拟机类型的平均带宽
    double refCst = CheapSchedule_MOELS();        //估算把所有任务分配给最便宜虚拟机的一个实例时的成本在，并作为参考成本
    double refMs = FastSchedule_MOELS(rank_b);    //为每个任务分配一个处理速度最快的虚拟机实例时的makespan ，并作为参考makespan
    vector<chromosome> population;
//    set<chromosome> iniPop;
    /*初始化种群*/
//    while( iniPop.size() < inputParameter.NumOfChromPerPop ){
//        chromosome Chrom;
//        initChrom(Chrom);
//        int nmb = iniPop.size();
//        Chrom.Weight = (double)nmb / (inputParameter.NumOfChromPerPop - 1); //初始化权重设置
//        Chrom.TskSchLst = GenerateAChrom_order(rank_b);     //基于rank为概率进行选择生成任务调度顺序，每次生成调度顺序不一样
//        GnrRscAlc_BenefitRatio_MOELS(Chrom); //根据把任务分配给使部分工作流（所有已分配任务）的归一化时间成本最小的实例的方法生成虚拟机分配方案并计算cost和makespan
//        iniPop.insert(Chrom);
//    }
//    /*对初始种群中的个体按帕累托支配关系和拥挤度距离进行排序并放到population中*/
//    vector<chromosome> temNDVec;
//    while(population.size() < inputParameter.NumOfChromPerPop){
//        GnrParetoFront_set(iniPop, temNDVec);                //计算fstPop的第一层非支配解（ParetoFront）并按拥挤距离排序存放到temNDVec中，-xy；
//        for (chromosome chrom: temNDVec) iniPop.erase(chrom);//在fstPop中去除第一层非支配解，即temNDVec中的所有个体-xy
//        population.insert(population.end(), temNDVec.begin(), temNDVec.end());
//    }
    for (int n = 0; n < inputParameter.NumOfChromPerPop ; ++n) {
        chromosome Chrom;
        initChrom(Chrom);
        Chrom.Weight = (double)n / (inputParameter.NumOfChromPerPop - 1); //初始化权重设置
        Chrom.TskSchLst = GenerateAChrom_order(rank_b);     //基于rank为概率进行选择生成任务调度顺序，每次生成调度顺序不一样
        GnrRscAlc_BenefitRatio_MOELS(Chrom, refCst, refMs); //根据把任务分配给使部分工作流（所有已分配任务）的归一化时间成本最小的实例的方法生成虚拟机分配方案并计算cost和makespan
        population.push_back(Chrom);
    }
    sort(population.begin(), population.end(), sortPopByMsCstAsc_fstMs); //在生成非支配解(ParetoFront)前必须对种群中的个体按makespan从小到大排序，相同makespan则按cost小到大排序
    vector<chromosome> temPop;
    while (temPop.size() < inputParameter.NumOfChromPerPop) {
        vector<chromosome> PF = GnrParetoFront(population);  //计算种群的第一层非支配解（ParetoFront）并在种群中去除第一层非支配解，-xy；
        sortParetoFrontByCrowdingDistance(PF);               //按拥挤距离排序从大到小排序
        temPop.insert(temPop.end(), PF.begin(), PF.end());
    }
    population = temPop;

    for(int i = 0; i < inputParameter.NumOfChromPerPop; ++i){
    }
    /*迭代进化*/
    while(1){
        ++iteration;
        runtime =double (clock() - start) / CLOCKS_PER_SEC;
        if (runtime > SchTime){
            SchTime = runtime;  break;
        }
        vector<chromosome> NewPopulation;
        /*cross*/
        for (int n = 0; n < inputParameter.NumOfChromPerPop; n += 2){
            int parent1 = selectChrom_Tournament(inputParameter.NumOfChromPerPop);
            int parent2 = selectChrom_Tournament(inputParameter.NumOfChromPerPop);
            while (parent1 == parent2)
                parent2 = selectChrom_Tournament(inputParameter.NumOfChromPerPop);
            chromosome chrom1 = population[parent1];
            chromosome chrom2 = population[parent2];
            if(rand()%1000/1000.0 <= inputParameter.CrossProb){
                if(rand()%1000/1000.0 < 0.5){
                    /*交叉方法1 -- 互换两父体中的权重*/
                    double Wight;
                    Wight = chrom1.Weight;
                    chrom1.Weight = chrom2.Weight;
                    chrom2.Weight = Wight;
                } else{
                    /*交叉方法2 -- 基于拓扑排序的任务调度顺序交叉和基于SBX的权重交叉*/
                    CrossoverTaskOrder(chrom1, chrom2, true);  //基于拓扑排序的任务调度顺序交叉(右交叉，前面保存不变)
                    CrossoverWight(chrom1, chrom2);            //基于SBX(模拟二进制交叉)的权重交叉
                }
            }
            NewPopulation.push_back(chrom1);
            NewPopulation.push_back(chrom2);
        }
        /*Mutation*/
        for(int n = 0; n < inputParameter.NumOfChromPerPop; ++n){
            if(rand()%1000/1000.0 <= inputParameter.MutationProb){
                if(rand()%1000/1000.0 < 0.5){
                    MutationWight(NewPopulation[n]);      //变异方法1 -- 基于PM(多项式变异)的权重变异
                } else {
                    MutationTaskOrder(NewPopulation[n]);  //变异方法2 -- 基于拓扑排序的任务调度顺序单点变异
                }
            }
        }
        for (int n = 0; n < inputParameter.NumOfChromPerPop; ++n) {        //生成虚拟机分配方案并解码
            GnrRscAlc_BenefitRatio_MOELS(NewPopulation[n], refCst, refMs); //根据把任务分配给使部分工作流（所有已分配任务）的归一化时间成本最小的实例的方法生成虚拟机分配方案并计算cost和makespan
        }

        population.insert(population.end(), NewPopulation.begin(), NewPopulation.end());  //两代合一代(NSGA-II采用的策略)
        sort(population.begin(), population.end(), sortPopByMsCstAsc_fstMs); //在生成非支配解(ParetoFront)前必须对种群中的个体按makespan从小到大排序，如果makespan相同则cost小的在前
        /*从合并种群中选取pareto非支配解集合，形成下代种群*/
        vector<chromosome> nxtPop;
        while(nxtPop.size() < inputParameter.NumOfChromPerPop){
            vector<chromosome> PF = GnrParetoFront(population);  //计算种群的第一层非支配解（ParetoFront）并在种群中去除第一层非支配解，-xy；
            sortParetoFrontByCrowdingDistance(PF);               //按拥挤距离排序从大到小排序
            nxtPop.insert(nxtPop.end(), PF.begin(), PF.end());
        }
        population.assign(nxtPop.begin(),nxtPop.begin()+inputParameter.NumOfChromPerPop);
    }

    sort(population.begin(), population.end(), sortPopByMsCstAsc_fstMs); //在生成非支配解(ParetoFront)前必须对种群中的个体按makespan从小到大排序
    vector<chromosome> ParetoFront = GnrParetoFront(population);  //找出Population的非支配解(ParetoFront),只求出第一层
    return ParetoFront;
}