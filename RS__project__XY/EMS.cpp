//
// Created by xieyi on 2023/4/18.
//

#include "EMS.h"
#include <fstream>
#include "config.h"
#include "tools.h"
#include "GenOper.h"
#include "ChrOper.h"

//原论文种始化时所有虚拟机实例共共享同一种类型，故虚拟机可得率为1.0
vector<chromosome> runEMS(string XMLFile, string ResAlcFile, double& SchTime, int& iteration){
    //记录程序执行时间
    clock_t start =clock();
    double runtime;
    ReadFile(XMLFile, ResAlcFile);//读取文件
    EMS_ConfigInputParameter();   //设置配置
    CalculateLevelList();          //计算任务层次值
    vector<double> ww(comConst.NumOfTasks, 0);
    vector<vector<double>> cc(comConst.NumOfTasks, vector<double>(comConst.NumOfTasks,0));
    vector<double> rank_b(comConst.NumOfTasks, 0);
    ClcAvrTskExcTime(ww);
    ClcAvrCmnTime(cc);
    Calculate_Rank_b(rank_b,ww,cc);
    vector<chromosome> population;
    chromosome ch;
    initChrom(ch);
    IndexSortByValueDesc(ch.TskSchLst,rank_b);
    ch.Weight = 1;
    GnrRscAlc_BenefitRatio_MOELS(ch,1,1); //由于完工时间权重系数=1，因此生成的调度为给定调度顺序下不考虑成本的最小时间调度方案
    //通过GnrRscAlc_BenefitRatio生成的个体中实例从0开始连续编号，未使用的实例的类型值为初始值-1，
    //EMS中要求未使用的实例的类型也要有一个合法的编码，因此需要把虚拟机实例类型等于-1（即未使用的实例）的类型随机变为有效的整数（类型编码）
    for(int j = *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end())+1; j < comConst.NumOfTasks; ++j){
        ch.TypeLst[j] = rand() % comConst.NumOfVMTypes;
    }
    population.push_back(ch);
    ch.Weight = 0;
    GnrRscAlc_BenefitRatio_MOELS(ch,1,1); //由于完工时间权重系数=0，因此生成的调度为给定调度顺序下不考虑时间的最小成本调度方案
    for(int j = *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end())+1; j < comConst.NumOfTasks; ++j){
        ch.TypeLst[j] = rand() % comConst.NumOfVMTypes; //把虚拟机实例类型等于-1（即未使用的实例）的类型随机变为有效的整数（类型编码）
    }
    population.push_back(ch);

    while (population.size() < inputParameter.NumOfChromPerPop){
        chromosome chrom ;
        initChrom(chrom);
        //1.生成调度顺序
        chrom.TskSchLst = GenerateTopoOrder_Ran();
        //2.生成实例编号
        if(rand() % 1000 / 1000.0 < 0.5){
            for(int i = 0; i < comConst.NumOfTasks; ++i) chrom.InsAlcLst[i] = 0;
        } else {
            for(int i = 0; i < comConst.NumOfTasks; ++i) chrom.InsAlcLst[i] = rand() % comConst.NumOfTasks;
        }
        //3.生成实例类型 --- 原文假设任一资源可以处理任一任务,即资源可得率为1，即EMC不支持资源可得率小于1的情况！！
        for(int j = 0; j < comConst.NumOfTasks; ++j)  {
            map<int,VMType>::iterator it = VMTypes.begin();
            advance(it,rand() % comConst.NumOfVMTypes);
            chrom.TypeLst[j] = (*it).first;
        }
        population.push_back(chrom);
    }
    //解码
    for(int i = 0; i < population.size(); ++i){
        DecodeChrom(population[i]);
    }
    sort(population.begin(), population.end(), sortPopByMsCstAsc_fstMs); //在生成非支配解(ParetoFront)前必须对种群中的个体按makespan从小到大排序，相同makespan则按cost小到大排序
    vector<chromosome> srtPop;   //NSGA-II中的个体是按非支配关系和拥挤距离排序的，即首先按非支配层次从小到大进行排序，相同层次的则根据拥挤距离从大到小排序。
    while(srtPop.size() < inputParameter.NumOfChromPerPop){
        vector<chromosome> temPF = GnrParetoFront(population); //计算population的第一层ParetoFront存放到temPF中，同时在population中去除第一层ParetoFront-xy；
        sortParetoFrontByCrowdingDistance(temPF);              //计算拥挤距离并对temPF中的个体按拥挤距离从大到小进行排序
        srtPop.insert(srtPop.end(), temPF.begin(), temPF.end());
    }
    population = srtPop;

    while (1){
        ++iteration;
        runtime =double (clock() - start) / CLOCKS_PER_SEC;
        if (runtime > SchTime){
            SchTime = runtime;  break;
        }
        vector<chromosome> NewPopulation;
        //交叉 -- 交叉率为 1 -- 两种交叉方式是否是同时进行的
        for (int n = 0; n < inputParameter.NumOfChromPerPop; n += 2){
            //NSGA-II中采用锦标赛法进行个体选择
            int parent1 = selectChrom_Tournament(inputParameter.NumOfChromPerPop);
            int parent2 = selectChrom_Tournament(inputParameter.NumOfChromPerPop);
            while (parent1 == parent2)
                parent2 = selectChrom_Tournament(inputParameter.NumOfChromPerPop);
            chromosome chrom1 = population[parent1];
            chromosome chrom2 = population[parent2];
            if (rand()%2 == 0)
                CrossoverTaskOrder(chrom1, chrom2, true);//任务调度顺序交叉
            else
                CrossoverIns_EMS(chrom1, chrom2);        //虚拟机实例交叉
            NewPopulation.push_back(chrom1);
            NewPopulation.push_back(chrom2);
        }
        //变异 -- 变异率为 1 / comConst.NumOfTasks(1/n)
        for(int n = 0; n < inputParameter.NumOfChromPerPop; ++n){
            if( rand()%1000/1000.0 < 0.1 )
                mutation_EMS(NewPopulation[n]);
        }
        //解码
        for(int i = 0; i < population.size(); ++i){
            DecodeChrom(NewPopulation[i]);
        }
        population.insert(population.end(), NewPopulation.begin(), NewPopulation.end()); //两代合一代（NSGA-II中采用）
        sort(population.begin(), population.end(), sortPopByMsCstAsc_fstMs);//在生成ParetoFront前必须对种群中的个体按makespan从小到大排序，如果makespan相同则cost小的在前
        vector<chromosome> nxtPop;
        while(nxtPop.size() < inputParameter.NumOfChromPerPop){
            vector<chromosome> PF = GnrParetoFront(population);  //计算种群的第一层非支配解（ParetoFront）并在种群mrgPop中去除第一层非支配解，-xy；
            sortParetoFrontByCrowdingDistance(PF);               //对PF中的个体按拥挤距离从大到小排序
            nxtPop.insert(nxtPop.end(), PF.begin(), PF.end());
        }
        population.assign(nxtPop.begin(), nxtPop.begin()+inputParameter.NumOfChromPerPop);
    }
    sort(population.begin(), population.end(), sortPopByMsCstAsc_fstMs); //在生成非支配解(ParetoFront)前必须对种群中的个体按makespan从小到大排序
    vector<chromosome> ParetoFront = GnrParetoFront(population);  //找出Population的非支配解(ParetoFront),只求出第一层
    return ParetoFront;
}
