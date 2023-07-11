#include <fstream>
#include "MDMOGA.h"
#include "config.h"
#include "tools.h"
#include "GenOper.h"
#include "ChrOper.h"

vector<chromosome> runMDMOGA(string XMLFile, string ResAlcFile, double& SchTime, int& iteration){
    //记录程序执行时间
    clock_t start =clock();
    double runtime;
    ReadFile(XMLFile, ResAlcFile);//读取文件
    MDMOGA_ConfigInputParameter();//设置参数配置

    CalculateLevelList();          //计算任务层次值
    CalculateDescendants();        //计算子孙任务
    CalculateAncestors();          //计算祖先任务
    vector<double> ww(comConst.NumOfTasks, 0);
    vector<vector<double>> cc(comConst.NumOfTasks, vector<double>(comConst.NumOfTasks,0));
    vector<double> rank_b(comConst.NumOfTasks, 0), rank_t(comConst.NumOfTasks, 0);
    ClcAvrTskExcTime(ww);
    ClcAvrCmnTime(cc);  //考虑无需传输文件及有资格处理任务的虚拟机类型可能不同的情况下计算两个任务之间文件的平均传输时间
    Calculate_Rank_b(rank_b, ww, cc);
    Calculate_Rank_t(rank_t, ww, cc);
    vector<int> TskSchLst_b(comConst.NumOfTasks);
    vector<int> TskSchLst_t(comConst.NumOfTasks);
    IndexSortByValueDesc(TskSchLst_b, rank_b);
    IndexSortByValueAsc(TskSchLst_t, rank_t);

    // /*估计最小预算和最小完工时间*/
    double refCst = CheapSchedule_xy(); //获得把所有任务分配给同一个虚拟机实例时处理费用最小的成本，即估计的最小预算
//    refCst = 2;                         //对Montage 25 1.0， refCst的设置影响比较大，大的好一点-xy
    double refMs = FastSchedule_xy(rank_b); //获得把所有任务分配给无需等待的速度最快的虚拟机实例时的完工时间/Deadline，即估计的最小完工时间
    //根据rank值产生一个任务调度顺序并根据效益比计算得出一个方案
    set<chromosome> fstPop, temPop;
    chromosome Chrom ;
    initChrom(Chrom);
//    IndexSortByValueDesc(Chrom.TskSchLst,rank_b);
    int nmbOfHMCMB = 100+ 2 * comConst.NumOfTasks;
    for (int n = 0; n < nmbOfHMCMB; ++n) {
        Chrom.Weight = n / double(nmbOfHMCMB -1 );     //根据效益比生成资源配置方案时需要权重系数
        Chrom.TskSchLst = TskSchLst_b;
        GnrRscAlc_BenefitRatio_MOELS(Chrom, refCst, refMs);
        GnrTwoDimCoding(Chrom);
        temPop.insert(Chrom);
        Chrom.TskSchLst = TskSchLst_t;
        GnrRscAlc_BenefitRatio_MOELS(Chrom, refCst, refMs);
        GnrTwoDimCoding(Chrom);
        temPop.insert(Chrom);
        Chrom.TskSchLst = GenerateTopoOrder_lev();
        GnrRscAlc_BenefitRatio_MOELS(Chrom, refCst, refMs);
        GnrTwoDimCoding(Chrom);
        temPop.insert(Chrom);
    }

    vector<chromosome> paretoFront = GnrParetoFront(temPop);
    fstPop.insert(paretoFront.begin(),paretoFront.end());

    /*生成初始种群*/
    int nmb = 0;
    int ndNmb = inputParameter.NumOfChromPerPop - fstPop.size();
    while( fstPop.size() < inputParameter.NumOfChromPerPop ){
        chromosome chrom;
        initChrom(chrom);
        if(nmb < ndNmb * 0.5){ //一半个体的资源配置由最大效益比的启发式方法生成，且权重系数从0变到1
            chrom.Weight = 1 - (nmb * 2.0 / ndNmb);
            chrom.TskSchLst = GenerateTopoOrder_Ran();   //根据拓扑排序生成调度顺序
            GnrRscAlc_BenefitRatio_MOELS(chrom, refCst, refMs);  //由最大效益比的启发式方法生成个体的资源配置
            GnrTwoDimCoding(chrom);
            fstPop.insert(chrom);
            nmb++;
        } else {
            chrom.TskSchLst = GenerateTopoOrder_lev();
            GnrRscAlc_Random(chrom);//随机产生
            GnrTwoDimCoding(chrom); //生成二维编码，
            GnrOneDimCoding(chrom); //再由二维编码生成一维编码，使实例编号连续，因为在类型变异中直接修改一维编码的资源类型时需要其连续！！-xy
            DecodeChrom(chrom);
            fstPop.insert(chrom);
        }
    }
    /*对初始种群中的个体按帕累托支配关系和拥挤度距离进行*/
    vector<chromosome> population;
    while(population.size() < inputParameter.NumOfChromPerPop){
        vector<chromosome> temPF = GnrParetoFront(fstPop); //计算fstPop的第一层非支配解（ParetoFront）存放到temPF中，同时在fstPop中去除第一层非支配解-xy；
        sortParetoFrontByCrowdingDistance(temPF);          //计算拥挤距离并按拥挤距离从大到小进行排序
        population.insert(population.end(), temPF.begin(), temPF.end());
    }

    vector<double> A(inputParameter.NumOfChromPerPop);
    CalSlctProb_Rank(1+1.0/inputParameter.NumOfChromPerPop, A); //计算基于排序值的累计选择概率，
    while(1){
        ++iteration;
        runtime =double (clock() - start) / CLOCKS_PER_SEC;
        if (runtime > SchTime){
            SchTime = runtime;  break;
        }
        vector<chromosome> NewPopulation;
        for(int n = 0; n < inputParameter.NumOfChromPerPop; ++n){
            double rnd = double(rand()%100)/100;
            if(rnd <= inputParameter.MutationProb) {
                int nmb = rand() % 3;
                int parent = ChromSelect(A);  //变异的个体也要选择-xy
                chromosome OrgChrom = population[parent];
                if( nmb == 0){                //类型变异
                    Mutation_Type(OrgChrom);
                    NewPopulation.push_back(OrgChrom);
                } else if (nmb == 1) {        //分组变异
                    if ( Mutation_Group(OrgChrom) )
                        NewPopulation.push_back(OrgChrom);
                    else
                        --n;
                } else {                     //调度顺序变异
                    if ( Mutation_TskOrd(OrgChrom) )
                        NewPopulation.push_back(OrgChrom);
                    else
                        --n;
                }
            } else {
                int parent1 = ChromSelect(A);
                int parent2 = parent1;
                while (parent1 == parent2) {
                    parent2 = ChromSelect(A);
                }
                chromosome NewChrom;
                initChrom(NewChrom);
                if ( Cross_TwoDim(population[parent1], population[parent2], NewChrom) )
                    NewPopulation.push_back(NewChrom);
                else
                    --n;
            }
        }
        for (int n = 0; n < inputParameter.NumOfChromPerPop; ++n) {
            DecodeChrom(NewPopulation[n]);     //decode
            RPI_ShfTsk2(NewPopulation[n]);
        }
//        IFBSI(NewPopulation[0]);
//        LBCRI(NewPopulation[0]);
        set<chromosome> mrgPop;  //两代合一代并去重；
        mrgPop.insert(population.begin(), population.end());
        mrgPop.insert(NewPopulation.begin(), NewPopulation.end());
        vector<chromosome> nxtPop;
        while(nxtPop.size() < inputParameter.NumOfChromPerPop){
            vector<chromosome> temPF = GnrParetoFront(mrgPop); //计算第一层非支配解（ParetoFront）存放到temPF中，同时去除第一层非支配解-xy；
            sortParetoFrontByCrowdingDistance(temPF);          //按拥挤距离排序
            nxtPop.insert(nxtPop.end(), temPF.begin(), temPF.end());
        }
        population.assign(nxtPop.begin(), nxtPop.begin()+inputParameter.NumOfChromPerPop);  //此时Population存的就是选择后按非支配和拥挤度排序的当代种群
    }

    set<chromosome> setPop;
    setPop.insert(population.begin(),population.end());
    vector<chromosome> ParetoFront = GnrParetoFront(setPop);
//    sort(ParetoFront.begin(), ParetoFront.end(), sortPopByCostDesc);
    return ParetoFront;
}