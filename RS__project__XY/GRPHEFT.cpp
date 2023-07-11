//
// Created by 11379 on 2021-08-16.
//

#include "GRPHEFT.h"
#include <fstream>
#include "config.h"
#include "tools.h"
#include "ChrOper.h"

vector<chromosome> runGRPHEFT(string XMLFile, string ResAlcFile, double& SchTime){
    //在选取生成的虚拟机实例集合Y时没有考虑实例对任务的支持情况，但在进行实例分配时考虑了实例对任务的支持情况，
    //因此本算法适用于资源可得率小于1的情况，但可能会出现很多情况下找不到解（return 0）
    //此外在选取生成的虚拟机实例集合Y时也没有考虑固定起租成本和最小起租时间，算法虽然能用于考虑固定起租成本和最小起租时间，但效果可能会变差！
    clock_t start =clock();        //记录程序执行时间
    ReadFile(XMLFile, ResAlcFile); //读取文件
    CalculateLevelList();          //计算任务层次值

    struct BdgMsCt{
        double budget, makespan, cost;
    };
    vector<BdgMsCt> BdgMsCt;      //存放预算及对应的工期和成本

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
    if (MinBdg[XMLFile] + PrecisionValue < MinUintCost ) {
        cout << endl << "The budget is less than the minimum unit cost of VM in GRPHEFT, and No feasible solution!";
        exit(0);
    }

    vector<int> TskSchLst(comConst.NumOfTasks);
    vector<double> ww(comConst.NumOfTasks, 0);
    vector<double> rank_b(comConst.NumOfTasks, 0);
    ClcAvrTskExcTime(ww);
    Calculate_Rank_b(rank_b,ww);
    IndexSortByValueDesc(TskSchLst,rank_b);  //内部好像可以优化-xy, 已优化-xy1
    set<chromosome> result;

    double stepValue = (MaxBdg[XMLFile] - MinBdg[XMLFile]) / comConst.NumOfTasks;
    double currValue = MinBdg[XMLFile] + 0 * stepValue;
    while (currValue - PrecisionValue < MaxBdg[XMLFile] - 0 * stepValue ){
        chromosome BestChrom;
        initChrom(BestChrom);
        BestChrom.MakeSpan = 999999999; BestChrom.Cost = 999999999;
        vector<Instance> LstY;
        //通过while循环确定最终方案
        for (int m = 0; m < TypeListBySrt.size(); ++m) {
            double Remaining_Budget = currValue;
            chromosome chrom;
            initChrom(chrom);
            chrom.TskSchLst = TskSchLst; //生成任务调度顺序
            vector<Instance> Y;          //选取的虚拟机实例集合
            map<int,int> NumOfInsUnused; //记录每种类型的实例的未使用数量
            for (int j = 0; j < TypeListBySrt.size(); ++j) {
                NumOfInsUnused[TypeListBySrt[j].type] = -9999999;
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
                    TemInst.ITL.insert(TypeListBySrt[j].BootTime); TemInst.ITL.insert(99999999 * 1.0);
                    TemInst.LET = 0, TemInst.LST = 99999999;
                    TemInst.state = 0;
                    Y.push_back(TemInst);
                }
                Remaining_Budget -= TypeListBySrt[j].UnitCost * NumOfInsUnused[TypeListBySrt[j].type];
                ++j;
            }

            //根据所选的虚拟机实例生成一种调度方案
//            Modified_HEFT_XY(chrom, Y, NumOfInsUnused);        //更正源代码有问题的方法-xy修改

//            //方法1：按完成时间从小到大分配给不增加成本的实例，如果不存在不增加成本的实例，则按完成时间从小到大分配给存在补偿实例的实例，否则返回无可行调度方案；
//            //即优先考虑不增加成本的，然后再考虑会增加成本的
//            if (Modified_HEFT_XY_1(chrom, Y, NumOfInsUnused) == 0) {
//                continue;
//            }

            //方法2：按完成时间从小到大分配给不增加成本或存在补偿实例的实例，否则返回无可行调度方案，即优先考虑时间最小，然后再考虑是否会增加成本
            //其中Modified_HEFT_XY_Imp相对于Modified_HEFT_XY_2增加了当前访问遍历的实例与最近访问遍历的实例是否类型相同且都是未被使用的判断，如果是则可以跳过当前的实例，以提高算法效率！
            if (Modified_HEFT_XY_Imp(chrom, Y, NumOfInsUnused) == 0) {
                continue;
            }

            if(chrom.Cost - PrecisionValue < currValue && chrom.MakeSpan - PrecisionValue < BestChrom.MakeSpan){
                BestChrom = chrom; LstY = Y;
            }
        }
        //记录预算和对应的结果
        BdgMsCt.push_back({currValue, BestChrom.MakeSpan, BestChrom.Cost});
        if (BestChrom.MakeSpan + PrecisionValue < 999999999) result.insert(BestChrom);
        currValue += stepValue;
    }
    SchTime = double (clock() - start) / CLOCKS_PER_SEC;
    //将原始数据结果输出到磁盘
    ofstream  fout("../result_GRPHEFT_Org.txt",ios::app);
    if (!fout) {
        cout << endl << "Open result_GRPHEFT_Org failure...\n";
        exit(1);
    }
    fout.setf(ios::fixed, ios::floatfield);
    fout.precision(3);
    fout << XMLFile << "; " << ResAlcFile;
    fout << endl << "工期： " ;
    for(int i = BdgMsCt.size() - 1 ; i >= 0; --i)  fout << BdgMsCt[i].makespan << " ";
    fout << endl << "成本： " ;
    for(int i = BdgMsCt.size() - 1 ; i >= 0; --i)  fout << BdgMsCt[i].cost << " ";
    fout << endl << "预算： " ;
    for(int i = BdgMsCt.size() - 1 ; i >= 0; --i)   fout << BdgMsCt[i].budget << " ";
    fout << endl;
    vector<chromosome> ParetoFront = GnrParetoFront(result);
    return ParetoFront;
}
