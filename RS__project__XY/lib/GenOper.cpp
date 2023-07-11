#include "tools.h"
#include "GenOper.h"
#include "ChrOper.h"


int ChromSelect(vector<double>& A) {
    double rnd = double(rand()%1000)/1000;
    for (int i = 0; i < A.size(); i++) {
        if ( rnd < A[i] ) {
            return i;
        }
    }
}

void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A) {
    for (int n = 0; n < A.size(); ++n) {
        A[n] = pow(RtOfSltPrb,A.size()-1-n) * (RtOfSltPrb - 1) / (pow(RtOfSltPrb, A.size()) - 1);
    }
    for (int n = 1; n < A.size(); ++n){
        A[n] = A[n] + A[n - 1];
    }
}

bool Cross_TwoDim(chromosome &chrom1 ,chromosome &chrom2, chromosome &newChrom){
    //随机生成一个整数
    //先判断Code.size()是否大于1，然后crossPoint = rand() % (chrom1.Code.size()-1) + 1; -xy3 //已改 -G
    if(chrom1.Code.size() > 1){
        int nmb = rand() % (chrom1.Code.size()-1) + 1;
        //set<list<int> > 转换成 vector<list<int> >
        list<list<int>> Code1 = chrom1.Code;
        list<list<int>> Code2 = chrom2.Code;
        list<list<int>> NewCode;                 //交叉后的输出结果
        list<list<int>>::iterator LstPt;
        //Step2 在父体1（Code1）中随机选择nmb个实例编号及其任务列表,把它们放到NewCode中；
        for(int j = 0; j < nmb; ++j){
            int b = rand() % Code1.size();       //随机Code1中的一个list
            LstPt = Code1.begin();
            advance(LstPt,b);
            NewCode.push_back(*LstPt);           //将此位置的list实例添加到NewCode中
            Code1.erase(LstPt);                  //删除Code1中此list实例
        }
        // Step3.1在父体2（Code2）中删除NewCode中已存在的元素
        list<int>::iterator LstNewPt, LstPt2;
        for(auto lstInNC: NewCode){                         //遍历NewCode中的每个list
            LstNewPt = lstInNC.begin();  LstNewPt++;        //从list的第二个元素开始
            for(; LstNewPt != lstInNC.end(); ++LstNewPt){   //遍历list中的每个元素
                for(list<list<int>>::iterator it = Code2.begin(); it != Code2.end(); ++it){   //遍历Code2容器中每个list
                    LstPt2 = (*it).begin(); LstPt2++;
                    list<int>::iterator temPt = find(LstPt2, (*it).end(), *LstNewPt);
                    if (temPt != (*it).end()){
                        (*it).erase(temPt);
                        break;
                    }
                }
            }
        }
        // Step3.2把Code2中的剩余元素添加到NewCode中
        for(list<list<int>>::iterator Pt = Code2.begin(); Pt != Code2.end(); ++Pt){
            if((*Pt).size() == 1)  continue;
            NewCode.push_back(*Pt);
        }
        //将交叉后的VecNewCode赋值到初始化后的chrom1中//重新定义一个新的染色体NewChrom，-xy3 //已改- G
        newChrom.Code = NewCode;
        return GnrOneDimCoding(newChrom);
    } else{
        newChrom = chrom1;
        return true;
    }
}

//分组变异
bool Mutation_Group(chromosome &chrom){
//    //set集合不允许修改里面的值，所以先set<list<int> > 转换成 vector<list<int>>-
//    vector<list<int>> VctCode;                       //转换后的VecCode1
//    VctCode.assign(chrom.Code.begin(), chrom.Code.end());
    list<int>::iterator LstPt;
    list<list<int>>::iterator dLstPt = chrom.Code.begin();
    int TaskId = rand() % comConst.NumOfTasks;    //随机选择一个任务
    //Step1：找出任务分配的资源实例
    int k = 0;
    for(; dLstPt!=chrom.Code.end(); dLstPt++){    //遍历vector中每一个list
        LstPt = (*dLstPt).begin();  LstPt++;      //从每个list第二个位置开始
        list<int>::iterator PtTskId = find(LstPt, (*dLstPt).end(), TaskId);
        if ( PtTskId != (*dLstPt).end() ) {       //找到任务在第k个list中
            (*dLstPt).erase(PtTskId);             //删除该任务
            break;
        }
        ++k;
    }
    //Step2：如实例k中只有任务TaskId则删除实例k，否则在实例k中删除任务TaskId，
    if((*dLstPt).size() == 1){                     //如果实例k中已经没有任务，（任务TaskId已被删除了）
        FndPosInsertTsk_xy(chrom.Code, TaskId, k); //找出可以处理任务TaskId的资源实例集合AINS，若AINS不为空从中随机选择一个插入该任务，否则租用一个新的插入
        chrom.Code.erase(dLstPt);                  //因为k已定，因此实例的删除必须在调用FndPosInsertTsk_xy后才进行，另外此时该实例已没有分配的任务因此要删除,
    } else {                                       //如果实例k中有两个及以上的任务
//        LstPt = (*dLstPt).begin(); LstPt++;
//        (*dLstPt).erase(find(LstPt,(*dLstPt).end(),TaskId)); //在实例k中删除任务TaskId
        if (rand() % 1000 / 1000.0 < 0.8) {
            FndPosInsertTsk_xy(chrom.Code,TaskId, k);
        } else {
            int nmb = rand() % Tasks[TaskId].ElgRsc.size();  //nmb为随机选择的位置
            set<int>::iterator pt = Tasks[TaskId].ElgRsc.begin();
            advance(pt,nmb);
            list<int> Instance;
            Instance.push_back(*pt);           //第一个位置存放的是虚拟机类型
            Instance.push_back(TaskId);        //第二个位置存放的是任务
            chrom.Code.push_back(Instance);    //将添加一个新的list到Code中
        }
    }
    return GnrOneDimCoding(chrom);
}

void FndPosInsertTsk_xy(list<list<int>>& Code, int TaskId, int k){
    list<int>::iterator LstPt;
    vector<int> AINS;
    for(int j = 0; j < Code.size(); ++j){        //Step3：找出可以处理任务TaskId的资源实例集合AINS(已存在&&可处理)
        if( j == k ) continue;                   //任务已存在的那个实例不添加到AINS中去
        list<list<int>>::iterator dLstPt = Code.begin(); advance(dLstPt,j);
        LstPt = (*dLstPt).begin();
        if (Tasks[TaskId].ElgRsc.find(*LstPt) != Tasks[TaskId].ElgRsc.end())  //查看是否在ElgRsc中
            AINS.push_back(j);                   //记录此实例在Code中的位置
    }
    if (!AINS.empty()) { //Step4：从AINS随机找一个实例，并在该实例中找出一个有效插入位置插入任务TaskId
        //4、Step：从AINS随机选择一个存储的数，并在合适位置插入TaskId，-xy3    //改过 -G
        int k2 = AINS[rand() % AINS.size()];
        //4.1 在选择的k1实例中找任务i的合法位置
        list<int>::iterator DscPoint;            //记录list容器中第一个子孙任务的位置指针，改名DscPoint,-xy3  //已改 -G
        list<int>::iterator AncPoint;            //记录list容器中最后一个祖父任务的位置指针，改名AncPoint,-xy3 //已改 -G
        list<list<int>>::iterator dLstPt = Code.begin(); advance(dLstPt,k2);
        DscPoint = (*dLstPt).begin();      DscPoint++;    //从第二个位置开始算
        for (; DscPoint != (*dLstPt).end(); ++DscPoint) { //从前向后找，找到第一个子孙任务
            if (Descendants[TaskId].find(*DscPoint) != Descendants[TaskId].end()){
                break;
            }
        }
        AncPoint = DscPoint; AncPoint--;  //如果DscPoint()为end(),则  AncPoint--  实际上指向最后一个任务（最后一个任务的后一个指针为end()）
        //AncPoint必在DscPoint前面，因此只需从DscPoint位置向前找,找到第一个祖先任务，如果找不到迭代器指向第一个位置
        for (; AncPoint != (*dLstPt).begin(); AncPoint--) {
            if(Ancestors[TaskId].find(*AncPoint) != Ancestors[TaskId].end()){
                break;
            }
        }
        int NumSite = distance(AncPoint, DscPoint);   //找出AncPoint和DscPoint中间有多少个位置，能否直接计算？-xy3 //已改 -G
        int a = (rand() % NumSite)+1;                 //找插入的位置
        advance(AncPoint, a);                         //表示从AncPoint指针开始，向后移动a个位置，此时AncPoint的位置就是要插入的位置
        (*dLstPt).insert(AncPoint,TaskId);            //在list的AncPoint位置插入TaskId
    } else { //Step6：随机选择一个可利用类型，租用一个新的虚拟机
        int nmb = rand() % Tasks[TaskId].ElgRsc.size(); //nmb为随机选择的位置
        set<int>::iterator pt = Tasks[TaskId].ElgRsc.begin();
        advance(pt,nmb);
        list<int> Instance;
        Instance.push_back(*pt);                     //第一个位置存放的是虚拟机类型
        Instance.push_back(TaskId);                  //第二个位置存放的是任务
        Code.push_back(Instance);                    //将添加一个新的list到VecCode容器中
    }
}

//调度顺序变异
bool Mutation_TskOrd(chromosome &chrom){
    list<int>::iterator pt;
    int TaskId = rand() % comConst.NumOfTasks;                 //随机选择一个任务
    //Step1：找出任务分配的资源实例
    int k = 0;
    list<list<int>>::iterator dLstPt = chrom.Code.begin();
    for(;  dLstPt != chrom.Code.end(); dLstPt++){                            //遍历vector中每一个list
        pt = (*dLstPt).begin();  pt++;                      //从每个list第二个位置开始
        pt = find(pt, (*dLstPt).end(), TaskId);               //直接在for循环中找，-xy3 //已修改 -G
        if(pt != (*dLstPt).end()){                            //找到任务在第k个list中
            break;
        }
    }
    //Step2：将任务插入到资源实例kt中合适的位置
    list<int>::iterator DscPoint;                   //记录list容器中第一个子孙任务的位置指针
    list<int>::iterator AncPoint;                   //记录list容器中最后一个祖先任务的位置指针
    DscPoint = pt;     DscPoint++;                 //指向了任务TaskId后边那个任务，如果TaskId为最后一个任务，则DscPoint为end()
    AncPoint = pt;     AncPoint--;                 //指向了任务TaskId前边那个任务，如果TaskId为第一个任务，则AncPoint为begin()（虚拟机类型）
    /* pt必须是原来位置的pt。 DscPoint AncPoint 迭代器不能停留在pt位置，否则同样会被删除掉 */
    (*dLstPt).erase(pt);                          //从第k个实例中删除任务TaskId, 删除之后，pt记录的位置还是任务TaskId的位置
    while (DscPoint != (*dLstPt).end() && find(Descendants[TaskId].begin(), Descendants[TaskId].end(),*DscPoint) == Descendants[TaskId].end())
        DscPoint++;
    while (AncPoint != (*dLstPt).begin() && find(Ancestors[TaskId].begin(), Ancestors[TaskId].end(),*AncPoint) == Ancestors[TaskId].end())
        AncPoint--;
    int Num = distance(AncPoint, DscPoint); //找出AncPoint和DscPoint中间有多少个位置，
    int a = (rand() % Num)+1;                        //找插入的位置
    advance(AncPoint, a);                            //表示从AncPoint指针开始，向后移动a个位置，此时AncPoint的位置就是要插入的位置
    (*dLstPt).insert(AncPoint,TaskId);              //在list的AncPoint位置插入TaskId
    return GnrOneDimCoding(chrom);
}

//类型变异
void Mutation_Type(chromosome &chrom){
    //1、Step：随机选择一个资源实例
    int k = rand() % chrom.Code.size();
    set<int> NoAvlTypeSt, AvlTypeSt = VMTypeSet;
    list<list<int>>::iterator lstPt = chrom.Code.begin();
    advance(lstPt,k);
    //Step2：找出能处理资源实例k的任务编号列表中所有任务的资源类型
    list<int>::iterator pt = (*lstPt).begin();
    AvlTypeSt.erase(*pt); //从AvlTypeSt中删除原来使用的实例类型，这样可以确保后面从中选择的类型与原来的必定不同；
    pt++;    //第二个位置为第一个任务
    for(; pt != (*lstPt).end(); ++pt){      //第k个list中每个任务的循环
        NoAvlTypeSt.insert(Tasks[*pt].NoElgRsc.begin(),Tasks[*pt].NoElgRsc.end());
    }
    AvlTypeSt.erase(NoAvlTypeSt.begin(), NoAvlTypeSt.end());

    //Step3 从AvlTypeSt中随机选择一个类型，赋值给随机选出的资源实例
    if(AvlTypeSt.size() > 0){
        int rnd = rand() % AvlTypeSt.size();
        set<int>::iterator it = AvlTypeSt.begin();
        advance(it, rnd);
        *((*lstPt).begin()) = *it;               //修改二维编码的中的资源类型
        chrom.TypeLst[k] = *((*lstPt).begin());  //（直接修改一维，不调用转一维函数)直接用k，这里要求一维中的实例编号连续-xy
    }
}

vector<int> GenerateAChrom_order(vector<double>& rank){
    vector<int> Order;
    vector<int> upr(comConst.NumOfTasks,0);
    vector<int> RTI;
    for (int i = 0; i < comConst.NumOfTasks; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i] == 0) {
            RTI.push_back(i);
        }
    }
    while (!RTI.empty()){
        int SlcTsk = chooseNextTask(rank, RTI);
        RTI.erase(find(RTI.begin(), RTI.end(), SlcTsk));
        for (int childId: Tasks[SlcTsk].children) {
            upr[childId] = upr[childId] - 1;
            if (upr[childId] == 0) {
                RTI.push_back(childId);
            }
        }
        Order.push_back(SlcTsk);
    }
    return Order;
}

int chooseNextTask(vector<double> rank, vector<int> RT){
    int sum = 0;
    for(int i = 0; i < RT.size(); ++i){
        sum += rank[RT[i]];
    }
    double randNum = rand() % 100 / 100.0;
    double ProbSum = 0;
    int chosenTask = -1;
    for(int i = 0; i < RT.size(); ++i){
        ProbSum += rank[RT[i]] / sum;
        if (ProbSum > randNum) {
            chosenTask = RT[i];
            break;
        }
    }
    return chosenTask;
}

chromosome Crossover_MOELS(chromosome &chrom1, chromosome &chrom2){
    int crossPoint = rand()%(comConst.NumOfTasks - 1) + 1;
    vector<int> tem2;
    tem2 = chrom2.TskSchLst;
    //从tem2中删除其值为chrom1交叉点前的元素
    for (int i = 0; i < crossPoint; i++)
        for (vector<int>::iterator it = tem2.begin(); it != tem2.end(); it++)
            if ((*it) == chrom1.TskSchLst[i]) {
                tem2.erase(it);
                break;
            }
    int n = 0;
    for (int i = crossPoint; i < comConst.NumOfTasks; i++) {
        chrom1.TskSchLst[i] = tem2[n];
        n++;
    }
    return chrom1;
}

void CrossoverTaskOrder(chromosome &chrom1, chromosome &chrom2, bool flag){
    int crossPoint = rand()%(comConst.NumOfTasks - 1) + 1;
    vector<int> tem1, tem2;
    tem1 = chrom1.TskSchLst;
    tem2 = chrom2.TskSchLst;
    if (flag) { //flag=true 是交叉点前保留不变
        //从tem1中删除其值为chrom2交叉点前的元素
        for (int i = 0; i < crossPoint; i++)
            for (vector<int>::iterator it = tem1.begin(); it != tem1.end(); it++)
                if ((*it) == chrom2.TskSchLst[i]) {
                    tem1.erase(it);
                    break;
                }
        //从tem2中删除其值为chrom1交叉点前的元素
        for (int i = 0; i < crossPoint; i++)
            for (vector<int>::iterator it = tem2.begin(); it != tem2.end(); it++)
                if ((*it) == chrom1.TskSchLst[i]) {
                    tem2.erase(it);
                    break;
                }
        int n = 0;
        for (int i = crossPoint; i < comConst.NumOfTasks; i++) {
            chrom1.TskSchLst[i] = tem2[n];
            chrom2.TskSchLst[i] = tem1[n];
            n++;
        }
    } else {  //flag=false 是交叉点后保留不变
        //从tem1中删除其值为chrom2 交叉点后的元素
        for (int i = crossPoint; i < comConst.NumOfTasks; i++)
            for (vector<int>::iterator it = tem1.begin(); it != tem1.end(); it++)
                if ((*it) == chrom2.TskSchLst[i]) {
                    tem1.erase(it);
                    break;
                }
        //从tem2中删除其值为chrom1交叉点后的元素
        for (int i = crossPoint; i < comConst.NumOfTasks; i++)
            for (vector<int>::iterator it = tem2.begin(); it != tem2.end(); it++)
                if ((*it) == chrom1.TskSchLst[i]) {
                    tem2.erase(it);
                    break;
                }
        int n = 0;
        for (int i = 0; i < crossPoint; i++) {
            chrom1.TskSchLst[i] = tem2[n];
            chrom2.TskSchLst[i] = tem1[n];
            n++;
        }
    }
}

void CrossoverTasklist(chromosome &chrom1, chromosome &chrom2){
    vector<int> tem1, tem2;
    tem1 = chrom1.TskSchLst;
    tem2 = chrom2.TskSchLst;
    int crossPoint = rand()%(comConst.NumOfTasks - 1) + 1;

    //从tem1中删除其值为chrom2交叉点前的元素
    for (int i = 0; i < crossPoint; ++i)
        for (vector<int>::iterator it = tem1.begin(); it != tem1.end(); it++)
            if ((*it) == chrom2.TskSchLst[i]) {
                tem1.erase(it);
                break;
            }
    //从tem2中删除其值为chrom1交叉点前的元素
    for (int i = 0; i < crossPoint; ++i)
        for (vector<int>::iterator it = tem2.begin(); it != tem2.end(); it++)
            if ((*it) == chrom1.TskSchLst[i]) {
                tem2.erase(it);
                break;
            }
    int n = 0;
    for (int i = crossPoint; i < comConst.NumOfTasks; ++i) {
        chrom1.TskSchLst[i] = tem2[n];
        chrom2.TskSchLst[i] = tem1[n];
        n++;
    }
}

void CrossoverInsAndType(chromosome &chrom1, chromosome &chrom2) {  //发明23的交叉方法
    double mrgPro = 0.8;          //lamda的值，实例合并概率
    int crossPoint = rand()%(comConst.NumOfTasks-1) + 1;
    chromosome chrom3, chrom4;
    initChrom(chrom3); initChrom(chrom4);
    chrom3.TskSchLst = chrom1.TskSchLst; chrom4.TskSchLst = chrom2.TskSchLst;
    chrom3.TypeLst = chrom1.TypeLst; chrom4.TypeLst = chrom2.TypeLst;
    for (int i = 0; i < crossPoint; ++i) {  //初始化交叉点crossPoint前（不包括crossPoint）的实例分配
        chrom3.InsAlcLst[i] = chrom1.InsAlcLst[i]; chrom4.InsAlcLst[i] = chrom2.InsAlcLst[i];
    }
    for (int j = crossPoint; j < comConst.NumOfTasks; ++j) {//初始化交叉点及其后的实例分配
        chrom3.InsAlcLst[j] = -1;  chrom4.InsAlcLst[j] = -1;
    }
    int cv1 = *max_element(chrom3.InsAlcLst.begin(), chrom3.InsAlcLst.begin() + crossPoint); //交叉点crossPoint前（不包括crossPoint）的实例编号的最大值
    int cv2 = *max_element(chrom4.InsAlcLst.begin(), chrom4.InsAlcLst.begin() + crossPoint);
//    int a = crossPoint; int b = crossPoint 不需要之间用k 或 crossPoint - 已修改xy
    for (int k = crossPoint; k < comConst.NumOfTasks; ++k) {
        if (chrom3.InsAlcLst[k] == -1) {
            vector<int> STI;    //个体的某一段中相同实例的任务编号-G
            vector<int> SII;    //个体的某一段中实例类型与某值相同的任务编号-G
            for (int i = k; i < comConst.NumOfTasks; ++i) {
                if (chrom2.InsAlcLst[k] == chrom2.InsAlcLst[i])
                    STI.push_back(i);//存放的是k及其以后的任务中与k有着相同实例的任务编号
            }
            for (int j = 0; j < crossPoint; ++j) {//k号任务的类型与前面存在的半段实例编号的类型是否相同
                if (chrom3.TypeLst[chrom3.InsAlcLst[j]] == chrom2.TypeLst[chrom2.InsAlcLst[k]])
                    SII.push_back(j);//存放的是任务编号，分配实例的类型与K的相同的任务编号
            }
            if (SII.size() != 0 && rand() % 1000 / 1000.0 < mrgPro ) {
                int rnd = rand() % SII.size();
                for (int i = 0; i < STI.size(); ++i) chrom3.InsAlcLst[STI[i]] = chrom3.InsAlcLst[SII[rnd]];
            } else {
                cv1 = cv1 + 1;
                for (int i = 0; i < STI.size(); ++i)  chrom3.InsAlcLst[STI[i]] = cv1;
                chrom3.TypeLst[cv1] = chrom2.TypeLst[chrom2.InsAlcLst[k]];
            }
        }
        if (chrom4.InsAlcLst[k] == -1) {
            vector<int> STI1, SII1;
            for (int i = k; i < comConst.NumOfTasks; ++i) {
                if (chrom1.InsAlcLst[k] == chrom1.InsAlcLst[i])   STI1.push_back(i);
            }
            for (int j = 0; j < crossPoint; ++j) {
                if (chrom4.TypeLst[chrom4.InsAlcLst[j]] == chrom1.TypeLst[chrom1.InsAlcLst[k]])   SII1.push_back(j);
            }
            if (SII1.size() != 0 && rand() % 1000 / 1000.0 < mrgPro) {
                int indexType = rand() % SII1.size();
                for (int i = 0; i < STI1.size(); ++i)  chrom4.InsAlcLst[STI1[i]] = chrom4.InsAlcLst[SII1[indexType]];
            } else {
                cv2 = cv2 + 1;
                for (int i = 0; i < STI1.size(); ++i)  chrom4.InsAlcLst[STI1[i]] = cv2;
                chrom4.TypeLst[cv2] = chrom1.TypeLst[chrom1.InsAlcLst[k]];
            }
        }
    }
    for (int k = cv1 + 1; k < comConst.NumOfTasks; ++k)  chrom3.TypeLst[k] = -1; //大于实例编码最大值的对应虚拟机类型变为-1，即未使用的实例的类型赋为-1
    for (int k = cv2 + 1; k < comConst.NumOfTasks; ++k)  chrom4.TypeLst[k] = -1; //大于实例编码最大值的对应虚拟机类型变为-1，即未使用的实例的类型赋为-1
    chrom1 = chrom3; chrom2 = chrom4;
}

void CrossoverWight(chromosome &chrom1, chromosome &chrom2){
    /* Simulated Binary Crossover (SBX，模拟二进制交叉) for MOELS */
    double Ri;
    double Ui = double(rand()%100)/100;
    if(Ui < 0.5){
        Ri = pow(2*Ui, 1/(10 + 1));       //原文中把交叉分布指数设为10
    } else{
        Ri = pow(1/(2*(1-Ui)), 1/(10 + 1));
    }
    chrom1.Weight = 0.5*((1 + Ri) * chrom1.Weight + (1 - Ri) * chrom2.Weight);
    chrom2.Weight = 0.5*((1 - Ri) * chrom1.Weight + (1 + Ri) * chrom2.Weight);
    //判断使其在定义域内
    if(chrom1.Weight > 1)
        chrom1.Weight = 1;
    else if(chrom1.Weight < 0)
        chrom1.Weight = 0;

    if(chrom2.Weight > 1)
        chrom2.Weight = 1;
    else if(chrom2.Weight < 0)
        chrom2.Weight = 0;
}

void CrossoverIns(chromosome &chrom1, chromosome &chrom2){
    chromosome chrom3, chrom4;
    initChrom(chrom3);
    initChrom(chrom4);
    chrom3.TskSchLst =chrom1.TskSchLst;
    chrom4.TskSchLst =chrom2.TskSchLst;
    chrom3.TypeLst =chrom1.TypeLst;
    chrom4.TypeLst =chrom2.TypeLst;
    int crosspoint = rand()%(comConst.NumOfTasks - 1) + 1;
    for (int i = crosspoint; i < comConst.NumOfTasks; ++i) {
        chrom3.InsAlcLst[i] = chrom1.InsAlcLst[i];
        chrom4.InsAlcLst[i] = chrom2.InsAlcLst[i];
    }
    for (int i = 0; i <crosspoint; ++i) {
        chrom3.InsAlcLst[i] = chrom2.InsAlcLst[i];
        chrom4.InsAlcLst[i] = chrom1.InsAlcLst[i];
    }
    vector<int >SI1, SI2;
    for (int j = crosspoint; j < comConst.NumOfTasks; ++j) {
        for (int i = 0; i < crosspoint; ++i) {
            if (chrom1.InsAlcLst[j] == chrom2.InsAlcLst[i]){
//                    if (SI1.size() == 0){
//                        SI1.push_back(chrom1.VMAllocationList[j]);
//                    } else{
//                        for (int k = 0; k < SI1.size(); ++k) {
//                            if (chrom1.VMAllocationList[j] == SI1[k]){
//                                continue;
//                            } else{
//                                SI1.push_back(chrom1.VMAllocationList[j]);
//                            }
//                        }
//                    }
//                    for (int l = 0; l < SI1.size(); ++l) {
//                        outfile<<endl<<"SI1"<<SI1[l]+1;
//                    }
                SI1.push_back(chrom1.InsAlcLst[j]);
                if(rand()%2 == 0){
                    chrom3.TypeLst[chrom1.InsAlcLst[j]] = chrom1.TypeLst[chrom1.InsAlcLst[j]];
                } else {
                    chrom3.TypeLst[chrom1.InsAlcLst[j]] = chrom2.TypeLst[chrom1.InsAlcLst[j]];
                }
            }

            if (chrom2.InsAlcLst[j] == chrom1.InsAlcLst[i]){
//                    if (SI2.size() == 0){
//                        SI2.push_back(chrom2.VMAllocationList[j]);
//                    } else{
//                        for (int k = 0; k < SI2.size(); ++k) {
//                            if (chrom2.VMAllocationList[j] == SI2[k]){
//                                continue;
//                            } else{
//                                SI2.push_back(chrom2.VMAllocationList[j]);
//                            }
//                        }
//                    }
                SI2.push_back(chrom2.InsAlcLst[j]);
                int randomType2 = rand()%2;
                if(randomType2 == 0){
                    chrom4.TypeLst[chrom2.InsAlcLst[j]] = chrom1.TypeLst[chrom2.InsAlcLst[j]];
                }
                if (randomType2 == 1){
                    chrom4.TypeLst[chrom2.InsAlcLst[j]] = chrom2.TypeLst[chrom2.InsAlcLst[j]];
                }
            }
        }
    }
    chrom1 = chrom3;
    chrom2 = chrom4;
}

void CrossoverIns_EMS(chromosome &chrom1, chromosome &chrom2){
    chromosome chrom3 = chrom1, chrom4 = chrom2;
    chrom3.MakeSpan = 0; chrom3.Cost = 0; chrom4.MakeSpan = 0; chrom4.Cost = 0;
    int crosspoint = rand()%(comConst.NumOfTasks - 1) + 1; //设置交叉点；
    for (int i = 0; i <crosspoint; ++i) {
        //交换交叉点前的实例，即交叉点前儿子chrom3的实例分配方案来自于母亲chrom2，而女儿chrom4的实例分配方案来自于父亲chrom1
        chrom3.InsAlcLst[i] = chrom2.InsAlcLst[i];
        chrom4.InsAlcLst[i] = chrom1.InsAlcLst[i];
    }
    //调整实例类型
    for (int i = 0; i < crosspoint; ++i) {//调整儿子chrom3的实例类型
        bool flg1 = true;
        for (int j = crosspoint; j < comConst.NumOfTasks; ++j) {
            if ( chrom2.InsAlcLst[i] == chrom1.InsAlcLst[j]){//存在使用相同实例编号的情况
                if (chrom1.TypeLst[chrom1.InsAlcLst[j]] != chrom2.TypeLst[chrom2.InsAlcLst[i]]) {
                    //如果两者的实例类型不相同，则随机选择一个, 即50%来自于母亲chrom2,50%来自于父亲chrom1（保持不变）
                    if(rand()%2 == 0) {
                        if (chrom3.InsAlcLst[j] != chrom1.InsAlcLst[j]) cout << "EMC INS crossover is woring 1 !"<<endl;
                        chrom3.TypeLst[chrom1.InsAlcLst[j]] = chrom2.TypeLst[chrom2.InsAlcLst[i]]; //来自于母亲chrom2
                    } else {
                        chrom3.TypeLst[chrom1.InsAlcLst[j]] = chrom1.TypeLst[chrom1.InsAlcLst[j]]; //来自于父亲chrom1，尽管初始化为父亲的类型，但由于中间可能会改变，因此此句不能省略，
                    }
                }
                flg1 = false;
                break;
            }
        }
        if (flg1 == true) {   //不存在使用相同实例编号的情况
            if (chrom3.InsAlcLst[i] != chrom2.InsAlcLst[i]) cout << "EMC INS crossover is woring 3 !"<<endl;
            chrom3.TypeLst[chrom2.InsAlcLst[i]] = chrom2.TypeLst[chrom2.InsAlcLst[i]]; //直接遗传，即来与母亲的实例类型相同
            if (rand()%1000/1000.0 < inputParameter.MutationProb) {  //以小的概率变异实例类型；
                map<int,VMType>::iterator it = VMTypes.begin();
                advance(it, rand() % comConst.NumOfVMTypes);
                chrom3.TypeLst[chrom2.InsAlcLst[i]] = (*it).first;
            }
        }
    }
    for (int i = 0; i < crosspoint; ++i) {  //调整女儿chrom4的实例类型
        bool flg2 = true;
        for (int j = crosspoint; j < comConst.NumOfTasks; ++j) {
            if ( chrom1.InsAlcLst[i] == chrom2.InsAlcLst[j]){ //存在使用相同实例编号的情况
                if (chrom2.TypeLst[chrom2.InsAlcLst[j]] != chrom1.TypeLst[chrom1.InsAlcLst[i]]) {
                    //如果两者的实例类型不相同，则随机选择一个, 即50%来自于父亲chrom1,50%来自于母亲chrom2（保持不变）
                    if(rand()%2 == 0) {
                        if (chrom4.InsAlcLst[j] != chrom2.InsAlcLst[j]) cout << "EMC INS crossover is woring 4 !"<<endl;
                        chrom4.TypeLst[chrom2.InsAlcLst[j]] = chrom1.TypeLst[chrom1.InsAlcLst[i]]; //来自于父亲chrom1
                    }else {
                        chrom4.TypeLst[chrom2.InsAlcLst[j]] = chrom2.TypeLst[chrom2.InsAlcLst[j]]; //来自于母亲chrom2，因此此句可以省略，
                    }
                }
                flg2 = false;
                break;
            }
        }
        if (flg2 == true) {   //不存在使用相同实例编号的情况
            if (chrom4.InsAlcLst[i] != chrom1.InsAlcLst[i]) cout << "EMC INS crossover is woring 6 !"<<endl;
            chrom4.TypeLst[chrom1.InsAlcLst[i]] = chrom1.TypeLst[chrom1.InsAlcLst[i]]; //直接遗传，即来与母亲的实例类型相同
            if (rand()%1000/1000.0 < inputParameter.MutationProb) {  //以小的概率变异实例类型；
                map<int,VMType>::iterator it = VMTypes.begin();
                advance(it,rand() % comConst.NumOfVMTypes);
                chrom4.TypeLst[chrom1.InsAlcLst[i]] = (*it).first;
            }
        }
    }
    chrom1 = chrom3; chrom2 = chrom4;
}

void MutationWight(chromosome &chrom1){
    /*Polynomial Mutation(PM，多项式变异)*/
    double Ri;
    double Ui = double(rand()%100)/100;
    if(Ui < 0.5){
        Ri = pow(2*Ui, (double)1/(20 + 1)) - 1;      // 原文中把变异分布指数设为20
    } else{
        Ri = 1 - pow(2*(1-Ui), (double)1/(20 + 1));  // 原来的"pow(1-(2*(1-Ui)), (double)1/(20 + 1))"有问题，修修正了-xy
    }
    chrom1.Weight += Ri;
    /*判断使其在定义域内*/
    if(chrom1.Weight > 1)
        chrom1.Weight = 1;
    else if(chrom1.Weight < 0)
        chrom1.Weight = 0;
}

void MutationTaskOrder(chromosome &chrom1){
    int index = rand()%comConst.NumOfTasks;
    int taskID = chrom1.TskSchLst[index];
    int str = index - 1, end = index + 1;
    while (str > -1 && (find(Tasks[taskID].parents.begin(), Tasks[taskID].parents.end(), chrom1.TskSchLst[str]) == Tasks[taskID].parents.end()))
        str--;
    str++;
    while (end < comConst.NumOfTasks && (find(Tasks[taskID].children.begin(), Tasks[taskID].children.end(), chrom1.TskSchLst[end]) == Tasks[taskID].children.end()))
        end++;
    if (end - str <= 1) {
        return;
    }
    int insertPoint = rand() % (end - str) + str;
    while (insertPoint ==index) {
        insertPoint = rand() % (end - str) + str;
    }
    if (insertPoint < index) {
        for (int i = index; i > insertPoint; i--) {
            chrom1.TskSchLst[i] = chrom1.TskSchLst[i - 1];
        }
        chrom1.TskSchLst[insertPoint] = taskID;
    } else {
        for (int i = index; i < insertPoint; i++) {
            chrom1.TskSchLst[i] = chrom1.TskSchLst[i + 1];
        }
        chrom1.TskSchLst[insertPoint] = taskID;
    }
}

void mutation_EMS(chromosome &chrom1){
    int rnd = rand() % 3;
    if (rnd == 0) {
        MutationTaskOrder(chrom1); //任务调度顺序变异
    } else if (rnd == 1) {         //虚拟机实例变异
        for (int i = 0; i < comConst.NumOfTasks; ++i) {
            if (rand()%1000/1000.0 < inputParameter.MutationProb) {
                chrom1.InsAlcLst[i] == rand()%comConst.NumOfTasks;
            }
        }
    } else {                       //虚拟机类型变异
        for (int j = 0; j < comConst.NumOfTasks; ++j) {
            if (rand()%1000/1000.0 < inputParameter.MutationProb) {
                map<int,VMType>::iterator it = VMTypes.begin();
                advance(it,rand() % comConst.NumOfVMTypes);
                chrom1.TypeLst[j] = (*it).first;
            }
        }
    }
}

void MutationIns(chromosome &chrom1){
    // /*for RSGA1 */
    int mtnId = rand() % (comConst.NumOfTasks-1) + 1; //由于第1个任务（0号任务）的实例编码固定为0，故从第2到最后一个任务中随机选择一个进行变异，
    double mrgPrb = 0.8;
    chromosome chrom2 = chrom1;
    int alpha = *max_element(chrom1.InsAlcLst.begin(), chrom1.InsAlcLst.begin()+mtnId) + 1; //找出变异点mtnId前实例编码的最大值(alpha)
//    int alpha1 = -1; //测试用-xy
//    for (int i = 0; i < mtnId; ++i) { //找到变异点前实例编码的最大值(alpha)
//        if (chrom1.InsAlcLst[i] > alpha1)   alpha1 = chrom1.InsAlcLst[i];
//    }
//    alpha1 = alpha1 + 1;
//    if (alpha != alpha1) {
//        cout<<endl<<"alpha is wrong"; exit(2);
//    }
    //随机生成0～1之间的值
    vector<int> SII;    //能处理变异点任务的所有虚拟机实例编号集合
    int MaxInsId = *max_element(chrom1.InsAlcLst.begin(),chrom1.InsAlcLst.end()); //虚拟机实例中最大的数
    for (int j = 0; j < MaxInsId; ++j) {  //找到除了InsAlcLst[mtnId]之外的能处理变异点任务的所有虚拟机实例的编号集合
        if(j == chrom1.InsAlcLst[mtnId]) continue;
        if(Tasks[mtnId].ElgRsc.find(chrom1.TypeLst[j]) != Tasks[mtnId].ElgRsc.end())  SII.push_back(j);
    }

    if(SII.size()!=0){   //S3
        if(rand() % 1000 / 1000.0 < mrgPrb) {    //S4:合并变异
            int insOfMrg = rand() % SII.size();  //取得要合并的实例编号
            if (SII[insOfMrg] > alpha) {        //S4.2:后合并操作 SII[insOfMrg]为k'
                chrom2.InsAlcLst[mtnId] = alpha; //S4.2.1
                for (int i = mtnId+1; i < comConst.NumOfTasks; ++i) {   //S4.2.1 把原来等于SII[insOfMrg]的值都改成alpha
                    if (chrom1.InsAlcLst[i] == SII[insOfMrg]) chrom2.InsAlcLst[i] = alpha;
                }
                chrom2.TypeLst[alpha] = chrom1.TypeLst[SII[insOfMrg]];
                if (chrom1.InsAlcLst[mtnId] == alpha) {    //S4.2.3
                    vector<int> tskOfEqlA;     //tskOfEqlA为所有虚拟机实例等于alpha的任务下标集合
                    for (int i = mtnId + 1; i < comConst.NumOfTasks; ++i) {
                        if (chrom1.InsAlcLst[i] == alpha) tskOfEqlA.push_back(i);
                    }
                    if (tskOfEqlA.size() != 0) {   //S4.2.4
                        int beta = -1;                //用来记录mtnId与tskOfEqlA[0]之间的最大值
                        for (int i = mtnId + 1 ; i < tskOfEqlA[0]; ++i) {
                            if (chrom1.InsAlcLst[i] > beta)  beta = chrom1.InsAlcLst[i];
                        }
                        if (beta >= SII[insOfMrg]) {     //beta >= k' S4.2.4.1
                            for (int i = mtnId + 1; i < comConst.NumOfTasks; ++i) {
                                if (chrom1.InsAlcLst[i] > SII[insOfMrg] && chrom1.InsAlcLst[i] <= beta){ //往前移，即-1
                                    chrom2.InsAlcLst[i] = chrom1.InsAlcLst[i] - 1;
                                }
                            }
                            for (int j = SII[insOfMrg]; j < beta; ++j)  chrom2.TypeLst[j] = chrom1.TypeLst[j+1]; //应该放到上面的循环外面！已改-xy
                            for (int i: tskOfEqlA)  chrom2.InsAlcLst[i] = beta;      //应该放到上面的循环外面！已改-xy
                            chrom2.TypeLst[beta] = chrom1.TypeLst[chrom1.InsAlcLst[mtnId]];  //应该放到上面的循环外面！已改-xy
                        }
                        if (alpha < beta && beta < SII[insOfMrg]) {//S4.2.4.2 大于beta小于SII[indexVmOfMrg]之间的虚拟机实例后移
                            for (int i = tskOfEqlA[0] + 1; i < comConst.NumOfTasks; ++i) {
                                if(chrom1.InsAlcLst[i] > beta && chrom1.InsAlcLst[i] < SII[insOfMrg]){ //实例后移,即+1
                                    chrom2.InsAlcLst[i]= chrom1.InsAlcLst[i] + 1;
                                }
                            }
                            for (int k = beta + 2; k <= SII[insOfMrg]; ++k)  chrom2.TypeLst[k]=chrom1.TypeLst[k-1]; //实例后移后修改相应实例的虚拟机类型
                            for (int i: tskOfEqlA)  chrom2.InsAlcLst[i] = beta + 1;
                            chrom2.TypeLst[beta + 1] = chrom1.TypeLst[chrom1.InsAlcLst[mtnId]];
                        }
                        if (beta < alpha) {       //S4.2.4.3  大于等于alpha小于SII[index]之间的虚拟机实例后移,即+1
                            for (int i = mtnId + 1; i < comConst.NumOfTasks; ++i) {
                                if(chrom1.InsAlcLst[i] >= alpha && chrom1.InsAlcLst[i] < SII[insOfMrg])
                                    chrom2.InsAlcLst[i] = chrom1.InsAlcLst[i] + 1;
                            }
                            for (int j = alpha + 1; j <= SII[insOfMrg]; ++j)  chrom2.TypeLst[j]=chrom1.TypeLst[j-1];
                        }
                    }
                    if (tskOfEqlA.size() == 0) {  //S4.2.5 大于SII[indexVM]的虚拟机实例前移,即-1
                        for (int i = mtnId + 1; i < comConst.NumOfTasks; ++i) {
                            if(chrom1.InsAlcLst[i] > SII[insOfMrg]) chrom2.InsAlcLst[i] = chrom1.InsAlcLst[i]-1;
                        }
                        for (int j = SII[insOfMrg]; j < MaxInsId; ++j)  chrom2.TypeLst[j]=chrom1.TypeLst[j+1];
                        chrom2.TypeLst[MaxInsId]= -1 ;
                    }
                }
                if (chrom1.InsAlcLst[mtnId] != alpha) { //S4.2.4.3 大于等于alpha小于SII[indexVm]的虚拟机实例后移,即+1
                    for (int i = mtnId + 1; i < comConst.NumOfTasks; ++i) {
                        if(chrom1.InsAlcLst[i] >= alpha && chrom1.InsAlcLst[i] < SII[insOfMrg])
                            chrom2.InsAlcLst[i] = chrom1.InsAlcLst[i] + 1;
                    }
                    for (int j = alpha + 1; j <= SII[insOfMrg]; ++j)  chrom2.TypeLst[j]=chrom1.TypeLst[j-1];
                }
            }
            if (SII[insOfMrg] <= alpha) {    //S4.3前合并操作
                chrom2.InsAlcLst[mtnId] = SII[insOfMrg]; // S4.3.1
                if(chrom1.InsAlcLst[mtnId] == alpha){    // S4.3.2
                    vector<int > tskOfEqlA;
                    for (int i = mtnId+1; i < comConst.NumOfTasks; ++i) {
                        if(chrom1.InsAlcLst[i] == alpha) tskOfEqlA.push_back(i);
                    }
                    if(tskOfEqlA.size()!=0){//S4.3.3
                        int beta = -1;   //找到最大的, 即记录mtnId与tskOfEqlA[0]之间的最大值
                        for (int i = mtnId + 1; i < tskOfEqlA[0]; ++i) {
                            if (chrom1.InsAlcLst[i] > beta) beta = chrom1.InsAlcLst[i];
                        }
                        if(beta > alpha){ //S4.3.4 大于alpha小于等于beta的虚拟机实例编号向前移1位,即-1
                            for (int i = mtnId + 1; i < comConst.NumOfTasks; ++i) {
                                if(chrom1.InsAlcLst[i] > alpha &&chrom1.InsAlcLst[i] <= beta){
                                    chrom2.InsAlcLst[i] = chrom1.InsAlcLst[i]-1;
                                }
                            }
                            for (int j = alpha; j < beta; ++j) chrom2.TypeLst[j]=chrom1.TypeLst[j + 1];
                            for (int i: tskOfEqlA)  chrom2.InsAlcLst[i] = beta; //-w
                            chrom2.TypeLst[beta] = chrom1.TypeLst[chrom1.InsAlcLst[mtnId]];
                        }
                    }
                    if(tskOfEqlA.size() == 0){    //S4.3.5 所有大于alpha的虚拟机实例编号向前移1位
                        for (int i = mtnId + 1; i < comConst.NumOfTasks; ++i) {
                            if (chrom1.InsAlcLst[i] > alpha)  chrom2.InsAlcLst[i] = chrom1.InsAlcLst[i] - 1;
                        }
                        for (int j = alpha ; j < MaxInsId; ++j)  chrom2.TypeLst[j] = chrom1.TypeLst[j+1];
                        chrom2.TypeLst[MaxInsId] = -1;
                    }
                }
            }
        } else {
            SplitMutation(chrom1,chrom2,mtnId,alpha,MaxInsId);
        }
    }else{
        SplitMutation(chrom1,chrom2,mtnId,alpha, MaxInsId);
    }
    chrom1=chrom2;
}

void SplitMutation(chromosome &chrom1, chromosome &chrom2, int mtnId, int alpha, int MaxInsId){
    chrom2.InsAlcLst[mtnId] = alpha;
    int rnd = rand() % Tasks[mtnId].ElgRsc.size();
    set<int>::iterator it = Tasks[mtnId].ElgRsc.begin();
    advance(it,rnd);
    chrom2.TypeLst[alpha] = *it;
    if (chrom1.InsAlcLst[mtnId] < alpha){  //由于chrom1.InsAlcLst[mutationIndex] <= alpha，也即chrom1.InsAlcLst[mutationIndex] != alpha
        for (int i = mtnId + 1; i < comConst.NumOfTasks; ++i) {//虚拟机实例编号大于等于alpha的编号往后移动1位//-w从0开始，-不用从0开始！-xy；
            if(chrom1.InsAlcLst[i] >= alpha )
                chrom2.InsAlcLst[i] = chrom1.InsAlcLst[i] + 1;
        }
        for (int j = alpha + 1; j <= MaxInsId + 1; ++j) {
            chrom2.TypeLst[j] = chrom1.TypeLst[j-1];
        }
    } else { //chrom1.InsAlcLst[mutationIndex] == alpha; S5.2
//        if (chrom1.InsAlcLst[mtnId] > alpha) { //测试用-xy
//            cout<<endl<<"code is wrong"; exit(2);
//        }
        vector<int> SAI;
        for (int i = mtnId + 1; i < comConst.NumOfTasks; ++i) {
            if(chrom1.InsAlcLst[i] == alpha)  SAI.push_back(i);
        }
        if(SAI.size() != 0){   //S5.3
            int max = *max_element(chrom1.InsAlcLst.begin(), chrom1.InsAlcLst.begin()+SAI[0]) + 1;//找到SAI[0]之前最大的虚拟机实例编号,并+1
            for (int j = 0; j < SAI.size(); ++j)  chrom2.InsAlcLst[SAI[j]] = max; //令所有虚拟机实例编号等于alpha重新编号为max
            chrom2.TypeLst[max] = chrom1.TypeLst[alpha];
            for (int j = SAI[0] + 1; j < comConst.NumOfTasks; ++j) {//虚拟机实例编号大于等于max的编号往后移动1位
                if(chrom1.InsAlcLst[j] >= max) chrom2.InsAlcLst[j] = chrom1.InsAlcLst[j] + 1;
            }
            for (int k = max + 1; k <= MaxInsId + 1; ++k)  chrom2.TypeLst[k] = chrom1.TypeLst[k-1];
        }
    }
}

void MutationType(chromosome &chrom1){
    int maxInsId = *max_element(chrom1.InsAlcLst.begin(),chrom1.InsAlcLst.end());
    int slcInsId = rand() % (maxInsId + 1); //从所有已使用的实例中随机选择一个实例
    vector<int> tskAsgToSlcIns;
    for (int i = 0; i < comConst.NumOfTasks; ++i) { //找出所有已分配给slcInsId的任务
        if (chrom1.InsAlcLst[i] == slcInsId) tskAsgToSlcIns.push_back(i);
    }
    set<int> avlVMTypes = VMTypeSet;
    for (int tsk: tskAsgToSlcIns) {  //找出所有能处理已分配给slcInsId的所有任务的虚拟机类型；
        for (int ty: Tasks[tsk].NoElgRsc) {
            avlVMTypes.erase(ty);
        }
    }
    if (avlVMTypes.size() == 0){ //测试用cs-xy
        cout<<endl<<"it is wrong in MutationType for avlVMTypes"; exit(2);
    }
    if (avlVMTypes.size() > 1) { //重新随机分配一个不同是实例类型
        set<int>::iterator it = avlVMTypes.begin();
        int rnd = rand() % avlVMTypes.size();
        advance(it, rnd);
        while (*it == chrom1.TypeLst[slcInsId]) {
            it = avlVMTypes.begin();
            rnd = rand() % avlVMTypes.size();
            advance(it, rnd);
        }
        chrom1.TypeLst[slcInsId] = *it;
    }
}

void MutationTask(chromosome &chrom1,int index){
    int taskID = chrom1.TskSchLst[index];
    int str = index - 1, end = index + 1;
    while (str > -1 && (find(Tasks[taskID].parents.begin(), Tasks[taskID].parents.end(), chrom1.TskSchLst[str]) ==
                        Tasks[taskID].parents.end()))
        str--;
    str++;

    while (end < comConst.NumOfTasks &&
           (find(Tasks[taskID].children.begin(), Tasks[taskID].children.end(), chrom1.TskSchLst[end]) ==
            Tasks[taskID].children.end()))
        end++;

    if (end - str <= 1) {
        return;
    }
    int insertPoint = rand() % (end - str) + str;
    while (insertPoint ==index) {
        insertPoint = rand() % (end - str) + str;
    }
    if (insertPoint < index) {
        for (int i = index; i > insertPoint; i--) {
            chrom1.TskSchLst[i] = chrom1.TskSchLst[i - 1];
        }
        chrom1.TskSchLst[insertPoint] = taskID;
    } else {
        for (int i = index; i < insertPoint; i++) {
            chrom1.TskSchLst[i] = chrom1.TskSchLst[i + 1];
        }
        chrom1.TskSchLst[insertPoint] = taskID;
    }
    chrom1.TypeLst = chrom1.TypeLst;
    chrom1.InsAlcLst = chrom1.InsAlcLst;
}

void EMS_DecodeChrom(chromosome& ch){//没有对makespan初始化，也没有计算fitness-xy
    vector<set<double> > ITL;
    for (int j = 0; j < *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end())+1; ++j) {
        set<double> a;
        a.insert(VMTypes[ch.TypeLst[j]].BootTime);//正向解码时初始化为roottime；反向解码时初始为0，但在计算makespan时每个资源实例的最大完成时间上加roottims；-xy2
        a.insert(InfInt * 1.0);
        ITL.push_back(a);
    }
    //解码
    for(int i = 0; i < comConst.NumOfTasks; ++i){
        int TaskIndex = ch.TskSchLst[i];
        int RscIndex = ch.InsAlcLst[TaskIndex];
        double ReadyTime = 0;
        for (int ParentTask: Tasks[TaskIndex].parents) {
            int ParentRsc = ch.InsAlcLst[ParentTask];
            double TransferTime = 0;
            if(RscIndex != ParentRsc) {
                TransferTime = ParChildTranFileSizeSum[ParentTask][TaskIndex]  / (XY_MIN(VMTypes[ch.TypeLst[RscIndex]].bw, VMTypes[ch.TypeLst[ParentRsc]].bw));
            }
            double sum = ch.EndTime[ParentTask] + TransferTime;
            if (ReadyTime < sum) {
                ReadyTime = sum;
            }
        }
        set<double>::iterator pre  = ITL[RscIndex].begin();
        set<double>::iterator post = ITL[RscIndex].begin();
        ++post;
        double ExecutionTime = Tasks[TaskIndex].length / VMTypes[ch.TypeLst[RscIndex]].ps;
        while(post != ITL[RscIndex].end()) {
            if((*post - *pre) >= ExecutionTime && ReadyTime <= (*post)-ExecutionTime) {
                ch.StartTime[TaskIndex] = XY_MAX(*pre, ReadyTime);
                break;
            } else {
                ++pre;
                ++pre;
                ++post;
                ++post;
            }
        }
        ch.EndTime[TaskIndex] = ch.StartTime[TaskIndex] + ExecutionTime;
        if (ch.MakeSpan < ch.EndTime[TaskIndex]) {
            ch.MakeSpan = ch.EndTime[TaskIndex];
        }
        //{update ITL}
        if(ITL[RscIndex].find(ch.StartTime[TaskIndex]) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(ch.StartTime[TaskIndex]);
        } else {
            ITL[RscIndex].insert(ch.StartTime[TaskIndex]);
        }

        if(ITL[RscIndex].find(ch.EndTime[TaskIndex]) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(ch.EndTime[TaskIndex]);
        } else {
            ITL[RscIndex].insert(ch.EndTime[TaskIndex]);
        }
    }
    //计算成本
    double cost =0;
    for (int k = 0; k <= *max_element(ch.InsAlcLst.begin(),ch.InsAlcLst.end()) ; ++k){
        double Hrt = 10000000000;
        double Rnt = 0;
        vector<int> HandleTask;

        for (int i = 0; i < comConst.NumOfTasks; ++i) {         //找到相同实例的任务
            if(ch.InsAlcLst[i] == k) HandleTask.push_back(i);
        }
        if(HandleTask.size() == 0){
            continue;
        } else {
            for (int i = 0; i < HandleTask.size(); ++i) {       //虚拟机开始租用时间的计算（Hrt）
                double hrt = ch.StartTime[HandleTask[i]] - VMTypes[ch.TypeLst[k]].BootTime;
                if (hrt < Hrt) Hrt = hrt;
            }
            for (int i = 0; i < HandleTask.size(); ++i) {       //虚拟机归还时间计算（Rnt）
                //Rnt时间的计算是该任务完成时间
                if (Rnt < ch.EndTime[HandleTask[i]]) Rnt = ch.EndTime[HandleTask[i]];
            }
            double Lt = Rnt - Hrt;                              //最晚归还时间减去最早开始时间为虚拟机的租用时间
            double SpanTime = XY_MAX(Lt, VMTypes[ch.TypeLst[k]].MinLeaseTime);
            cost += VMTypes[ch.TypeLst[k]].FixedCost + (ceil((SpanTime) / VMTypes[ch.TypeLst[k]].ChargeUnit)) * VMTypes[ch.TypeLst[k]].UnitCost;
        }
    }
    ch.Cost = cost;
}
