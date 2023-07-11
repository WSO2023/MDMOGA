//
// Created by xieyi on 2023/3/29.
//

#ifndef FRAME_CHROPER_H
#define FRAME_CHROPER_H
#include "common.h"
void initChrom(chromosome &chrom);
void GnrRscAlc_Random(chromosome& chrom);
//void GnrRscAlc_BenefitRatio(chromosome &chrom);
//void GnrRscAlc_BenefitRatio_xy(chromosome &chrom, double refCst, double refMs);
void GnrRscAlc_BenefitRatio_MOELS(chromosome &chrom, double refCst, double refMs);
void GnrRscAlc_BenefitRatioNew(chromosome &chrom);
void ClcAvrTskExcTime(vector<double>& w);
void ClcAvrCmnTime_AB(vector<vector<double>>& c);
void ClcAvrCmnTime(vector<vector<double>>& c);
void Calculate_Rank_t(vector<double>& RankList, vector<double>& w, vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& RankList, vector<double>& w, vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& RankList, vector<double>& w);
vector<chromosome> GnrChrsByGRPHEFT(double estMinBdg, double estMaxBdg, vector<int>& TskSchLstByRnk_b);

void Modified_HEFT_XY(chromosome &chrom, vector<Instance> &Y, map<int,int> &NumOfInsMrg);
int Modified_HEFT_XY_1(chromosome &chrom, vector<Instance> &Y, map<int,int> &NumOfInsMrg);
int Modified_HEFT_XY_2(chromosome &chrom, vector<Instance> &Y, map<int,int> &NumOfInsMrg);
int Modified_HEFT_XY_Imp(chromosome &chrom, vector<Instance> &Y, map<int,int> &NumOfInsMrg);
void Modified_HEFT_f(chromosome &chrom, vector<VMType> &Y);

void GenericParticleByRank(chromosome& chrom, vector<double>& Rank);
void GnrTypeAlcPartByRnd(chromosome& chrom);
void RepairTypeAlcPart(chromosome& ch);
int FindNearestRscId(int TaskId, double value);
void RepairPriorityAndGnrSchOrd(chromosome& chrom);
void UpdateParticle(chromosome &ch,chromosome &Pbest, chromosome &Gbest, double &runtime, double &SchTime);
void UpdateParticle_gui(chromosome &chrom, chromosome &Pbest, chromosome &Gbest, int t);
double FindTaskReadyTime(int tsk, int ins, chromosome &chrom);
void  DecodeChrom(chromosome& ch);
void  DecodeChrom_NI(chromosome& ch);

double FindIdleTimeSlot(set<double>& ITLofRscId, double& ExeTime, double& ReadyTime);
void UpdateITL(set<double>& ITLofRscId,double& StartTime,double& EndTime);

void DecodeParticle(chromosome& chrom, vector<double>& tskLFT);
void GnrTwoDimCoding(chromosome &chrom);
bool GnrOneDimCoding(chromosome &chrom);
void RecodeInsAndType(chromosome &chrom);
int IFBDI(chromosome& ch);
int LBCAI(chromosome& chrom);
int RPI_ShfTsk1(chromosome& ch);
int RPI_ShfTsk2(chromosome& ch);
int RPI_SplIns1(chromosome& ch);
int RPI_SplIns2(chromosome& ch);

void InitProModelOfOdr1(vector<vector<double>>& PMO, vector<int>& NumOfAncestors, vector<int>& NumOfNonDescendants, vector<double>& Rank_b);
void InitProModelOfOdr2(vector<vector<double>>& PMO);
void InitProModelOfGrp(vector<vector<double>>& PMG);
void InitProModelOfRsc1(vector<map<int, double>>& PMR);
void InitProModelOfRsc2(vector<map<int, double>>& PMR);
chromosome GnrTskLstOfChr1_prp(vector<vector<double>>& PMO, vector<double>& eta_TSO);
chromosome GnrTskLstOfChr2_prp(vector<vector<double>>& PMO, vector<double>& eta_TSO);
void GnrInsLstOfChr(chromosome& chrom, vector<vector<double>>& PMG);
void GnrTypeLstOfChr1(chromosome& chrom, vector<map<int, double>>& PMR);
void GnrTypeLstOfChr2(chromosome& chrom, vector<map<int, double>>& PMR);
void UpdatePMO1(vector<vector<double>>& PMO, vector<chromosome>& elitPop);
void UpdatePMO2(vector<vector<double>>& PMO, vector<chromosome>& elitPop);
void UpdatePMG(vector<vector<double>>& PMG, vector<chromosome>& elitPop);
void UpdatePMR1(vector<map<int, double>>& PMR, vector<chromosome>& elitPop);
void UpdatePMR2(vector<map<int, double>>& PMR, vector<chromosome>& elitPop);

//void FindIdleTimeSlot(vector<set<double>>& ITL,int& RscIndex,double& StartTime,double& ExeTime,double& ReadyTime);
//void UpdateITL(vector<set<double>>& ITL,int& RscIndex,double& StartTime,double& EndTime);

#endif //FRAME_CHROPER_H
