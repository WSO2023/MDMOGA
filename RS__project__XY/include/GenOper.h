//
// Created by xieyi on 2023/3/29.
//

#ifndef FRAME_GENOPER_H
#include "common.h"
#define FRAME_GENOPER_H

int ChromSelect(vector<double>& A);
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A);
bool Cross_TwoDim(chromosome &chrom1 ,chromosome &chrom2, chromosome &newChrom);
bool Mutation_Group(chromosome &chrom);
void FndPosInsertTsk_xy(list<list<int>> &Code, int TaskId, int k);
bool Mutation_TskOrd(chromosome &chrom);
void Mutation_Type(chromosome &chrom);

vector<int> GenerateAChrom_order(vector<double>& rank);
int chooseNextTask(vector<double> rank, vector<int> RT);
void CrossoverTaskOrder(chromosome &chrom1, chromosome &chrom2, bool flag);
void CrossoverTasklist(chromosome &chrom1, chromosome &chrom2);
void CrossoverInsAndType(chromosome &chrom1, chromosome &chrom2);
void CrossoverWight(chromosome &chrom1, chromosome &chrom2);
void MutationWight(chromosome &chrom1);
void CrossoverIns(chromosome &chrom1, chromosome &chrom2);
void CrossoverIns_EMS(chromosome &chrom1, chromosome &chrom2);
void mutation_EMS(chromosome &chrom1);
void MutationTaskOrder(chromosome &chrom1);
void MutationTask(chromosome &chrom1,int index);
void MutationIns(chromosome &chrom1);
void SplitMutation(chromosome &chrom1,chromosome &chrom2,int mutationIndex,int alpha, int MaxInsId);
void MutationType(chromosome &chrom1);

//void RepairPriorityAndGnrSchOrd(chromosome& chrom);
//void UpdateParticle(chromosome &ch,chromosome &Pbest, chromosome &Gbest, double &runtime, double &SchTime);
void EMS_DecodeChrom(chromosome& ch);


#endif //FRAME_GENOPER_H
