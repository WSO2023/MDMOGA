//
// Created by 11379 on 2020-10-30.
//

#ifndef FRAME_TOOLS_H
#include "common.h"
#define FRAME_TOOLS_H

void WriteOrgResultToTxt (string AlgName, string MdlName, pair<double, vector<chromosome>>& result);
void WriteOrgResultToTxt (string AlgName, string MdlName, int iterations, pair<double, vector<chromosome>>& result);
map<string, double> CalucateHV(const map<string, pair<double,vector<chromosome>>>& result);
map<string, double> CalucateHV(const map<string, vector<pair<double,double>>>& result);
map<string, double> CalucateIGD(const map<string, pair<double,vector<chromosome>>>& result);
map<string, double> CalucateIGD(const map<string, vector<pair<double,double>>>& result);
map<string, map<string, double>> CalucateCM(const map<string, pair<double,vector<chromosome>>>& result);
map<string, map<string, double>> CalucateCM(const map<string, vector<pair<double,double>>>& result);

bool sortPopOnDistanceByDesc(chromosome& a, chromosome& b);
bool SortValueByDescend(pair<int,double>& a, pair<int,double>& b);
bool SortValueByAscend(pair<double,double>& a, pair<double,double>& b);
void IndexSortByValueAsc(vector<int>& ind, vector<double>& fitness);
void IndexSortByValueOnAsc(vector<int>& ind, vector<int>& value);
void IndexSortByValueDesc(vector<int>& ind, vector<double>& fitness);
void IndexSortByValueDesc(vector<int>& ind, vector<int>& tskLst, vector<double>& value);
bool sortPopByFitAsc(chromosome& a, chromosome& b);
bool sortPopByCostDesc(chromosome& a, chromosome& b);
bool sortPopByCostAsc(chromosome& a, chromosome& b);
bool sortPopByMakespanDesc(chromosome& a, chromosome& b);
bool sortPopByMakespanAsc(chromosome& a, chromosome& b);
bool sortPopByMsCstAsc_fstMs(chromosome& a, chromosome& b);
void CalculateLevelList();
void CalculateDescendants();
void CalculateAncestors();
vector<int> GenerateTopoOrder_Ran();
vector<int> GenerateTopoOrder_lev();
vector<int> GenerateTopoOrder_lev_bt(vector<double >& rank_bt);
vector<double> CalculateDeadline(double Deadline);
vector<double> calculateTskLFT(double Deadline);
double CheapSchedule();
double CheapSchedule(chromosome &ch);
double CheapSchedule_xy();
double CheapSchedule_xy(chromosome& ch);
double CheapSchedule_MOELS();
double FastSchedule();
double FastSchedule(chromosome &ch);
double FastSchedule_xy(vector<double>& rank_b);
double FastSchedule_xy(vector<double>& rank_b, chromosome& chrom);
double FastSchedule_MOELS(vector<double>& rank_b);
bool BetterParticle(chromosome& ch1, chromosome& ch2, double& Deadline);
chromosome BetterPartcle(chromosome& p1, chromosome& p2, double& TotalDeadline);
double estMaxMakespan_HPSO();
double estMinMakespan_HPSO();
void CalucateAdjMatrix();
void UTtoRT(set<int>& UT, vector<int>& RT, vector<vector<int>>& linkMatrix);
double RandomDouble(int start, int end);

void selection_Tournament(vector<chromosome>& population, int& parent_1, int& parent_2);
int selectChrom_Tournament(int& NumOfChromPerPop);

vector<chromosome> ParetoFront(vector<chromosome>& population);
void SelectParetoFront(vector<chromosome>& population, vector<chromosome>& NDVec);
vector<chromosome> GnrParetoFront(set<chromosome>& population);
vector<pair<double, double>> GnrParetoFront(set<pair<double, double>>& solutions);
vector<chromosome> GnrParetoFront(list<chromosome>& pop);
vector<chromosome> GnrParetoFront(vector<chromosome>& pop);
void sortParetoFrontByCrowdingDistance (vector<chromosome>& PFVec);
double clcCstOfInsLeased (double strTimeOfUse, double endTimeOfUse, int VMType);



#endif //FRAME_TOOLS_H
