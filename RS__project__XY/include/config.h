//
// Created by 11379 on 2020-10-30.
//

#ifndef FRAME_CONFIG_H
#define FRAME_CONFIG_H
#include "common.h"

void DeleteFirstLineInFile(string fileName);
int ReadID(string id);
void ReadFile(string XmlFile,string RscAlcFile);
void ClearALL();
void TaskCombine();
void HPSO_ConfigInputParameter();
void NEDA_ConfigInputParameter();
void RSGA1_ConfigInputParameter();
void MDMOGA_ConfigInputParameter();
void MOELS_ConfigInputParameter();
void EMS_ConfigInputParameter();

#endif //FRAME_CONFIG_H
