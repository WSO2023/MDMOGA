#include <iostream>
#include <config.h>
#include <fstream>
#include "tools.h"
#include "MDMOGA.h"
#include "MOELS.h"
#include "GRPHEFT.h"
#include "EMS.h"
#include "MOHEFT.h"
#include "HMCMB.h"

using namespace std;

int main() {
    srand((int) time(0));
    //{set the runtime (termination) time of the algorithm}
    map<string, double> SchTime;
    SchTime["Montage25_0.4"] = 1.166;
    SchTime["Montage25_0.7"] = 1.553;
    SchTime["Montage25_1.0"] = 2.953;
    SchTime["Montage50_0.4"] = 2.999;
    SchTime["Montage50_0.7"] = 4.400;
    SchTime["Montage50_1.0"] = 19.229;
    SchTime["Montage100_0.4"] =9.100;
    SchTime["Montage100_0.7"] =15.142;
    SchTime["Montage100_1.0"] =99.766;
    SchTime["CyberShake30_0.4"] = 1.151;
    SchTime["CyberShake30_0.7"] = 2.279;
    SchTime["CyberShake30_1.0"] = 3.777;
    SchTime["CyberShake50_0.4"] = 3.021;
    SchTime["CyberShake50_0.7"] = 5.694;
    SchTime["CyberShake50_1.0"] = 18.747;
    SchTime["CyberShake100_0.4"] =9.627;
    SchTime["CyberShake100_0.7"] =15.455;
    SchTime["CyberShake100_1.0"] =36.899;
    SchTime["Epigenomics24_0.4"] = 0.880;
    SchTime["Epigenomics24_0.7"] = 1.202;
    SchTime["Epigenomics24_1.0"] = 3.058;
    SchTime["Epigenomics47_0.4"] = 2.217;
    SchTime["Epigenomics47_0.7"] = 3.158;
    SchTime["Epigenomics47_1.0"] = 38.069;
    SchTime["Epigenomics100_0.4"] =4.199;
    SchTime["Epigenomics100_0.7"] =7.431;
    SchTime["Epigenomics100_1.0"] =87.597;
    SchTime["Inspiral30_0.4"] = 1.257;
    SchTime["Inspiral30_0.7"] = 1.689;
    SchTime["Inspiral30_1.0"] = 12.641;
    SchTime["Inspiral50_0.4"] = 3.040;
    SchTime["Inspiral50_0.7"] = 4.092;
    SchTime["Inspiral50_1.0"] = 42.18;
    SchTime["Inspiral100_0.4"] =7.067;
    SchTime["Inspiral100_0.7"] =9.676;
    SchTime["Inspiral100_1.0"] =98.897;
    SchTime["Sipht29_1.0"] =3.617;
    SchTime["Sipht58_1.0"] =17.72;
    SchTime["Sipht97_1.0"] =67.29;

    MinBdg["Montage_25.xml"] = 0.08;
    MinBdg["Montage_50.xml"] = 0.08;
    MinBdg["Montage_100.xml"] = 0.08;
    MinBdg["CyberShake_30.xml"] = 0.26;
    MinBdg["CyberShake_50.xml"] = 0.48;
    MinBdg["CyberShake_100.xml"] = 0.76;
    MinBdg["Epigenomics_24.xml"] = 6.3;
    MinBdg["Epigenomics_47.xml"] = 14.7;
    MinBdg["Epigenomics_100.xml"] = 140;
    MinBdg["Inspiral_30.xml"] = 2.8;
    MinBdg["Inspiral_50.xml"] = 5.3;
    MinBdg["Inspiral_100.xml"] = 9.2;
    MinBdg["Sipht_29.xml"] = 2;
    MinBdg["Sipht_58.xml"] = 4.3;
    MinBdg["Sipht_97.xml"] = 6.5;
    MaxBdg["Montage_25.xml"] = 4.68;
    MaxBdg["Montage_50.xml"] = 14.56;
    MaxBdg["Montage_100.xml"] = 32.24;
    MaxBdg["CyberShake_30.xml"] = 4.16;
    MaxBdg["CyberShake_50.xml"] = 4.68;
    MaxBdg["CyberShake_100.xml"] = 0.4;
    MaxBdg["Epigenomics_24.xml"] = 45.24;
    MaxBdg["Epigenomics_47.xml"] = 105.04;
    MaxBdg["Epigenomics_100.xml"] = 944.72;
    MaxBdg["Inspiral_30.xml"] = 21.84;
    MaxBdg["Inspiral_50.xml"] = 39;
    MaxBdg["Inspiral_100.xml"] = 68.12;
    MaxBdg["Sipht_29.xml"] = 36.33;
    MaxBdg["Sipht_58.xml"] = 38.56;
    MaxBdg["Sipht_97.xml"] = 47.4;

    //{clear "result"}
    ofstream outfile;
    outfile.open("../result.txt", ios::out);
    outfile.close();
    outfile.open("../result_SchTime.txt", ios::out);
    outfile.close();
    outfile.open("../result_EMS.txt", ios::out);
    outfile.close();
    outfile.open("../result_MOELS.txt", ios::out);
    outfile.close();
    outfile.open("../result_MDMOGA.txt", ios::out);
    outfile.close();
    outfile.open("../result_GRPHEFT.txt", ios::out);
    outfile.close();
    outfile.open("../result_GRPHEFT_Org.txt", ios::out);
    outfile.close();
    outfile.open("../result_MOHEFT.txt", ios::out);
    outfile.close();
    outfile.open("../result_HMCMB.txt", ios::out);
    outfile.close();
    outfile.open("../result_Alg.txt", ios::out);  //用于输出帕累托前沿
    outfile.close();
    int stop_target =0;
    string Model, NumOfTask, RscAvlRatio;
    bool flg1 = true, flg2 = true;
    do {
        string StrLine;
        ifstream iFile("../fileList.txt");
        if (!iFile) {
            cout << endl << "filelist open failed!\n";
            exit(0);
        }
        getline(iFile, StrLine);
        if (StrLine.size() < 1) {
            cout << endl << "Empty input file" << endl;
            exit(0);
        }
        iFile.close();
        string XmlFile, RscAlcFile;
        istringstream is(StrLine);
        is >> Model >> NumOfTask >> RscAvlRatio;  //NumOfTask, RscAvlRatio
        XmlFile = Model + "_" + NumOfTask + ".xml";
        RscAlcFile = NumOfTask + "_" + RscAvlRatio + "_0.txt";
        string MdlName =  Model + "_" + NumOfTask + "_" + RscAvlRatio;
        cout << endl << Model << " " << NumOfTask << " " << RscAvlRatio;
        map<string, pair<double,vector<chromosome>>> result;
        map<string, pair<double,vector<chromosome>>> tem_result;
        vector<chromosome> paretoFront;

        double EMS_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int EMS_Iteration = 0;
        paretoFront = runEMS(XmlFile, RscAlcFile, EMS_SchTime, EMS_Iteration);
        result["EMS"] = {EMS_SchTime, paretoFront};
        WriteOrgResultToTxt("EMS", MdlName, EMS_Iteration, result["EMS"]);
        ClearALL();

        double MOELS_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int MOELS_Iteration = 0;
        paretoFront = runMOELS(XmlFile, RscAlcFile, MOELS_SchTime, MOELS_Iteration);
        result["MOELS"] = {MOELS_SchTime, paretoFront};
        WriteOrgResultToTxt("MOELS", MdlName, MOELS_Iteration, result["MOELS"]);
        ClearALL();

//        double NEDA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        int NEDA_Iteration = 0;
//        paretoFront = runNEDA(XmlFile, RscAlcFile, NEDA_SchTime, NEDA_Iteration);
//        result["NEDA"] = {NEDA_SchTime, paretoFront};;
//        WriteOrgResultToTxt("NEDA", MdlName, NEDA_Iteration, result["NEDA"]);
//        ClearALL();

//        double RSGA1_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        int RSGA1_Iteration = 0;
//        paretoFront = runRSGA1(XmlFile, RscAlcFile, RSGA1_SchTime, RSGA1_Iteration);
//        result["RSGA1"] = {RSGA1_SchTime, paretoFront};
//        WriteOrgResultToTxt("RSGA1", MdlName, RSGA1_Iteration, result["RSGA1"]);
//        ClearALL();
//
        double MDMOGA_SchTime  = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
        int MDMOGA_Iteration = 0;
        paretoFront = runMDMOGA(XmlFile, RscAlcFile, MDMOGA_SchTime, MDMOGA_Iteration);
        result["MDMOGA"] = {MDMOGA_SchTime, paretoFront};
        WriteOrgResultToTxt("MDMOGA", MdlName, MDMOGA_Iteration, result["MDMOGA"]);
        ClearALL();
//
        double MOHEFT_SchTime = 0;
        paretoFront = runMOHEFT(XmlFile,RscAlcFile,MOHEFT_SchTime);
        result["MOHEFT"] = {MOHEFT_SchTime,paretoFront};
        WriteOrgResultToTxt("MOHEFT", MdlName, result["MOHEFT"]);
        ClearALL();
//
        double HMCMB_SchTime = 0;
        paretoFront = runHMCMB(XmlFile,RscAlcFile,HMCMB_SchTime);
        result["HMCMB"] = {HMCMB_SchTime,paretoFront };
        WriteOrgResultToTxt("HMCMB", MdlName, result["HMCMB"]);
        ClearALL();


        if (!result.empty()) { //由于HPSO, GRPHEFT的需要，从已获得的结果中估算最大最小预算和截止时间
            MinBdg[XmlFile] = 9999999999, MaxBdg[XmlFile] = -1;
            MinDT[XmlFile] = 9999999999, MaxDT[XmlFile] = -1;
            for (const auto & alg: result) { //由于result中是按完工时间升序（成本降序）排序的，因此只需考虑第一个和最后一个值即可
                if (alg.second.second[0].MakeSpan + PrecisionValue < MinDT[XmlFile]) MinDT[XmlFile] = alg.second.second[0].MakeSpan;
                if (alg.second.second[alg.second.second.size()-1].MakeSpan - PrecisionValue > MaxDT[XmlFile]) MaxDT[XmlFile] = alg.second.second[alg.second.second.size()-1].MakeSpan;
                if (alg.second.second[0].Cost - PrecisionValue > MaxBdg[XmlFile]) MaxBdg[XmlFile] = alg.second.second[0].Cost;
                if (alg.second.second[alg.second.second.size()-1].Cost + PrecisionValue < MinBdg[XmlFile]) MinBdg[XmlFile] = alg.second.second[alg.second.second.size()-1].Cost;
            }
        }

        double GRPHEFT_SchTime = 0;
        paretoFront = runGRPHEFT(XmlFile, RscAlcFile, GRPHEFT_SchTime);
        result["GRPHEFT"] = {GRPHEFT_SchTime, paretoFront};
        WriteOrgResultToTxt("GRPHEFT", MdlName, result["GRPHEFT"]);
        ClearALL();

//        double HPSO_SchTime = SchTime[Model + NumOfTask + "_" + RscAvlRatio];
//        int HPSO_Iteration = 0;
//        result["HPSO"] = runHPSO(XmlFile, RscAlcFile, HPSO_SchTime, HPSO_Iteration);
//        WriteOrgResultToTxt("HPSO", MdlName, HPSO_SchTime, HPSO_Iteration, result["HPSO"]);
//        ClearALL();

        /*把算法的调度/运行时间写入到文件result_SchTime */;
        outfile.open("../result_SchTime.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(3);
        if (flg1) { //输出标题行
            outfile << "case ";
            for (auto slt: result) outfile << " " + slt.first;
            outfile << endl;
            flg1 = false;
        }
        outfile << Model + "_" + NumOfTask + "_" + RscAvlRatio + " ";
        for (auto slt: result) outfile << " " << slt.second.first;
        outfile << endl;
        outfile.close();
        string temmur = Model + NumOfTask + "_" + RscAvlRatio;

        if(stop_target !=4) {
            for (auto & tem_al : result) {
                if(tem_result[tem_al.first].second.empty()){
                    tem_result[tem_al.first] = tem_al.second;
                }else{
                    set<chromosome> mrgPop;
                    mrgPop.insert(tem_result[tem_al.first].second.begin(),tem_result[tem_al.first].second.end());
                    mrgPop.insert(tem_al.second.second.begin(),tem_al.second.second.end());
                    vector<chromosome> refPF = GnrParetoFront(mrgPop); //生成参考帕累托解集refPF；
                    tem_result[tem_al.first].second = refPF;
                }
            }
            stop_target++;
        }else{
            stop_target =0;
            for (auto & tem_al : result) {
                if(tem_result[tem_al.first].second.empty()){
                    tem_result[tem_al.first] = tem_al.second;
                }else{
                    set<chromosome> mrgPop;
                    mrgPop.insert(tem_result[tem_al.first].second.begin(),tem_result[tem_al.first].second.end());
                    mrgPop.insert(tem_al.second.second.begin(),tem_al.second.second.end());
                    vector<chromosome> refPF = GnrParetoFront(mrgPop); //生成参考帕累托解集refPF；
                    tem_result[tem_al.first].second = refPF;
                }
            }
            outfile.open("../result_Alg.txt", ios::app);
            for (auto& tem_al: tem_result) {
                outfile << tem_al.first <<" "<< temmur <<endl;
                outfile<<"cost ";
                for (auto & tem_ch: tem_al.second.second) {
                    outfile<< tem_ch.Cost<<" ";
                }
                outfile<<endl;
                outfile<<"makespan ";
                for (auto & tem_ch: tem_al.second.second) {
                    outfile<< tem_ch.MakeSpan<<" ";
                }
                outfile<<endl;
            }
            outfile.close();
            tem_result.clear();
        }
        //clculate HV, IGD, CM
        map<string, double> HVs = CalucateHV(result);
        map<string, double> IGDs = CalucateIGD(result);
        map<string, map<string, double>> CMs = CalucateCM(result);
        //HV, IGD, CM are written into the file
        outfile.open("../result.txt", ios::app);
        if (!outfile) {
            cout << "Open the file failure...\n";
            exit(0);
        }
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(3);
        if (flg2) { //输出标题行
            outfile << "case ";
            for (auto HV: HVs) outfile << " HV_" + HV.first;
            for (auto IGD: IGDs) outfile << " IGD_" + IGD.first;
            for (auto frt: CMs) {
                for (auto scn: frt.second)  outfile <<" CM_"+ frt.first + "," + scn.first;;
            }
            outfile << endl;
            flg2 = false;
        }
        outfile << Model + "_" + NumOfTask + "_" + RscAvlRatio + " ";
        for (auto HV: HVs) outfile << " " << HV.second;
        for (auto IGD: IGDs) outfile << " " << IGD.second;
        for (auto frt: CMs) {
            for (auto scn: frt.second)  outfile << " "<< scn.second;
        }
        outfile << endl;
        outfile.close();
        DeleteFirstLineInFile("../fileList.txt");//delete the first line in the file
    } while (1);
}
