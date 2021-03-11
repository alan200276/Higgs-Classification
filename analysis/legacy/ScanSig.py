# Import useful libraries
# Import useful libraries
import importlib
import numpy as np
import substructure
import matplotlib.pyplot as plt
import time
import pandas as pd
import matplotlib
matplotlib.rc('font', size=15)
import os

# Import local libraries
import csv_decoder
import save_and_load
import function as fn


def PT(px, py):
    return np.sqrt(px*px + py*py)

def M(e, px, py, pz):
    return np.sqrt(e*e - (px*px + py*py + pz*pz))

my_data_file_path = './BDT_Model_for_test/Data_train_w_pre.csv'
# read the data and store data in DataFrame
Data_test = pd.read_csv(my_data_file_path)

list = ["higgs_pt","higgs_eta","higgs_m",
           "non_higgs_leading_pt","non_higgs_leading_eta","non_higgs_leading_m",
           "non_higgs_second_pt","non_higgs_second_eta","non_higgs_second_m",
           "girth","CIJS","SIJS","DeltaEta",
            "weight","isHJ","isVBF","isVH","eventindex","prediction"]

Data_test["PTH"],Data_test["evtweight"] = Data_test["higgs_pt"],Data_test["weight"]



def Z(s,b):
    z = np.sqrt(2*((s+b)*np.log(1+s/b)-s))
#     z = s/np.sqrt(b)
    return z


Luminosity = 300 #fb-1 run3
significance = [] 
PTCUT = np.linspace(400,1250,86)

zh_ratio = 0.661816317958867
wh_ratio = 0.338183682041133
ggF_Factor, VBF_Factor, VH_Factor, ttH_Factor = 229066./229000., 374654./229000., 278766./229000., 259114./229000.



for ptcut in PTCUT:
    N = 100
    cut = np.zeros([N+1])
    NggH, NVBF, NVH, NttH, Nother = np.zeros([N+1]), np.zeros([N+1]), np.zeros([N+1]), np.zeros([N+1]), np.zeros([N+1])
    cut[0] = 0
    ggF = Data_test[(Data_test["isGGH"]==1) & (Data_test["prediction"] > 0) & (Data_test["PTH"] > ptcut)]    #prediction  #BDT
    VBF = Data_test[(Data_test["isVBF"]==1) & (Data_test["prediction"] > 0) & (Data_test["PTH"] > ptcut)]
    VH = Data_test[(Data_test["isVH"]==1) & (Data_test["prediction"] > 0) & (Data_test["PTH"] > ptcut)]
    ttH = Data_test[(Data_test["isttH"]==1) & (Data_test["prediction"] > 0) & (Data_test["PTH"] > ptcut)]
    Other = Data_test[(Data_test["isGGH"]!=1) & (Data_test["prediction"] > 0) & (Data_test["PTH"] > ptcut)]
    
    ggh_BDTcut = fn.Weight(ggF,"ggF",PTmin=ptcut)
    vbf_BDTcut = fn.Weight(VBF,"VBF",PTmin=ptcut)
    vh_BDTcut = fn.Weight(VH,"VH",PTmin=ptcut)
    tth_BDTcut = fn.Weight(ttH,"ttH",PTmin=ptcut)

    NggH[0] = (np.array(ggh_BDTcut[2])*0.5824*ggF_Factor)[0]*Luminosity*1000 
    NVBF[0] = (np.array(vbf_BDTcut[2])*0.5824*VBF_Factor)[0]*Luminosity*1000 
    NVH[0] =  (np.array(vh_BDTcut[2])*0.5824*(wh_ratio*0.6741+zh_ratio*0.6991)*VH_Factor)[0]*Luminosity*1000 
    NttH[0] = (np.array(tth_BDTcut[2])*0.5824*ttH_Factor)[0]*Luminosity*1000 
    Nother[0] = len(Other)

    for i in range(N):
        cut[i+1] = cut[i] + 1./N
        ggF = Data_test[(Data_test["isGGH"]==1) & (Data_test["prediction"] > cut[i+1]) & (Data_test["PTH"] > ptcut)]
        VBF = Data_test[(Data_test["isVBF"]==1) & (Data_test["prediction"] > cut[i+1]) & (Data_test["PTH"] > ptcut)]
        VH = Data_test[(Data_test["isVH"]==1) & (Data_test["prediction"] > cut[i+1]) & (Data_test["PTH"] > ptcut)]
        ttH = Data_test[(Data_test["isttH"]==1) & (Data_test["prediction"] > cut[i+1]) & (Data_test["PTH"] > ptcut)]
        Other = Data_test[(Data_test["isGGH"]!=1) & (Data_test["prediction"] > cut[i+1]) & (Data_test["PTH"] > ptcut)]
        
        
        ggh_BDTcut = fn.Weight(ggF,"ggF",PTmin=ptcut)
        vbf_BDTcut = fn.Weight(VBF,"VBF",PTmin=ptcut)
        vh_BDTcut = fn.Weight(VH,"VH",PTmin=ptcut)
        tth_BDTcut = fn.Weight(ttH,"ttH",PTmin=ptcut)

        NggH[i+1] = (np.array(ggh_BDTcut[2])*0.5824*ggF_Factor)[0]*Luminosity*1000 
        NVBF[i+1] = (np.array(vbf_BDTcut[2])*0.5824*VBF_Factor)[0]*Luminosity*1000 
        NVH[i+1] = (np.array(vh_BDTcut[2])*0.5824*(wh_ratio*0.6741+zh_ratio*0.6991)*VH_Factor)[0]*Luminosity*1000 
        NttH[i+1] = (np.array(tth_BDTcut[2])*0.5824*ttH_Factor)[0]*Luminosity*1000 

        if NggH[i+1] <= 0:
            N = i
            break
        elif NVBF[i+1]+NVH[i+1]+NttH[i+1] ==0:
            N = i
            break
    
    sig = np.array(Z(NggH[:N],(NVBF[:N]+NVH[:N]+NttH[:N])))
    
    if len(sig) == 0:
        break
        
    cut_opt = cut[np.where(sig == max(sig))][0]
    print(np.where(sig == max(sig)))
    print(NggH[np.where(sig == max(sig))],NVBF[np.where(sig == max(sig))])
    print(cut_opt)

    max_sig = max(sig)
    significance.append(max_sig)
    
np.save("./BDT_Model_for_test/significancescan",np.array(significance))
# np.save("./2CNN_Model_0/significancescan",np.array(significance))