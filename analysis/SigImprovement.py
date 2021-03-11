# Import useful libraries
# Import useful libraries
import importlib
import numpy as np
import time
import pandas as pd
import os
import sys
# Import local libraries
import function as fn

# import csv_decoder
# import save_and_load
# import substructure
# import matplotlib.pyplot as plt

try :
    prediction_file = sys.argv[1]
    save_file = sys.argv[2]
    axis = sys.argv[3]  #0:"ggF" 1:"VBF" 2:"VH" 3:"ttH" 
    
#     if int(axis) == 0:
#         process = "isGGH"
#     elif int(axis) == 1:
#         process = "isVBF"
#     elif int(axis) == 2:
#         process = "isVH"
#     elif int(axis) == 3:
#         process = "isttH"
    
except:
    print("**************Check Input Again!!!!**************")



def PT(px, py):
    return np.sqrt(px*px + py*py)

def M(e, px, py, pz):
    return np.sqrt(e*e - (px*px + py*py + pz*pz))

my_data_file_path = './BDT_Model_for_test/Data_train.csv'
# read the data and store data in DataFrame
Data_test = pd.read_csv(my_data_file_path)


prediction = np.load(prediction_file)

list = ["higgs_pt","higgs_eta","higgs_m",
           "non_higgs_leading_pt","non_higgs_leading_eta","non_higgs_leading_m",
            "mjj",
           "non_higgs_second_pt","non_higgs_second_eta","non_higgs_second_m",
           "girth","CIJS","SIJS","DeltaEta",
            "weight","isHJ","isVBF","isVH","eventindex","prediction"]

Data_test["PTH"],Data_test["evtweight"] = Data_test["higgs_pt"],Data_test["weight"]
Data_test["pre_out1"],Data_test["pre_out2"] = prediction[:,0], prediction[:,1]
Data_test["pre_out3"],Data_test["pre_out4"] = prediction[:,2], prediction[:,3]

# Data_test["pre_GBDT_out1"],Data_test["pre_GBDT_out2"] = pre_GBDT[:,0], pre_GBDT[:,1]
# Data_test["pre_GBDT_out3"],Data_test["pre_GBDT_out4"] = pre_GBDT[:,2], pre_GBDT[:,3]

# Data_test["pre_2CNN_out1"],Data_test["pre_2CNN_out2"] = pre_2CNN[:,0], pre_2CNN[:,1]
# Data_test["pre_2CNN_out3"],Data_test["pre_2CNN_out4"] = pre_2CNN[:,2], pre_2CNN[:,3]



def Z(s,b):
    z = np.sqrt(2*((s+b)*np.log(1+s/b)-s))
#     z = s/np.sqrt(b)
    return z


Luminosity = 300 #fb-1 run3
# significance = [] 
# PTCUT = np.linspace(400,1200,81)

ptcut = 400


zh_ratio = 0.661816317958867
wh_ratio = 0.338183682041133
ggF_Factor, VBF_Factor, VH_Factor, ttH_Factor = 229066./229000., 374654./229000., 278766./229000., 259114./229000.


print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()
############################################################################################################################################################



# for ptcut in PTCUT:
N = 100
cut = np.zeros([N+1])
NggH, NVBF, NVH, NttH, Nother = np.zeros([N+1]), np.zeros([N+1]), np.zeros([N+1]), np.zeros([N+1]), np.zeros([N+1])
cut[0] = 0
ggF = Data_test[(Data_test["isGGH"]==1) & (Data_test["pre_out1"] > 0) & (Data_test["PTH"] > ptcut)]    #prediction  #BDT
VBF = Data_test[(Data_test["isVBF"]==1) & (Data_test["pre_out2"] > 0) & (Data_test["PTH"] > ptcut)]
VH = Data_test[(Data_test["isVH"]==1) & (Data_test["pre_out3"] > 0) & (Data_test["PTH"] > ptcut)]
ttH = Data_test[(Data_test["isttH"]==1) & (Data_test["pre_out4"] > 0) & (Data_test["PTH"] > ptcut)]
# Other = Data_test[(Data_test[process]!=1) & (Data_test["pre_out"+str(int(axis)+1)] > 0) & (Data_test["PTH"] > ptcut)] 

ggh_cut = fn.Weight(ggF,"ggF",PTmin=ptcut, index=2)
vbf_cut = fn.Weight(VBF,"VBF",PTmin=ptcut, index=2)
vh_cut = fn.Weight(VH,"VH",PTmin=ptcut, index=2)
tth_cut = fn.Weight(ttH,"ttH",PTmin=ptcut, index=2)


NggH[0] = (np.array(ggh_cut[2])*0.5824*ggF_Factor)[0]*Luminosity*1000 
NVBF[0] = (np.array(vbf_cut[2])*0.5824*VBF_Factor)[0]*Luminosity*1000 
NVH[0] =  (np.array(vh_cut[2])*0.5824*(wh_ratio*0.6741+zh_ratio*0.6991)*VH_Factor)[0]*Luminosity*1000 
NttH[0] = (np.array(tth_cut[2])*0.5824*ttH_Factor)[0]*Luminosity*1000 
# Nother[0] = len(Other)

for i in range(N):
    cut[i+1] = cut[i] + 1./N
    ggF = Data_test[(Data_test["isGGH"]==1) & (Data_test["pre_out"+str(int(axis)+1)] > cut[i+1]) & (Data_test["PTH"] > ptcut)]
    VBF = Data_test[(Data_test["isVBF"]==1) & (Data_test["pre_out"+str(int(axis)+1)] > cut[i+1]) & (Data_test["PTH"] > ptcut)]
    VH = Data_test[(Data_test["isVH"]==1) & (Data_test["pre_out"+str(int(axis)+1)] > cut[i+1]) & (Data_test["PTH"] > ptcut)]
    ttH = Data_test[(Data_test["isttH"]==1) & (Data_test["pre_out"+str(int(axis)+1)] > cut[i+1]) & (Data_test["PTH"] > ptcut)]
#     Other = Data_test[(Data_test[process]!=1) & (Data_test["pre_out"+str(axis+1)] > cut[i+1]) & (Data_test["PTH"] > ptcut)]
    
    ggh_cut = fn.Weight(ggF,"ggF",PTmin=ptcut, index=2)
    vbf_cut = fn.Weight(VBF,"VBF",PTmin=ptcut, index=2)
    vh_cut = fn.Weight(VH,"VH",PTmin=ptcut, index=2)
    tth_cut = fn.Weight(ttH,"ttH",PTmin=ptcut, index=2)
    
    NggH[i+1] = (np.array(ggh_cut[2])*0.5824*ggF_Factor)[0]*Luminosity*1000 
    NVBF[i+1] = (np.array(vbf_cut[2])*0.5824*VBF_Factor)[0]*Luminosity*1000 
    NVH[i+1] = (np.array(vh_cut[2])*0.5824*(wh_ratio*0.6741+zh_ratio*0.6991)*VH_Factor)[0]*Luminosity*1000 
    NttH[i+1] = (np.array(tth_cut[2])*0.5824*ttH_Factor)[0]*Luminosity*1000 

#     if len(ggF) == 0:
#         N = i
#         break
    if int(axis) == 0:
        if NggH[i+1] <= 0:
            N = i
            break
        elif NVBF[i+1]+NVH[i+1]+NttH[i+1] ==0:
            N = i
            break
    elif int(axis) == 1:        
        if NVBF[i+1] <= 0:
            N = i
            break
        elif NggH[i+1]+NVH[i+1]+NttH[i+1] ==0:
            N = i
            break
    elif int(axis) == 2:   
        if NVH[i+1] <= 0:
            N = i
            break
        elif NggH[i+1]+NVBF[i+1]+NttH[i+1] ==0:
            N = i
            break
    elif int(axis) == 3:    
        if NttH[i+1] <= 0:
            N = i
            break
        elif NggH[i+1]+NVBF[i+1]+NVH[i+1] ==0:
            N = i
            break

        
        
    
        
        
if int(axis) == 0:
    sig = np.array(Z(NggH[:N],(NVBF[:N]+NVH[:N]+NttH[:N]))/Z(NggH[0],(NVBF[0]+NVH[0]+NttH[0])))
    
    print(np.where(sig == max(sig)))
    index = np.where(sig == max(sig))
    cut_opt = NggH[index]/NggH[0]
    print("cut_opt",cut_opt)
    print("ggF efficiency: ", NggH[index]/NggH[0])
    print("Others rejection rate: ",  1/((NVBF[index]+NVH[index]+NttH[index])/(NVBF[0]+NVH[0]+NttH[0])))
    rejection_rate = 1/((NVBF[index]+NVH[index]+NttH[index])/(NVBF[0]+NVH[0]+NttH[0]))
    max_sig = max(sig)
    
elif int(axis) == 1:
    sig = np.array(Z(NVBF[:N],(NggH[:N]+NVH[:N]+NttH[:N]))/Z(NggH[0],(NggH[0]+NVH[0]+NttH[0])))
    
    print(np.where(sig == max(sig)))
    index = np.where(sig == max(sig))
    cut_opt = NVBF[index]/NVBF[0]
    print("cut_opt",cut_opt)
    print("ggF efficiency: ", NVBF[index]/NVBF[0])
    print("Others rejection rate: ", 1/((NggH[index]+NVH[index]+NttH[index])/(NggH[0]+NVH[0]+NttH[0])))
    rejection_rate = 1/((NggH[index]+NVH[index]+NttH[index])/(NggH[0]+NVH[0]+NttH[0]))
    max_sig = max(sig)
    
elif int(axis) == 2:
    sig = np.array(Z(NVH[:N],(NggH[:N]+NVBF[:N]+NttH[:N]))/Z(NVH[0],(NggH[0]+NVBF[0]+NttH[0])))    
    
    print(np.where(sig == max(sig)))
    index = np.where(sig == max(sig))
    cut_opt = NVH[index]/NVH[0]
    print("cut_opt",cut_opt)
    print("ggF efficiency: ", NVH[index]/NVH[0])
    print("Others rejection rate: ", 1/((NggH[index]+NVBF[index]+NttH[index])/(NggH[0]+NVBF[0]+NttH[0])))
    rejection_rate = 1/((NggH[index]+NVBF[index]+NttH[index])/(NggH[0]+NVBF[0]+NttH[0]))
    max_sig = max(sig)
    
elif int(axis) == 3:
    sig = np.array(Z(NttH[:N],(NggH[:N]+NVBF[:N]+NVH[:N]))/Z(NttH[0],(NggH[0]+NVBF[0]+NVH[0])))   
    
    print(np.where(sig == max(sig)))
    index = np.where(sig == max(sig))
    cut_opt = NttH[index]/NttH[0]
    print("cut_opt",cut_opt)
    print("ggF efficiency: ", NttH[index]/NttH[0])
    print("Others rejection rate: ", 1/((NggH[index]+NVBF[index]+NVH[index])/(NggH[0]+NVBF[0]+NVH[0])))
    rejection_rate = 1/((NggH[index]+NVBF[index]+NVH[index])/(NggH[0]+NVBF[0]+NVH[0]))
    max_sig = max(sig)
    
# sig = np.array(Z(NggH[:N],(NVBF[:N]+NVH[:N]+NttH[:N])))
# sig = np.array(Z(NVBF[:N],(NggH[:N]+NVH[:N]+NttH[:N])))
# sig = np.array(Z(NVH[:N],(NggH[:N]+NVBF[:N]+NttH[:N])))

# cut_opt = cut[np.where(sig == max(sig))][0]
# print(np.where(sig == max(sig)))
# print(NggH[np.where(sig == max(sig))],NVBF[np.where(sig == max(sig))])


# index = np.where(sig == max(sig))
# cut_opt = NggH[index]/NggH[0]
# print("cut_opt",cut_opt)
# print("ggF efficiency: ", NggH[index]/NggH[0])
# print("Others rejection rate: ", 1-(NVBF[index]+NVH[index]+NttH[index])/(NVBF[0]+NVH[0]+NttH[0]))

# max_sig = max(sig)
    
# significance = np.array([NggH[:N]/NggH[0],sig])
if int(axis) == 0:  
    np.savez_compressed("./Significance/sig_impro"+str(save_file)+"_ggF", 
                        ggF_eff = NggH[:N]/NggH[0], 
                        sig_impro = sig,
                        rejection_rate = rejection_rate
                        )
    
elif int(axis) == 1:
    np.savez_compressed("./Significance/sig_impro"+str(save_file)+"_VBF", 
                        VBF_eff = NVBF[:N]/NVBF[0], 
                        sig_impro = sig,
                        rejection_rate = rejection_rate
                       )

elif int(axis) == 2:
    np.savez_compressed("./Significance/sig_impro"+str(save_file)+"_VH", 
                        VH_eff = NVH[:N]/NVH[0], 
                        sig_impro = sig,
                        rejection_rate = rejection_rate
                       )

elif int(axis) == 3:
    np.savez_compressed("./Significance/sig_impro"+str(save_file)+"_ttH", 
                        ttH_eff = NttH[:N]/NttH[0], 
                        sig_impro = sig,
                        rejection_rate = rejection_rate
                       )
    
############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime consumption : {:.4f} min for GBDT\033[0;m".format(totaltime/60.))