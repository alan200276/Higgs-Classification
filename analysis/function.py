import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def PT(px, py):
    return np.sqrt(px*px + py*py)

def M(e, px, py, pz):
    return np.sqrt(e*e - (px*px + py*py + pz*pz))


def Weight(process, name, PTmin = 400 , PTmax = 2000, index = 1 ):
    def reweight_bin(LHCHXSWG_process,LHCHXSWG_EW):
        LHCHXSWG_process_bin = []
        for i in range(1,len(LHCHXSWG_process)):
            LHCHXSWG_process_bin.append(LHCHXSWG_process[i-1]*(1-LHCHXSWG_EW[i]/100.)-LHCHXSWG_process[i]*(1-LHCHXSWG_EW[i]/100.))
        LHCHXSWG_process_bin = np.array(LHCHXSWG_process_bin)
        return LHCHXSWG_process_bin
    
    
    # 400,410,420,...,600,650,700,750,800,....,1250 GeV
    LHCHXSWG_ggF = np.array([33.30,29.34,25.95,22.97,20.39,18.08,16.01,14.27,12.77,11.39,
                    10.17,9.11,8.15,7.29,6.52,5.87,5.29,4.76,4.29,3.86,3.48,
                    2.13,1.32,0.84,0.54,0.36,0.24,0.16,0.11,0.0733,0.0504,0.0349,0.0243,0.0168])
    ggh_bin_weight = np.load("./reweighting/ggh_bin_weight_"+str(index)+".npy")
    LHCHXSWG_ggF_bin = []
    for i in range(1,len(LHCHXSWG_ggF)):
        LHCHXSWG_ggF_bin.append(LHCHXSWG_ggF[i-1]-LHCHXSWG_ggF[i])
#         LHCHXSWG_ggF_bin.append(LHCHXSWG_ggF[-1])
    LHCHXSWG_ggF_bin = np.array(LHCHXSWG_ggF_bin)
    
    
    # 400,450,500,550,600,650,700,750,800 GeV
    LHCHXSWG_VBF = np.array([14.23,8.06,4.75,2.90,1.82,1.17,0.77,0.51,0.35])
    EW_VBF = np.array([17.80,19.43,21.05,22.34,23.73,25.03,26.29,27.35,28.42])
    vbf_bin_weight = np.load("./reweighting/vbf_bin_weight_"+str(index)+".npy")
    
    LHCHXSWG_VH = np.array([11.16,6.87,4.39,2.87,1.91,1.30,0.90,0.62,0.44])
    EW_VH = np.array([19.05,20.83,22.50,24.07,25.56,26.98,28.30,29.60,30.83])
    vh_bin_weight = np.load("./reweighting/vh_bin_weight_"+str(index)+".npy")
    
    LHCHXSWG_ttH = np.array([6.89,4.24,2.66,1.76,1.11,0.72,0.47,0.32,0.22])
    EW_ttH = np.array([6.95,7.75,8.94,9.11,9.91,10.67,11.37,11.94,12.51])
    tth_bin_weight = np.load("./reweighting/tth_bin_weight_"+str(index)+".npy")
    
    if name == "ggF":
        LHCHXSWG_process_bin = LHCHXSWG_ggF_bin
        for_reweight = ggh_bin_weight
    elif name == "VBF":
        LHCHXSWG_process_bin = reweight_bin(LHCHXSWG_VBF,EW_VBF)
        for_reweight = vbf_bin_weight
    elif name == "VH":
        LHCHXSWG_process_bin = reweight_bin(LHCHXSWG_VH,EW_VH)
        for_reweight = vh_bin_weight
    elif name == "ttH": 
        LHCHXSWG_process_bin = reweight_bin(LHCHXSWG_ttH,EW_ttH)
        for_reweight = tth_bin_weight
        
    weight = []
    totalweight, uncertainty, uncertainty_perbin = [], [], []
    
    region1 = np.linspace(400,600,int((600-400)/10+1))
    region2 = np.linspace(650,1250,int((1250-650)/50+1))
    PTCUT_ggF = np.concatenate((region1,region2))
    LHCHXSWG_pt = np.linspace(400,800,int((800-400)/50+1))
    PTCUT = np.linspace(PTmin,PTmax,int((PTmax-PTmin)/10+1))
#     PTCUT = np.linspace(400,1250,int((1250-400)/10+1))
    
    if name == "ggF":
        k = 0

        for i in range(len(PTCUT)-1):
            if len(process[(process["PTH"] >= PTCUT[i]) & (process["PTH"] < PTCUT[i+1]) ]) == 0:
                length = i
                break
            else:
                bin_weights = process[(process["PTH"] >= PTCUT[i]) & (process["PTH"] < PTCUT[i+1]) ]["evtweight"]/300000
                reweight_factor = LHCHXSWG_process_bin[k]*0.001/for_reweight[k]

                length = i+1
            if PTCUT_ggF[k] < 1200:
                if PTCUT[i] == PTCUT_ggF[k+1]:
                    k += 1
            elif PTCUT_ggF[k] >= 1250:
                reweight_factor = 0.050843596565854896
                
            weight.append(sum(bin_weights*reweight_factor))
            uncertainty_perbin.append(np.sqrt(np.sum((bin_weights*reweight_factor)**2)))
                

#         for i in range(len(PTCUT[:length+1])):
#             totalweight.append(sum(weight[i:]))
# #             print(sum(weight[i:]))
#             uncertainty.append(np.sqrt(np.sum(np.array(uncertainty_perbin)[i:]**2)))
        for i in range(len(PTCUT[:length+1])):
            if PTCUT[i] < 1450:
                totalweight.append(sum(weight[i:np.where(PTCUT==1450)[0][0]]))
                uncertainty.append(np.sqrt(np.sum(np.array(uncertainty_perbin)[i:np.where(PTCUT==1450)[0][0]]**2)))
            elif PTCUT[i] >= 1450:
                totalweight.append(sum(weight[i:]))
                uncertainty.append(np.sqrt(np.sum(np.array(uncertainty_perbin)[i:]**2)))
    else:
        k = 0
        for i in range(len(PTCUT)-1):
            if len(process[(process["PTH"] >= PTCUT[i]) & (process["PTH"] < PTCUT[i+1]) ]) == 0:
                length = i
                break
            else:
                bin_weights = process[(process["PTH"] >= PTCUT[i]) & (process["PTH"] < PTCUT[i+1]) ]["evtweight"]/300000
                reweight_factor = LHCHXSWG_process_bin[k]*0.001/for_reweight[k]

                weight.append(sum(bin_weights*reweight_factor))
                uncertainty_perbin.append(np.sqrt(np.sum((bin_weights*reweight_factor)**2)))

                length = i+1
            if LHCHXSWG_pt[k] < 750:
                if PTCUT[i] == LHCHXSWG_pt[k+1]:
                    k += 1
            
        for i in range(len(PTCUT[:length+1])):
            totalweight.append(sum(weight[i:]))
            uncertainty.append(np.sqrt(np.sum(np.array(uncertainty_perbin)[i:]**2)))
        
            
    return PTCUT, weight, totalweight, uncertainty_perbin, uncertainty

def error_propagation(xs_1,xs_2,xs_3,xs_4,error_1,error_2,error_3,error_4):
    xs_1,xs_2,xs_3,xs_4 = np.array(xs_1),np.array(xs_2),np.array(xs_3),np.array(xs_4)
    totalerror, div_error_1, div_error_2, div_error_3, div_error_4 = [],[],[],[],[]
    totaleXS = []
    length = min(len(error_1),len(error_2),len(error_3),len(error_4))
    for i in range(length):
        totalerror.append(np.sqrt(np.array(error_1)[i]**2+np.array(error_2)[i]**2+np.array(error_3)[i]**2+np.array(error_4)[i]**2))
        totaleXS.append(xs_1[i]+xs_2[i]+xs_3[i]+xs_4[i])
    totalerror, totaleXS = np.array(totalerror), np.array(totaleXS)
    
    for i in range(length):
        div_error_1.append((xs_1[i]/totaleXS[i])*np.sqrt((error_1[i]/xs_1[i])**2+(totalerror[i]/totaleXS[i])**2))
        div_error_2.append((xs_2[i]/totaleXS[i])*np.sqrt((error_2[i]/xs_2[i])**2+(totalerror[i]/totaleXS[i])**2))
        div_error_3.append((xs_3[i]/totaleXS[i])*np.sqrt((error_3[i]/xs_3[i])**2+(totalerror[i]/totaleXS[i])**2))
        div_error_4.append((xs_4[i]/totaleXS[i])*np.sqrt((error_4[i]/xs_4[i])**2+(totalerror[i]/totaleXS[i])**2))
    
    return totalerror, div_error_1, div_error_2, div_error_3, div_error_4


def DrawCumulativeXection(ggh_weight, vbf_weight, vh_weight, tth_weight, 
                          ggh_factor=1, vbf_factor=1, vh_factor=1, tth_factor=1):
    fig, ax = plt.subplots(1,1, figsize=(9,6))

    plt.fill_between(ggh_weight[0][:len(ggh_weight[2])], np.array(ggh_weight[2])*ggh_factor
                     -np.array(ggh_weight[4])*ggh_factor, 
                     np.array(ggh_weight[2])*ggh_factor
                     +np.array(ggh_weight[4])*ggh_factor, color = 'green', alpha =0.5,linewidth=1)

    plt.fill_between(vbf_weight[0][:len(vbf_weight[2])], np.array(vbf_weight[2])*vbf_factor
                     -np.array(vbf_weight[4])*vbf_factor, 
                     np.array(vbf_weight[2])*vbf_factor
                     +np.array(vbf_weight[4])*vbf_factor, color = 'red', alpha =0.5,linewidth=1)

    plt.fill_between(vh_weight[0][:len(vh_weight[2])], np.array(vh_weight[2])*vh_factor
                     -np.array(vh_weight[4])*vh_factor, 
                     np.array(vh_weight[2])*vh_factor
                     +np.array(vh_weight[4])*vh_factor, color = 'blue', alpha =0.5,linewidth=1)

    plt.fill_between(tth_weight[0][:len(tth_weight[2])], np.array(tth_weight[2])*tth_factor
                     -np.array(tth_weight[4])*tth_factor, 
                     np.array(tth_weight[2])*tth_factor
                     +np.array(tth_weight[4])*tth_factor, color = 'magenta', alpha =0.5,linewidth=1)

    plt.plot(ggh_weight[0][:len(ggh_weight[2])],np.array(ggh_weight[2])*ggh_factor, color = "green",linewidth=3,label="ggF") #weighted
    plt.plot(vbf_weight[0][:len(vbf_weight[2])],np.array(vbf_weight[2])*vbf_factor, color = "red",linewidth=3,label="VBF") #weighted
    plt.plot(vh_weight[0][:len(vh_weight[2])],np.array(vh_weight[2])*vh_factor, color = "blue",linewidth=3,label="VH") #weighted
    plt.plot(tth_weight[0][:len(tth_weight[2])],np.array(tth_weight[2])*tth_factor, color = "magenta",linewidth=3,label="ttH") #weighted

    plt.yscale("log")
    plt.xlim((400,1250))
    plt.ylim((1E-6,1E-1))
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.yaxis.set_ticks_position('both')
    plt.tick_params(axis='y', which='both', labelleft='on',labelright=0,direction="in")
    plt.tick_params(axis='x', direction="in")
    plt.grid(True)
    plt.xlabel('$P^H_{t}$ [GeV]', fontsize=25,horizontalalignment='right',x=1)
    plt.ylabel("$\sum(p^H_t)$ [pb] ",fontsize=25,horizontalalignment='right',y=1)
    plt.legend(ncol=1,fontsize=20)
    plt.tight_layout()


    # plt.savefig("./Higgs_Pt/Cumulative_Xection", transparent=True)
    plt.show()
    
    
def DrawFractional(ggh_weight, vbf_weight, vh_weight, tth_weight, 
                   ggh_factor=1, vbf_factor=1, vh_factor=1, tth_factor=1):
    
    fig, ax = plt.subplots(1,1, figsize=(9,6))

    length = min(len(ggh_weight[2]),len(vbf_weight[2]),len(vh_weight[2]),len(tth_weight[2]))-1

    error = error_propagation(np.array(ggh_weight[2][:length])*ggh_factor,
                              np.array(vbf_weight[2][:length])*vbf_factor,
                              np.array(vh_weight[2][:length])*vh_factor ,
                              np.array(tth_weight[2][:length])*tth_factor,
                              np.array(ggh_weight[4][:length])*ggh_factor,
                              np.array(vbf_weight[4][:length])*vbf_factor,
                              np.array(vh_weight[4][:length])*vh_factor ,
                              np.array(tth_weight[4][:length])*tth_factor
                             )

    total_Xection = []
    vbf_rate, vh_rate, tth_rate = [], [], []


    for i in range(length): 
        total_Xection.append(np.array(ggh_weight[2])[i]*ggh_factor +\
                    np.array(vbf_weight[2])[i]*vbf_factor + \
                    np.array(vh_weight[2])[i]*vh_factor + \
                    np.array(tth_weight[2])[i]*tth_factor)

    total_Xection = np.array(total_Xection)

    plt.fill_between(ggh_weight[0][:length], np.array(ggh_weight[2])[:length]/total_Xection*ggh_factor 
                     -np.array(error[1]), 
                     np.array(ggh_weight[2])[:length]/total_Xection*ggh_factor 
                     +np.array(error[1]), color = 'green', alpha =0.5,linewidth=1)

    plt.fill_between(vbf_weight[0][:length], np.array(vbf_weight[2])[:length]/total_Xection*vbf_factor
                     -np.array(error[2]), 
                     np.array(vbf_weight[2])[:length]/total_Xection*vbf_factor
                     +np.array(error[2]), color = 'red', alpha =0.5,linewidth=1)

    plt.fill_between(vh_weight[0][:length], np.array(vh_weight[2])[:length]/total_Xection*vh_factor
                     -np.array(error[3]), 
                     np.array(vh_weight[2])[:length]/total_Xection*vh_factor
                     +np.array(error[3]), color = 'blue', alpha =0.5,linewidth=1)

    plt.fill_between(tth_weight[0][:length], np.array(tth_weight[2])[:length]/total_Xection*tth_factor
                     -np.array(error[4]), 
                     np.array(tth_weight[2])[:length]/total_Xection*tth_factor
                     +np.array(error[4]), color = 'magenta', alpha =0.5,linewidth=1)
    
    plt.plot(ggh_weight[0][:length],np.array(ggh_weight[2])[:length]/total_Xection*ggh_factor , color = "green",linewidth=3,label="ggF") #weighted
    plt.plot(vbf_weight[0][:length],np.array(vbf_weight[2])[:length]/total_Xection*vbf_factor, color = "red",linewidth=3,label="VBF") #weighted
    plt.plot(vh_weight[0][:length],np.array(vh_weight[2])[:length]/total_Xection*vh_factor, color = "blue",linewidth=3,label="VH") #weighted
    plt.plot(tth_weight[0][:length],np.array(tth_weight[2])[:length]/total_Xection*tth_factor, color = "magenta",linewidth=3,label="ttH") #weighted

    # plt.yscale("log")
    plt.xlim((400,1250))
    plt.ylim((0,1))
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.yaxis.set_ticks_position('both')
    plt.tick_params(axis='y', which='both', labelleft='on',labelright=0,direction="in")
    plt.tick_params(axis='x', direction="in")
    plt.grid(True)
    plt.xlabel('$P^H_{t}$ [GeV]', fontsize=25,horizontalalignment='right',x=1)
    plt.ylabel("Fractional Contribution",fontsize=20,horizontalalignment='right',y=1)
    plt.legend(ncol=1,fontsize=20)
    plt.tight_layout()

    # plt.savefig("./Higgs_Pt/CumulativeXection_w_BDTcut_Ratio", transparent=True)
    plt.show()