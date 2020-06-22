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



importlib.reload(save_and_load)
# time counter
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()
############################################################################################################################################################
# save_and_load.load(folder_name, process_name)
ggh_event_list,ggh_mass_list, ggh_higgs_list, ggh_weight_list,ggh_image_list, ggh_recluster_images = save_and_load.load_numpy("./numpy_file_3","ggH")
vbf_event_list,vbf_mass_list, vbf_higgs_list, vbf_weight_list, vbf_image_list, vbf_recluster_images = save_and_load.load_numpy("./numpy_file_3","VBF")
vh_event_list,vh_mass_list, vh_higgs_list, vh_weight_list, vh_image_list, vh_recluster_images = save_and_load.load_numpy("./numpy_file_3","VH")
tth_event_list,tth_mass_list, tth_higgs_list, tth_weight_list, tth_image_list, tth_recluster_images = save_and_load.load_numpy("./numpy_file_3","ttH")

############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime Cost : {:.4f} min\033[0;m".format(totaltime/60.))




Train_data_file_path = './BDT_Model/Data_train.csv'
Test_data_file_path = './BDT_Model/Data_test.csv'
# Val_data_file_path = './BDT_Model/Data_val.csv'
# # read the data and store data in DataFrame
Data_train = pd.read_csv(Train_data_file_path)
Data_test = pd.read_csv(Test_data_file_path)
# Data_val = pd.read_csv(Val_data_file_path)


ggH_train = Data_train[Data_train["isGGH"]==1]
VBF_train = Data_train[Data_train["isVBF"]==1]
VH_train = Data_train[Data_train["isVH"]==1]
ttH_train = Data_train[Data_train["isttH"]==1]

ggH_test = Data_test[Data_test["isGGH"]==1]
VBF_test = Data_test[Data_test["isVBF"]==1]
VH_test = Data_test[Data_test["isVH"]==1]
ttH_test = Data_test[Data_test["isttH"]==1]

# ggH_val = Data_val[Data_val["isGGH"]==1]
# VBF_val = Data_val[Data_val["isVBF"]==1]
# VH_val = Data_val[Data_val["isVH"]==1]
# ttH_val = Data_val[Data_val["isttH"]==1]


ggh_image_train, ggh_jet_train = ggh_image_list[ggH_train["eventindex"]], ggh_recluster_images[ggH_train["eventindex"]]
ggh_image_test, ggh_jet_test = ggh_image_list[ggH_test["eventindex"]], ggh_recluster_images[ggH_test["eventindex"]]
# ggh_image_val, ggh_jet_val = ggh_image_list[ggH_val["eventindex"]], ggh_recluster_images[ggH_val["eventindex"]]

vbf_image_train, vbf_jet_train = vbf_image_list[VBF_train["eventindex"]], vbf_recluster_images[VBF_train["eventindex"]]
vbf_image_test, vbf_jet_test = vbf_image_list[VBF_test["eventindex"]], vbf_recluster_images[VBF_test["eventindex"]]
# vbf_image_val, vbf_jet_val = vbf_image_list[VBF_val["eventindex"]], vbf_recluster_images[VBF_val["eventindex"]]

vh_image_train, vh_jet_train = vh_image_list[VH_train["eventindex"]], vh_recluster_images[VH_train["eventindex"]]
vh_image_test, vh_jet_test = vh_image_list[VH_test["eventindex"]], vh_recluster_images[VH_test["eventindex"]]
# vh_image_val, vh_jet_val = vh_image_list[VH_val["eventindex"]], vh_recluster_images[VH_val["eventindex"]]

tth_image_train, tth_jet_train = tth_image_list[ttH_train["eventindex"]], tth_recluster_images[ttH_train["eventindex"]]
tth_image_test, tth_jet_test = tth_image_list[ttH_test["eventindex"]], tth_recluster_images[ttH_test["eventindex"]]
# tth_image_val, tth_jet_val = tth_image_list[ttH_val["eventindex"]], tth_recluster_images[ttH_val["eventindex"]]



print("ggH Training: Event Image {}, Jet Image {} ".format(ggh_image_train.shape,ggh_jet_train.shape))
print("ggH Test: Event Image {}, Jet Image {} ".format(ggh_image_test.shape,ggh_jet_test.shape))
# print("ggH Val.: Event Image {}, Jet Image {} ".format(ggh_image_val.shape,ggh_jet_val.shape))
print("VBF Training: Event Image {}, Jet Image {} ".format(vbf_image_train.shape,vbf_jet_train.shape))
print("VBF Test: Event Image {}, Jet Image {} ".format(vbf_image_test.shape,vbf_jet_test.shape))
# print("VBF Val.: Event Image {}, Jet Image {} ".format(vbf_image_val.shape,vbf_jet_val.shape))
print("VH Training: Event Image {}, Jet Image {} ".format(vh_image_train.shape,vh_jet_train.shape))
print("VH Test: Event Image {}, Jet Image {} ".format(vh_image_test.shape,vh_jet_test.shape))
# print("VH Val.: Event Image {}, Jet Image {} ".format(vh_image_val.shape,vh_jet_val.shape))
print("ttH Training: Event Image {}, Jet Image {} ".format(tth_image_train.shape,tth_jet_train.shape))
print("ttH Test: Event Image {}, Jet Image {} ".format(tth_image_test.shape,tth_jet_test.shape))
# print("ttH Val.: Event Image {}, Jet Image {} ".format(tth_image_val.shape,tth_jet_val.shape))

x_train_eve = np.concatenate((ggh_image_train, vbf_image_train))
x_train_eve = np.concatenate((x_train_eve, vh_image_train))
x_train_eve = np.concatenate((x_train_eve, tth_image_train))
x_test_eve = np.concatenate((ggh_image_test, vbf_image_test))
x_test_eve = np.concatenate((x_test_eve, vh_image_test))
x_test_eve = np.concatenate((x_test_eve, tth_image_test))
# x_val_eve = np.concatenate((ggh_image_val, vbf_image_val))
# x_val_eve = np.concatenate((x_val_eve, vh_image_val))
# x_val_eve = np.concatenate((x_val_eve, tth_image_val))

x_train_jet = np.concatenate((ggh_jet_train, vbf_jet_train))
x_train_jet = np.concatenate((x_train_jet, vh_jet_train))
x_train_jet = np.concatenate((x_train_jet, tth_jet_train))
x_test_jet = np.concatenate((ggh_jet_test, vbf_jet_test))
x_test_jet = np.concatenate((x_test_jet, vh_jet_test))
x_test_jet = np.concatenate((x_test_jet, tth_jet_test))
# x_val_jet = np.concatenate((ggh_jet_val, vbf_jet_val))
# x_val_jet = np.concatenate((x_val_jet, vh_jet_val))
# x_val_jet = np.concatenate((x_val_jet, tth_jet_val))

# # sample_weight = np.concatenate((ggh_weight_list[HJ_train["eventindex"]], vbf_weight_list[VBF_train["eventindex"]]))
# # sample_weight = np.concatenate((sample_weight,vh_weight_list[VH_train["eventindex"]]))

y_train = np.concatenate((np.full(len(ggh_image_train), 0), np.full(len(vbf_image_train), 1)))
y_train = np.concatenate((y_train, np.full(len(vh_image_train), 2)))
y_train = np.concatenate((y_train, np.full(len(tth_image_train), 3)))
y_test = np.concatenate((np.full(len(ggh_image_test), 0), np.full(len(vbf_image_test), 1)))
y_test = np.concatenate((y_test, np.full(len(vh_image_test), 2)))
y_test = np.concatenate((y_test, np.full(len(tth_image_test), 3)))
# y_val = np.concatenate((np.full(len(ggh_image_val), 0), np.full(len(vbf_image_val), 1)))
# y_val = np.concatenate((y_val, np.full(len(vh_image_val), 2)))
# y_val = np.concatenate((y_val, np.full(len(tth_image_val), 3)))

# print("Event: Training {} Test {} Val {}".format(x_train_eve.shape,x_test_eve.shape,x_val_eve.shape))
# print("Jet: Training {} Test {} Val {}".format(x_train_jet.shape,x_test_jet.shape,x_val_jet.shape))
# print("Target: Training {} Test {} Val {}".format(y_train.shape,y_test.shape,y_val.shape))
# input_shape = x_train_eve[0].shape

k = 0 
print(len(x_train_eve),len(x_train_jet))
for i in range(len(x_train_eve)):
    if i%25000 == 0 :
        k += 1
        
        os.mkdir("./2CNN_Model/EventTrain_"+str(k))
        os.mkdir("./2CNN_Model/JetTrain_"+str(k))
        
        event_filepath = "./2CNN_Model/EventTrain_"+str(k)
        jet_filepath = "./2CNN_Model/JetTrain_"+str(k)
        print("k= ", k)
        print("event_filepath= ",event_filepath)
        print("jet_filepath= ", jet_filepath)
        
    np.savez_compressed(event_filepath+"/x_train_eve_"+str(i+1)+".npz", x_train_eve[i])
    np.savez_compressed(jet_filepath+"/x_train_jet_"+str(i+1)+".npz", x_train_jet[i])

k = 0  
dataframe = pd.DataFrame()
path_event_image, path_jet_image = [],[]
for i in range(len(y_train)):
    if i%25000 == 0 :
        k += 1
        event_filepath = "EventTrain_"+str(k)
        jet_filepath = "JetTrain_"+str(k)
        
    path_event_image.append(event_filepath+"/x_train_eve_"+str(i+1)+".npz")
    path_jet_image.append(jet_filepath+"/x_train_jet_"+str(i+1)+".npz")
    
dataframe["EventImageTrain"] = path_event_image
dataframe["JetImageTrain"] = path_jet_image
dataframe["Y"] = y_train
dataframe.to_csv("./2CNN_Model/Train_dict.csv",index = 0)
    

    
k = 0 
print(len(x_test_eve),len(x_test_jet))
for i in range(len(x_test_eve)):
    if i%25000 == 0 :
        k += 1
        
        os.mkdir("./2CNN_Model/EventTest_"+str(k))
        os.mkdir("./2CNN_Model/JetTest_"+str(k))
        
        event_filepath = "./2CNN_Model/EventTest_"+str(k)
        jet_filepath = "./2CNN_Model/JetTest_"+str(k)
        print("k= ", k)
        print("event_filepath= ",event_filepath)
        print("jet_filepath= ", jet_filepath)
        
    np.savez_compressed(event_filepath+"/x_test_eve_"+str(i+1)+".npz", x_test_eve[i])
    np.savez_compressed(jet_filepath+"/x_test_jet_"+str(i+1)+".npz", x_test_jet[i])

    
k = 0     
dataframe = pd.DataFrame()
path_event_image, path_jet_image = [],[]
for i in range(len(y_test)):
    if i%25000 == 0 :
        k += 1
        event_filepath = "EventTest_"+str(k)
        jet_filepath = "JetTest_"+str(k)
        print("k= ", k)
        print("event_filepath= ",event_filepath)
        print("jet_filepath= ", jet_filepath)
        
    path_event_image.append(event_filepath+"/x_test_eve_"+str(i+1)+".npz")
    path_jet_image.append(jet_filepath+"/x_test_jet_"+str(i+1)+".npz")
    
dataframe["EventImageTest"] = path_event_image
dataframe["JetImageTest"] = path_jet_image
dataframe["Y"] = y_test
dataframe.to_csv("./2CNN_Model/Test_dict.csv",index = 0)