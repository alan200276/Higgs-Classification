#!/usr/bin/env python
# encoding: utf-8

# Import local libraries
import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import importlib

import csv_decoder
import save_and_load
import substructure



importlib.reload(csv_decoder)
# time counter
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()
############################################################################################################################################################

# Jet and event image setting
width = 40
height = 40

# Reading in custom showered files;
# This read produces event_list (collection of raw vectors) and event images

# This will be used to test saving mechanisms.


## Loading events passed through ghost-association algorithm
csv_test_path = "/Storage/alan/FrankShower"

# print('Loading ggH events')
# hj_event_list,hj_mass_list,hj_image_list, hj_higgs_list, hj_weight_list, num_hj_files = \
#     csv_decoder.load_events(path= csv_test_path + "/ggh_test_2/",\
#                 contains=".csv",pt_cut=1, width=width, height=height)

print('Loading VBF events')
vbf_event_list,vbf_mass_list, vbf_image_list, vbf_higgs_list, vbf_weight_list, num_vbf_files = \
    csv_decoder.load_events(path= csv_test_path + "/vbf_test_2/",\
                contains=".csv",pt_cut=1, width=width, height=height)

# print('Loading VH events')
# vh_event_list,vh_mass_list, vh_image_list, vh_higgs_list, vh_weight_list, num_vh_files = \
#     csv_decoder.load_events(path= csv_test_path + "/vh_test_2/",\
#                 contains=".csv",pt_cut=1, width=width, height=height)

# print('Loading ttH events')
# tth_event_list,tth_mass_list, tth_image_list, tth_higgs_list, tth_weight_list, num_tth_files = \
#     csv_decoder.load_events(path= csv_test_path + "/tth_test_2/",\
#                 contains=".csv",pt_cut=1, width=width, height=height)


# np.save("./numpy_file_4/vbf_event_list",vbf_event_list)
# np.save("./numpy_file_4/vh_event_list",vh_event_list)
# np.save("./numpy_file_4/tth_event_list",tth_event_list)

############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime consumption : {:.4f} min\033[0;m".format(totaltime/60.))




importlib.reload(csv_decoder)
# time counter
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()
###########################################################################################################################################################

print("Jet Clustering")

# hj_clustered = csv_decoder.cluster_event(hj_event_list)
vbf_clustered = csv_decoder.cluster_event(vbf_event_list)
# vh_clustered = csv_decoder.cluster_event(vh_event_list)
# tth_clustered = csv_decoder.cluster_event(tth_event_list)

# hj_non_higgs_jets, hj_non_higgs_fat_jets, hj_higgs_jet, hj_higgs_fat_jet, hj_non_higgs_jet_list = csv_decoder.recluster_event(hj_clustered, hj_higgs_list)
vbf_non_higgs_jets, vbf_non_higgs_fat_jets, vbf_higgs_jet, vbf_higgs_fat_jet, vbf_non_higgs_jet_list = csv_decoder.recluster_event(vbf_clustered, vbf_higgs_list)
# vh_non_higgs_jets, vh_non_higgs_fat_jets, vh_higgs_jet, vh_higgs_fat_jet, vh_non_higgs_jet_list = csv_decoder.recluster_event(vh_clustered, vh_higgs_list)
# tth_non_higgs_jets, tth_non_higgs_fat_jets, tth_higgs_jet, ttH_higgs_fat_jet, tth_non_higgs_jet_list = csv_decoder.recluster_event(tth_clustered, tth_higgs_list)

# print("# of HJ",len(hj_clustered))
print("# of VBF",len(vbf_clustered))
# print("# of VH",len(vh_clustered))
# print("# of ttH",len(tth_clustered))




############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime consumption : {:.4f} min\033[0;m".format(totaltime/60.))




importlib.reload(csv_decoder)
# time counter
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()
############################################################################################################################################################


# Produce jet images, the zero-center and normalize
print('Producing jet images')
# hj_recluster_images = csv_decoder.return_jet_image_list(hj_event_list,hj_non_higgs_jets,0.8, width=width, height=height)
vbf_recluster_images = csv_decoder.return_jet_image_list(vbf_event_list,vbf_non_higgs_jets,0.8, width=width, height=height)
# vh_recluster_images = csv_decoder.return_jet_image_list(vh_event_list,vh_non_higgs_jets,0.8, width=width, height=height)
# tth_recluster_images = csv_decoder.return_jet_image_list(tth_event_list,tth_non_higgs_jets,0.8, width=width, height=height)

# Zero centering and normalizing
# # hj_image_list, vbf_image_list, vh_image_list, tth_image_list = csv_decoder.zero_center_and_normalize((hj_image_list,vbf_image_list,vh_image_list,tth_image_list))
# # hj_recluster_images, vbf_recluster_images, vh_recluster_images, tth_recluster_images  = csv_decoder.zero_center_and_normalize(( hj_recluster_images,vbf_recluster_images,vh_recluster_images,tth_recluster_images))

vbf_image_list = csv_decoder.zero_center_and_normalize(vbf_image_list)
vbf_recluster_images = csv_decoder.zero_center_and_normalize(vbf_recluster_images)

############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime consumption : {:.4f} min\033[0;m".format(totaltime/60.))



importlib.reload(save_and_load)
# time counter
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()
############################################################################################################################################################
print("Saving Data")
# ## save_and_load.save(folder_name, process_name, event_list, mass_list, higgs, weight, image_list):
# save_and_load.save_numpy("./numpy_file_4","ggH",hj_event_list,hj_mass_list, hj_higgs_list, hj_weight_list,hj_image_list,hj_recluster_images)
save_and_load.save_numpy("./numpy_file_4","VBF",vbf_event_list,vbf_mass_list, vbf_higgs_list, vbf_weight_list, vbf_image_list,vbf_recluster_images)
# save_and_load.save_numpy("./numpy_file_4","VH",vh_event_list,vh_mass_list, vh_higgs_list, vh_weight_list, vh_image_list,vh_recluster_images)
# save_and_load.save_numpy("./numpy_file_4","ttH",tth_event_list,tth_mass_list, tth_higgs_list, tth_weight_list, tth_image_list,tth_recluster_images)

############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime Cost : {:.4f} min\033[0;m".format(totaltime/60.))




print("Finish")


