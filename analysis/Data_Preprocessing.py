# Import local libraries
from sklearn.utils import shuffle
import csv_decoder  #self-defined
import save_and_load #self-defined
import substructure #self-defined
import importlib
import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import os

from matplotlib.ticker import MaxNLocator
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, BoundaryNorm
from matplotlib.collections import LineCollection
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
import copy


importlib.reload(save_and_load)
# time counter
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()
############################################################################################################################################################

file_path = "../numpy_file_train"

# save_and_load.load(folder_name, process_name)
ggh_event_list,ggh_mass_list, ggh_higgs_list, ggh_weight_list = save_and_load.load_numpy(file_path,"ggH")
vbf_event_list,vbf_mass_list, vbf_higgs_list, vbf_weight_list = save_and_load.load_numpy(file_path,"VBF")
vh_event_list,vh_mass_list, vh_higgs_list, vh_weight_list = save_and_load.load_numpy(file_path,"VH")
tth_event_list,tth_mass_list, tth_higgs_list, tth_weight_list = save_and_load.load_numpy(file_path,"ttH")

############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime consumption : {:.4f} min\033[0;m".format(totaltime/60.))


importlib.reload(csv_decoder)
# time counter
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()
############################################################################################################################################################


print("Jet Clustering")

ggh_clustered = csv_decoder.cluster_event(ggh_event_list)
vbf_clustered = csv_decoder.cluster_event(vbf_event_list)
vh_clustered = csv_decoder.cluster_event(vh_event_list)
tth_clustered = csv_decoder.cluster_event(tth_event_list)

ggh_non_higgs_jets, ggh_non_higgs_fat_jets, ggh_higgs_jet, ggh_higgs_fat_jet, ggh_non_higgs_jet_list = csv_decoder.recluster_event(ggh_clustered, ggh_higgs_list)
vbf_non_higgs_jets, vbf_non_higgs_fat_jets, vbf_higgs_jet, vbf_higgs_fat_jet, vbf_non_higgs_jet_list = csv_decoder.recluster_event(vbf_clustered, vbf_higgs_list)
vh_non_higgs_jets, vh_non_higgs_fat_jets, vh_higgs_jet, vh_higgs_fat_jet, vh_non_higgs_jet_list = csv_decoder.recluster_event(vh_clustered, vh_higgs_list)
tth_non_higgs_jets, tth_non_higgs_fat_jets, tth_higgs_jet, tth_higgs_fat_jet, tth_non_higgs_jet_list = csv_decoder.recluster_event(tth_clustered, tth_higgs_list)

print("# of HJ",len(ggh_clustered))
print("# of VBF",len(vbf_clustered))
print("# of VH",len(vh_clustered))
print("# of ttH",len(tth_clustered))

############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime consumption : {:.4f} min\033[0;m".format(totaltime/60.))

# importlib.reload(csv_decoder)
# # time counter
# print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
# ticks_1 = time.time()
# ############################################################################################################################################################

# ggh_event_image_rotated = csv_decoder.return_image_list(csv_decoder.Rotate_Event_List(ggh_event_list,ggh_higgs_list), width=40, height=40)
# vbf_event_image_rotated = csv_decoder.return_image_list(csv_decoder.Rotate_Event_List(vbf_event_list,vbf_higgs_list), width=40, height=40)
# vh_event_image_rotated = csv_decoder.return_image_list(csv_decoder.Rotate_Event_List(vh_event_list,vh_higgs_list), width=40, height=40)
# tth_event_image_rotated = csv_decoder.return_image_list(csv_decoder.Rotate_Event_List(tth_event_list,tth_higgs_list), width=40, height=40)

# ggh_event_images_rotated_zcno, \
# vbf_event_images_rotated_zcno, \
# vh_event_images_rotated_zcno, \
# tth_event_images_rotated_zcno = csv_decoder.zero_center_and_normalize((copy.deepcopy(ggh_event_image_rotated),\
#                                                                        copy.deepcopy(vbf_event_image_rotated),\
#                                                                        copy.deepcopy(vh_event_image_rotated), \
#                                                                        copy.deepcopy(tth_event_image_rotated)))

# np.savez_compressed("./Test_Images/ggh_event_images_rotated_zcno_1n1c1c", ggh_event_images_rotated_zcno)
# np.savez_compressed("./Test_Images/vbf_event_images_rotated_zcno_1n1c1c", vbf_event_images_rotated_zcno)
# np.savez_compressed("./Test_Images/vh_event_images_rotated_zcno_1n1c1c", vh_event_images_rotated_zcno)
# np.savez_compressed("./Test_Images/tth_event_images_rotated_zcno_1n1c1c", tth_event_images_rotated_zcno)


# ###########################################################################################################################################################
# ticks_2 = time.time()
# totaltime =  ticks_2 - ticks_1
# print("\033[3;33mTime consumption : {:.4f} min\033[0;m".format(totaltime/60.))


importlib.reload(csv_decoder)

# Jet and event image setting
width = 40
height = 40

# time counter
print(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
ticks_1 = time.time()
############################################################################################################################################################

# Produce jet images, the zero-center and normalize
print('Producing jet images')
ggh_jet_images = csv_decoder.return_jet_image_list(ggh_event_list,ggh_non_higgs_jets,0.8, width=width, height=height)
vbf_jet_images = csv_decoder.return_jet_image_list(vbf_event_list,vbf_non_higgs_jets,0.8, width=width, height=height)
vh_jet_images = csv_decoder.return_jet_image_list(vh_event_list,vh_non_higgs_jets,0.8, width=width, height=height)
tth_jet_images = csv_decoder.return_jet_image_list(tth_event_list,tth_non_higgs_jets,0.8, width=width, height=height)

# Zero centering and normalizing
ggh_jet_images_zcno, \
vbf_jet_images_zcno, \
vh_jet_images_zcno, \
tth_jet_images_zcno  = csv_decoder.zero_center_and_normalize(( copy.deepcopy(ggh_jet_images), \
                                                               copy.deepcopy(vbf_jet_images), \
                                                               copy.deepcopy(vh_jet_images), \
                                                               copy.deepcopy(tth_jet_images)))

np.savez_compressed("./Train_Images/ggh_jet_images_zcno_1n1c1c", ggh_jet_images_zcno)
np.savez_compressed("./Train_Images/vbf_jet_images_zcno_1n1c1c", vbf_jet_images_zcno)
np.savez_compressed("./Train_Images/vh_jet_images_zcno_1n1c1c", vh_jet_images_zcno)
np.savez_compressed("./Train_Images/tth_jet_images_zcno_1n1c1c", tth_jet_images_zcno)

############################################################################################################################################################
ticks_2 = time.time()
totaltime =  ticks_2 - ticks_1
print("\033[3;33mTime consumption : {:.4f} min\033[0;m".format(totaltime/60.))
