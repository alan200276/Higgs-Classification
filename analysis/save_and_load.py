import numpy as np
import pickle
import csv_decoder
import os

def save(folder_name, process_name, event_list, mass_list, higgs, weight,image_list,recluster_images):
	if not os.path.exists(folder_name):
		os.makedirs(folder_name)
	if not os.path.exists(folder_name+'/'+ process_name + '_events/' ):
		os.makedirs(folder_name+'/'+ process_name + '_events/' )
		
	# Saving each event into an individual npy file (safer)
	for i in range(len(event_list)):
		np.save(folder_name+'/'+ process_name + '_events/' + str(i) + '.npy', event_list[i])

	# Saving all arrays with np.save (faster)

	np.save(folder_name+'/'+process_name+'_mass_list.npy', mass_list)
	np.save(folder_name+'/'+process_name+'_higgs.npy', higgs)
	np.save(folder_name+'/'+process_name+'_weight.npy', weight)
	np.save(folder_name+'/'+process_name+'_image_list.npy', image_list)
	np.save(folder_name+'/'+process_name+'_recluster_images.npy', recluster_images) 
	
def load(folder_name, process_name):

	# Loading all arrays with np.save (faster)
	new_mass_list = np.load(folder_name+'/'+process_name+'_mass_list.npy', allow_pickle=True)
	new_weight = np.load(folder_name+'/'+process_name+'_weight.npy', allow_pickle=True)
	new_higgs = np.load(folder_name+'/'+process_name+'_higgs.npy', allow_pickle=True)
	new_image_list = np.load(folder_name+'/'+process_name+'_image_list.npy', allow_pickle=True)
	new_recluster_images = np.load(folder_name+'/'+process_name+'_recluster_images.npy', allow_pickle=True)
	
	new_event_list = []
	for i in range(len(new_mass_list)):
		new_event_list.append(np.load(folder_name+'/'+ process_name + '_events/' + str(i) + '.npy', allow_pickle=True))


	return new_event_list, new_mass_list,\
		new_higgs, new_weight, new_image_list, \
		new_recluster_images


def load_numpy(folder_name, process_name):

	# Loading all arrays with np.save (faster)
	new_mass_list = np.load(folder_name+'/'+process_name+'_mass_list.npy', allow_pickle=True)
	new_weight = np.load(folder_name+'/'+process_name+'_weight.npy', allow_pickle=True)
	new_higgs = np.load(folder_name+'/'+process_name+'_higgs.npy', allow_pickle=True)
	new_image_list = np.load(folder_name+'/'+process_name+'_image_list.npy', allow_pickle=True)
	new_recluster_images = np.load(folder_name+'/'+process_name+'_recluster_images.npy', allow_pickle=True)
	new_event_list = np.load(folder_name+'/'+process_name+'_event_list.npy', allow_pickle=True)

	return new_event_list, new_mass_list,\
		new_higgs, new_weight, new_image_list, \
		new_recluster_images

def save_numpy(folder_name, process_name, event_list, mass_list, higgs, weight,image_list,recluster_images):
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    if not os.path.exists(folder_name+'/'+ process_name + '_events/' ):
        os.makedirs(folder_name+'/'+ process_name + '_events/' )
        
    # Saving all arrays with np.save (faster)
    np.save(folder_name+'/'+process_name+'_mass_list.npy', mass_list)
    np.save(folder_name+'/'+process_name+'_higgs.npy', higgs)
    np.save(folder_name+'/'+process_name+'_weight.npy', weight)
    np.save(folder_name+'/'+process_name+'_image_list.npy', image_list)
    np.save(folder_name+'/'+process_name+'_recluster_images.npy', recluster_images)
    np.save(folder_name+'/'+process_name+'_event_list.npy', event_list)