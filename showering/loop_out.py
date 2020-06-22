import numpy as np
import sys
import os


N_truth = []
N_selected = []
total = 0

path = sys.argv[1]


N_truth = []
N_selected = []
total = 0

for j, filename in enumerate(os.listdir(path)):
    if os.path.isfile(os.path.join(path,filename)) and (".log" in filename):
#         print('Currently reading: '+ str(filename))
        
        fp = open(path + str(filename), 'r')

        lines = fp.readlines()
        length = len(lines)
        N_truth.append(np.array(lines[length - 1].strip().split()[0],dtype=float))
        N_selected.append(np.array(lines[length - 2].strip().split()[0],dtype=float))
        total += 100000.

        fp.close()
    
    
N_truth = np.array(N_truth)
N_selected= np.array(N_selected)
total_truth = np.sum(N_truth)
std_truth = np.std(N_truth)/np.sqrt(10)
total_selected = np.sum(N_selected)
std_selected = np.std(N_selected)/np.sqrt(10)
    
    
print("Total Events before merging: {:.0f}".format(total))
print("Number of events after merging: {:.0f} +-{:.3f}%".format(total_truth,std_truth/total_truth*100))
print("Number of events after selected: {:.0f} +-{:.3f}%".format(total_selected,std_selected/total_selected*100))
print("Merging efficiency: {:.4f} +-{:.4f} %".format(np.average(N_truth/100000.),np.std(N_truth/100000.)/np.sqrt(10)/(np.average(N_truth/100000.))*100))
print("Selection efficiency: {:.4f} +-{:.4f} %".format(np.average(N_selected/N_truth), np.std(N_selected/N_truth)/np.sqrt(10)/(total_selected/total_truth)*100))

print("Merging X Selection efficiency: {:.4f} +-{:.4f} %".format(total_selected/total, np.std(N_selected/100000.)/np.sqrt(10)/(total_selected/total)*100))


