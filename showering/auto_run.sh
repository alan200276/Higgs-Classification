#!/bin/bash

date


this_path="/Storage/alan/MC_Samples"
store_path="/u2/alan-2/higgsproduction"


echo "Start Running"

# i=3
# while [ $i != 13 ]
i=4
while [ $i != 7 ]
do
   echo i=$i
   
    if [ $i -lt 10 ]; then
        a=0$i
        echo a=$a
        
    elif [ $i -gt 9 ]; then
        a=$i
        echo a=$a
        
    fi  
    date +"%Y %b %m"
#     date +"%r"
    
#     echo "GGF"
#     echo "Ghost"
#     ./Ghost_Associated_B_hadron.exe ggh_boost.cmnd ./ggh_test_2/ggh_1 $this_path/higgsproduction_1/ggh/Events/run_"$a"/events.lhe ./ggh_truth_2/ggh_1 > ./ggh_test_2/ggh_1_out_"$a".log
    
#     ./Ghost_Associated_B_hadron.exe ggh_boost.cmnd ./ggh_test_2/ggh_2 $this_path/higgsproduction_2/ggh/Events/run_"$a"/events.lhe ./ggh_truth_2/ggh_2 > ./ggh_test_2/ggh_2_out_"$a".log
    
#     ./Ghost_Associated_B_hadron.exe ggh_boost.cmnd ./ggh_test_2/ggh_3 $this_path/higgsproduction_3/ggh/Events/run_"$a"/events.lhe ./ggh_truth_2/ggh_3 > ./ggh_test_2/ggh_3_out_"$a".log

    date +"%r"
    
    echo "VBF"
    echo "Ghost"
    ./Ghost_Associated_B_hadron.exe vbf_boost.cmnd ./vbf_test_2/vbf_1 $this_path/higgsproduction_1/vbf/Events/run_"$a"/events.lhe ./vbf_truth_2/vbf_1 > ./vbf_test_2/vbf_1_out_"$a".log
    
    ./Ghost_Associated_B_hadron.exe vbf_boost.cmnd ./vbf_test_2/vbf_2 $this_path/higgsproduction_2/vbf/Events/run_"$a"/events.lhe ./vbf_truth_2/vbf_2 > ./vbf_test_2/vbf_2_out_"$a".log
    
    ./Ghost_Associated_B_hadron.exe vbf_boost.cmnd ./vbf_test_2/vbf_3 $this_path/higgsproduction_3/vbf/Events/run_"$a"/events.lhe ./vbf_truth_2/vbf_3 > ./vbf_test_2/vbf_3_out_"$a".log
    
    
#     date +"%r"
#     echo "VH"
#     echo "Ghost"
#     ./Ghost_Associated_B_hadron.exe vh_boost.cmnd ./vh_test_2/vh_1 $this_path/higgsproduction_1/vh/Events/run_"$a"/events.lhe ./vh_truth_2/vh_1 > ./vh_test_2/vh_1_out_"$a".log
    
#     ./Ghost_Associated_B_hadron.exe vh_boost.cmnd ./vh_test_2/vh_2 $this_path/higgsproduction_2/vh/Events/run_"$a"/events.lhe ./vh_truth_2/vh_2  > ./vh_test_2/vh_2_out_"$a".log
    
#     ./Ghost_Associated_B_hadron.exe vh_boost.cmnd ./vh_test_2/vh_3 $this_path/higgsproduction_3/vh/Events/run_"$a"/events.lhe ./vh_truth_2/vh_3  > ./vh_test_2/vh_3_out_"$a".log


#     date +"%r"
#     echo "ttH"
#     echo "Ghost"
#     ./Ghost_Associated_B_hadron.exe tth_boost.cmnd ./tth_test_2/tth_1 $this_path/higgsproduction_1/tth/Events/run_"$a"/events.lhe ./tth_truth_2/tth_1 > ./tth_test_2/tth_1_out_"$a".log
    
#     ./Ghost_Associated_B_hadron.exe tth_boost.cmnd ./tth_test_2/tth_2 $this_path/higgsproduction_2/tth/Events/run_"$a"/events.lhe ./tth_truth_2/tth_2 > ./tth_test_2/tth_2_out_"$a".log
    
#     ./Ghost_Associated_B_hadron.exe tth_boost.cmnd ./tth_test_2/tth_3 $this_path/higgsproduction_3/tth/Events/run_"$a"/events.lhe ./tth_truth_2/tth_3 > ./tth_test_2/tth_3_out_"$a".log



#     date +"%Y %b %m"
#     date +"%r"
#     echo "GGF"
#     echo "Frank"
#     ./myexample.exe ggh_boost.cmnd ./ggh_test_2/normal_2/ggh $this_path/ggh/Events/run_"$a"/events.lhe ./ggh_truth_2/normal_2/ggh > ./ggh_test_2/ggh_out_"$a"
    
#     date +"%r"
#     echo "Ghost"
#     ./myexample_ghost.exe ggh_boost.cmnd ./ggh_test_2/ghost_2/ggh_ghost $this_path/ggh/Events/run_"$a"/events.lhe ./ggh_truth_2/ghost_2/ggh_ghost > ./ggh_test_2/ggh_out_ghost_"$a"

#     date +"%Y %b %m"
#     date +"%r"
#     echo "VBF"
#     echo "Frank"
#     ./myexample.exe vbf_boost.cmnd ./vbf_test_2/normal_2/vbf $this_path/vbf/Events/run_"$a"/events.lhe ./vbf_truth_2/normal_2/vbf > ./vbf_test_2/vbf_out_"$a"
    
#     date +"%r"
#     echo "Ghost"
#     ./myexample_ghost.exe vbf_boost.cmnd ./vbf_test_2/ghost_2/vbf_ghost $this_path/vbf/Events/run_"$a"/events.lhe ./vbf_truth_2/ghost_2/vbf_ghost > ./vbf_test_2/vbf_out_ghost_"$a"


#     date +"%Y %b %m"
#     date +"%r"
#     echo "VH"
#     echo "Frank"
#     ./myexample.exe vh_boost.cmnd ./vh_test_2/normal_2/vh $this_path/vh/Events/run_"$a"/events.lhe ./vh_truth_2/normal_2/vh > ./vh_test_2/vh_out_"$a"
    
#     date +"%r"
#     echo "Ghost"
#     ./myexample_ghost.exe vh_boost.cmnd ./vh_test_2/ghost_2/vh_ghost $this_path/vh/Events/run_"$a"/events.lhe ./vh_truth_2/ghost_2/vh_ghost > ./vh_test_2/vh_out_ghost_"$a"


#     date +"%Y %b %m"
#     date +"%r"
#     echo "ttH"
#     echo "Frank"
#     ./myexample.exe tth_boost.cmnd ./tth_test_2/normal_2/tth $this_path/tth/Events/run_"$a"/events.lhe ./tth_truth_2/normal_2/tth > ./tth_test_2/tth_out_"$a"
    
#     date +"%r"
#     echo "Ghost"
#     ./myexample_ghost.exe tth_boost.cmnd ./tth_test_2/ghost_2/tth_ghost $this_path/tth/Events/run_"$a"/events.lhe ./tth_truth_2/ghost_2/tth_ghost > ./tth_test_2/tth_out_ghost_"$a"

   date +"%Y %b %m"
   date +"%r"
   i=$(($i+1))

done

echo "Finish"

date
