#!/bin/bash

date

this_path="/u4/alan/higgsproduction_1" 
#store_path="/u2/alan-2/higgsproduction"

echo "Start Running"


i=1
while [ $i != 11 ]
do
   echo i=$i

   date +"%Y %b %m"
   date +"%r"
   echo "GGH"
   ./bin/mg5_aMC $this_path/run_txt/ggh.txt > $this_path/outinfo/ggh/GGH_out_"$i"
   #rm -rf  $this_path/ggh/Events/run*/*.hepmc.gz
   #rsync -arz  $this_path/ggh/Events $store_path/ggh


   date +"%r"
   echo "VBFH"
   ./bin/mg5_aMC $this_path/run_txt/vbf.txt > $this_path/outinfo/vbf/VBFH_out_"$i"
   #rm -rf  $this_path/vbf/Events/run*/*.hepmc.gz
   #rsync -arz  $this_path/vbf/Events $store_path/vbf
   
   date +"%r"
   echo "VH"
   ./bin/mg5_aMC $this_path/run_txt/vh.txt > $this_path/outinfo/vh/VH_out_"$i"
   #rm -rf  $this_path/vh/Events/run*/*.hepmc.gz
   #rsync -arz  $this_path/vh/Events $store_path/vh


   date +"%r"
   echo "ttH"
   ./bin/mg5_aMC $this_path/run_txt/tth.txt > $this_path/outinfo/tth/ttH_out_"$i"
   #rm -rf  $this_path/tth/Events/run*/*.hepmc.gz
   #rsync -arz  $this_path/tth/Events $store_path/tth

   date +"%Y %b %m"
   date +"%r"
   i=$(($i+1))

done

echo "Finish"

date
