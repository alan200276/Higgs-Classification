#!/bin/bash

setup_PYTHIA() {
    export PYTHIA8LOCATION=/root/pythia8244
    export LD_LIBRARY_PATH=/root/pythia8244/lib/:$LD_LIBRARY_PATH
    export PYTHIA8=/root/pythia8244
}


setup_fastjet() {
    export FASTJETLOCATION=/root/fastjet_3_2_1
    export LD_LIBRARY_PATH=/root/fastjet_3_2_1/lib/:$LD_LIBRARY_PATH
}


setup_PYTHIA
setup_fastjet



#export LD_LIBRARY_PATH=/home/alan/root/rootinstall_2/lib/:$LD_LIBRARY_PATH
