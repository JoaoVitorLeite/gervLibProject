#!/bin/bash

dataset_train="../datasets/cities_norm.csv"
separator_train=","
distance_function="EUCLIDEAN"
pivot_type=("RANDOM" "GNAT" "CONVEX" "KMEDOIDS" "MAXSEPARATED" "MAXVARIANCE" "SELECTION" "PCA" "SSS" "HFI" "IS" "WDR")
sample_size=(1.0 1.0 1.0 0.2 1.0 0.3 0.3 0.2 0.2 0.2 0.5 0.5)
seed=($(shuf -i 0-500000 -n 13))
path_save_results="../results/"

cp pivot_expt.cpp main.cpp
#cd ../gervLib_Project/
#make clean
#/usr/lib/qt5/bin/qmake -makefile gervLib.pro
Qt/6.5.0/gcc_64/bin/qmake -makefile gervLib.pro 
make
#cd ../gervLib/bin
cd bin/

for((i=0; i<12; i++));
do
    ./gervLib -DATASET_TRAIN ${dataset_train} -DATASET_TRAIN_SEPARATOR ${separator_train} -DISTANCE_FUNCTION ${distance_function} -PIVOT_TYPE ${pivot_type[$i]} -SAMPLE_SIZE_PIVOT ${sample_size[$i]} -SEED ${seed[$i]} -REP 1 -PATH_SAVE_RESULTS ${path_save_results}
done


