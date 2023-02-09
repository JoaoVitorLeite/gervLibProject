#!/bin/bash

index="OMNIKDTREE"
dataset_train="../datasets/train_nasa.csv"
separator_train=","
dataset_test="../datasets/test_nasa.csv"
separator_test=","
distance_function="EUCLIDEAN"
pivot_type=("RANDOM" "GNAT" "CONVEX" "KMEDOIDS" "MAXSEPARETED" "MAXVARIANCE" "SELECTION" "PCA" "SSS" "FFT" "HFI" "IS" "WDR")
sample_size=(1.0 1.0 1.0 0.2 1.0 0.3 0.3 0.01 0.03 1.0 0.1 0.5 0.5)
num_pivots=2
seed=($(shuf -i 0-500000 -n 13))
k_max=10
num_query=1000
num_per_leaf=360 #only omni,kdtree
path_save_results="../results/"

mkdir -p omni/omni_files/
mkdir -p kdtree/kdtree_files/
mkdir -p results/

#make clean
cp query.cpp main.cpp
#cp query_str.cpp main.cpp
cd ../gervLib_Project/
make clean
#/usr/lib/qt5/bin/qmake -makefile gervLib.pro
Qt/6.1.0/gcc_64/bin/qmake -makefile gervLib.pro 
#cd bin/
make
cd ../gervLib/bin

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# for((i=0; i<13; i++));
# do
#     nohup ./gervLib -INDEX ${index} -DATASET_TRAIN ${dataset_train} -DATASET_TRAIN_SEPARATOR ${separator_train} -DATASET_TEST ${dataset_test} -DATASET_TEST_SEPARATOR ${separator_test} -DISTANCE_FUNCTION ${distance_function} -PIVOT_TYPE ${pivot_type[$i]} -SAMPLE_SIZE_PIVOT ${sample_size[$i]} -NUM_PIVOTS ${num_pivots} -SEED ${seed[$i]} -K_MAX ${k_max} -PATH_SAVE_RESULTS ${path_save_results} -NUM_PER_LEAF ${num_per_leaf} &
# done
# 


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for((i=0; i<9; i++));
do
    nohup ./gervLib -INDEX ${index} -DATASET_TRAIN ${dataset_train} -DATASET_TRAIN_SEPARATOR ${separator_train} -DATASET_TEST ${dataset_test} -DATASET_TEST_SEPARATOR ${separator_test} -DISTANCE_FUNCTION ${distance_function} -PIVOT_TYPE ${pivot_type[$i]} -SAMPLE_SIZE_PIVOT ${sample_size[$i]} -NUM_PIVOTS ${num_pivots} -SEED ${seed[$i]} -K_MAX ${k_max} -PATH_SAVE_RESULTS ${path_save_results} -NUM_PER_LEAF ${num_per_leaf} &
done

wait

for((i=9; i<13; i++));
do
    nohup ./gervLib -INDEX ${index} -DATASET_TRAIN ${dataset_train} -DATASET_TRAIN_SEPARATOR ${separator_train} -DATASET_TEST ${dataset_test} -DATASET_TEST_SEPARATOR ${separator_test} -DISTANCE_FUNCTION ${distance_function} -PIVOT_TYPE ${pivot_type[$i]} -SAMPLE_SIZE_PIVOT ${sample_size[$i]} -NUM_PIVOTS ${num_pivots} -SEED ${seed[$i]} -K_MAX ${k_max} -PATH_SAVE_RESULTS ${path_save_results} -NUM_PER_LEAF ${num_per_leaf} &
done
