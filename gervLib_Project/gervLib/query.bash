#!/bin/bash

index="OMNIKDTREE"
dataset_train="../datasets/train_cities_norm.csv"
separator_train=","
dataset_test="../datasets/test_cities_norm.csv"
separator_test=","
distance_function="EUCLIDEAN"
pivot_type=("BPP" "CONVEX" "GNAT" "SSS" "MAXSEPARATED" "MAXVARIANCE" "PCA" "IS" "HFI" "WDR" "SELECTION" "KMEDOIDS" "RANDOM")
sample_size=(1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0)
num_pivots=2
seed=($(shuf -i 0-500000 -n 13))
k_max=100
rep=10
page_size=2220
num_bins=256
num_per_leaf=55
path_save_results="../results/tcc_expt/consultas/cities/omni/"

mkdir -p omni/omni_files/
mkdir -p kdtree/kdtree_files/
mkdir -p pm_tree/pm_files/
mkdir -p spb_tree/spb_files/
mkdir -p results/

cp query.cpp main.cpp

sed -i "s/const int LC = -1/const int LC = "${num_per_leaf}"/g" main.cpp
sed -i '20s/.*/static const size_t PAGE_SIZE = '${page_size}';/' config_spb.h
sed -i '22s/.*/static const size_t MIN_BTREE_LEAF_NUM = '${num_per_leaf}';/' config_spb.h
#sed -i "s/static const size_t PAGE_SIZE = -1/static const size_t PAGE_SIZE = "${page_size}"/g" config_spb.h
#sed -i "s/static const size_t MIN_BTREE_LEAF_NUM = -1/static const size_t MIN_BTREE_LEAF_NUM = "${num_per_leaf}"/g" config_spb.h

/usr/lib/qt5/bin/qmake -makefile gervLib.pro
#Qt/6.5.0/gcc_64/bin/qmake -makefile gervLib.pro 
make
cd bin/

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------

for((i=0; i<9; i++));
do
    nohup ./gervLib -INDEX ${index} -DATASET_TRAIN ${dataset_train} -DATASET_TRAIN_SEPARATOR ${separator_train} -DATASET_TEST ${dataset_test} -DATASET_TEST_SEPARATOR ${separator_test} -DISTANCE_FUNCTION ${distance_function} -PIVOT_TYPE ${pivot_type[$i]} -SAMPLE_SIZE_PIVOT ${sample_size[$i]} -NUM_PIVOTS ${num_pivots} -SEED ${seed[$i]} -K_MAX ${k_max} -PATH_SAVE_RESULTS ${path_save_results} -NUM_PER_LEAF ${num_per_leaf} -REP ${rep} -PAGE_SIZE ${page_size} -NUM_BINS ${num_bins} &
done

wait

for((i=9; i<13; i++));
do
    nohup ./gervLib -INDEX ${index} -DATASET_TRAIN ${dataset_train} -DATASET_TRAIN_SEPARATOR ${separator_train} -DATASET_TEST ${dataset_test} -DATASET_TEST_SEPARATOR ${separator_test} -DISTANCE_FUNCTION ${distance_function} -PIVOT_TYPE ${pivot_type[$i]} -SAMPLE_SIZE_PIVOT ${sample_size[$i]} -NUM_PIVOTS ${num_pivots} -SEED ${seed[$i]} -K_MAX ${k_max} -PATH_SAVE_RESULTS ${path_save_results} -NUM_PER_LEAF ${num_per_leaf} -REP ${rep} -PAGE_SIZE ${page_size} -NUM_BINS ${num_bins} &
done



