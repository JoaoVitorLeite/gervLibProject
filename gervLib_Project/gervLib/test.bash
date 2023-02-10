#!/bin/bash

num_per_leaf=340

sed -i "s/const int LC = -1/const int LC = "${num_per_leaf}"/g" main.cpp


