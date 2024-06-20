#!/bin/bash

Name1=obs
Name2=PolyChain
for((i=1;i<=1;i++))
do  
    g++ $Name1.cpp $Name2.cpp -o $Name1.o -O2 -std=c++11
    OMP_NUM_THREADS=1 ./$Name1.o 
done