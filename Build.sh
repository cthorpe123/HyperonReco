#!/bin/bash

mkdir -p lib

INC=-I${HYP_TOP}/Core

g++ `${ROOTSYS}/bin/root-config --cflags` -Wall -c -fPIC FittedV.cxx -o lib/FittedV.o
g++  `${ROOTSYS}/bin/root-config --cflags` -Wall -c -fPIC VFitter.cxx -o lib/VFitter.o
g++ `${ROOTSYS}/bin/root-config --cflags` -Wall -c -fPIC HoughTransformer.cxx -o lib/HoughTransformer.o 
g++ `${ROOTSYS}/bin/root-config --cflags` -Wall -c -fPIC ${INC} FitOrganiser.cxx -o lib/FitOrganiser.o 
g++ -shared `${ROOTSYS}/bin/root-config --libs` lib/HoughTransformer.o lib/FittedV.o lib/VFitter.o lib/FitOrganiser.o -o lib/libVFitter.so
