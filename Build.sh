#!/bin/bash

g++ `${ROOTSYS}/bin/root-config --cflags` -Wall -c -fPIC FittedV.cxx -o FittedV.o
g++  `${ROOTSYS}/bin/root-config --cflags` -Wall -c -fPIC VFitter.cxx -o VFitter.o
g++ `${ROOTSYS}/bin/root-config --cflags` -Wall -c -fPIC HoughTransformer.cxx -o HoughTransformer.o 
g++ `${ROOTSYS}/bin/root-config --cflags` -Wall -c -fPIC FitOrganiser.cxx -o FitOrganiser.o 
g++ -shared `${ROOTSYS}/bin/root-config --libs` HoughTransformer.o FittedV.o VFitter.o FitOrganiser.o -o libVFitter.so
