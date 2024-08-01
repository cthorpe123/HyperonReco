#!/bin/bash

g++ `${ROOTSYS}/bin/root-config --cflags` -c -fPIC FittedV.cxx -o FittedV.o
g++  `${ROOTSYS}/bin/root-config --cflags` -c -fPIC VFitter.cxx -o VFitter.o
g++ `${ROOTSYS}/bin/root-config --cflags` -c -fPIC HoughTransformer.cxx -o HoughTransformer.o 
g++ -shared `${ROOTSYS}/bin/root-config --libs` HoughTransformer.o FittedV.o VFitter.o -o libVFitter.so
