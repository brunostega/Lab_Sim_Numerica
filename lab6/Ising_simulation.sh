#!/bin/bash                                                                     

echo "Ising with no magnetic Field h=0"

cp input_h_off.dat input.dat
echo "Ising with Metropolis sampling"
make Metropolis
cp outputMetropolis/*.* outputMetropolis/H_field_off

echo "Ising with Gibbs sampling"
make Gibbs
cp outputGibbs/*.* outputGibbs/H_field_off


echo "Ising with magnetic Field h=0.2"

cp input_h_on.dat input.dat
echo "Ising with Metropolis sampling"
make Metropolis
cp outputMetropolis/*.* outputMetropolis/H_field_on

echo "Ising with Gibbs sampling"
make Gibbs
cp outputGibbs/*.* outputGibbs/H_field_on




