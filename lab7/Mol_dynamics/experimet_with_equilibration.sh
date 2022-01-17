#!/bin/bash

echo "Solid Argon with equilibration phase"
make Solid

echo "Liquid Argon with equilibration phase"
make Liquid

echo "Gas Argon with equilibration phase"
make Gas

make clean
