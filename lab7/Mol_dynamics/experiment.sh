#!/bin/bash

echo "Solid Argon already equilibrised"
make Solid_termalized

echo "Liquid Argon already equilibrised"
make Liquid_termalized

echo "Gas Argon already equilibrised"
make Gas_termalized

make clean