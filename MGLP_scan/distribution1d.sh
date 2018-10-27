#!/bin/bash
alias sed=gsed
for ((ir=0;ir<=100;ir++))
do
    mkdir Datar$[ir*1+0]
    cd Datar$[ir*1+0]
    cp ../input.deck ./
    sed -i  s/ratio_out=100/ratio_out=$[ir*1+0]/g ./input.deck
    cd ../
    echo distribution in Datab$[ir*1+0] is ok!
done
