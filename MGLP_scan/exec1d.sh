#!/bin/bash
for ((ir=0;ir<=100;ir++))
do
    date >> finish.deck
    mv  Datar$[ir*1+0] Data
    ./bin/gztrace
    mv Data Datar$[ir*1+0]
    echo Datar$[ir*1+0] is finished! >> finish.deck
done
