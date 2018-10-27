#!/bin/bash
for ((ib=0;ib<=40;ib++))
do
  for ((ip=0;ip<=40;ip++))
  do
    date >> finish.deck
    mv  Datab$[ib*1000+0]p$[ip*6+10] Data
    mkdir C-Rad
    ./a.out
    mv C-Rad Data
    mv Data Datab$[ib*1000+0]p$[ip*6+10]
    echo post Datab$[ib*1000+0]p$[ip*6+10] ! >> finish.deck
  done
done
#mkdir txt
#python3 para_txt2d.py
