#!/bin/bash
for ((ib=0;ib<=40;ib++))
do
  for ((ip=0;ip<=40;ip++))
  do
    date >> finish.deck
    mv  Datab$[ib*1000+0]p$[ip*6+10] Data
    ./bin/gztrace
    mv Data Datab$[ib*1000+0]p$[ip*6+10]
    echo Datab$[ib*1000+0]p$[ip*6+10] is finished! >> finish.deck
  done
done
#mkdir txt
#python3 para_txt2d.py
