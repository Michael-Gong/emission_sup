#!/bin/bash
for ((ib=0;ib<=40;ib++))
do
  for ((ip=0;ip<=40;ip++))
  do
    mkdir Datab$[ib*1000+0]p$[ip*6+10]
    cd Datab$[ib*1000+0]p$[ip*6+10]
    cp ../input.deck ./
    sed -i  s/a0_out=0/a0_out=$[ib*1000+0]/g ./input.deck
    sed -i  s/py0=10/py0=$[ip*6+10]/g ./input.deck
   # sed -i  s/ratio_out=40/ratio_out=$[ir*1+0]/g ./input.deck
   # sed -i '.original'  s/a0_out=1/a0_out=$[ib*1+0]/g ./input.deck
    cd ../
    echo distribution in Datab$[ib*1000+0]p$[ip*6+10] is ok!
  done
done
