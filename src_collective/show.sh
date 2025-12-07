#!/bin/sh


for i in `seq 1 9` ;do
  cat ${1}/res_0_${i} | grep "Packet" | tail -1 | cut -d " " -f 5
done
cat ${1}/res_1_0 | grep "Packet" | tail -1 | cut -d " " -f 5


