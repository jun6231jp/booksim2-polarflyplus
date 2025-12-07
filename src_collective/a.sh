#!/bin/sh


#$1: file

line=`cat $1 | grep "##step" | wc -l`

for i in `seq 1 ${line}`
do
 node=`cat ${1} | grep "##step" | cut -d ":" -f 1 | head -$i | tail -1 | cut -d " " -f 2` 
 for j in `seq 2 1000`
 do
  step=`cat ${1} | grep "##step" | cut -d ":" -f 2 | head -$i | tail -1 | cut -d " " -f $j`
  txstep=`cat ${1} | grep "##txstep" | cut -d ":" -f 2 | head -$i | tail -1 | cut -d " " -f $j`
  rxstep=`cat ${1} | grep "##rxstep" | cut -d ":" -f 2 | head -$i | tail -1 | cut -d " " -f $j`
  if [ ${step} -eq 0 ]; then
	  echo "${node} step $((${j}-1)) : ${step} ${rxstep} ${txstep}"
    break
  fi
 done 
done
