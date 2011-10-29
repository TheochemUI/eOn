#!/bin/sh
akmc_rate=`awk '$1 ~ /[0-9]/ {rate+= $(NF-1)} END {print rate}' akmc/states/0/processtable`
pr_rate=`awk '$1 ~ /[0-9]/ {time+=$4} END {print (NR-1)/time}' pr/states/0/processtable`
printf "akmc rate: %.3e\n" $akmc_rate
printf "pr   rate: %.3e\n" $pr_rate
