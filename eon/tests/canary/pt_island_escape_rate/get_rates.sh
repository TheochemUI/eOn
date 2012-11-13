#!/bin/sh
akmc_rate=`awk '$1 ~ /[0-9]/ {rate+= $(NF-1)} END {print rate}' akmc/states/0/processtable`
pr_rate=`awk '$1 ~ /[0-9]/ {time+=$4} END {print (NR-1)/time}' pr/states/0/processtable`
pr_hd_rate=`awk '$1 ~ /[0-9]/ {time+=$4} END {print (NR-1)/time}' pr_hd/states/0/processtable`
pr_hd_nose_rate=`awk '$1 ~ /[0-9]/ {time+=$4} END {print (NR-1)/time}' pr_hd_nose/states/0/processtable`
pr_hd_andersen_tweak=`awk '$1 ~ /[0-9]/ {time+=$4} END {print (NR-1)/time}' pr_hd_andersen_tweak/states/0/processtable`
printf "akmc                 rate: %.3e\n" $akmc_rate
printf "pr                   rate: %.3e\n" $pr_rate
printf "pr_hd                rate: %.3e\n" $pr_hd_rate
printf "pr_hd_nose           rate: %.3e\n" $pr_hd_nose_rate
printf "pr_hd_andersen_tweak rate: %.3e\n" $pr_hd_andersen_tweak
