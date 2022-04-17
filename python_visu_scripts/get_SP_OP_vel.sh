#!/bin/bash

while read name;
do
 awk '{if(NR==9 || NR==29 || NR==46 || NR==57) print $2,$3,$4,$5,$6,$7,sqrt($5*$5+$6*$6+$7*$7)}' "../"$name"/point_values.txt" > "../"$name"/point_values_SP_OP.txt"
 awk '{if(NR==1) v=$7; if (NR==2) printf "%-10s %1.3f\n",  "SP vel = ", (v+$7)/2; if(NR==3) v=$7; if (NR==4) printf "%-10s %1.3f\n", "OP vel = ", (v+$7)/2}' "../"$name"/point_values_SP_OP.txt" > "../"$name"/average_vel_SP_OP.txt"
 cat "../"$name"/average_vel_SP_OP.txt"
done < list_of_model_names.txt
