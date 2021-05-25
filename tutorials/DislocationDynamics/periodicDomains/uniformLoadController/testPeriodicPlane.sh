#!/bin/bash
counter=1
while [ $counter -le 50 ]
do
clear;
./microstructureGenerator
echo $?
((counter++))
done
echo All done
