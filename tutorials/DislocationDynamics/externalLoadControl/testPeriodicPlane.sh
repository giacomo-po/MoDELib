#!/bin/bash
counter=1
while [ $counter -le 50 ]
do
clear;./microstructureGenerator
((counter++))
done
echo All done
