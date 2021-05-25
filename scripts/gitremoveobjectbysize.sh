#!/bin/bash

while read line
  do
    set $line          # assigns words in line to positional parameters
    if [[ "$5" == "tutorials/DislocationDynamics/"*"T_"*".txt" 
       || "$5" == "tutorials/DislocationDynamics/"*"N_"*".txt" 
       || "$5" == "tutorials/DislocationDynamics/"*"evl/evl_"*".txt" 
       || "$5" == "tutorials/DislocationDynamics/"*"evl/evl_"*".bin" 
       || "$5" == "tutorials/DislocationDynamics/"*"evl/ddAux_"*".txt" 
       || "$5" == "tutorials/DislocationDynamics/"*"evl/ddAux_"*".bin" 
       || "$5" == "tutorials/"*"/mesh/mesh.neigh" 
       || "$5" == "tutorials/"*"/mesh/mesh.ele" 
       || "$5" == "tutorials/"*"DDparallel" 
       || "$5" == "tutorials/"*"DDserial" 
       || "$5" == "tutorials/"*"DDcode" 
       || "$5" == *"piGenerator" 
       || "$5" == *"straightGenerator" 
       || "$5" == *"frGenerator" 
       || "$5" == *"dipoleGenerator" 
       || "$5" == *"junctionGenerator"
       || "$5" == *"pairGenerator" 
       || "$5" == "tutorials/"*".off" 
       || "$5" == "tutorials/"*"tetgen" 
       || "$5" == "tutorials/ParticleInteraction/ChargedParticles/P_0compare.txt"
       || "$5" == "test/logNormalDist/probability.txt" ]]; then
  		echo Removing "$4 $5"
  		git filter-branch -f --tree-filter "rm -f $5" -- --all
	fi
  done < sortedFilesBySize.txt
