#!/bin/bash
#
# $1: BioPax input dir RELATIVE TO CURRENT DIRECTORY


for filename in `pwd`/$1/*.owl.gz
do
  echo Converting $filename...
  ./biopax2sbml.sh $filename yes "${filename%.*}".sbml.zip interactive >> infos.txt 2>> errors.txt
done;
