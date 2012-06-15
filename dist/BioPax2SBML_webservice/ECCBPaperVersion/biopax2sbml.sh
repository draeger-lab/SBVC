#!/bin/bash

# concatenate all, but first argument (output file)
args=""
i=1
for arg in "$@"; do
  if [ $i -gt 2 ]; then
    args=$args" "$arg
  fi
  let i=i+1
done

# remove commas
args=$(echo $args | tr ',' ' ')
infile=$2
cp $infile $infile.xml

# run BioPax2SBML
GALAXY_DIR=`pwd`"/../../.."
java -jar $GALAXY_DIR"/tools/ra_tools/biopax2sbml/biopax2sbml_0.1.jar" $args --biopax-input-file=$infile.owl> /dev/null 2>&1

# copy output file to Galaxy
rm $infile.xml
outfile=${3##*=}
mv $outfile $1
exit 0
