#!/bin/bash
#
# $1: BioPax input file
# $2: SBML output file
# $3: (real) input file name

# run BioPax2SBML
# echo $1
# echo $2
# echo $3
#exit
GALAXY_DIR=`pwd`"/../../../.."
applicationPath=$GALAXY_DIR"/tools/ra_tools/BioPax2SBML";# BioPax2SBML Verzeichnis

# Make a path for the outfile
dataPath=$2;
dataPath="${dataPath%.*}_files"
if [ ! -d $dataPath ]; then
  mkdir $dataPath
fi

outputFile=$dataPath
outputFile="${dataPath}/${3}.sbml.xml"

cd "$applicationPath"
java -Xms128m -Xmx2G -jar biopax2sbml_1.0.jar --input=$1 --output=$outputFile --format=SBML_QUAL --log-level=SEVERE --gui=false > /dev/null 2>&1



cd $dataPath
#Get a single file as result
returnFile="biopax2sbml.zip"
if [ `ls -1 | wc -l` -gt 0 ]; then
  zip -9 -o -r -q $returnFile * > /dev/null 2>&1
else
  returnFile=`ls -1`
fi

#Copy this file for the user
mv $returnFile $2


