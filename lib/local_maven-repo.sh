#!/bin/bash
###############################################
# This script installs the required java libraries not available via maven repositories from the 
# lib folder in the local maven repository
# 
# Must be run/updated when dependencies change.
###############################################

JSBML_VERSION=2498
SYSBIO_VERSION=1390
PAXTOOLS_VERSION=4.3.0
SBGN_VERSION=0.2
KEGGTRANSLATOR_VERSION=2.6


# lib directory
DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# JSBML
mvn install:install-file -DgroupId=org.sbml -DartifactId=JSBML-incl-libs -Dversion=${JSBML_VERSION} -Dfile=$DIR/JSBML-${JSBML_VERSION}-incl-libs.jar -Dpackaging=jar -DgeneratePom=true -DlocalRepositoryPath=$DIR -DcreateChecksum=true

# SysBio
mvn install:install-file -DgroupId=de.zbit -DartifactId=SysBio -Dversion=${SYSBIO_VERSION} -Dfile=$DIR/SysBio-${SYSBIO_VERSION}.jar -Dpackaging=jar -DgeneratePom=true -DlocalRepositoryPath=$DIR -DcreateChecksum=true

# SBGN
mvn install:install-file -Dfile=${DIR}/org.sbgn.jar -DgroupId=org -DartifactId=sbgn -Dversion=${SBGN_VERSION} -Dpackaging=jar -DlocalRepositoryPath=${DIR}/

# PaxTools
mvn install:install-file -Dfile=${DIR}/paxtools-${PAXTOOLS_VERSION}.jar -DgroupId=org.biopax.paxtools -DartifactId=paxtools-archetype -Dversion=${PAXTOOLS_VERSION} -Dpackaging=jar -DlocalRepositoryPath=${DIR}/

# KEGGtranslator
mvn install:install-file -Dfile=${DIR}/KEGGtranslator-${KEGGTRANSLATOR_VERSION}.jar -DgroupId=de.zbit -DartifactId=KEGGtranslator -Dversion=${KEGGTRANSLATOR_VERSION} -Dpackaging=jar -DlocalRepositoryPath=${DIR}/
