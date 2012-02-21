/*
 * $Id$
 * $URL$
 * ---------------------------------------------------------------------
 * This file is part of the SysBio API library.
 *
 * Copyright (C) 2009-2012 by the University of Tuebingen, Germany.
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation. A copy of the license
 * agreement is provided in the file named "LICENSE.txt" included with
 * this software distribution and also available online as
 * <http://www.gnu.org/licenses/lgpl-3.0-standalone.html>.
 * ---------------------------------------------------------------------
 */
package de.zbit.biopax;

import java.io.FileNotFoundException;
import java.util.Map.Entry;
import java.util.logging.Level;

import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.model.level3.PhysicalEntity;

import de.zbit.biopax.BioPaxPathwayHolder;
import de.zbit.biopax.BioPax2KGML;
import de.zbit.parser.Species;
import de.zbit.util.logging.LogUtil;

/**
 * Class to test the BioCartaTools
 * 
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class BioPax2KGMLTest {
  
  private BioPax2KGML bc = null;
  
  public BioPax2KGMLTest(){
    bc = new BioPax2KGML();
  }
  

  private void testCreateKGMLsFromBioCartaModel(String file) {
    // test for pathway gene ids
    //    bioCartafile Level 3(can be downloaded from http://pid.nci.nih.gov/download.shtml)
    Model m = bc.getModel(file);
    bc.createKGMLsFromModel(m);
  }
  
  private void testCreateKGMLsFromBioPaxFile(String file) {
    Species species = new Species("Homo sapiens", "_HUMAN", "human", "hsa",9606);    
    // test for pathway gene ids
    //    bioCartafile Level 3(can be downloaded from http://pid.nci.nih.gov/download.shtml)
    Model m = bc.getModel(file);
    bc.createKGMLsFromBioPaxFile(m, species, "BIOMD0000000201", 201, "Goldbeter2008_Somite_Segmentation_Clock_Notch_Wnt_FGF");
  }
  
  /**
   * method to test the {@link BioPax2KGML#getPathwaysWithGeneID(String, Model)}
   * Be carefull this method uses a {@link BioPax2KGML#getModel(String)} call where a local BioCarta file
   * of level 3 is needed. The file could be downloaded from http://pid.nci.nih.gov/download.shtml
   */
  private void testGetPathwaysWithGeneID() {
    String species = "human";    
    // test for pathway gene ids
    //    bioCartafile Level 3(can be downloaded from http://pid.nci.nih.gov/download.shtml)
    Model m = bc.getModel("C:/Users/buechel/Downloads/BioCarta.bp3.owl");
    for (BioPaxPathwayHolder pw : bc.getPathwaysWithEntrezGeneID(species, m)) {
      System.out.println(pw.getRDFid() + "\t" + pw.getName());
    }
  }
  
  /**
   * @param args
   * @throws FileNotFoundException
   */
  public static void main(String[] args) throws FileNotFoundException {
    LogUtil.initializeLogging(Level.CONFIG);
    
    BioPax2KGMLTest bft = new BioPax2KGMLTest();
    
    // Biocarta file with several pathways in one file
    String file = "C:/Users/buechel/Downloads/BioCarta.bp3.owl";
    bft.testCreateKGMLsFromBioCartaModel(file);
   
//    // Reactome file with several pathways in one file
//    String file = "C:/Users/buechel/Downloads/ReactomePathways/Homo sapiens.owl";
//    bft.testCreateKGMLsFromBioCartaModel(file);
    
//    // Singel Biomodels file of one pathway
//    String file = "C:/Users/buechel/Dropbox/Uni/BioPax-SBML-Projekt/BIOMD0000000201-biopax3.owl";
//    bft.testCreateKGMLsFromBioPaxFile("C:/Users/buechel/Dropbox/Uni/BioPax-SBML-Projekt/BIOMD0000000201-biopax3.owl");

    if(true)return;
   
    bft.testGetPathwaysWithGeneID();

  }

}
