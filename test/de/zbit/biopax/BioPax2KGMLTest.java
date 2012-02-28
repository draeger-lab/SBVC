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
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.FileHandler;
import java.util.logging.Filter;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

import org.biopax.paxtools.model.Model;

import de.zbit.kegg.KGMLWriter;
import de.zbit.kegg.parser.KeggParser;
import de.zbit.kegg.parser.pathway.Pathway;
import de.zbit.parser.Species;
import de.zbit.util.logging.LogUtil;
import de.zbit.util.logging.OneLineFormatter;

/**
 * Class to test the BioCartaTools
 * 
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class BioPax2KGMLTest {
  
  public static final Logger log = Logger.getLogger(BioPax2KGMLTest.class.getName());
  
  private BioPax2KGML bc = null;
  
  public BioPax2KGMLTest(){
    bc = new BioPax2KGML();
  }
  

  private void testCreateKGMLsFromBioCartaModel(String file) {
    // test for pathway gene ids
    //    bioCartafile Level 3(can be downloaded from http://pid.nci.nih.gov/download.shtml)
    Model m = BioPaxConverter.getModel(file);
    bc.createKGMLsFromModel(m);
  }
  
  private void testCreateKGMLsFromBioPaxFile(String file) {
    Species species = new Species("Homo sapiens", "_HUMAN", "human", "hsa",9606);    
    // test for pathway gene ids
    //    bioCartafile Level 3(can be downloaded from http://pid.nci.nih.gov/download.shtml)
    Model m = BioPaxConverter.getModel(file);
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
    Model m = BioPaxConverter.getModel("C:/Users/buechel/Downloads/BioCarta.bp3.owl");
    for (BioPaxPathwayHolder pw : bc.getPathwaysWithEntrezGeneID(species, m)) {
      System.out.println(pw.getRDFid() + "\t" + pw.getName());
    }
  }
  
  private static void createExtendedKGML() {
    FileHandler h = null;
    try {
      h = new FileHandler("relationsAdded.txt");
    } catch (SecurityException e1) {
      e1.printStackTrace();
    } catch (IOException e1) {
      e1.printStackTrace();
    }
    h.setFormatter(new OneLineFormatter());
    h.setFilter(new Filter() {

      @Override
      public boolean isLoggable(LogRecord record) {
        if ((record.getSourceClassName().equals(BioPax2KGML.class.getName()) || record
            .getLoggerName().equals(BioPax2KGML.class.getName()))
            && record.getLevel().equals(Level.INFO)
            && record.getSourceMethodName().equals("addRelationsToPathway")) {

          return true;
        }
        return false;
      }
    });
    LogUtil.addHandler(h, LogUtil.getInitializedPackages());

    Model m = BioPaxConverter.getModel("C:/Users/buechel/Downloads/BioCarta.bp3.owl");
    List<String> fileList = new ArrayList<String>();
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05014.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa04115.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05210.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05215.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05152.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa04210.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa04110.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05219.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05222.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05200.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05016.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa04722.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05010.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05220.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05416.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05223.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05145.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05212.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa04380.xml");
    fileList.add("C:/Users/buechel/Desktop/KGML_hsa/non-metabolic/hsa/hsa05160.xml");

    for (String filename : fileList) {
      List<Pathway> pathways = null;
      try {
        pathways = KeggParser.parse(filename);
      } catch (Exception e) {
        log.log(Level.SEVERE, "Doof.");
      }
      BioPax2KGML bp2k = new BioPax2KGML();
      for (Pathway p : pathways) {
        bp2k.addRelationsToPathway(p, m);
        String fn = filename.replace(".xml", "_extended.xml");
        KGMLWriter.writeKGML(p, fn);  
      }    
    }    
  }
  
  /**
   * @param args
   * @throws FileNotFoundException
   */
  public static void main(String[] args) throws FileNotFoundException {
    LogUtil.initializeLogging(Level.INFO);
    
    createExtendedKGML();
    
    if(true)return;
    
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
