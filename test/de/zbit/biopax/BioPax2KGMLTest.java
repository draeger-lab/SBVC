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

import org.biopax.paxtools.model.BioPAXLevel;
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
  private BioPaxL32KGML bc3 = null;
  private BioPaxL22KGML bc2 = null;
  
  public BioPax2KGMLTest(){
    bc3 = new BioPaxL32KGML();
    bc2 = new BioPaxL22KGML();
  }
  

  private void testCreateKGMLsFromBioCartaModel(String file) {
    Model m = BioPax2KGML.getModel(file);
    if (m.getLevel().equals(BioPAXLevel.L2)) {
      bc2.createKGMLsFromModel(m);
    } else if (m.getLevel().equals(BioPAXLevel.L3)){
      bc3.createKGMLsFromModel(m);
    }
  }
  
  private void testCreateKGMLsFromBioPaxFile(String file) {
    Species species = new Species("Homo sapiens", "_HUMAN", "human", "hsa",9606);    
    // test for pathway gene ids    
    Model m = BioPax2KGML.getModel(file);
    if (m.getLevel().equals(BioPAXLevel.L2)) {
      bc2.createKGMLsFromBioPaxFile(m, species, "BIOMD0000000201", 201, "alk1_2pathway");
    } else if (m.getLevel().equals(BioPAXLevel.L3)){
      bc3.createKGMLsFromBioPaxFile(m, species, "BIOMD0000000201", 201, "alk1_2pathway");
    }
  }
  
  /**
   * method to test the {@link BioPaxL32KGML#getPathwaysWithGeneID(String, Model)}
   * Be carefull this method uses a {@link BioPaxL32KGML#getModel(String)} call where a local BioCarta file
   * of level 3 is needed. The file could be downloaded from http://pid.nci.nih.gov/download.shtml
   */
  private void testGetPathwaysWithGeneID(String file) {
    String species = "human";    
    Model m = BioPax2KGML.getModel(file);
    if (m.getLevel().equals(BioPAXLevel.L2)) {
//      for (BioPaxPathwayHolder pw : bc2.getPathwaysWithEntrezGeneID(species, m)) {
//        System.out.println(pw.getRDFid() + "\t" + pw.getName());
//      }
      System.out.println("up to known not implemented for level2");
    } else if (m.getLevel().equals(BioPAXLevel.L3)){
      for (BioPaxPathwayHolder pw : bc3.getPathwaysWithEntrezGeneID(species, m)) {
        System.out.println(pw.getRDFid() + "\t" + pw.getName());
      }
    }
    
  }
  
  private static void createExtendedKGML(List<String> fileList, String file) {
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

      /* (non-Javadoc)
       * @see java.util.logging.Filter#isLoggable(java.util.logging.LogRecord)
       */
      public boolean isLoggable(LogRecord record) {
        if ((record.getSourceClassName().equals(BioPaxL32KGML.class.getName()) || record
            .getLoggerName().equals(BioPaxL32KGML.class.getName()))
            && record.getLevel().equals(Level.INFO)
            && record.getSourceMethodName().equals("addRelationsToPathway")) {

          return true;
        }
        return false;
      }
    });
    LogUtil.addHandler(h, LogUtil.getInitializedPackages());

    Model m = BioPax2KGML.getModel(file);    

    for (String filename : fileList) {
      List<Pathway> pathways = null;
      try {
        pathways = KeggParser.parse(filename);
      } catch (Exception e) {
        log.log(Level.SEVERE, "Doof.");
      }
      BioPaxL32KGML bp2k = new BioPaxL32KGML();
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

    String keggFolder = System.getenv("KEGG_FOLDER");
    String bioCartaFile = System.getenv("BIOCARTA_FILE");
    String fileFolder = System.getenv("FILE_FOLDER");
    
    BioPax2KGMLTest bft = new BioPax2KGMLTest();    
//    bft.testCreateKGMLsFromBioPaxFile(fileFolder + "alk1_2pathway_changed.owl");
    bft.testCreateKGMLsFromBioCartaModel(fileFolder + "BioCarta.bp2.owl");
    
    if(true) return;
    
    List<String> fileList = new ArrayList<String>();
    fileList.add(keggFolder + "hsa05014.xml");
    fileList.add(keggFolder + "hsa04115.xml");
    fileList.add(keggFolder + "hsa05210.xml");
    fileList.add(keggFolder + "hsa05215.xml");
    fileList.add(keggFolder + "hsa05152.xml");
    fileList.add(keggFolder + "hsa04210.xml");
    fileList.add(keggFolder + "hsa04110.xml");
    fileList.add(keggFolder + "hsa05219.xml");
    fileList.add(keggFolder + "hsa05222.xml");
    fileList.add(keggFolder + "hsa05200.xml");
    fileList.add(keggFolder + "hsa05016.xml");
    fileList.add(keggFolder + "hsa04722.xml");
    fileList.add(keggFolder + "hsa05010.xml");
    fileList.add(keggFolder + "hsa05220.xml");
    fileList.add(keggFolder + "hsa05416.xml");
    fileList.add(keggFolder + "hsa05223.xml");
    fileList.add(keggFolder + "hsa05145.xml");
    fileList.add(keggFolder + "hsa05212.xml");
    fileList.add(keggFolder + "hsa04380.xml");
    fileList.add(keggFolder + "hsa05160.xml");
    createExtendedKGML(fileList, bioCartaFile);
    
    if(true)return;
    
   
    bft.testCreateKGMLsFromBioCartaModel(bioCartaFile);

    if(true)return;
   
    bft.testGetPathwaysWithGeneID(bioCartaFile);

  }  

}
