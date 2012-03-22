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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Set;
import java.util.logging.FileHandler;
import java.util.logging.Filter;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level2.pathway;
import org.biopax.paxtools.model.level2.unificationXref;
import org.biopax.paxtools.model.level3.UnificationXref;

import de.zbit.kegg.KGMLWriter;
import de.zbit.kegg.parser.KeggParser;
import de.zbit.kegg.parser.pathway.Pathway;
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
  

  private void testCreateKGMLsFromBioCartaModel(String file, String destinationFolder, 
      boolean singleMode, boolean writeEntryExtended) {   
    BioPax2KGML.createKGMLsFromModel(file, destinationFolder, singleMode, writeEntryExtended);    
  }  
  
  private void testCreateKGMLsFromDirectory(String fileFolder, String destinationFolder, 
      boolean singleMode, boolean writeEntryExtended) {
    File f = new File(fileFolder);
    if (f.isDirectory()) {
      String[] files = f.list();
      for (String file : files) {
        BioPax2KGML.createKGMLsFromModel(fileFolder + file, destinationFolder, singleMode, writeEntryExtended);    
      }      
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
    if (m!=null && m.getLevel().equals(BioPAXLevel.L2)) {
//      for (BioPaxPathwayHolder pw : bc2.getPathwaysWithEntrezGeneID(species, m)) {
//        System.out.println(pw.getRDFid() + "\t" + pw.getName());
//      }
      System.out.println("up to known not implemented for level2");
    } else if (m!=null && m.getLevel().equals(BioPAXLevel.L3)){
      BioPaxL32KGML bc3 = new BioPaxL32KGML();
      for (BioPaxPathwayHolder pw : bc3.getPathwaysWithEntrezGeneID(species, m)) {
        System.out.println(pw.getRDFid() + "\t" + pw.getName());
      }
    } else if (m==null){
      log.log(Level.SEVERE,"Model = null.");
    }
    
  }
  
  /**
   * for extending existing KGML files with further information for the relations
   * @param fileList
   * @param file
   * @param writeEntryExtended
   */
  private static void createExtendedKGML(List<String> fileList, String file, boolean writeEntryExtended) {
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
    if(m!=null){
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
          KGMLWriter.writeKGML(p, fn, false);  
        }    
      }  
    } else {
      log.log(Level.SEVERE, "Could not continue, because the model is null.");
    }

      
  }
  

  private void parseAndWritePathway(String file, String destinationFolder, String pwName) {
    Model m = BioPax2KGML.getModel(file);
    if (m!=null){
      Pathway keggPW = null;
      if(m.getLevel().equals(BioPAXLevel.L2)){
        BioPaxL22KGML b22 = new BioPaxL22KGML();
        pathway pw = b22.getPathwayByName(m, pwName);      
        if(pw!=null)
          keggPW = b22.parsePathway(m, pw, b22.determineSpecies(pw.getORGANISM()));
      } else if(m.getLevel().equals(BioPAXLevel.L3)){
        BioPaxL32KGML b23 = new BioPaxL32KGML();
        org.biopax.paxtools.model.level3.Pathway pw = b23.getPathwayByName(m, pwName);
        if(pw!=null)
          keggPW = b23.parsePathway(m, pw, b23.determineSpecies(pw.getOrganism()));
      }
      
      if(keggPW != null){
        KGMLWriter.writeKGML(keggPW, destinationFolder + KGMLWriter.createFileName(keggPW), false);
      } else {
        System.out.println("keggPW is null.");
      }
    } else {
      log.log(Level.SEVERE, "Could not continue, because the model is null.");
    }
  }

  /**
   * print the pathwy list of a model
   * @param file
   */
  private void printPathwayList(String file) {
    Model m = BioPax2KGML.getModel(file);
    List<String> pws = null;
    if(m.getLevel().equals(BioPAXLevel.L2)){
      BioPaxL22KGML b22 = new BioPaxL22KGML();
      pws = b22.getListOfPathways(m);
    } else if(m.getLevel().equals(BioPAXLevel.L3)){
      BioPaxL32KGML b23 = new BioPaxL32KGML();
      pws = b23.getPathways(m);
    }
    
    if(pws!=null && pws.size()>0){
      for (String string : pws) {
        System.out.println(string);
      }
    }
    
  }  
  
/**
 * simple check to find out which database abbreviations are used
 * @param file
 */
  private void logUnificationXRefs(String file) {
   Model m = BioPax2KGML.getModel(file);
    if (m != null && m.getLevel().equals(BioPAXLevel.L2)) {
      Set<unificationXref> refs = m.getObjects(unificationXref.class);
      for (unificationXref ur : refs) {
        log.info(ur.getDB() + " - " + ur.getID());
      }
    } else if (m != null && m.getLevel().equals(BioPAXLevel.L3)) {
      Set<UnificationXref> refs = m.getObjects(UnificationXref.class);
      for (UnificationXref ur : refs) {
        log.info(ur.getDb() + " - " + ur.getId());
      }
    }
   
    
  }

  
  /**
   * @param args
   * @throws FileNotFoundException
   */
  public static void main(String[] args) throws FileNotFoundException {
    LogUtil.initializeLogging(Level.INFO);
    FileHandler h = null;
    try {
      h = new FileHandler("log.txt");
    } catch (SecurityException e1) {
      e1.printStackTrace();
    } catch (IOException e1) {
      e1.printStackTrace();
    }
    h.setFormatter(new OneLineFormatter());    
    LogUtil.addHandler(h, LogUtil.getInitializedPackages());
    
    FileHandler h2 = null;
    try {
      h2 = new FileHandler("unificationRefs.txt");
    } catch (SecurityException e1) {
      e1.printStackTrace();
    } catch (IOException e1) {
      e1.printStackTrace();
    }
    h2.setFormatter(new OneLineFormatter());
    h2.setFilter(new Filter() {

      /* (non-Javadoc)
       * @see java.util.logging.Filter#isLoggable(java.util.logging.LogRecord)
       */
      public boolean isLoggable(LogRecord record) {
        if ((record.getSourceClassName().equals(BioPax2KGMLTest.class.getName()) || record
            .getLoggerName().equals(BioPax2KGMLTest.class.getName()))
            && record.getLevel().equals(Level.INFO)) {

          return true;
        }
        return false;
      }
    });
    LogUtil.addHandler(h, LogUtil.getInitializedPackages());
    
    String fileFolder = System.getenv("FILE_FOLDER");
    
    BioPax2KGMLTest bft = new BioPax2KGMLTest();
    
    
    // Database testing of the parser
    
//    // INOH Pathway
//    String subINOH_SigL2 = "INOH Pathway/Signal_Level2/";
//    String subINOH_SigL3 = "INOH Pathway/Signal_Level3/";
//    String subINOH_MetL2 = "INOH Pathway/Metabolic_Level2/";
//    String subINOH_MetL3 = "INOH Pathway/Metabolic_Level3/";
//    
//    bft.testCreateKGMLsFromDirectory(fileFolder + subINOH_SigL2, fileFolder + "INOH Pathway/res_s2/", true, false);    
//    bft.testCreateKGMLsFromDirectory(fileFolder + subINOH_SigL3, fileFolder + "INOH Pathway/res_s3/", true, false);
//    bft.testCreateKGMLsFromDirectory(fileFolder + subINOH_MetL2, fileFolder + "INOH Pathway/res_m2/", true, false);
//    bft.testCreateKGMLsFromDirectory(fileFolder + subINOH_MetL3, fileFolder + "INOH Pathway/res_m3/", true, false);
//    
//    bft.logUnificationXRefs(fileFolder + subINOH_SigL2 + "GPCR_Adenosine_A2A_receptor.owl");
//    bft.logUnificationXRefs(fileFolder + subINOH_SigL3 + "Wnt_CAEEL.owl");
//    bft.logUnificationXRefs(fileFolder + subINOH_MetL2 + "Phenylalanine_degradation.owl");
//    bft.logUnificationXRefs(fileFolder + subINOH_MetL3 + "Steroids_metabolism.owl");

//    // NetPathway
//    String subNetPath = "NetPath Pathway/";
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subNetPath + 
//        "NetPath_3.owl", fileFolder + subNetPath, true, false);
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subNetPath + 
//        "NetPath_4.owl", fileFolder + subNetPath, true, false);
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subNetPath + 
//        "NetPath_11.owl", fileFolder + subNetPath, true, false);
//    
//    bft.logUnificationXRefs(fileFolder + subNetPath + "NetPath_3.owl");
//    bft.logUnificationXRefs(fileFolder + subNetPath + "NetPath_4.owl");
//    bft.logUnificationXRefs(fileFolder + subNetPath + "NetPath_11.owl");
    
//    // Panther
//    String subPanther = "Panther/owl/";
////    bft.testCreateKGMLsFromDirectory(fileFolder + subPanther, fileFolder + "Panther/res/", true, false);
//   
//    bft.logUnificationXRefs(fileFolder + subPanther + "Aminobutyrate_degradation.owl");
//    bft.logUnificationXRefs(fileFolder + subPanther + 
//        "Heterotrimeric_G_protein_signaling_pathway_Gi_and_Gs_mediated_pathway.owl");
//    bft.logUnificationXRefs(fileFolder + subPanther + "Oxidative_stress_response.owl");
//    bft.logUnificationXRefs(fileFolder + subPanther + "Wnt_signaling_pathway.owl");
    
//    // PathwayCommons
//    String subPathwayCommons = "PathwayCommons/";
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPathwayCommons + 
//        "biogrid.owl", fileFolder + subPathwayCommons, true, false);
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPathwayCommons + 
//        "cell-map.owl", fileFolder + subPathwayCommons, true, false);
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPathwayCommons + 
//        "hprd.owl", fileFolder + subPathwayCommons, true, false);
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPathwayCommons + 
//        "humancyc.owl", fileFolder + subPathwayCommons, true, false);
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPathwayCommons + 
//        "imid.owl", fileFolder + subPathwayCommons, true, false);
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPathwayCommons + 
//        "intact.owl", fileFolder + subPathwayCommons, true, false);
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPathwayCommons + 
//        "mint.owl", fileFolder + subPathwayCommons, true, false);
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPathwayCommons + 
//        "nci-nature.owl", fileFolder + subPathwayCommons, true, false);
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPathwayCommons + 
//        "reactome.owl", fileFolder + subPathwayCommons, true, false);
//  //TODO: just if I have time left overs
    
//    // PhosphoSitePlus
//    String subPhosphpSitePlus = "PhosphpSitePlus/";
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPhosphpSitePlus + 
//        "Kinase_substrates.owl", fileFolder + subPhosphpSitePlus, true, false);
//    bft.logUnificationXRefs(fileFolder + subPhosphpSitePlus + 
//        "Kinase_substrates.owl");
   
    // PID 
    String subPID = "PID_Pathways/";
    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPID + 
        "ceramidepathway.owl", 
        fileFolder + subPID, true, false);
//    
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPID + 
//        "alk1_2pathway_changed.owl", fileFolder + subPID, true, false);    
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPID + 
//        "parkinsonspathway.owl", fileFolder + subPID, true, false);    
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPID + 
//        "BioCarta.bp2.owl", fileFolder + subPID, false, false);    
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subPID + 
//        "BioCarta.bp3.owl", fileFolder + subPID, false, false);
//    
//    bft.logUnificationXRefs(fileFolder + subPID + "BioCarta.bp3.owl");
//    bft.logUnificationXRefs(fileFolder + subPID + "BioCarta.bp2.owl");
    
//    // Reactome
//    String subReactomeL2 = "ReactomePathways/Level2/";
//    String subReactomeL3 = "ReactomePathways/Level3/";
//    
//    bft.parseAndWritePathway(fileFolder + subReactomeL2 + "Homo sapiens.owl", 
//        fileFolder + subReactomeL2, "Sema4D mediated inhibition of cell attachment and migration");
//    bft.parseAndWritePathway(fileFolder + subReactomeL3 + "Homo sapiens.owl", 
//        fileFolder + subReactomeL3, "Sema4D mediated inhibition of cell attachment and migration");
//    
//    bft.logUnificationXRefs(fileFolder + subReactomeL2 + "Homo sapiens.owl");
//    bft.logUnificationXRefs(fileFolder + subReactomeL3 + "Homo sapiens.owl");
    
 
//    //SPIKE
//    String subSpike = "SPIKE/Spike2_LatestSpikeDB.owl";
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subSpike, 
//        fileFolder + "SPIKE/", true, false); 
//    bft.logUnificationXRefs(fileFolder + subSpike);
    

    // WikiPathways
//    String subWikiPath = "WikiPathway/WP1984_46412.owl";
//    bft.testCreateKGMLsFromBioCartaModel(fileFolder + subWikiPath, fileFolder + 
//        "WikiPathway/", true, false);    
//    bft.logUnificationXRefs(fileFolder + subWikiPath);

  }

}
