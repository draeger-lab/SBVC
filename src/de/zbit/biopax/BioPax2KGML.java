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
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biopax.paxtools.io.BioPAXIOHandler;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level2.pathway;
import org.biopax.paxtools.model.level3.Pathway;
import org.biopax.paxtools.util.IllegalBioPAXArgumentException;

import de.zbit.io.FileTools;
import de.zbit.mapper.GeneID2KeggIDMapper;
import de.zbit.mapper.GeneSymbol2GeneIDMapper;
import de.zbit.util.Species;
import de.zbit.util.Utils;

/**
 * This class is a base class to convert BioPax files and contains all methods
 * which are use for LEVEL 2 and LEVEL 3 converting
 * 
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public abstract class BioPax2KGML {

  public static final Logger log = Logger.getLogger(BioPax2KGML.class.getName());

  /**
   * default folder name for the KGMLs "pws"
   */
  static String defaultFolderName = "pws";
  static List<Species> allSpecies = null;

  static {
    // TODO: Replace this with a flat-text file in resources to support more
    // organisms.
    allSpecies = new ArrayList<Species>(4);
    allSpecies.add(new Species("Homo sapiens", "_HUMAN", "Human", "hsa", 9606));
    allSpecies.add(new Species("Mus musculus", "_MOUSE", "Mouse", "mmu", 10090));
    allSpecies.add(new Species("Rattus norvegicus", "_RAT", "Rat", "rno", 10116));
    allSpecies.add(new Species("Enterococcus faecalis", "_ENTFA", "Enterococcus", "efa", 226185));
  }

  /**
   * This variable must be set to false for normal biopax2kgml conversion.
   * 
   * If it is true, an existing pathway will be augmented with relations, and
   * relation subtypes. NO reactions and NO entries are added!!!
   */
  boolean augmentOriginalKEGGpathway = false;
  boolean addSelfReactions = false;
  int newAddedRelations = 0;
  int selfRelation = 0;
  int addedSubTypes = 0;

  /**
   * number which is used to determine a pathway id, if it is not possible to
   * exclude the id from the BioCarta file
   */
  int keggPathwayNumberCounter = 100000;

  /**
   * this variable is used to determine the kegg id of an entry
   */
  int keggEntryID = 0;

  /**
   * this variable is used to determine the kegg id of an entry
   */
  int keggUnknownNo = 0;

  /**
   * this variable is used to determine the kegg reaction id
   */
  int keggReactionID = 0;

  /**
   * default organism for KEGG parsing - "hsa"
   */
  String organism = "hsa";

  /**
   * undefined, if we have no gene id to set the kegg name of an entry we use
   * this name
   */
  String keggUnknownName = "unknown";

  /**
   * mapper to map gene symbols to gene ids
   */
  protected static GeneSymbol2GeneIDMapper geneSymbolMapper = null;

  /**
   * mapper to map gene ids to KEGG ids
   */
  protected static GeneID2KeggIDMapper geneIDKEGGmapper = null;

  /**
   * transforms a set to a map. The key is a RDFId and the value the
   * corresponding object
   * 
   * @param <T>
   * @param set
   *          to convert
   * @return the converted map
   */
  protected static <T extends BioPAXElement> Map<String, T> getMapFromSet(Set<T> set) {
    Map<String, T> map = new HashMap<String, T>();
    for (T elem : set) {
      map.put(elem.getRDFId(), elem);
    }

    return map;
  }

  /**
   * 
   * @return the new KEGG unknown "unknownx", whereas x is set to the
   *         {@link BioPax2KGML#keggUnknownNo}.{@link BioPax2KGML#keggUnknownNo}
   *         is incremented after this step
   */
  protected String getKEGGUnkownName() {
    return keggUnknownName + String.valueOf(++keggUnknownNo);
  }

  /**
   * 
   * @return the new KEGG reaction name "rn:unknownx", whereas x is set to the
   *         {@link BioPax2KGML#keggReactionID}.
   *         {@link BioPax2KGML#keggReactionID} is augmented after this step
   */
  protected String getReactionName() {
    return "rn:unknown" + String.valueOf(++keggReactionID);
  }

  /**
   * mapps an entered gene id to a kegg id, if this is not possible the species
   * abbreviation:geneID is returned
   * 
   * @param mapper
   * @return
   */
  protected String mapGeneIDToKEGGID(Integer geneID, Species species) {
    String keggName = null;
    try {
      keggName = geneIDKEGGmapper.map(geneID);
    } catch (Exception e) {
      log.log(Level.WARNING, "Could not map geneid: '" + geneID.toString() + "' to a KEGG id, "
          + "'speciesAbbreviation:geneID will be used instead.", e);
    }

    if (keggName == null) {
      keggName = species.getKeggAbbr() + ":" + geneID.toString();
    }

    return keggName;
  }

  /**
   * The rdfID is in the format: http://pid.nci.nih.gov/biopaxpid_9717
   * 
   * From this id the number is excluded and used as pathway number, if this is
   * not possible the {@link BioPaxL22KGML#keggPathwayNumberCounter} is used and
   * incremented
   * 
   * @param rdfId
   * @return
   */
  protected int getKeggPathwayNumber(String rdfId) {
    int posUnderscore = rdfId.indexOf('_');
    if (posUnderscore > -1 && posUnderscore <= rdfId.length()) {
      try {
        return Integer.parseInt(rdfId.substring(posUnderscore + 1));
      } catch (Exception e) {
        return keggPathwayNumberCounter++;
      }
    }

    return keggPathwayNumberCounter++;
  }

  /**
   * creates out of the pathwayName the KEGG no which is need to describe the
   * pathway
   * 
   * @param pathwayName
   * @return
   */
  protected int determineKEGGPathwayNumber(String pathwayName) {
    // TODO: better idea instead of using hashCode???
    return pathwayName.hashCode();
  }

  /**
   * Converts the inputStream of an owl file containing BioPAX entries
   * 
   * @param io
   * @return
   */
  public static Model getModel(InputStream io) {
    BioPAXIOHandler handler = new SimpleIOHandler();
    Model m = null;
    try {
      m = handler.convertFromOWL(io);
    } catch (IllegalBioPAXArgumentException e) {
      log.log(Level.SEVERE, "Could not read model!", e);
    }
    return m;
  }

  /**
   * The {@link BioPax2KGML}{@link #geneSymbolMapper} and
   * {@link BioPax2KGML#geneIDKEGGmapper} are initialized for the entered
   * species
   * 
   * @param species
   */
  public static void initalizeMappers(Species species) {
    try {
      geneSymbolMapper = new GeneSymbol2GeneIDMapper(species.getCommonName());
    } catch (IOException e) {
      log.log(Level.SEVERE, "Could not initalize mapper for species '" + species.toString() + "'!",
          e);
      System.exit(1);
    }

    try {
      geneIDKEGGmapper = new GeneID2KeggIDMapper(species);
    } catch (IOException e) {
      log.log(Level.SEVERE, "Error while initializing gene id to KEGG ID mapper for species '"
          + species.toString() + "'.", e);
      System.exit(1);
    }
  }

  /**
   * @return a unique {@link BioPaxL22KGML#keggEntryID}.
   */
  protected int getKeggEntryID() {
    keggEntryID++;
    return keggEntryID;
  }

  /**
   * determines the link for the pathway image
   * 
   * @param species
   * @param pathway
   * @param keggPW
   */
  public void addImageLinkToKEGGpathway(Species species, String pathwayName,
      de.zbit.kegg.parser.pathway.Pathway keggPW) {
    String linkName = pathwayName;
    if (!linkName.equals("") && linkName.contains("pathway")) {
      linkName = linkName.replace("pathway", "Pathway");

      if (species.getKeggAbbr().equals("hsa")) {
        keggPW.setLink("http://www.biocarta.com/pathfiles/h_" + linkName + ".asp");
        keggPW.setImage("http://www.biocarta.com/pathfiles/h_" + linkName + ".gif");
      } else if (species.getKeggAbbr().equals("mmu")) {
        keggPW.setLink("http://www.biocarta.com/pathfiles/m_" + linkName + ".asp");
        keggPW.setImage("http://www.biocarta.com/pathfiles/m_" + linkName + ".gif");
      }
    }
  }

  /**
   * Creates a folder depending on the {@link BioPax2KGML#defaultFolderName} and
   * the {@link BioPAXLevel}
   * 
   * @param level
   * @return the folderName
   */
  protected static String createDefaultFolder(BioPAXLevel level) {
    String folderName = defaultFolderName + level.toString() + "/";
    if (!new File(folderName).exists()) {
      boolean success = (new File(folderName)).mkdir();
      if (success) {
        log.log(Level.SEVERE, "Could not create directory '" + folderName + "'");
        System.exit(1);
      }
    }

    return folderName;
  }

  /**
   * Calls the method {@link BioPax2KGML#getModel(InputStream) for an entered
   * owl file}
   * 
   * @param file
   * @return Model
   */
  public static Model getModel(String file) {
    InputStream io = null;
    try {
      io = new FileInputStream(new File(file));
    } catch (FileNotFoundException e) {
      log.log(Level.SEVERE, "Could not parse file: " + file + ".", e);
    }

    log.log(Level.CONFIG, "Model sucessfully created");
    return getModel(io);
  }

  /**
   * Creates for an entered {@link Model} the corresponding KEGG pathways
   * 
   * @param m
   */
  public static void createKGMLsFromModel(String fileName, String destinationFolder,
      boolean singleMode, boolean writeEntryExtended) {
    Model m = BioPax2KGML.getModel(fileName);
    if (m!=null){
      File f = new File(fileName);
      if (destinationFolder == null || destinationFolder.isEmpty()) {
        destinationFolder = createDefaultFolder(m.getLevel());
      }

      if (m.getLevel().equals(BioPAXLevel.L2)) {
        BioPaxL22KGML bp = new BioPaxL22KGML();
        Set<pathway> pathways = m.getObjects(pathway.class);
        // if (pathways!=null && pathways.size()>0) {
        if (!singleMode) {
          bp.createKGMLsForPathways(m, destinationFolder, pathways, false);
        } else {
          bp.createKGMLForBioPaxFile(m, FileTools.removeFileExtension(f.getName()),
              destinationFolder, false);
        }
      } else if (m.getLevel().equals(BioPAXLevel.L3)) {
        BioPaxL32KGML bp = new BioPaxL32KGML();
        Set<Pathway> pathways = m.getObjects(Pathway.class);
        // if (pathways!=null && pathways.size()>0) {
        if (!singleMode) {
          bp.createKGMLsForPathways(m, destinationFolder, pathways, writeEntryExtended);
        } else {
          bp.createKGMLForBioPaxFile(m, FileTools.removeFileExtension(f.getName()),
              destinationFolder, writeEntryExtended);
        }
      } else {
        log.log(Level.SEVERE, "Unkown BioPax Level '" + m.getLevel().toString()
            + "' is not supported.");
        System.exit(1);
      }
    } else {
      log.log(Level.SEVERE, "Could not continue, because the model is null.");
    }
    
  }

  /**
   * @param
   * @return name of the set without blanks
   */
  protected String getNameWithoutBlanks(Set<String> names) {
    String name = "";

    if (names != null && names.size() > 0) {
      List<String> names2 = Utils.getListOfCollection(names);
      for (int i = names2.size() - 1; i > 0; i--) {
        name = names2.get(i);
        if (name.length() > 0 && !name.contains(" "))
          return name;
      }
    }

    return name;
  }
}
