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

import de.zbit.mapper.GeneID2KeggIDMapper;
import de.zbit.mapper.GeneSymbol2GeneIDMapper;
import de.zbit.parser.Species;

/**
 * This class is a base class to convert BioPax elements
 * 
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public abstract class BioPax2KGML {

  public static final Logger log = Logger.getLogger(BioPax2KGML.class.getName());

  static List<Species> allSpecies = null;
  
  static {
    // TODO: Replace this with a flat-text file in resources to support more
    // organisms.
    allSpecies = new ArrayList<Species>(3);
    allSpecies.add(new Species("Homo sapiens", "_HUMAN", "Human", "hsa", 9606));
    allSpecies.add(new Species("Mus musculus", "_MOUSE", "Mouse", "mmu", 10090));
    allSpecies.add(new Species("Rattus norvegicus", "_RAT", "Rat", "rno", 10116));
    allSpecies.add(new Species("Enterococcus faecalis", "_ENTFA", "Enterococcus", "efa", 226185));
  }
  
  /**
   * This variable must be set to false for normal biopax2kgml conversion. 
   * 
   * If it is true, an existing pathway will be augmented with relations, and relation subtypes.
   * NO reactions and NO entries are added!!!
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
   * Converts the inputStream of an owl file containing BioPAX entries
   * 
   * @param io
   * @return
   */
  public static Model getModel(InputStream io) {
    BioPAXIOHandler handler = new SimpleIOHandler();
    return handler.convertFromOWL(io);
  }  
  
  
  /**
   * The {@link BioPaxL32KGML}{@link #geneSymbolMapper} and
   * {@link BioPaxL32KGML#geneIDKEGGmapper} are initalized for the entered species
   * 
   * @param species
   */
  public static void initalizeMappers(Species species) {
    try {
      geneSymbolMapper = new GeneSymbol2GeneIDMapper(species.getCommonName());
    } catch (IOException e) {
      log.log(Level.SEVERE, "Could not initalize mapper for species '" + species.toString() 
          + "'!", e);
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
   * 
   * @return a unique {@link BioPaxL22KGML#keggEntryID}.
   */
  protected int getKeggEntryID() {
    keggEntryID++;
    return keggEntryID;
  }
  
  /**
   * Calls the method {@link BioPaxL32KGML#getModel(InputStream) for an entered
   * owl file}
   * 
   * @param file
   * @return
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

  public void createKGMLsFromModel(Model m) {
    String folder = "pws" + m.getLevel().toString() + "/";
    if (m.getLevel().equals(BioPAXLevel.L2)){
      BioPaxL22KGML bp = new BioPaxL22KGML();
      bp.createKGMLsFromModel(m, folder);
    } else if (m.getLevel().equals(BioPAXLevel.L3)){
      BioPaxL32KGML bp = new BioPaxL32KGML();
      bp.createKGMLsFromModel(m, folder);
    }  
  }
  

  
  /**
   * methods for model and file parsing
   * 
   */
  
  abstract public void createKGMLsFromModel(Model m, String folder);

  
}
