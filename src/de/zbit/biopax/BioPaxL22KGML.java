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

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level2.ControlType;
import org.biopax.paxtools.model.level2.InteractionParticipant;
import org.biopax.paxtools.model.level2.bioSource;
import org.biopax.paxtools.model.level2.biochemicalReaction;
import org.biopax.paxtools.model.level2.catalysis;
import org.biopax.paxtools.model.level2.complex;
import org.biopax.paxtools.model.level2.complexAssembly;
import org.biopax.paxtools.model.level2.control;
import org.biopax.paxtools.model.level2.conversion;
import org.biopax.paxtools.model.level2.entity;
import org.biopax.paxtools.model.level2.interaction;
import org.biopax.paxtools.model.level2.modulation;
import org.biopax.paxtools.model.level2.pathway;
import org.biopax.paxtools.model.level2.pathwayComponent;
import org.biopax.paxtools.model.level2.physicalEntity;
import org.biopax.paxtools.model.level2.physicalEntityParticipant;
import org.biopax.paxtools.model.level2.process;
import org.biopax.paxtools.model.level2.protein;
import org.biopax.paxtools.model.level2.relationshipXref;
import org.biopax.paxtools.model.level2.rna;
import org.biopax.paxtools.model.level2.sequenceEntity;
import org.biopax.paxtools.model.level2.smallMolecule;
import org.biopax.paxtools.model.level2.transport;
import org.biopax.paxtools.model.level2.transportWithBiochemicalReaction;
import org.biopax.paxtools.model.level3.Complex;

import de.zbit.kegg.KGMLWriter;
import de.zbit.kegg.parser.pathway.Entry;
import de.zbit.kegg.parser.pathway.EntryType;
import de.zbit.kegg.parser.pathway.Graphics;
import de.zbit.kegg.parser.pathway.Pathway;
import de.zbit.kegg.parser.pathway.Reaction;
import de.zbit.kegg.parser.pathway.ReactionComponent;
import de.zbit.kegg.parser.pathway.ReactionType;
import de.zbit.kegg.parser.pathway.Relation;
import de.zbit.kegg.parser.pathway.RelationType;
import de.zbit.kegg.parser.pathway.SubType;
import de.zbit.kegg.parser.pathway.ext.GeneType;
import de.zbit.mapper.GeneSymbol2GeneIDMapper;
import de.zbit.parser.Species;
import de.zbit.util.ProgressBar;
import de.zbit.util.Utils;

/**
 * This class works with PaxTools. It is used to fetch information out of a
 * level 2 BioCarta files. Example files could be downloaded from
 * http://pid.nci.nih.gov/download.shtml
 * 
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class BioPaxL22KGML extends BioPax2KGML {
  
  public static final Logger log = Logger.getLogger(BioPaxL22KGML.class.getName());
  
  /**
   * parse an BioPax file by calling 
   * {@link BioPaxL22KGML#createKGMLForBioPaxFile(Model, String, String)}
   * the folder is determined with 
   * {@link BioPax2KGML#createDefaultFolder(org.biopax.paxtools.model.BioPAXLevel)}
   * @param m
   * @param fileName
   */
  public void createKGMLForBioPaxFile(Model m, String fileName) {
    String folder = createDefaultFolder(m.getLevel());
    createKGMLForBioPaxFile(m, fileName, folder);
  }
  /**
   * The methods parse a BioPax file which contains no <bp>Pathway: ....</bp> tag
   * 
   * @param m
   * @param pathwayName used as pathway name and title 
   * @param folder where the KGML is saved
   */
  public void createKGMLForBioPaxFile(Model m, String pathwayName, String folder) {
    // determine the organism
    String defaultSpecies = "Homo sapiens", newSpecies = defaultSpecies;
    Species species = new Species(defaultSpecies, "_HUMAN", "human", "hsa", 9606);
    
    Set<bioSource> orgs = m.getObjects(bioSource.class);
    if (orgs != null && orgs.size() > 0) {
      bioSource org = orgs.iterator().next();
      if (org != null && org.getNAME() != null) {
        newSpecies = org.getNAME();
        if (!newSpecies.equals(defaultSpecies)) {
          species = Species.search(allSpecies, newSpecies, Species.COMMON_NAME);
        }
      }
      log.info("Determined pathway species '" + species.getCommonName() + "'.");
    } else {
      log.info("No pathway species could be determined, default species '" 
          + species.getCommonName() + "' is used.");
    }
    
    initalizeMappers(species);    

    int pathNo = determineKEGGPathwayNumber(pathwayName);
    de.zbit.kegg.parser.pathway.Pathway keggPW = new de.zbit.kegg.parser.pathway.Pathway(
        pathwayName, species.getKeggAbbr(), pathNo, pathwayName);

    for (entity entity : m.getObjects(entity.class)) {
      parseEntity(entity, keggPW, m, species);
    }

    String fileName = folder + KGMLWriter.createFileName(keggPW);
    KGMLWriter.writeKGML(keggPW, fileName);
  } 

  
  /**
   * This method creates for each pathway in the model a KGML file with the
   * pathway name and saves the pathways in an default folder see 
   * {@link BioPax2KGML#createDefaultFolder(org.biopax.paxtools.model.BioPAXLevel)}
   * 
   * @param m
   */
  public void createKGMLsForPathways(Model m, Set<pathway> pathways) {
    String folder = createDefaultFolder(m.getLevel());
    createKGMLsForPathways(m, folder, pathways);
  }
  
  /**
   * This method creates for each pathway in the model a KGML file with the
   * pathway name
   * 
   * @param m
   */
  public void createKGMLsForPathways(Model m, String folder, Set<pathway> pathways) {
    log.info("Creating KGML files for pathways.");
    Collection<Pathway> keggPWs = parsePathways(m, pathways);
    
    for (Pathway keggPW : keggPWs) {
      String fileName = folder + KGMLWriter.createFileName(keggPW);
      KGMLWriter.writeKGML(keggPW, fileName);
    }    
  }

  /**
   * this method parses the biopax pathway by firstly determining the pathway species and
   * then parsing the single pathway
   * @param m
   * @param pathways
   * @return
   */
  private Collection<Pathway> parsePathways(Model m, Set<pathway> pathways) {
    String oldSpecies = "Homo sapiens", newSpecies = oldSpecies;
    Species species = new Species(oldSpecies, "_HUMAN", "human", "hsa", 9606);
    initalizeMappers(species);
    
    Collection<Pathway> keggPWs = new ArrayList<Pathway>();

    for (pathway pathway : pathways) {
      // determine the pathway organism - it's done here to save time, while initializing the mappers
      bioSource pwOrg = pathway.getORGANISM();
      if (pwOrg != null) {          
        if (pwOrg.getNAME() != null){
          newSpecies = pwOrg.getNAME();
          if (!newSpecies.equals(oldSpecies)){
            species = Species.search(allSpecies, newSpecies, Species.COMMON_NAME);
            initalizeMappers(species);
            oldSpecies = newSpecies;
          }
        }
      }
      keggPWs.add(parsePathway(m, pathway, species));
    }
    return keggPWs;
  }
  
  /**
   * parses the biopax pathway
   * @param m
   * @param pathway
   * @param species
   * @return
   */
  private Pathway parsePathway(Model m, pathway pathway, Species species) {
    // create the pathway
    int number = getKeggPathwayNumber(pathway.getRDFId());
    de.zbit.kegg.parser.pathway.Pathway keggPW = new de.zbit.kegg.parser.pathway.Pathway(
        "biocarta:" + number, species.getKeggAbbr(), number,
        pathway.getNAME());

    addImageLinkToKEGGpathway(species, pathway.getNAME(), keggPW);

    for (pathwayComponent interaction : pathway.getPATHWAY_COMPONENTS()) {
//      if (pathwayStep.class.isAssignableFrom(interaction.getClass())){
//        
//      } else if (interaction.class.isAssignableFrom(interaction.getClass())){
        parseInteraction((interaction) interaction, keggPW, m, species);  
//      } else if (pathway.class.isAssignableFrom(interaction.getClass())){
//        parsePathway(m, (pathway) interaction, species);  
//      } else {
//        log.log(Level.SEVERE, "Could not parse: '" + interaction.getModelInterface() + "'.");
//      }      
    }

    return keggPW;
  }
  /**
   * deteremines the gene ids of the elements in a pathway
   * 
   * This method is not so clean should be rewritten, becuase in the method
   * {@link BioPaxL22KGML#getEntrezGeneIDsForPathways(List, String, Model)}
   * complexes are not treated right
   * 
   * @param pathways
   * @param species
   * @param m
   * @return
   */
  public List<BioPaxPathwayHolder> getEntrezGeneIDsForPathways(
      List<BioPaxPathwayHolder> pathways, String species, Model m) {
    log.info("Start parsing gene ids.");
    ProgressBar bar = new ProgressBar(pathways.size());

    Map<String, relationshipXref> xrefs = getMapFromSet(m.getObjects(relationshipXref.class));
    for (BioPaxPathwayHolder pw : pathways) {
      // if (pw.getRDFid().equals("http://pid.nci.nih.gov/biopaxpid_9796"))
      // {//TODO: is necessary to uncomment!!!!
      log.log(Level.FINER, "Pathway: " + pw.getPathwayName() + ": " + pw.getNoOfEntities());
      Set<entity> pwEntities = new HashSet<entity>();
      for (BioPAXElement entity : pw.entities) {
        pwEntities.addAll(getEntitiesWithName((entity)entity));
        if (!(Pathway.class.isAssignableFrom(entity.getClass())))// entity
                                                                 // instanceof
                                                                 // Pathway ||
                                                                 // entity
                                                                 // instanceof
                                                                 // PathwayImpl))
          log.log(Level.FINER, "--Input: " + entity.getRDFId() + "\t" + entity.getModelInterface());
      }

      Map<entity, Integer> geneIDs = getEntrezGeneIDs(pwEntities, species, xrefs);
      for (java.util.Map.Entry<entity, Integer> entity : geneIDs.entrySet()) {
        log.log(Level.FINER, "----res: " + entity.getKey() + " " + entity.getValue());
        pw.addGeneID(entity.getValue());
      }
      // }//TODO: is necessary to uncomment!!!!
      bar.DisplayBar();
    }

    return pathways;
  }
  
  /**
   * The method returns the smallest entity having a name, i.e. a gene symbol,
   * which could be parsed
   * 
   * @param entity
   * @return Collection containing {@link entity}s having a name and are not
   *         instance of a complex or ComplexAssembly
   */
  protected Collection<? extends entity> getEntitiesWithName(entity entity) {
    Set<entity> resEntities = new HashSet<entity>();
    String name = entity.getNAME();

    if (!name.isEmpty() && !(pathway.class.isAssignableFrom(entity.getClass()))) {
      if (complex.class.isAssignableFrom(entity.getClass())) {
        complex c = (complex) entity;
        for (physicalEntityParticipant pe : c.getCOMPONENTS()) {
          resEntities.addAll(getEntitiesWithName(pe.getPHYSICAL_ENTITY()));
        }
      } else if (complexAssembly.class.isAssignableFrom(entity.getClass())) {
        complexAssembly c = (complexAssembly) entity;
        for (InteractionParticipant pe : c.getPARTICIPANTS()) {
          resEntities.addAll(getEntitiesWithName(((physicalEntityParticipant)pe).getPHYSICAL_ENTITY()));
        }

      } else {
        resEntities.add(entity);
      }
    } 
    return resEntities;
  }

  /**
   * This method maps to all gene symbols of the entered entities the
   * corresponding gene id It's an advantage to preprocess the entities to
   * exclusively having entities with a name call therefore the method
   * {@link BioPaxL22KGML#getEntitiesWithName(Eetity)} This method also considers
   * 
   * @link {@link Complex} classes. NOTE: the method is not so clean and should
   *       be rewritten!!
   * 
   * @param entities
   * @param species
   * @param xrefs
   * @return
   */
  private Map<entity, Integer> getEntrezGeneIDs(Set<entity> entities, String species,
      Map<String, relationshipXref> xrefs) {
    Map<entity, Integer> geneIDs = new HashMap<entity, Integer>();

    // searching for gene ids
    Map<String, entity> geneSymbols = new HashMap<String, entity>();
    for (entity entity : entities) {
      Integer geneID = null;
      String id = entity.getRDFId();

      // gene id in xrefs?
      if (xrefs.containsKey(id)) {
        geneID = getEntrezGeneIDFromDBxref(xrefs.get(id));
      }

      if (geneID == null) {
        // we have to search the gene id with the gene symbol, adding symbol to
        // the gene symbol set
        String names = entity.getNAME();
        if (names != null && !names.isEmpty()) {
          
            if (names.contains("/")) {
              String[] split = names.split("/");
              for (String rs : split) {
                if (!geneSymbols.containsKey(rs))
                  geneSymbols.put(rs, entity);
              }
            } else if (!geneSymbols.containsKey(names)) {
              geneSymbols.put(names, entity);
            }
          
        }
      } else if (!geneIDs.containsKey(geneID)) {
        // gene id found directly
        log.log(Level.FINER, "found: " + entity.getRDFId() + " " + geneID);
        geneIDs.put(entity, geneID);
      }
    }

    // getting the gene ids which could not be found directly
    if (geneSymbols.size() > 0) {
      geneIDs.putAll(getEntrezGeneIDOverGeneSymbol(geneSymbols));
    }

    return geneIDs;
  }
  
  /**
   * firstly set {@link BioPaxL22KGML#augmentOriginalKEGGpathway} to true, 
   * secondly a model from biocarta is created, and
   * thirdly all relations are added to p if entry1 and entry2 of the relation are in p too. 
   * 
   * @param p
   * @return
   */
  public de.zbit.kegg.parser.pathway.Pathway addRelationsToPathway(
      de.zbit.kegg.parser.pathway.Pathway p, Model m) {
    augmentOriginalKEGGpathway = true;
    int relationsBegin = p.getRelations().size();
    newAddedRelations = 0;
    selfRelation = 0;
    addedSubTypes = 0;
    Set<pathway> pathways = m.getObjects(pathway.class);
    Species species = Species.search(allSpecies, p.getOrg(), Species.KEGG_ABBR);
    if(species != null){
      initalizeMappers(species);
      
      for (pathway pathway : pathways) {    
        // determine the pathway organism
        bioSource pwOrg = pathway.getORGANISM();  
        if (pwOrg != null) {
          if (pwOrg.getNAME() != null) {
            String pathwaySpecies = pwOrg.getNAME();
            if (pathwaySpecies.equals(species.getScientificName())) {
              for (pathwayComponent interaction : pathway.getPATHWAY_COMPONENTS()) {
                parseInteraction((interaction) interaction, p, m, species);
              }      
            } else {
              log.log(Level.WARNING, "No additional information available for species '" + 
                  species.getScientificName() + "'.");
            }
          }          
        }
        
      }
    } else {
      log.log(Level.SEVERE, "It was not possible to initialize the pathway species '" +
          p.getOrg()+ "'.");
    }
    
    log.log(Level.INFO, (p.getRelations().size()-relationsBegin) + "|" + newAddedRelations + " new relations are added to pathway '" + p.getName() 
        + "', (" + (newAddedRelations-selfRelation) + " relations with different kegg identifiers, "
        + selfRelation + " self relations), and " 
        + addedSubTypes + " subtypes are added to existing relations.");
    
    return p;
  }


  /**
   * parse a BioPax entity element
   * 
   * @param entity
   * @param keggPW
   * @param mapper
   * @param m
   * @param species
   * @return
   */
  private Entry parseEntity(entity entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species) {
    Entry keggEntry = null;
    if (physicalEntity.class.isAssignableFrom(entity.getClass())) {
      keggEntry = parsePhysicalEntity((physicalEntity) entity, keggPW, m, species);
    } else if (pathway.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.map, null, ",", null);
    } else if (interaction.class.isAssignableFrom(entity.getClass())) {
      parseInteraction((interaction)entity, keggPW, m, species);
    } else {
      log.severe("Unknonw entity type: " + entity.getModelInterface() + "-" + entity.getRDFId());
      System.exit(1);
    }

    return keggEntry;
  }

  /**
   * parse a BioPax PhysicalEntity element
   * 
   * @param entity
   * @param keggPW
   * @param m
   * @param species
   * @return
   */
  private Entry parsePhysicalEntity(physicalEntity entity,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    Entry keggEntry = null;

    if (complex.class.isAssignableFrom(entity.getClass())) {
      List<Integer> components = createComplexComponentList(((complex) entity).getCOMPONENTS(),
          keggPW, m, species);
      keggEntry = createKEGGEntry((entity) entity, keggPW, m, species, EntryType.group, null, "/",
          components);
    } else if (sequenceEntity.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, GeneType.dna, ",",
          null);
    } else if (protein.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, GeneType.protein,
          ",", null);
    } else if (rna.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, GeneType.rna, ",",
          null);
    } else if (smallMolecule.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.compound, GeneType.unknown,
          ",", null);
    } else {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.other, GeneType.unknown,
          ",", null);
    }

    return keggEntry;
  }

  /**
   * Creates a list of all complex entities. Each entity is checked if it
   * already exists in the KEGG pathway and if not it is created. to
   * 
   * @param set
   * @param keggPW
   * @param m
   * @param species
   * @return
   */
  private List<Integer> createComplexComponentList(Set<physicalEntityParticipant> set,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    List<Integer> components = new ArrayList<Integer>();

    for (physicalEntityParticipant physicalEntity : set) {
      Entry keggEntry = parsePhysicalEntity(physicalEntity.getPHYSICAL_ENTITY(), keggPW, m, species);
      if (keggEntry != null)
        components.add(keggEntry.getId());
    }

    return components;
  }


  /**
   * This method maps the gene id to the gene symbol of the entered entity. This
   * method is not able to map complexes call for that the method
   * {@link BioPaxL22KGML#getEntrezGeneIDsOfComplex(Complex, Map)} It's an
   * advantage to preprocess the entities to exclusively having entities with a
   * name call therefore the method
   * {@link BioPaxL22KGML#getEntitiesWithName(entity)}
   * 
   * The default for a gene id which is not found is -1
   * 
   * @param entity
   * @param xrefs
   * @return
   */
  protected int getEntrezGeneID(entity entity, Map<String, relationshipXref> xrefs) {
    Integer geneID = -1;
    String id = entity.getRDFId();

    // gene id in xrefs?
    if (xrefs.containsKey(id)) {
      geneID = getEntrezGeneIDFromDBxref(xrefs.get(id));
    }

    if (geneID != -1) {
      // gene id found directly
      log.log(Level.FINER, "found: " + entity.getRDFId() + " " + geneID);
      return geneID;
    } else {
      // we have to search the gene id with the gene symbol, adding symbol to
      // the gene symbol set
      String name = entity.getNAME();
      if (name != null) {
        Map<String, entity> geneSymbols = new HashMap<String, entity>();
        
          geneSymbols.put(name, entity);

        if (getEntrezGeneIDOverGeneSymbol(geneSymbols).get(entity) != null)
          geneID = getEntrezGeneIDOverGeneSymbol(geneSymbols).get(entity).intValue();
      }
    }

    return geneID;
  }
  
  /**
   * 
   * @param xref
   * @return Integer of the entered xref or null
   */
  public Integer getEntrezGeneIDFromDBxref(relationshipXref xref) {
    if (xref.getDB().equals("LL")) {
      log.info(xref.getRDFId() + "|" + xref.getID());
      return Integer.parseInt(xref.getID());
    } else {
      return null;
    }
  }  

  /**
   * Maps the entered gene symbol names to a geneID
   * 
   * @param geneSymbolMapper
   *          {@link GeneSymbol2GeneIDMapper}
   * @param nameSet
   *          with names of one {@link BioPAXElement} element
   * @return the gene id o
   */
  protected Map<entity, Integer> getEntrezGeneIDOverGeneSymbol(Map<String, entity> nameSet) {
    log.finest("getGeneIDOverGeneSymbol");
    Map<entity, Integer> geneIDs = new HashMap<entity, Integer>();
    Integer geneID = null;

    for (java.util.Map.Entry<String, entity> symbol : nameSet.entrySet()) {
      try {
        geneID = geneSymbolMapper.map(symbol.getKey());
      } catch (Exception e) {
        log.log(Level.WARNING, "Error while mapping name: " + symbol.getKey() + ".", e);
      }

      if (geneID != null) {
        log.log(Level.FINER, "----- found! Geneid: " + geneID + " " + symbol.getValue().getRDFId()
            + " " + symbol.getKey());
        geneIDs.put(symbol.getValue(), geneID);
        break;
      } else if (symbol.getKey().contains("-")) {
        log.log(Level.FINER,
            "recall for symbol: " + symbol.getValue().getRDFId() + " " + symbol.getKey());
        Map<String, entity> set = new HashMap<String, entity>();
        set.put(symbol.getKey().replace("-", ""), symbol.getValue());
        geneIDs.putAll(getEntrezGeneIDOverGeneSymbol(set));
      } else if (symbol.getKey().contains(" ")) {
        log.log(Level.FINER,
            "recall for symbol: " + symbol.getValue().getRDFId() + " " + symbol.getKey());
        Map<String, entity> set = new HashMap<String, entity>();
        set.put(symbol.getKey().replace(" ", "_"), symbol.getValue());
        geneIDs.putAll(getEntrezGeneIDOverGeneSymbol(set));
      } else {
        log.log(Level.FINER,
            "----- not found " + symbol.getValue().getRDFId() + " " + symbol.getKey());
      }
    }

    return geneIDs;
  }
  
  /**
   * Adds the created {@link Entry} to the
   * {@link BioPaxL22KGML#bc2KeggEntry} map and to the
   * {@link de.zbit.kegg.parser.pathway.Pathway}
   * 
   * @param entity
   * @param keggPW
   * @param mapper
   * @param m
   * @return
   */
  protected Entry createKEGGEntry(entity entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species, EntryType eType, GeneType gType, String graphNameSeparator,
      List<Integer> components) {
    Entry keggEntry;

    // creating new KEGG entry
    String keggname = null;
    Integer geneID = getEntrezGeneID(entity, getMapFromSet(m.getObjects(relationshipXref.class)));
    if (geneID != -1) {
      keggname = mapGeneIDToKEGGID(geneID, species);
    } else {
      keggname = getKEGGUnkownName();
    }

    String graphName = "";
    String names = entity.getNAME();
    if (names != null) {
      
      names = names.trim();
      names = names.replace(" ", "_");

    }
    graphName = names;

    Graphics graphics = null;
    if (eType.equals(EntryType.map))
      graphics = Graphics.createGraphicsForPathwayReference(graphName);
    else if (eType.equals(EntryType.compound))
      graphics = Graphics.createGraphicsForCompound(graphName);
    else
      graphics = new Graphics(graphName);

    keggEntry = new Entry(keggPW, getKeggEntryID(), keggname, eType, graphics);    

    if (components != null) {      
      keggEntry.addComponents(components);
    } 

    

    // checking if entry already exists
    if (!augmentOriginalKEGGpathway){
      Collection<de.zbit.kegg.parser.pathway.Entry> entries = keggPW.getEntries();
      if (entries != null && entries.size() > 0) {
        for (de.zbit.kegg.parser.pathway.Entry entry : entries) {        
            // important to ignore id, because this can differ from file to file
            if (entry.equalsWithoutIDNameReactionComparison(keggEntry)) {            
              keggEntry = entry;
              return keggEntry;
            }        
        }
      }
      // add entry to pathway
      keggPW.addEntry(keggEntry);      
    } else {
      keggEntry=null;
      if (!keggname.startsWith(keggUnknownName)){
        // Search an existing kegg entry, that contains this keggname
        Collection<de.zbit.kegg.parser.pathway.Entry> entries = keggPW.getEntriesForName(keggname);
        keggEntry = de.zbit.kegg.parser.pathway.Pathway.getBestMatchingEntry(keggname, entries);        
      }
    }        
    
    return keggEntry;
  }

  /**
   * parse a BioPax Interaction element
   * 
   * @param entity
   * @param keggPW
   * @param mapper
   * @param m
   * @param species
   */
  protected void parseInteraction(interaction entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species) {
    if (control.class.isAssignableFrom(entity.getClass())) {
      parseControl((control) entity, keggPW, m, species);
    } else if (conversion.class.isAssignableFrom(entity.getClass())) {
      parseConversion((conversion) entity, keggPW, m, species);
    } else {      
      createKEGGEntry(entity, keggPW, m, species);
    }
  }

  /**
   * parse a BioPax Control element
   * 
   * @param entity
   * @param keggPW
   * @param mapper
   * @param m
   * @param species
   */
  private void parseControl(control entity, de.zbit.kegg.parser.pathway.Pathway keggPW, Model m,
      Species species) {
    if (catalysis.class.isAssignableFrom(entity.getClass())) {
      createKEGGReactionRelation(((catalysis) entity).getCONTROLLER(),
          ((catalysis) entity).getCONTROLLED(), getSubtype(((catalysis) entity).getCONTROL_TYPE()),
          keggPW, m, species);
    } else if (modulation.class.isAssignableFrom(entity.getClass())) {
      createKEGGReactionRelation(((modulation) entity).getCONTROLLER(),
          ((modulation) entity).getCONTROLLED(),
          getSubtype(((modulation) entity).getCONTROL_TYPE()), keggPW, m, species);
    } else {
      createKEGGReactionRelation(((control) entity).getCONTROLLER(),
          ((control) entity).getCONTROLLED(), getSubtype(((control) entity).getCONTROL_TYPE()),
          keggPW, m, species);
    } 
  }

  /**
   * Returns the subtypes for a specific ControlType
   * 
   * @param cType
   * @return
   */
  private SubType getSubtype(ControlType cType) {   
    switch (cType) {
      case ACTIVATION:
        return (new SubType(SubType.ACTIVATION));
      case ACTIVATION_ALLOSTERIC:
        return (new SubType(SubType.ACTIVATION));
      case ACTIVATION_NONALLOSTERIC:
        return (new SubType(SubType.ACTIVATION));
      case ACTIVATION_UNKMECH:
        return (new SubType(SubType.ACTIVATION));
      case INHIBITION:
        return (new SubType(SubType.INHIBITION));
      case INHIBITION_ALLOSTERIC:
        return (new SubType(SubType.INHIBITION));
      case INHIBITION_COMPETITIVE:
        return (new SubType(SubType.INHIBITION));
      case INHIBITION_IRREVERSIBLE:
        return (new SubType(SubType.INHIBITION));
      case INHIBITION_NONCOMPETITIVE:
        return (new SubType(SubType.INHIBITION));
      case INHIBITION_OTHER:
        return (new SubType(SubType.INHIBITION));
      case INHIBITION_UNCOMPETITIVE:
        return (new SubType(SubType.INHIBITION));
      case INHIBITION_UNKMECH:
        return (new SubType(SubType.INHIBITION));
      default:
        log.log(Level.SEVERE, "Unkown ControlType: '" + cType.toString() + "'.");
        System.exit(1);
        return null;
    }
  }

  /**
   * For a list of controllers and controlled elements the corresponding KEGG
   * reactions and relations are created
   * 
   * @param controllers
   * @param controlleds
   * @param subtypes
   * @param keggPW
   * @param m
   * @param species
   * @return
   */
  private Entry createKEGGReactionRelation(Set<physicalEntityParticipant> controllers,
      Set<process> controlleds, SubType subtype,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    Entry keggEntry1 = null;
    RelationType relType = null;

    if (controllers.size() >= 1) {
      for (physicalEntityParticipant controller : controllers) {
        keggEntry1 = parseEntity(controller.getPHYSICAL_ENTITY(), keggPW, m, species);
        
        if(keggEntry1.getType().equals(EntryType.map)){
          relType = RelationType.maplink;
        } else {
          relType = RelationType.PPrel;
        }

        if (keggEntry1 != null && controlleds.size() > 0) {
          for (process process : controlleds) {
            if (conversion.class.isAssignableFrom(process.getClass())) {
              conversion con = (conversion) process;
              if (biochemicalReaction.class.isAssignableFrom(con.getClass())
                  || complexAssembly.class.isAssignableFrom(con.getClass())
                  || transportWithBiochemicalReaction.class.isAssignableFrom(con.getClass())) {
                if (!augmentOriginalKEGGpathway) {
                  Reaction r = null;
                  try {
                    r = createKEGGReaction(((biochemicalReaction) con).getLEFT(),
                        ((biochemicalReaction) con).getRIGHT(), keggPW, m, species);

                  } catch (ClassCastException e) {
                    try {
                      r = createKEGGReaction(((complexAssembly) con).getLEFT(),
                          ((complexAssembly) con).getRIGHT(), keggPW, m, species);
                    } catch (ClassCastException e2) {
                      r = createKEGGReaction(((transportWithBiochemicalReaction) con).getLEFT(),
                          ((transportWithBiochemicalReaction) con).getRIGHT(), keggPW, m, species);
                    }
                  }

                  if (relType.equals(RelationType.maplink) && r!=null) {
                    for (ReactionComponent rc : r.getSubstrates()) {
                      createKEGGRelation(keggPW, keggEntry1.getId(), rc.getId(), relType, subtype);
                    }
                  } else if (relType.equals(RelationType.PPrel) && r!=null) {
                    keggEntry1.appendReaction(r.getName());
                  }  
                }
              } else if (transport.class.isAssignableFrom(con.getClass())) {
                List<Relation> rels = createKEGGRelations(((transport) con).getLEFT(),
                    ((transport) con).getRIGHT(), keggPW, m, species, RelationType.PPrel, 
                    new SubType(SubType.STATE_CHANGE));                
                for (Relation rel : rels) {
                  if (rel !=null)
                    createKEGGRelation(keggPW, keggEntry1.getId(), rel.getEntry2(), relType, subtype);
                }
              } else {
                log.severe("Not programmed case: controlled interface '" + con.getModelInterface()
                    + "'");
                System.exit(1);
              }
            } else if (pathway.class.isAssignableFrom(process.getClass())) {
              Entry keggEntry2 = createKEGGEntry((pathway) process, keggPW, m, species,
                  EntryType.map, null, ",", null);
              if (keggEntry2 !=null)
                createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(),
                  relType, subtype);
            } else if (interaction.class.isAssignableFrom(process.getClass())) {
              interaction inter = ((interaction)process); 
              if (!inter.getNAME().isEmpty()){                
                Entry keggEntry2 = createKEGGEntry(inter, keggPW, m, species);
                if (keggEntry2 !=null)
                  createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(),
                    relType, subtype);
              }
            } else {
              log.severe("Process: " + process.getModelInterface() + "-This should not happen!");
              System.exit(1);
            }
            // ControlType (0 or 1) - up to now ignored

            // Cofactor = PhysicalEntity (0 or more) - up to now ignored

            // CatalysisDirection - up to now ignored
          }
        } else {
          break;
        }
      }
    }

    return keggEntry1;
  }

  /**
   * converts an interaction if it could not be mapped to another subclass like control or 
   * conversion
   * 
   * @param inter
   * @param keggPW
   * @param m
   * @param species
   * @return
   */
  private Entry createKEGGEntry(interaction inter, Pathway keggPW, Model m, Species species) {
    Entry keggEntry1 = null;
    List<InteractionParticipant> participants = Utils.getListOfCollection(inter.getPARTICIPANTS());
    if (participants.size() > 1) {
      for (int i = 0; i < participants.size(); i++) {
        if (pathway.class.isAssignableFrom(participants.get(i).getClass())) {
          keggEntry1 = parseEntity(((pathway) participants.get(i)), keggPW, m, species);
        } else if (physicalEntityParticipant.class.isAssignableFrom(participants.get(i).getClass())) {
          keggEntry1 = parseEntity(
              ((physicalEntityParticipant) participants.get(i)).getPHYSICAL_ENTITY(), keggPW, m,
              species);
        } else {
          log.log(Level.SEVERE, "1 This should not happen: '"
              + participants.get(i).getModelInterface() + "'.");
          System.exit(1);
        }
        for (int j = 1; j < participants.size(); j++) {
          Entry keggEntry2 = null;
          if (pathway.class.isAssignableFrom(participants.get(j).getClass())) {
            keggEntry2 = parseEntity(((pathway) participants.get(j)), keggPW, m, species);
          } else if (physicalEntityParticipant.class.isAssignableFrom(participants.get(j)
              .getClass())) {
            keggEntry2 = parseEntity(
                ((physicalEntityParticipant) participants.get(j)).getPHYSICAL_ENTITY(), keggPW, m,
                species);
          } else {
            log.log(Level.SEVERE, "2 This should not happen: '"
                + participants.get(j).getModelInterface() + "'.");
            System.exit(1);
          }

          createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), RelationType.maplink,
              null);
        }

      }
    } else if (participants.size() > 0) {
      if (pathway.class.isAssignableFrom(participants.get(0).getClass())) {
        keggEntry1 = parseEntity(((pathway) participants.get(0)), keggPW, m, species);
      } else if (physicalEntityParticipant.class.isAssignableFrom(participants.get(0).getClass())) {
        keggEntry1 = parseEntity(
            ((physicalEntityParticipant) participants.get(0)).getPHYSICAL_ENTITY(), keggPW, m,
            species);
      } else {
        log.log(Level.SEVERE, "3 - This should not happen: '"
            + participants.get(0).getModelInterface() + "'.");
        System.exit(1);
      }
    } else {
      // creating new KEGG entry
      String keggname = getKEGGUnkownName();

      String graphName = "";
      String names = inter.getNAME();
      if (names != null) {

        names = names.trim();
        names = names.replace(" ", "_");

      }
      graphName = names;

      Graphics graphics = Graphics.createGraphicsForPathwayReference(graphName);
      ;
      EntryType eType = EntryType.map;

      keggEntry1 = new Entry(keggPW, getKeggEntryID(), keggname, eType, graphics);

      // checking if entry already exists
      if (!augmentOriginalKEGGpathway) {
        Collection<de.zbit.kegg.parser.pathway.Entry> entries = keggPW.getEntries();
        if (entries != null && entries.size() > 0) {
          for (de.zbit.kegg.parser.pathway.Entry entry : entries) {
            // important to ignore id, because this can differ from file to file
            if (entry.equalsWithoutIDNameReactionComparison(keggEntry1)) {
              keggEntry1 = entry;
              return keggEntry1;
            }
          }
        }
        // add entry to pathway
        keggPW.addEntry(keggEntry1);
      } else {
        keggEntry1 = null;
        if (!keggname.startsWith(keggUnknownName)) {
          // Search an existing kegg entry, that contains this keggname
          Collection<de.zbit.kegg.parser.pathway.Entry> entries = keggPW
              .getEntriesForName(keggname);
          keggEntry1 = de.zbit.kegg.parser.pathway.Pathway.getBestMatchingEntry(keggname, entries);
        }
      }
    }

    return keggEntry1;
  }


  /**
   * parse a BioPax Conversion element
   * 
   * @param entity
   * @param keggPW
   * @param species
   * @param m
   */
  private void parseConversion(interaction entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species) {
    if (complexAssembly.class.isAssignableFrom(entity.getClass())) {
      if  (!augmentOriginalKEGGpathway)
        createKEGGReaction(((complexAssembly) entity).getLEFT(),
          ((complexAssembly) entity).getRIGHT(), keggPW, m, species);
    } else if (biochemicalReaction.class.isAssignableFrom(entity.getClass())) {
      if  (!augmentOriginalKEGGpathway)
        createKEGGReaction(((biochemicalReaction) entity).getLEFT(),
          ((biochemicalReaction) entity).getRIGHT(), keggPW, m, species);
    } else if (transport.class.isAssignableFrom(entity.getClass())) {
      createKEGGRelations(((transport) entity).getLEFT(), ((transport) entity).getRIGHT(), keggPW,
          m, species, RelationType.PPrel, new SubType(SubType.STATE_CHANGE));
    } else if (transportWithBiochemicalReaction.class.isAssignableFrom(entity.getClass())) {
      if  (!augmentOriginalKEGGpathway)
        // BiochemicalReaction br = (TransportWithBiochemicalReaction) entity;
        // deltaG, deltaH, deltaS, ec, and KEQ are ignored
        createKEGGReaction(((biochemicalReaction) entity).getLEFT(),
          ((biochemicalReaction) entity).getRIGHT(), keggPW, m, species);
    } else {
      log.warning("Unknown kind of Conversion: " + entity.getModelInterface() + "-"
          + entity.getRDFId());
    }
  }

  /**
   * This method is called to create for two PhysicalEntitiy sets relations
   * 
   * @param set
   * @param set2
   * @param keggPW
   * @param m
   * @param species
   * @param type
   * @return
   */
  private List<Relation> createKEGGRelations(Set<physicalEntityParticipant> set, Set<physicalEntityParticipant> set2,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species, RelationType type, 
      SubType subType) {

    List<Relation> relations = new ArrayList<Relation>();

    for (physicalEntityParticipant left : set) {     
      Entry keggEntry1 = parseEntity(left.getPHYSICAL_ENTITY(), keggPW, m, species);
      
      if (keggEntry1 !=null){
        for (physicalEntityParticipant right : set2) {
          Entry keggEntry2 = parsePhysicalEntity(right.getPHYSICAL_ENTITY(), keggPW, m, species);
          if (keggEntry1 !=null){
            Relation r = createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), type, subType);
            if (!relations.contains(r))
              relations.add(r);  
          }          
        }  
      }    
    }

    return relations;
  }

  /**
   * This method first checks if a relation already exist for the given
   * variables, if yes the method returns the existing relation, otherwise it
   * creates a new KEGG relation and adds it to the pathway
   * 
   * @param keggPW
   * @param keggEntry1Id
   * @param keggEntry2Id
   * @param type
   * @param subTypes
   * @return
   */
  private Relation createKEGGRelation(de.zbit.kegg.parser.pathway.Pathway keggPW, int keggEntry1Id,
      int keggEntry2Id, RelationType type, SubType subType) {
    ArrayList<Relation> existingRels = keggPW.getRelations();
    Relation r = null;

    // Check if it already exists and only create novel relations.
    if (existingRels.size() > 0) {
      for (Relation rel : existingRels) {
        boolean relExists = true;
        if ((rel.getEntry1() == keggEntry1Id && rel.getEntry2() == keggEntry2Id)) {
          relExists &= rel.isSetType() == (type != null);
          if (relExists && type != null)
            relExists &= (rel.getType().equals(type));
         
          if (relExists && subType !=null) {
            r = rel;
            boolean added = r.addSubtype(subType);
            if (augmentOriginalKEGGpathway && added)
              addedSubTypes++;
            return r;
          }          
        }
      }

      r = new Relation(keggEntry1Id, keggEntry2Id, type, subType);
    } else {
      r = new Relation(keggEntry1Id, keggEntry2Id, type, subType);
    }
    // If we are here, the relation r is NOVEL AND NOT CURRENTLY IN THE PATHWAY

    // Add the relation to the pathway
    
    
    if (augmentOriginalKEGGpathway){
      if(keggPW.getEntryForId(keggEntry1Id)!=null && keggPW.getEntryForId(keggEntry2Id)!=null) {
        // Only add relations if nodes for the relation are present.
        if (keggEntry1Id != keggEntry2Id){
          keggPW.addRelation(r);
          newAddedRelations++; 
        }        
        
      }
    } else {
      keggPW.addRelation(r);
    }

    return r;
  }

  /**
   * Checks if the reaction already exists in the kegg pathway. If yes the
   * existing relaction is returned, otherwise a new reaction is created and
   * added to the pathway
   * 
   * @param entity
   * @param list
   * @param keggPW
   * @param m
   */
  private Reaction createKEGGReaction(Set<physicalEntityParticipant> lefts, Set<physicalEntityParticipant> rights,
      Pathway keggPW, Model m, Species species) {
    List<ReactionComponent> products = new ArrayList<ReactionComponent>();
    List<ReactionComponent> substrates = new ArrayList<ReactionComponent>();

    for (physicalEntityParticipant left : lefts) {
      Entry keggEntry =  parsePhysicalEntity(left.getPHYSICAL_ENTITY(), keggPW, m, species);
         
      if (keggEntry != null) {
        ReactionComponent rc = new ReactionComponent(keggEntry.getId(), keggEntry.getName());
        products.add(rc);
  
      }
    }
    
    for (physicalEntityParticipant right : rights) {
      Entry keggEntry = parsePhysicalEntity(right.getPHYSICAL_ENTITY(), keggPW, m, species);
      if (keggEntry != null) {
        ReactionComponent rc = new ReactionComponent(keggEntry.getId(), keggEntry.getName());        
        substrates.add(rc);
  
      }
    }

    Reaction r = null;
    boolean reactionExists = false;
    for (Reaction existingReact : keggPW.getReactions()) {
      List<ReactionComponent> existingProds = existingReact.getProducts();
      List<ReactionComponent> extistingSubs = existingReact.getSubstrates();

      if (existingReact.getType().equals(species) && existingProds.size() == products.size()
          && extistingSubs.size() == substrates.size()) {
        boolean allReactantsIn = true;

        for (ReactionComponent prod : products) {
          if (!existingProds.contains(prod)) {
            allReactantsIn = false;
            break;
          }
        }

        if (allReactantsIn) {
          for (ReactionComponent sub : substrates) {
            if (!extistingSubs.contains(sub)) {
              allReactantsIn = false;
              break;
            }
          }
        }

        if (allReactantsIn) {
          reactionExists = true;
          r = existingReact;
          break;
        }
      }
    }

    if (!reactionExists) {
      r = new Reaction(keggPW, getReactionName(), ReactionType.other);
      r.addProducts(products);
      r.addSubstrates(substrates);
      keggPW.addReaction(r);
    }

    return r;
  }
}
