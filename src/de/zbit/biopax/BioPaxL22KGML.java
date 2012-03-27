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
import org.biopax.paxtools.model.level2.dataSource;
import org.biopax.paxtools.model.level2.entity;
import org.biopax.paxtools.model.level2.interaction;
import org.biopax.paxtools.model.level2.modulation;
import org.biopax.paxtools.model.level2.openControlledVocabulary;
import org.biopax.paxtools.model.level2.pathway;
import org.biopax.paxtools.model.level2.pathwayComponent;
import org.biopax.paxtools.model.level2.pathwayStep;
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
import org.biopax.paxtools.model.level2.unificationXref;
import org.biopax.paxtools.model.level2.xref;

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
import de.zbit.kegg.parser.pathway.ext.EntryExtended;
import de.zbit.kegg.parser.pathway.ext.EntryTypeExtended;
import de.zbit.util.DatabaseIdentifiers;
import de.zbit.util.DatabaseIdentifiers.IdentifierDatabases;
import de.zbit.util.SortedArrayList;
import de.zbit.util.Species;
import de.zbit.util.Utils;
import de.zbit.util.progressbar.ProgressBar;

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
  public void createKGMLForBioPaxFile(Model m, String comment, String fileName, boolean writeEntryExtended) {
    String folder = createDefaultFolder(m.getLevel());
    createKGMLForBioPaxFile(m, comment, fileName, folder, writeEntryExtended);
  }
  /**
   * The methods parse a BioPax file which contains no <bp>Pathway: ....</bp> tag
   * 
   * @param m
   * @param pathwayName used as pathway name and title 
   * @param folder where the KGML is saved
   */
  public void createKGMLForBioPaxFile(Model m, String comment, String pathwayName, String folder, boolean writeEntryExtended) {
    // determine the organism
    Species species = new Species("Homo sapiens", "_HUMAN", "human", "hsa", 9606);
    
    Set<bioSource> orgs = m.getObjects(bioSource.class);
    if (orgs != null && orgs.size() > 0) {
      bioSource org = orgs.iterator().next();
      species = determineSpecies(org);      
    } else {
      log.info("No pathway species could be determined, default species '"
          + species.getCommonName() + "' is used.");
    }
    
    initalizeMappers(species);    

    int pathNo = determineKEGGPathwayNumber(pathwayName);
    de.zbit.kegg.parser.pathway.Pathway keggPW = new de.zbit.kegg.parser.pathway.Pathway(
        pathwayName, species.getKeggAbbr(), pathNo, pathwayName);
    keggPW.setComment(comment);
    
    log.info("Converting pathway '" + pathwayName + "'.");
    
    for (entity entity : m.getObjects(entity.class)) {
      parseEntity(entity, keggPW, m, species);
    }

    String fileName = folder + KGMLWriter.createFileName(keggPW);
    KGMLWriter.writeKGML(keggPW, fileName, writeEntryExtended);
  } 

  
  /**
   * This method creates for each pathway in the model a KGML file with the
   * pathway name and saves the pathways in an default folder see 
   * {@link BioPax2KGML#createDefaultFolder(org.biopax.paxtools.model.BioPAXLevel)}
   * 
   * @param m
   */
  public void createKGMLsForPathways(Model m, String comment, Set<pathway> pathways, boolean writeEntryExtended) {
    String folder = createDefaultFolder(m.getLevel());
    createKGMLsForPathways(m, comment, folder, pathways, writeEntryExtended);
  }
  
  /**
   * This method creates for each pathway in the model a KGML file with the
   * pathway name
   * 
   * @param m
   */
  public void createKGMLsForPathways(Model m, String comment, String folder, Set<pathway> pathways, boolean writeEntryExtended) {
    log.info("Creating KGML files for pathways.");
    Collection<Pathway> keggPWs = parsePathways(m, comment, pathways);
    
    for (Pathway keggPW : keggPWs) {
      String fileName = folder + KGMLWriter.createFileName(keggPW);
      KGMLWriter.writeKGML(keggPW, fileName, writeEntryExtended);
    }    
  }

  /**
   * this method parses the biopax pathway by firstly determining the pathway species and
   * then parsing the single pathway
   * @param m
   * @param pathways
   * @return
   */
  protected Collection<Pathway> parsePathways(Model m, String comment, Set<pathway> pathways) {    
    Species oldSpecies = new Species("Homo sapiens", "_HUMAN", "human", "hsa", 9606);
    initalizeMappers(oldSpecies);
    
    Collection<Pathway> keggPWs = new ArrayList<Pathway>();

    for (pathway pathway : pathways) {
      // determine the pathway organism - it's done here to save time, while initializing the mappers
      Species newSpecies = determineSpecies(pathway.getORGANISM());
      if(!newSpecies.equals(oldSpecies)){
        initalizeMappers(newSpecies);
        oldSpecies = newSpecies;
      }
      keggPWs.add(parsePathway(m, comment, pathway, oldSpecies));
    }
    return keggPWs;
  }
  
  /**
   * determines the species of the pathway and returns {@link Species}
   * the default species it "Homo sapiens"
   * @param pathway
   * @return
   */
  protected Species determineSpecies(bioSource pwOrg) {
    Species species = new Species("Homo sapiens", "_HUMAN", "human", "hsa", 9606), 
    detSpecies = null;

    if (pwOrg != null) {
      unificationXref ref = pwOrg.getTAXON_XREF();
      if (ref != null) {
        if (ref.getDB().toLowerCase().equals(DatabaseIdentifiers.IdentifierDatabases.NCBI_Taxonomy.toString().toLowerCase())) {
          detSpecies = Species.search(allSpecies, ref.getID(), Species.NCBI_TAX_ID);
        }
      }
      if (pwOrg.getNAME() != null) {
        String newSpecies = pwOrg.getNAME();
        detSpecies = Species.search(allSpecies, newSpecies, Species.COMMON_NAME);        
      } 
    } 

    if(detSpecies != null){
      log.info("Determined pathway species '" + species.getCommonName() + "'.");
      return detSpecies;
    }
    
    log.info("No pathway species could be determined, default species '"
        + species.getCommonName() + "' is used.");
    return species;
  }
  
  /**
   * parses the biopax pathway
   * @param m
   * @param pathway
   * @param species
   * @return
   */
  Pathway parsePathway(Model m, String comment, pathway pathway, Species species) {
    if(geneSymbolMapper==null && geneIDKEGGmapper==null){
      initalizeMappers(species);
    }
    
    // create the pathway
    int number = getKeggPathwayNumber(pathway.getRDFId());
    String sourceDB = getSourceDB(pathway.getDATA_SOURCE());
    String pwName = pathway.getNAME();
    de.zbit.kegg.parser.pathway.Pathway keggPW = new de.zbit.kegg.parser.pathway.Pathway(
        sourceDB + String.valueOf(number), species.getKeggAbbr(), number,
        pwName);
    keggPW.setComment(comment);
    
    log.info("Converting pathway '" + pwName + "'.");
    
    if(!sourceDB.isEmpty())
      keggPW.setLink(sourceDB);

    //TODO: it is not possible to define which image link to set, perhaps using datasource..., 
    // but too much databases to be conform for each
//    addImageLinkToKEGGpathway(species, pathway.getNAME(), keggPW);
    

    for (pathwayComponent interaction : pathway.getPATHWAY_COMPONENTS()) {
      if (pathwayStep.class.isAssignableFrom(interaction.getClass())){
        parsePathwayStep((pathwayStep)interaction, keggPW, m, species);
      } else if (interaction.class.isAssignableFrom(interaction.getClass())){
        parseInteraction((interaction) interaction, keggPW, m, species);  
      } else {
        log.log(Level.SEVERE, "Could not parse: '" + interaction.getModelInterface() + "'.");
      }      
    }

    return keggPW;
  }
  
  /**
   * sets the source of the data if available to the class
   * @param sources
   * @param keggPW
   */
  private String getSourceDB(Set<dataSource> sources) {
    String sourceDB = "";
    if (sources!=null && sources.size()>0){
      for (dataSource source : sources) {
        if(source.getCOMMENT()!=null && source.getCOMMENT().size()>0){
          sourceDB = source.getCOMMENT().iterator().next();
        }      
      }
    }
    
    return sourceDB;
  }
  /**
   * 
   * @param interaction
   * @param keggPW
   * @param m
   * @param species
   */
  private void parsePathwayStep(pathwayStep pwStep, Pathway keggPW, Model m, Species species) {
    Set<process> interactions = pwStep.getSTEP_INTERACTIONS();
    for (process process : interactions) {
      parseInteraction((interaction)process, keggPW, m, species);
    }
    Set<pathwayStep> nextSteps = pwStep.getNEXT_STEP();
    for (pathwayStep pathwayStep : nextSteps) {
      parsePathwayStep(pathwayStep, keggPW, m, species);
    }
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

//      Map<entity, Integer> geneIDs = getEntrezGeneIDs(pwEntities, species, xrefs);
//      for (java.util.Map.Entry<entity, Integer> entity : geneIDs.entrySet()) {
//        log.log(Level.FINER, "----res: " + entity.getKey() + " " + entity.getValue());
//        pw.addGeneID(entity.getValue());
//      }
      //TODO rewrite this method
      // }//TODO: is necessary to uncomment!!!!
      bar.DisplayBar();
    }

    return pathways;
  }
  
  /**
   * returns a list of all pathways containing pathway components
   * @param m
   * @return
   */
  protected List<String> getListOfPathways(Model m){
    List<String> pws = new SortedArrayList<String>();
    
    Set<pathway> list = m.getObjects(pathway.class);
    for (pathway pw : list) {
      if (pw.getPATHWAY_COMPONENTS().size()>0)
        pws.add(pw.getNAME());
    }
    
    return pws;
  }
  
  /**
   * 
   * @param m
   * @param name
   * @return the BioPaxPathway with the specific name
   */
  protected pathway getPathwayByName(Model m, String name){
    pathway pw = null;
    Set<pathway> list = m.getObjects(pathway.class);
    
    for (pathway p : list) {
      if(p.getNAME().equals(name))
        return p;
    }
    
    return pw;
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
  private EntryExtended parseEntity(entity entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species) {
    EntryExtended keggEntry = null;
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
  private EntryExtended parsePhysicalEntity(physicalEntity entity,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    EntryExtended keggEntry = null;

    if (complex.class.isAssignableFrom(entity.getClass())) {
      List<Integer> components = createComplexComponentList(((complex) entity).getCOMPONENTS(),
          keggPW, m, species);
      keggEntry = createKEGGEntry((entity) entity, keggPW, m, species, EntryType.group, null, "/",
          components);
    } else if (sequenceEntity.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, EntryTypeExtended.dna, ",",
          null);
    } else if (protein.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.enzyme, EntryTypeExtended.protein,
          ",", null);
    } else if (rna.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, EntryTypeExtended.rna, ",",
          null);
    } else if (smallMolecule.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.compound, EntryTypeExtended.unknown,
          ",", null);
    } else {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.other, EntryTypeExtended.unknown,
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
      EntryExtended keggEntry = parsePhysicalEntity(physicalEntity.getPHYSICAL_ENTITY(), keggPW, m, species);
      if (keggEntry != null)
        components.add(keggEntry.getId());
    }

    return components;
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
  protected EntryExtended createKEGGEntry(entity entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species, EntryType eType, EntryTypeExtended gType, String graphNameSeparator,
      List<Integer> components) {
    EntryExtended keggEntry;

    // get all availabe database identifiers of the entity
    Map<IdentifierDatabases, Collection<String>> identifiers = 
      getDatabaseIdentifiers(entity, eType, gType);
   
    // determine graph name and gene symbols   
    String graphName = "";
    String names = entity.getNAME();
    if (names != null) {      
      names = names.trim();
      names = names.replace(" ", "_");
    }
    graphName = names;
        
    // determine gene id
    Integer geneID = -1;
    Collection<String> geneIDs = identifiers.get(IdentifierDatabases.EntrezGene);   
    if (geneIDs!=null && geneIDs.size()>0){
      geneID = Integer.parseInt((Utils.collectionToList(geneIDs).get(0)));  
    } else {
      // we have to search the gene id with the gene symbol, adding symbol to
      // the gene symbol set  
      Collection<String> geneSymbols = identifiers.get(IdentifierDatabases.GeneSymbol);
      if(geneSymbols!=null && geneSymbols.size()>0)
        geneID = getEntrezGeneIDForGeneSymbol(geneSymbols).intValue();   
    }
    
    // determine kegg name
    String keggname = null;
    if (geneID != -1) {
      keggname = mapGeneIDToKEGGID(geneID, species);
    } else {
      keggname = getKEGGUnkownName();
    }

    // create graphics
    Graphics graphics = null;
    if (eType.equals(EntryType.map))
      graphics = Graphics.createGraphicsForPathwayReference(graphName);
    else if (eType.equals(EntryType.compound))
      graphics = Graphics.createGraphicsForCompound(graphName);
    else
      graphics = new Graphics(graphName);

    keggEntry = new EntryExtended(keggPW, getKeggEntryID(), keggname, eType, gType, graphics);    
    
    
    // set further information to the entry
    keggEntry.addDatabaseIdentifiers(identifiers);   
    
    
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
              keggEntry = (EntryExtended) entry;
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
        keggEntry = (EntryExtended) de.zbit.kegg.parser.pathway.Pathway.getBestMatchingEntry(keggname, entries);        
      }
    }        
    
    return keggEntry;
  }

  /**
   * add if available further entity information
   * 
   * @param keggEntry
   * @param entity
   */
  private Map<IdentifierDatabases, Collection<String>> getDatabaseIdentifiers(entity entity,
      EntryType eType, EntryTypeExtended gType) {    
    Map<IdentifierDatabases, Collection<String>> map = 
      new HashMap<DatabaseIdentifiers.IdentifierDatabases, Collection<String>>();
    // data source
    Set<dataSource> ds = entity.getDATA_SOURCE();
    Set<xref> xrefs = entity.getXREF();
    if (ds.size() != 0) {
      for (dataSource d : ds) {
        if (d.getNAME().size()>0){
          String db = Utils.collectionToList(d.getNAME()).get(0);
          IdentifierDatabases dbIdentifier = DatabaseIdentifiers.getIdentifier(db);
          if (dbIdentifier != null){
            if (Utils.collectionToList(d.getCOMMENT())!=null &&
                Utils.collectionToList(d.getCOMMENT()).size()>0) {
              String comment = Utils.collectionToList(d.getCOMMENT()).get(0);
              if(!comment.isEmpty())
                Utils.addToMapOfSets(map, dbIdentifier, comment);  
            }            
          }
        }
      }
    } 
    
    // xrefs
    if (xrefs.size() != 0) {
      for (xref d : xrefs) {
        if (!d.getDB().isEmpty()){          
          IdentifierDatabases dbIdentifier = DatabaseIdentifiers.getIdentifier(d.getDB());
          if (dbIdentifier != null) {
            if (Utils.collectionToList(d.getCOMMENT())!=null &&
                Utils.collectionToList(d.getCOMMENT()).size()>0) {
              String comment = Utils.collectionToList(d.getCOMMENT()).get(0);
              if(!comment.isEmpty())
                Utils.addToMapOfSets(map, dbIdentifier, comment);
            }
          } else if (d.getDB().equalsIgnoreCase("LL")) { // special case in PID database files
            if (Utils.collectionToList(d.getCOMMENT())!=null &&
                Utils.collectionToList(d.getCOMMENT()).size()>0) {
              String comment = Utils.collectionToList(d.getCOMMENT()).get(0);
              if (!comment.isEmpty())
                Utils.addToMapOfSets(map, IdentifierDatabases.EntrezGene, comment);              
            }
          }
        }
      }
    }
    
    // gene symbols are assigend depending on the EntryTypes
    // EntryTypeExtended: protein, dna_region, rna_region, dna, rna, unknown;
    // EntryType:         ortholog, enzyme, reaction, gene, group, compound, map, genes, other
    if(!eType.equals(EntryType.map)){
      String names = entity.getNAME();
      if (names != null) {      
        names = names.trim();
        Utils.addToMapOfSets(map, IdentifierDatabases.GeneSymbol, names);
      }  
    }
  
    
    return map;
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
    if (cType!=null){
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
    
    return null;    
  }
  
  private SubType getSubtype(Set<openControlledVocabulary> iTypes) {
    for (openControlledVocabulary iType : iTypes) {
      String type = iType.getRDFId(); 
      if (type!=null)
        if(type.toLowerCase().endsWith("activation")){
          return new SubType(SubType.ACTIVATION);
        } else if(type.toLowerCase().endsWith("inhibition")){
          return new SubType(SubType.INHIBITION);
        } else if(type.toLowerCase().endsWith("transcription")){
          return new SubType(SubType.EXPRESSION);
        } else if(type.toLowerCase().endsWith("translation")){
          return new SubType(SubType.EXPRESSION);
        } else if(type.toLowerCase().endsWith("molecular_interaction")){
          return new SubType(SubType.BINDING);
        } else if(type.toLowerCase().endsWith("hedgehog_cleavage_and_lipidation")){
          return new SubType(SubType.INDIRECT_EFFECT);
        } else {
          log.info("--- Type --- " + type);
          return new SubType(SubType.STATE_CHANGE);
        }      
    }
    return null;
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
  private EntryExtended createKEGGReactionRelation(Set<physicalEntityParticipant> controllers,
      Set<process> controlleds, SubType subtype,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    EntryExtended keggEntry1 = null;
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
              } else if (conversion.class.isAssignableFrom(con.getClass())){
                  List<Relation> rels =  createKEGGRelations(con.getLEFT(), con.getRIGHT(), keggPW, m, species, 
                    RelationType.PPrel, getSubtype(con.getINTERACTION_TYPE()));
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
              EntryExtended keggEntry2 = createKEGGEntry((pathway) process, keggPW, m, species,
                  EntryType.map, null, ",", null);
              if (keggEntry2 !=null)
                createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(),
                  relType, subtype);
            } else if (interaction.class.isAssignableFrom(process.getClass())) {
              interaction inter = ((interaction)process); 
              if (!inter.getNAME().isEmpty()){                
                EntryExtended keggEntry2 = createKEGGEntry(inter, keggPW, m, species);
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
  private EntryExtended createKEGGEntry(interaction inter, Pathway keggPW, Model m, Species species) {
    EntryExtended keggEntry1 = null;
    List<InteractionParticipant> participants = Utils.iterableToList(inter.getPARTICIPANTS());
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
          EntryExtended keggEntry2 = null;
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

          createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), RelationType.other,
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

      keggEntry1 = new EntryExtended(keggPW, getKeggEntryID(), keggname, eType, graphics);

      // checking if entry already exists
      if (!augmentOriginalKEGGpathway) {
        Collection<de.zbit.kegg.parser.pathway.Entry> entries = keggPW.getEntries();
        if (entries != null && entries.size() > 0) {
          for (de.zbit.kegg.parser.pathway.Entry entry : entries) {
            // important to ignore id, because this can differ from file to file
            if (entry.equalsWithoutIDNameReactionComparison(keggEntry1)) {
              keggEntry1 = (EntryExtended) entry;
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
          keggEntry1 = (EntryExtended) de.zbit.kegg.parser.pathway.Pathway.getBestMatchingEntry(keggname, entries);
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
    } else if (conversion.class.isAssignableFrom(entity.getClass())){
      createKEGGRelations(((conversion)entity).getLEFT(), ((conversion)entity).getRIGHT(), keggPW, 
          m, species, RelationType.PPrel, getSubtype(((conversion)entity).getINTERACTION_TYPE()));
    } else {
      log.log(Level.SEVERE, "Unknown kind of Conversion: " + entity.getModelInterface());
      System.exit(1);
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
      EntryExtended keggEntry1 = parseEntity(left.getPHYSICAL_ENTITY(), keggPW, m, species);
      
      if (keggEntry1 !=null){
        for (physicalEntityParticipant right : set2) {
          EntryExtended keggEntry2 = parsePhysicalEntity(right.getPHYSICAL_ENTITY(), keggPW, m, species);
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
      EntryExtended keggEntry =  parsePhysicalEntity(left.getPHYSICAL_ENTITY(), keggPW, m, species);
         
      if (keggEntry != null) {
        ReactionComponent rc = new ReactionComponent(keggEntry.getId(), keggEntry.getName());
        products.add(rc);
  
      }
    }
    
    for (physicalEntityParticipant right : rights) {
      EntryExtended keggEntry = parsePhysicalEntity(right.getPHYSICAL_ENTITY(), keggPW, m, species);
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

      if (existingProds.size() == products.size()
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
