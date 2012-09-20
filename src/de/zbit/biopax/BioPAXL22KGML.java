/*
 * $Id$
 * $URL$
 * ---------------------------------------------------------------------
 * This file is part of SBVC, the systems biology visualizer and
 * converter. This tools is able to read a plethora of systems biology
 * file formats and convert them to an internal data structure.
 * These files can then be visualized, either using a simple graph
 * (KEGG-style) or using the SBGN-PD layout and rendering constraints.
 * Some currently supported IO formats are SBML (+qual, +layout), KGML,
 * BioPAX, SBGN, etc. Please visit the project homepage at
 * <http://www.cogsys.cs.uni-tuebingen.de/software/SBVC> to obtain the
 * latest version of SBVC.
 *
 * Copyright (C) 2012 by the University of Tuebingen, Germany.
 *
 * SBVC is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation. A copy of the license
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
import java.util.Iterator;
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
import org.biopax.paxtools.model.level3.BioSource;

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
import de.zbit.util.DatabaseIdentifierTools;
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
public class BioPAXL22KGML extends BioPAX2KGML {
  
  public static final Logger log = Logger.getLogger(BioPAXL22KGML.class.getName());
 
  /**
   * The methods parse a BioPax file which contains no <bp>Pathway: ....</bp> tag
   * 
   * @param m
   * @param pathwayName used as pathway name and title 
   * @param folder where the KGML is saved
   */
  public de.zbit.kegg.parser.pathway.Pathway createPathwayFromBioPaxFile(Model m, String comment, 
      String pathwayName, Species species) {
    // determine the organism
    if (species != null)
      initalizeMappers(species);

    
    de.zbit.kegg.parser.pathway.Pathway keggPW = createPathwayInstance(comment, m, species, pathwayName, null);
    log.info("Converting pathway '" + keggPW.getTitle() + "'.");
    
    for (entity entity : m.getObjects(entity.class)) {
      parseEntity(entity, keggPW, m, species, null);
    }

    return keggPW;
  }

  
  /**
   * determines the pathway species
   * @param m
   * @return
   */
  public static Collection<Species> getSpecies(Model m) {
    Collection<Species> species = new HashSet<Species>();
    Set<bioSource> orgs = m.getObjects(bioSource.class);
    if (orgs != null && orgs.size() > 0) {
      Iterator<bioSource> it = orgs.iterator();
      while (it.hasNext()){
        bioSource org = it.next();
        Species s = determineSpecies(org);
        if (s!=null)
          species.add(s);  
      }
      
    } else {
      log.info("No specific pathway species could be determined.");
    }
    return species;
  } 
  
  /**
   * this method parses the biopax pathway by firstly determining the pathway species and
   * then parsing the single pathway
   * @param m
   * @param pathways
   * @return
   */
  protected Collection<Pathway> createPathways(Model m, String comment, Set<pathway> pathways, 
      Species species) {    
    if (species!=null) {
      initalizeMappers(species);
    }
    
    Collection<Pathway> keggPWs = new ArrayList<Pathway>();

    for (pathway pathway : pathways) {
      // determine the pathway organism - it's done here to save time, while initializing the mappers
      Species newSpecies = determineSpecies(pathway.getORGANISM());
      if(newSpecies!=null && !newSpecies.equals(species)){
        initalizeMappers(newSpecies);
        species = newSpecies;
      }
      keggPWs.add(createPathway(m, comment, pathway, species));
    }
    return keggPWs;
  }
  
  /**
   * determines the species of the pathway and returns {@link Species}
   * the default species is null
   * @param pathway
   * @return
   */
  protected static Species determineSpecies(bioSource pwOrg) {
    Species detSpecies = null;

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
      log.info("Determined pathway species '" + detSpecies.getCommonName() + "'.");
      return detSpecies;
    }
    
    log.info("No pathway species could be determined.");
    return detSpecies;
  }
  
  /**
   * parses the biopax pathway
   * @param m
   * @param pathway
   * @param species
   * @return
   */
  Pathway createPathway(Model m, String comment, pathway pathway, Species species) {
    if(species!=null && geneSymbolMapper==null && geneIDKEGGmapper==null){
      initalizeMappers(species);
    }
    
    // get the pathway references
    Map<IdentifierDatabases, Collection<String>> map = 
      new HashMap<DatabaseIdentifiers.IdentifierDatabases, Collection<String>>();
    addXrefToDatabaseMap(map, pathway.getXREF());
    
    de.zbit.kegg.parser.pathway.Pathway keggPW = createPathwayInstance(comment, pathway, species, null, map);
    log.info("Converting pathway '" + keggPW.getTitle() + "'.");

    //TODO: it is not possible to define which image link to set, perhaps using datasource..., 
    // but too much databases to be conform for each
//    addImageLinkToKEGGpathway(species, pathway.getNAME(), keggPW);
    

    for (pathwayComponent pathComp : pathway.getPATHWAY_COMPONENTS()) {
      if (pathwayStep.class.isAssignableFrom(pathComp.getClass())){
        parsePathwayStep((pathwayStep)pathComp, keggPW, m, species);
      } else if (interaction.class.isAssignableFrom(pathComp.getClass())){
        parseInteraction((interaction) pathComp, keggPW, m, species);
      } else if (pathway.class.isAssignableFrom(pathComp.getClass())){
        createKEGGEntry((pathway) pathComp, keggPW, m, species, 
            EntryType.map, null, ",", null, null, false);  
      } else {
        log.log(Level.SEVERE, "Could not parse: '" + pathComp.getModelInterface() + "'.");
      }      
    }

    return keggPW;
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
   * {@link BioPAXL22KGML#getEntrezGeneIDsForPathways(List, String, Model)}
   * complexes are not treated right
   * 
   * @param pathways
   * @param species
   * @param m
   * @return
   */
  public List<BioPAXPathwayHolder> getEntrezGeneIDsForPathways(
      List<BioPAXPathwayHolder> pathways, String species, Model m) {
    log.info("Start parsing gene ids.");
    ProgressBar bar = new ProgressBar(pathways.size());

    Map<String, relationshipXref> xrefs = getMapFromSet(m.getObjects(relationshipXref.class));
    for (BioPAXPathwayHolder pw : pathways) {
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
  public static List<String> getListOfPathways(Model m){
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
  public pathway getPathwayByName(Model m, String name){
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
    String name = entity==null?null:entity.getNAME();

    if (name!=null && !name.isEmpty() && !(pathway.class.isAssignableFrom(entity.getClass()))) {
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
   * firstly set {@link BioPAXL22KGML#augmentOriginalKEGGpathway} to true, 
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
      Model m, Species species, openControlledVocabulary cv) {
    if (entity==null) return null;
    EntryExtended keggEntry = null;
    if (physicalEntity.class.isAssignableFrom(entity.getClass())) {
      keggEntry = parsePhysicalEntity((physicalEntity) entity, keggPW, m, species, cv, false);
    } else if (pathway.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.map, null, ",", null, cv, false);
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
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species, 
      openControlledVocabulary cv, boolean complexComponent) {
    EntryExtended keggEntry = null;

    if (complex.class.isAssignableFrom(entity.getClass())) {
      List<Integer> components = createComplexComponentList(((complex) entity).getCOMPONENTS(),
          keggPW, m, species);
      keggEntry = createKEGGEntry((entity) entity, keggPW, m, species, EntryType.group, null, "/",
          components, cv, complexComponent);
    } else if (sequenceEntity.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.other, EntryTypeExtended.gene, ",",
          null, cv, complexComponent);
    } else if (protein.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, EntryTypeExtended.protein,
          ",", null, cv, complexComponent);
    } else if (rna.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.other, EntryTypeExtended.rna, ",",
          null, cv, complexComponent);
    } else if (smallMolecule.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.compound, EntryTypeExtended.unknown,
          ",", null, cv, complexComponent);
    } else {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.other, EntryTypeExtended.unknown,
          ",", null, cv, complexComponent);
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
      if (physicalEntity.getPHYSICAL_ENTITY()==null) continue;
      EntryExtended keggEntry = parsePhysicalEntity(physicalEntity.getPHYSICAL_ENTITY(), keggPW, m, 
          species, physicalEntity.getCELLULAR_LOCATION(), true);
      if (keggEntry != null)
        components.add(keggEntry.getId());
    }

    return components;
  }

  /**
   * Adds the created {@link Entry} to the
   * {@link BioPAXL22KGML#bc2KeggEntry} map and to the
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
      List<Integer> components, openControlledVocabulary cv, boolean complexComponent) {
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
    
    if (entity.getSHORT_NAME()!=null && !entity.getSHORT_NAME().isEmpty())
        graphName = entity.getSHORT_NAME();
    else 
      graphName = names;
        

    String keggname = BioPAX2KGML.getKEGGName(identifiers, species);


    // create graphics
    Graphics graphics = null;
    if (eType.equals(EntryType.map)) {
      graphics = Graphics.createGraphicsForPathwayReference(graphName);
    } else if (eType.equals(EntryType.compound)) {
      graphics = Graphics.createGraphicsForCompound(graphName);
    } else if (eType.equals(EntryType.group) || eType.equals(EntryType.genes)) {
      graphics = Graphics.createGraphicsForGroupOrComplex(graphName);
    } else if (eType.equals(EntryType.gene) || eType.equals(EntryType.other)) {
      graphics = Graphics.createGraphicsForProtein(graphName);
    } else if (eType.equals(EntryType.compound)) {
      graphics = Graphics.createGraphicsForCompound(graphName);
    } else {
      graphics = new Graphics(graphName);
    }

    keggEntry = new EntryExtended(keggPW, getKeggEntryID(), keggname, eType, gType, graphics);    
    
    
    // set further information to the entry
    keggEntry.addDatabaseIdentifiers(identifiers);   
    
    if (components != null) {      
      keggEntry.addComponents(components);
    } 
    
    if (cv!=null && cv.getTERM().size()>0) {
      keggEntry.setCompartment(cv.getTERM().iterator().next());
    }

    // checking if entry already exists
    if (!augmentOriginalKEGGpathway){
      if (!complexComponent){
        Collection<de.zbit.kegg.parser.pathway.Entry> entries = keggPW.getEntries();
        if (entries != null && entries.size() > 0) {
          for (de.zbit.kegg.parser.pathway.Entry entry : entries) {        
              // important to ignore id, because this can differ from file to file
              if (((EntryExtended)entry).equalsWithoutIDNameReactionComparison(keggEntry)) {            
                keggEntry = (EntryExtended) entry;
                return keggEntry;
              }        
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

    // xrefs
    addXrefToDatabaseMap(map, entity.getXREF());
    
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
   * @param map
   * @param xref
   */
  private void addXrefToDatabaseMap(Map<IdentifierDatabases, Collection<String>> map, Set<xref> xrefs) {
    if (xrefs.size() != 0) {
      for (xref d : xrefs) {
        if (d!=null){
          if(d.getDB()!=null && !d.getDB().isEmpty()){
            IdentifierDatabases dbIdentifier = null;
            if (!d.getDB().equalsIgnoreCase("kegg")) {
              dbIdentifier = DatabaseIdentifiers.getDatabase(d.getDB());
            }
            if (dbIdentifier != null) {
              String id = d.getID();
              if (id!=null && !id.isEmpty()){
                Utils.addToMapOfSets(map, dbIdentifier, id);
              }
            } else if (d.getDB().equalsIgnoreCase("LL")) { // special case in PID database files
              if (d.getID()!=null) {
                String id = d.getID();
                if (!id.isEmpty())
                  Utils.addToMapOfSets(map, IdentifierDatabases.EntrezGene, id);              
              }
            } else if (d.getDB().equalsIgnoreCase("KEGG")) { // Infer correct KEGG db
              if (d.getID()!=null) {
                String id = d.getID();
                if (!id.isEmpty()){
                  IdentifierDatabases db = DatabaseIdentifierTools.getKEGGdbFromID(id);
                  if (db!=null) {
                    Utils.addToMapOfSets(map, db, id);
                  }
                }
              }
            }
          }
        }
      }
    }
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
    } else if (interaction.class.isAssignableFrom(entity.getClass())) {
      parseInteractionFromInteraction(entity, keggPW, m, species, null);
    } else {      
      log.severe("Unknonw entity type: " + entity.getModelInterface() + "-" + entity.getRDFId());
      System.exit(1);
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
          keggPW, m, species, ((catalysis)entity).getXREF());
    } else if (modulation.class.isAssignableFrom(entity.getClass())) {
      createKEGGReactionRelation(((modulation) entity).getCONTROLLER(),
          ((modulation) entity).getCONTROLLED(),
          getSubtype(((modulation) entity).getCONTROL_TYPE()), keggPW, m, species, 
          ((modulation)entity).getXREF());
    } else {
      createKEGGReactionRelation(((control) entity).getCONTROLLER(),
          ((control) entity).getCONTROLLED(), getSubtype(((control) entity).getCONTROL_TYPE()),
          keggPW, m, species, ((control)entity).getXREF());
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
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species, Set<xref> xrefs) {
    EntryExtended keggEntry1 = null;
    RelationType relType = null;

    if (controllers.size() >= 1) {
      for (physicalEntityParticipant controller : controllers) {
        if (controller.getPHYSICAL_ENTITY()==null) {
          continue;
        }
        keggEntry1 = parseEntity(controller.getPHYSICAL_ENTITY(), keggPW, m, species, 
            controller.getCELLULAR_LOCATION());
        
        if(keggEntry1!=null && keggEntry1.getType().equals(EntryType.map)){
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
                        ((biochemicalReaction) con).getRIGHT(), keggPW, m, species, xrefs);

                  } catch (ClassCastException e) {
                    try {
                      r = createKEGGReaction(((complexAssembly) con).getLEFT(),
                          ((complexAssembly) con).getRIGHT(), keggPW, m, species, xrefs);
                    } catch (ClassCastException e2) {
                      r = createKEGGReaction(((transportWithBiochemicalReaction) con).getLEFT(),
                          ((transportWithBiochemicalReaction) con).getRIGHT(), keggPW, m, species, xrefs);
                    }
                  }

                  if (relType.equals(RelationType.maplink) && r!=null) {
                    for (ReactionComponent rc : r.getSubstrates()) {
                      createKEGGRelation(keggPW, keggEntry1.getId(), rc.getId(), relType, subtype, xrefs);
                    }
                  } else if (relType.equals(RelationType.PPrel) && r!=null) {
                    keggEntry1.appendReaction(r.getName());
                  }  
                }
              } else if (transport.class.isAssignableFrom(con.getClass())) {
                Reaction r = createKEGGReaction(((transport) con).getLEFT(),
                    ((transport) con).getRIGHT(), keggPW, m, species, ((transport)con).getXREF());                
                if (relType.equals(RelationType.maplink) && r!=null) {
                  for (ReactionComponent rc : r.getSubstrates()) {
                    createKEGGRelation(keggPW, keggEntry1.getId(), rc.getId(), relType, subtype, xrefs);
                  }
                } else if (relType.equals(RelationType.PPrel) && r!=null) {
                  keggEntry1.appendReaction(r.getName());
                }
              } else if (conversion.class.isAssignableFrom(con.getClass())){
                  List<Relation> rels =  createKEGGRelations(con.getLEFT(), con.getRIGHT(), keggPW, m, species, 
                    RelationType.PPrel, getSubtype(con.getINTERACTION_TYPE()), xrefs);
                  for (Relation rel : rels) {
                    if (rel !=null)
                      createKEGGRelation(keggPW, keggEntry1.getId(), rel.getEntry2(), relType, subtype, xrefs);
                  }
              } else {
                log.severe("Not programmed case: controlled interface '" + con.getModelInterface()
                    + "'");
                System.exit(1);
              }
            } else if (pathway.class.isAssignableFrom(process.getClass())) {
              EntryExtended keggEntry2 = createKEGGEntry((pathway) process, keggPW, m, species,
                  EntryType.map, null, ",", null, null, false);
              if (keggEntry2 !=null)
                createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(),
                  relType, subtype, xrefs);
            } else if (interaction.class.isAssignableFrom(process.getClass())) {                            
              parseInteractionFromInteraction( ((interaction)process), keggPW, m, species, 
                  keggEntry1);
            } else {
              log.severe("Process: " + process.getModelInterface() + "-This should not happen!");
              System.exit(1);
            }
            // ControlType (0 or 1) - up to now ignored

            // Cofactor = PhysicalEntity (0 or more) - up to now ignored

            // CatalysisDirection - up to now ignored
          }
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
  private void parseInteractionFromInteraction(interaction inter, Pathway keggPW, 
      Model m, Species species, EntryExtended baseEntry) {
    EntryExtended keggEntry1 = null;
    List<InteractionParticipant> participants = Utils.iterableToList(inter.getPARTICIPANTS());
    if (participants.size() > 1) {
      for (int i = 0; i < participants.size(); i++) {
        if (pathway.class.isAssignableFrom(participants.get(i).getClass())) {
          keggEntry1 = parseEntity(((pathway) participants.get(i)), keggPW, m, species, null);
        } else if (physicalEntityParticipant.class.isAssignableFrom(participants.get(i).getClass())) {
          keggEntry1 = parseEntity(
              ((physicalEntityParticipant) participants.get(i)).getPHYSICAL_ENTITY(), keggPW, m,
              species, ((physicalEntityParticipant) participants.get(i)).getCELLULAR_LOCATION());
        } else {
          log.log(Level.SEVERE, "1 This should not happen: '"
              + participants.get(i).getModelInterface() + "'.");
          System.exit(1);
        }
        for (int j = 1; j < participants.size(); j++) {
          EntryExtended keggEntry2 = null;
          if (pathway.class.isAssignableFrom(participants.get(j).getClass())) {
            keggEntry2 = parseEntity(((pathway) participants.get(j)), keggPW, m, species, null);
          } else if (physicalEntityParticipant.class.isAssignableFrom(participants.get(j)
              .getClass())) {
            keggEntry2 = parseEntity(
                ((physicalEntityParticipant) participants.get(j)).getPHYSICAL_ENTITY(), keggPW, m,
                species, ((physicalEntityParticipant) participants.get(j)).getCELLULAR_LOCATION());
          } else {
            log.log(Level.SEVERE, "2 This should not happen: '"
                + participants.get(j).getModelInterface() + "'.");
            System.exit(1);
          }

          if (keggEntry1!=null && keggEntry2!=null) {
            createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), RelationType.other, 
                null, inter.getXREF());
          }
          if (baseEntry!=null && keggEntry1!=null) {
            createKEGGRelation(keggPW, baseEntry.getId(), keggEntry1.getId(), RelationType.other, 
                null, inter.getXREF());
          }
        }

      }
    } else if (participants.size() > 0) {
      if (pathway.class.isAssignableFrom(participants.get(0).getClass())) {
        keggEntry1 = parseEntity(((pathway) participants.get(0)), keggPW, m, species, null);
      } else if (physicalEntityParticipant.class.isAssignableFrom(participants.get(0).getClass())) {
        keggEntry1 = parseEntity(
            ((physicalEntityParticipant) participants.get(0)).getPHYSICAL_ENTITY(), keggPW, m,
            species, ((physicalEntityParticipant) participants.get(0)).getCELLULAR_LOCATION());
      } else {
        log.log(Level.SEVERE, "3 - This should not happen: '"
            + participants.get(0).getModelInterface() + "'.");
        System.exit(1);
      }
      if (baseEntry!=null && keggEntry1!=null) {
        createKEGGRelation(keggPW, baseEntry.getId(), keggEntry1.getId(), RelationType.other, 
            null, inter.getXREF());
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
      
      if (baseEntry!=null && keggEntry1!=null) {
        createKEGGRelation(keggPW, baseEntry.getId(), keggEntry1.getId(), RelationType.other, 
            null, inter.getXREF());
      }
    }
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
        createKEGGReaction(((complexAssembly) entity).getLEFT(), ((complexAssembly) entity).getRIGHT(), 
            keggPW, m, species, ((complexAssembly)entity).getXREF());
    } else if (biochemicalReaction.class.isAssignableFrom(entity.getClass())) {
      if  (!augmentOriginalKEGGpathway)
        createKEGGReaction(((biochemicalReaction) entity).getLEFT(),
          ((biochemicalReaction) entity).getRIGHT(), keggPW, m, species,
          ((biochemicalReaction) entity).getXREF());
    } else if (transport.class.isAssignableFrom(entity.getClass())) {
      createKEGGReaction(((transport) entity).getLEFT(), ((transport) entity).getRIGHT(), keggPW,
          m, species, ((transport) entity).getXREF());
    } else if (transportWithBiochemicalReaction.class.isAssignableFrom(entity.getClass())) {
      if  (!augmentOriginalKEGGpathway)
        // BiochemicalReaction br = (TransportWithBiochemicalReaction) entity;
        // deltaG, deltaH, deltaS, ec, and KEQ are ignored
        createKEGGReaction(((biochemicalReaction) entity).getLEFT(),
          ((biochemicalReaction) entity).getRIGHT(), keggPW, m, species,
          ((biochemicalReaction) entity).getXREF());
    } else if (conversion.class.isAssignableFrom(entity.getClass())){
      createKEGGRelations(((conversion)entity).getLEFT(), ((conversion)entity).getRIGHT(), keggPW, 
          m, species, RelationType.PPrel, getSubtype(((conversion)entity).getINTERACTION_TYPE()),
          ((conversion)entity).getXREF());
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
      SubType subType, Set<xref> xrefs) {

    List<Relation> relations = new ArrayList<Relation>();

    for (physicalEntityParticipant left : set) {
      if (left.getPHYSICAL_ENTITY()==null) continue;
      EntryExtended keggEntry1 = parsePhysicalEntity(left.getPHYSICAL_ENTITY(), keggPW, m, 
          species, left.getCELLULAR_LOCATION(), false);
      
      if (keggEntry1 !=null){
        for (physicalEntityParticipant right : set2) {
          if (right.getPHYSICAL_ENTITY()==null) continue;
          EntryExtended keggEntry2 = parsePhysicalEntity(right.getPHYSICAL_ENTITY(), keggPW, m, 
              species, right.getCELLULAR_LOCATION(), false);
          if (keggEntry1 !=null){
            Relation r = createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), 
                type, subType, xrefs);
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
      int keggEntry2Id, RelationType type, SubType subType, Set<xref> xrefs) {
    ArrayList<Relation> existingRels = keggPW.getRelations();
    Relation r = null;

    // get the references
    Map<IdentifierDatabases, Collection<String>> map = 
      new HashMap<DatabaseIdentifiers.IdentifierDatabases, Collection<String>>();
    addXrefToDatabaseMap(map, xrefs);
    
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
            if (augmentOriginalKEGGpathway && added){
              addedSubTypes++;
              r.setSource("KEGG_AND_BIOCARTA");
            }
            return r;
          }          
        }
      }

      r = new Relation(keggEntry1Id, keggEntry2Id, type, subType);
    } else {
      r = new Relation(keggEntry1Id, keggEntry2Id, type, subType);
    }
    
    r.addDatabaseIdentifiers(map);
    // If we are here, the relation r is NOVEL AND NOT CURRENTLY IN THE PATHWAY

    // Add the relation to the pathway
    
    
    if (augmentOriginalKEGGpathway){
      if(keggPW.getEntryForId(keggEntry1Id)!=null && keggPW.getEntryForId(keggEntry2Id)!=null) {
        // Only add relations if nodes for the relation are present.
        if (keggEntry1Id != keggEntry2Id){
          r.setSource("BIOCARTA");
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
      Pathway keggPW, Model m, Species species, Set<xref> xrefs) {
    List<ReactionComponent> products = new ArrayList<ReactionComponent>();
    List<ReactionComponent> substrates = new ArrayList<ReactionComponent>();

    for (physicalEntityParticipant left : lefts) {
      if (left.getPHYSICAL_ENTITY()==null) continue;
      EntryExtended keggEntry =  parsePhysicalEntity(left.getPHYSICAL_ENTITY(), keggPW, m, 
          species, left.getCELLULAR_LOCATION(), false);
         
      if (keggEntry != null) {
        ReactionComponent rc = new ReactionComponent(keggEntry.getId(), keggEntry.getName());
        Integer stoich = (int) left.getSTOICHIOMETRIC_COEFFICIENT();
        if (stoich >0 && stoich<Integer.MAX_VALUE) {
          rc.setStoichiometry(stoich);
        }
        products.add(rc);
      }
    }
    
    for (physicalEntityParticipant right : rights) {
      if (right.getPHYSICAL_ENTITY()==null) continue;
      EntryExtended keggEntry = parsePhysicalEntity(right.getPHYSICAL_ENTITY(), keggPW, m, 
          species, right.getCELLULAR_LOCATION(), false);
      if (keggEntry != null) {
        ReactionComponent rc = new ReactionComponent(keggEntry.getId(), keggEntry.getName());   
        Integer stoich = (int)right.getSTOICHIOMETRIC_COEFFICIENT();
        if (stoich >0 && stoich<Integer.MAX_VALUE) {
          rc.setStoichiometry(stoich);
        }
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
      // get the references
      Map<IdentifierDatabases, Collection<String>> map = 
        new HashMap<DatabaseIdentifiers.IdentifierDatabases, Collection<String>>();
      addXrefToDatabaseMap(map, xrefs);
      
      r = new Reaction(keggPW, getReactionName(), ReactionType.other);
      r.addProducts(products);
      r.addSubstrates(substrates);
      r.addDatabaseIdentifiers(map);
      keggPW.addReaction(r);
    }

    return r;
  }
}
