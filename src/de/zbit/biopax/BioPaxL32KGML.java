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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
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
import org.biopax.paxtools.model.level3.BioSource;
import org.biopax.paxtools.model.level3.BiochemicalReaction;
import org.biopax.paxtools.model.level3.Catalysis;
import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.model.level3.ComplexAssembly;
import org.biopax.paxtools.model.level3.Control;
import org.biopax.paxtools.model.level3.ControlType;
import org.biopax.paxtools.model.level3.Controller;
import org.biopax.paxtools.model.level3.Conversion;
import org.biopax.paxtools.model.level3.ConversionDirectionType;
import org.biopax.paxtools.model.level3.Degradation;
import org.biopax.paxtools.model.level3.Dna;
import org.biopax.paxtools.model.level3.DnaRegion;
import org.biopax.paxtools.model.level3.Entity;
import org.biopax.paxtools.model.level3.Gene;
import org.biopax.paxtools.model.level3.GeneticInteraction;
import org.biopax.paxtools.model.level3.Interaction;
import org.biopax.paxtools.model.level3.InteractionVocabulary;
import org.biopax.paxtools.model.level3.Modulation;
import org.biopax.paxtools.model.level3.MolecularInteraction;
import org.biopax.paxtools.model.level3.Pathway;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.model.level3.Protein;
import org.biopax.paxtools.model.level3.Provenance;
import org.biopax.paxtools.model.level3.RelationshipXref;
import org.biopax.paxtools.model.level3.Rna;
import org.biopax.paxtools.model.level3.RnaRegion;
import org.biopax.paxtools.model.level3.SmallMolecule;
import org.biopax.paxtools.model.level3.Stoichiometry;
import org.biopax.paxtools.model.level3.TemplateReaction;
import org.biopax.paxtools.model.level3.TemplateReactionRegulation;
import org.biopax.paxtools.model.level3.Transport;
import org.biopax.paxtools.model.level3.TransportWithBiochemicalReaction;
import org.biopax.paxtools.model.level3.Xref;

import de.zbit.kegg.KGMLWriter;
import de.zbit.kegg.parser.pathway.Entry;
import de.zbit.kegg.parser.pathway.EntryType;
import de.zbit.kegg.parser.pathway.Graphics;
import de.zbit.kegg.parser.pathway.Reaction;
import de.zbit.kegg.parser.pathway.ReactionComponent;
import de.zbit.kegg.parser.pathway.ReactionType;
import de.zbit.kegg.parser.pathway.Relation;
import de.zbit.kegg.parser.pathway.RelationType;
import de.zbit.kegg.parser.pathway.SubType;
import de.zbit.kegg.parser.pathway.ext.EntryExtended;
import de.zbit.kegg.parser.pathway.ext.GeneType;
import de.zbit.mapper.GeneSymbol2GeneIDMapper;
import de.zbit.util.DatabaseIdentifiers;
import de.zbit.util.SortedArrayList;
import de.zbit.util.Species;
import de.zbit.util.Utils;
import de.zbit.util.progressbar.ProgressBar;

/**
 * This class works with PaxTools. It is used to fetch information out of a
 * level 3 BioCarta files. Example files could be downloaded from
 * http://pid.nci.nih.gov/download.shtml
 * 
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class BioPaxL32KGML extends BioPax2KGML {
  
  public static final Logger log = Logger.getLogger(BioPaxL32KGML.class.getName());

  /**
   * This method creates for each pathway in the model a KGML file with the
   * pathway name and saves the pathways in an default folder see 
   * {@link BioPax2KGML#createDefaultFolder(org.biopax.paxtools.model.BioPAXLevel)}
   * 
   * @param m
   */
  public void createKGMLsForPathways(Model m, Set<Pathway> pathways, boolean writeEntryExtended) {
    String folder = createDefaultFolder(m.getLevel());
    createKGMLsForPathways(m, folder, pathways, writeEntryExtended);
  }
  
  /**
   * This method creates for each pathway in the model a KGML file with the
   * pathway name
   * 
   * @param m
   */
  public void createKGMLsForPathways(Model m, String folder, Set<Pathway> pathways, boolean writeEntryExtended) {
    log.info("Creating for each pathway a KGML file.");

    Collection<de.zbit.kegg.parser.pathway.Pathway> keggPWs = parsePathways(m, pathways);

    for (de.zbit.kegg.parser.pathway.Pathway keggPW : keggPWs) {
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
  private Collection<de.zbit.kegg.parser.pathway.Pathway> parsePathways(Model m, Set<Pathway> pathways) {
    Species oldSpecies = new Species("Homo sapiens", "_HUMAN", "human", "hsa", 9606);
    initalizeMappers(oldSpecies);

    Collection<de.zbit.kegg.parser.pathway.Pathway> keggPWs = 
      new ArrayList<de.zbit.kegg.parser.pathway.Pathway>();

    for (Pathway pathway : pathways) {
      // determine the pathway organism
      Species newSpecies = determineSpecies(pathway.getOrganism());
      if(!newSpecies.equals(oldSpecies)){
        initalizeMappers(newSpecies);
        oldSpecies = newSpecies;
      }
      
      keggPWs.add(parsePathway(m, pathway, oldSpecies));
    }
    
    return keggPWs;
  }

  /**
   * determines the species of the pathway and returns {@link Species}
   * the default species it "Homo sapiens"
   * @param pathway
   * @return
   */
  protected Species determineSpecies(BioSource pwOrg){
    Species species = new Species("Homo sapiens", "_HUMAN", "human", "hsa", 9606),
    detSpecies = null;

    if (pwOrg != null) {
      Set<Xref> references = pwOrg.getXref();
      if(references!=null && references.size()>0){
        for (Xref xref : references) {
          if(xref.getDb().toLowerCase().equals(
              DatabaseIdentifiers.IdentifierDatabases.NCBI_Taxonomy.toString().toLowerCase())){
            detSpecies = Species.search(allSpecies, xref.getId(), Species.NCBI_TAX_ID);
          }
        }
      }
      
      if (pwOrg.getName() != null) {
        String newSpecies = Utils.iterableToList(pwOrg.getName()).get(0);
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
   * returns a list of all pathways containing pathway components
   * @param m
   * @return
   */
  protected List<String> getPathways(Model m){
    List<String> pws = new SortedArrayList<String>();
    
    Set<Pathway> list = m.getObjects(Pathway.class);
    for (Pathway pw : list) {
      if(pw.getPathwayComponent().size()>0)
        pws.add(getPathwayName(pw));
    }
    
    return pws;
  }
  
  /**
   * 
   * @param m
   * @param name
   * @return the BioPaxPathway with the specific name
   */
  protected Pathway getPathwayByName(Model m, String name){
    Pathway pw = null;
    Set<Pathway> list = m.getObjects(Pathway.class);
    
    for (Pathway p : list) {
      if(getPathwayName(p).equals(name))
        return p;
    }
    
    return pw;
  }  
  
  /**
   * parses the biopax pathway
   * @param m
   * @param pathway
   * @param species
   * @return
   */
  protected de.zbit.kegg.parser.pathway.Pathway parsePathway(Model m, Pathway pathway, Species species) {
    if(geneSymbolMapper==null && geneIDKEGGmapper==null){
      initalizeMappers(species);
    }
       
    //create the pathway
    int number = getKeggPathwayNumber(pathway.getRDFId());
    String sourceDB = getSourceDB(pathway.getDataSource());
    String pwName = getPathwayName(pathway);
 
    de.zbit.kegg.parser.pathway.Pathway keggPW = new de.zbit.kegg.parser.pathway.Pathway(
        sourceDB + String.valueOf(number), species.getKeggAbbr(), number,
        pwName);
    log.info("Converting pathway '" + pwName + "'.");
    if (!sourceDB.isEmpty())
      keggPW.setLink(sourceDB);

    //TODO: it is not possible to define which link to set, perhaps using datasource..., 
    // but too much databases to be conform for each
//    addImageLinkToKEGGpathway(species, getNameWithoutBlanks(pathway.getName()), keggPW);

    for (org.biopax.paxtools.model.level3.Process interaction : pathway.getPathwayComponent()) {
      parseInteraction((Interaction) interaction, keggPW, m, species);
    }
    
    return keggPW;
  }
  
  /**
   * sets the source of the data if available to the class
   * @param sources
   * @param keggPW
   */
  private String getSourceDB(Set<Provenance> sources) {
    String source = "";
    if (sources!=null && sources.size()>0){
      for (Provenance p : sources) {
        if(p.getComment()!=null && p.getComment().size()>0){
          source = p.getComment().iterator().next();          
        }
      }
    }
    
    return source;
  }

  /**
   * firstly set {@link BioPaxL32KGML#augmentOriginalKEGGpathway} to true, 
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
    Set<Pathway> pathways = m.getObjects(Pathway.class);
    Species species = Species.search(allSpecies, p.getOrg(), Species.KEGG_ABBR);
    if(species != null){
      initalizeMappers(species);
      
      for (Pathway pathway : pathways) {    
        // determine the pathway organism
        BioSource pwOrg = pathway.getOrganism();  
        if (pwOrg != null) {
          if (pwOrg.getName() != null) {
            String pathwaySpecies = Utils.iterableToList(pwOrg.getName()).get(0);
            if (pathwaySpecies.equals(species.getScientificName())) {
              for (org.biopax.paxtools.model.level3.Process interaction : pathway.getPathwayComponent()) {
                parseInteraction((Interaction) interaction, p, m, species);
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
   * parse an BioPax file by calling 
   * {@link BioPaxL22KGML#createKGMLForBioPaxFile(Model, String, String)}
   * the folder is determined with 
   * {@link BioPax2KGML#createDefaultFolder(org.biopax.paxtools.model.BioPAXLevel)}
   * @param m
   * @param fileName
   */
  public void createKGMLForBioPaxFile(Model m, String fileName, boolean writeEntryExtended) {
    String folder = createDefaultFolder(m.getLevel());
    createKGMLForBioPaxFile(m, fileName, folder, writeEntryExtended);
  }
  
  /**
   * The methods parse a BioPax file which contains no <bp>Pathway: ....</bp> tag
   * 
   * @param m
   * @param fileName used as pathway name and title 
   * @param folder where the KGML is saved
   */
  public void createKGMLForBioPaxFile(Model m, String pathwayName, String folder, boolean writeEntryExtended) {
    // determine the organism
    String defaultSpecies = "Homo sapiens", newSpecies = defaultSpecies;
 // determine the organism
    Species species = new Species("Homo sapiens", "_HUMAN", "human", "hsa", 9606);
    
    Set<BioSource> orgs = m.getObjects(BioSource.class);
    if (orgs != null && orgs.size() > 0) {
      BioSource org = orgs.iterator().next();
      species = determineSpecies(org);      
    } else {
      log.info("No pathway species could be determined, default species '"
          + species.getCommonName() + "' is used.");
    }
    initalizeMappers(species);    
    
    int pathNo = determineKEGGPathwayNumber(pathwayName);

    de.zbit.kegg.parser.pathway.Pathway keggPW = new de.zbit.kegg.parser.pathway.Pathway(
        pathwayName, species.getKeggAbbr(), pathNo, pathwayName);

    log.info("Converting pathway '" + pathwayName + "'.");


    //TODO: it is not possible to define which link to set, perhaps using datasource..., 
    // but too much databases to be conform for each
//    addImageLinkToKEGGpathway(species, getNameWithoutBlanks(pathway.getName()), keggPW);
    
    for (Entity entity : m.getObjects(Entity.class)) {
      parseEntity(entity, keggPW, m, species);
    }

    String fileName = folder + KGMLWriter.createFileName(keggPW);
    KGMLWriter.writeKGML(keggPW, fileName, writeEntryExtended);
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
  private EntryExtended parseEntity(Entity entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species) {
    EntryExtended keggEntry = null;
    if (PhysicalEntity.class.isAssignableFrom(entity.getClass())) {
      keggEntry = parsePhysicalEntity((PhysicalEntity) entity, keggPW, m, species);
    } else if (Interaction.class.isAssignableFrom(entity.getClass())) {
      parseInteraction((Interaction) entity, keggPW, m, species);
    } else if (Pathway.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.map, null, ",", null);
    } else if (Gene.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, null, ",", null);
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
  private EntryExtended parsePhysicalEntity(PhysicalEntity entity,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    EntryExtended keggEntry = null;

    if (Complex.class.isAssignableFrom(entity.getClass())) {
      List<Integer> components = createComplexComponentList(((Complex) entity).getComponent(),
          keggPW, m, species);
      keggEntry = createKEGGEntry((Entity) entity, keggPW, m, species, EntryType.group, null, "/",
          components);
    } else if (Dna.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, GeneType.dna, ",",
          null);
    } else if (DnaRegion.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, GeneType.dna_region,
          ",", null);
    } else if (Protein.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, GeneType.protein,
          ",", null);
    } else if (Rna.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, GeneType.rna, ",",
          null);
    } else if (RnaRegion.class.isAssignableFrom(entity.getClass())) {
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, GeneType.rna_region,
          ",", null);
    } else if (SmallMolecule.class.isAssignableFrom(entity.getClass())) {
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
   * @param complexEntries
   * @param keggPW
   * @param m
   * @param species
   * @return
   */
  private List<Integer> createComplexComponentList(Set<PhysicalEntity> complexEntries,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    List<Integer> components = new ArrayList<Integer>();

    for (PhysicalEntity physicalEntity : complexEntries) {
      EntryExtended keggEntry = parsePhysicalEntity(physicalEntity, keggPW, m, species);
      if (keggEntry != null)
        components.add(keggEntry.getId());
    }

    return components;
  }

  /**
   * This method maps the gene id to the gene symbol of the entered entity. This
   * method is not able to map complexes call for that the method
   * {@link BioPaxL32KGML#getEntrezGeneIDsOfComplex(Complex, Map)} It's an
   * advantage to preprocess the entities to exclusively having entities with a
   * name call therefore the method
   * {@link BioPaxL32KGML#getEntitiesWithName(Entity)}
   * 
   * The default for a gene id which is not found is -1
   * 
   * @param entity
   * @param xrefs
   * @return
   */
  protected static int getEntrezGeneID(Entity entity, Map<String, RelationshipXref> xrefs) {
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
      Set<String> names = entity.getName();
      if (names != null && names.size() > 0) {
        Map<String, Entity> geneSymbols = new HashMap<String, Entity>();
        for (String name : names) {
          geneSymbols.put(name, entity);
        }
        if (getEntrezGeneIDOverGeneSymbol(geneSymbols).get(entity) != null)
          geneID = getEntrezGeneIDOverGeneSymbol(geneSymbols).get(entity).intValue();
      }
    }

    return geneID;
  }
  
  /**
   * This method returns a list of all biocarta pathways with the containing
   * gene IDs
   * 
   * @param species
   * @param model
   * @return
   */
  public List<BioPaxPathwayHolder> getPathwaysWithEntrezGeneID(String species, Model m) {
    return getEntrezGeneIDsForPathways(getPathwayEntities(species, m), species, m);
  }
  
  /**
   * creates a list with {@link BioPaxPathwayHolder} and the entities of the
   * pathways
   * 
   * @param species
   * @param m
   * @return
   */
  private List<BioPaxPathwayHolder> getPathwayEntities(String species, Model m) {
    List<BioPaxPathwayHolder> pathways = new ArrayList<BioPaxPathwayHolder>();

    for (Entity entity : m.getObjects(Entity.class)) {
      for (Interaction string : entity.getParticipantOf()) {
        for (Pathway pw : string.getPathwayComponentOf()) {
          BioPaxPathwayHolder helper = new BioPaxPathwayHolder(pw.getRDFId(),
              getPathwayName(pw));
          int index = pathways.indexOf(helper);
          if (index > -1) {
            helper = pathways.get(index);
          } else {
            pathways.add(helper);
          }
          helper.addEntity(entity);
        }
      }
    }

    return pathways;
  }
  
  /**
   * creates out of the pathway set the corresponding name of the pathway, the
   * name with blanks is preferred
   * 
   * @param pathway
   * @return
   */
  private String getPathwayName(Entity pathway) {
    Set<String> names = pathway.getName();
    String name = "";

    if (names != null && names.size() > 0) {
      List<String> names2 = Utils.iterableToList(names);
      if (names2.size() == 1) {
        name = names2.get(0);
      } else if (names2.size() == 2) {
        String h = names2.get(0);
        if (h.contains(" ")) {
          name = h;
        } else {
          name = names2.get(1);
        }
      } else {
        StringBuffer sb = new StringBuffer();
        for (String n : names2) {
          sb.append(n + ";");
        }
        name = sb.toString();
      }
    }
    return name;
  }
  
  /**
   * method to create a pathay file out of the whole owl file
   * 
   * @param species
   * @param m
   * @param pathway
   * @throws IOException
   */
  private void writePathwayOwlFile(Pathway pathway, Model m, Species species) throws IOException {
    String name = getPathwayName(pathway);
    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(name + ".owl")));

    for (Entity entity : m.getObjects(Entity.class)) {
      for (Interaction string : entity.getParticipantOf()) {
        for (Pathway pw : string.getPathwayComponentOf()) {
          if (pw.equals(pathway)) {
            bw.append(entity.getModelInterface() + "\t" + entity.getRDFId() + "\n");
          }
        }
      }
    }

    bw.close();
  }

  /**
   * This method returns a list of all model-pathways with the containing gene
   * IDs
   * 
   * @param species
   * @param model
   * @return
   */
  public List<BioPaxPathwayHolder> getPathwaysWithGeneID(String species, Model m) {
    List<BioPaxPathwayHolder> pathways = new ArrayList<BioPaxPathwayHolder>();

    for (Entity entity : m.getObjects(Entity.class)) {
      for (Interaction string : entity.getParticipantOf()) {
        for (Pathway pw : string.getPathwayComponentOf()) {
          BioPaxPathwayHolder helper = new BioPaxPathwayHolder(pw.getRDFId());
          int index = pathways.indexOf(helper);
          if (index > -1) {
            helper = pathways.get(index);
          } else {
            pathways.add(helper);
          }
          helper.addEntity(entity);
        }
      }
    }

    Map<String, RelationshipXref> xrefs = getMapFromSet(m.getObjects(RelationshipXref.class));
    for (BioPaxPathwayHolder pw : pathways) {
      if (pw.getPathwayName().equals("http://pid.nci.nih.gov/biopaxpid_9796")) {
        log.log(Level.FINER, "Pathway: " + pw.getPathwayName() + ": " + pw.getNoOfEntities());
        Set<BioPAXElement> pwEntities = new HashSet<BioPAXElement>();
        for (BioPAXElement entity : pw.entities) {
          pwEntities.addAll(getEntitiesWithName((Entity)entity));
          if (!(Pathway.class.isAssignableFrom(entity.getClass())))
            log.log(Level.FINER,
                "--Input: " + entity.getRDFId() + "\t" + entity.getModelInterface());
        }

        Map<Entity, Integer> geneIDs = getEntrezGeneIDs(pwEntities, species, xrefs);
        for (java.util.Map.Entry<Entity, Integer> entity : geneIDs.entrySet()) {
          log.log(Level.FINER, "----res: " + entity.getValue() + " " + entity.getValue());
          pw.addGeneID(entity.getValue());
        }
      }
    }

    return pathways;
  }
  
  /**
   * 
   * @param xref
   * @return Integer of the entered xref or null
   */
  public static Integer getEntrezGeneIDFromDBxref(RelationshipXref xref) {
    if (xref.getDb().equals("LL")) {
      log.info(xref.getRDFId() + "|" + xref.getId());
      return Integer.parseInt(xref.getId());
    } else {
      return null;
    }
  }
  
  /**
   * Adds the created {@link Entry} to the
   * {@link BioPaxL32KGML#bc2KeggEntry} map and to the
   * {@link de.zbit.kegg.parser.pathway.Pathway}
   * 
   * @param entity
   * @param keggPW
   * @param mapper
   * @param m
   * @return
   */
  private EntryExtended createKEGGEntry(Entity entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species, EntryType eType, GeneType gType, String graphNameSeparator,
      List<Integer> components) {
    EntryExtended keggEntry;

    // creating new KEGG entry
    String keggname = null;
    Integer geneID = getEntrezGeneID(entity, getMapFromSet(m.getObjects(RelationshipXref.class)));
    if (geneID != -1) {
      keggname = mapGeneIDToKEGGID(geneID, species);
    } else {
      keggname = getKEGGUnkownName();
    }

    String graphName = "";
    Set<String> names = entity.getName();
    StringBuffer name = new StringBuffer();
    if (names != null && names.size() > 0) {
      for (String n : names) {
        n = n.trim();
        n = n.replace(" ", "_");
        name.append(n);
        name.append(graphNameSeparator);
        name.append(" ");
      }
      name.delete(name.length() - 2, name.length() - 1);
    }
    graphName = name.toString().trim();

    Graphics graphics = null;
    if (eType.equals(EntryType.map))
      graphics = Graphics.createGraphicsForPathwayReference(graphName);
    else if (eType.equals(EntryType.compound))
      graphics = Graphics.createGraphicsForCompound(graphName);
    else
      graphics = new Graphics(graphName);

    keggEntry = new EntryExtended(keggPW, getKeggEntryID(), keggname, eType, gType, graphics);
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
   * parse a BioPax Interaction element
   * 
   * @param entity
   * @param keggPW
   * @param mapper
   * @param m
   * @param species
   */
  private void parseInteraction(Interaction entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species) {
    if (Control.class.isAssignableFrom(entity.getClass())) {
      parseControl((Control) entity, keggPW, m, species);
    } else if (Conversion.class.isAssignableFrom(entity.getClass())) {
      parseConversion((Conversion) entity, keggPW, m, species);
    } else if (GeneticInteraction.class.isAssignableFrom(entity.getClass())) {
      createKEGGRelationForParticipantList(
          Utils.iterableToList(((GeneticInteraction) entity).getParticipant()), keggPW, m,
          species, RelationType.GErel, new SubType(SubType.ASSOCIATION));
    } else if (MolecularInteraction.class.isAssignableFrom(entity.getClass())) {
      createKEGGRelationForParticipantList(
          Utils.iterableToList(((MolecularInteraction) entity).getParticipant()), keggPW, m,
          species, RelationType.PPrel, new SubType(SubType.INDIRECT_EFFECT));
    } else if (TemplateReaction.class.isAssignableFrom(entity.getClass())) {
      createKEGGRelationForTemplateReaction((TemplateReaction) entity, keggPW, m, species);
    } else {
      createKEGGEntry(entity, keggPW, m, species);
    }
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
  private EntryExtended createKEGGEntry(Interaction inter, 
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    EntryExtended keggEntry1 = null;
    List<Entity> participants = Utils.iterableToList(inter.getParticipant());
    if (participants.size() > 1) {
      for (int i = 0; i < participants.size(); i++) {
        if (Pathway.class.isAssignableFrom(participants.get(i).getClass())) {
          keggEntry1 = parseEntity(((Pathway) participants.get(i)), keggPW, m, species);
        } else if (PhysicalEntity.class.isAssignableFrom(participants.get(i).getClass())) {
          keggEntry1 = parseEntity(
              ((PhysicalEntity) participants.get(i)), keggPW, m,
              species);
        } else {
          log.log(Level.SEVERE, "1 This should not happen: '"
              + participants.get(i).getModelInterface() + "'.");
          System.exit(1);
        }
        for (int j = 1; j < participants.size(); j++) {
          EntryExtended keggEntry2 = null;
          if (Pathway.class.isAssignableFrom(participants.get(j).getClass())) {
            keggEntry2 = parseEntity(((Pathway) participants.get(j)), keggPW, m, species);
          } else if (PhysicalEntity.class.isAssignableFrom(participants.get(j)
              .getClass())) {
            keggEntry2 = parseEntity(
                ((PhysicalEntity) participants.get(j)), keggPW, m,
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
      if (Pathway.class.isAssignableFrom(participants.get(0).getClass())) {
        keggEntry1 = parseEntity(((Pathway) participants.get(0)), keggPW, m, species);
      } else if (PhysicalEntity.class.isAssignableFrom(participants.get(0).getClass())) {
        keggEntry1 = parseEntity(
            ((PhysicalEntity) participants.get(0)), keggPW, m,
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
      Set<String> names = inter.getName();
      StringBuffer name = new StringBuffer();
      if (names != null && names.size() > 0) {
        for (String n : names) {
          n = n.trim();
          n = n.replace(" ", "_");
          name.append(n);
          name.append(" ");
        }
        name.delete(name.length() - 2, name.length() - 1);
      }
      graphName = name.toString().trim();

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
   * defines the interaction between two or more genes
   * 
   * @param entity
   * @param keggPW
   * @param m
   * @param species
   */
  private void createKEGGRelationForParticipantList(List<Entity> participants,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species, RelationType relType,
      SubType subType) {
    // interactionType (0 or 1)

    // Participant (2 or more)
    if (participants.size() == 1) {
      createKEGGEntry(participants.get(0), keggPW, m, species, EntryType.gene, null, ",", null);      
    } else if (participants.size() > 1) {
      for (int i = 0; i < participants.size() - 1; i++) {
        for (int j = i + 1; j < participants.size(); j++) {
          EntryExtended keggEntry1 = parsePhysicalEntity((PhysicalEntity)participants.get(i), keggPW, m, species);
          EntryExtended keggEntry2 = parsePhysicalEntity((PhysicalEntity)participants.get(j), keggPW, m, species);
          if (keggEntry1!=null && keggEntry2!=null)
            createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), relType, subType);
        }
      }

      // interactionScore (0 or more) - up to now ignored
      // phenotype (1) - up to now ignored
    }
  }

  /**
   * describes the production of a protein, mRNA, DNA
   * 
   * @param entity
   * @param keggPW
   * @param m
   * @param species
   */
  private Relation createKEGGRelationForTemplateReaction(TemplateReaction entity,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    Relation rel = null;
    for (PhysicalEntity product : entity.getProduct()) {
      // XXX: Not sure if this is completely right to model this as an relation to itself
      EntryExtended keggEntry = parsePhysicalEntity(product, keggPW, m, species);
      
      if (keggEntry !=null)
        rel = createKEGGRelation(keggPW, keggEntry.getId(), keggEntry.getId(), RelationType.GErel,
          new SubType(SubType.EXPRESSION));
    }

    return rel;
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
  private void parseControl(Control entity, de.zbit.kegg.parser.pathway.Pathway keggPW, Model m,
      Species species) {
    if (Catalysis.class.isAssignableFrom(entity.getClass())) {
      createKEGGReactionRelation(((Catalysis) entity).getController(),
          ((Catalysis) entity).getControlled(), getSubtype(((Catalysis) entity).getControlType()),
          keggPW, m, species);
    } else if (TemplateReactionRegulation.class.isAssignableFrom(entity.getClass())) {
      createKEGGReactionRelation(((TemplateReactionRegulation) entity).getController(),
          ((TemplateReactionRegulation) entity).getControlled(),
          getSubtype(((TemplateReactionRegulation) entity).getControlType()), keggPW, m, species);
    } else if (Modulation.class.isAssignableFrom(entity.getClass())) {
      createKEGGReactionRelation(((Modulation) entity).getController(),
          ((Modulation) entity).getControlled(),
          getSubtype(((Modulation) entity).getControlType()), keggPW, m, species);
    } else {
      createKEGGReactionRelation(((Control) entity).getController(),
          ((Control) entity).getControlled(), getSubtype(((Control) entity).getControlType()),
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
    if(cType!=null){
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

  /**
   * deteremines the gene ids of the elements in a pathway
   * 
   * This method is not so clean should be rewritten, becuase in the method
   * {@link BioPaxL32KGML#getEntrezGeneIDsForPathways(List, String, Model)}
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

    Map<String, RelationshipXref> xrefs = getMapFromSet(m.getObjects(RelationshipXref.class));
    for (BioPaxPathwayHolder pw : pathways) {
      // if (pw.getRDFid().equals("http://pid.nci.nih.gov/biopaxpid_9796"))
      // {//TODO: is necessary to uncomment!!!!
      log.log(Level.FINER, "Pathway: " + pw.getPathwayName() + ": " + pw.getNoOfEntities());
      Set<BioPAXElement> pwEntities = new HashSet<BioPAXElement>();
      for (BioPAXElement entity : pw.entities) {
        pwEntities.addAll(getEntitiesWithName((Entity) entity));
        if (!(Pathway.class.isAssignableFrom(entity.getClass())))
          log.log(Level.FINER, "--Input: " + entity.getRDFId() + "\t" + entity.getModelInterface());
      }

      Map<Entity, Integer> geneIDs = getEntrezGeneIDs(pwEntities, species, xrefs);
      for (java.util.Map.Entry<Entity, Integer> entity : geneIDs.entrySet()) {
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
   * @return Collection containing {@link Entity}s having a name and are not
   *         instance of a complex or ComplexAssembly
   */
  protected static Collection<? extends BioPAXElement> getEntitiesWithName(Entity entity) {
    Set<BioPAXElement> resEntities = new HashSet<BioPAXElement>();
    Set<String> name = entity.getName();

    if (name.size() > 0 && !(Pathway.class.isAssignableFrom(entity.getClass()))) {
      if (Complex.class.isAssignableFrom(entity.getClass())) {
        Complex c = (Complex) entity;
        for (PhysicalEntity pe : c.getComponent()) {
          resEntities.addAll(getEntitiesWithName(pe));
        }
      } else if (ComplexAssembly.class.isAssignableFrom(entity.getClass())) {
        ComplexAssembly c = (ComplexAssembly) entity;
        for (Stoichiometry pe : c.getParticipantStoichiometry()) {
          resEntities.addAll(getEntitiesWithName((Entity) pe));
        }

      } else {
        resEntities.add(entity);
      }
    } else if (entity.getParticipantOf().size() > 0
        && !(Pathway.class.isAssignableFrom(entity.getClass()))) {
      for (Entity entity2 : entity.getParticipantOf()) {
        resEntities.addAll(getEntitiesWithName(entity2));
      }
    }
    return resEntities;
  }

  /**
   * This method maps to all gene symbols of the entered entities the
   * corresponding gene id It's an advantage to preprocess the entities to
   * exclusively having entities with a name call therefore the method
   * {@link BioPaxL32KGML#getEntitiesWithName(Entity)} This method also considers
   * 
   * @link {@link Complex} classes. NOTE: the method is not so clean and should
   *       be rewritten!!
   * 
   * @param pwEntities
   * @param species
   * @param xrefs
   * @return
   */
  private static Map<Entity, Integer> getEntrezGeneIDs(Set<BioPAXElement> pwEntities, String species,
      Map<String, RelationshipXref> xrefs) {
    Map<Entity, Integer> geneIDs = new HashMap<Entity, Integer>();

    // searching for gene ids
    Map<String, Entity> geneSymbols = new HashMap<String, Entity>();
    for (BioPAXElement ent : pwEntities) {
      Integer geneID = null;
      Entity entity = (Entity)ent;
      String id = entity.getRDFId();

      // gene id in xrefs?
      if (xrefs.containsKey(id)) {
        geneID = getEntrezGeneIDFromDBxref(xrefs.get(id));
      }

      if (geneID == null) {
        // we have to search the gene id with the gene symbol, adding symbol to
        // the gene symbol set
        Set<String> names = entity.getName();
        if (names != null && names.size() > 0) {
          for (String name : names) {
            if (name.contains("/")) {
              String[] split = name.split("/");
              for (String rs : split) {
                if (!geneSymbols.containsKey(rs))
                  geneSymbols.put(rs, entity);
              }
            } else if (!geneSymbols.containsKey(name)) {
              geneSymbols.put(name, entity);
            }
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
   * Maps the entered gene symbol names to a geneID
   * 
   * @param geneSymbolMapper
   *          {@link GeneSymbol2GeneIDMapper}
   * @param nameSet
   *          with names of one {@link BioPAXElement} element
   * @return the gene id o
   */
  private static Map<Entity, Integer> getEntrezGeneIDOverGeneSymbol(Map<String, Entity> nameSet) {
    log.finest("getGeneIDOverGeneSymbol");
    Map<Entity, Integer> geneIDs = new HashMap<Entity, Integer>();
    Integer geneID = null;

    for (java.util.Map.Entry<String, Entity> symbol : nameSet.entrySet()) {
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
        Map<String, Entity> set = new HashMap<String, Entity>();
        set.put(symbol.getKey().replace("-", ""), symbol.getValue());
        geneIDs.putAll(getEntrezGeneIDOverGeneSymbol(set));
      } else if (symbol.getKey().contains(" ")) {
        log.log(Level.FINER,
            "recall for symbol: " + symbol.getValue().getRDFId() + " " + symbol.getKey());
        Map<String, Entity> set = new HashMap<String, Entity>();
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
  private EntryExtended createKEGGReactionRelation(Set<Controller> controllers,
      Set<org.biopax.paxtools.model.level3.Process> controlleds, SubType subtype,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    EntryExtended keggEntry1 = null;
    RelationType relType = null;

    if (controllers.size() >= 1) {
      for (Controller controller : controllers) {
        if (PhysicalEntity.class.isAssignableFrom(controller.getClass())) {
          keggEntry1 = parsePhysicalEntity((PhysicalEntity) controller, keggPW, m, species);
          relType = RelationType.PPrel;
        } else if (Pathway.class.isAssignableFrom(controller.getClass())) {
          keggEntry1 = createKEGGEntry((Entity) controller, keggPW, m, species, EntryType.map,
              null, ",", null);
          relType = RelationType.maplink;
        } else {
          log.severe("Controller: " + controller.getModelInterface() + "-This should not happen!");
          System.exit(1);
        }

        if (keggEntry1 != null && controlleds.size() > 0) {
          for (org.biopax.paxtools.model.level3.Process process : controlleds) {
            if (Conversion.class.isAssignableFrom(process.getClass())) {
              Conversion con = (Conversion) process;
              if (BiochemicalReaction.class.isAssignableFrom(con.getClass())
                  || ComplexAssembly.class.isAssignableFrom(con.getClass()) 
                  || TransportWithBiochemicalReaction.class.isAssignableFrom(con.getClass())) {
                if (!augmentOriginalKEGGpathway) {
                  Reaction r = null;
                  try {
                    r = createKEGGReaction(((BiochemicalReaction) con).getLeft(),
                        ((BiochemicalReaction) con).getRight(),
                        Utils.iterableToList(((BiochemicalReaction) con)
                            .getParticipantStoichiometry()), keggPW, m, species,
                            getReactionType(con.getConversionDirection()));

                  } catch (ClassCastException e) {
                    try {
                      r = createKEGGReaction(((ComplexAssembly) con).getLeft(),
                          ((ComplexAssembly) con).getRight(),
                          Utils.iterableToList(((ComplexAssembly) con)
                              .getParticipantStoichiometry()), keggPW, m, species,
                              getReactionType(con.getConversionDirection()));
                    } catch (ClassCastException e2) {
                      r = createKEGGReaction(((TransportWithBiochemicalReaction) con).getLeft(),
                          ((TransportWithBiochemicalReaction) con).getRight(),
                          Utils.iterableToList(((TransportWithBiochemicalReaction) con)
                              .getParticipantStoichiometry()), keggPW, m, species,
                              getReactionType(con.getConversionDirection()));
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
              } else if (Transport.class.isAssignableFrom(con.getClass())) {
                List<Relation> rels = createKEGGRelations(((Transport) con).getLeft(),
                    ((Transport) con).getRight(), keggPW, m, species, RelationType.PPrel, 
                    new SubType(SubType.STATE_CHANGE));                
                for (Relation rel : rels) {
                  if (rel !=null)
                    createKEGGRelation(keggPW, keggEntry1.getId(), rel.getEntry2(), relType, subtype);
                }
              } else if (Degradation.class.isAssignableFrom(con.getClass())) {
                List<Relation> rels = createKEGGRelations(((Degradation) con).getLeft(),
                    ((Degradation) con).getRight(), keggPW, m, species, RelationType.PPrel, 
                    new SubType(SubType.STATE_CHANGE));                
                for (Relation rel : rels) {
                  if (rel !=null)
                    createKEGGRelation(keggPW, keggEntry1.getId(), rel.getEntry2(), relType, subtype);
                }
              }  else if (Conversion.class.isAssignableFrom(con.getClass())){
                  List<Relation> rels =  createKEGGRelations(con.getLeft(), con.getRight(), keggPW, m, species, 
                    RelationType.PPrel, getSubtype(con.getInteractionType()));
                  for (Relation rel : rels) {
                    if (rel !=null)
                      createKEGGRelation(keggPW, keggEntry1.getId(), rel.getEntry2(), relType, subtype);
                  }
              } else {
                log.severe("Not programmed case: controlled interface '" + con.getModelInterface()
                    + "'");
                System.exit(1);
              }
            } else if (Pathway.class.isAssignableFrom(process.getClass())) {
              EntryExtended keggEntry2 = createKEGGEntry((Pathway) process, keggPW, m, species,
                  EntryType.map, null, ",", null);
              if (keggEntry2 !=null)
                createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(),
                  relType, subtype);
            } else if (TemplateReaction.class.isAssignableFrom(process.getClass())) {
              Relation rel = createKEGGRelationForTemplateReaction((TemplateReaction) process,
                  keggPW, m, species);
              if (rel != null) {
                createKEGGRelation(keggPW, keggEntry1.getId(), rel.getEntry2(), relType,
                    subtype);
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
   * returns teh appropriate KEGG ReactionType to the ConverstionDirectionType
   * @param conversionDirection
   * @return
   */
  private ReactionType getReactionType(ConversionDirectionType conversionDirection) {
    if (conversionDirection != null){
      switch (conversionDirection) {
        case LEFT_TO_RIGHT:
          return ReactionType.irreversible;
        case RIGHT_TO_LEFT:
          return ReactionType.irreversible;
        case REVERSIBLE:
          return ReactionType.reversible;
        default:
          return ReactionType.other;
      }  
    }
    return ReactionType.other;
  }

  /**
   * parse a BioPax Conversion element
   * 
   * @param entity
   * @param keggPW
   * @param species
   * @param m
   */
  private void parseConversion(Interaction entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species) {
    if (ComplexAssembly.class.isAssignableFrom(entity.getClass())) {
      if  (!augmentOriginalKEGGpathway)
        createKEGGReaction(((ComplexAssembly) entity).getLeft(),
          ((ComplexAssembly) entity).getRight(),
          Utils.iterableToList(((ComplexAssembly) entity).getParticipantStoichiometry()),
          keggPW, m, species, getReactionType(((ComplexAssembly) entity).getConversionDirection()));
    } else if (BiochemicalReaction.class.isAssignableFrom(entity.getClass())) {
      if  (!augmentOriginalKEGGpathway)
        createKEGGReaction(((BiochemicalReaction) entity).getLeft(),
          ((BiochemicalReaction) entity).getRight(),
          Utils.iterableToList(((BiochemicalReaction) entity).getParticipantStoichiometry()),
          keggPW, m, species, getReactionType(((BiochemicalReaction) entity).getConversionDirection()));
    } else if (Degradation.class.isAssignableFrom(entity.getClass())) {
      createKEGGRelations(((Degradation) entity).getLeft(), ((Degradation) entity).getRight(),
          keggPW, m, species, RelationType.PPrel, new SubType(SubType.STATE_CHANGE));
    } else if (Transport.class.isAssignableFrom(entity.getClass())) {
      createKEGGRelations(((Transport) entity).getLeft(), ((Transport) entity).getRight(), keggPW,
          m, species, RelationType.PPrel, new SubType(SubType.STATE_CHANGE));
    } else if (TransportWithBiochemicalReaction.class.isAssignableFrom(entity.getClass())) {
      if  (!augmentOriginalKEGGpathway)
        // BiochemicalReaction br = (TransportWithBiochemicalReaction) entity;
        // deltaG, deltaH, deltaS, ec, and KEQ are ignored
        createKEGGReaction(((BiochemicalReaction) entity).getLeft(),
          ((BiochemicalReaction) entity).getRight(),
          Utils.iterableToList(((BiochemicalReaction) entity).getParticipantStoichiometry()),
          keggPW, m, species, getReactionType(((BiochemicalReaction) entity).getConversionDirection()));
    } else if (Conversion.class.isAssignableFrom(entity.getClass())){
      createKEGGRelations(((Conversion)entity).getLeft(), ((Conversion)entity).getRight(), keggPW, 
          m, species, RelationType.PPrel, getSubtype(((Conversion)entity).getInteractionType()));
    } else {
      log.warning("Unknown kind of Conversion: " + entity.getModelInterface() + "-"
          + entity.getRDFId());
    }
  }

  private SubType getSubtype(Set<InteractionVocabulary> iTypes) {
    for (InteractionVocabulary iType : iTypes) {
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
   * This method is called to create for two PhysicalEntitiy sets relations
   * 
   * @param lefts
   * @param rights
   * @param keggPW
   * @param m
   * @param species
   * @param type
   * @return
   */
  private List<Relation> createKEGGRelations(Set<PhysicalEntity> lefts, Set<PhysicalEntity> rights,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species, RelationType type, 
      SubType subType) {

    List<Relation> relations = new ArrayList<Relation>();

    for (PhysicalEntity left : lefts) {
      EntryExtended keggEntry1 = parsePhysicalEntity(left, keggPW, m, species);
      if (keggEntry1 !=null){
        for (PhysicalEntity right : rights) {
          EntryExtended keggEntry2 = parsePhysicalEntity(right, keggPW, m, species);
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
         
          if (relExists) {
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
   * @param keggPW
   * @param m
   * @param species
   */
  private Reaction createKEGGReaction(Set<PhysicalEntity> lefts, Set<PhysicalEntity> rights,
      List<Stoichiometry> stoichiometry, de.zbit.kegg.parser.pathway.Pathway keggPW, Model m,
      Species species, ReactionType rType) {
    List<ReactionComponent> products = new ArrayList<ReactionComponent>();
    List<ReactionComponent> substrates = new ArrayList<ReactionComponent>();

    for (PhysicalEntity left : lefts) {
      EntryExtended keggEntry = parsePhysicalEntity(left, keggPW, m, species);
      if (keggEntry != null) {
        ReactionComponent rc = new ReactionComponent(keggEntry.getId(), keggEntry.getName());

        for (Stoichiometry stoichi : stoichiometry) {
          if (stoichi.getPhysicalEntity().equals(left)) {
            rc.setStoichiometry((int) stoichi.getStoichiometricCoefficient());
            break;
          }
        }
        substrates.add(rc);  
      }
    }

    for (PhysicalEntity right : rights) {
      EntryExtended keggEntry = parsePhysicalEntity(right, keggPW, m, species);
      if (keggEntry != null) {
        ReactionComponent rc = new ReactionComponent(keggEntry.getId(), keggEntry.getName());

        for (Stoichiometry stoichi : stoichiometry) {
          if (stoichi.getPhysicalEntity().equals(right)) {
            rc.setStoichiometry((int) stoichi.getStoichiometricCoefficient());
            break;
          }
        }
        products.add(rc);  
      }
    }

    Reaction r = null;
    boolean reactionExists = false;
    for (Reaction existingReact : keggPW.getReactions()) {
      List<ReactionComponent> existingProds = existingReact.getProducts();
      List<ReactionComponent> extistingSubs = existingReact.getSubstrates();

      if (existingReact.getType().equals(rType) && existingProds.size() == products.size()
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
