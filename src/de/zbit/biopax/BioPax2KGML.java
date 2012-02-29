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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

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
import org.biopax.paxtools.model.level3.Modulation;
import org.biopax.paxtools.model.level3.MolecularInteraction;
import org.biopax.paxtools.model.level3.Pathway;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.model.level3.Protein;
import org.biopax.paxtools.model.level3.RelationshipXref;
import org.biopax.paxtools.model.level3.Rna;
import org.biopax.paxtools.model.level3.RnaRegion;
import org.biopax.paxtools.model.level3.SmallMolecule;
import org.biopax.paxtools.model.level3.Stoichiometry;
import org.biopax.paxtools.model.level3.TemplateReaction;
import org.biopax.paxtools.model.level3.TemplateReactionRegulation;
import org.biopax.paxtools.model.level3.Transport;
import org.biopax.paxtools.model.level3.TransportWithBiochemicalReaction;

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
import de.zbit.kegg.parser.pathway.ext.GeneType;
import de.zbit.mapper.GeneID2KeggIDMapper;
import de.zbit.mapper.GeneSymbol2GeneIDMapper;
import de.zbit.parser.Species;
import de.zbit.util.Utils;

/**
 * This class works with PaxTools. It is used to fetch information out of a
 * level 3 BioCarta file This file could be downloaded from
 * http://pid.nci.nih.gov/download.shtml
 * 
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class BioPax2KGML extends BioPaxConverter {
  
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
   * This method creates for each pathway in the model a KGML file with the
   * pathway name
   * 
   * @param m
   */
  public void createKGMLsFromModel(Model m) {
    log.info("Creating for each pathway a KGML file.");
    Set<Pathway> pathways = m.getObjects(Pathway.class);

    String oldSpecies = "Homo sapiens", newSpecies = oldSpecies;
    Species species = new Species(oldSpecies, "_HUMAN", "human", "hsa", 9606);
    initalizeMappers(species);

    int i = 0;
    for (Pathway pathway : pathways) {

//      if (++i == 3) {
        // try {
        // writePathwayOwlFile(pathway, m, species);
        // } catch (IOException e) {
        // log.log(Level.WARNING, "File writing was not successful!", e);
        // }
//        System.out.println(i);
        
        // create the pathway
        int number = getKeggPathwayNumber(pathway.getRDFId());
        de.zbit.kegg.parser.pathway.Pathway keggPW = new de.zbit.kegg.parser.pathway.Pathway(
            "biocarta:" + number, getKEGGOrganism(pathway.getOrganism()), number,
            pathway.getStandardName());

        // determine the pathway organism
        BioSource pwOrg = pathway.getOrganism();
        if (pwOrg != null) {          
          if (pwOrg.getName() != null){
            newSpecies = Utils.getListOfCollection(pwOrg.getName()).get(0);
            if (!newSpecies.equals(oldSpecies)){
              species = Species.search(allSpecies, newSpecies, Species.COMMON_NAME);
              initalizeMappers(species);
              oldSpecies = newSpecies;
            }
          }
        }

        addImageLinkToKEGGpathway(species, pathway, keggPW);

        for (org.biopax.paxtools.model.level3.Process interaction : pathway.getPathwayComponent()) {
          parseInteraction((Interaction) interaction, keggPW, m, species);
        }


        KGMLWriter.writeKGML(keggPW);

//        break;
//      }
    }
  }


  /**
   * firstly set {@link BioPax2KGML#augmentOriginalKEGGpathway} to true, 
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
            String pathwaySpecies = Utils.getListOfCollection(pwOrg.getName()).get(0);
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
   * determine the link for the pathway
   * 
   * @param species
   * @param pathway
   * @param keggPW
   */
  public void addImageLinkToKEGGpathway(Species species, Pathway pathway,
      de.zbit.kegg.parser.pathway.Pathway keggPW) {
    String linkName = getNameWithoutBlanks(pathway.getName());
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
   * The methods should be used if the owl file exclusively contains one pathway
   * and the elements in the BioPax file are not explicitly denoted to the
   * pathway (i.e. the file does not contain <bp>Pathway: ....</bp> parts
   * 
   * @param m
   * @param species
   * @param pathwayName
   * @param pathNo
   * @param pathwayStandardName
   */
  public void createKGMLsFromBioPaxFile(Model m, Species species, String pathwayName, int pathNo,
      String pathwayStandardName) {
    initalizeMappers(species);

    de.zbit.kegg.parser.pathway.Pathway keggPW = new de.zbit.kegg.parser.pathway.Pathway(
        pathwayName, species.getKeggAbbr(), pathNo, pathwayStandardName);

    for (Entity entity : m.getObjects(Entity.class)) {
      parseEntity(entity, keggPW, m, species);
    }

    KGMLWriter.writeKGML(keggPW);
  }

  /**
   * returns name of the set without blanks 
   * 
   * @param 
   * @return name
   */
  private String getNameWithoutBlanks(Set<String> names) {
    String link = "";

    if (names != null && names.size() > 0) {
      List<String> names2 = Utils.getListOfCollection(names);
      for (int i = names2.size() - 1; i > 0; i--) {
        link = names2.get(i);
        if (link.length()>0 && !link.contains(" "))
          return link;
      }
    }

    return link;
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
  private Entry parseEntity(Entity entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species) {
    Entry keggEntry = null;
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
  private Entry parsePhysicalEntity(PhysicalEntity entity,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    Entry keggEntry = null;

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
      Entry keggEntry = parsePhysicalEntity(physicalEntity, keggPW, m, species);
      if (keggEntry != null)
        components.add(keggEntry.getId());
    }

    return components;
  }

  /**
   * mapps an entered gene id to a kegg id, if this is not possible the species
   * abbreviation:geneID is returned
   * 
   * @param mapper
   * @return
   */
  private String mapGeneIDToKEGGID(Integer geneID, Species species) {
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
   * Adds the created {@link Entry} to the
   * {@link BioPax2KGML#bc2KeggEntry} map and to the
   * {@link de.zbit.kegg.parser.pathway.Pathway}
   * 
   * @param entity
   * @param keggPW
   * @param mapper
   * @param m
   * @return
   */
  private Entry createKEGGEntry(Entity entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species, EntryType eType, GeneType gType, String graphNameSeparator,
      List<Integer> components) {
    Entry keggEntry;

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
  private void parseInteraction(Interaction entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species) {
    if (Control.class.isAssignableFrom(entity.getClass())) {
      parseControl((Control) entity, keggPW, m, species);
    } else if (Conversion.class.isAssignableFrom(entity.getClass())) {
      parseConversion((Conversion) entity, keggPW, m, species);
    } else if (GeneticInteraction.class.isAssignableFrom(entity.getClass())) {
      createKEGGRelationForParticipantList(
          Utils.getListOfCollection(((GeneticInteraction) entity).getParticipant()), keggPW, m,
          species, RelationType.GErel, new SubType(SubType.ASSOCIATION));
    } else if (MolecularInteraction.class.isAssignableFrom(entity.getClass())) {
      createKEGGRelationForParticipantList(
          Utils.getListOfCollection(((MolecularInteraction) entity).getParticipant()), keggPW, m,
          species, RelationType.PPrel, new SubType(SubType.INDIRECT_EFFECT));
    } else if (TemplateReaction.class.isAssignableFrom(entity.getClass())) {
      createKEGGRelationForTemplateReaction((TemplateReaction) entity, keggPW, m, species);
    } else {
      log.log(Level.SEVERE, "Unknown interaction type: " + entity.getModelInterface() + ".");
      System.exit(1);
    }
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
    for (int i = 0; i <= participants.size() - 1; i++) {
      for (int j = i + 1; j <= participants.size(); j++) {
        Entry keggEntry1 = createKEGGEntry((Gene) participants.get(i), keggPW, m, species,
            EntryType.gene, null, ",", null);
        Entry keggEntry2 = createKEGGEntry((Gene) participants.get(j), keggPW, m, species,
            EntryType.gene, null, ",", null);
        if (keggEntry1!=null && keggEntry2!=null)
          createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), relType, subType);
      }
    }

    // interactionScore (0 or more) - up to now ignored

    // phenotype (1) - up to now ignored
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
      // XXX: Not sure if this is completely right to model this as an relation
      // to itself
      Entry keggEntry = parsePhysicalEntity(product, keggPW, m, species);
      
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
    } else if (TemplateReaction.class.isAssignableFrom(entity.getClass())) {
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
  private Entry createKEGGReactionRelation(Set<Controller> controllers,
      Set<org.biopax.paxtools.model.level3.Process> controlleds, SubType subtype,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    Entry keggEntry1 = null;
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
                  || ComplexAssembly.class.isAssignableFrom(con.getClass())) {
                if (!augmentOriginalKEGGpathway) {
                  Reaction r = null;
                  try {
                    r = createKEGGReaction(((BiochemicalReaction) con).getLeft(),
                        ((BiochemicalReaction) con).getRight(),
                        Utils.getListOfCollection(((BiochemicalReaction) con)
                            .getParticipantStoichiometry()), keggPW, m, species,
                            getReactionType(con.getConversionDirection()));

                  } catch (ClassCastException e) {
                    r = createKEGGReaction(((ComplexAssembly) con).getLeft(),
                        ((ComplexAssembly) con).getRight(),
                        Utils.getListOfCollection(((ComplexAssembly) con)
                            .getParticipantStoichiometry()), keggPW, m, species,
                            getReactionType(con.getConversionDirection()));
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
              } else {
                log.severe("Not programmed case: controlled interface '" + con.getModelInterface()
                    + "'");
                System.exit(1);
              }
            } else if (Pathway.class.isAssignableFrom(process.getClass())) {
              Entry keggEntry2 = createKEGGEntry((Pathway) process, keggPW, m, species,
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
          Utils.getListOfCollection(((ComplexAssembly) entity).getParticipantStoichiometry()),
          keggPW, m, species, getReactionType(((ComplexAssembly) entity).getConversionDirection()));
    } else if (BiochemicalReaction.class.isAssignableFrom(entity.getClass())) {
      if  (!augmentOriginalKEGGpathway)
        createKEGGReaction(((BiochemicalReaction) entity).getLeft(),
          ((BiochemicalReaction) entity).getRight(),
          Utils.getListOfCollection(((BiochemicalReaction) entity).getParticipantStoichiometry()),
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
          Utils.getListOfCollection(((BiochemicalReaction) entity).getParticipantStoichiometry()),
          keggPW, m, species, getReactionType(((BiochemicalReaction) entity).getConversionDirection()));
    } else {
      log.warning("Unknown kind of Conversion: " + entity.getModelInterface() + "-"
          + entity.getRDFId());
    }
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
      Entry keggEntry1 = parsePhysicalEntity(left, keggPW, m, species);
      if (keggEntry1 !=null){
        for (PhysicalEntity right : rights) {
          Entry keggEntry2 = parsePhysicalEntity(right, keggPW, m, species);
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
      Entry keggEntry = parsePhysicalEntity(left, keggPW, m, species);
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
      Entry keggEntry = parsePhysicalEntity(right, keggPW, m, species);
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

  /**
   * 
   * @return the new KEGG reaction name "rn:unknownx", whereas x is set to the
   *         {@link BioPax2KGML#keggReactionID}.
   *         {@link BioPax2KGML#keggReactionID} is augmented after this step
   */
  private String getReactionName() {
    return "rn:unknown" + String.valueOf(++keggReactionID);
  }

  /**
   * 
   * @return the new KEGG unknonw "unknownx", whereas x is set to the
   *         {@link BioPax2KGML#keggUnknownNo}.{@link BioPax2KGML#keggUnknownNo}
   *         is incremented after this step
   */
  private String getKEGGUnkownName() {
    return keggUnknownName + String.valueOf(++keggUnknownNo);
  }

  /**
   * The rdfID is in the format: http://pid.nci.nih.gov/biopaxpid_9717
   * 
   * From this id the number is excluded and used as pathway number, if this is
   * not possible the {@link BioPax2KGML#keggPathwayNumberCounter} is used and
   * incremented
   * 
   * @param rdfId
   * @return
   */
  private int getKeggPathwayNumber(String rdfId) {
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
   * 
   * @return a unique {@link BioPax2KGML#keggEntryID}.
   */
  private int getKeggEntryID() {
    keggEntryID++;
    return keggEntryID;
  }

  /**
   * returns for a given species the KEGG abbreviation, up to now exclusively
   * homo sapiens and mus musculus are encoded. The default value is
   * {@link BioPax2KGML#organism}
   * 
   * @param org
   * @return
   */
  private String getKEGGOrganism(BioSource org) {
    if (org != null && org.getStandardName() != null) {
      if (org.getStandardName().equals("Homo sapiens")) {
        return "hsa";
      } else if (org.getStandardName().equals("Mus musculus")) {
        return "mmu";
      }
    }
    return organism;
  }

  /**
   * The {@link BioPax2KGML}{@link #geneSymbolMapper} and
   * {@link BioPax2KGML#geneIDKEGGmapper} are initalized for the entered species
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
}
