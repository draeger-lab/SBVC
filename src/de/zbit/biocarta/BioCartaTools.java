/*
 * $Id$
 * $URL$
 * ---------------------------------------------------------------------
 * This file is part of the SysBio API library.
 *
 * Copyright (C) 2011 by the University of Tuebingen, Germany.
 *
 * This library is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation. A copy of the license
 * agreement is provided in the file named "LICENSE.txt" included with
 * this software distribution and also available online as
 * <http://www.gnu.org/licenses/lgpl-3.0-standalone.html>.
 * ---------------------------------------------------------------------
 */
package de.zbit.biocarta;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.xml.stream.XMLStreamException;

import org.biopax.paxtools.io.BioPAXIOHandler;
import org.biopax.paxtools.io.SimpleIOHandler;
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
import org.biopax.paxtools.model.level3.UnificationXref;
import org.biopax.paxtools.model.level3.Xref;

import de.zbit.kegg.KGMLWriter;
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
import de.zbit.mapper.GeneID2KeggIDMapper;
import de.zbit.mapper.GeneSymbol2GeneIDMapper;
import de.zbit.parser.Species;
import de.zbit.util.ProgressBar;
import de.zbit.util.StringUtil;
import de.zbit.util.Utils;

/**
 * This class works with PaxTools. It is used to fetch information out of a level 3 BioCarta file
 * This file could be downloaded from http://pid.nci.nih.gov/download.shtml  
 * 
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class BioCartaTools {
  
  static List<Species> allSpecies=null;

  static {
    // TODO: Replace this with a flat-text file in resources to support more organisms.
    allSpecies = new ArrayList<Species>(3);
    allSpecies.add( new Species("Homo sapiens", "_HUMAN", "Human", "hsa", 9606) );
    allSpecies.add( new Species("Mus musculus", "_MOUSE", "Mouse", "mmu", 10090) );
    allSpecies.add( new Species("Rattus norvegicus", "_RAT", "Rat", "rno", 10116) );
  }
    
  /**
   * number which is used to determine a pathway id, if it is not possible
   * to exclude the id from the BioCarta file
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
  String defaultOrganism = "hsa";
  
  /**
   * undefined, if we have no gene id to set the kegg name of an entry we use this name
   */
  String keggUnknownName = "unknown";
  
  /**
   * mapper to map gene symbols to gene ids
   */
  private static GeneSymbol2GeneIDMapper geneSymbolMapper = null;
  
  /**
   * mapper to map gene ids to KEGG ids
   */
  private static GeneID2KeggIDMapper geneIDKEGGmapper = null;
  
  public static final Logger log = Logger.getLogger(BioCartaTools.class
      .getName());

  public BioCartaTools() {
  }

  public Model getModel(String file) {
    InputStream io = null;
    try {
      io = new FileInputStream(new File(file));
    } catch (FileNotFoundException e) {
      log.log(Level.SEVERE, "Could not parse file: " + file + ".", e);
    }
    
    log.log(Level.INFO, "Model sucessfully created");
    return getModel(io);
  }
  
  /**
   * return the parsed file
   * @param io
   * @return
   */
  public Model getModel(InputStream io) {
    BioPAXIOHandler handler = new SimpleIOHandler();
    //TODO: check if level is level 3 or higher, because this class is constructed to work for relations 
    return handler.convertFromOWL(io);
  }
    
  public void createKGMLsFromBioCartaModel(Model m) {
    Set<Pathway> pathways = m.getObjects(Pathway.class);
    Map<String, UnificationXref> xrefs = getMapFromSet(m.getObjects(UnificationXref.class));

    int oldTax = 9606, newTax = 9606;
    Species species = new Species("Homo sapiens", "_HUMAN", "human", "hsa",9606);
    initalizeMappers(species);
    
    int i = 0;
    for (Pathway pathway : pathways) {
            
      if (++i == 3) {
        // try {
        // writePathwayOwlFile(pathway, m, species);
        // } catch (IOException e) {
        // log.log(Level.WARNING, "File writing was not successful!", e);
        // }

        String name = getPathwayName(pathway);
        int number = getKeggPathwayNumber(pathway.getRDFId());
        de.zbit.kegg.parser.pathway.Pathway keggPW = new de.zbit.kegg.parser.pathway.Pathway(
            "biocarta:" + number, getKEGGOrganism(pathway.getOrganism()), number,
            pathway.getStandardName());
        
        BioSource organism = pathway.getOrganism();
        if (organism!=null){
          Set<Xref> s = organism.getXref();
          UnificationXref org = xrefs.get(s);
          if (org != null && org.getDb().equals("TAXONOMY")){
            newTax = Integer.parseInt(org.getId());
          }
          
          if(oldTax != newTax) {
            oldTax = newTax;
            species = Species.search(allSpecies, Integer.toString(newTax), Species.NCBI_TAX_ID);
            initalizeMappers(species);
          }
        }
        

        String linkName = getPathwayLinkName(pathway);
        if (!linkName.equals("") && linkName.contains("pathway")){
          linkName = linkName.replace("pathway", "Pathway");
          
          if (species.getKeggAbbr().equals("hsa")){          
            keggPW.setLink("http://www.biocarta.com/pathfiles/h_" + linkName + ".asp");
            keggPW.setImage("http://www.biocarta.com/pathfiles/h_" + linkName + ".gif");
          } else if (species.getKeggAbbr().equals("mmu")){
            keggPW.setLink("http://www.biocarta.com/pathfiles/m_" + linkName + ".asp");
            keggPW.setImage("http://www.biocarta.com/pathfiles/m_" + linkName + ".gif");
          }
        }
        

        for (org.biopax.paxtools.model.level3.Process interaction : pathway.getPathwayComponent()) {
          createKEGGReactionRelationForInteraction((Interaction) interaction, keggPW, m, species);
        }


        ArrayList<de.zbit.kegg.parser.pathway.Entry> entries = keggPW.getEntries();
        
        int counter = 0;
        int counterMacro = 0;
        int counterMaps = 0;
        for (de.zbit.kegg.parser.pathway.Entry entry : entries) {
//          System.out.println(entry.getName() + " - " + entry.getType());
          if (entry.getType().equals(EntryType.group)) counter++;
          if (entry.getType().equals(EntryType.map)) counterMaps++;
          if (entry.getType().equals(EntryType.compound)) counterMacro++;
        }
        System.out.println("Entries.size(): " + entries.size() + " - " + counter + " are groups, "
            + " " + counterMacro + " are marcromolecules, "
            + " " + counterMaps + " are maps.");
        
        System.out.println("Reactions.size(): " + keggPW.getReactions().size());
        System.out.println("Relations.size(): " + keggPW.getRelations().size());
        
        if (entries.size() > 0) {
          String fileName = name + ".xml";
          fileName = fileName.replace(" ", "_");  
          writeKGML(keggPW, "pws/" + StringUtil.removeAllNonFileSystemCharacters(fileName));
        }
        break;
      }
    }
  }  
  
  /**
   * method to create a pathay file out of the whole owl file
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
        for (Pathway pw : string.getPathwayComponentOf()){
          if(pw.equals(pathway)){
            bw.append(entity.getModelInterface() + "\t" + entity.getRDFId() + "\n");
          }
        }
      }
    }    

    
    bw.close();
  }

  public void createKGMLsFromBioPaxFile(Model m, Species species, String pathwayName, int pathNo, 
      String pathwayStandardName) {
    initalizeMappers(species);
    
    de.zbit.kegg.parser.pathway.Pathway keggPW = new de.zbit.kegg.parser.pathway.Pathway
      (pathwayName, species.getKeggAbbr(), pathNo, pathwayStandardName);
    
    for (Entity entity : m.getObjects(Entity.class)) {
      createKEGGEntryReactionRelation(entity, keggPW, m, species);
    }    

    
    String fileName  = keggPW.getName() + ".xml";
    fileName = fileName.replace(" ", "_");

    writeKGML(keggPW, "pws/" + StringUtil.removeAllNonFileSystemCharacters(fileName)); 
  }  
  
  /**
   * creates out of the pathway set the corresponding name of the pathway, the name with blanks is 
   * prefered
   * @param pathway
   * @return 
   */
  private String getPathwayName(Entity pathway) {
    Set<String> names = pathway.getName();
    String name = "";

    if (names != null && names.size() > 0) {
      List<String> names2 = Utils.getListOfCollection(names);
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
   * returns the name without links 
   * prefered
   * @param pathway
   * @return 
   */
  private String getPathwayLinkName(Entity pathway) {
    Set<String> names = pathway.getName();
    
    String link = "";

    if (names != null && names.size() > 0) {
      List<String> names2 = Utils.getListOfCollection(names);
      for (int i=names2.size()-1; i>0; i--){
        link = names2.get(i); 
        if (!link.contains(" "))
          return link;
      }
    }
    
    return link;
  }

  /**
   * if the fileName is not set it will be set automatically to the pathway name. The file will
   * be saved in the current folder
   * 
   * @param keggPW
   * @param fileName
   */
  private void writeKGML(de.zbit.kegg.parser.pathway.Pathway keggPW, String fileName) {
    if (fileName == null){
      fileName = keggPW.getName();
      if (fileName == null) fileName = Integer.toString(keggPW.hashCode());
      fileName = StringUtil.removeAllNonFileSystemCharacters(fileName);
    }
    if(!fileName.toLowerCase().endsWith(".xml")){
      fileName += ".xml";
    }
    
    try {
      KGMLWriter.writeKGMLFile(keggPW, fileName);
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (XMLStreamException e) {
      log.log(Level.WARNING, "Could not write file for pathway: '" + keggPW.getName() + "'.", e);
    }    
  }
  
  /**
   * 
   * @param entity
   * @param keggPW
   * @param mapper
   * @param m
   * @param species
   * @return 
   * @return
   */
  private EntryExtended createKEGGEntryReactionRelation(Entity entity, de.zbit.kegg.parser.pathway.Pathway keggPW,
      Model m, Species species) {
    EntryExtended keggEntry = null;    
    if (PhysicalEntity.class.isAssignableFrom(entity.getClass())){
      keggEntry = createKEGGEntryForPhysicalEntity((PhysicalEntity)entity, keggPW, m, species);      
    } else if (Interaction.class.isAssignableFrom(entity.getClass())){
      createKEGGReactionRelationForInteraction((Interaction) entity, keggPW, m, species);
    } else if (Pathway.class.isAssignableFrom(entity.getClass())){
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.map, null, ",", null);      
    } else if (Gene.class.isAssignableFrom(entity.getClass())){
      keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, null, ",", null);
    } else {
      log.severe("Unknonw entity type: " + entity.getModelInterface() + "-" + entity.getRDFId());
      System.exit(1);
    }
    
    return keggEntry;
  }

  /**
   * creates a {@link EntryExtended} for a given {@link PhysicalEntity} and adds the pair to the
   * {@link BioCartaTools#bc2KeggEntry} map.
   * The keggEntry is automatically added to the {@link de.zbit.kegg.parser.pathway.Pathway}
   * 
   * @param entity
   * @param keggPW
   */
  private de.zbit.kegg.parser.pathway.ext.EntryExtended createKEGGEntryForPhysicalEntity
  (PhysicalEntity entity, de.zbit.kegg.parser.pathway.Pathway keggPW, 
      Model m, Species species) {
    de.zbit.kegg.parser.pathway.ext.EntryExtended keggEntry = null;
    
    if (Complex.class.isAssignableFrom(entity.getClass())){
      List<Integer> components = createComponentsList(((Complex)entity).getComponent(),
          keggPW, m, species);      
      keggEntry = createKEGGEntry((Entity)entity, keggPW, m, species, EntryType.group, null, 
          "/", components);
    } else if (Dna.class.isAssignableFrom(entity.getClass())) {
        keggEntry = createKEGGEntry(entity, keggPW, m, species, EntryType.gene, GeneType.dna, 
            ",", null);
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
  
  private List<Integer> createComponentsList(Set<PhysicalEntity> complexEntries, 
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    List<Integer> components = new ArrayList<Integer>();
    
    for (PhysicalEntity physicalEntity : complexEntries) {         
      EntryExtended keggEntry = createKEGGEntryForPhysicalEntity(physicalEntity, keggPW, m, species);          
      components.add(keggEntry.getId());        
    }
    
    return components;
  }

  /**
   * mapps an entered gene id to a kegg id, if this is not possible the species abbreviation:geneID is returned
   * @param mapper
   * @return
   */
  private String mapGeneIDToKEGGID(Integer geneID, Species species) {
    String keggName = null;
    try {
      keggName = geneIDKEGGmapper.map(geneID);
    } catch (Exception e) {
      log.log(Level.WARNING, "Could not map geneid: '" + geneID.toString() + "' to a KEGG id, " +
      		"'speciesAbbreviation:geneID will be used instead.", e);
    }
    
    if (keggName==null) {
      keggName = species.getKeggAbbr() + ":" + geneID.toString();
    }
    
    return keggName;
  }

  /**
   * Adds the created {@link EntryExtended} to the {@link BioCartaTools#bc2KeggEntry} map and
   * to the {@link de.zbit.kegg.parser.pathway.Pathway}
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
    de.zbit.kegg.parser.pathway.ext.EntryExtended keggEntry; 
    
    // creating new KEGG entry
    String keggname = null;    
    Integer geneID = getEntrezGeneID(entity, getMapFromSet(m.getObjects(RelationshipXref.class)));
    if(geneID != -1){  
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
      name.delete(name.length()-2, name.length()-1);
    }      
    graphName = name.toString().trim();
    
    Graphics graphics = null;
    if(eType.equals(EntryType.map))
      graphics = Graphics.createGraphicsForPathwayReference(graphName);
    else if (eType.equals(EntryType.compound))
      graphics = Graphics.createGraphicsForCompound(graphName);
    else 
      graphics = new Graphics(graphName);
    
    
    
    if(components != null){
      keggEntry = new EntryExtended(keggPW, getKeggEntryID(), keggname, eType, graphics);
      keggEntry.addComponents(components);
    } else {
      keggEntry = new EntryExtended(keggPW, getKeggEntryID(), keggname, eType, gType, graphics);  
    }
      
    
    
    // checking if entry already exists
    Collection<de.zbit.kegg.parser.pathway.Entry> entries = keggPW.getEntries();    
    if (entries !=null && entries.size()>0) {
      for (de.zbit.kegg.parser.pathway.Entry entry : entries) {
        EntryExtended e = (EntryExtended)entry;
        if(e.equalsWithoutIDNameReactionComparison(keggEntry)){
          if (e.isSetReaction()){
            e.appendReaction(keggEntry.getReactionString());
          }
          keggEntry = e;
          System.out.println("entry already exists");
          return keggEntry;
        } 
      }      
    }
    
    // add entry to pathway
    keggPW.addEntry(keggEntry);
    
    return keggEntry;
  }  
    
  
  public int controls=0;
  /**
   * 
   * @param entity
   * @param keggPW
   * @param mapper
   * @param m
   * @param species
   */
  private void createKEGGReactionRelationForInteraction(Interaction entity, de.zbit.kegg.parser.pathway.Pathway keggPW, 
      Model m, Species species) {
    if (Control.class.isAssignableFrom(entity.getClass())){
      controls++;
      createKEGGReactionForControl((Control)entity, keggPW, m, species);
    } else if (Conversion.class.isAssignableFrom(entity.getClass())){
      createKEGGReactionForConversion((Conversion)entity, keggPW, m, species);     
    } else if (GeneticInteraction.class.isAssignableFrom(entity.getClass())){      
      createKEGGRelationForGeneticInteraction((GeneticInteraction)entity, keggPW, m, species);
    } else if (MolecularInteraction.class.isAssignableFrom(entity.getClass())){
      createKEGGRelationForMolecularInteraction((MolecularInteraction)entity, keggPW, m, species);
    } else if (TemplateReaction.class.isAssignableFrom(entity.getClass())){
      createKEGGRelationForTemplateReaction((TemplateReaction)entity, keggPW, m, species);      
    } else {
      log.log(Level.SEVERE, "Unknown interaction type: " + entity.getModelInterface() + ".");
      System.exit(1);
    }
  }

  /**
   * defines the interaction between two or more genes
   * @param entity
   * @param keggPW
   * @param m
   * @param species
   */
  private void createKEGGRelationForGeneticInteraction(GeneticInteraction entity,
    de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    log.finer("GeneticInteraction: " + entity.getModelInterface() + "-" + entity.getRDFId());
    List<Entity> participants = Utils.getListOfCollection(entity.getParticipant());
    
    // interactionType (0 or 1)
    //    RelationType type = determineRelationType(entity.getInteractionType());
        
    // Participant (2 or more)
    for (int i=0; i<=participants.size()-1; i++){ 
      for (int j=i+1; j<=participants.size(); j++){
        EntryExtended keggEntry1 = createKEGGEntry((Gene)participants.get(i), keggPW, m, species, 
            EntryType.gene, null, ",", null);
        EntryExtended keggEntry2 = createKEGGEntry((Gene)participants.get(j), keggPW, m, species, 
            EntryType.gene, null, ",", null);
        
        createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), RelationType.other, null);
        
        log.config(" 3 " + entity.getRDFId() + "-" + entity.getModelInterface() + " : " 
          + keggEntry1.getId() + "-" + keggEntry1.getName()+ "-" + keggEntry2.getId() + "-" + keggEntry2.getName()
          );
      }      
    }      
    
    // interactionScore (0 or more) - up to now ignored
    
    // phenotype (1) - up to now ignored
  }

  /**
  * describes an interaction involving molecular contact between but the exact mechanism is not known
  * @param entity
  * @param keggPW
  * @param m
  * @param species
  */
 private void createKEGGRelationForMolecularInteraction(MolecularInteraction entity,
     de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
   log.finer("MolecularInteraction: " + entity.getModelInterface() + "-" + entity.getRDFId());
   List<Entity> participants = Utils.getListOfCollection(entity.getParticipant());
       
   // Participant (2 or more)
   for (int i=0; i<=participants.size()-1; i++){ 
     for (int j=i+1; j<=participants.size(); j++){
       EntryExtended keggEntry1 = createKEGGEntry((Gene)participants.get(i), keggPW, m, species, 
           EntryType.gene, null, ",", null);
       EntryExtended keggEntry2 = createKEGGEntry((Gene)participants.get(j), keggPW, m, species, 
           EntryType.gene, null, ",", null);
       
       createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), RelationType.PPrel, null);
       log.config(" 3 " + entity.getRDFId() + "-" + entity.getModelInterface() + " : " 
         + keggEntry1.getId() + "-" + keggEntry1.getName() + "-" + keggEntry2.getId() + "-" + keggEntry2.getName() 
         );
     }      
   }  
   
 }

  /**
   * describes the production of a protein, mRNA, DNA
   * @param entity
   * @param keggPW
   * @param m
   * @param species
   */
  private Relation createKEGGRelationForTemplateReaction(TemplateReaction entity,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {    
    log.finer("TemplateReaction: " + entity.getRDFId() + "-" + entity.getModelInterface());
    Relation rel = null;
    for (PhysicalEntity product : entity.getProduct()) {
      log.info("TemplateReaction: product-" + product.getRDFId() + "-" + product.getModelInterface());
      // XXX: Not sure if this is completely right to model this as an relation to itself
      EntryExtended keggEntry = createKEGGEntryForPhysicalEntity(product, keggPW, m, species);
      rel = createKEGGRelation(keggPW, keggEntry.getId(), keggEntry.getId(), RelationType.GErel, null);
         
      log.config(" 3 " + entity.getRDFId() + "-" + entity.getModelInterface() + " : " 
        + keggEntry.getId() + "-" + keggEntry.getName() + "-" + keggEntry.getId() + "-" + keggEntry.getName() 
        );
    }
    
    return rel; 
  }

  
  /**
   * 
   * @param entity
   * @param keggPW
   * @param mapper
   * @param m
   * @param species
   */
  private void createKEGGReactionForControl(Control entity,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {   
    if(Catalysis.class.isAssignableFrom(entity.getClass())){
      createKEGGReactionForCatalysis((Catalysis)entity, keggPW, m, species);
    } else if (TemplateReaction.class.isAssignableFrom(entity.getClass())){
      createKEGGReactionForTemplateReactionRegulation((TemplateReactionRegulation)entity, keggPW, m, species);
    } else if (Modulation.class.isAssignableFrom(entity.getClass())){
      createKEGGReactionForModulation((Modulation)entity, keggPW, m, species);
    } else {
      createKEGGRelationForControl((Control)entity, keggPW, m, species);
    }
  }

  private EntryExtended createKEGGRelationForControl(Control control,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    EntryExtended keggEntry1 = null, keggEntry2 = null;
    RelationType relType = null;
    Set<Controller> controllers = control.getController();
    Set<org.biopax.paxtools.model.level3.Process> controlleds = control.getControlled();
    // Controller = PhysicalEntity, Pathway (0 or 1)
    // Controlled = Conversion (0 or 1)    

    if(controllers.size()>=1){
      for (Controller controller : control.getController()) {
        if (PhysicalEntity.class.isAssignableFrom(controller.getClass())){
          keggEntry1 = createKEGGEntryForPhysicalEntity((PhysicalEntity)controller, keggPW, m, species);
          relType = RelationType.ECrel;
        } else if (Pathway.class.isAssignableFrom(controller.getClass())){
          keggEntry1 = createKEGGEntry((Entity)controller, keggPW, m, species, EntryType.map,
              null, ",", null);
          relType = RelationType.maplink;
        } else {
          log.severe("Controller: " + controller.getModelInterface() + "-This should not happen!");
          System.exit(1);
        }
        if(keggEntry1 != null && controlleds.size()>0) {
          for (org.biopax.paxtools.model.level3.Process process : controlleds) {
            if (Conversion.class.isAssignableFrom(process.getClass())){
              Conversion con = (Conversion) process;
              if(BiochemicalReaction.class.isAssignableFrom(con.getClass()) && !relType.equals(RelationType.maplink)){
                Reaction r = createKEGGReactionForBiochemicalReaction(((BiochemicalReaction)con).getLeft(),
                    ((BiochemicalReaction)con).getRight(), 
                    Utils.getListOfCollection(((BiochemicalReaction)con).getParticipantStoichiometry()),keggPW, m, species);
                keggEntry1.appendReaction(r.getName()); 
              } else if(Transport.class.isAssignableFrom(con.getClass())){
                List<Relation> rels = createKEGGRelations(((Transport)con).getLeft(),
                    ((Transport)con).getRight(), keggPW, m, species, RelationType.other);
                for (Relation rel : rels) {
                  createKEGGRelation(keggPW, keggEntry1.getId(), rel.getEntry2(), RelationType.other, 
                      (getSubtypes(control.getControlType())));  
                }                
              }
            } else if (Pathway.class.isAssignableFrom(process.getClass())){
              keggEntry2 = createKEGGEntry((Entity)process, keggPW, m, species, EntryType.map,
                  null, ",", null);
              createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), RelationType.other, 
                  getSubtypes(control.getControlType()));              
            } else if (TemplateReaction.class.isAssignableFrom(process.getClass())){
              Relation rel = createKEGGRelationForTemplateReaction((TemplateReaction) process, keggPW, m, species);
              if (rel != null) {
                createKEGGRelation(keggPW, keggEntry1.getId(), rel.getEntry2(), 
                    RelationType.other, getSubtypes(control.getControlType()));
              }              
            } else {
              log.severe("Not programmed case: Catalysis with controlled interface '" + 
                  process.getModelInterface() + "'");
              System.exit(1);
            }
          }
          
          if (keggEntry2 != null) {
            createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), relType, null);
          }
        } else {
          break;
        }
      }
    }
    
    return keggEntry1;
  }
  
  private List<SubType> getSubtypes(ControlType cType) {
    List<SubType> types = new ArrayList<SubType>();
    
    switch (cType) {
      case ACTIVATION:
        types.add(new SubType(SubType.ACTIVATION));
        break;
      case ACTIVATION_ALLOSTERIC:
        types.add(new SubType(SubType.ACTIVATION));
        break;
      case ACTIVATION_NONALLOSTERIC:
        types.add(new SubType(SubType.ACTIVATION));
        break;
      case ACTIVATION_UNKMECH:
        types.add(new SubType(SubType.ACTIVATION));
        break;
      case INHIBITION: 
        types.add(new SubType(SubType.INHIBITION));
        break;
      case INHIBITION_ALLOSTERIC:
        types.add(new SubType(SubType.INHIBITION));
        break;
      case INHIBITION_COMPETITIVE:
        types.add(new SubType(SubType.INHIBITION));
        break;
      case INHIBITION_IRREVERSIBLE:
        types.add(new SubType(SubType.INHIBITION));
        break;
      case INHIBITION_NONCOMPETITIVE:
        types.add(new SubType(SubType.INHIBITION));
        break;
      case INHIBITION_OTHER:
        types.add(new SubType(SubType.INHIBITION));
        break;
      case INHIBITION_UNCOMPETITIVE:
        types.add(new SubType(SubType.INHIBITION));
        break;
      case INHIBITION_UNKMECH:
        types.add(new SubType(SubType.INHIBITION));
        break;
      default:
        types.add(new SubType(SubType.INDIRECT_EFFECT));
    }
    
    return types;
  }

  private EntryExtended createKEGGReactionForCatalysis(Catalysis catalysis,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    EntryExtended keggEntry1 = null;
    RelationType relType = null;
    Set<Controller> controllers = catalysis.getController();
    Set<org.biopax.paxtools.model.level3.Process> controlleds = catalysis.getControlled();
    // Controller = PhysicalEntity, Pathway (0 or 1)
    // Controlled = Conversion (0 or 1)    

    if(controllers.size()>=1){
      for (Controller controller : catalysis.getController()) {
        if (PhysicalEntity.class.isAssignableFrom(controller.getClass())){
          keggEntry1 = createKEGGEntryForPhysicalEntity((PhysicalEntity)controller, keggPW, m, species);
          relType = RelationType.ECrel;
        } else if (Pathway.class.isAssignableFrom(controller.getClass())) {
          keggEntry1 = createKEGGEntry((Entity)controller, keggPW, m, species, EntryType.map,
              null, ",", null);
          relType = RelationType.maplink;
        } else {
          log.severe("Controller: " + controller.getModelInterface() + "-This should not happen!");
          System.exit(1);
        }
        
        if(keggEntry1 != null && controlleds.size()>0) {
          for (org.biopax.paxtools.model.level3.Process process : controlleds) {
            if (Conversion.class.isAssignableFrom(process.getClass())){
              Conversion con = (Conversion) process;
              if(BiochemicalReaction.class.isAssignableFrom(con.getClass()) ||
                  ComplexAssembly.class.isAssignableFrom(con.getClass())){
                Reaction r = null;
                try{
                  r = createKEGGReactionForBiochemicalReaction(((BiochemicalReaction)con).getLeft(),
                      ((BiochemicalReaction)con).getRight(), 
                      Utils.getListOfCollection(((BiochemicalReaction)con).getParticipantStoichiometry()),keggPW, m, species);
                    
                } catch (ClassCastException e){
                  r = createKEGGReactionForBiochemicalReaction(((ComplexAssembly)con).getLeft(),
                      ((ComplexAssembly)con).getRight(), 
                      Utils.getListOfCollection(((ComplexAssembly)con).getParticipantStoichiometry()),keggPW, m, species);
                  
                }
                
                if(relType.equals(RelationType.maplink)){
                  for (ReactionComponent rc : r.getSubstrates()){
                    createKEGGRelation(keggPW, keggEntry1.getId(), rc.getId(), relType, getSubtypes(catalysis.getControlType()));
                  }
                } else if (relType.equals(RelationType.ECrel)){
                  keggEntry1.appendReaction(r.getName());  
                }                 
                log.config(" 3 controller: " + keggEntry1.getId() + "- Controlled: " + r.getName());
              } else if (Transport.class.isAssignableFrom(con.getClass())){
                List<Relation> rels = createKEGGRelations(((Transport)con).getLeft(), ((Transport)con).getRight(), 
                    keggPW, m, species, RelationType.other);
                for (Relation rel : rels) {
                  createKEGGRelation(keggPW, keggEntry1.getId(), rel.getEntry2(), relType, getSubtypes(catalysis.getControlType()));  
                }                                
              } else {
                log.severe("Not programmed case: Catalysis with controlled interface '" + 
                    con.getModelInterface() + "'");
                System.exit(1);
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
  
  private void createKEGGReactionForModulation(Modulation modulation,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    Set<Controller> controllers = modulation.getController();
    Set<org.biopax.paxtools.model.level3.Process> controlleds = modulation.getControlled();    
    // Controller = physical entity (0 or 1)
    // Controlled = Catalysis (0 or 1)
    EntryExtended keggEntry1 = null, keggEntry2 = null;
    if (controllers.size()>0){
      for (Controller entity : controllers) {
        if (PhysicalEntity.class.isAssignableFrom(entity.getClass())){
          if (controlleds.size()>0) {
            for (org.biopax.paxtools.model.level3.Process process : controlleds) {
              if (Catalysis.class.isAssignableFrom(process.getClass())){
                keggEntry1 = createKEGGEntryForPhysicalEntity((PhysicalEntity) entity, keggPW, m, species);
                keggEntry2 = createKEGGReactionForCatalysis((Catalysis) process, keggPW, m, species);
                createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), RelationType.other, null);
              } else {
                log.severe("Process: " + process.getModelInterface() + "This should not happen!");
                System.exit(1);
              }
            } 
          }
        } else {
          log.severe("Controller: " + entity.getModelInterface() + "-This should not happen!");
          System.exit(1);
        }
      }
    }
    log.config(" 3 " + modulation.getRDFId() + "-" + modulation.getModelInterface() + " : " 
        + (keggEntry1!=null ? (keggEntry1.getId() + "-" + keggEntry1.getName()) : "") 
        + (keggEntry2!=null ? (keggEntry2.getId() + "-" + keggEntry2.getName()) : ""));
  }

  private void createKEGGReactionForTemplateReactionRegulation(TemplateReactionRegulation tmpReaction,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    Set<Controller> controllers = tmpReaction.getController();
    Set<org.biopax.paxtools.model.level3.Process> controlleds = tmpReaction.getControlled();
    // Controller = physical entity (0 or 1)
    // Controlled = TemplateReaction (0 or 1)
    
    if (controllers.size()>0){
      for (Controller entity : controllers) {
        if (PhysicalEntity.class.isAssignableFrom(entity.getClass())){
          if (controlleds.size()>0) {
            for (org.biopax.paxtools.model.level3.Process tempReaction : controlleds) {
              if (TemplateReaction.class.isAssignableFrom(tempReaction.getClass())){ 
                EntryExtended keggEntry1 = createKEGGEntryForPhysicalEntity((PhysicalEntity)entity, keggPW, m, species);
                Relation r = createKEGGRelationForTemplateReaction((TemplateReaction)tempReaction, keggPW, m, species);
                if (r!=null)
                  createKEGGRelation(keggPW, keggEntry1.getId(), r.getEntry2(), RelationType.other, null); 
              } else {
                log.severe("Process: " + tempReaction.getModelInterface() + "This should not happen!");
                System.exit(1);
              }
            } 
          }
        } else {
          log.severe("Controller: " + entity.getModelInterface() + "-This should not happen!");
          System.exit(1);
        }
      }
    }    
  }

  /**
   * 
   * @param entity
   * @param keggPW
   * @param species
   * @param m
   */
  private void createKEGGReactionForConversion(Interaction entity,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    
    // left  = Physical Entity (0 or more)
    
    // right = Physical Entity (0 or more)
    
    // participantStoichometry = = Physical Entity (0 or more)
    
    // spontaneous     
    if (ComplexAssembly.class.isAssignableFrom(entity.getClass())){
      createKEGGReactionForBiochemicalReaction(((ComplexAssembly)entity).getLeft(),
          ((ComplexAssembly)entity).getRight(), 
          Utils.getListOfCollection(((ComplexAssembly)entity).getParticipantStoichiometry()), keggPW, m, species);
    } else if (BiochemicalReaction.class.isAssignableFrom(entity.getClass())){
      createKEGGReactionForBiochemicalReaction(((BiochemicalReaction)entity).getLeft(),
          ((BiochemicalReaction)entity).getRight(), 
          Utils.getListOfCollection(((BiochemicalReaction)entity).getParticipantStoichiometry()), keggPW, m, species);
    } else if (Degradation.class.isAssignableFrom(entity.getClass())){
      createKEGGRelations(((Transport)entity).getLeft(), ((Transport)entity).getRight(), 
          keggPW, m, species, RelationType.other);  
    } else if (Transport.class.isAssignableFrom(entity.getClass())){
      createKEGGRelations(((Transport)entity).getLeft(), ((Transport)entity).getRight(), 
          keggPW, m, species, RelationType.other);
    } else if (TransportWithBiochemicalReaction.class.isAssignableFrom(entity.getClass())){
      BiochemicalReaction br = (TransportWithBiochemicalReaction) entity;
      log.info("deltaG: " + br.getDeltaG() + ", deltaH: " + br.getDeltaH()
          + ", deltaS: " + br.getDeltaS() + ", ec: " + br.getECNumber()
          + ", KEQ: " + br.getKEQ());
      createKEGGReactionForBiochemicalReaction(((BiochemicalReaction)entity).getLeft(),
          ((BiochemicalReaction)entity).getRight(), 
          Utils.getListOfCollection(((BiochemicalReaction)entity).getParticipantStoichiometry()), keggPW, m, species);
    }  else {
      log.warning("Unknown kind of Conversion: " + entity.getModelInterface() + "-"  + entity.getRDFId());
    }   
  }

 

  private List<Relation> createKEGGRelations(Set<PhysicalEntity> lefts, Set<PhysicalEntity> rights,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m,
      Species species, RelationType type) {

    List<Relation> relations = new ArrayList<Relation>();

    for (PhysicalEntity left : lefts) {
      EntryExtended keggEntry1 = createKEGGEntryForPhysicalEntity(left, keggPW, m, species);   
      for (PhysicalEntity right : rights) {
        EntryExtended keggEntry2 = createKEGGEntryForPhysicalEntity(right, keggPW, m, species);   
        Relation r = createKEGGRelation(keggPW, keggEntry1.getId(), keggEntry2.getId(), type, null);
        if (!relations.contains(r))
          relations.add(r);
      }  
    }
    
    return relations;
  }

  private Relation createKEGGRelation(de.zbit.kegg.parser.pathway.Pathway keggPW, 
      int keggEntry1Id, int keggEntry2Id, RelationType type, List<SubType> subTypes) {
    ArrayList<Relation> existingRels = keggPW.getRelations();
    Relation r = null; 
    
    // Check if it already exists and only create novel relations.
    if(existingRels.size()>0){
      for (Relation rel : existingRels) {
        boolean relExists = true;
        if((rel.getEntry1() == keggEntry1Id &&
            rel.getEntry2() == keggEntry2Id)){
          
          relExists &= rel.isSetType()== (type!=null);
          if(relExists && type!=null)
            relExists &= (rel.getType().equals(type));
            
          relExists &= rel.isSetSubTypes()==(subTypes!=null);
          if(relExists && (subTypes!=null)) 
            relExists &= (rel.getSubtypes().equals(subTypes));
          
          if (relExists) {
            r = rel;
            return r;
          }
        }
      }
      
      r = new Relation(keggEntry1Id, keggEntry2Id, type);
    } else {
      r = new Relation(keggEntry1Id, keggEntry2Id, type);
    }
    
    // Add the relation to the pathway
    keggPW.addRelation(r);
    
    return r;
  }

  /**
   * 
   * @param entity
   * @param keggPW
   * @param m
   * @param species
   */
  private Reaction createKEGGReactionForBiochemicalReaction(Set<PhysicalEntity> lefts, Set<PhysicalEntity> rights, 
      List<Stoichiometry> stoichiometry, de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {   
    List<ReactionComponent> products = new ArrayList<ReactionComponent>();
    List<ReactionComponent> substrates = new ArrayList<ReactionComponent>();
    ReactionType rType = ReactionType.other;
    // TODO: Reaction type must be set!!! How to distinguish between irreversible or reversible

    for (PhysicalEntity left : lefts) {
      EntryExtended keggEntry = createKEGGEntryForPhysicalEntity(left, keggPW, m, species);   
      ReactionComponent rc = new ReactionComponent(keggEntry.getId(), keggEntry.getName());
      
      for (Stoichiometry stoichi : stoichiometry) {
        if (stoichi.getPhysicalEntity().equals(left)) {
          rc.setStoichiometry((int) stoichi.getStoichiometricCoefficient());
          break;
        }
      }
      substrates.add(rc);
    }

    for (PhysicalEntity right : rights) {
      EntryExtended keggEntry = createKEGGEntryForPhysicalEntity(right, keggPW, m, species);   
      ReactionComponent rc = new ReactionComponent(keggEntry.getId(), keggEntry.getName());
      
      for (Stoichiometry stoichi : stoichiometry) {
        if (stoichi.getPhysicalEntity().equals(right)) {
          rc.setStoichiometry((int) stoichi.getStoichiometricCoefficient());
          break;
        }
      }
      products.add(rc);
    }
     
    Reaction r = null;          
    boolean reactionExists = false;
    for (Reaction existingReact : keggPW.getReactions()) {
      List<ReactionComponent> existingProds = existingReact.getProducts();
      List<ReactionComponent> extistingSubs = existingReact.getSubstrates();
      
      if (existingReact.getType().equals(rType) &&
          existingProds.size() == products.size()
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
    
      
    if(!reactionExists){
      r = new Reaction(keggPW, getReactionName(), ReactionType.other);
      r.addProducts(products);
      r.addSubstrates(substrates);
      keggPW.addReaction(r);
    }
    
    for (ReactionComponent reactionComponent : products) {
      log.log(Level.CONFIG, " 3 " + r.getName() + " Left: " + reactionComponent.getName() + "(" +reactionComponent.getId() + ")");  
    }
    for (ReactionComponent reactionComponent : substrates) {
      log.log(Level.CONFIG, " 3 " + r.getName() + " Right: " +reactionComponent.getName() + "(" +reactionComponent.getId() + ")");  
    }
    
    
    
    return r;
  }

  /**
   * 
   * @return the new KEGG reaction name "rn:unknownx", whereas x is set to the 
   * {@link BioCartaTools#keggReactionID}.{@link BioCartaTools#keggReactionID} is augmented after
   * this step 
   */
  private String getReactionName() {    
    return "rn:unknown" + String.valueOf(++keggReactionID);
  }
  
  /**
   * 
   * @return the new KEGG reaction name "rn:unknownx", whereas x is set to the 
   * {@link BioCartaTools#keggReactionID}.{@link BioCartaTools#keggReactionID} is augmented after
   * this step 
   */
  private String getKEGGUnkownName() {    
    return keggUnknownName + String.valueOf(++keggUnknownNo);
  }

  /**
   * The rdfID is in the format: http://pid.nci.nih.gov/biopaxpid_9717 
   *
   * From this id the number is excluded and used as pathway number, if this is
   * not possible the {@link BioCartaTools#keggPathwayNumberCounter} is used and
   * incremented
   * 
   * @param rdfId
   * @return
   */
  private int getKeggPathwayNumber(String rdfId) {
    int posUnderscore = rdfId.indexOf('_');
    if(posUnderscore>-1 && posUnderscore<=rdfId.length()){
      try{
        return Integer.parseInt(rdfId.substring(posUnderscore+1));
      } catch (Exception e) {
        return keggPathwayNumberCounter++;
      }
    }
 
    return keggPathwayNumberCounter++;
  }

  /**
   * The rdfID is in the format: http://pid.nci.nih.gov/biopaxpid_9717 
   *
   * From this id the number is excluded and used as pathway number
   * @param rdfId
   * @return
   */
  private int getKeggEntryID() {
    keggEntryID++;
     return keggEntryID;
  }
  
  /**
   * returns for a given species the KEGG abbreviation, up to now exclusively 
   * homo sapiens and mus musculus are encoded. The default value is {@link BioCartaTools#defaultOrganism}
   * @param organism
   * @return
   */
  private String getKEGGOrganism(BioSource organism) {
    if(organism != null && organism.getStandardName()!=null){
      if(organism.getStandardName().equals("Homo sapiens")){
        return "hsa";
      } else if (organism.getStandardName().equals("Mus musculus")){
        return "mmu";
      }  
    }    
    return defaultOrganism;  
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
 * This method maps to all gene symbols of the entered entities the corresponding gene id
 * It's an advantage to preprocess the entities to exclusively having entities with a name
 * call therefore the method {@link BioCartaTools#getEntitiesWithName(Entity)}  This method also considers 
 * @link {@link Complex} classes but the method is not so clean, should be rewritten!!
 *    
 * @param entities 
 * @param species
 * @param xrefs
 * @return
 */
  private static Map<Entity, Integer> getEntrezGeneIDs(Set<Entity> entities, String species, Map<String, RelationshipXref> xrefs) {
    Map<Entity, Integer> geneIDs = new HashMap<Entity, Integer>();    
    
    //searching for gene ids
    Map<String, Entity> geneSymbols = new HashMap<String, Entity>();
    for (Entity entity : entities) {
      Integer geneID = null;
      String id = entity.getRDFId();
      
      // gene id in xrefs?
      if(xrefs.containsKey(id)) {            
         geneID = getEntrezGeneIDFromDBxref(xrefs.get(id));           
      } 
      
      if (geneID == null) {
        // we have to search the gene id with the gene symbol, adding symbol to the gene symbol set
        Set<String> names = entity.getName();
        if (names != null && names.size()>0) {
          for (String name : names) {
            if(name.contains("/")) {
              String[] split = name.split("/");
              for (String rs : split) {
                if(!geneSymbols.containsKey(rs))
                  geneSymbols.put(rs, entity);
              }
            } else if (!geneSymbols.containsKey(name)){
              geneSymbols.put(name, entity);
            }          
          }     
        }           
      } else if(!geneIDs.containsKey(geneID)) {
        // gene id found directly
        log.log(Level.FINER, "found: " + entity.getRDFId() + " " + geneID);
        geneIDs.put(entity, geneID);
      }
    }
    
    // getting the gene ids which could not be found directly
    if(geneSymbols.size()>0){      
      geneIDs.putAll(getEntrezGeneIDOverGeneSymbol(geneSymbols));
    } 
    
    return geneIDs;
  }

  /**
   * This method maps the gene id to the gene symbol of the entered entity. This method is not able to 
   * map complexes call for that the method {@link BioCartaTools#getEntrezGeneIDsOfComplex(Complex, Map)}
   * It's an advantage to preprocess the entities to exclusively having entities with a name
   * call therefore the method {@link BioCartaTools#getEntitiesWithName(Entity)}  
   * 
   * The default for a gene id which is not found is -1
   *    
   * @param entity 
   * @param xrefs
   * @return
   */
  private static int getEntrezGeneID(Entity entity, Map<String, RelationshipXref> xrefs) {
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
  
  public void initalizeMappers(Species species){    
    try {
      geneSymbolMapper = new GeneSymbol2GeneIDMapper(species.getCommonName());
    } catch (IOException e) {
      log.log(Level.SEVERE, "Could not initalize mapper!", e);
    }
    
    try {
      geneIDKEGGmapper = new GeneID2KeggIDMapper(species);
    } catch (IOException e) {
      log.log(Level.SEVERE, "Error while initializing gene id to KEGG ID mapper.", e);
    }
  }
  /**
   * Maps the entered gene symbol names to a geneID
   * 
   * @param geneSymbolMapper {@link GeneSymbol2GeneIDMapper} 
   * @param nameSet with names of one {@link BioPAXElement} element 
   * @return the gene id o
   */
  private static Map<Entity, Integer> getEntrezGeneIDOverGeneSymbol(Map<String, Entity> nameSet) {
    log.finest("getGeneIDOverGeneSymbol");
    Map<Entity, Integer> geneIDs = new HashMap<Entity, Integer>();
    Integer geneID = null;

      
    for (Entry<String, Entity> symbol : nameSet.entrySet()) {
      try {
        geneID = geneSymbolMapper.map(symbol.getKey()); 
      } catch (Exception e) {
        log.log(Level.WARNING, "Error while mapping name: " + symbol.getKey() + ".", e);
      }
      
      if(geneID!=null){
        log.log(Level.FINER, "----- found! Geneid: " + geneID + " "  + symbol.getValue().getRDFId() + " "+ symbol.getKey());
        geneIDs.put(symbol.getValue(), geneID);
        break;
      } else if (symbol.getKey().contains("-")) {
        log.log(Level.FINER, "recall for symbol: " + symbol.getValue().getRDFId() + " " + symbol.getKey());
        Map<String, Entity> set = new HashMap<String, Entity>();
        set.put(symbol.getKey().replace("-", ""), symbol.getValue());        
        geneIDs.putAll(getEntrezGeneIDOverGeneSymbol(set));
      } else if (symbol.getKey().contains(" ")) {
        log.log(Level.FINER, "recall for symbol: " + symbol.getValue().getRDFId() + " " + symbol.getKey());
        Map<String, Entity> set = new HashMap<String, Entity>();
        set.put(symbol.getKey().replace(" ", "_"), symbol.getValue());      
        geneIDs.putAll(getEntrezGeneIDOverGeneSymbol(set));
      }
      else {
        log.log(Level.FINER, "----- not found " + symbol.getValue().getRDFId() + " " + symbol.getKey());
      }
    }
    
    return geneIDs;
  }


   /**
    * transforms a set to a map. The key is a RDFId and the value the corresponding object
   * @param <T>
   * @param set to convert
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
   * The method returns the smallest entity having a name, i.e. a gene symbol, which could 
   * be parsed
   * 
   * @param entity
   * @return Collection containing {@link Entity}s having a name and are not instance of
   * a complex or ComplexAssembly
   */
  private static Collection<? extends Entity> getEntitiesWithName(Entity entity) {
    Set<Entity> resEntities = new HashSet<Entity>();
    Set<String> name = entity.getName();

    if(name.size()>0 && !(Pathway.class.isAssignableFrom(entity.getClass()))){
      if(Complex.class.isAssignableFrom(entity.getClass())){
        Complex c = (Complex)entity;
        for (PhysicalEntity pe : c.getComponent()) {
          resEntities.addAll(getEntitiesWithName(pe));
        }
      } else if(ComplexAssembly.class.isAssignableFrom(entity.getClass())){
        ComplexAssembly c = (ComplexAssembly)entity;
        for (Stoichiometry pe : c.getParticipantStoichiometry()) {
          resEntities.addAll(getEntitiesWithName((Entity) pe));
        }
        
      } else {
        resEntities.add(entity);
      }
    } else if (entity.getParticipantOf().size()>0 && !(Pathway.class.isAssignableFrom(entity.getClass()))){      
      for (Entity entity2 : entity.getParticipantOf()) {
        resEntities.addAll(getEntitiesWithName(entity2));
      }      
    } 
    return resEntities;
  }
  
  /**
   * This method returns a list of all biocarta pathways with the containing gene IDs
   * 
   * @param species
   * @param model
   * @return
   */
  public List<BioCartaPathwayHolder> getPathwaysWithEntrezGeneID(String species, Model m){     
    return getEntrezGeneIDsForPathways(getPathwayEntities(species, m), species, m);
  }

  /**
   * creates a list with {@link BioCartaPathwayHolder} and the entities of the pathways
   * @param species
   * @param m
   * @return
   */
  private List<BioCartaPathwayHolder> getPathwayEntities(String species, Model m) {
    List<BioCartaPathwayHolder> pathways = new ArrayList<BioCartaPathwayHolder>();
    
    for (Entity entity : m.getObjects(Entity.class)) {
      for (Interaction string : entity.getParticipantOf()) {
        for (Pathway pw : string.getPathwayComponentOf()){
          BioCartaPathwayHolder helper = new BioCartaPathwayHolder(pw.getRDFId(), getPathwayName(pw));
          int index = pathways.indexOf(helper);
          if(index >-1){
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
   * deteremines the gene ids of the elements in a pathway
   * 
   * This method is not so clean should be rewritten, becuase in the method
   * {@link BioCartaTools#getEntrezGeneIDsForPathways(List, String, Model)}
   * complexes are not treated right
   * 
   * @param pathways
   * @param species
   * @param m
   * @return
   */
  public List<BioCartaPathwayHolder> getEntrezGeneIDsForPathways(List<BioCartaPathwayHolder> pathways, String species, Model m){
    log.info("Start parsing gene ids.");
    ProgressBar bar = new ProgressBar(pathways.size());
    
    Map<String, RelationshipXref> xrefs = getMapFromSet(m.getObjects(RelationshipXref.class));
    for (BioCartaPathwayHolder pw : pathways) {
//      if (pw.getRDFid().equals("http://pid.nci.nih.gov/biopaxpid_9796")) {//TODO: is necessary to uncomment!!!!
        log.log(Level.FINER, "Pathway: " + pw.getPathwayName() + ": " + pw.getNoOfEntities());
        Set<Entity> pwEntities = new HashSet<Entity>();
        for (Entity entity : pw.entities) {
          pwEntities.addAll(getEntitiesWithName(entity));
          if(!(Pathway.class.isAssignableFrom(entity.getClass())))//entity instanceof Pathway || entity instanceof PathwayImpl))
            log.log(Level.FINER, "--Input: " + entity.getRDFId() + "\t" + entity.getModelInterface());
        }

        Map<Entity, Integer> geneIDs = getEntrezGeneIDs(pwEntities, species, xrefs);
        for (Entry<Entity, Integer> entity : geneIDs.entrySet()) {
          log.log(Level.FINER, "----res: " + entity.getKey() + " " + entity.getValue());
          pw.addGeneID(entity.getValue());
        }        
//      }//TODO: is necessary to uncomment!!!!
        bar.DisplayBar();
    }
    
    return pathways;
  }
  
  /**
   * This method returns a list of all biocarta pathways with the containing gene IDs
   * 
   * @param species
   * @param model
   * @return
   */
  public List<BioCartaPathwayHolder> getPathwaysWithGeneID(String species, Model m){    
    List<BioCartaPathwayHolder> pathways = new ArrayList<BioCartaPathwayHolder>();
    
    for (Entity entity : m.getObjects(Entity.class)) {
      for (Interaction string : entity.getParticipantOf()) {
        for (Pathway pw : string.getPathwayComponentOf()){
          BioCartaPathwayHolder helper = new BioCartaPathwayHolder(pw.getRDFId());
          int index = pathways.indexOf(helper);
          if(index >-1){
            helper = pathways.get(index);
          } else {
            pathways.add(helper);
          }
          helper.addEntity(entity); 
        }
      }
    }
    
    Map<String, RelationshipXref> xrefs = getMapFromSet(m.getObjects(RelationshipXref.class));
    for (BioCartaPathwayHolder pw : pathways) {
      if (pw.getPathwayName().equals("http://pid.nci.nih.gov/biopaxpid_9796")) {
        log.log(Level.FINER, "Pathway: " + pw.getPathwayName() + ": " + pw.getNoOfEntities());
        Set<Entity> pwEntities = new HashSet<Entity>();
        for (Entity entity : pw.entities) {
          pwEntities.addAll(getEntitiesWithName(entity));
          if(!(Pathway.class.isAssignableFrom(entity.getClass())))
            log.log(Level.FINER, "--Input: " + entity.getRDFId() + "\t" + entity.getModelInterface());
        }

        Map<Entity, Integer> geneIDs = getEntrezGeneIDs(pwEntities, species, xrefs);
        for (Entry<Entity, Integer> entity : geneIDs.entrySet()) {
          log.log(Level.FINER, "----res: " + entity.getValue() + " " + entity.getValue());
          pw.addGeneID(entity.getValue());
        }        
      }
    }
    
    return pathways;
  }
}
