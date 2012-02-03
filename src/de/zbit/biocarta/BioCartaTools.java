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
import org.biopax.paxtools.model.level3.CatalysisDirectionType;
import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.model.level3.ComplexAssembly;
import org.biopax.paxtools.model.level3.Control;
import org.biopax.paxtools.model.level3.Controller;
import org.biopax.paxtools.model.level3.Conversion;
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
import de.zbit.kegg.parser.pathway.EntryType;
import de.zbit.kegg.parser.pathway.Graphics;
import de.zbit.kegg.parser.pathway.Reaction;
import de.zbit.kegg.parser.pathway.ReactionComponent;
import de.zbit.kegg.parser.pathway.ReactionType;
import de.zbit.kegg.parser.pathway.Relation;
import de.zbit.kegg.parser.pathway.RelationType;
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
  String keggUndefinedName = "undefined";
  
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

  public void createKGMLsFromBioCartaModel(Model m, Species species) {
    initalizeMappers(species);
    
    Set<Pathway> pathways = m.getObjects(Pathway.class);
    int i = 0;
    for (Pathway pathway : pathways) {
      if (++i == 3) {
        // try {
        // writePathwayOwlFile(pathway, m, species);
        // } catch (IOException e) {
        // log.log(Level.WARNING, "File writing was not successful!", e);
        // }

        String name = getName(pathway);
        int number = getKeggPathwayNumber(pathway.getRDFId());
        de.zbit.kegg.parser.pathway.Pathway keggPW = new de.zbit.kegg.parser.pathway.Pathway(
            name, getKEGGOrganism(pathway.getOrganism()), number,
            pathway.getStandardName());

        for (org.biopax.paxtools.model.level3.Process interaction : pathway.getPathwayComponent()) {
          createKEGGReactionRelationForInteraction((Interaction) interaction, keggPW, m, species);
        }

        String fileName = keggPW.getName() + ".xml";
        fileName = fileName.replace(" ", "_");

        writeKGML(keggPW, "pws/" + StringUtil.removeAllNonFileSystemCharacters(fileName));

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
    Boolean convertingSuccessfull = false;
    String name = getName(pathway);
    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(name + ".owl")));
    
    int number  = getKeggPathwayNumber(pathway.getRDFId());
    de.zbit.kegg.parser.pathway.Pathway keggPW = new de.zbit.kegg.parser.pathway.Pathway
      (name, getKEGGOrganism(pathway.getOrganism()), number, pathway.getStandardName());
    EntryExtended keggEntry = null;
    
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
   * creates out of the pathway set the corresponding name of the pathway
   * @param pathway
   * @return 
   */
  protected String getName(Entity pathway) {
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
   * calls @link {@link BioCartaTools#writeKGML(de.zbit.kegg.parser.pathway.Pathway, String)
   * @param keggPW
   */
  private void writeKGML(de.zbit.kegg.parser.pathway.Pathway keggPW) {
    writeKGML(keggPW, null);     
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
//    log.config(" 1 " + entity.getRDFId() + "-" + entity.getModelInterface() + " : " + keggEntry.getId() + "-" + keggEntry.getName());
    if (PhysicalEntity.class.isAssignableFrom(entity.getClass())){
      keggEntry = createKEGGEntryForPhysicalEntity((PhysicalEntity)entity, keggPW, m, species);      
    } else if (Interaction.class.isAssignableFrom(entity.getClass())){
      createKEGGReactionRelationForInteraction((Interaction) entity, keggPW, m, species);
    } else if (Pathway.class.isAssignableFrom(entity.getClass())){
      keggEntry = createKEGGEntryForPathway((Pathway)entity, keggPW, m, species);      
    } else if (Gene.class.isAssignableFrom(entity.getClass())){
      keggEntry = createKEGGEntryForGene((Gene)entity, keggPW, m, species);
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
    
    String graphName = null;
    String keggname = null;
    GeneType gType = null;
    EntryType eType = null;
    if (Complex.class.isAssignableFrom(entity.getClass())){
      log.finer("Complex: " + entity.getModelInterface());
      keggEntry = createKEGGEntryForComplex((Complex)entity, keggPW, m, species);
    } else {
      Integer geneID = getEntrezGeneID(entity, getMapFromSet(m.getObjects(RelationshipXref.class)));      
      if (geneID != -1) {
        keggname = mapGeneIDToKEGGID(geneID, species);
      } else {
        Set<String> names = entity.getName();
        StringBuffer name = new StringBuffer(); 
        if (names != null && names.size() > 0) {
          for (String n : names) {
            n = n.trim();
            n = n.replace(" ", "_");
            name.append(n);
            name.append(" ");
          }
        }      
        graphName = name.toString().trim();
        keggname = keggUndefinedName;
      }
      
      Class<?> entityClass = entity.getClass();
      if (Dna.class.isAssignableFrom(entityClass)) {
        log.finer("Dna: " + entity.getModelInterface() + " " + keggname);
        eType = EntryType.gene;
        gType = GeneType.dna;
      } else if (DnaRegion.class.isAssignableFrom(entityClass)) {
        log.finer("DnaRegion: " + entity.getModelInterface() + " " + keggname);
        eType = EntryType.gene;
        gType = GeneType.dna_region;
      } else if (Protein.class.isAssignableFrom(entityClass)) {
        log.finer("Protein: " + entity.getModelInterface() + " " + keggname);
        eType = EntryType.gene;
        gType = GeneType.protein;
      } else if (Rna.class.isAssignableFrom(entityClass)) {
        log.finer("Rna: " + entity.getModelInterface() + " " + keggname);
        eType = EntryType.gene;
        gType = GeneType.rna;
      } else if (RnaRegion.class.isAssignableFrom(entityClass)) {
        log.finer("RnaRegion: " + entity.getModelInterface() + " " + keggname);
        eType = EntryType.gene;
        gType = GeneType.rna_region;
      } else if (SmallMolecule.class.isAssignableFrom(entityClass)) {
        log.finer("SmallMolecule: " + entity.getModelInterface() + " " + keggname);
        eType = EntryType.compound;
        gType = GeneType.unknown;
      } else {     
        log.finer("GenericPhysicalEntity: " + entity.getModelInterface() + " " + keggname);
        eType = EntryType.other;
        gType = GeneType.unknown;
      }     

      Collection<de.zbit.kegg.parser.pathway.Entry> entries = keggPW.getEntriesForName(keggname);
      boolean alreadyAddedToPathway = false;
      if (entries != null && entries.size() > 0) {
        for (de.zbit.kegg.parser.pathway.Entry entry : entries) {
          EntryExtended e = (EntryExtended) entry;
          if (e.getGeneType() != null && e.getGeneType().equals(gType)
              && (e.getType() != null && e.getType().equals(eType))) {
            // TODO: compounds test?
            if (keggname.equals(keggUndefinedName)) {
              if (e.isSetGraphics()
                  && e.getGraphics().getName().equals(graphName)) {
                keggEntry = e;
                alreadyAddedToPathway = true;
                break;
              }
            } else {
              keggEntry = e;
              alreadyAddedToPathway = true;
              break;
            }
          }
        }
      }

      if (!alreadyAddedToPathway) {
        keggEntry = new EntryExtended(keggPW, getKeggEntryID(), keggname, eType, gType);
        if (graphName != null) {
          keggEntry.addGraphics(new Graphics(graphName));
        }
      }

      log.config(" 2 " + entity.getRDFId() + "-" + entity.getModelInterface()
          + " : " + keggEntry.getId() + "-" + keggEntry.getName());
    }

    log.finer("Entity: " + entity.getRDFId() + " - keggEntryID: "
        + keggEntry.getId() + ".");
    return keggEntry;
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
   * creates an {@link EntryExtended} but does NOT add it to the {@link de.zbit.kegg.parser.pathway.Pathway}
   * 
   * @param complex
   * @param keggPW
   * @param mapper
   * @param m
   * @param species
   * @return
   */
  private de.zbit.kegg.parser.pathway.ext.EntryExtended createKEGGEntryForComplex(Complex complex,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    
    de.zbit.kegg.parser.pathway.ext.EntryExtended keggEntry = null, keggEntry2 = null;    
    Set<PhysicalEntity> components = complex.getComponent();
    
    boolean alreadyAddedToPathway = false;
    Collection<de.zbit.kegg.parser.pathway.Entry> entries = keggPW.getEntriesForName(keggUndefinedName);
    if (entries !=null && entries.size()>0) {      
      for (de.zbit.kegg.parser.pathway.Entry entry : entries) {
        EntryExtended entryExtended = (EntryExtended)entry;
        
        if (entryExtended.getType()!=null && entryExtended.getType().equals(EntryType.group)){
          List<Integer> eComponents = entryExtended.getComponents();
          if(eComponents.size()==components.size()){
            boolean allFound = true;
            for (Integer id : eComponents) {
              if(keggPW.getEntryForId(id.intValue())==null){
                allFound = false;
                break;
              }
            }
            
            if (allFound) {
              alreadyAddedToPathway = true;
              keggEntry = entryExtended;
              break;
            }
          }
        }
      }      
    }
    
    if(!alreadyAddedToPathway){
      keggEntry = new EntryExtended(keggPW, getKeggEntryID(), keggUndefinedName, EntryType.group);
      for (PhysicalEntity physicalEntity : components) {         
        keggEntry2 = createKEGGEntryForPhysicalEntity(physicalEntity, keggPW, m, species);          
        keggEntry.addComponent(keggEntry2.getId());        
      }
    }
    
    log.config(" 2 " + complex.getRDFId() + "-" + complex.getModelInterface() + " : " + keggEntry.getId() + "-" + keggEntry.getName());
    return keggEntry;
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
  private EntryExtended createKEGGEntryForGene(Gene entity, de.zbit.kegg.parser.pathway.Pathway keggPW, 
      Model m, Species species) {
    de.zbit.kegg.parser.pathway.ext.EntryExtended keggEntry = null; 
    
    String keggname = null;
    
    Integer geneID = getEntrezGeneID(entity, getMapFromSet(m.getObjects(RelationshipXref.class)));
    if(geneID != -1){  
      keggname = mapGeneIDToKEGGID(geneID, species);
    } else {
      keggname = keggUndefinedName;
    }
          
    Collection<de.zbit.kegg.parser.pathway.Entry> entries = keggPW.getEntriesForName(keggname);
    boolean alreadyAddedToPathway = false;
    if (entries !=null && entries.size()>0) {
      for (de.zbit.kegg.parser.pathway.Entry entry : entries) {
        EntryExtended e = (EntryExtended)entry;
        if(e.getType()!=null && e.getType().equals(EntryType.gene)){
          keggEntry = e;
          alreadyAddedToPathway = true;          
          break;
        } 
      }      
    }
    
    if(!alreadyAddedToPathway){
      keggEntry = new EntryExtended(keggPW, getKeggEntryID(), keggname, EntryType.gene);
    }
    
    log.config(" 2 " + entity.getRDFId() + "-" + entity.getModelInterface() + " : " + keggEntry.getId() + "-" + keggEntry.getName());
    return keggEntry;
  }
  

  /**
   * 
   * Creates an {@link EntryExtended} and adds it to the {@link BioCartaTools#bc2KeggEntry} map
   * and to the {@link de.zbit.kegg.parser.pathway.Pathway}
   * 
   * @param entity
   * @param keggPW
   * @param mapper
   * @param m
   */
  private de.zbit.kegg.parser.pathway.ext.EntryExtended createKEGGEntryForPathway(Pathway entity,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    
    String keggname = keggUndefinedName;
    String graphName = getName(entity);
    int number  = getKeggPathwayNumber(entity.getRDFId());
   
    EntryExtended keggEntry = null;
    
    Collection<de.zbit.kegg.parser.pathway.Entry> entries = keggPW.getEntriesForName(keggname);
    boolean alreadyAddedToPathway = false;
    if (entries !=null && entries.size()>0) {
      for (de.zbit.kegg.parser.pathway.Entry entry : entries) {
        EntryExtended e = (EntryExtended)entry;
        if(e.getType()!=null && e.getType().equals(EntryType.map)){
          if(e.isSetGraphics() && e.getGraphics().getName().equals(graphName)){
            keggEntry = e;
            alreadyAddedToPathway = true;          
            break;
          } else if (!e.isSetGraphics()) {
            keggEntry = e;
            alreadyAddedToPathway = true;          
            break;
          }            
        }
      }      
    }
    
    if(!alreadyAddedToPathway){
      keggEntry = new EntryExtended(keggPW, getKeggEntryID(), keggname, EntryType.map);
      if(graphName!=null)
        keggEntry.addGraphics(new Graphics(graphName));
    }
    log.config(" 2 " + entity.getRDFId() + "-" + entity.getModelInterface() + " : " + keggEntry.getId() + "-" + keggEntry.getName());
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
    log.config(" 3 " + entity.getModelInterface() + "-" + entity.getRDFId());

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
        EntryExtended keggEntry1 = createKEGGEntryForGene((Gene)participants.get(i), keggPW, m, species);
        EntryExtended keggEntry2 = createKEGGEntryForGene((Gene)participants.get(j), keggPW, m, species);
        
        new Relation(keggPW, keggEntry1.getId(), keggEntry2.getId(), RelationType.other);
        log.config(" 3 " + entity.getRDFId() + "-" + entity.getModelInterface() + " : " 
          + keggEntry1.getId() + "-" + keggEntry1.getName()+ "-" + keggEntry2.getId() + "-" + keggEntry2.getName()
          );
      }      
    }      
    
    // interactionScore (0 or more) - up to now ignored
    
    // phenotype (1) - up to now ignored
  }
  
  /**
   * determine depending on the MI identifier the KEGG {@link RelationType}
   * The default is {@link RelationType#other}
   * @param interactionType
   * @return {@link RelationType}
   */
  private RelationType determineRelationType(Set<InteractionVocabulary> interactionType) {
    if(interactionType.size()>0){
      for (InteractionVocabulary iv : interactionType) {
        //TODO: mapper
//        distinguish between RelationType.maplink, RelationType.PCrel, RelationType.GErel, RelationType.PPrel,
//        RelationType.ECrel
      }
    }
    return RelationType.other;
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
       EntryExtended keggEntry1 = createKEGGEntryForGene((Gene)participants.get(i), keggPW, m, species);
       EntryExtended keggEntry2 = createKEGGEntryForGene((Gene)participants.get(j), keggPW, m, species);
       
       new Relation(keggPW, keggEntry1.getId(), keggEntry2.getId(), RelationType.PPrel);
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
      rel = new Relation(keggPW, keggEntry.getId(), keggEntry.getId(), RelationType.GErel);      
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
          keggEntry1 = createKEGGEntryForPathway((Pathway)controller, keggPW, m, species);
          relType = RelationType.maplink;
        } else {
          log.severe("Controller: " + controller.getModelInterface() + "-This should not happen!");
          System.exit(1);
        }
        if(keggEntry1 != null && controlleds.size()>0) {
          for (org.biopax.paxtools.model.level3.Process process : controlleds) {
            if (Conversion.class.isAssignableFrom(process.getClass())){
              Conversion con = (Conversion) process;
              if(BiochemicalReaction.class.isAssignableFrom(con.getClass())){
                Reaction r = createKEGGReactionForBiochemicalReaction(((BiochemicalReaction)con).getLeft(),
                    ((BiochemicalReaction)con).getRight(), 
                    Utils.getListOfCollection(((BiochemicalReaction)con).getParticipantStoichiometry()),keggPW, m, species);
                keggEntry1.appendReaction(r.getName()); 
              } else if(Transport.class.isAssignableFrom(con.getClass())){
                Reaction r = createKEGGReactionForBiochemicalReaction(((Transport)con).getLeft(),
                    ((Transport)con).getRight(), 
                    Utils.getListOfCollection(((Transport)con).getParticipantStoichiometry()),keggPW, m, species);
                keggEntry1.appendReaction(r.getName()); 
              }
            } else if (Pathway.class.isAssignableFrom(process.getClass())){
              keggEntry1 = createKEGGEntryForPathway((Pathway)process, keggPW, m, species);
            } else {
              log.severe("Process: " + process.getModelInterface() + "-This should not happen!");
              System.exit(1);
            }
          }
          
          if (keggEntry2 != null) {
            
            new Relation(keggPW, keggEntry1.getId(), keggEntry2.getId(), relType);
          }
//          log.config(" 3 " + control.getRDFId() + "-" + control.getModelInterface() + " : " 
//              + keggEntry1.getId() + "-" + keggEntry1.getName()+ keggEntry2.getId() + "-" + keggEntry2.getName());
        } else {
          break;
        }
      }
    }
    
    return keggEntry1;
  }
  
  private EntryExtended createKEGGReactionForCatalysis(Catalysis catalysis,
      de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    EntryExtended keggEntry1 = null, keggEntry2 = null;
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
          keggEntry1 = createKEGGEntryForPathway((Pathway)controller, keggPW, m, species);
          relType = RelationType.maplink;
        } else {
          log.severe("Controller: " + controller.getModelInterface() + "-This should not happen!");
          System.exit(1);
        }
        if(keggEntry1 != null && controlleds.size()>0) {
          for (org.biopax.paxtools.model.level3.Process process : controlleds) {
            if (Conversion.class.isAssignableFrom(process.getClass())){
              Conversion con = (Conversion) process;
              if(BiochemicalReaction.class.isAssignableFrom(con.getClass())){
                Reaction r = createKEGGReactionForBiochemicalReaction(((BiochemicalReaction)con).getLeft(),
                    ((BiochemicalReaction)con).getRight(), 
                    Utils.getListOfCollection(((BiochemicalReaction)con).getParticipantStoichiometry()),keggPW, m, species);
                keggEntry1.appendReaction(r.getName()); 
              }
            } else {
              log.severe("Process: " + process.getModelInterface() + "-This should not happen!");
              System.exit(1);
            }
          }
          
          if (keggEntry2 != null) {
            // ControlType (0 or 1)

            // Cofactor = PhysicalEntity (0 or more)

            // CatalysisDirection

            if (catalysis.getCatalysisDirection() != null) {
              CatalysisDirectionType catType = catalysis
                  .getCatalysisDirection();
              if (catType == CatalysisDirectionType.LEFT_TO_RIGHT) {
                // TODO:
                log.fine("Catalysis: left-to-right");
              } else if (catType == CatalysisDirectionType.RIGHT_TO_LEFT) {
                // TODO:
                log.fine("Catalysis: right-to-left");
              } else {
                // TODO: unknown catType???
                log.fine("Catalysis: unknown");
              }
            }
            if (catalysis.getCofactor() != null) {
              // TODO: up to now there was no cofactor in catalysis, but what to
              // do with this thing?
              for (PhysicalEntity physEnt : catalysis.getCofactor()) {
                log.fine("Catalysis: cofactor=" + physEnt.getModelInterface());
              }
            }
            new Relation(keggPW, keggEntry1.getId(), keggEntry2.getId(), relType);
          }
//          log.config(" 3 " + catalysis.getRDFId() + "-" + catalysis.getModelInterface() + " : " 
//              + keggEntry1.getId() + "-" + keggEntry1.getName()+ keggEntry2.getId() + "-" + keggEntry2.getName());
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
    
    if (controllers.size()>0){
      for (Controller entity : controllers) {
        if (PhysicalEntity.class.isAssignableFrom(entity.getClass())){
          if (controlleds.size()>0) {
            for (org.biopax.paxtools.model.level3.Process process : controlleds) {
              if (Catalysis.class.isAssignableFrom(process.getClass())){
                //TODO: 
//                does catalaysis exist?
//                if (no) {
//                  create catalysis
//                } 
                EntryExtended keggEntry1 = createKEGGEntryForPhysicalEntity((PhysicalEntity) entity, keggPW, m, species);
                EntryExtended keggEntry2 = createKEGGReactionForCatalysis((Catalysis) process, keggPW, m, species);
                new Relation(keggPW, keggEntry1.getId(), keggEntry2.getId(), RelationType.other);
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
//        + keggEntry.getId() + "-" + keggEntry1.getName()+ keggEntry2.getId() + "-" + keggEntry2.getName()
        );
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
                //TODO: 
//                does TemplateReaction exist?
//                if (no) {
//                  create TemplateReaction
//                } 
                EntryExtended keggEntry1 = createKEGGEntryForPhysicalEntity((PhysicalEntity)entity, keggPW, m, species);
                Relation r = createKEGGRelationForTemplateReaction((TemplateReaction)tempReaction, keggPW, m, species);
                new Relation(keggPW, keggEntry1.getId(), r.getEntry1(), RelationType.other); //TODO: to simply fetch the entry1 from
                                                                                             // the relation is wrong!!!
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
    log.config(" 3 " + tmpReaction.getRDFId() + "-" + tmpReaction.getModelInterface() + " : " 
//      + keggEntry.getId() + "-" + keggEntry1.getName()+ keggEntry2.getId() + "-" + keggEntry2.getName()
      );
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
    log.config(" 3 " + entity.getModelInterface() + "-" + entity.getRDFId());        
    if (ComplexAssembly.class.isAssignableFrom(entity.getClass())){
      createKEGGReactionForBiochemicalReaction(((ComplexAssembly)entity).getLeft(),
          ((ComplexAssembly)entity).getRight(), 
          Utils.getListOfCollection(((ComplexAssembly)entity).getParticipantStoichiometry()), keggPW, m, species);
    } else if (BiochemicalReaction.class.isAssignableFrom(entity.getClass())){
      createKEGGReactionForBiochemicalReaction(((BiochemicalReaction)entity).getLeft(),
          ((BiochemicalReaction)entity).getRight(), 
          Utils.getListOfCollection(((BiochemicalReaction)entity).getParticipantStoichiometry()), keggPW, m, species);
    } else if (Degradation.class.isAssignableFrom(entity.getClass())){
      // TODO: implement      
    } else if (Transport.class.isAssignableFrom(entity.getClass())){
      createKEGGReactionForBiochemicalReaction(((Transport)entity).getLeft(),
          ((Transport)entity).getRight(), 
          Utils.getListOfCollection(((Transport)entity).getParticipantStoichiometry()), keggPW, m, species);
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

 

  /**
   * 
   * @param entity
   * @param keggPW
   * @param m
   * @param species
   */
  private Reaction createKEGGReactionForBiochemicalReaction(Set<PhysicalEntity> lefts, Set<PhysicalEntity> rights, 
      List<Stoichiometry> stoichiometry, de.zbit.kegg.parser.pathway.Pathway keggPW, Model m, Species species) {
    String rName = getReactionName();
    
    List<ReactionComponent> products = new ArrayList<ReactionComponent>();
    List<ReactionComponent> substrates = new ArrayList<ReactionComponent>();
    ReactionType rType = ReactionType.other;
    // TODO: Reaction type must be set!!! How to distinguish between irreversible or reversible

//    sb.append(name);
    for (PhysicalEntity left : lefts) {
      EntryExtended keggEntry = createKEGGEntryForPhysicalEntity(left, keggPW, m, species);   
      ReactionComponent rc = new ReactionComponent(keggEntry.getId(), keggEntry.getName());
      
      for (Stoichiometry stoichi : stoichiometry) {
        if (stoichi.getPhysicalEntity().equals(left)) {
          rc.setStoichiometry((int) stoichi.getStoichiometricCoefficient());
          break;
        }
      }
//      sb.append("product: " + left.getRDFId() + left.getDisplayName() + "-" + keggEntry.getName() + ",");
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
//      sb.append("substrate: " + right.getRDFId() + right.getDisplayName() + "-" + keggEntry.getName() + ",");
      products.add(rc);
    }
     
    Reaction r = null;            
    boolean reactionExists = false;
    for (Reaction react : keggPW.getReactions()) {
      List<ReactionComponent> reactProds = react.getProducts();
      List<ReactionComponent> reactSubs = react.getSubstrates();
      
      if (react.getType().equals(rType) &&
          reactProds.size() == products.size()
          && reactSubs.size() == substrates.size()) {
        boolean rcMissing = false;
        for (ReactionComponent sub : substrates) {
          if (!reactSubs.contains(sub)) { 
            rcMissing = true;
            break;
          }
        }
        if (!rcMissing) { 
          for (ReactionComponent prod : products) {
            if (!reactProds.contains(prod)) {
              rcMissing = true;
              break;
            }
          }
          if (!rcMissing) {
            reactionExists = true;
            break;
          }
        }
        
      }
    }
    
    if(!reactionExists){
      r = new Reaction(keggPW, rName, ReactionType.other);
      r.addProducts(products);
      r.addSubstrates(substrates);
    } else{
      keggReactionID--;
    }
//    log.info(sb.toString());
    
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


  private void createRelation(Control ent,
      de.zbit.kegg.parser.pathway.Pathway keggPW, GeneID2KeggIDMapper mapper,
      Model m, Species species) {
    Set<Controller> controllers = ent.getController();
    Set<org.biopax.paxtools.model.level3.Process> controlleds = ent.getControlled();
    
    // create a relation
    for (Controller controller : controllers) {
      
      for (org.biopax.paxtools.model.level3.Process process : controlleds) {
        if((PhysicalEntity.class.isAssignableFrom(controller.getClass())) && 
        Pathway.class.isAssignableFrom(process.getClass())){  
          //case 1 controller = physical entity, controlled = pathway
          //TODO: easy
          log.fine("Create Relation physEnt-pathway.");
          log.fine("  ---  " + controller.getRDFId() +"|" + controller.getModelInterface() + "-" 
                             + process.getRDFId()    +"|" + process.getModelInterface());
        } else if(Pathway.class.isAssignableFrom(controller.getClass()) &&  
            Pathway.class.isAssignableFrom(process.getClass())){           
          //case 2 controller = pathway, controlled = Pathway
          //TODO: easy
          log.fine("Create Relation pathway-pathway.");
          log.fine("  ---  " + controller.getRDFId() +"|" + controller.getModelInterface() + "-" 
              + process.getRDFId()    +"|" + process.getModelInterface());
        } else if (Pathway.class.isAssignableFrom(controller.getClass()) &&  
            Interaction.class.isAssignableFrom(process.getClass())){  
          //case 3 controller = pathway, controlled = interaction
          //TODO: if it is a template regulation-> cool we just can use the product
          //      what to do with another Control object???
          log.fine("Create Relation pathway-interaction.");
          log.fine("  ---  " + controller.getRDFId() +"|" + controller.getModelInterface() + "-" 
              + process.getRDFId()    +"|" + process.getModelInterface());
        } else {
          log.log(Level.WARNING, "Something went wrong: " + controller.getRDFId() +"|" + controller.getModelInterface() + "-" 
              + process.getRDFId()    +"|" + process.getModelInterface());
        }
      }

      if (Catalysis.class.isAssignableFrom(ent.getClass())){
        Catalysis catalysis =  (Catalysis) ent;
        if (catalysis.getCatalysisDirection()!= null) {
          CatalysisDirectionType catType = catalysis.getCatalysisDirection();
          if (catType == CatalysisDirectionType.LEFT_TO_RIGHT){
            //TODO:
            log.fine("Catalysis: left-to-right");
          } else if (catType == CatalysisDirectionType.RIGHT_TO_LEFT) {
            //TODO:
            log.fine("Catalysis: right-to-left");
          } else {
            //TODO: unknown catType???
            log.fine("Catalysis: unknown");
          }
        }
        if (catalysis.getCofactor() != null) {
          //TODO: up to now there was no cofactor in catalysis, but what to do with this thing?
          for (PhysicalEntity physEnt : catalysis.getCofactor()) {
            log.fine("Catalysis: cofactor=" + physEnt.getModelInterface());
          }          
        }
      } 
    }    
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
    //TODO: Check if this is valid for all BioPax formats!!!
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
      GeneSymbol2GeneIDMapper mapper = null;
      try {
        mapper = new GeneSymbol2GeneIDMapper(species);
      } catch (IOException e) {
        log.log(Level.SEVERE, "Could not initalize mapper!", e);
        e.printStackTrace();
      } 
      
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
          BioCartaPathwayHolder helper = new BioCartaPathwayHolder(pw.getRDFId(), getName(pw));
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
