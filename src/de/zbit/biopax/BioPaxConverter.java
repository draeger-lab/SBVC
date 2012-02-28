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

import org.biopax.paxtools.io.BioPAXIOHandler;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.model.level3.ComplexAssembly;
import org.biopax.paxtools.model.level3.Entity;
import org.biopax.paxtools.model.level3.Interaction;
import org.biopax.paxtools.model.level3.Pathway;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.model.level3.RelationshipXref;
import org.biopax.paxtools.model.level3.Stoichiometry;

import de.zbit.mapper.GeneID2KeggIDMapper;
import de.zbit.mapper.GeneSymbol2GeneIDMapper;
import de.zbit.parser.Species;
import de.zbit.util.ProgressBar;
import de.zbit.util.Utils;

/**
 * This class is a base class to convert BioPax elements
 * 
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class BioPaxConverter {

  public static final Logger log = Logger.getLogger(BioPaxConverter.class.getName());

  /**
   * mapper to map gene symbols to gene ids
   */
  protected static GeneSymbol2GeneIDMapper geneSymbolMapper = null;

  /**
   * mapper to map gene ids to KEGG ids
   */
  protected static GeneID2KeggIDMapper geneIDKEGGmapper = null;

  /**
   * This method maps the gene id to the gene symbol of the entered entity. This
   * method is not able to map complexes call for that the method
   * {@link BioPax2KGML#getEntrezGeneIDsOfComplex(Complex, Map)} It's an
   * advantage to preprocess the entities to exclusively having entities with a
   * name call therefore the method
   * {@link BioPax2KGML#getEntitiesWithName(Entity)}
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
   * This method maps to all gene symbols of the entered entities the
   * corresponding gene id It's an advantage to preprocess the entities to
   * exclusively having entities with a name call therefore the method
   * {@link BioPax2KGML#getEntitiesWithName(Entity)} This method also considers
   * 
   * @link {@link Complex} classes. NOTE: the method is not so clean and should
   *       be rewritten!!
   * 
   * @param entities
   * @param species
   * @param xrefs
   * @return
   */
  private static Map<Entity, Integer> getEntrezGeneIDs(Set<Entity> entities, String species,
      Map<String, RelationshipXref> xrefs) {
    Map<Entity, Integer> geneIDs = new HashMap<Entity, Integer>();

    // searching for gene ids
    Map<String, Entity> geneSymbols = new HashMap<String, Entity>();
    for (Entity entity : entities) {
      Integer geneID = null;
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

    for (Entry<String, Entity> symbol : nameSet.entrySet()) {
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
   * deteremines the gene ids of the elements in a pathway
   * 
   * This method is not so clean should be rewritten, becuase in the method
   * {@link BioPax2KGML#getEntrezGeneIDsForPathways(List, String, Model)}
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
      Set<Entity> pwEntities = new HashSet<Entity>();
      for (Entity entity : pw.entities) {
        pwEntities.addAll(getEntitiesWithName(entity));
        if (!(Pathway.class.isAssignableFrom(entity.getClass())))// entity
                                                                 // instanceof
                                                                 // Pathway ||
                                                                 // entity
                                                                 // instanceof
                                                                 // PathwayImpl))
          log.log(Level.FINER, "--Input: " + entity.getRDFId() + "\t" + entity.getModelInterface());
      }

      Map<Entity, Integer> geneIDs = getEntrezGeneIDs(pwEntities, species, xrefs);
      for (Entry<Entity, Integer> entity : geneIDs.entrySet()) {
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
  private static Collection<? extends Entity> getEntitiesWithName(Entity entity) {
    Set<Entity> resEntities = new HashSet<Entity>();
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
        Set<Entity> pwEntities = new HashSet<Entity>();
        for (Entity entity : pw.entities) {
          pwEntities.addAll(getEntitiesWithName(entity));
          if (!(Pathway.class.isAssignableFrom(entity.getClass())))
            log.log(Level.FINER,
                "--Input: " + entity.getRDFId() + "\t" + entity.getModelInterface());
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

  /**
   * Converts the inputStream of an owl file containing BioPAX entries
   * 
   * @param io
   * @return
   */
  public static Model getModel(InputStream io) {
    BioPAXIOHandler handler = new SimpleIOHandler();
    // TODO: check if level is level 3 or higher, because this class is
    // constructed to work for relations
    return handler.convertFromOWL(io);
  }

  /**
   * Calls the method {@link BioPax2KGML#getModel(InputStream) for an entered
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
  
}
