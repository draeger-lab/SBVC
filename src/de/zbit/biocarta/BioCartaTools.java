package de.zbit.biocarta;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
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
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level2.ControlType;
import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.model.level3.ComplexAssembly;
import org.biopax.paxtools.model.level3.Control;
import org.biopax.paxtools.model.level3.Controller;
import org.biopax.paxtools.model.level3.Entity;
import org.biopax.paxtools.model.level3.Interaction;
import org.biopax.paxtools.model.level3.Pathway;
import org.biopax.paxtools.model.level3.PhysicalEntity;
import org.biopax.paxtools.model.level3.Process;
import org.biopax.paxtools.model.level3.RelationshipXref;
import org.biopax.paxtools.model.level3.Stoichiometry;
import org.biopax.paxtools.model.level3.Xref;

import de.zbit.mapper.GeneSymbol2GeneIDMapper;
import de.zbit.util.logging.LogUtil;

public class BioCartaTools {

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
    
    return getModel(io);
  }
  
  public Model getModel(InputStream io) {
    BioPAXIOHandler handler = new SimpleIOHandler();//TODO: check level
    return handler.convertFromOWL(io);
  }

  public Set<Pathway> getPathways(Model m) {
    return m.getObjects(Pathway.class);

  }

  public void getEntryGeneIDsOfPathwayElements() {

  }

  public static void getControllEntities(Control control, ControlType ctype) {
    log.info(control.getName().toString());

    if (ctype == null) {
      log.info("ControlType: " + control.getControlType());
      getControlerName(control);
      getControlledNames(control);
    } else if (control.getControlType().equals(ctype)) {
      log.info("ControlType: " + control.getControlType());
      getControlerName(control);
      getControlledNames(control);
    }
  }

  private static void getControlledNames(Control control) {
    log.info("Controlled: ");
    for (Process p : control.getControlled()) {
      if (p instanceof Interaction) {
        Interaction i = (Interaction) p;
        for (Entity ent : i.getParticipant()) {
          log.info("\t" + ent.getName());
          for (Xref xref : ent.getXref()) {
            log.info("Ref: " + xref.getDb() + ":" + xref.getId());
          }
        }
      }
    }
  }

  private static void getControlerName(Control control) {
    for (Controller ctrl : control.getController()) {
      log.info("\t" + ctrl.getName());
    }
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
 * call therefore the method {@link BioCartaTools#getEntitiesWithName(Entity)}  
 *    
 * @param entities 
 * @param species
 * @param xrefs
 * @return
 */
  private static Map<Integer, Entity> getEntrezGeneIDs(Set<Entity> entities, String species, Map<String, RelationshipXref> xrefs) {
    Map<Integer, Entity> geneIDs = new HashMap<Integer, Entity>();    
    
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
        geneIDs.put(geneID, entity);
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
      
      geneIDs.putAll(getGeneIDOverGeneSymbol(mapper, geneSymbols));
    } 
    
    return geneIDs;
  }

  /**
   * Maps the entered gene symbol names to a geneID
   * 
   * @param mapper {@link GeneSymbol2GeneIDMapper} 
   * @param nameSet with names of one {@link BioPAXElement} element 
   * @return the gene id o
   */
  private static Map<Integer, Entity> getGeneIDOverGeneSymbol(GeneSymbol2GeneIDMapper mapper, Map<String, Entity> nameSet) {
    log.finest("getGeneIDOverGeneSymbol");
    Map<Integer, Entity> geneIDs = new HashMap<Integer, Entity>();
    Integer geneID = null;
    
      
    for (Entry<String, Entity> symbol : nameSet.entrySet()) {
  
      try {
        geneID = mapper.map(symbol.getKey());        
      } catch (Exception e) {
        log.log(Level.WARNING, "Error while mapping name: " + symbol.getKey() + ".", e);
      }
      
      if(geneID!=null){
        log.log(Level.FINER, "----- found! Geneid: " + geneID + " "  + symbol.getValue().getRDFId() + " "+ symbol.getKey());
        geneIDs.put(geneID, symbol.getValue());
      } else if (symbol.getKey().contains("-")) {
        log.log(Level.FINER, "recall for symbol: " + symbol.getValue().getRDFId() + " " + symbol.getKey());
        Map<String, Entity> set = new HashMap<String, Entity>();
        set.put(symbol.getKey().replace("-", ""), symbol.getValue());        
        geneIDs.putAll(getGeneIDOverGeneSymbol(mapper, set));
      } else if (symbol.getKey().contains(" ")) {
        log.log(Level.FINER, "recall for symbol: " + symbol.getValue().getRDFId() + " " + symbol.getKey());
        Map<String, Entity> set = new HashMap<String, Entity>();
        set.put(symbol.getKey().replace(" ", "_"), symbol.getValue());      
        geneIDs.putAll(getGeneIDOverGeneSymbol(mapper, set));
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
  private static <T extends BioPAXElement> Map<String, T> getMapFromSet(Set<T> set) {
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

    if(name.size()>0 && !(entity instanceof Pathway)){
      if(entity instanceof Complex){
        Complex c = (Complex)entity;
        for (PhysicalEntity pe : c.getComponent()) {
          resEntities.addAll(getEntitiesWithName(pe));
        }
      } else if(entity instanceof ComplexAssembly){
        ComplexAssembly c = (ComplexAssembly)entity;
        for (Stoichiometry pe : c.getParticipantStoichiometry()) {
          resEntities.addAll(getEntitiesWithName((Entity) pe));
        }
        
      } else {
        resEntities.add(entity);
      }
    } else if (entity.getParticipantOf().size()>0 && !(entity instanceof Pathway)){      
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
          if(!(entity instanceof Pathway))
            log.log(Level.FINER, "--Input: " + entity.getRDFId() + "\t" + entity.getModelInterface());
        }

        Map<Integer, Entity> geneIDs = getEntrezGeneIDs(pwEntities, species, xrefs);
        for (Entry<Integer, Entity> entity : geneIDs.entrySet()) {
          log.log(Level.FINER, "----res: " + entity.getKey() + " " + entity.getValue());
          pw.addGeneID(entity.getKey());
        }        
      }
    }
    
    return pathways;
  }

  /**
   * @param args
   * @throws FileNotFoundException
   */
  public static void main(String[] args) throws FileNotFoundException {
    String species = "human";    
    LogUtil.initializeLogging(Level.FINER);
    
    BioCartaTools bf = new BioCartaTools();    
    
    
    // test for pathway gene ids
    //    bioCartafile Level 3(can be downloaded from http://pid.nci.nih.gov/download.shtml)
    Model m = bf.getModel("C:/Users/buechel/Downloads/BioCarta.bp3.owl");
    bf.getPathwaysWithGeneID(species, m);
   
  }


}
