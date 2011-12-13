package de.zbit.biocarta;
import java.util.HashSet;
import java.util.Set;

import org.biopax.paxtools.model.level3.Entity;


public class BioCartaPathwayHolder {

  Set<Entity> entities = null;
  Set<Integer> geneIDs = null;
  String name;
  
  BioCartaPathwayHolder(String name){
    this.name = name;
  }
  
  public boolean addEntity(Entity ent){
    if(entitiesIsSet()){
      if (!entities.contains(ent)) {
        entities.add(ent);
        return true;
      }
    }
    return false;
  }

  private boolean entitiesIsSet() {
   if(entities == null)
    return createEntities();
   else
     return true;
  }

  private boolean createEntities() {
    entities = new HashSet<Entity>();
    return true;
  }
  
  public String getPathwayName(){
    return name;
  }

  public int getNoOfEntities() {
    return entitiesIsSet()? entities.size() : 0;
  }
  
  @Override
  public boolean equals(Object obj) {
    if (obj instanceof BioCartaPathwayHolder){
      if(((BioCartaPathwayHolder) obj).getPathwayName().equals(this.getPathwayName())){
        return true;
      }
      else 
        return false;
    } else 
      return false;
  }
  
  @Override
  public int hashCode() {   
    return getPathwayName().hashCode()*13;
  }

  public boolean addGeneID(Integer key) {
    if(geneIDsIsSet()){
      return geneIDs.add(key);
    }
    return false;
  }

  private boolean geneIDsIsSet() {
    if(geneIDs == null)
      return createGeneIDs();
     else
       return true;
  }

  private boolean createGeneIDs() {
    geneIDs = new HashSet<Integer>();
    return true;
  }

  public Set<Integer> getGeneIDs() {
    return geneIDs;
  }

  public String getName() {    
    return name;
  }
}
