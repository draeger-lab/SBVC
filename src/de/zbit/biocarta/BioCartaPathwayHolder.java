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
import java.util.HashSet;
import java.util.Set;

import org.biopax.paxtools.model.level3.Entity;

/**
 * This class stores one BioCarta pathway, its BioCarta id, the standard name, 
 * and the contained entities and gene ids in {@link Set} objects 
 * 
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
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
