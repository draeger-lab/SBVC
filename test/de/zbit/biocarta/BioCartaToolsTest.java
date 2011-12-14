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

import java.io.FileNotFoundException;
import java.util.logging.Level;

import org.biopax.paxtools.model.Model;

import de.zbit.util.logging.LogUtil;

/**
 * Class to test the BioCartaTools
 * 
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class BioCartaToolsTest {
  
  private BioCartaTools bc = null;
  
  public BioCartaToolsTest(){
    bc = new BioCartaTools();
  }
  
  /**
   * method to test the {@link BioCartaTools#getPathwaysWithGeneID(String, Model)}
   * Be carefull this method uses a {@link BioCartaTools#getModel(String)} call where a local BioCarta file
   * of level 3 is needed. The file could be downloaded from http://pid.nci.nih.gov/download.shtml
   */
  private void testGetPathwaysWithGeneID() {
    String species = "human";    
    // test for pathway gene ids
    //    bioCartafile Level 3(can be downloaded from http://pid.nci.nih.gov/download.shtml)
    Model m = bc.getModel("C:/Users/buechel/Downloads/BioCarta.bp3.owl");
    bc.getPathwaysWithGeneID(species, m);
  }
  
  /**
   * @param args
   * @throws FileNotFoundException
   */
  public static void main(String[] args) throws FileNotFoundException {
    LogUtil.initializeLogging(Level.FINER);
    
    BioCartaToolsTest bft = new BioCartaToolsTest();
    
    bft.testGetPathwaysWithGeneID();

  }

  
}
