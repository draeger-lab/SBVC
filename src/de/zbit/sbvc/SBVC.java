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
 * BioPax, SBGN, etc. Please visit the project homepage at
 * <http://www.cogsys.cs.uni-tuebingen.de/software/SBVC> to obtain the
 * latest version of SBVC.
 *
 * Copyright (C) 2012-2012 by the University of Tuebingen, Germany.
 *
 * SBVC is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation. A copy of the license
 * agreement is provided in the file named "LICENSE.txt" included with
 * this software distribution and also available online as
 * <http://www.gnu.org/licenses/lgpl-3.0-standalone.html>.
 * ---------------------------------------------------------------------
 */
package de.zbit.sbvc;

import java.awt.Window;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import de.zbit.AppConf;
import de.zbit.Launcher;
import de.zbit.biopax.BioPax2KGML;
import de.zbit.kegg.KEGGtranslatorOptions.NODE_NAMING;
import de.zbit.kegg.KGMLWriter;
import de.zbit.kegg.Translator;
import de.zbit.kegg.api.cache.KeggInfoManagement;
import de.zbit.kegg.io.AbstractKEGGtranslator;
import de.zbit.kegg.io.KEGG2SBMLqual;
import de.zbit.sbvc.io.SBVCIOOptions;
import de.zbit.util.Utils;
import de.zbit.util.prefs.KeyProvider;
import de.zbit.util.prefs.SBProperties;

/**
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class SBVC extends Launcher{
  
  /**
   * 
   */
  private static final long serialVersionUID = 5858547705130260355L;
  
  /**
   * The {@link Logger} for this class.
   */
  private static final transient Logger log = Logger.getLogger(Translator.class.getName());

  public SBVC(String[] args) {
    super(args);
  }

  public SBVC() {
    // TODO Auto-generated constructor stub
  }

  @Override
  public void commandLineMode(AppConf appConf) {
    SBProperties props = appConf.getCmdArgs();
    String folderName = SBVCIOOptions.OUTPUT.getValue(props).getPath();
    String input = SBVCIOOptions.INPUT.getValue(props).getPath();
    boolean splitMode = (boolean)(SBVCIOOptions.SPLIT_MODE.getValue(props)); 
    
    parseInputFile(input, folderName, splitMode);
  }

  public void parseInputFile(String input, String folderName, boolean splitMode) {
    // getting the KEGG Pathways of the model
    Collection<de.zbit.kegg.parser.pathway.Pathway> keggPWs = 
      BioPax2KGML.createPathwaysFromModel(input, folderName, splitMode, false);
    
    // translation to sbml
    KEGG2SBMLqual k2s = null;
    if (new File(Translator.cacheFileName).exists()
        && new File(Translator.cacheFileName).length() > 1) {
      KeggInfoManagement manager = null;
      try {
        manager = (KeggInfoManagement) KeggInfoManagement.loadFromFilesystem(Translator.cacheFileName);
        k2s = new KEGG2SBMLqual(manager);
      } catch (IOException e) {
        k2s = new KEGG2SBMLqual();
      }      
    } else {
      k2s = new KEGG2SBMLqual();
    }
    
    // necessary that both reactions and relations are written to the file
    k2s.setConsiderReactions(true);
    k2s.setAddCellDesignerAnnots(false);
    k2s.setNameToAssign(NODE_NAMING.INTELLIGENT);
    k2s.setRemoveOrphans(false);
    k2s.setAutocompleteReactions(false);
    k2s.setRemoveWhiteNodes(false);
    k2s.setShowFormulaForCompounds(false);
    k2s.setRemovePathwayReferences(false);
    k2s.setAddLayoutExtension(false);
    k2s.setCheckAtomBalance(false);
    
    for (de.zbit.kegg.parser.pathway.Pathway p : keggPWs) {
      k2s.translate(p, Utils.ensureSlash(folderName) + KGMLWriter.createFileName(p));
    }
    
    // Remember already queried objects (save cache)
    if (AbstractKEGGtranslator.getKeggInfoManager().hasChanged()) {
      KeggInfoManagement.saveToFilesystem(Translator.cacheFileName, AbstractKEGGtranslator.getKeggInfoManager());
    }
  }

  @Override
  public List<Class<? extends KeyProvider>> getCmdLineOptions() {
    //XXX: we can add further IOOptions later, for example those of KEGG translator
    List<Class<? extends KeyProvider>> configList = new ArrayList<Class<? extends KeyProvider>>(3);
    configList.add(SBVCIOOptions.class);
    return configList;
  }

  @Override
  public List<Class<? extends KeyProvider>> getInteractiveOptions() {
    // TODO Auto-generated method stub
    return null;
  }

  @Override
  public URL getURLlicenseFile() {
    URL url = null;
    try {
      url = new URL("http://www.gnu.org/licenses/lgpl-3.0-standalone.html");
    } catch (MalformedURLException exc) {
      log.log(Level.FINE, exc.getLocalizedMessage(), exc);
    }
    return url;
  }

  @Override
  public URL getURLOnlineUpdate() {
    // TODO Auto-generated method stub
    return null;
  }

  @Override
  public String getVersionNumber() {
    return "0.1";
  }

  @Override
  public short getYearOfProgramRelease() {    
    return (short) 2012;
  }

  @Override
  public short getYearWhenProjectWasStarted() {
    return (short) 2012;
  }

  @Override
  public Window initGUI(AppConf appConf) {
    // TODO Auto-generated method stub
    return null;
  }
  
  /**
   * 
   * @param args
   */
  public static void main(String args[]) {
    // --input files/KGMLsamplefiles/hsa00010.xml --format GraphML --output test.txt
    new SBVC(args);
  }

}
