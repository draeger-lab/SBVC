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
 * BioPAX, SBGN, etc. Please visit the project homepage at
 * <http://www.cogsys.cs.uni-tuebingen.de/software/SBVC> to obtain the
 * latest version of SBVC.
 *
 * Copyright (C) 2012-2014 by the University of Tuebingen, Germany.
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

import java.util.logging.Level;
import java.util.logging.Logger;

import de.zbit.util.Species;
import de.zbit.util.logging.LogUtil;

/**
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class SBVCTest {

  public static final Logger log = Logger.getLogger(SBVCTest.class.getName());
  /**
   * @param args
   */
  public static void main(String[] args) {
    LogUtil.initializeLogging(Level.FINE);
    SBVC sbvc = new SBVC();
    
    sbvc.convertBioPAXToSBML("C:/Users/buechel/Downloads/pathway-biopax.owl",
        "C:/Users/buechel/Downloads/res/", 
        null);
    
    if(true)return;
    sbvc.convertBioPAXToSBML("C:/Users/buechel/Downloads/PID_Pathways/BioCarta.bp2.owl",
        "C:/Users/buechel/Downloads/PID_Pathways/BioCarta Level 2/", new Species("Homo sapiens", "_HUMAN", "human", "hsa",9606));
    
    sbvc.convertBioPAXToSBML("C:/Users/buechel/Downloads/PID_Pathways/BioCarta.bp3.owl",
        "C:/Users/buechel/Downloads/PID_Pathways/BioCarta Level 3/", new Species("Homo sapiens", "_HUMAN", "human", "hsa",9606));
    
    
    sbvc.convertBioPAXToSBML("C:/Users/buechel/Downloads/PID_Pathways/NCI-Nature_Curated.bp2.owl",
        "C:/Users/buechel/Downloads/PID_Pathways/NCI-Nature Curated Level 2/", new Species("Homo sapiens", "_HUMAN", "human", "hsa",9606));
    
    sbvc.convertBioPAXToSBML("C:/Users/buechel/Downloads/PID_Pathways/NCI-Nature_Curated.bp3.owl",
        "C:/Users/buechel/Downloads/PID_Pathways/NCI-Nature Curated Level 3/", new Species("Homo sapiens", "_HUMAN", "human", "hsa",9606));
    
    sbvc.convertBioPAXToSBML("C:/Users/buechel/Downloads/PID_Pathways/Reactome.bp2.owl",
        "C:/Users/buechel/Downloads/PID_Pathways/Reactome Level 2/", new Species("Homo sapiens", "_HUMAN", "human", "hsa",9606));
    
    sbvc.convertBioPAXToSBML("C:/Users/buechel/Downloads/PID_Pathways/Reactome.bp3.owl",
        "C:/Users/buechel/Downloads/PID_Pathways/Reactome Level 3/", new Species("Homo sapiens", "_HUMAN", "human", "hsa",9606));
    
    
//    sbvc.convertBioPAXToSBML("C:/Users/buechel/Desktop/workspace/SBVC/doc/Bioinformatics_paper/Comparision/Homarus americanus.owl",
//        "C:/Users/buechel/Desktop/workspace/SBVC/doc/Bioinformatics_paper/Comparision/");
////    
//    sbvc.convertBioPAXToSBML("C:/Users/buechel/Desktop/workspace/SBVC/doc/Bioinformatics_paper/Comparision/Human immunodeficiency virus 1.owl",
//    "C:/Users/buechel/Desktop/workspace/SBVC/doc/Bioinformatics_paper/Comparision/");
//
//    
//    sbvc.convertBioPAXToSBML("C:/Users/buechel/Desktop/workspace/SBVC/doc/Bioinformatics_paper/Comparision/ceramidepathway_level2.owl",
//    "C:/Users/buechel/Desktop/workspace/SBVC/doc/Bioinformatics_paper/Comparision/");
    
  }

}
