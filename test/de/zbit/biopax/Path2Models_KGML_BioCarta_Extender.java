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
package de.zbit.biopax;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.logging.FileHandler;
import java.util.logging.Filter;
import java.util.logging.Level;
import java.util.logging.LogRecord;
import java.util.logging.Logger;

import org.biopax.paxtools.model.Model;

import com.ibm.icu.text.SimpleDateFormat;

import de.zbit.io.DirectoryParser;
import de.zbit.kegg.KGMLWriter;
import de.zbit.kegg.Translator;
import de.zbit.kegg.io.BatchKEGGtranslator;
import de.zbit.kegg.io.KEGGtranslator;
import de.zbit.kegg.io.KEGGtranslatorIOOptions.Format;
import de.zbit.kegg.parser.KeggParser;
import de.zbit.kegg.parser.pathway.Pathway;
import de.zbit.util.Utils;
import de.zbit.util.logging.LogUtil;
import de.zbit.util.logging.OneLineFormatter;

/**
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class Path2Models_KGML_BioCarta_Extender {

  public static final Logger log = Logger.getLogger(Path2Models_KGML_BioCarta_Extender.class.getName());
  
  /**
   * for extending existing KGML files with further information for the relations
   * @param fileList (just the filename)
   * @param informationFile containing new relation information, i.e. BioCarta data
   * @param writeEntryExtended 
   */
  private static void createExtendedKGML(List<String> fileList, String informationFile,
      String sourceDir, String destDir, boolean writeEntryExtended, boolean createInstantSBML) {
    log.info("fileList consists of " + fileList.size() + " entries, informationFile= '" + 
        informationFile + "', sourceDir= '" + sourceDir + "', destDir= '" + destDir + "', " +
            "writeExtended= '" + writeEntryExtended + "'.");
    FileHandler h = null;
    try {
      h = new FileHandler("relationsAdded.txt");
    } catch (SecurityException e1) {
      e1.printStackTrace();
    } catch (IOException e1) {
      e1.printStackTrace();
    }
    h.setFormatter(new OneLineFormatter());
    h.setFilter(new Filter() {

      /* (non-Javadoc)
       * @see java.util.logging.Filter#isLoggable(java.util.logging.LogRecord)
       */
      public boolean isLoggable(LogRecord record) {
        if ((record.getSourceClassName().equals(BioPAXL32KGML.class.getName()) || record
            .getLoggerName().equals(BioPAXL32KGML.class.getName()))
            && record.getLevel().equals(Level.INFO)
            && record.getSourceMethodName().equals("addRelationsToPathway")) {

          return true;
        }
        return false;
      }
    });
    LogUtil.addHandler(h, LogUtil.getInitializedPackages());

    Model m = BioPAX2KGML.getModel(informationFile);    
    if(m!=null){
      
      // For direct translations 2 SBML-qual get and adjust KEGGtranslator
      KEGGtranslator<?> t = null;
      if (createInstantSBML) {
        t = BatchKEGGtranslator.getTranslator(Format.SBML_QUAL, Translator.getManager());
        
        // Adjust options for p2m
        Translator.adjustForPath2Models();
        
        // Pretend to be an original KEGGtranslator instance 
        Translator copyValues = new Translator(null);
        System.setProperty("app.name", copyValues.getAppName());
        System.setProperty("app.version", copyValues.getVersionNumber());
        if (copyValues.getCitation(true)!=null) {
          System.setProperty("app.citation.html", copyValues.getCitation(true));
        }
        if (copyValues.getCitation(false)!=null) {
          System.setProperty("app.citation", copyValues.getCitation(false));
        }
        
      }
      
      for (String filename : fileList) {
        List<Pathway> pathways = null;
        String pw = Utils.ensureSlash(sourceDir) + filename;
        log.info("Parsing pathway '" + pw + "'.");
        try {
          pathways = KeggParser.parse(pw);
        } catch (Exception e) {
          log.log(Level.SEVERE, "Parsing of KEGG pathway '" + pw + "' was not successful.");
        }
        BioPAXL32KGML bp2k = new BioPAXL32KGML();
        if (pathways != null && !pathways.isEmpty()) {
          for (Pathway p : pathways) {
            bp2k.addRelationsToPathway(p, m);
            if (bp2k.getNewAddedRelations()>0 || bp2k.getAddedSubTypes()>0){
              p.setAdditionalText("Pathway information was augmented with BioCarta information " +
              "(2010-08-10)");
              String fn = filename.replace(".xml", "_extended.xml");
              KGMLWriter.writeKGML(p, Utils.ensureSlash(destDir) + fn, writeEntryExtended);
              // As of today, the KGML Pathway data structure is NOT SERIALIZABLE.
              //SerializableTools.saveObject(filename+".dat", p);
              
              // Convert directly to Qual
              if (t!=null) {
                t.translate(p, Utils.ensureSlash(destDir) + filename.replace(".xml", ".sbml.xml"));
              }
            }
          }  
        } else {
          log.severe("Pathways are null, application exits.");
          System.exit(1);
        }
            
      }  
    } else {
      log.log(Level.SEVERE, "Could not continue, because the model is null.");
      System.exit(1);
    }
  }

  
  /**
   * @param args
   */
  public static void main(String[] args) {
    LogUtil.initializeLogging(Level.INFO);
    FileHandler h = null;
    try {
      h = new FileHandler("log.txt");
    } catch (SecurityException e1) {
      e1.printStackTrace();
    } catch (IOException e1) {
      e1.printStackTrace();
    }
    h.setFormatter(new OneLineFormatter());    
    LogUtil.addHandler(h, LogUtil.getInitializedPackages());
    
    FileHandler h2 = null;
    try {
      h2 = new FileHandler("unificationRefs.txt");
    } catch (SecurityException e1) {
      e1.printStackTrace();
    } catch (IOException e1) {
      e1.printStackTrace();
    }
    h2.setFormatter(new OneLineFormatter());
    h2.setFilter(new Filter() {

      /* (non-Javadoc)
       * @see java.util.logging.Filter#isLoggable(java.util.logging.LogRecord)
       */
      public boolean isLoggable(LogRecord record) {
        if ((record.getSourceClassName().equals(Path2Models_KGML_BioCarta_Extender.class.getName()) || record
            .getLoggerName().equals(Path2Models_KGML_BioCarta_Extender.class.getName()))
            && record.getLevel().equals(Level.INFO)) {

          return true;
        }
        return false;
      }
    });
    LogUtil.addHandler(h, LogUtil.getInitializedPackages());
    
  
    List<String> fileList = new ArrayList<String>();
   
    // augment kgmls with BioCarta information
    String KGMLinputDir = "W:/non-metabolic/organisms/hsa"; // "C:/Tools/hsa/";
    KGMLinputDir = "C:/Tools/hsa/";
    String outputFolder = "V:/";
    outputFolder = "C:/Tools/p2m";
    SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd"); 
    
    DirectoryParser d = new DirectoryParser(KGMLinputDir, "xml");
    while (d.hasNext()) {
      fileList.add(d.next());
    }
    
    Path2Models_KGML_BioCarta_Extender.createExtendedKGML(fileList, Utils.ensureSlash(outputFolder) + "BioCarta.bp3_utf8.owl", 
      KGMLinputDir, 
      String.format("%sPath2models/%s-KEGG_NON-METABOLIC_Extended", Utils.ensureSlash(outputFolder), sdf.format(new Date())), true, true);
    
    
  }

}
