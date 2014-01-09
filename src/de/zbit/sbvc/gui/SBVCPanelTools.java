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
package de.zbit.sbvc.gui;

import java.awt.event.ActionListener;
import java.io.File;
import java.net.MalformedURLException;

import de.zbit.graph.gui.TranslatorGraphLayerPanel;
import de.zbit.graph.gui.TranslatorPanel;
import de.zbit.graph.gui.options.GraphBackgroundImageProvider;
import de.zbit.graph.gui.options.TranslatorPanelOptions;
import de.zbit.gui.GUITools;
import de.zbit.kegg.ext.KEGGTranslatorPanelOptions;
import de.zbit.kegg.gui.TranslatorBioPAXPanel;
import de.zbit.kegg.gui.TranslatorGraphPanel;
import de.zbit.kegg.gui.TranslatorSBGNPanel;
import de.zbit.kegg.gui.TranslatorSBMLPanel;
import de.zbit.kegg.gui.TranslatorUI;
import de.zbit.kegg.io.KEGGtranslatorIOOptions;
import de.zbit.kegg.io.KEGGtranslatorIOOptions.Format;
import de.zbit.util.prefs.SBPreferences;

/**
 * @author Finja B&uuml;chel
 * @version $Rev$
 */
public class SBVCPanelTools {

  
  /**
   * 
   * @param inputFile
   * @param outputFormat
   * @param translationResult
   */
  public static TranslatorPanel<?> createPanel(final File inputFile, final Format outputFormat, ActionListener translationResult) {
    TranslatorPanel<?> panel = null;
    
    
    de.zbit.kegg.io.KEGGtranslatorIOOptions.Format keggFormat = 
      KEGGtranslatorIOOptions.Format.valueOf(outputFormat.toString());

    
    switch (outputFormat) {
      case SBML: case SBML_QUAL: case SBML_CORE_AND_QUAL: case SBML_L2V4: case SBML_L3V1: /*case LaTeX: */
        panel = new TranslatorSBMLPanel(inputFile, keggFormat, translationResult);
        break;
        
      case GraphML: case GML: case JPG: case GIF: case YGF: case TGF: /* case SVG:*/
        panel = new TranslatorGraphPanel(inputFile, keggFormat, translationResult);
        break;
        
      case SBGN:
        panel = new TranslatorSBGNPanel(inputFile, translationResult);
        break;
        
      case BioPAX_level2: case BioPAX_level3: case SIF:
        panel = new TranslatorBioPAXPanel(inputFile, keggFormat, translationResult);
        break;
        
      default:
        GUITools.showErrorMessage(null, "Unknown output Format: '" + outputFormat + "'.");
    } 
    
    if (panel!=null && (panel instanceof TranslatorGraphLayerPanel)) {
      setupBackgroundImage((TranslatorGraphLayerPanel<?>) panel);
    }
    return panel;
  }
  
  public static TranslatorPanel<?> createPanel(final String pathwayID, final Format outputFormat, ActionListener translationResult) {
    TranslatorPanel<?> panel = null;
    
    de.zbit.kegg.io.KEGGtranslatorIOOptions.Format keggFormat = 
      KEGGtranslatorIOOptions.Format.valueOf(outputFormat.toString());
    
    switch (outputFormat) {
      case SBML: case SBML_QUAL: case SBML_CORE_AND_QUAL: case SBML_L2V4: case SBML_L3V1: /*case LaTeX: */
        panel = new TranslatorSBMLPanel(pathwayID, keggFormat, translationResult);
        break;
        
      case GraphML: case GML: case JPG: case GIF: case YGF: case TGF: /* case SVG:*/
        panel = new TranslatorGraphPanel(pathwayID, keggFormat, translationResult);
        break;
        
      case SBGN:
        panel = new TranslatorSBGNPanel(pathwayID, translationResult);
        break;
        
      case BioPAX_level2: case BioPAX_level3: case SIF:
        panel = new TranslatorBioPAXPanel(pathwayID, keggFormat, translationResult);
        break;
        
      default:
        GUITools.showErrorMessage(null, "Unknwon output Format: '" + outputFormat + "'.");
        return null;
    }
    
    if (panel!=null && (panel instanceof TranslatorGraphLayerPanel)) {
      setupBackgroundImage((TranslatorGraphLayerPanel<?>) panel);
    }
    return panel;
  }
  
  /**
   * Setup the background image as set in the preferences
   * @param pane the pane to add the background image
   * @param translator the translator used for translation
   * @param prefs might be null, else, prefs object for {@link TranslatorPanelOptions}
   * @throws MalformedURLException
   */
  public static void setupBackgroundImage(TranslatorGraphLayerPanel<?> panel) {
    SBPreferences prefs = SBPreferences.getPreferencesFor(KEGGTranslatorPanelOptions.class);
    GraphBackgroundImageProvider provider = null;
    
    if (KEGGTranslatorPanelOptions.SHOW_LOGO_IN_GRAPH_BACKGROUND.getValue(prefs)) {
      provider = GraphBackgroundImageProvider.Factory.createStaticImageProvider(TranslatorUI.getWatermarkLogoResource());
      
    } else if (KEGGTranslatorPanelOptions.SHOW_KEGG_PICTURE_IN_GRAPH_BACKGROUND.getValue(prefs)) {
        Integer brighten = (KEGGTranslatorPanelOptions.BRIGHTEN_KEGG_BACKGROUND_IMAGE.getValue(prefs));
        if (brighten==null || brighten<0) brighten = 0;
        boolean greyscale = (KEGGTranslatorPanelOptions.GREYSCALE_KEGG_BACKGROUND_IMAGE.getValue(prefs));
        provider = GraphBackgroundImageProvider.Factory.createDynamicTranslatorImageProvider(brighten, greyscale);
    }
    
    // Setup the provider
    panel.setBackgroundImageProvider(provider);
  }
  
}
