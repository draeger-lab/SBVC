/*
 * $Id: SBVC.java 190 2015-10-20 11:17:21Z roemer $
 * $URL: https://rarepos.cs.uni-tuebingen.de/svn-path/SBVC/trunk/src/de/zbit/sbvc/SBVC.java $
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

import static de.zbit.util.Utils.getMessage;

import java.awt.Window;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import de.zbit.AppConf;
import de.zbit.Launcher;
import de.zbit.biopax.BioPAX2KGML;
import de.zbit.cache.InfoManagement;
import de.zbit.gui.GUIOptions;
import de.zbit.io.FileTools;
import de.zbit.kegg.KEGGtranslatorOptions;
import de.zbit.kegg.KEGGtranslatorOptions.NODE_NAMING;
import de.zbit.kegg.KGMLWriter;
import de.zbit.kegg.Translator;
import de.zbit.kegg.api.cache.KeggInfoManagement;
import de.zbit.kegg.gui.TranslatorUI;
import de.zbit.kegg.io.AbstractKEGGtranslator;
import de.zbit.kegg.io.KEGG2SBMLqual;
import de.zbit.resources.Resource;
import de.zbit.sbvc.gui.SBVCUI;
import de.zbit.sbvc.io.SBVCIOOptions;
import de.zbit.util.Species;
import de.zbit.util.Utils;
import de.zbit.util.logging.LogUtil;
import de.zbit.util.prefs.KeyProvider;
import de.zbit.util.prefs.SBProperties;

/**
 * @author Finja B&uuml;chel
 * @author Clemens Wrzodek
 * @version $Rev: 190 $
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

	public SBVC() {
		// this(new String[0]);
	}

	public SBVC(String[] args) {
		super(args);
	}

	@Override
	public void commandLineMode(AppConf appConf) {
		SBProperties props = appConf.getCmdArgs();
		String folderName = SBVCIOOptions.OUTPUT.getValue(props).getPath();
		String input = SBVCIOOptions.INPUT.getValue(props).getPath();
		Integer speciesInput = SBVCIOOptions.SPECIES.getValue(props);

		Species species = null;
		if (speciesInput != null) {
			int tax = speciesInput.intValue();
			if (tax == 9606) {
				species = new Species("Homo sapiens", "_HUMAN", "human", "hsa",9606);
			} else if (tax == 10090) {
				species = new Species("Mus musculus", "_MOUSE", "mouse", "mmu", 10090);
			} else {
				species = Species.getSpeciesWithTaxonomyIDInList(speciesInput, new BufferedReader(
						new InputStreamReader(Resource.class.getResourceAsStream("speclist.txt"))), false);
			}
		}

		if(species == null){
			log.warning("Taxonomy id '" + speciesInput + "' could not be assigned to a species.");
		} else {
			log.info("Entered species '" + species.getCommonName() + "' was identified.");
		}

		convertBioPAXToSBML(input, folderName, species);
	}

	/**
	 * 
	 * @param input <code>BioPAX</code> file.
	 * @param outputFolderName Result folder in which the created file should be put.
	 * @param splitMode {@code true} if separate files for each pathway should
	 * be created. {@code false} to create a single output file.
	 */
	public void convertBioPAXToSBML(String input, String outputFolderName, Species species) {
		// getting the KEGG Pathways of the model
		Collection<de.zbit.kegg.parser.pathway.Pathway> keggPWs =
				BioPAX2KGML.createPathwaysFromModel(input, species);

		// translation to sbml
		KEGG2SBMLqual k2s = null;
		if (new File(Translator.cacheFileName).exists()
				&& new File(Translator.cacheFileName).length() > 1) {
			KeggInfoManagement manager = null;
			try {
				manager = (KeggInfoManagement) InfoManagement.loadFromFilesystem(Translator.cacheFileName);
				k2s = new KEGG2SBMLqual(manager);
			} catch (IOException e) {
				k2s = new KEGG2SBMLqual();
			}
		} else {
			k2s = new KEGG2SBMLqual();
		}

		//    // original version
		//    // necessary that both reactions and relations are written to the file
		//    k2s.setConsiderReactions(true);
		//    k2s.setAddCellDesignerAnnots(false);
		//    k2s.setNameToAssign(NODE_NAMING.INTELLIGENT);
		//    k2s.setRemoveOrphans(false);
		//    k2s.setAutocompleteReactions(false);
		//    k2s.setRemoveWhiteNodes(false);
		//    k2s.setShowFormulaForCompounds(false);
		//    k2s.setRemovePathwayReferences(false);
		//    k2s.setAddLayoutExtension(false);
		//    k2s.setCheckAtomBalance(false);

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
		k2s.setUseGroupsExtension(false);
		k2s.setCheckAtomBalance(false);

		// Translate and write output
		if (keggPWs.size()==1 && !new File(outputFolderName).isDirectory()) {
			// Single in and single out file, no split mode.
			de.zbit.kegg.parser.pathway.Pathway p = keggPWs.iterator().next();
			String title = FileTools.removeFileExtension(outputFolderName);
			if (title.toLowerCase().endsWith(".sbml")) {
				title = title.substring(0, title.length()-5);
			}
			p.setTitle(title);
			k2s.translate(p, outputFolderName);

		} else {
			// We had multiple biopax pathway objects in input
			// SBML does not permit multiple models => write one file per model
			for (de.zbit.kegg.parser.pathway.Pathway p : keggPWs) {
				k2s.translate(p, Utils.ensureSlash(outputFolderName) + KGMLWriter.createFileName(p));
			}

		}

		// Remember already queried objects (save cache)
		if (AbstractKEGGtranslator.getKeggInfoManager().hasChanged()) {
			KeggInfoManagement.saveToFilesystem(Translator.cacheFileName, AbstractKEGGtranslator.getKeggInfoManager());
		}
	}

	@Override
	public List<Class<? extends KeyProvider>> getCmdLineOptions() {
		/*
		 * XXX: we can add further IOOptions later, for example those of KEGG translator
		 * TODO: getPublicationXref() in KEGG2BioPAX - find a way to insert other publications?
		 */
		List<Class<? extends KeyProvider>> configList = new ArrayList<Class<? extends KeyProvider>>(3);
		configList.add(SBVCIOOptions.class);
		configList.add(GUIOptions.class);
		configList.add(KEGGtranslatorOptions.class);
		return configList;
	}

	@Override
	public List<Class<? extends KeyProvider>> getInteractiveOptions() {
		// TODO Auto-generated method stub
		return null;
	}

	/* (non-Javadoc)
	 * @see de.zbit.Launcher#getAppName()
	 */
	@Override
	public String getAppName() {
		//    return "System Biology Visualizer and Converter";
		return "BioPAX2SBML"; // XXX: This is for the BioPAX2SBML converter
	}

	/* (non-Javadoc)
	 * @see de.zbit.Launcher#getCitation(boolean)
	 */
	@Override
	public String getCitation(boolean HTMLstyle) {
		// TODO Return a citation string, as soon as this is published.
		return super.getCitation(HTMLstyle);
	}

	/* (non-Javadoc)
	 * @see de.zbit.Launcher#addCopyrightToSplashScreen()
	 */
	@Override
	protected boolean addCopyrightToSplashScreen() {
		return false;
	}

	@Override
	public URL getURLlicenseFile() {
		URL url = null;
		try {
			url = new URL("http://www.gnu.org/licenses/lgpl-3.0-standalone.html");
		} catch (MalformedURLException exc) {
			log.log(Level.FINE, getMessage(exc), exc);
		}
		return url;
	}

	@Override
	public URL getURLOnlineUpdate() {
		try {
			return new URL("http://www.cogsys.cs.uni-tuebingen.de/software/SBVC/downloads/");
		} catch (MalformedURLException e) {
			log.log(Level.FINE, e.getLocalizedMessage(), e);
		}
		return null;
	}

	@Override
	public String getVersionNumber() {
		return "1.0.2";
	}

	@Override
	public short getYearOfProgramRelease() {
		return (short) 2016;
	}

	/* (non-Javadoc)
	 * @see de.zbit.Launcher#getOrganization()
	 */
	@Override
	public String getOrganization() {
		return "by the individual authors";
	}

	/* (non-Javadoc)
	 * @see de.zbit.Launcher#getInstitute()
	 */
	@Override
	public String getInstitute() {
		return "See https://github.com/cogsys-tuebingen/SBVC/";
	}

	@Override
	public short getYearWhenProjectWasStarted() {
		return (short) 2012;
	}

	@Override
	public Window initGUI(AppConf appConf) {
		return new SBVCUI(appConf);
	}

	/**
	 * 
	 * @param args
	 */
	public static void main(String args[]) {
		LogUtil.initializeLogging(Level.INFO);
		// "Merge" other applications with this one
		// Must be done first, because option defaults are changed
		integrateIntoKEGGtranslator();
		GUIOptions.GUI.setDefaultValue(Boolean.FALSE);
		// Make an instance of this application
		new SBVC(args);
	}

	/**
	 * This method changes some default values of KEGGtranslator,
	 * option visibilities, logos, etc. to look like it would be
	 * the {@link SBVC} application.
	 */
	public static void integrateIntoKEGGtranslator() {
		// Set default values for KEGGtranslator
		KEGGtranslatorOptions.REMOVE_ORPHANS.setDefaultValue(false);
		KEGGtranslatorOptions.REMOVE_WHITE_GENE_NODES.setDefaultValue(false);
		KEGGtranslatorOptions.AUTOCOMPLETE_REACTIONS.setDefaultValue(false);
		KEGGtranslatorOptions.AUTOCOMPLETE_REACTIONS.setVisible(false);

		KEGGtranslatorOptions.REMOVE_PATHWAY_REFERENCES.setDefaultValue(false);

		// TODO: Modify all options further to fit the needs of SBVC.
		// Also set the real values (not only the default values) e.g. for autocomplete r.

		TranslatorUI.watermarkLogoResource = "img/logo.png";
		// TODO: BackgroundImageProviders? on Graphs...
		//TranslatorGraphLayerPanel.optionClass = KEGGTranslatorPanelOptions.class;
	}

}
