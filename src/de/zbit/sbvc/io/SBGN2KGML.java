/*
 * $Id: package-info.java 139 2012-08-24 15:19:03Z draeger $
 * $URL: https://rarepos.cs.uni-tuebingen.de/svn-path/SBVC/trunk/src/de/zbit/sbvc/io/package-info.java $
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
 * Copyright (C) 2012-2013 by the University of Tuebingen, Germany.
 *
 * SBVC is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation. A copy of the license
 * agreement is provided in the file named "LICENSE.txt" included with
 * this software distribution and also available online as
 * <http://www.gnu.org/licenses/lgpl-3.0-standalone.html>.
 * ---------------------------------------------------------------------
 */
package de.zbit.sbvc.io;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.xml.bind.JAXBException;

import org.sbgn.SbgnUtil;
import org.sbgn.bindings.Arc;
import org.sbgn.bindings.Glyph;
import org.sbgn.bindings.Port;
import org.sbgn.bindings.Sbgn;
import org.xml.sax.SAXException;

import de.zbit.graph.io.def.SBGNProperties.ArcType;
import de.zbit.graph.io.def.SBGNProperties.GlyphType;
import de.zbit.kegg.KGMLWriter;
import de.zbit.kegg.parser.pathway.Entry;
import de.zbit.kegg.parser.pathway.EntryType;
import de.zbit.kegg.parser.pathway.Graphics;
import de.zbit.kegg.parser.pathway.Pathway;
import de.zbit.kegg.parser.pathway.Reaction;
import de.zbit.kegg.parser.pathway.ReactionComponent;
import de.zbit.kegg.parser.pathway.ReactionType;
import de.zbit.kegg.parser.pathway.Relation;
import de.zbit.kegg.parser.pathway.RelationType;
import de.zbit.kegg.parser.pathway.SubType;
import de.zbit.kegg.parser.pathway.ext.EntryExtended;
import de.zbit.sbvc.io.helper.SBGN2KGMLHelper;

/**
 * Class for converting {@link Sbgn} into a {@link Pathway}
 *
 * @author 	Manuel Ruff
 * @date 	2012-06-03
 * @version $Rev: 165$
 * @since	Revision 99
 *
 */
public class SBGN2KGML {

	public static final Logger log = Logger.getLogger(SBGN2KGML.class.getName());						// logger for errors, warnings, etc.
	
	public SBGN2KGMLHelper helper = new SBGN2KGMLHelper();												// helper class to get the proper EntryType and ArcType
	
	
	HashMap<String, Entry> entryLookup = new HashMap<String, Entry>();									// mapping from the glyph id onto the entries
	HashMap<String, String> clonemarkerLookup = new HashMap<String, String>();							// mapping from the clonemarker id onto the kegg id
	HashMap<String, Glyph> processGlyphLookup = new HashMap<String, Glyph>();							// mapping from the process glyph id onto the glyph itself
	
	HashMap<String, ArrayList<String>> portSubstrates = new HashMap<String, ArrayList<String>>();		// mapping from the port id onto the corresponding substrate ids
	HashMap<String, ArrayList<String>> portProducts = new HashMap<String, ArrayList<String>>();			// mapping from the port id onto the corresponding product ids
	HashMap<String, ArrayList<String>> reactionModifiers = new HashMap<String, ArrayList<String>>();	// mapping from the process glyph id onto the reactionmodifiers
	
	ArrayList<Arc> arcsWithoutPorts = new ArrayList<Arc>();												// storage of all arcs without ports (relations)
	
	private int entryID = 0;																			// incremental number for all entries
	private int reactionID = 0;																			// incremental number for all reactions
	private int skipCount = 0;																			// counts the clonemarker occurences to guarantee ongoing names
	
	private boolean considerGlyphs = true;																
	private boolean considerArcs = true;

	/**
	 * Method for reading a {@link Sbgn}-file
	 * @param filename		the name of the file as {@link String}
	 * @return {@link Sbgn}	returns a {@link Sbgn}-object
	 * @exception JAXBException
	 */
	public Sbgn read(String filename) {

		// initialize the objects
		File f = new File(filename);
		Sbgn sbgn = null;
		
		// try reading the file
		try {
			sbgn = SbgnUtil.readFromFile(f);
			
		} catch (JAXBException e) {
			log.log(Level.SEVERE, String.format("Couldn't read in the the file: %s", f.getAbsoluteFile()), e);
		}
		
		return sbgn;
	}

	/**
	 * Method for validating a {@link Sbgn}-file
	 * @param filename	the name of the file as {@link String}
	 * @return {@link boolean}	<ul>
	 * 								<li> {@code true} if the file is valid</li>
	 * 								<li> {@code false} otherwise</li>
	 * 							</ul>
	 * @exception JAXBException
	 * @exception SAXException
	 * @exception IOException
	 */
	public boolean validateSBGN(String filename) {
		
		// initialize the objects
		File f = new File(filename);
		boolean isValid = false;
		
		// try validating the file
		try {
			isValid = SbgnUtil.isValid(f);
		} catch (JAXBException e) {
			log.log(Level.SEVERE, String.format("Couldn't read in the the file: %s", f.getAbsoluteFile()), e);
		} catch (SAXException e) {
			log.log(Level.SEVERE, String.format("Invalid element: %s", f.getAbsoluteFile()), e);
		} catch (IOException e) {
			log.log(Level.SEVERE, String.format("Couldn't create or find the file: %s", f.getAbsoluteFile()), e);
		}
		
		return isValid;
	}

	/**
	 * Method for translating a {@link Sbgn}-file into a {@link Pathway}
	 * @param {@link Sbgn}		(see {@link SBGN2KGML#read(String)})
	 * @return {@link Pathway}	(see {@link KGMLWriter})
	 */
	protected Pathway translate(Sbgn sbgn) {

		// create a new pathway
		Pathway p = new Pathway("unknown", "unknown", 10000, "unknown_title");

		// translate the glyphs
		if(considerGlyphs)
			handleAllGlyphs(sbgn, p);
		
		// translate the arcs
		if(considerArcs)
			handleAllArcs(sbgn, p);
		
		
		// add entries to pathway
		addAllEntriesToPathway(p);
		
		return p;
	}
	
	/**
	 * Method for adding all {@link Entry}'s to the {@link Pathway}
	 * @param p	{@link Pathway}
	 */
	private void addAllEntriesToPathway(Pathway p) {

		// create an iterator
		Iterator<Entry> entryIterator = entryLookup.values().iterator();
		
		// add the entries to the pathway
		while(entryIterator.hasNext()){
			Entry e = entryIterator.next();
			p.addEntry(e);
		}
		
	}

	/**
	 * Method to convert all {@link Glyph}'s into {@link Pathway#getEntries()}
	 * @param sbgn	{@link Sbgn}
	 * @param p		{@link Pathway}
	 */
	protected void handleAllGlyphs(Sbgn sbgn, Pathway p){
		
		for (Glyph g : sbgn.getMap().getGlyph()) {
			
			// create a mapping for the gylpsh and store the informations in the hashmaps
			createMappingForProcessGlyphsOnly(g);
			
			// no label means process glyph
			if(g.getLabel() != null){
				
				// create a new entry with an ongoing id
				Entry e = new EntryExtended(p, ++entryID, "unknown_entry:"+(entryID-skipCount));

				// convert the GlyphType into an EntryType
				e.setType(helper.getEntryTypeFromGlyphType(GlyphType.valueOfString(g.getClazz())));
				
				// create the graphics
				Graphics gr = null;
				
				// determine the graphic type of the entry
				if (e.getType().compareTo(EntryType.compound) == 0)
					gr = Graphics.createGraphicsForCompound(g.getLabel().getText());
				else if (e.getType().compareTo(EntryType.map) == 0)
					gr = Graphics.createGraphicsForPathwayReference(g.getLabel().getText());
				else if (e.getType().compareTo(EntryType.group) == 0)
				  gr = Graphics.createGraphicsForGroupOrComplex(g.getLabel().getText());
				else if (e.getType().compareTo(EntryType.gene) == 0)
				  gr = Graphics.createGraphicsForProtein(g.getLabel().getText());
				else
				  gr = new Graphics(g.getLabel().getText());
				
				if(gr != null){
					gr.setDefaults(e.getType());
					gr.setX((int) g.getBbox().getX());
					gr.setY((int) g.getBbox().getY());
					gr.setHeight((int) g.getBbox().getH());
					gr.setWidth((int) g.getBbox().getW());
				}
	
				// add the graphics to the entry
				e.addGraphics(gr);
				
				// map the glyph onto the entries and the glyphs
				e.setCustom(g);
				entryLookup.put(g.getId(), e);
				
				// same entries have the same name
				// TODO: maybe decrease the number if a clone is found to get ongoing numbering!
				if(g.getClone() != null) {
					if(clonemarkerLookup.containsKey(g.getLabel().getText())) {
						e.setName(clonemarkerLookup.get(g.getLabel().getText()));
						skipCount++;
					} else {
						clonemarkerLookup.put(g.getLabel().getText(), e.getName());
					}
				} else {
					clonemarkerLookup.put(g.getLabel().getText(), e.getName());
				}
			}
		}
		
	}
	
	/**
	 * Method to convert all the {@link Arc}'s into {@link Pathway#getReactions()} and {@link Pathway#getRelations()}
	 * @param sbgn	{@link Sbgn}
	 * @param p		{@link Pathway}
	 */
	protected void handleAllArcs(Sbgn sbgn, Pathway p){
		
		// create a mapping a the arcs and ports and store the informations in the hashmaps
		createMappingForAllPortsAndArcs(sbgn);
		
		/////////////////////////////////////
		// create edges for the relations //
		///////////////////////////////////
		for(Arc awp : arcsWithoutPorts) {
			
			// get source and target glyphs
			Glyph sourceGlyph = (Glyph) awp.getSource();
			Glyph targetGlyph = (Glyph) awp.getTarget();
			
			// get the corresponding entries
			int sourceEntryID = entryLookup.get(sourceGlyph.getId()).getId();
			int targetEntryID = entryLookup.get(targetGlyph.getId()).getId();
			
			// create a relation
			// TODO: maybe change the relation type accordingly
			Relation relation = new Relation(sourceEntryID, targetEntryID, RelationType.other);
			
			// get the arctype
			ArcType atype = ArcType.valueOfString(awp.getClazz());
			
			// look up the subtype
			SubType subtype = helper.getSubTypeFromArcType(atype);
			
			// set the subtype only if there is a subtype
			if(subtype != null){
				relation.addSubtype(subtype);
			}
			
			// add the relation to the pathway
			p.addRelation(relation);
			
		}

		/////////////////////////////////////
		// create edges for the reactions //
		///////////////////////////////////
		
		// create an iterator for all process glyphs
		Iterator<Glyph> processGlyphIterator = processGlyphLookup.values().iterator();
		
		while(processGlyphIterator.hasNext()){
			
			// initialize
			Glyph processGlyph = processGlyphIterator.next();
			ArrayList<String> substrates = new ArrayList<String>();
			ArrayList<String> products = new ArrayList<String>();
			
			// collect the ids for the ports of the process glyph and determine if the reactants are substrates or products
			for(Port port : processGlyph.getPort()){
				if(portSubstrates.containsKey(port.getId()))
					substrates.add(port.getId());
				if(portProducts.containsKey(port.getId()))
					products.add(port.getId());
			}
			
			// create a new reaction
			Reaction r = null;
			
			// if there are no substrates or products  it is a reversible reaction
			if(substrates.isEmpty() || products.isEmpty()){
				r = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.reversible);
				
				// we need to set the substrates and products even if there are none of one
				if(!products.isEmpty()){
					for(int i = 0; i < products.size(); i++){
						for(String s : portProducts.get(products.get(i))){
							if(i == 0)
								r.addSubstrate(new ReactionComponent(entryLookup.get(s)));
							else
								r.addProduct(new ReactionComponent(entryLookup.get(s)));
						}
					}
				}
				if(!substrates.isEmpty()){
					for(int i = 0; i < substrates.size(); i++){
						for(String s : portSubstrates.get(substrates.get(i))){
							if(i == 0)
								r.addSubstrate(new ReactionComponent(entryLookup.get(s)));
							else
								r.addProduct(new ReactionComponent(entryLookup.get(s)));
						}
					}
				}
			}
			// otherwise it is an irreversible reaction
			else{
				r = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.irreversible);
			
				// add the products and substrates to the reaction
				for(int i = 0; i < products.size(); i++){
					for(String s : portProducts.get(products.get(i))){
						r.addProduct(new ReactionComponent(entryLookup.get(s)));
					}
				}
				
				for(int i = 0; i < substrates.size(); i++){
					for(String s : portSubstrates.get(substrates.get(i))){
						r.addSubstrate(new ReactionComponent(entryLookup.get(s)));
					}
				}
			}
			
			// if there are reaction modifiers i.e. enzymes update the reaction where they are involved
			if(reactionModifiers.containsKey(processGlyph.getId())){
				ArrayList<String> rmIDs = reactionModifiers.get(processGlyph.getId());
				for(String rm : rmIDs){
					Entry e = entryLookup.get(rm);
					if(e.isSetReaction())
						e.setReaction(e.getReactionString() + " " + r.getName());
					else
						entryLookup.get(rm).setReaction(r.getName());
				}
			}
			
			// add the reaction to the pathway
			p.addReaction(r);
		}
			
	}
	
	/**
	 * Method to create a mapping for all {@link Port}'s onto their linked {@link Glyph}'s
	 * and gather all {@link Arc}'s without any {@link Port}'s
	 * @param sbgn	{@link Sbgn}
	 */
	protected void createMappingForAllPortsAndArcs(Sbgn sbgn) {
		
		for (Arc arc : sbgn.getMap().getArc()) {
		
			// initialize
			ArrayList<String> collectedGlyphNamesForOnePort = new ArrayList<String>();
			ArrayList<String> collectedReactionModifiers = new ArrayList<String>();
			Port port = null;
			boolean isProduct = false;
			boolean isSubstrate = false;
			
			// check if source or target is a port and if its already gathered
			// update the glyphname list
			if(arc.getSource().getClass().equals(Port.class)) {
				port = (Port) arc.getSource();
				isProduct = true;
				if(portProducts.containsKey(port.getId()))
						collectedGlyphNamesForOnePort = portProducts.get(port.getId());
				collectedGlyphNamesForOnePort.add(((Glyph) arc.getTarget()).getId());
			}
			
			if(arc.getTarget().getClass().equals(Port.class)) {
				port = (Port) arc.getTarget();
				isSubstrate = true;
				if(portSubstrates.containsKey(port.getId()))
						collectedGlyphNamesForOnePort = portSubstrates.get(port.getId());
				collectedGlyphNamesForOnePort.add(((Glyph) arc.getSource()).getId());
			}
			
			if(isProduct)
				portProducts.put(port.getId(), collectedGlyphNamesForOnePort);
			
			if(isSubstrate)
				portSubstrates.put(port.getId(), collectedGlyphNamesForOnePort);
			
			// glyphs without ports are relations and stored in arcsWithoutPorts
			if(port == null){
				if(((Glyph) arc.getTarget()).getPort().isEmpty() && ((Glyph) arc.getSource()).getPort().isEmpty()){
					arcsWithoutPorts.add(arc);
				} else if(!((Glyph) arc.getTarget()).getPort().isEmpty()){
					if(reactionModifiers.containsKey(((Glyph) arc.getTarget()).getId()))
							collectedReactionModifiers = reactionModifiers.get(((Glyph) arc.getTarget()).getId());
					collectedReactionModifiers.add(((Glyph) arc.getSource()).getId());
					reactionModifiers.put(((Glyph) arc.getTarget()).getId(), collectedReactionModifiers);
				}
			}
				
		}
	}
	
	/**
	 * Method for writing a {@link Pathway} into a file
	 * @param p			{@link Pathway}
	 * @param filename	{@link String}
	 */
	protected void saveToFile(Pathway p, String filename){
		KGMLWriter.writeKGML(p, filename, false);
	}
	
	/**
	 * Method to create a mapping only for a Process-{@link Glyph}'s
	 * @param g	{@link Glyph}
	 */
	protected void createMappingForProcessGlyphsOnly(Glyph g){
		// glyphs with portlists are process glyphs
		if(!g.getPort().isEmpty())
			processGlyphLookup.put(g.getId(), g);		
	}
	
	
	public static void main(String args[]) {
//		String filename = "glycolysis.sbgn";
//		String filename = "mapk_cascade.sbgn";
//		String filename = "files/SBGNExamples/test.sbgn";
		String filename = "test.sbgn";
		SBGN2KGML sbgn2kgml = new SBGN2KGML();
		Sbgn sbgn = sbgn2kgml.read(filename);
		
		System.out.println(sbgn2kgml.validateSBGN(filename));
		
		Pathway p = sbgn2kgml.translate(sbgn);
		sbgn2kgml.saveToFile(p, "PathwayTest.xml");
	}

}
