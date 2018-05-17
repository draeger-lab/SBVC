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
package de.zbit.sbvc.io;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import de.zbit.kegg.KGMLWriter;
import de.zbit.kegg.parser.pathway.Entry;
import de.zbit.kegg.parser.pathway.EntryType;
import de.zbit.kegg.parser.pathway.Graphics;
import de.zbit.kegg.parser.pathway.Pathway;
import de.zbit.kegg.parser.pathway.Reaction;
import de.zbit.kegg.parser.pathway.ReactionComponent;
import de.zbit.kegg.parser.pathway.ReactionType;
import de.zbit.kegg.parser.pathway.Relation;
import de.zbit.kegg.parser.pathway.SubType;
import de.zbit.kegg.parser.pathway.ext.EntryExtended;
import de.zbit.sbvc.io.helper.SIFPathway;
import de.zbit.sbvc.io.helper.SIFProperties;
import de.zbit.sbvc.io.helper.SIFProperties.InteractionType;
import de.zbit.sbvc.io.helper.SIFRelation;

/**
 * Class for converting Simple Interaction Files into a {@link Pathway}
 *
 * @author 	Manuel Ruff
 * @date 	2012-08-08
 * @version $Rev: 138$
 * @since	Revision 135
 *
 */
public class SIF2KGML {
	
	public static final Logger log = Logger.getLogger(SIF2KGML.class.getName());					// logger for errors, warnings, etc.
	
	private int entryID = 0;																		// incremental number for all entries
	private int reactionID = 0;																		// incremental number for all reactions
	
	/**
	 * Method for reading in a Simple Interaction File
	 * This Method tries to determine if the SIF-File is tab- or space-separated
	 * @param filename	{@link String}
	 * @return {@Link ArrayList<SifRelation>)
	 * @throws IOException 
	 */
	public static SIFPathway readSIF(String filename) throws IOException {
		
		if(!filename.endsWith(".sif"))
			log.log(Level.SEVERE, String.format("The file %s doesn't end with .sif", filename));
		
		Scanner scanner = null;
		
		// take the filename and set it as Pathwayname
		SIFPathway sif = new SIFPathway(filename.replaceAll("\\.sif", "").replaceAll(".*/", "").replaceAll(".*\\\\", ""));
		
		// manage the SIFRelations in an arraylist
		ArrayList<SIFRelation> relations = new ArrayList<SIFRelation>();

		// try reading in the sif file line by line
		try {
			scanner = new Scanner(new FileReader(filename)).useDelimiter("\n");
			
			// initialize
			String line;
			int counter = 0;
			
			while(scanner.hasNext()) {
				
				// create objects to store the information
				String[] tabsplitted;
				String[] spacesplitted;
				String substrate = null;
				String product = null;
				InteractionType iType = null;
				
				// read
				line = scanner.nextLine();
				counter++;
				tabsplitted = line.split("\t");
				spacesplitted = line.split(" ");
				
				// check if the file fulfills the SIF patterns
				if(tabsplitted.length == 3) {
					substrate = tabsplitted[0];
					product = tabsplitted[2];
					iType = InteractionType.valueOf(tabsplitted[1]);
				} else if (spacesplitted.length == 3) {
					substrate = spacesplitted[0];
					product = spacesplitted[2];
					iType = InteractionType.valueOf(spacesplitted[1]);
				}
				
				// store the information in the relation objects
				if(substrate != null && product != null && iType != null ) {
					SIFRelation relation = new SIFRelation(substrate, iType, product);
					relations.add(relation);
				} else {
					log.log(Level.WARNING, String.format("The line %s in the file %s doesn't have the proper Simple Interaction Format", new Object[]{counter, filename}));
				}
				
			}
		} catch (FileNotFoundException e) {
			log.log(Level.SEVERE, String.format("Couldn't find the file " + filename));
		}
		
		// add the relations to the sif pathway
		sif.setRelations(relations);
		
		// close the scanner
		scanner.close();
		
		return sif;
	}
	
	/**
	 * Method for translating a {@link SIFPathway} into a {@link Pathway}
	 * @param {@link SIFPathway}	(see {@link SIF2KGML#readSIF(String)})
	 * @return {@link Pathway}	(see {@link KGMLWriter})
	 */
	protected Pathway translate(SIFPathway sif) {

		// create a new pathway
		Pathway p = new Pathway("unknown", "unknown", 10000, sif.getName());
		
		// Entry lookup
		HashMap<String, Entry> entryLookup = new HashMap<String, Entry>();
		
		for(SIFRelation relation : sif.getRelations()) {
						
			// initialize entries
			Entry source = null;
			Entry target = null;
			
			// check if entries already exists
			if(entryLookup.containsKey(relation.getSource()))
				source = entryLookup.get(relation.getSource());
			
			if(entryLookup.containsKey(relation.getTarget()))
				target = entryLookup.get(relation.getTarget());
			
			// if not create them
			if(source == null) {
				source = new EntryExtended(p, ++entryID, relation.getSource());
				entryLookup.put(relation.getSource(), source);
			}
		
			if(target == null) {
				target = new EntryExtended(p, ++entryID, relation.getTarget());
				entryLookup.put(relation.getTarget(), target);
			}
			
			// set EntryTypes
			source.setType(SIFProperties.getEntryTypeFromInteractionType(relation.getInteractionType(), true));
			target.setType(SIFProperties.getEntryTypeFromInteractionType(relation.getInteractionType(), false));
			
			// create Graphics
			Graphics sourceGraphic = null;
			Graphics targetGraphic = null;
			
			// create reaction and relation
			Reaction rea = null;
			Relation rel = null;
			
			if(relation.getInteractionType() != null){
				switch(relation.getInteractionType()){
					case pp:
						sourceGraphic = Graphics.createGraphicsForProtein(source.getName());
						targetGraphic = Graphics.createGraphicsForProtein(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						rel.addSubtype(new SubType(SubType.BINDING_ASSOCIATION));
						break;
					case pd:
						sourceGraphic = Graphics.createGraphicsForProtein(source.getName());
						targetGraphic = new Graphics(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						rel.addSubtype(new SubType(SubType.EXPRESSION));
						break;
					case pr:
						sourceGraphic = Graphics.createGraphicsForProtein(source.getName());
						targetGraphic = Graphics.createGraphicsForPathwayReference(source.getName());
						rea = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.irreversible);
						rea.addSubstrate(new ReactionComponent(source));
						rea.addProduct(new ReactionComponent(target));
						break;
					case rc:
						sourceGraphic = Graphics.createGraphicsForPathwayReference(source.getName());
						targetGraphic = Graphics.createGraphicsForCompound(target.getName());
						rea = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.irreversible);
						rea.addSubstrate(new ReactionComponent(source));
						rea.addProduct(new ReactionComponent(target));
						break;
					case cr:
						sourceGraphic = Graphics.createGraphicsForCompound(source.getName());
						targetGraphic = Graphics.createGraphicsForPathwayReference(target.getName());
						rea = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.irreversible);
						rea.addSubstrate(new ReactionComponent(source));
						rea.addProduct(new ReactionComponent(target));
						break;
					case gl:
						sourceGraphic = new Graphics(source.getName());
						targetGraphic = new Graphics(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						rel.addSubtype(new SubType(SubType.MISSING_INTERACTION));
						break;
					case pm:
						sourceGraphic = Graphics.createGraphicsForProtein(source.getName());
						targetGraphic = Graphics.createGraphicsForCompound(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						rel.addSubtype(new SubType(SubType.MISSING_INTERACTION));
						break;
					case mp:
						sourceGraphic = Graphics.createGraphicsForCompound(source.getName());
						targetGraphic = Graphics.createGraphicsForProtein(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						rel.addSubtype(new SubType(SubType.MISSING_INTERACTION));
						break;
					case CO_CONTROL:
						sourceGraphic = new Graphics(source.getName());
						targetGraphic = new Graphics(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						rel.addSubtype(new SubType(SubType.ASSOCIATION));
						break;
					case COMPONENT_OF:
						sourceGraphic = new Graphics(source.getName());
						targetGraphic = new Graphics(target.getName());
						target.addComponent(source.getId());
						target.setType(EntryType.reaction); // not sure if this is correct
						break;
					case IN_SAME_COMPONENT:
						sourceGraphic = new Graphics(source.getName());
						targetGraphic = new Graphics(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						rel.addSubtype(new SubType(SubType.BINDING_ASSOCIATION));
						break;
					case METABOLIC_CATALYSIS:
						sourceGraphic = new Graphics(source.getName());
						targetGraphic = new Graphics(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						rel.addSubtype(new SubType(SubType.MISSING_INTERACTION));
						break;
					case REACTS_WITH:
						sourceGraphic = new Graphics(source.getName());
						targetGraphic = new Graphics(target.getName());
						rea = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.irreversible);
						rea.addSubstrate(new ReactionComponent(source));
						rea.addProduct(new ReactionComponent(target));
						break;
					case SEQUENTIAL_CATALYSIS:
						sourceGraphic = new Graphics(source.getName());
						targetGraphic = new Graphics(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						rel.addSubtype(new SubType(SubType.INDIRECT_EFFECT));
						break;
					case STATE_CHANGE:
						sourceGraphic = new Graphics(source.getName());
						targetGraphic = new Graphics(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						rel.addSubtype(new SubType(SubType.STATE_CHANGE));
						break;
					case INTERACTS_WITH:
						sourceGraphic = new Graphics(source.getName());
						targetGraphic = new Graphics(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						rel.addSubtype(new SubType(SubType.MISSING_INTERACTION));
						break;
					default:
						sourceGraphic = new Graphics(source.getName());
						targetGraphic = new Graphics(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
				}
				
				// set the graphic defaults
				sourceGraphic.setDefaults(source.getType());
				targetGraphic.setDefaults(target.getType());
				
				// add graphics to the entries
				source.addGraphics(sourceGraphic);
				target.addGraphics(targetGraphic);
				
				// check if all went well
				if(rea != null)
					p.addReaction(rea);
				if(rel != null){
					rel.setEntry1(source.getId());
					rel.setEntry2(target.getId());
					p.addRelation(rel);
				}
				
				// add the entries to the pathway
				p.addEntry(source);
				p.addEntry(target);
				
			}
			
		}
		
		return p;
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
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String filename = "hsa00010.sif";
//		String filename = "hsa04010.sif";
		SIF2KGML sif2kgml = new SIF2KGML();
		SIFPathway sifpathway = SIF2KGML.readSIF(filename);
		// debug
//		for(SIFRelation sif : sifpathway.getRelations())
//			System.out.println(sif.getSource() + " " + sif.getInteractionType() + " " + sif.getTarget());
		Pathway p = sif2kgml.translate(sifpathway);
		sif2kgml.saveToFile(p, "SIFPathwayTest.xml");
	}
}
