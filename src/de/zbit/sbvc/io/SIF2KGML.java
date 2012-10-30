package de.zbit.sbvc.io;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import de.zbit.kegg.KGMLWriter;
import de.zbit.kegg.parser.pathway.Entry;
import de.zbit.kegg.parser.pathway.EntryType;
import de.zbit.kegg.parser.pathway.Graphics;
import de.zbit.kegg.parser.pathway.Pathway;
import de.zbit.kegg.parser.pathway.Reaction;
import de.zbit.kegg.parser.pathway.ReactionType;
import de.zbit.kegg.parser.pathway.Relation;
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
 * @since	$Rev: 135$
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
		
		for(SIFRelation relation : sif.getRelations()) {
			
			// create entries
			Entry source = new EntryExtended(p, ++entryID, relation.getSource());
			Entry target = new EntryExtended(p, ++entryID, relation.getTarget());
			
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
						break;
					case pd:
						sourceGraphic = Graphics.createGraphicsForProtein(source.getName());
						targetGraphic = new Graphics(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						break;
					case pr:
						sourceGraphic = Graphics.createGraphicsForProtein(source.getName());
						targetGraphic = Graphics.createGraphicsForPathwayReference(source.getName());
						rea = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.irreversible);
						source.setReaction(rea.getName());
						target.setReaction(rea.getName());
						break;
					case rc:
						sourceGraphic = Graphics.createGraphicsForPathwayReference(source.getName());
						targetGraphic = Graphics.createGraphicsForCompound(target.getName());
						rea = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.irreversible);
						source.setReaction(rea.getName());
						break;
					case cr:
						sourceGraphic = Graphics.createGraphicsForCompound(source.getName());
						targetGraphic = Graphics.createGraphicsForPathwayReference(target.getName());
						rea = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.irreversible);
						target.setReaction(rea.getName());
						break;
					case pm:
						sourceGraphic = Graphics.createGraphicsForProtein(source.getName());
						targetGraphic = Graphics.createGraphicsForCompound(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						break;
					case mp:
						sourceGraphic = Graphics.createGraphicsForCompound(source.getName());
						targetGraphic = Graphics.createGraphicsForProtein(target.getName());
						rel = new Relation(source.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
						break;
					case CO_CONTROL:
					case COMPONENT_OF:
					case IN_SAME_COMPONENT:
					case METABOLIC_CATALYSIS:
					case REACTS_WITH:
					case SEQUENTIAL_CATALYSIS:
					case STATE_CHANGE:
//						rel = new Relation(target.getId(), target.getId(), SIFProperties.getRelationTypeFromInteractionType(relation.getInteractionType()));
//						rel.addSubtype(source);
					case INTERACTS_WITH:
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
//		String filename = "Test.sif";
		String filename = "hsa00010.sif";
		SIF2KGML sif2kgml = new SIF2KGML();
		SIFPathway sifpathway = SIF2KGML.readSIF(filename);
		System.out.println(sifpathway.getName());
		System.out.println(sifpathway.getRelations().size());
		for(SIFRelation sif : sifpathway.getRelations())
			System.out.println(sif.getSource() + " " + sif.getInteractionType() + " " + sif.getTarget());
//		Pathway p = sif2kgml.translate(sifpathway);
//		sif2kgml.saveToFile(p, "SIFPathwayTest.xml");
	}
}
