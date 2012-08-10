package de.zbit.sbvc.io;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import de.zbit.io.csv.CSVReader;
import de.zbit.kegg.KGMLWriter;
import de.zbit.kegg.parser.pathway.Entry;
import de.zbit.kegg.parser.pathway.EntryType;
import de.zbit.kegg.parser.pathway.Graphics;
import de.zbit.kegg.parser.pathway.Pathway;
import de.zbit.kegg.parser.pathway.Reaction;
import de.zbit.kegg.parser.pathway.ReactionType;
import de.zbit.kegg.parser.pathway.Relation;
import de.zbit.kegg.parser.pathway.ext.EntryExtended;

/**
 * Class for converting Simple Interaction Files into a {@link Pathway}
 *
 * @author 	Manuel Ruff
 * @date 	2012-08-08
 * @version $Rev: 135$
 * @since	$Rev: 135$
 *
 */

public class SIF2KGML {
	
	public static final Logger log = Logger.getLogger(SIF2KGML.class.getName());					// logger for errors, warnings, etc.
	
	private int entryID = 0;																		// entries id
	private int reactionID = 0;
	
	/**
	 * Method for reading in a Simple Interaction File
	 * @param filename	{@link String}
	 * @return {@Link ArrayList<SifRelation>)
	 * @throws IOException 
	 */
	public static SIFPathway readSIF(String filename) throws IOException {
		
		if(!filename.endsWith(".sif"))
			log.log(Level.SEVERE, String.format("The file %s doesn't end with .sif", filename));
		
		Scanner scanner = null;
		// take the only the filename and set it as Pathwayname
		SIFPathway sif = new SIFPathway(filename.replaceAll("\\.sif", "").replaceAll(".*/", "").replaceAll(".*\\\\", ""));
		ArrayList<SIFRelation> relations = new ArrayList<SIFRelation>();

		try {
			scanner = new Scanner(new FileReader(filename)).useDelimiter("\n");
			
			String line;
			String[] splitted;
			int counter = 0;
			
			while(scanner.hasNext()) {
				line = scanner.nextLine();
				counter++;
				splitted = line.split("\t");
				// check if the file fulfills the patterns
				if(splitted.length != 3)
					splitted = line.split(" ");
				
				if(splitted.length != 3)
					log.log(Level.WARNING, String.format("The line %s in the file %s doesn't have the form <nodeA> <relationship> <nodeB>", new Object[]{counter, filename}));
				
				relations.add(new SIFRelation(splitted[0], SIFProperties.InteractionType.valueOf(splitted[1]), splitted[2]));
				
				}
		} catch (FileNotFoundException e) {
			
		}
		return sif;
	}
	
	/**
	 * Method for translating a {@link SIFPathway} a {@link Pathway}
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
			source.setType(EntryType.compound);
			target.setType(EntryType.compound);
			
			//create Graphics
			Graphics sourceGraphic = null;
			Graphics targetGraphic = null;
			
			// since all are compounds createGraphicsForCompound
			sourceGraphic = Graphics.createGraphicsForCompound(source.getName());
			targetGraphic = Graphics.createGraphicsForCompound(target.getName());
			
			// set the graphic defaults
			sourceGraphic.setDefaults(source.getType());
			targetGraphic.setDefaults(target.getType());

			// add graphics to the entries
			source.addGraphics(sourceGraphic);
			target.addGraphics(targetGraphic);
			
			Reaction rea = null;
			Relation rel = null;
			
			switch(relation.getInteractionType()){
				case pr:
					rea = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.irreversible);
					source.setReaction(rea.getName());
					target.setReaction(rea.getName());
					break;
				case rc:
					rea = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.irreversible);
					source.setReaction(rea.getName());
					target.setReaction(rea.getName());
					break;
				case cr:
					rea = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.irreversible);
					source.setReaction(rea.getName());
					target.setReaction(rea.getName());
					break;
				default:
					rel = new Relation(source.getId(), target.getId(), null);
			}
			
			if(rea != null)
				p.addReaction(rea);
			if(rel != null)
				p.addRelation(rel);
			
			p.addEntry(source);
			p.addEntry(target);
			
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
	 */
	public static void main(String[] args) {
		String filename = "Test.sif";
		SIF2KGML sif2kgml = new SIF2KGML();
		SIFPathway sifpathway = sif2kgml.readSIF(filename);
		Pathway p = sif2kgml.translate(sifpathway);
		sif2kgml.saveToFile(p, "SIFPathwayTest.xml");
	}
}
