package de.zbit.sbvc.io;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumSet;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.xml.bind.JAXBException;

import org.sbgn.SbgnUtil;
import org.sbgn.bindings.Arc;
import org.sbgn.bindings.Glyph;
import org.sbgn.bindings.Port;
import org.sbgn.bindings.Sbgn;
import org.xml.sax.SAXException;

import de.zbit.graph.io.def.SBGNProperties.GlyphType;
import de.zbit.kegg.KGMLWriter;
import de.zbit.kegg.parser.pathway.Entry;
import de.zbit.kegg.parser.pathway.EntryType;
import de.zbit.kegg.parser.pathway.Graphics;
import de.zbit.kegg.parser.pathway.Pathway;
import de.zbit.kegg.parser.pathway.Relation;
import de.zbit.kegg.parser.pathway.RelationType;
import de.zbit.kegg.parser.pathway.ext.EntryExtended;

/**
 * Class for converting {@link Sbgn} into a {@link Pathway}
 *
 * @author 	Manuel Ruff
 * @date 	2012-06-03
 * @version $Rev: 99$
 * @since	$Rev: 99$
 *
 */
public class SBGN2KGML {

	public static final Logger log = Logger.getLogger(SBGN2KGML.class.getName());	// logger for errors, warnings, etc.
	
	HashMap<String, GlyphType> glyphTypeLookup = new HashMap<String, GlyphType>();	// mapping from the glyphs clazz onto the GlyphType
	HashMap<String, Glyph> glyphLookup = new HashMap<String, Glyph>();				// mapping from the glyph ids onto the glyph itself
	HashMap<String, Entry> entryLookup = new HashMap<String, Entry>();				// mapping from the glyph ids onto the entries
	HashMap<String, Entry> clonemarkerLookup = new HashMap<String, Entry>();		// mapping from the clonemarker name onto the entry
	HashMap<String, Glyph> processGlyphLookup = new HashMap<String, Glyph>();		// mapping from the process glyph id's onto the glyph itself
	HashMap<String, ArrayList<String>> portLinkageCollectorLookup = new HashMap<String, ArrayList<String>>();	// mapping collection for all in or outgoing glyphs of a port
	private int id = 0;	// entries id
	
	private boolean considerGlyphs = true;
	private boolean considerArcs = true;
	
	/**
	 * Constructor
	 */
	public SBGN2KGML(){
		initialize();
	}
	
	/**
	 * Method for initializing {@link SBGN2KGML}. It simply creates a mapping for the {@link Glyph#getClazz()} onto {@link GlyphType}
	 */
	protected void initialize(){
		// for all GlyphTypes
		for(GlyphType g : EnumSet.allOf(GlyphType.class)){
			// generate a lookup from String to GlyphType
			glyphTypeLookup.put(g.toString(), g);
		}
	}
	
	/**
	 * Method to get the corresponding {@link EntryType} from the {@link GlyphType}
	 * @param gtype	{@link GlyphType}
	 * @return		{@link EntryType}
	 */
	// TODO: this mapping is kinda wired, because it is surjective
	protected EntryType getEntryTypeFromGlyphType(GlyphType gtype){
		switch(gtype){
			case simple_chemical:
				return EntryType.compound;
			case macromolecule:
				return EntryType.enzyme;
			case macromolecule_multimer:
				return EntryType.genes;
			default:
				return EntryType.other;
		}
	}

	/**
	 * Method for reading a {@link Sbgn}-file
	 * @param filename		the name of the file as {@link String}
	 * @return {@link Sbgn}	returns a {@link Sbgn}-object
	 * @exception JAXBException
	 */
	public Sbgn read(String filename) {
		// initialize a file
		File f = new File(filename);
		// initilize the sbgn object
		Sbgn sbgn = null;
		try {
			// read in the file into the sbgn object
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
	 * 								<li> <code>true</code> if the file is valid</li>
	 * 								<li> <code>false</code> otherwise</li>
	 * 							</ul>
	 * @exception JAXBException
	 * @exception SAXException
	 * @exception IOException
	 */
	public boolean validateSBGN(String filename) {
		File f = new File(filename);
		boolean isValid = false;
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
	 * Method for translating a {@link Sbgn}-file into a KGML-file
	 * @param {@link Sbgn}		(see {@link SBGN2KGML#read(String)})
	 * @return {@link Pathway}	(see {@link KGMLWriter})
	 */
	protected Pathway translate(Sbgn sbgn) {

		// create a new pathway
		// TODO: the number of the pathway should be a five-digit number so what?
		Pathway p = new Pathway("unknown", "unknown", 10000, "unknown title");

		// translate the glyphs
		if(considerGlyphs)
			handleAllGlyphs(sbgn, p);
		
		// translate the arcs
		if(considerArcs)
			handleAllArcs(sbgn, p);
		
		return p;
	}
	
	/**
	 * Method to convert all {@link Glyph}'s into {@link Pathway#getEntries()}
	 * @param sbgn	{@link Sbgn}
	 * @param p		{@link Pathway}
	 */
	protected void handleAllGlyphs(Sbgn sbgn, Pathway p){
		
		// for every glyph
		for (Glyph g : sbgn.getMap().getGlyph()) {
			
			// map all the ports and subglyphs
			createMappingForGlyph(g);
			
			// no label means process glyph
			if(g.getLabel() != null){
				
				// create a new entry
				Entry e = new EntryExtended(p, ++id, "unknown:"+id);

				// convert the glyphtype into entrytype
				e.setType(getEntryTypeFromGlyphType(GlyphType.valueOfString(g.getClazz())));
				
//				if (e instanceof EntryExtended) {
//				  // TODO: Set extended type
//				  ((EntryExtended)e).setGeneType(null);
//				}
				
				// create the graphics
				// TODO: eventually use different creators
//				Graphics.createGraphicsForCompound(name);
//				Graphics.createGraphicsForPathwayReference(name);
//				else
				Graphics gr = new Graphics(g.getLabel().getText());
				
				gr.setDefaults(e.getType());
				gr.setX((int) g.getBbox().getX());
				gr.setY((int) g.getBbox().getY());
				gr.setHeight((int) g.getBbox().getH());
				gr.setWidth((int) g.getBbox().getW());
	
				// add the graphics to the entry
				e.addGraphics(gr);
				
				// map the glyph onto the entries
				entryLookup.put(g.getId(), e);
				
				// same entries have the same name
				if(g.getClone() != null)
					if(clonemarkerLookup.containsKey(g.getLabel().getText()))
						e.setName(clonemarkerLookup.get(g.getLabel().getText()).getName());
				else
					clonemarkerLookup.put(g.getLabel().getText(), e);
				
				// add the entry to the pathway
				p.addEntry(e);
			}
		}
		
	}
	
	/**
	 * Method to convert all the {@link Arc}'s into {@link Pathway#getReactions()} and {@link Pathway#getRelations()}
	 * @param sbgn	{@link Sbgn}
	 * @param p		{@link Pathway}
	 */
	protected void handleAllArcs(Sbgn sbgn, Pathway p){
		
		// TODO determine the RelationType
		// TODO NOW it connects all together - thats wrong
		// TODO add enzyme support
		
		createMappingForAllPorts(sbgn);
		
		// for every arc
		for (Arc arc : sbgn.getMap().getArc()) {
			
			Port firstPort = null;
			Port secondPort = null;
			Entry source = null;
			Entry target = null;
			ArrayList<Relation> relations = new ArrayList<Relation>();
			
			if(arc.getTarget().getClass().equals(Port.class)) {
				firstPort = (Port) arc.getTarget();
			}
			
			if(firstPort == null && arc.getSource().getClass().equals(Glyph.class) && arc.getTarget().getClass().equals(Glyph.class)) {
				source = entryLookup.get(((Glyph) arc.getSource()).getId());
				target = entryLookup.get(((Glyph) arc.getTarget()).getId());
				
				if(source != null && target != null) {
					Relation relation = new Relation(source.getId(), target.getId(), RelationType.other);
					relations.add(relation);
				}
			}
			
			if(firstPort != null && arc.getSource().getClass().equals(Glyph.class)) {
				Glyph processGlyph = processGlyphLookup.get(firstPort.getId());
				
				for (Port currentPort : processGlyph.getPort()) {
					if(currentPort.getId() != firstPort.getId())
						secondPort = currentPort;
				}
					
				ArrayList<String> collectedGlyphs = portLinkageCollectorLookup.get(secondPort.getId());
						
				for (String s : collectedGlyphs) {
					Relation relation = new Relation(entryLookup.get(((Glyph) arc.getSource()).getId()).getId(), entryLookup.get(s).getId(), RelationType.other);
					relations.add(relation);
				}
			}
			
			for(Relation relation : relations)
				p.addRelation(relation);

		}
	}
	
	/**
	 * Method to create a mapping for all {@link Port}'s and their linked {@link Glyph}'s from a {@link Sbgn}
	 * @param sbgn	{@link Sbgn}
	 */
	protected void createMappingForAllPorts(Sbgn sbgn) {
		
		for (Arc arc : sbgn.getMap().getArc()) {
			
			ArrayList<String> collectedGlyphs = new ArrayList<String>();
			Port port = null;
			
			if(arc.getSource().getClass().equals(Port.class)) {
				port = (Port) arc.getSource();
				if(portLinkageCollectorLookup.containsKey(port.getId()))
						collectedGlyphs = portLinkageCollectorLookup.get(port.getId());
				collectedGlyphs.add(((Glyph) arc.getTarget()).getId());
			}
			
			if(arc.getTarget().getClass().equals(Port.class)) {
				port = (Port) arc.getTarget();
				if(portLinkageCollectorLookup.containsKey(port.getId()))
						collectedGlyphs = portLinkageCollectorLookup.get(port.getId());
				collectedGlyphs.add(((Glyph) arc.getSource()).getId());
			}
			
			if(port != null)
				portLinkageCollectorLookup.put(port.getId(), collectedGlyphs);
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
	 * Method to create a mapping for a {@link Glyph} which will be saved in {@link SBGN2KGML#glyphLookup} and {@link SBGN2KGML#portLookup}
	 * @param g	{@link Glyph}
	 */
	protected void createMappingForGlyph(Glyph g){
		// put the original glyph into the glyph lookup table
		if(g.getLabel() != null)
			glyphLookup.put(g.getId(), g);
		else {
			processGlyphLookup.put(g.getId(), g);

			for (Port port : g.getPort())
				processGlyphLookup.put(port.getId(), g);
		}
	}
	
	
	public static void main(String args[]) {
		String filename = "glycolysis.sbgn";
		SBGN2KGML sbgn2kgml = new SBGN2KGML();
		sbgn2kgml.initialize();
		Sbgn sbgn = sbgn2kgml.read(filename);
		
		System.out.println(sbgn2kgml.validateSBGN(filename));
		
		Pathway p = sbgn2kgml.translate(sbgn);
		sbgn2kgml.saveToFile(p, "PathwayTest.xml");
	}

}
