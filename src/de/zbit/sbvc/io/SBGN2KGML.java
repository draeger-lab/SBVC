package de.zbit.sbvc.io;

import java.io.File;
import java.io.IOException;
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
	HashMap<String, Entry> handleGlyphs = new HashMap<String, Entry>();				// mapping from the clonemarker names onto the entries
	HashMap<String, Glyph> glyphLookup = new HashMap<String, Glyph>();				// mapping from the glyph id's onto the glyph itself
	HashMap<String, Port> portLookup = new HashMap<String, Port>();					// mapping from the port id's onto the ports itself
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
		Pathway p = new Pathway("path:unknown", "unknown organism", 10000);

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
				Entry e = new Entry(p, ++id, "unknown:"+id);
				// TODO: get the proper EntryType from the Glyphtype
				e.setType(getEntryTypeFromGlyphType(GlyphType.valueOfString(g.getClazz())));
				
				// create the graphics
				Graphics gr = new Graphics(g.getLabel().getText());
				gr.setX((int) g.getBbox().getX());
				gr.setY((int) g.getBbox().getY());
				gr.setHeight((int) g.getBbox().getH());
				gr.setWidth((int) g.getBbox().getW());
				gr.setDefaults(e.getType());
	
				// add the graphics to the entry
				e.addGraphics(gr);
				
				// check if the glyph is a clonemarker and handle it
				if(handleGlyphs.containsKey(g.getLabel().getText())){
					e.setName(handleGlyphs.get(g.getLabel().getText()).getName());
				} else {
					handleGlyphs.put(g.getLabel().getText(), e);
				}
				
				// add the entry to the pathway
				p.addEntry(e);
			} else {
				//TODO: do something with the process glyphs etc.S
			}
		}
		
	}
	
	/**
	 * Method to convert all the {@link Arc}'s into {@link Pathway#getReactions()} and {@link Pathway#getRelations()}
	 * @param sbgn	{@link Sbgn}
	 * @param p		{@link Pathway}
	 */
	protected void handleAllArcs(Sbgn sbgn, Pathway p){
		
		// for every arc
		for (Arc arc : sbgn.getMap().getArc()) {

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
		glyphLookup.put(g.getId(), g);

		// put all subglyphs (that should be state variables for example) into the glyph lookup table
		for(Glyph subGlyph : g.getGlyph())
			glyphLookup.put(subGlyph.getId(), subGlyph);
		
		// put all ports into the port lookup table
		for(Port subPort : g.getPort())
			portLookup.put(subPort.getId(), subPort);
	}
	
	
	public static void main(String args[]) {
//		String filename = "C:/Users/manu/Desktop/HiWi/SBGN2KGML/Examples/glycolysis.sbgn";
		SBGN2KGML sbgn2kgml = new SBGN2KGML();
		Glyph g = new Glyph();
		g.setClazz(GlyphType.simple_chemical.toString());
		EntryType test = sbgn2kgml.getEntryTypeFromGlyphType(GlyphType.valueOfString(g.getClazz()));
		System.out.println(test.name());
//		sbgn2kgml.initialize();
//		Sbgn sbgn = sbgn2kgml.read(filename);
//		if(sbgn2kgml.validateSBGN(filename))
//			System.out.println(true);
//		else
//			System.out.println(false);
//		Pathway p = sbgn2kgml.translate(sbgn);
//		sbgn2kgml.saveToFile(p, "PathwayTest.xml");
//		Sbgn sbgn = sbgn2kgml.read(filename);
//		Pathway p = sbgn2kgml.translate(sbgn);
//		
//		KGMLWriter.writeKGML(p, "test.xml", false);
		
		
//		System.out.println(SBGNProperties.GlyphType.valueOfString("simple chemical").getClass());
	}

}
