package de.zbit.sbvc.io;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
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
import de.zbit.kegg.parser.pathway.ext.EntryExtended;

/**
 * Class for converting {@link Sbgn} into a {@link Pathway}
 *
 * @author 	Manuel Ruff
 * @date 	2012-06-03
 * @version $Rev: 136$
 * @since	$Rev: 99$
 *
 */
public class SBGN2KGML {

	public static final Logger log = Logger.getLogger(SBGN2KGML.class.getName());					// logger for errors, warnings, etc.
	
	public SBGN2KGMLHelper helper = new SBGN2KGMLHelper();											// helper class to get the proper EntryType and ArcType
	
	
	HashMap<String, Entry> entryLookup = new HashMap<String, Entry>();								// mapping from the glyph id onto the entries
	HashMap<String, String> clonemarkerLookup = new HashMap<String, String>();						// mapping from the clonemarker name onto the keggid
	HashMap<String, Glyph> processGlyphLookup = new HashMap<String, Glyph>();						// mapping from the process glyph id onto the glyph itself
	
	HashMap<String, ArrayList<String>> portSubstrates = new HashMap<String, ArrayList<String>>();
	HashMap<String, ArrayList<String>> portProducts = new HashMap<String, ArrayList<String>>();
	HashMap<String, ArrayList<String>> reactionModifiers = new HashMap<String, ArrayList<String>>();	// mapping from the process glyph id onto the reactionmodifiers
	
	ArrayList<Arc> arcsWithoutPorts = new ArrayList<Arc>();											// mapping for arcs without ports (relations)
	
	private int entryID = 0;																		// entries id
	private int reactionID = 0;
	private int skipCount = 0;																		// counts the clonemarker occurences to guarantee ongoing names
	
	private boolean considerGlyphs = true;
	private boolean considerArcs = true;
	
	/**
	 * Constructor
	 */
	public SBGN2KGML(){}

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
	 * 								<li> <code>true</code> if the file is valid</li>
	 * 								<li> <code>false</code> otherwise</li>
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
		
		Iterator<Entry> entryIterator = entryLookup.values().iterator();
		
		while(entryIterator.hasNext()){
			Entry e = entryIterator.next();
			p.addEntry(e);
		}
		
		return p;
	}
	
	/**
	 * Method to convert all {@link Glyph}'s into {@link Pathway#getEntries()}
	 * @param sbgn	{@link Sbgn}
	 * @param p		{@link Pathway}
	 */
	protected void handleAllGlyphs(Sbgn sbgn, Pathway p){
		
		for (Glyph g : sbgn.getMap().getGlyph()) {
			
			createMappingForProcessGlyphsOnly(g);
			
			// no label means process glyph
			if(g.getLabel() != null){
				
				Entry e = new EntryExtended(p, ++entryID, "unknown_entry:"+(entryID-skipCount));

				// convert the GlyphType into an EntryType
				e.setType(helper.getEntryTypeFromGlyphType(GlyphType.valueOfString(g.getClazz())));
				
//				if (e instanceof EntryExtended) {
//				  // TODO: Set extended type
//				  ((EntryExtended)e).setGeneType(null);
//				}
				
				// create the graphics
				// TODO: eventually use different creators
				Graphics gr = null;
				
				if(e.getType().compareTo(EntryType.compound) == 0)
					gr = Graphics.createGraphicsForCompound(g.getLabel().getText());
				else if(e.getType().compareTo(EntryType.map) == 0)
					gr = Graphics.createGraphicsForPathwayReference(g.getLabel().getText());
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
				}
				else {
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
		
		createMappingForAllPortsAndArcs(sbgn);
			
		// create edges for the relations
		for(Arc awp : arcsWithoutPorts) {
			
			Glyph sourceGlyph = (Glyph) awp.getSource();
			Glyph targetGlyph = (Glyph) awp.getTarget();
			
			int sourceEntryID = entryLookup.get(sourceGlyph.getId()).getId();
			int targetEntryID = entryLookup.get(targetGlyph.getId()).getId();
			Relation relation = new Relation(sourceEntryID, targetEntryID, helper.getRelationTypeFromArcType(ArcType.valueOfString(awp.getClazz())));
			p.addRelation(relation);
			
		}

		// create edges for the reactions
		Iterator<Glyph> processGlyphIterator = processGlyphLookup.values().iterator();
		
		while(processGlyphIterator.hasNext()){
			
			Glyph processGlyph = processGlyphIterator.next();
			ArrayList<String> substrates = new ArrayList<String>();
			ArrayList<String> products = new ArrayList<String>();
			
			for(Port port : processGlyph.getPort()){
				if(portSubstrates.containsKey(port.getId()))
					substrates.add(port.getId());
				if(portProducts.containsKey(port.getId()))
					products.add(port.getId());
			}
			
			Reaction r = null;
			
			// TODO
			if(substrates.isEmpty() || products.isEmpty()){
				r = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.reversible);
				
				if(!products.isEmpty()){
					r.addSubstrate(new ReactionComponent(entryLookup.get(portProducts.get(products.get(0)).get(0))));
					r.addProduct(new ReactionComponent(entryLookup.get(portProducts.get(products.get(1)).get(0))));
				}
				if(!substrates.isEmpty()){
					r.addSubstrate(new ReactionComponent(entryLookup.get(portSubstrates.get(substrates.get(0)).get(0))));
					r.addProduct(new ReactionComponent(entryLookup.get(portSubstrates.get(substrates.get(1)).get(0))));
				}
			}
			else{
				r = new Reaction(p, "unknown_reaction:"+(++reactionID), ReactionType.irreversible);
			
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
			
			for(ReactionComponent rc : r.getSubstrates()){
				Entry e = rc.getCorrespondingEntry();
				if(e.isSetReaction())
					e.setReaction(e.getReactionString() + " " + r.getName());
				else
					e.setReaction(r.getName());
			}
			
			for(ReactionComponent rc : r.getProducts()){
				Entry e = rc.getCorrespondingEntry();
				if(e.isSetReaction())
					e.setReaction(e.getReactionString() + " " + r.getName());
				else
					e.setReaction(r.getName());
			}
			
			if(reactionModifiers.containsKey(processGlyph.getId())){
				ArrayList<String> rmIDs = reactionModifiers.get(processGlyph.getId());
				for(String rm : rmIDs){
					entryLookup.get(rm).setReaction(r.getName());
				}
			}
			
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
		String filename = "glycolysis.sbgn";
//		String filename = "Test.sbgn";
		SBGN2KGML sbgn2kgml = new SBGN2KGML();
		Sbgn sbgn = sbgn2kgml.read(filename);
		
		System.out.println(sbgn2kgml.validateSBGN(filename));
		
		Pathway p = sbgn2kgml.translate(sbgn);
		sbgn2kgml.saveToFile(p, "PathwayTest.xml");
	}

}
