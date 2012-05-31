package de.zbit.sbvc.io;

import java.io.File;
import java.io.IOException;

import javax.xml.bind.JAXBException;

import org.sbgn.SbgnUtil;
import org.sbgn.bindings.Arc;
import org.sbgn.bindings.Glyph;
import org.sbgn.bindings.Sbgn;
import org.xml.sax.SAXException;

import de.zbit.kegg.parser.pathway.Entry;
import de.zbit.kegg.parser.pathway.Graphics;
import de.zbit.kegg.parser.pathway.Pathway;

public class SBGN2KGML {

	  public static void main(String args[]) {
		  String filename = "glycolysis.sbgn";
		  SBGN2KGML sbgn2kgml = new SBGN2KGML();
		  Sbgn sbgn = sbgn2kgml.read(filename);
	  }
	  
	  public Sbgn read(String filename){
		  File f = new File(filename);
		  Sbgn sbgn = null;
		  try {
			sbgn = SbgnUtil.readFromFile(f);
		} catch (JAXBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return sbgn;
	  }
	  
	  public boolean sbgnValidation(String filename){
		  File f = new File(filename);
		  boolean isValid = false;
		  try {
			isValid = SbgnUtil.isValid(f);
		} catch (JAXBException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (SAXException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return isValid;
	  }
	  
	  protected Pathway translate(Sbgn sbgn){
		  
		  // create a new pathway
		  // TODO: how to get the parameter of this object?
		  Pathway p = new Pathway("name", "org", 0);
		  
		  // for every glyph
		  for(Glyph g : sbgn.getMap().getGlyph()){
			  
			  // create the graphics
			  Graphics gr = new Graphics(e);
			  gr.setX((int) g.getBbox().getX());
			  gr.setY((int) g.getBbox().getY());
			  gr.setHeight((int) g.getBbox().getH());
			  gr.setWidth((int) g.getBbox().getW());
			  
			  // create an entry
			  // TODO: determine the entrytype
			  Entry e = new Entry(p, 0, g.getLabel().toString(), null, gr);
			  
			  // add the entry to the pathway
			  p.addEntry(e);
		  }
		  
		  // for every arc
		  for(Arc arc : sbgn.getMap().getArc()){
			  
		  }
		  
		  return null;
	  }
	
}
