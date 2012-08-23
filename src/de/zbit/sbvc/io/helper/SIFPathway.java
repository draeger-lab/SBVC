package de.zbit.sbvc.io.helper;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Class for storing a SIF-Pathway
 *
 * @author 	Manuel Ruff
 * @date 	2012-08-08
 * @version $Rev: 138$
 * @since	$Rev: 135$
 *
 */

public class SIFPathway {
	
	private String name;														// name of the Pathway (mainly this is the raw filename)
	private ArrayList<SIFRelation> relations = new ArrayList<SIFRelation>();	// storage of all relations of a Pathway
	
	
	/**
	 * Constructor
	 * @param name	{@link String}
	 */
	public SIFPathway(String name){
		this.name = name;
	}

	//////////////////////////
	// Getters and Setters //
	////////////////////////
	
	/**
	 * Method for getting the {@link SIFPathway#relations} element
	 * @return relations {@link ArrayList<SIFRelation>}
	 */
	public ArrayList<SIFRelation> getRelations() {
		return relations;
	}

	/**
	 * Method for setting the {@link SIFPathway#relations} element
	 * @param relations {@link ArrayList<SIFRelation>}
	 */
	public void addRelations(ArrayList<SIFRelation> relations) {
		this.relations = relations;
	}

	/**
	 * Method for getting the {@link SIFPathway#name} element
	 * @return name {@link String}
	 */
	public String getName() {
		return name;
	}

}
