package de.zbit.sbvc.io;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Class for storing a SIF-Pathway
 *
 * @author 	Manuel Ruff
 * @date 	2012-08-08
 * @version $Rev: 135$
 * @since	$Rev: 135$
 *
 */

public class SIFPathway {
	
	private String name;
	private ArrayList<SIFRelation> relations = new ArrayList<SIFRelation>();
	
	public SIFPathway(String name){
		this.name = name;
	}

	public ArrayList<SIFRelation> getRelations() {
		return relations;
	}

	public void setRelations(ArrayList<SIFRelation> relations) {
		this.relations = relations;
	}

	public String getName() {
		return name;
	}

}
