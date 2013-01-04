/*
 * $Id$
 * $URL$
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
package de.zbit.sbvc.io.helper;

import java.util.ArrayList;

/**
 * Class for storing a SIF-Pathway
 *
 * @author 	Manuel Ruff
 * @date 	2012-08-08
 * @version $Rev: 138$
 * @since	Revision 135
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
	public void setRelations(ArrayList<SIFRelation> relations) {
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
