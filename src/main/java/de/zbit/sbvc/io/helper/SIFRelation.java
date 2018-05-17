/*
 * $Id: SIFRelation.java 180 2014-01-09 22:21:09Z draeger $
 * $URL: https://rarepos.cs.uni-tuebingen.de/svn-path/SBVC/trunk/src/de/zbit/sbvc/io/helper/SIFRelation.java $
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
package de.zbit.sbvc.io.helper;

import de.zbit.sbvc.io.helper.SIFProperties.InteractionType;

/**
 * Class for storing a SIF Relation
 *
 * @author 	Manuel Ruff
 * @author Clemens Wrzodek
 * @date 	2012-08-08
 * @version $Rev: 138$
 * @since	Revision 135
 */
public class SIFRelation {
	
	private String source;						// the source or substrate
	private InteractionType interactionType;	// interaction or relation type between the source and target
	private String target;						// the target or product
	
	public SIFRelation(String source, InteractionType interactionType, String target){
		this.source = source;
		this.interactionType = interactionType;
		this.target = target;
	}

	//////////////////////////
	// Getters and Setters //
	////////////////////////
	
	/**
	 * Method for getting the {@link SIFRelation#source} element
	 * @return source {@link String}
	 */
	public String getSource() {
		return source;
	}

	/**
	 * Method for getting the {@link SIFRelation#interactionType} element
	 * @return interactionType {@link InteractionType}
	 */
	public InteractionType getInteractionType() {
		return interactionType;
	}

	/**
	 * Method for getting the {@link SIFRelation#target} element
	 * @return target {@link String}
	 */
	public String getTarget() {
		return target;
	}

}
