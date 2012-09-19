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
 * Copyright (C) 2012 by the University of Tuebingen, Germany.
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

import java.util.logging.Logger;

import de.zbit.graph.io.def.SBGNProperties.GlyphType;
import de.zbit.kegg.parser.pathway.EntryType;
import de.zbit.kegg.parser.pathway.Pathway;
import de.zbit.kegg.parser.pathway.RelationType;

/**
 * Helper-Class for converting Simple Interaction Files into a {@link Pathway}
 *
 * @author 	Manuel Ruff
 * @author Clemens Wrzodek
 * @date 	2012-08-08
 * @version $Rev: 138$
 * @since	$Rev: 135$
 *
 */

public class SIFProperties {
	
	/**
	 * InteractionTypes of the Simple Interaction Format (SIF):<ul>
	 * <li>pp (protein-protein interaction)</li>
	 * <li>pd (protein->dna)</li>
	 * <li>pr (protein->reaction)</li>
	 * <li>rc (reaction->compound)</li>
	 * <li>cr (compound->reaction)</li>
	 * <li>gl (genetic lethal relationship)</li>
	 * <li>pm (protein-metabolite interaction)</li>
	 * <li>mp (metabolite-protein interaction)</li>
	 * </ul>
	 */
	public static enum InteractionType
	{
		pp,
		pd,
		pr,
		rc,
		cr,
		gl,
		pm,
		mp;
	}
	
	/**
	 * Method to get the corresponding {@link EntryType} from the {@link InteractionType} for the first or second component
	 * @param iType {@link InteractionType}
	 * @param firstComponent {@link Boolean}
	 * @return		{@link EntryType}
	 */
	public static EntryType getEntryTypeFromInteractionType(InteractionType iType, boolean firstComponent){
		
		if(firstComponent){
			switch(iType) {
			case pp:
				return EntryType.gene;
			case pd:
				return EntryType.gene;
			case pr:
				return EntryType.gene;
			case rc:
				return EntryType.reaction;
			case cr:
				return EntryType.compound;
			case gl:
				return EntryType.other;
			case pm:
				return EntryType.gene;
			case mp:
				return EntryType.other;
			default:
				return EntryType.other;
			}
		} else {
			switch(iType) {
			case pp:
				return EntryType.gene;
			case pd:
				return EntryType.other;
			case pr:
				return EntryType.reaction;
			case rc:
				return EntryType.compound;
			case cr:
				return EntryType.reaction;
			case gl:
				return EntryType.other;
			case pm:
				return EntryType.other;
			case mp:
				return EntryType.gene;
			default:
				return EntryType.other;
			}
		}
		
	}

	/**
	 * Method to get the corresponding {@link RelationType} from the {@link InteractionType}
	 * @param gtype	{@link InteractionType}
	 * @return		{@link RelationType}
	 */
	public static RelationType getRelationTypeFromInteractionType(InteractionType iType){
		
		switch(iType) {
			case pp:
				return RelationType.PPrel;
			case pd:
				return RelationType.GErel;
//			case gl:
//				return null;
//			case pm:
//				return null;
//			case mp:
//				return null;
			default:
				return RelationType.other;
		}
		
	}
}
