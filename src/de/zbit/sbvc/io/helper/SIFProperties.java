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
 * @since	Revision 135
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
	 * 
	 * <li>COMPONENT_OF</li>
	 * <li>CO_CONTROL</li>
	 * <li>INTERACTS_WITH</li>
	 * <li>IN_SAME_COMPONENT</li>
	 * <li>METABOLIC_CATALYSIS</li>
	 * <li>REACTS_WITH</li>
	 * <li>SEQUENTIAL_CATALYSIS</li>
	 * <li>STATE_CHANGE</li>
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
		mp,
		COMPONENT_OF,
		CO_CONTROL,
		INTERACTS_WITH,
		IN_SAME_COMPONENT,
		METABOLIC_CATALYSIS,
		REACTS_WITH,
		SEQUENTIAL_CATALYSIS,
		STATE_CHANGE;
		
		/**
		 * Check if the {@link InteractionType} belongs to the two letter code
		 * @return {@link Boolean}
		 */
		public boolean isShortInteractionType(){
			return (   this.equals(pp) 
					|| this.equals(pd) 
					|| this.equals(pr) 
					|| this.equals(rc) 
					|| this.equals(cr) 
					|| this.equals(gl) 
					|| this.equals(pm) 
					|| this.equals(mp));
		}
		
		/**
		 * Check if the {@link InteractionType} belongs to the long version code
		 * @return {@link Boolean}
		 */
		public boolean isLongInteractionType(){
			return (   this.equals(COMPONENT_OF) 
					|| this.equals(CO_CONTROL) 
					|| this.equals(INTERACTS_WITH) 
					|| this.equals(IN_SAME_COMPONENT) 
					|| this.equals(METABOLIC_CATALYSIS) 
					|| this.equals(REACTS_WITH) 
					|| this.equals(SEQUENTIAL_CATALYSIS) 
					|| this.equals(STATE_CHANGE));
		}
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
				case pd:
				case pr:
				case pm:
					return EntryType.gene;
				case rc:
					return EntryType.reaction;
				case cr:
	      case mp:
					return EntryType.compound;
				default:
					return EntryType.other;
			}
		} else {
			switch(iType) {
				case pp:
				case mp:
					return EntryType.gene;
				case pr:
				case cr:
					return EntryType.reaction;
				case rc:
				case pm:
					return EntryType.compound;
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
//			case COMPONENT_OF:
//				return null;
//			case CO_CONTROL:
//				return null;
//			case INTERACTS_WITH:
//				return null;
//			case IN_SAME_COMPONENT:
//				return null;
//			case METABOLIC_CATALYSIS:
//				return null;
//			case REACTS_WITH:
//				return null;
//			case SEQUENTIAL_CATALYSIS:
//				return null;
//			case STATE_CHANGE:
//				return null;
			default:
				return RelationType.other;
		}
		
	}
	
}
