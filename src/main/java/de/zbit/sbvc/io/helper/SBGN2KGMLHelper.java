/*
 * $Id: SBGN2KGMLHelper.java 180 2014-01-09 22:21:09Z draeger $
 * $URL: https://rarepos.cs.uni-tuebingen.de/svn-path/SBVC/trunk/src/de/zbit/sbvc/io/helper/SBGN2KGMLHelper.java $
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

import java.util.logging.Logger;

import de.zbit.graph.io.def.SBGNProperties.ArcType;
import de.zbit.graph.io.def.SBGNProperties.GlyphType;
import de.zbit.kegg.parser.pathway.EntryType;
import de.zbit.kegg.parser.pathway.SubType;
import de.zbit.sbvc.io.SBGN2KGML;

/**
 * Helper Class for {@link SBGN2KGML}
 * 
 * @author 	Manuel Ruff
 * @date 	2012-06-21
 * @version $Rev: 138$
 * @since	Revision 114
 *
 */
public class SBGN2KGMLHelper {
	
	public static final Logger log = Logger.getLogger(SBGN2KGMLHelper.class.getName());			// logger for errors, warnings, etc.
	
	/**
	 * Method to get the corresponding {@link EntryType} from the {@link GlyphType}
	 * @param gtype	{@link GlyphType}
	 * @return		{@link EntryType}
	 */
	public EntryType getEntryTypeFromGlyphType(GlyphType gtype){
		
		switch(gtype){
			case interaction:
				return EntryType.reaction;
			case macromolecule:
				return EntryType.enzyme;
			case macromolecule_multimer:
				return EntryType.genes;
			case nucleic_acid_feature:
				return EntryType.gene;
			case simple_chemical:
				return EntryType.compound;
			case submap:
				return EntryType.map;
			default:
				return EntryType.other;
		}
		
	}
	
	/**
	 * Method to get the corresponding {@link SubType} from the {@link ArcType}
	 * @param atype	{@link ArcType}
	 * @return		{@link SubType}
	 */
	public SubType getSubTypeFromArcType(ArcType atype){
		
		// uncommenting the unused cases results in wired results
		
		switch(atype){
			case absolute_inhibition:
				return new SubType(SubType.INHIBITION);
			case absolute_stimulation:
				return new SubType(SubType.ACTIVATION);
//			case assignment:
			case catalysis:
				return new SubType(SubType.ACTIVATION);
//			case consumption:
//			case equivalence_arc:
			case inhibition:
				return new SubType(SubType.INHIBITION);
			case interaction:
				return new SubType(SubType.MISSING_INTERACTION);
//			case logic_arc:
//			case modulation:
//			case necessary_stimulation:
			case negative_influence:
				return new SubType(SubType.INHIBITION);
			case positive_influence:
				return new SubType(SubType.ACTIVATION);
//			case production:
			case stimulation:
				return new SubType(SubType.ACTIVATION);
//			case unknown_influence:
			default:
				return null;
				
		}
	}
	
}
