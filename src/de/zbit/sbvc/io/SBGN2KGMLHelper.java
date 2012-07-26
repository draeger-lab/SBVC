package de.zbit.sbvc.io;

import java.util.logging.Logger;

import de.zbit.graph.io.def.SBGNProperties.ArcType;
import de.zbit.graph.io.def.SBGNProperties.GlyphType;
import de.zbit.kegg.parser.pathway.EntryType;
import de.zbit.kegg.parser.pathway.RelationType;
import de.zbit.kegg.parser.pathway.ext.EntryTypeExtended;

/**
 * Helper Class for {@link SBGN2KGML}
 * 
 * @author 	Manuel Ruff
 * @date 	2012-06-21
 * @version $Rev: 114$
 * @since	$Rev: 114$
 *
 */
public class SBGN2KGMLHelper {
	
	public static final Logger log = Logger.getLogger(SBGN2KGMLHelper.class.getName());			// logger for errors, warnings, etc.
	
	/**
	 * Constructor
	 */
	public SBGN2KGMLHelper(){}
	
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
	 * Method to get the corresponding {@link RelationType} from the {@link ArcType}
	 * @param atype	{@link ArcType}
	 * @return		{@link RelationType}
	 */
	public RelationType getRelationTypeFromArcType(ArcType atype){
		
		switch(atype){
			default:
				return RelationType.other;
		}
	}
	
	//#######################
	//# Getters and Setters #
	//#######################
	
}
