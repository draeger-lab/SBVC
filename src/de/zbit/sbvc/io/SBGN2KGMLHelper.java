package de.zbit.sbvc.io;

import java.util.logging.Logger;

import de.zbit.graph.io.def.SBGNProperties.ArcType;
import de.zbit.graph.io.def.SBGNProperties.GlyphType;
import de.zbit.kegg.parser.pathway.EntryType;
import de.zbit.kegg.parser.pathway.Relation;
import de.zbit.kegg.parser.pathway.RelationType;
import de.zbit.kegg.parser.pathway.SubType;
import de.zbit.kegg.parser.pathway.ext.EntryTypeExtended;

/**
 * Helper Class for {@link SBGN2KGML}
 * 
 * @author 	Manuel Ruff
 * @date 	2012-06-21
 * @version $Rev: 136$
 * @since	$Rev: 114$
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
	 * Method to get the corresponding {@link RelationType} from the {@link ArcType}
	 * @param atype	{@link ArcType}
	 * @return		{@link RelationType}
	 */
	public RelationType getRelationTypeFromArcType(ArcType atype){
		// TODO: RelationType:
//	  RelationType.
		switch(atype){
			case absolute_inhibition:
			  return SubType.INHIBITION;
      case absolute_stimulation:
        
      case assignment:
      case catalysis:
        // TODO: das ist dann ein enzym, welches eine reaktion katalysisert.
        // d.h. hier sollte entry.setReaction() gesetzt werden.
        // Wenn es eine Relation ist, sollte SubType.ACTIVATION verwendet werden.
      case consumption:
      case equivalence_arc:
      case inhibition:
        return SubType.INHIBITION;
      case interaction:
        return SubType.MISSING_INTERACTION;
      case logic_arc:
      case modulation:
      case necessary_stimulation:
      case negative_influence:
        return SubType.INHIBITION;
      case positive_influence:
        return SubType.ACTIVATION;
      case production:
      case stimulation:
        return SubType.ACTIVATION;
      case unknown_influence:
      default:
				return RelationType.other;
		}
	}
	
	//#######################
	//# Getters and Setters #
	//#######################
	
}
