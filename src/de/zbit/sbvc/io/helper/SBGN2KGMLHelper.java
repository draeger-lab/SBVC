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
