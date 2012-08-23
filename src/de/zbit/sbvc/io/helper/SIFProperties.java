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
 * @date 	2012-08-08
 * @version $Rev: 138$
 * @since	$Rev: 135$
 *
 */

public class SIFProperties {
	
	/**
	 * InteractionTypes of the Simple Interaction Format
	 * pp (protein-protein interaction)
	 * pd (protein->dna)
	 * pr (protein->reaction)
	 * rc (reaction->compound)
	 * cr (compound->reaction)
	 * gl (genetic lethal relationship)
	 * pm (protein-metabolite interaction)
	 * mp (metabolite-protein interaction)
	 *
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
