package de.zbit.sbvc.io;

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
 * @version $Rev: 135$
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

	public RelationType getRelationTypeFromInteractionType(InteractionType iType){
		
		switch(iType) {
			case pp:
				return RelationType.PPrel;
			case pd:
				return RelationType.GErel;
			case pr:
				return null;
			case rc:
				return null;
			case cr:
				return null;
			case gl:
				return null;
			case pm:
				return null;
			case mp:
				return null;
			default:
				return RelationType.other;
		}
		
	}
}
