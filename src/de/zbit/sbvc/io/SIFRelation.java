package de.zbit.sbvc.io;

import de.zbit.sbvc.io.SIFProperties.InteractionType;

/**
 * Class for storing a SIF Relation
 *
 * @author 	Manuel Ruff
 * @date 	2012-08-08
 * @version $Rev: 135$
 * @since	$Rev: 135$
 *
 */

public class SIFRelation {
	
	private String source;
	private InteractionType interactionType;
	private String target;
	
	public SIFRelation(String source, InteractionType interactionType, String target){
		this.source = source;
		this.interactionType = interactionType;
		this.target = target;
	}

	public String getSource() {
		return source;
	}

	public InteractionType getInteractionType() {
		return interactionType;
	}

	public String getTarget() {
		return target;
	}

}
