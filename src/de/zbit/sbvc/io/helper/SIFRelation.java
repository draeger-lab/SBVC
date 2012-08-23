package de.zbit.sbvc.io.helper;

import de.zbit.sbvc.io.helper.SIFProperties.InteractionType;

/**
 * Class for storing a SIF Relation
 *
 * @author 	Manuel Ruff
 * @date 	2012-08-08
 * @version $Rev: 138$
 * @since	$Rev: 135$
 *
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
	 * @return interactionType {@link InteractionType#}
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
