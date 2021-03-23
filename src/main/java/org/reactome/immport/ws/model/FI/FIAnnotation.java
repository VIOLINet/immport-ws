package org.reactome.immport.ws.model.FI;

/**
 * 
 * @author brunsont
 *
 */
public class FIAnnotation {
	private String annotation;
	private String reverseAnnotation;
	private String direction;
	private int interactionId;
	private double score;
	
	public FIAnnotation() {
		super();
	}

	public FIAnnotation(String annotation, String direction, int interactionId, double score) {
		super();
		this.annotation = annotation;
		this.direction = direction;
		this.interactionId = interactionId;
		this.score = score;
	}

	public FIAnnotation(String annotation, String reverseAnnotation, String direction, int interactionId,
			double score) {
		super();
		this.annotation = annotation;
		this.reverseAnnotation = reverseAnnotation;
		this.direction = direction;
		this.interactionId = interactionId;
		this.score = score;
	}

	public String getAnnotation() {
		return annotation;
	}

	public void setAnnotation(String annotation) {
		this.annotation = annotation;
	}

	public String getReverseAnnotation() {
		return reverseAnnotation;
	}

	public void setReverseAnnotation(String reverseAnnotation) {
		this.reverseAnnotation = reverseAnnotation;
	}

	public String getDirection() {
		return direction;
	}

	public void setDirection(String direction) {
		this.direction = direction;
	}

	public int getInteractionId() {
		return interactionId;
	}

	public void setInteractionId(int interactionId) {
		this.interactionId = interactionId;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}
	
}
