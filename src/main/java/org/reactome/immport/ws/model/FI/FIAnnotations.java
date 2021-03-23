package org.reactome.immport.ws.model.FI;

import java.util.List;

/**
 * 
 * @author brunsont
 *
 */
public class FIAnnotations {
	
	private List<FIAnnotation> fiAnnotation;
	
	public FIAnnotations() {
		super();
	}

	public FIAnnotations(List<FIAnnotation> fiAnnotation) {
		super();
		this.fiAnnotation = fiAnnotation;
	}

	public List<FIAnnotation> getFiAnnotation() {
		return fiAnnotation;
	}

	public void setFiAnnotation(List<FIAnnotation> fiAnnotation) {
		this.fiAnnotation = fiAnnotation;
	} 
	
}
