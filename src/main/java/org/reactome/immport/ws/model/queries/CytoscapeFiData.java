package org.reactome.immport.ws.model.queries;

import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonInclude.Include;

/**
 * 
 * @author brunsont
 *
 */
public class CytoscapeFiData {
	
	String id;
	@JsonInclude(Include.NON_NULL)
	String name;
	@JsonInclude(Include.NON_NULL)
	String taret;
	@JsonInclude(Include.NON_NULL)
	String source;
	
	public CytoscapeFiData(String id, String name, String taret, String source) {
		this.id = id;
		this.name = name;
		this.taret = taret;
		this.source = source;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getTaret() {
		return taret;
	}

	public void setTaret(String taret) {
		this.taret = taret;
	}

	public String getSource() {
		return source;
	}

	public void setSource(String source) {
		this.source = source;
	}
}
