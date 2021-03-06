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
	String source;
	@JsonInclude(Include.NON_NULL)
	String target;
	@JsonInclude(Include.NON_NULL)
	String clusterColor;
	@JsonInclude(Include.NON_NULL)
	String direction;
	
	public CytoscapeFiData() {}
	
	public CytoscapeFiData(String id, String name, String source, String target) {
		this.id = id;
		this.name = name;
		this.source = source;
		this.target = target;
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
	
	public String getSource() {
		return source;
	}

	public void setSource(String source) {
		this.source = source;
	}
	public String getTarget() {
		return target;
	}

	public void setTarget(String target) {
		this.target = target;
	}

	public void setClusterColor(String clusterColor) {
		this.clusterColor = clusterColor;
	}

	public String getClusterColor() {
		return clusterColor;
	}

	public String getDirection() {
		return direction;
	}

	public void setDirection(String direction) {
		this.direction = direction;
	}
	
}
