package org.reactome.immport.ws.model.queries;

/**
 * 
 * @author brunsont
 *
 */
public class CytoscapeFI {
	String group;
	CytoscapeFiData data;
	
	public CytoscapeFI() {}
	
	public CytoscapeFI(String group, CytoscapeFiData data) {
		this.group = group;
		this.data = data;
	}

	public String getGroup() {
		return group;
	}

	public void setGroup(String group) {
		this.group = group;
	}

	public CytoscapeFiData getData() {
		return data;
	}

	public void setData(CytoscapeFiData data) {
		this.data = data;
	}
}
