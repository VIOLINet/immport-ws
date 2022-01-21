package org.reactome.immport.ws.model;

import java.util.ArrayList;
import java.util.List;

import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonInclude.Include;

/**
 * This class is used to model Reactome pathways for the pathway hierarchy.
 * @author wug
 *
 */
@JsonInclude(Include.NON_NULL)
public class ReactomePathway {
	private String stId;
	private String name;
	private String species;
	private String type;
	private Boolean diagram;
	private String topPathway;
	// List used for order. However, the order is not there in the JSON output from
	// Reactome's content API.
	private List<ReactomePathway> children;
	
	public ReactomePathway() {
	}

	public String getSpecies() {
		return species;
	}

	public void setSpecies(String species) {
		this.species = species;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public Boolean getDiagram() {
		return diagram;
	}

	public void setDiagram(Boolean diagram) {
		this.diagram = diagram;
	}

	public String getStId() {
		return stId;
	}

	public void setStId(String stId) {
		this.stId = stId;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getTopPathway() {
		return topPathway;
	}

	public void setTopPathway(String topPathway) {
		this.topPathway = topPathway;
	}

	public List<ReactomePathway> getChildren() {
		return children;
	}

	public void setChildren(List<ReactomePathway> children) {
		this.children = children;
	}
	
	public void addChild(ReactomePathway pathway) {
		if (children == null)
			children = new ArrayList<>();
		children.add(pathway);
	}

}
