package org.reactome.immport.ws.model.analysis;

import java.util.List;

public class BiosampleAnalyses {
	List<BiosampleAnalysis> biosampleAnalyses;

	public BiosampleAnalyses() {
		super();
	}
	
	public BiosampleAnalyses(List<BiosampleAnalysis> biosampleAnalyses) {
		super();
		this.biosampleAnalyses = biosampleAnalyses;
	}

	public List<BiosampleAnalysis> getBiosampleAnalyses() {
		return biosampleAnalyses;
	}

	public void setBiosampleAnalyses(List<BiosampleAnalysis> biosampleAnalyses) {
		this.biosampleAnalyses = biosampleAnalyses;
	}
}
