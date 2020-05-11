package org.reactome.immport.ws.service;

import org.springframework.beans.factory.annotation.Autowired;

public class ReactomeAnalysisConfig {
    
    private String reactomeAnalysisURL;
    private String reactomeFIServiceURL;
    private String[] moduleColors;
    
    public ReactomeAnalysisConfig() {
        
    }

    public String getReactomeAnalysisURL() {
        return reactomeAnalysisURL;
    }

    public void setReactomeAnalysisURL(String reactomeAnalysisURL) {
        this.reactomeAnalysisURL = reactomeAnalysisURL;
    }

    public String getReactomeFIServiceURL() {
        return reactomeFIServiceURL;
    }

    public void setReactomeFIServiceURL(String reactomeFIServiceURL) {
        this.reactomeFIServiceURL = reactomeFIServiceURL;
    }

	public String[] getModuleColors() {
		return moduleColors;
	}

	public void setModuleColors(String[] moduleColors) {
		this.moduleColors = moduleColors;
	}

}
