package org.reactome.immport.ws.service;

import org.springframework.beans.factory.annotation.Autowired;

public class ReactomeAnalysisConfig {
    
    private String reactomeAnalysisURL;
    private String reactomeFIServiceURL;
    
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

}
