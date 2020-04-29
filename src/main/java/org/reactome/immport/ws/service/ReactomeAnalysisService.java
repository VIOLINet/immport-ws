package org.reactome.immport.ws.service;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class ReactomeAnalysisService {
    private final Logger logger = Logger.getLogger(ReactomeAnalysisService.class);
    
    @Autowired
    private ReactomeAnalysisConfig config;
    
    public ReactomeAnalysisService() {
    }
    
    public String constructFINetwork() {
        System.out.println("Reactome FI Url: " + config.getReactomeFIServiceURL());
        return "{\"FIs\": [\"EGF\\tEGFR\"]}";
    }
    
    public String doPathwayEnrichmentAnalysis() {
        System.out.println("Reactome URL: " + config.getReactomeAnalysisURL());
        return "{\n" + 
                "  \"summary\": {\n" + 
                "    \"token\": \"MjAyMDA0MjkwMzUyMDBfNDU3OQ%3D%3D\",\n" + 
                "    \"projection\": false,\n" + 
                "    \"interactors\": false,\n" + 
                "    \"type\": \"OVERREPRESENTATION\",\n" + 
                "    \"sampleName\": \"\",\n" + 
                "    \"text\": true,\n" + 
                "    \"includeDisease\": true\n" + 
                "  }}";
    }

}
