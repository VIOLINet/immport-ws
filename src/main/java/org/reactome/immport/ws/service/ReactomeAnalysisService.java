package org.reactome.immport.ws.service;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Collection;

import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.HttpStatus;
import org.apache.commons.httpclient.methods.PostMethod;
import org.apache.commons.httpclient.methods.StringRequestEntity;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class ReactomeAnalysisService {
    private final Logger logger = Logger.getLogger(ReactomeAnalysisService.class);
    
    private final String HTTP_POST = "Post";
    
    @Autowired
    private ReactomeAnalysisConfig config;
    
    public ReactomeAnalysisService() {
    }
    
    /**
     * Query Reactome FI service 
     * @param genes
     * @return
     */
    public String constructFINetwork(Collection<String> genes) {
    	String fiText = "";
    	try {
    		fiText = callHttp(config.getReactomeFIServiceURL() + "/network/buildNetwork", HTTP_POST, genes.toString());
    	} catch (IOException e) {
    		e.printStackTrace();
    	}
        System.out.println("Reactome FI Url: " + config.getReactomeFIServiceURL());
        return fiText;
    }
    
    /**
     * Query analysis service for a set of genes
     * @param genes
     * @return
     */
    public String doPathwayEnrichmentAnalysis(Collection<String> genes) {
    	String analysisText = "";
    	try {
    		analysisText = callHttp(config.getReactomeAnalysisURL(), HTTP_POST, genes.toString());
		} catch (IOException e) {
			e.printStackTrace();
		}
        System.out.println("Reactome URL: " + config.getReactomeAnalysisURL());
        return analysisText;
    }
    
    private String callHttp(String url, String type, String query) throws IOException {
    	PostMethod method = null;
    	HttpClient client = null;
    	
    	method = new PostMethod(url);
    	method.setRequestEntity(new StringRequestEntity(query, "text/plain", "UTF-8"));
    	method.setRequestHeader("Accept", "application/json");
    	
    	client = new HttpClient();
    	
    	int responseCode = client.executeMethod(method);
    	if(responseCode == HttpStatus.SC_OK) {
    		return method.getResponseBodyAsString();
    	}
    	else return "";
    	
    }
}
