package org.reactome.immport.ws.service;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.HttpStatus;
import org.apache.commons.httpclient.methods.PostMethod;
import org.apache.commons.httpclient.methods.StringRequestEntity;
import org.apache.log4j.Logger;
import org.reactome.immport.ws.model.queries.CytoscapeFI;
import org.reactome.immport.ws.model.queries.CytoscapeFiData;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;


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
    public String constructFINetwork(Set<String> genes) {
    	String rtn = "";
    	try {
    		String fiText = callHttp(config.getReactomeFIServiceURL() + "/network/queryFIs", HTTP_POST, String.join(",", genes));
    		rtn = convertToCyJson(fiText);
    	} catch (IOException e) {
    		return "";
    	}
        System.out.println("Reactome FI Url: " + config.getReactomeFIServiceURL());
        return rtn;
    }
    
    private String convertToCyJson(String fiText) throws IOException {
    	ObjectMapper objectMapper = new ObjectMapper();
    	JsonNode fis = objectMapper.readTree(fiText);
    	List<CytoscapeFI> rtn = new ArrayList<>();
    	
    	if(!fis.get("interaction").isArray()) return "";
    	
    	int counter = 0;
    	for(final JsonNode node : fis.get("interaction")) {
    		rtn.add(new CytoscapeFI("nodes", new CytoscapeFiData(node.get("firstProtein").get("name").asText(),node.get("firstProtein").get("name").asText(),null,null)));
    		rtn.add(new CytoscapeFI("nodes", new CytoscapeFiData(node.get("secondProtein").get("name").asText(),node.get("secondProtein").get("name").asText(),null,null)));
    		rtn.add(new CytoscapeFI("edges", new CytoscapeFiData(counter+"", null, node.get("firstProtein").get("name").asText(),node.get("secondProtein").get("name").asText())));
    		counter++;
    	}
    	
		return objectMapper.writeValueAsString(rtn);
	}

	/**
     * Query analysis service for a set of genes
     * @param genes
     * @return
     */
    public String doPathwayEnrichmentAnalysis(Set<String> genes) {
    	String analysisText = "";
    	try {
    		analysisText = callHttp(config.getReactomeAnalysisURL(), HTTP_POST, String.join(",", genes));
		} catch (IOException e) {
			e.printStackTrace();
		}
        System.out.println("Reactome URL: " + config.getReactomeAnalysisURL());
        return analysisText;
    }

	/**
     * Perform HTTP call to passed in URL. Currently only supports post requests
     * with passed in String query. Returns a string of the response body,
     * or an empty string on request error.
     * @param url
     * @param type
     * @param query
     * @return
     * @throws IOException
     */
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
