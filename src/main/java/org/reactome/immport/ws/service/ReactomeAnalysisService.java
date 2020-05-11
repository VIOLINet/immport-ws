package org.reactome.immport.ws.service;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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
    		ObjectMapper mapper = new ObjectMapper();
    		rtn = mapper.writeValueAsString(convertToCyJson(fiText, genes));
    	} catch (IOException e) {
    		return "";
    	}
        return rtn;
    }
    
    private List<CytoscapeFI> convertToCyJson(String fiText, Set<String> genes) throws IOException {
    	ObjectMapper objectMapper = new ObjectMapper();
    	JsonNode fis = objectMapper.readTree(fiText);
    	
    	if(!fis.get("interaction").isArray()) return new ArrayList<>();
    	    	
    	List<CytoscapeFI> rtn = new ArrayList<>();
    	genes.forEach(x -> {
    		rtn.add(new CytoscapeFI("nodes", new CytoscapeFiData(x,x,null,null)));
    	});
    	
    	List<CytoscapeFI> edges = new ArrayList<>();
    	for(final JsonNode node : fis.get("interaction")) {
    		edges.add(new CytoscapeFI("edges", new CytoscapeFiData("e"+edges.size(), null, node.get("firstProtein").get("name").asText(),node.get("secondProtein").get("name").asText())));
    	}
    	
    	rtn.addAll(edges);
		return rtn;
	}
    
    public List<CytoscapeFI> constructClusteredFINetwork(List<CytoscapeFI> network) {
    	Set<String> genes = new HashSet<>();
		List<String> fisToCluster = new ArrayList<>();
		for(CytoscapeFI fi : network) {
			if(fi.getGroup().contains("edges"))
				fisToCluster.add(fi.getData().getId() + "\t" + fi.getData().getSource() + "\t" + fi.getData().getTarget());
			else
				genes.add(fi.getData().getName());
		}
		try {
			String response = callHttp(config.getReactomeFIServiceURL()+"/network/cluster", HTTP_POST, String.join("\n", fisToCluster) + "\n");
			return addClustering(network, response);
		} catch(IOException e) {
			return new ArrayList<>();
		}
	}
    
	private List<CytoscapeFI> addClustering(List<CytoscapeFI> network, String clusterResponse) throws IOException {
		Map<String, String> clusterMap = makeClusterMap(clusterResponse);
		
		for(CytoscapeFI obj : network) {
			if(!obj.getGroup().contains("nodes")) continue;
			obj.getData().setClusterColor(clusterMap.get(obj.getData().getName()));
		}
		
		return network;
	}

	/**
	 * Makes map of cluster call to FI service from gene to a hex color based on cluster number
	 * @param clusterResponse
	 * @return
	 * @throws IOException 
	 */
	private Map<String, String> makeClusterMap(String clusterResponse) throws IOException {
		ObjectMapper objectMapper = new ObjectMapper();
    	JsonNode fis = objectMapper.readTree(clusterResponse);
    	if(!fis.get("geneClusterPairs").isArray()) return null;
    	
		Map<String, String> geneToColorMap = new HashMap<>();
		
		for(JsonNode node : fis.get("geneClusterPairs")) {
			geneToColorMap.put(node.get("geneId").toString().replaceAll("^\"|\"$", ""), node.get("cluster").toString().replaceAll("^\"|\"$", ""));
		}
		
		return geneToColorMap;
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
    	method.setRequestEntity(new StringRequestEntity(query, "text/plain", null));
    	method.setRequestHeader("Accept", "application/json");
    	    	
    	client = new HttpClient();
    	
    	int responseCode = client.executeMethod(method);
    	if(responseCode == HttpStatus.SC_OK) {
    		return method.getResponseBodyAsString();
    	}
    	else return "Response code: " + responseCode;
    	
    }
}
