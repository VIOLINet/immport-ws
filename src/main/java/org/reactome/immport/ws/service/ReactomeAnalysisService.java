package org.reactome.immport.ws.service;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
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
    	String rtn = null;
    	try {
    		String fiText = callHttp(config.getReactomeFIServiceURL() + "/network/queryFIs", HTTP_POST, String.join(",", genes));
    		//return List<CytoscapeFI>. Spring will convert on its own
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
    
    public Map<String, String> constructClusteredFINetwork(List<CytoscapeFI> network) {
		List<String> fisToCluster = new ArrayList<>();
		for(CytoscapeFI fi : network) {
			if(fi.getGroup().contains("edges"))
				fisToCluster.add(fi.getData().getId() + "\t" + fi.getData().getSource() + "\t" + fi.getData().getTarget());
		}
		try {
			String response = callHttp(config.getReactomeFIServiceURL()+"/network/cluster", HTTP_POST, String.join("\n", fisToCluster) + "\n");
			return makeClusterMap(response);
		} catch(IOException e) {
			return new HashMap<>();
		}
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
    	String[] colorList = config.getModuleColors();
    	if(!fis.get("geneClusterPairs").isArray()) return null;
    	
		Map<String, String> geneToColorMap = new HashMap<>();
		
		for(JsonNode node : fis.get("geneClusterPairs")) {
			geneToColorMap.put(node.get("geneId").asText(), 
							   colorList.length > node.get("cluster").asInt() ?
							   "rgb(" + colorList[node.get("cluster").asInt()] + ")":
							   null);
		}
		
		return geneToColorMap;
	}

	/**
     * Query analysis service for a set of genes
     * @param genes
     * @return
     */
    public String doPathwayEnrichmentAnalysis(Set<String> genes) {
    	String analysisText = null;
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
    	else return "Reactome service responded with code: " + responseCode;
    	
    }

	public String analyzeBiosamples(String nums) {
		
		String response = "";
		
		try {
			String json = new String(Files.readAllBytes(Paths.get("/Users/brunsont/git/immport-ws/src/main/resources/de_analysis/selections.json")));
			response = callHttp("http://127.0.01:8087/doDiffExpAnalysis", "POST", json);
			System.out.println(response);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		return response;
	}
}
