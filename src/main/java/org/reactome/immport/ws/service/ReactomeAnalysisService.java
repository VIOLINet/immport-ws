package org.reactome.immport.ws.service;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.HttpStatus;
import org.apache.commons.httpclient.methods.PostMethod;
import org.apache.commons.httpclient.methods.StringRequestEntity;
import org.apache.log4j.Logger;
import org.reactome.immport.ws.model.ReactomePathway;
import org.reactome.immport.ws.model.FI.FIAnnotation;
import org.reactome.immport.ws.model.FI.FIAnnotations;
import org.reactome.immport.ws.model.analysis.BiosampleAnalysis;
import org.reactome.immport.ws.model.queries.CytoscapeFI;
import org.reactome.immport.ws.model.queries.CytoscapeFiData;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;

@Component
public class ReactomeAnalysisService {
    private final Logger logger = Logger.getLogger(ReactomeAnalysisService.class);
    
    private final String HTTP_POST = "Post";
        
    @Autowired
    private ReactomeAnalysisConfig config;
    @Autowired
    private RScriptService rservice;
    // Want to catch the pathway list assuming it will not change
    private List<ReactomePathway> pathways;
    
    public ReactomeAnalysisService() {
    }
    
	/**
	 * Generate a list of pathways that are ordered based on their locations in the pathway
	 * hierarchy so that similar pathways are grouped together. This method uses a width-first-search
	 * algorithm to get the list.
	 * @return
	 * @throws Exception
	 */
    @SuppressWarnings("unchecked")
    public synchronized List<ReactomePathway> getHierarchicalOrderedPathways() {
    	if (pathways != null)
    		return pathways;
    	try {
    		URL url = new URL(config.getReactomePathwayHierarchyURL());
    		ObjectMapper mapper = new ObjectMapper();
    		List<ReactomePathway> topPathways = mapper.readValue(url,
    				new TypeReference<List<ReactomePathway>>() {
    		});
    		List<ReactomePathway> pathways = new ArrayList<>();
    		// Some pathways may be listed under multiple topics. For this application,
    		// we need to list them once only by choosing whatever they occur first.
    		// Use this set to keep tracking pathways based on names.
    		Set<String> listedPathways = new HashSet<>();
    		for (ReactomePathway topPathway : topPathways) {
    			traversePathway(topPathway, topPathway, pathways, listedPathways);
    		}
    		// Get rid of some values that don't need at the front-end
    		pathways.forEach(p -> {
    			p.setChildren(null);
    			p.setSpecies(null);
    			p.setType(null);
    			p.setDiagram(null);
    		});
    		this.pathways = pathways;
    		return pathways;
    	}
    	catch(Exception e) {
    		logger.error(e.getMessage(), e);
    	}
    	return Collections.EMPTY_LIST; // Just return an empty list.
    }
	
	private void traversePathway(ReactomePathway topPathway,
	                             ReactomePathway currentPathway,
	                             List<ReactomePathway> list,
	                             Set<String> listedPathways) {
		if (!currentPathway.getType().equals("Pathway") &&
			!currentPathway.getType().equals("TopLevelPathway"))
			return; // Check pathways only
		if (listedPathways.contains(currentPathway.getName()))
			return; // Added to the list already
		currentPathway.setTopPathway(topPathway.getName());
		list.add(currentPathway);
		listedPathways.add(currentPathway.getName());
		if (currentPathway.getChildren() == null)
			return; // No need to go down.
		for (ReactomePathway child : currentPathway.getChildren())
			traversePathway(topPathway, child, list, listedPathways);
	}
    
    /**
     * Query Reactome FI service 
     * @param genes
     * @return
     */
    public List<CytoscapeFI> constructFINetwork(Set<String> genes) {
    	try {
    		String fiText = callHttp(config.getReactomeFIServiceURL() + "/network/queryFIs", HTTP_POST, String.join(",", genes));

    		//Add all genes to list of fisToReturn after converting to cytoscape objects
    		List<CytoscapeFI> fisToReturn = getGeneCytoscapeFIs(genes);
    		//Add all edges after getting annotations
    		fisToReturn.addAll(annotateCytoscapeFIEdges(fiText));
    		
       		return fisToReturn;
    	} catch (IOException e) {
    		logger.error(e);
    		return new ArrayList<>();
    	}
    }
    
    private List<CytoscapeFI> getGeneCytoscapeFIs(Set<String> genes) {
    	List<CytoscapeFI> rtn = new ArrayList<>();
    	genes.forEach(x -> {
    		rtn.add(new CytoscapeFI("nodes", new CytoscapeFiData(x,x,null,null)));
    	});
    	
    	return rtn;
	}

	private Collection<CytoscapeFI> annotateCytoscapeFIEdges(String fiText) throws IOException {
    	ObjectMapper objectMapper = new ObjectMapper();
    	JsonNode fis = objectMapper.readTree(fiText);
    	
    	//null check in case there are no edges in the network
    	JsonNode interactions = fis.get("interaction");
    	if(interactions == null) return new ArrayList<>();
    	
    	StringBuilder annotatedFiPostString = new StringBuilder();
    	Map<Integer, CytoscapeFI> idToCyEdge = new HashMap<>();
    	
    	//loop over interactions and create two things
    	//annotatedFiPostString is "id\tprotein1\tprotein2\n" for annotation service
    	//idToCyEdge is key value from id to CytoscapeFI so annotation direction can easily be added
    	for(int i = 0; i< interactions.size(); i++) {
    		JsonNode node = interactions.get(i);
    		String[] proteins = {node.get("firstProtein").get("name").asText(), node.get("secondProtein").get("name").asText()};
    		Arrays.sort(proteins); //sort so always alphabetical
    		//append to post data for getting annotations
    		annotatedFiPostString.append(i + "\t" + proteins[0] + "\t" + proteins[1] + "\n");
    		//make Cytoscape edge to add annotation direction to later
    		idToCyEdge.put(i, new CytoscapeFI("edges", new CytoscapeFiData(proteins[0] + "&&" + proteins[1], null, proteins[0], proteins[1])));
    	}
    	
		String annotatedFiText = callHttp(config.getReactomeFIServiceURL()+"/network/annotate", HTTP_POST, annotatedFiPostString.toString());
		
		//convert annotation return to map of id to object to make adding direction to Cytoscape edge easier
		Map<Integer, FIAnnotation> annotationMap = convertAnnotatedFiTextToMap(annotatedFiText);
		
		//add direction to edges
		idToCyEdge.forEach((id, edge) -> {
			edge.getData().setDirection(annotationMap.get(id).getDirection());
		});
		
		//return just CytoscapeFI list from map values
    	return idToCyEdge.values();
	}

	private Map<Integer, FIAnnotation> convertAnnotatedFiTextToMap(String annotatedFiText) throws IOException {
		ObjectMapper mapper = new ObjectMapper();
		List<FIAnnotation> annotations = mapper.readValue(annotatedFiText, FIAnnotations.class).getFiAnnotation();
		
		return annotations.stream().collect(Collectors.toMap(FIAnnotation::getInteractionId, x -> x));
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
			logger.error(e);
			return "";
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
    	PostMethod method = new PostMethod(url);
    	HttpClient client = new HttpClient();
    	
    	method.setRequestEntity(new StringRequestEntity(query, "text/plain", "UTF-8"));
    	method.setRequestHeader("Accept", "application/json");
    	
    	int responseCode = client.executeMethod(method);
    	if(responseCode == HttpStatus.SC_OK) {
    		return method.getResponseBodyAsString();
    	}
    	else return "Reactome service responded with code: " + responseCode;
    	
    }

	public List<BiosampleAnalysis> analyzeBiosamples(String jsonText) {
//	    System.out.println(jsonText);
		
		String response = "";
 		try {			
			//make call to plumber  server for R script
 		    String url = "http://localhost:" + rservice.getPort() + "/doDiffExpAnalysis";
 		    // For local test with RStudio
// 		    String url = "http://localhost:8087/doDiffExpAnalysis";
			PostMethod method = new PostMethod(url);
			method.addParameter("selection.json", jsonText);
			HttpClient client = new HttpClient();
			client.executeMethod(method);
			
			//get body of response
			response =  method.getResponseBodyAsString();
			response = response.replaceAll("\\\\", "");
			response = response.substring(2, response.length()-2);
			return structureBiosampleAnalysis(response);

			
		} catch (IOException e) {
			logger.error(e);
			return new ArrayList<>();
		}
	}

	private List<BiosampleAnalysis> structureBiosampleAnalysis(String response) throws IOException {
		ObjectMapper mapper = new ObjectMapper();
		List<BiosampleAnalysis> analysisObjs = mapper.readValue(response, mapper.getTypeFactory().constructCollectionType(List.class, BiosampleAnalysis.class));
		return analysisObjs;
	}
}
