package org.reactome.immport.ws.service;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.Collection;

import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.HttpMethod;
import org.apache.commons.httpclient.HttpStatus;
import org.apache.commons.httpclient.methods.PostMethod;
import org.apache.commons.httpclient.methods.RequestEntity;
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
    	System.out.println(fiText);
        System.out.println("Reactome FI Url: " + config.getReactomeFIServiceURL());
        return "{\"FIs\": [\"EGF\\tEGFR\"]}";
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
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	System.out.println(analysisText);
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
    
    private String callHttp(String url, String type, String query) throws IOException {
    	PostMethod method = null;
    	HttpClient client = null;
    	
    	method = new PostMethod(url);
    	method.setRequestEntity(new StringRequestEntity(query, "text/plain", "UTF-8"));
    	method.setRequestHeader("Accept", "application/json");
    	
    	client = new HttpClient();
    	
    	int responseCode = client.executeMethod(method);
    	if(responseCode == HttpStatus.SC_OK) {
    		InputStream is = method.getResponseBodyAsStream();
    		return readMethodReturn(is);
    	}
    	else return "";
    	
    }
    protected String readMethodReturn(InputStream is) throws IOException {
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader reader = new BufferedReader(isr);
        StringBuilder builder = new StringBuilder();
        String line = null;
        while ((line = reader.readLine()) != null)
            builder.append(line).append("\n");
        reader.close();
        isr.close();
        is.close();
        // Remove the last new line
        String rtn = builder.toString();
        // Just in case an empty string is returned
        if (rtn.length() == 0)
            return rtn;
        return rtn.substring(0, rtn.length() - 1);
    }
}
