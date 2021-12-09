package org.reactome.immport.ws.test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import org.apache.commons.httpclient.HttpClient;
import org.apache.commons.httpclient.HttpMethod;
import org.apache.commons.httpclient.HttpStatus;
import org.apache.commons.httpclient.methods.GetMethod;
import org.apache.commons.httpclient.methods.PostMethod;
import org.apache.commons.httpclient.methods.RequestEntity;
import org.apache.commons.httpclient.methods.StringRequestEntity;
import org.junit.Test;
import org.reactome.immport.ws.model.requests.GSMForVOs;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;

public class WSTests {
    protected final String HOST_URL = "http://localhost:8076/reactome-immport-ws";
    protected final String HTTP_POST = "Post";
    protected final String HTTP_GET = "Get";

    
    @Test
    public void testGetStudy() throws Exception {
        String url = HOST_URL + "/study/SDY212";
        testGetMethod(url);
    }
    
    private void testGetMethod(String url) throws Exception {
        System.out.println(url);
        String rtn = callHttp(url, HTTP_GET, null);
        outputJSON(rtn);
    }
    
    @Test
    public void testRepository() throws Exception {
        String url = HOST_URL + "/repository/GEO";
        testGetMethod(url);
    }
    
    @Test
    public void testQueryStudiesForVO() throws Exception {
    	String url = HOST_URL + "/study/vaccine/VO_0000642";
    	String rtn = callHttp(url, HTTP_GET, null);
    	outputJSON(rtn);
    }
    
    @Test
    public void testQueryStudiesForVOs() throws Exception {
    	String url = HOST_URL + "/study/vaccine";
    	String query = "VO_0000044,VO_0000642,VO_0004809,VO_0000045,VO_0000046,VO_0000047";
    	String rtn = callHttp(url, HTTP_POST, query);
    	outputJSON(rtn);
    }
    
    @Test
    public void testQueryTimesForVO() throws Exception {
    	String url  = HOST_URL + "/collectionTimes/vaccine/VO_0000642";
    	String rtn = callHttp(url,HTTP_GET, null);
    	outputJSON(rtn);
    }
    
    @Test
    public void testQueryTimesForVOs() throws Exception {
    	String url  = HOST_URL + "/collectionTimes/vaccine";
    	String query = "VO_0004809,VO_0000047,VO_0000642,VO_0000044,VO_0000045"+"\n"+"male,female";
    	String rtn = callHttp(url, HTTP_POST, query);
    	outputJSON(rtn);
    }
    
    @Test
    public void testAnalyzeBiosamples() throws Exception {
    	String url = HOST_URL + "/analysis/geneExpression";
    	String queryJson = "{\"GSMids\":[\"GSM1934006\",\"GSM1934008\",\"GSM1934009\",\"GSM1934011\",\"GSM1934012\",\"GSM1934014\",\"GSM1934015\",\"GSM1934016\",\"GSM1934017\",\"GSM1934019\",\"GSM1934020\",\"GSM1934022\",\"GSM1934023\",\"GSM1934025\",\"GSM1934026\",\"GSM1934028\",\"GSM1934029\",\"GSM1934031\",\"GSM1934032\",\"GSM1934034\",\"GSM1934035\",\"GSM1934037\",\"GSM1934038\",\"GSM1934040\",\"GSM1934041\",\"GSM1934043\",\"GSM1934044\",\"GSM1934046\",\"GSM1934047\",\"GSM1934049\",\"GSM1934050\",\"GSM1934052\",\"GSM1934053\",\"GSM1934055\",\"GSM1934056\",\"GSM1934058\",\"GSM1934059\",\"GSM1934061\",\"GSM1934062\",\"GSM1934064\",\"GSM1934065\",\"GSM1934067\",\"GSM1934068\",\"GSM1934070\",\"GSM1934071\",\"GSM1934073\",\"GSM1934074\",\"GSM1934076\",\"GSM1934077\",\"GSM1934079\",\"GSM1934080\",\"GSM1934082\",\"GSM1934083\",\"GSM1934085\",\"GSM1934086\",\"GSM1934088\",\"GSM733942\",\"GSM733950\",\"GSM733965\",\"GSM733944\",\"GSM733947\",\"GSM733953\",\"GSM733956\",\"GSM733959\",\"GSM733962\",\"GSM733967\",\"GSM733970\",\"GSM733973\",\"GSM733976\",\"GSM733979\",\"GSM733982\",\"GSM733985\",\"GSM733987\",\"GSM733990\",\"GSM733993\",\"GSM733996\",\"GSM733999\",\"GSM734001\",\"GSM734004\",\"GSM734007\",\"GSM734010\",\"GSM734013\",\"GSM734016\",\"GSM734019\",\"GSM733943\",\"GSM733952\",\"GSM733946\",\"GSM733949\",\"GSM733955\",\"GSM733958\",\"GSM733961\",\"GSM733964\",\"GSM733969\",\"GSM733972\",\"GSM733975\",\"GSM733978\",\"GSM733981\",\"GSM733984\",\"GSM733989\",\"GSM733992\",\"GSM733995\",\"GSM733998\",\"GSM734000\",\"GSM734003\",\"GSM734006\",\"GSM734009\",\"GSM734012\",\"GSM734015\",\"GSM734018\",\"GSM734021\",\"GSM734051\",\"GSM734043\",\"GSM734068\",\"GSM734075\",\"GSM734091\",\"GSM734055\",\"GSM734047\",\"GSM734071\",\"GSM734079\",\"GSM734095\",\"GSM734053\",\"GSM734045\",\"GSM734070\",\"GSM734077\",\"GSM734093\",\"GSM734057\",\"GSM734049\",\"GSM734073\",\"GSM734081\",\"GSM734097\",\"GSM734052\",\"GSM734044\",\"GSM734069\",\"GSM734076\",\"GSM734092\",\"GSM734056\",\"GSM734048\",\"GSM734072\",\"GSM734080\",\"GSM734096\",\"GSM734050\",\"GSM734042\",\"GSM734066\",\"GSM734074\",\"GSM734090\",\"GSM734054\",\"GSM734046\",\"GSM734067\",\"GSM734078\",\"GSM734094\",\"GSM733816\",\"GSM733819\",\"GSM733822\",\"GSM733825\",\"GSM733828\",\"GSM733831\",\"GSM733834\",\"GSM733837\",\"GSM733840\",\"GSM733818\",\"GSM733821\",\"GSM733824\",\"GSM733827\",\"GSM733830\",\"GSM733833\",\"GSM733836\",\"GSM733839\",\"GSM733842\"],\"studyCohort\":{\"genders\":[\"Female\",\"Male\"],\"ages\":[\"0.0 - 46.0\",\"0.0 - 47.0\",\"21.0 - 47.0\"],\"races\":[\"Asian\",\"Black or African American\",\"White\"]},\"resultSetName\":\"Untitled_1\",\"modelTime\":false,\"analysisGroups\":{\"group1\":[0],\"group2\":[7]},\"studyVariables\":[],\"platformCorrection\":true,\"usePairedData\":true,\"variableGenes\":false}";
    	String rtn = callHttp(url,HTTP_POST, queryJson);
    	outputJSON(rtn);
    }
    
    @Test
    public void testPathwayAnalysisForGenes() throws Exception{
    	String url = HOST_URL + "/analysis/pathways";
    	List<String> genes = new ArrayList<>(Arrays.asList("LRRN3", "SERPINE2", "ACTA2", "CCR7", "SYT11", "ADRB2"));
    	ObjectMapper mapper = new ObjectMapper();
    	String json = mapper.writeValueAsString(genes);
    	String rtn = callHttp(url, HTTP_POST, json);
    	outputJSON(rtn);
    }
    
    @Test
    public void testGetBiosampleMetadata() throws Exception {
    	String url = HOST_URL + "/metadata/biosamples";
    	String rtn = callHttp(url, HTTP_GET, null);
    	System.out.println(rtn);
    }
    
    private String outputJSON(String json) throws JsonProcessingException, IOException {
        ObjectMapper mapper = new ObjectMapper();
        Object obj = mapper.readValue(json, Object.class);
        String rtn = mapper.writerWithDefaultPrettyPrinter().writeValueAsString(obj);
        System.out.println(rtn);
        return rtn;
    }

    protected String callHttp(String url,
                              String type,
                              String query) throws IOException {
        HttpMethod method = null;
        HttpClient client = null;
        if (type.equals(HTTP_POST)) {
            method = new PostMethod(url);
            client = initializeHTTPClient((PostMethod) method, query);
        } else {
            method = new GetMethod(url); // Default
            client = new HttpClient();
        }
        method.setRequestHeader("Accept", "application/json");
        method.setRequestHeader("content-type", "application/json");
        int responseCode = client.executeMethod(method);
        if (responseCode == HttpStatus.SC_OK) {
            InputStream is = method.getResponseBodyAsStream();
            return readMethodReturn(is);
        } else {
            System.err.println("Error from server: " + method.getResponseBodyAsString());
            System.out.println("Response code: " + responseCode);
            throw new IllegalStateException(method.getResponseBodyAsString());
        }
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

    private HttpClient initializeHTTPClient(PostMethod post, String query) throws UnsupportedEncodingException {
        RequestEntity entity = new StringRequestEntity(query, "text/plain", "UTF-8");
        //        RequestEntity entity = new StringRequestEntity(query, "application/XML", "UTF-8");
        post.setRequestEntity(entity);
        //        post.setRequestHeader("Accept", "application/JSON, application/XML, text/plain");
        post.setRequestHeader("Accept", "application/json");
        HttpClient client = new HttpClient();
        return client;
    }

}
