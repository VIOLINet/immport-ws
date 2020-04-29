package org.reactome.immport.ws.test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;

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
    protected final String HOST_URL = "http://localhost:8076/immportws";
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
