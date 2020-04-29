package org.reactome.immport.ws.test;

import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;

import java.io.IOException;

import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.mock.web.MockHttpServletResponse;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;
import org.springframework.test.context.web.WebAppConfiguration;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.test.web.servlet.MvcResult;
import org.springframework.test.web.servlet.ResultActions;
import org.springframework.test.web.servlet.request.MockMvcRequestBuilders;
import org.springframework.test.web.servlet.setup.MockMvcBuilders;
import org.springframework.web.context.WebApplicationContext;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;

/**
 * This class is used to test RESTful APIs.
 * @author wug
 *
 */
@RunWith(SpringJUnit4ClassRunner.class)
@WebAppConfiguration
@ContextConfiguration("test-servlet-context.xml")
public class SpringMVCTests {
    
    @Autowired
    private WebApplicationContext wac;
    private MockMvc mockMVC;
    
    public SpringMVCTests() {
    }
    
    @Before
    public void setup() {
        this.mockMVC = MockMvcBuilders.webAppContextSetup(this.wac).build();
    }
    
    @Test
    public void getStudy() throws Exception {
        String url = "/study/SDY212";
        testGet(url);
    }
    
    @Test
    public void queryStudyForVO() throws Exception {
        String url = "/study/vaccine/VO_0000044";
        testGet(url);
    }
    
    @Test
    public void queryStudyForVOviaPost() throws Exception {
        String url = "/study/vaccine";
        String voIds = "VO_0000044,VO_0000642,VO_0004809,VO_0000045,VO_0000046,VO_0000047";
        testPost(url, voIds);
    }
    
    @Test
    public void getAnalysisResults() throws Exception {
        String url = "/analysis/pathways";
        testGet(url);
    }
    
    @Test
    public void getFINetwork() throws Exception {
        String url = "/analysis/fi_network";
        testGet(url);
    }

    private void testGet(String url) throws Exception {
        ResultActions actions = this.mockMVC.perform(get(url));
        outputResult(actions);
    }
    
    private void testPost(String url, String body) throws Exception {
        ResultActions actions = mockMVC.perform(MockMvcRequestBuilders.post(url).content(body));
        outputResult(actions);
    }

    private void outputResult(ResultActions actions) throws Exception {
        MvcResult result = actions.andReturn();
        MockHttpServletResponse obj = result.getResponse();
        outputJSON(obj.getContentAsString());
    }
    
    @Test
    public void getRepository() throws Exception {
        String url = "/repository/GEO";
        testGet(url);
    }
    
    private String outputJSON(String json) throws JsonProcessingException, IOException {
//        System.out.println(json);
        ObjectMapper mapper = new ObjectMapper();
        Object obj = mapper.readValue(json, Object.class);
        String rtn = mapper.writerWithDefaultPrettyPrinter().writeValueAsString(obj);
        System.out.println(rtn);
        return rtn;
    }

}
