package org.reactome.immport.ws.web;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.reactome.immport.ws.model.PublicRepository;
import org.reactome.immport.ws.model.Study;
import org.reactome.immport.ws.model.queries.BioSampleCollectionTime;
import org.reactome.immport.ws.model.queries.VOToGSM;
import org.reactome.immport.ws.model.requests.GSMForVOs;
import org.reactome.immport.ws.service.ImmportDAO;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.web.bind.annotation.CrossOrigin;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.ModelAttribute;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RestController;

@RestController
@CrossOrigin
public class ImmportWSController {

    @Autowired
    private ImmportDAO studyDAO;
    
    public ImmportWSController() {
    }
    
    @GetMapping("study/{accession}")
    public Study getStudy(@PathVariable("accession") String accession) {
        return studyDAO.loadStudy(accession);
    }
    
    @GetMapping("study/vaccine/{vo_term}")
    public List<Study> queryStudiesForVO(@PathVariable("vo_term") String voId) {
        return studyDAO.queryStudiesForVO(voId);
    }
    
    /**
     * voIds should be delimited by "," in the request body text.
     * @param voIds
     * @return
     */
    @PostMapping("study/vaccine")
    public List<Study> queryStudiesForVOs(@RequestBody String voIds) {
        Collection<String> c = Arrays.asList(voIds.split(","));
        return studyDAO.queryStudiesForVO(c);
    }
    
    @GetMapping("repository/{name}")
    public PublicRepository getRepository(@PathVariable("name") String name) {
        return studyDAO.loadPublicRepository(name);
    }
    
    /**
     * Structure:
     * Line 1: comma separated list of voIds
     * Line 2: comma separated list of genders
     * @param voId
     * @return
     */
    @PostMapping("expSample/vaccine")
    @CrossOrigin
    public List<VOToGSM> queryGSMDataForVOs(@RequestBody GSMForVOs gsmForVOs){
    	if(gsmForVOs == null)
    		return new ArrayList<>();
    	
//    	Collection<String> voIds = Arrays.asList(lines[0].split(","));
//    	Collection<String> genderList = Arrays.asList(lines[1].split(","));
    	return studyDAO.queryGSMDataForVO(gsmForVOs);
    }
    
    @GetMapping("collectionTimes/vaccine/{voId}")
    public List<BioSampleCollectionTime> queryTimesCollectedForVO(@PathVariable("voId") String voId){
    	return studyDAO.queryStudyTimeCollectedForVO(voId);
    }
    
    /**
     * voIds should be delimited by "," in the request body text.
     * @param text
     * @return
     */
    @PostMapping("collectionTimes/vaccine")
    public List<BioSampleCollectionTime> queryTimesCollectedForVOs(@RequestBody String text){
    	Collection<String> voIds = Arrays.asList(text.split(","));
    	if(voIds == null || voIds.size() == 0)
    		return new ArrayList<>();
    	return studyDAO.queryStudyTimeCollectedForVO(voIds);
    }
    
}
