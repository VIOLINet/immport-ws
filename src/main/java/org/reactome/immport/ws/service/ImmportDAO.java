package org.reactome.immport.ws.service;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.reactome.immport.ws.model.Experiment;
import org.reactome.immport.ws.model.PublicRepository;
import org.reactome.immport.ws.model.Study;
import org.reactome.immport.ws.model.queries.BioSampleCollectionTime;
import org.reactome.immport.ws.model.queries.VOToGSM;
import org.reactome.immport.ws.model.requests.GSMForVOs;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Repository;
import org.springframework.transaction.annotation.Transactional;

@Repository
public class ImmportDAO {
    @Autowired
    private SessionFactory sessionFactory;
    
    private static boolean comparisonData = false;
    
    public ImmportDAO() {
    }
    
    @Transactional(readOnly = true)
    public Study loadStudy(String accession) {
        Study study = load(accession, Study.class);
        return study;
    }
    
    @Transactional(readOnly = true)
    public Experiment loadExperiment(String accession) {
        return load(accession, Experiment.class);
    }
    
    @Transactional(readOnly = true)
    public List<Study> queryStudiesForVO(String voId) {
        return queryStudiesForVO(Collections.singleton(voId));
    }
    
    @Transactional(readOnly = true)
    public List<Study> queryStudiesForVO(Collection<String> voIds) {
        String queryText = "SELECT DISTINCT s FROM Study s " +
                           "INNER JOIN s.arms a " + 
                           "INNER JOIN a.exposures e " + 
                           "WHERE e.exposureMaterialId in :voIds";
        Session session = sessionFactory.getCurrentSession();
        List<Study> studies = session.createQuery(queryText, Study.class)
                                     .setParameter("voIds", voIds)
                                     .getResultList();
        return studies;
    }

    @Transactional(readOnly = true)
    public List<VOToGSM> queryGSMDataForVO(GSMForVOs gSMForVOs){
    	String queryText = "SELECT ge.repositoryAccession, ie.exposureMaterialId, sub.gender, " +
    					   "bs.studyTimeCollected, bs.studyTimeCollectedUnit, bs.study.accession, ge.repository.name FROM Subject sub " + 
    					   "INNER JOIN sub.immuneExposures ie " +
    					   "INNER JOIN sub.biosamples bs " + 
    					   "INNER JOIN bs.expSamples es " +
    					   "INNER JOIN es.geneExpressions ge " +
    					   "WHERE ie.exposureMaterialId in :voIds " + 
    					   "AND sub.gender in :sjGender";
		Session session = sessionFactory.getCurrentSession();
		List<Object[]> repositoryAccessions = session.createQuery(queryText, Object[].class)
						                           .setParameter("voIds", gSMForVOs.getVoIds())
						                           .setParameter("sjGender", gSMForVOs.getGenderList())
						                           .list();
		List<VOToGSM> result = new ArrayList<>();
		for(Object[] obj : repositoryAccessions) {
			String date = obj[3].toString() + " " + obj[4];
			if(obj[0].toString().equals("Not available yet") || !gSMForVOs.getTimes().contains(date)) continue;
			result.add(new VOToGSM(obj[0].toString(), obj[1].toString(), obj[2].toString(), obj[3].toString() + " " + obj[4].toString(), obj[5].toString(), obj[6].toString()));
		}
		
    	return result;
    }
    
    public Set<String> getTestGeneSymbols(){
    	Set<String> geneIDs = new HashSet<>();
    	
    	if(comparisonData) return getTestComparisonData();
    	
    	InputStream is = getClass().getClassLoader().getResourceAsStream("test.tsv");
    	Scanner scanner = new Scanner(is);
    	String line = scanner.nextLine();
    	while((line = scanner.nextLine()) != null) {
    		String gene = line.split("\t")[2].replaceAll("^\"|\"$", "");
    		if(gene.isEmpty() || gene.contains("///") || gene.contains("@")) continue;  
    		geneIDs.add(gene);
    		if(!scanner.hasNextLine()) break;
    		
    	}
    	scanner.close();
    	comparisonData = !comparisonData;
    	return geneIDs;
    }
    
    private Set<String> getTestComparisonData() {
		Set<String> geneIDs = new HashSet<>();
		
		InputStream is = getClass().getClassLoader().getResourceAsStream("comparison.txt");
		Scanner scanner =  new Scanner(is);
		String line = scanner.nextLine();
		while((line = scanner.nextLine()) != null) {
			geneIDs.add(line);
			if(!scanner.hasNextLine()) break;
		}
		scanner.close();
		comparisonData = !comparisonData;
		return geneIDs;
	}

	@Transactional(readOnly = true)
    public List<BioSampleCollectionTime> queryStudyTimeCollectedForVO(String voId){
    	return queryStudyTimeCollectedForVO(Collections.singleton(voId));
    }
    
    @Transactional(readOnly = true)
    public List<BioSampleCollectionTime> queryStudyTimeCollectedForVO(Collection<String> voIds){
    	String queryText = "SELECT DISTINCT bs.studyTimeCollected, bs.studyTimeCollectedUnit FROM BioSample bs " +
    					   "INNER JOIN bs.subject sub " + 
    					   "INNER JOIN sub.immuneExposures ie " +
    					   "WHERE ie.exposureMaterialId in :voIds";
    	Session session = sessionFactory.getCurrentSession();
    	List<Object[]> times = session.createQuery(queryText, Object[].class)
    								  .setParameter("voIds", voIds)
    								  .list();
    	
    	List<BioSampleCollectionTime> rtn = new ArrayList<>();
    	for(Object[] obj : times) {
    		rtn.add(new BioSampleCollectionTime(Double.parseDouble(obj[0].toString()), obj[1].toString()));
    	}
    	
    	return rtn;
    }
    
    
    @Transactional(readOnly = true) // Apparently the transaction should not be moved to the load method.
    public PublicRepository loadPublicRepository(String name) {
        return load(name, PublicRepository.class);
    }
    
    private <T> T load(String accession, Class<T> cls) {
        Session session = sessionFactory.getCurrentSession();
        return session.get(cls, accession);
    }
    
}
