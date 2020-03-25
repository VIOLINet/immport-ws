package org.reactome.immport.ws.service;

import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.reactome.immport.ws.model.Experiment;
import org.reactome.immport.ws.model.PublicRepository;
import org.reactome.immport.ws.model.Study;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Repository;
import org.springframework.transaction.annotation.Transactional;

@Repository
public class ImmportDAO {
    @Autowired
    private SessionFactory sessionFactory;
    
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
                                     .list();
        return studies;
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
