/*
 * JBoss, Home of Professional Open Source
 * Copyright 2014, Red Hat, Inc. and/or its affiliates, and individual
 * contributors by the @authors tag. See the copyright.txt in the
 * distribution for a full listing of individual contributors.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.reactome.immport.ws.test;

import java.beans.BeanInfo;
import java.beans.Introspector;
import java.beans.PropertyDescriptor;
import java.lang.reflect.Method;
import java.util.List;
import java.util.Set;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.boot.Metadata;
import org.hibernate.boot.MetadataSources;
import org.hibernate.boot.registry.StandardServiceRegistry;
import org.hibernate.boot.registry.StandardServiceRegistryBuilder;
import org.junit.Test;
import org.reactome.immport.ws.model.AdverseEvent;
import org.reactome.immport.ws.model.ArmOrCohort;
import org.reactome.immport.ws.model.ArmToSubject;
import org.reactome.immport.ws.model.BioSample;
import org.reactome.immport.ws.model.Experiment;
import org.reactome.immport.ws.model.ExpSample;
import org.reactome.immport.ws.model.ImmuneExposure;
import org.reactome.immport.ws.model.Intervention;
import org.reactome.immport.ws.model.SampleGeneExpression;
import org.reactome.immport.ws.model.Study;
import org.reactome.immport.ws.model.Subject;
import org.reactome.immport.ws.model.Treatment;

public class HibernateTests {

    public HibernateTests() {
    }
    
    @Test
    public void testStudy() throws Exception {
        SessionFactory sf = createSessionFactory();
        Session session = sf.openSession();
        
        Study study = session.get(Study.class, "SDY212");
        checkObject(study);
        
        Set<Experiment> experiments = study.getExperiments();
        System.out.println("Total experiments: " + experiments.size());
        // This experiment should have gene expression
        Experiment experiment = experiments.stream()
                .filter(e -> e.getAccession().equals("EXP13387"))
                .findAny()
                .get();
        checkObject(experiment);
        
        Set<ExpSample> expSamples = experiment.getExpSamples();
        System.out.println("Total expSamples: " + expSamples.size());
        ExpSample expSample = expSamples.stream()
                .filter(s -> s.getGeneExpressions() != null)
                .findAny()
                .get();
        checkObject(expSample);
        
        Set<SampleGeneExpression> geneExpressions = expSample.getGeneExpressions();
        System.out.println("Total gene expression: " + (geneExpressions == null ? 0 : geneExpressions.size()));
        SampleGeneExpression geneExp = geneExpressions.stream().findAny().get();
        checkObject(geneExp);
        
        Set<ArmOrCohort> arms = study.getArms();
        System.out.println("Total arms: " + arms.size());
        ArmOrCohort arm = arms.stream().findAny().get();
        
        Set<ArmToSubject> armToSubjects = arm.getArmToSubjects();
        System.out.println("Total ArmToSubjects: " + armToSubjects.size());
        ArmToSubject armToSubject = armToSubjects.stream().findAny().get();
        checkObject(armToSubject);
        
        Subject subject = armToSubject.getSubject();
        checkObject(subject);
        
        Set<BioSample> biosamples = subject.getBiosamples();
        System.out.println("Total biosamples: " + biosamples.size());
        BioSample biosample = biosamples.stream().findAny().get();
        checkObject(biosample);
        
        Set<ImmuneExposure> exposures = subject.getImmuneExposures();
        System.out.println("Total immune exposures: " + exposures.size());
        ImmuneExposure exposure = exposures.stream().findAny().get();
        checkObject(exposure);
        
        session.close();
    }
    
    @Test
    public void testArmToSubject() throws Exception {
        SessionFactory sf = createSessionFactory();
        Session session = sf.openSession();
        List<ArmToSubject> armToSubjects = session.createQuery("FROM " + ArmToSubject.class.getSimpleName(), ArmToSubject.class)
                                         .list();
        System.out.println("Total ArmToSubject: " + armToSubjects.size());
        session.close();
    }
    
    /**
     * Run this method for testing the Hibernate mapping.
     * @throws Exception
     */
    @Test
    public void testLoad() throws Exception {
        SessionFactory sf = createSessionFactory();
        Session session = sf.openSession();
        // Test a study
        Study study = session.get(Study.class, "SDY212");
        checkObject(study);
        
        // Test an experiment
        Experiment experiment = session.get(Experiment.class,"EXP10074");
        checkObject(experiment);
        
        // Check ExperimentSample
        ExpSample expSample = session.get(ExpSample.class, "ES100008");
        checkObject(expSample);
        
        // Check SampleGeneExpression: Note the id should be Long
        SampleGeneExpression geneExpression = session.get(SampleGeneExpression.class, 1000L);
        checkObject(geneExpression);
        
        // Check ArmOrCohort
        ArmOrCohort arm = session.get(ArmOrCohort.class, "ARM3368");
        checkObject(arm);
        
        // Check Subject
        Subject subject = session.get(Subject.class, "SUB73456");
        checkObject(subject);
        
        // Check BioSample
        BioSample bioSample = session.get(BioSample.class, "BS02978");
        checkObject(bioSample);
        
        // Check Treatment
        Treatment treatment = session.get(Treatment.class, "TRT1017");
        checkObject(treatment);
        
        AdverseEvent adverseEvent = session.get(AdverseEvent.class, "AE13477");
        checkObject(adverseEvent);
        
        Intervention intervention = session.get(Intervention.class, "SM100009");
        checkObject(intervention);
        
        ImmuneExposure exposure = session.get(ImmuneExposure.class, "IM194");
        checkObject(exposure);
        
        session.close();
        sf.close();
    }
    
    private SessionFactory createSessionFactory() {
        StandardServiceRegistry standardRegistry = new StandardServiceRegistryBuilder().configure("hibernate.cfg.xml").build();
        Metadata metaData = new MetadataSources(standardRegistry).getMetadataBuilder().build();
        SessionFactory sessionFactory = metaData.getSessionFactoryBuilder().build();
        return sessionFactory;
    }
    
    private void checkObject(Object obj) throws Exception {
        System.out.println("\nCheck " + obj.getClass().getName() + ": ");
        BeanInfo info = Introspector.getBeanInfo(obj.getClass());
        PropertyDescriptor[] propDesces = info.getPropertyDescriptors();
        for (PropertyDescriptor propDesc : propDesces) {
            Method method = propDesc.getReadMethod();
            System.out.println(propDesc.getDisplayName() + ": " + method.invoke(obj));
        }
    }
    
}