package org.reactome.immport.ws.main;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.boot.Metadata;
import org.hibernate.boot.MetadataSources;
import org.hibernate.boot.registry.StandardServiceRegistry;
import org.hibernate.boot.registry.StandardServiceRegistryBuilder;
import org.junit.Test;
import org.reactome.immport.ws.config.AppConfig;
import org.reactome.immport.ws.model.BioSample;
import org.reactome.immport.ws.model.ExpSample;
import org.reactome.immport.ws.model.ImmuneExposure;
import org.reactome.immport.ws.model.SampleGeneExpression;
import org.reactome.immport.ws.model.Study;
import org.reactome.immport.ws.model.Subject;
import org.springframework.context.annotation.AnnotationConfigApplicationContext;

import com.fasterxml.jackson.databind.ObjectMapper;

/**
 * This class is used to handle the meta data file.
 * @author wug
 *
 */
public class BioSampleMetaDataProcessor {
    private final String SAMPLE_META_FILE = "src/main/resources/data/biosample_metadata.csv";
    
    public BioSampleMetaDataProcessor() {
    }
    
    /**
     * This method is used to dump gene expression data related to immune exposure.
     * @throws Exception
     */
    @Test
    public void dumpGeneExpression() throws Exception {
        SessionFactory sf = createSessionFactory();
        Session session = sf.openSession();
        
        String fileName = "output/ImmuneExposureGeneExpression_060921.txt";
        PrintWriter pr = new PrintWriter(fileName);
        StringBuilder builder = new StringBuilder();
        builder.append("Study Id\t"
                  + "Study Title\t"
                  + "Study Brief Desc\t"
                  + "Study Min Age\t"
                  + "Study Max Age\t"
                  + "Subject Id\t"
                  + "Subject Gender\t"
                  + "Subject Race\t"
                  + "Exposure Material Id\t"
                  + "Exposure Material Reported\t"
                  + "Biosample Id\t"
                  + "Biosample Study Time Collected\t"
                  + "Biosample Study Time Collected Units\t"
                  + "Biosample Time t0 Event\t"
                  + "Biosample Type\t"
                  + "Biosample Subtype\t"
                  + "ExpSample Id\t"
                  + "Expsample Repository Name\t"
                  + "Expsample Repository Accession\t"
                  + "Expsample To Multiple GSM Flag");
        pr.println(builder.toString());
        builder.setLength(0);
        List<ExpSample> expSamples = session.createQuery("FROM " + ExpSample.class.getName() + " es WHERE es.accession != null", 
                                                                   ExpSample.class)
                                            .getResultList();
        System.out.println("Total ImmuneExposure: " + expSamples.size());
        for (ExpSample expsample : expSamples) {
            Set<SampleGeneExpression> expressions = expsample.getGeneExpressions();
            if (expressions == null || expressions.size() == 0)
                continue;
            boolean flag = expressions.size() > 1;
            Set<BioSample> biosamples = expsample.getBioSamples();
            for (BioSample biosample : biosamples) {
                Subject subject = biosample.getSubject();
                Study study = biosample.getStudy();
                String t0event = biosample.getStudyTimeT0Event();
                if (t0event != null && t0event.equals("Other"))
                    t0event = biosample.getStudyTimeT0EventSpecify();
                Set<ImmuneExposure> immuneExposures = subject.getImmuneExposures();
                if (immuneExposures == null || immuneExposures.size() == 0)
                    continue;
                Set<String> exposureTexts = normalizeImmuneExposures(immuneExposures);
                for (String exposureText : exposureTexts) {
                    for (SampleGeneExpression expression : expressions) {
                        builder.append(study.getAccession()).append("\t")
                               .append(study.getBriefTitle()).append("\t")
                               .append(study.getBriefDescription()).append("\t")
                               .append(study.getMinimumAge()).append("\t")
                               .append(study.getMaximumAge()).append("\t")
                               .append(subject.getAccession()).append("\t")
                               .append(subject.getGender()).append("\t")
                               .append(subject.getRace()).append("\t")
                               .append(exposureText).append("\t")
                               .append(biosample.getAccession()).append("\t")
                               .append(biosample.getStudyTimeCollected()).append("\t")
                               .append(biosample.getStudyTimeCollectedUnit()).append("\t")
                               .append(t0event).append("\t")
                               .append(biosample.getType()).append("\t")
                               .append(biosample.getSubtype()).append("\t")
                               .append(expsample.getAccession()).append("\t")
                               .append(expression.getRepository().getName()).append("\t")
                               .append(expression.getRepositoryAccession()).append("\t")
                               .append(flag);
                        pr.println(builder.toString());
                        builder.setLength(0);
                    }
                }
            }
        }
        pr.close();
        session.close();
        sf.close();
    }
    
    private SessionFactory createSessionFactory() {
        StandardServiceRegistry standardRegistry = new StandardServiceRegistryBuilder().configure("hibernate.cfg.xml").build();
        Metadata metaData = new MetadataSources(standardRegistry).getMetadataBuilder().build();
        SessionFactory sessionFactory = metaData.getSessionFactoryBuilder().build();
        return sessionFactory;
    }
    
    /**
     * The same immune exposure may be used in different studies and therefore has different immune exposures ids,
     * resulting in multiple objects. This method is used to extract the common information among these objects.
     * @param exposures
     * @return
     */
    private Set<String> normalizeImmuneExposures(Set<ImmuneExposure> exposures) {
        return exposures.stream()
                        .filter(exposure -> exposure.getExposureMaterialId() != null)
                        .map(exposure -> exposure.getExposureMaterialId() + "\t" + exposure.getExposureMaterialReported())
                        .collect(Collectors.toSet());
    }
    
    /**
     * Generate some sample JSON test file for R script.
     * @throws IOException
     */
    @Test
    public void generateJSON() throws IOException {
        String study = "SDY269";
        Set<String> times = Stream.of("0", "7").collect(Collectors.toSet());
        Set<String> platforms = Stream.of("[HT_HG-U133A] Affymetrix HT Human Genome U133A Array").collect(Collectors.toSet());
        Set<String> vaccines = Stream.of("FluMist").collect(Collectors.toSet());
        Set<String> typeSubtypes = Stream.of("PBMC|CD19 + B cells isolated from PBMCs").collect(Collectors.toSet());
        List<String> gsmList = Files.lines(Paths.get(SAMPLE_META_FILE))
                                    .skip(1)
                                    .map(line -> line.split(","))
                                    .filter(tokens -> tokens[5].equals(study))
                                    .filter(tokens -> times.contains(tokens[12]))
                                    .filter(tokens -> platforms.contains(tokens[4]))
                                    .filter(tokens -> vaccines.contains(tokens[10]))
                                    .filter(tokens -> typeSubtypes.contains(tokens[22]))
                                    .map(tokens -> tokens[0])
                                    .collect(Collectors.toList());
        System.out.println("Samples: " + gsmList.size());
        ObjectMapper mapper = new ObjectMapper();
        mapper.writeValue(System.out, gsmList);
    }
    
    /**
     * This method is used to add the type and subtype information to the existing sample meta file
     * @throws Exception
     */
    @Test
    public void addTypeInfo() throws Exception {
        Map<String, String> gsmToType = new HashMap<>();
        Map<String, String> gsmToSubType = new HashMap<>();
        queryTypeInfo(gsmToType, gsmToSubType);
        addTypeInfo(gsmToType, gsmToSubType);
    }
    
    private void addTypeInfo(Map<String, String> gsmToType,
                             Map<String, String> gsmToSubtype) throws Exception {
        String dir = "src/main/resources/data/";
        String srcFile = dir + "biosample_metadata3.csv";
        String destFile = dir + "biosample_metadata4_052621.csv";
        FileReader fr = new FileReader(srcFile);
        BufferedReader br = new BufferedReader(fr);
        PrintWriter pr = new PrintWriter(destFile);
        String line = br.readLine();
        pr.println(line + ",type,subtype,type_subtype");
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split(",");
            String gsmId = tokens[0];
            String type = gsmToType.get(gsmId);
            if (type == null)
                type = "NA";
            String subType = gsmToSubtype.get(gsmId);
            if (subType == null)
                subType =  "NA";
            line = line + "," + type + "," + subType + "," + type + " : " + subType;
            pr.println(line);
        }
        br.close();
        fr.close();
        pr.close();
    }

    private void queryTypeInfo(Map<String, String> gsmToType, Map<String, String> gsmToSubType) {
        // Collection type and subtype information for GSM ids
        AnnotationConfigApplicationContext context = new AnnotationConfigApplicationContext(AppConfig.class);
        SessionFactory sf = context.getBean(SessionFactory.class);
        Session session = sf.openSession();
        List<BioSample> biosamples = session.createQuery("FROM " + BioSample.class.getName(), BioSample.class).getResultList();
        System.out.println("Total BioSamples: " + biosamples.size());
        biosamples.forEach(biosample -> {
            String type = biosample.getType();
            String subtype = biosample.getSubtype();
            if (type == null && subtype == null)
                return; // Nothing to do
            Set<ExpSample> expSamples = biosample.getExpSamples();
            if (expSamples == null || expSamples.size() == 0)
                return;
            for (ExpSample expSample : expSamples) {
                Set<SampleGeneExpression> geneExpression = expSample.getGeneExpressions();
                if (geneExpression == null || geneExpression.size() == 0)
                    continue;
                geneExpression.forEach(exp -> {
                    if (type != null)
                        gsmToType.put(exp.getRepositoryAccession(), type);
                    if (subtype != null)
                        gsmToSubType.put(exp.getRepositoryAccession(), subtype);
                });
            }
        });
        session.close();
        sf.close();
        context.close();
        // Do a quick check
        Set<String> types = gsmToType.values().stream().collect(Collectors.toSet());
        System.out.println("Total types: " + types.size());
        types.stream().sorted().forEach(System.out::println);
        Set<String> subtypes = gsmToSubType.values().stream().collect(Collectors.toSet());
        System.out.println("Total sub-types: " + subtypes.size());
        subtypes.stream().sorted().forEach(System.out::println);
    }

}
