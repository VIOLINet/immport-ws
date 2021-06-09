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
import org.junit.Test;
import org.reactome.immport.ws.config.AppConfig;
import org.reactome.immport.ws.model.BioSample;
import org.reactome.immport.ws.model.ExpSample;
import org.reactome.immport.ws.model.SampleGeneExpression;
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
