package org.reactome.immport.ws.main;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
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

import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.StringColumn;
import tech.tablesaw.api.Table;
import tech.tablesaw.columns.Column;
import tech.tablesaw.io.csv.CsvReadOptions;
import tech.tablesaw.io.csv.CsvWriteOptions;
import tech.tablesaw.selection.Selection;

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
     * This method is used to generate some stats from the meta file.
     * @throws IOException
     */
    @Test
    public void calculateStats() throws IOException {
    	// As of May 9, 2022, SAMPLE_META_FILE should point to output/ImmuneExposureGeneExpression_020922.csv
    	String src = SAMPLE_META_FILE;
        FileReader reader = new FileReader(src);
        Iterable<CSVRecord> records = CSVFormat.DEFAULT.withFirstRecordAsHeader().parse(reader);
        // Studies, subjects, min age and max age, vaccine ids, biosamples, min days and max days, type and subtype,
    	// biosamples, gsm, gse, gpl
    	Set<String> studies = new HashSet<String>();
    	Set<String> subjects = new HashSet<>();
    	double minAge = Double.MAX_VALUE;
    	double maxAge = Double.MIN_VALUE;
    	Set<String> vaccines = new HashSet<>();
    	double minDays = Double.MAX_VALUE;
    	double maxDays = Double.MIN_VALUE;
    	Set<String> types = new HashSet<>();
    	Set<String> subtypes = new HashSet<>();
    	Set<String> biosamples = new HashSet<>();
    	Set<String> gsms = new HashSet<>();
    	Set<String> gses = new HashSet<>();
    	Set<String> gpls = new HashSet<>();
    	Set<String> races = new HashSet<>();
    	int totalRows = 0;
    	// For individual vaccine
    	Map<String, VaccineRelatedStats> vaccine2info = new HashMap<>();
    	for (CSVRecord record : records) {
    		totalRows ++;
    		studies.add(record.get("immport_study_accession"));
    		subjects.add(record.get("immport_subject_accession"));
    		double minAge1 = Double.parseDouble(record.get("study_min_age"));
    		if (minAge1 < minAge) minAge = minAge1;
    		double maxAge1 = Double.parseDouble(record.get("study_max_age"));
    		if (maxAge1 > maxAge) maxAge = maxAge1;
    		vaccines.add(record.get("immport_immune_exposure_material_id"));
    		double day = Double.parseDouble(record.get("immport_vaccination_time"));
    		if (day < minDays) minDays = day;
    		if (day > maxDays) maxDays = day;
    		types.add(record.get("biosample_type"));
    		subtypes.add(record.get("subtype"));
    		biosamples.add(record.get("immport_biosample_accession"));
    		gsms.add(record.get("gsm"));
    		gses.add(record.get("gse"));
    		gpls.add(record.get("gpl"));
    		races.add(record.get("race"));
    		// Collect the vaccine specific information
    		VaccineRelatedStats vaccineStats = vaccine2info.get(record.get("vaccine"));
    		if (vaccineStats == null) {
    			vaccineStats = new VaccineRelatedStats();
    			vaccineStats.vaccine = record.get("vaccine");
    			vaccineStats.vaccineId = record.get("immport_immune_exposure_material_id");
    			vaccine2info.put(vaccineStats.vaccine, vaccineStats);
    		}
    		vaccineStats.studies.add(record.get("immport_study_accession"));
    		vaccineStats.subjects.add(record.get("immport_subject_accession"));
    		if (minAge1 < vaccineStats.minAge) vaccineStats.minAge = minAge1;
    		if (maxAge1 > vaccineStats.maxAge) vaccineStats.maxAge = maxAge1;
    		if (day < vaccineStats.minDays) vaccineStats.minDays = day;
    		if (day > vaccineStats.maxDays) vaccineStats.maxDays = day;
    		vaccineStats.cellTypes.add(record.get("biosample_type"));
    		vaccineStats.cellSubTypes.add(record.get("subtype"));
    		vaccineStats.biosamples.add(record.get("immport_biosample_accession"));
    		vaccineStats.gsms.add(record.get("gsm"));
    		vaccineStats.gses.add(record.get("gse"));
    		vaccineStats.gpls.add(record.get("gpl"));
    		vaccineStats.races.add(record.get("race"));
    	}
    	System.out.println("Total studies: " + studies.size());
    	System.out.println("Total subjects: " + subjects.size());
    	System.out.println("Min age: " + minAge);
    	System.out.println("Max age: " + maxAge);
    	System.out.println("Total vaccines: " + vaccines.size());
    	System.out.println("Min day: " + minDays);
    	System.out.println("Max day: " + maxDays);
    	System.out.println("Total types: " + types.size());
    	System.out.println("Total subtypes: " + subtypes.size());
    	System.out.println("Total biosamples: " + biosamples.size());
    	System.out.println("Total GSMs: " + gsms.size());
    	System.out.println("Total GSEs: " + gses.size());
    	System.out.println("Total GPLs: " + gpls.size());
    	System.out.println("Total races: " + races.size());
    	System.out.println("Total rows: " + totalRows);
    	
    	System.out.println("\nVaccine based stats:");
    	List<String> vaccineList = vaccine2info.keySet().stream().sorted().collect(Collectors.toList());
    	StringBuilder builder = new StringBuilder();
    	String headers = "Vaccine\tVO_ID\tCategories\t" + 
    	                 "Races\tminAge\tMaxAge\tminDays\t" +
    			         "maxDays\tCell_Types\tCell_Subtypes\t" + 
    	                 "Biosamples\tGSMs\tGSEs\tGPLs\tStudies";
    	System.out.println(headers);
    	for (String vaccine : vaccineList) {
    		VaccineRelatedStats stats = vaccine2info.get(vaccine);
    		builder.append(vaccine).append("\t");
    		builder.append(stats.vaccineId).append("\t");
    		builder.append("\t"); // To be filled manually
    		builder.append(stats.races.size()).append("\t");
    		builder.append(stats.minAge).append("\t");
    		builder.append(stats.maxAge).append("\t");
    		builder.append(stats.minDays).append("\t");
    		builder.append(stats.maxDays).append("\t");
    		builder.append(stats.cellTypes.size()).append("\t");
    		builder.append(stats.cellSubTypes.size()).append("\t");
    		builder.append(stats.biosamples.size()).append("\t");
    		builder.append(stats.gsms.size()).append("\t");
    		builder.append(stats.gses.size()).append("\t");
    		builder.append(stats.gpls.size()).append("\t");
    		builder.append(stats.studies.size()).append("\t");
    		System.out.println(builder.toString());
    		builder.setLength(0);
    	}
    }
    
    /**
     * There are some rows use hours. For these rows, we need to convert them as days by dividing 24.
     * @throws IOException
     */
    @Test
    public void convertHoursToDays() throws IOException {
    	String src = "output/ImmuneExposureGeneExpression_091621.csv";
    	String target = "output/ImmuneExposureGeneExpression_020922.csv";
    	
    	CsvReadOptions.Builder builder = CsvReadOptions.builder(src).separator(',');
    	CsvReadOptions options = builder.build();
    	Table table = Table.read().usingOptions(options);
    	System.out.println("Total rows: " + table.rowCount());
    	System.out.println("Total cols: " + table.columnCount());
    	
    	// Check how many time units are used
    	StringColumn timeUnitCol = table.stringColumn("immport_vaccination_time_unit");
    	Set<String> timeUnits = new HashSet<>();
    	DoubleColumn timeCol = table.doubleColumn("immport_vaccination_time");
    	int totalUpdate = 0;
    	for (int i = 0; i < table.rowCount(); i++) {
    		String timeUnit = timeUnitCol.getString(i);
    		timeUnits.add(timeUnit);
    		if (timeUnit.equals("Hours")) {
    			Double time = timeCol.getDouble(i);
    			Double day = time / 24.0d;
    			// Convert it into two decimal
    			// This is just rounded since it is not easy for negative days.
    			day = ((int)(day * 100)) / 100.0d;
    			timeUnitCol.set(i, "Days");
    			timeCol.set(i, day);
    			totalUpdate ++;
    		}
    	}
    	System.out.println("Total time units: " + String.join(", ", timeUnits));
    	System.out.println("Total updated rows: " + totalUpdate);

    	// Save the updated table
    	CsvWriteOptions.Builder writerBuilder = CsvWriteOptions.builder(target)
    			.separator(',');
    	table.write().usingOptions(writerBuilder.build());
    }
    
    /**
     * Call this method to update the headers so that they are compatible to JavaScript Vue code.
     * This method also changes the delimit from tab to common.
     * @throws IOException
     */
    @Test
    public void updateHeaders() throws IOException {
    	// The following is used to map from the immport table headers to the headers used by Nasim
    	// and then coded into Vue.
    	Map<String, String> headerMap = new HashMap<String, String>();
    	headerMap.put("Expsample_Repository_Accession", "gsm");
    	headerMap.put("GPL_Accession", "gpl");
    	headerMap.put("GSE_Accession", "gse");
    	headerMap.put("GPL_Title", "platform_desc");
    	headerMap.put("Study_Id", "immport_study_accession");
    	headerMap.put("Subject_Id", "immport_subject_accession");
    	headerMap.put("Subject_Gender", "gender");
    	headerMap.put("ExpSample_Id", "immport_biosample_accession");
    	headerMap.put("Exposure_Material_Id", "immport_immune_exposure_material_id");
    	headerMap.put("Exposure_Material_Reported", "vaccine");
    	headerMap.put("Biosample_Study_Time_Collected", "immport_vaccination_time");
    	headerMap.put("Biosample_Study_Time_Collected_Units", "immport_vaccination_time_unit");
    	headerMap.put("Subject_Race", "race");
    	headerMap.put("Biosample_Time_t0_Event", "day_0_def");
    	headerMap.put("type", "Biosample_Type");
    	headerMap.put("Biosample_Subtype", "subtype");
    	
    	String srcFile = "output/ImmuneExposureGeneExpression_091421.txt";
    	String targetFile = "output/ImmuneExposureGeneExpression_091621.csv";
    	FileReader fr = new FileReader(srcFile);
    	BufferedReader br = new BufferedReader(fr);
    	PrintWriter pr = new PrintWriter(targetFile);
    	String line = br.readLine();
    	line = line.replace("\t", ",").replace(" ", "_");
    	String[] tokens = line.split(",");
    	List<String> list = new ArrayList<String>();
    	for (String token : tokens) {
    		String mapped = headerMap.get(token);
    		if (mapped == null)
    			list.add(token.toLowerCase());
    		else
    			list.add(mapped);
    	}
    	list.add("type_subtype");
    	list.add("age_group");
    	pr.println(String.join(",", list));
    	// For generating type_subtype
    	int typeCol = 14;
    	int subtypeCol = 15;
    	int minAgeCol = 3;
    	int maxAgeCol = 4;
    	
    	while ((line = br.readLine()) != null) {
    		tokens = line.split("\t");
    		list.clear();
    		for (String token : tokens) {
    			if (token.contains(","))
    				token = "\"" + token + "\"";
    			list.add(token);
    		}
    		list.add(tokens[typeCol] + " : " + 
    		         (tokens[subtypeCol].length() == 0 ? "NA" : tokens[15]));
    		// For age group
    		list.add(tokens[minAgeCol] + " - " + tokens[maxAgeCol]);
    		pr.println(String.join(",", list));
    	}
    	pr.close();
    	br.close();
    }
    
    /**
     * We don't want to include GSE datasets having only one time point since they cannot be used for
     * differential expression analysis even though they may be merged with other datasets. However,
     * integration may always be an issue during differential expression analysis. 
     * @throws IOException
     */
    @Test
    public void filterOutGSEsWithOnlyOneTimePoint() throws IOException {
    	String src = "output/ImmuneExposureGeneExpression_082521.txt";
    	String target = "output/ImmuneExposureGeneExpression_090921.txt";
    	CsvReadOptions.Builder builder = CsvReadOptions.builder(src)
    			.separator('\t');
    	CsvReadOptions options = builder.build();
    	Table table = Table.read().usingOptions(options);
    	System.out.println("Total rows: " + table.rowCount());
    	System.out.println("Total cols: " + table.columnCount());
    	StringColumn gseCol = table.stringColumn("GSE Accession");
    	Column<?> timeCol = table.column("Biosample Study Time Collected");
    	Map<String, Set<String>> gse2times = new HashMap<>();
    	for (int i = 0; i < table.rowCount(); i++) {
    		String gse = gseCol.getString(i);
    		String time = timeCol.getString(i);
    		gse2times.compute(gse, (key, set) -> {
    			if (set == null)
    				set = new HashSet<>();
    			set.add(time);
    			return set;
    		});
    	}
    	System.out.println("\nThe following GSE datasets will not be included:");
//    	gse2times.forEach((g, t) -> System.out.println(g + ": " + t));
    	Set<String> excludedGSEs = gse2times.keySet()
    							           .stream()
    							           .filter(g -> gse2times.get(g).size() == 1)
    							           .collect(Collectors.toSet());
    	// See file, 07_2021_Result_Summary.xlsx, for the reason why the following is excluded
    	excludedGSEs.add("GSE65440");
    	excludedGSEs.stream().sorted().forEach(System.out::println);
    	// Filter on the table directly
    	Selection selection = gseCol.isNotIn(excludedGSEs);
    	System.out.println("Total rows to be included: " + selection.size());
    	Table filteredTable = table.where(selection);
    	System.out.println("After filtering: " + filteredTable.rowCount());
    	CsvWriteOptions.Builder writerBuilder = CsvWriteOptions.builder(target)
    			.separator('\t');
    	filteredTable.write().usingOptions(writerBuilder.build());
    }
    
    /**
     * Add batch annotation to the meta file. The majority of rows will add GSE ids as their
     * batches. However, for some of them, a hard-coded batch annotation was used. This annotation
     * was based on an external annotation file 07_2021_Result_Summary.xlsx. 
     * @throws IOException
     */
    @Test
    public void addBatchAnnotations() throws IOException {
        // The following batch annotation is based on file 07_2021_Result_Summary.xlsx, which was
        // manually edited based on PCA analysis. 
        String batch_annotations = "GSE30101,GSM744835-GSM744888,GSM744942-GSM745143,GSM745214-GSM745344\n"
                                 + "GSE59714,GSM1441196-GSM1441267,GSM1441830-GSM1441985";
        String[] lines = batch_annotations.split("\n");
        Map<String, List<Integer[]>> gse2batches = new HashMap<>();
        for (String line : lines) {
            List<Integer[]> batches = new ArrayList<>();
            String[] tokens = line.split(",");
            for (int i = 1; i < tokens.length; i++) {
                String[] tokens1 = tokens[i].split("-");
                Integer min = new Integer(tokens1[0].substring(3));
                Integer max = new Integer(tokens1[1].substring(3));
                batches.add(new Integer[]{min, max});
            }
            gse2batches.put(tokens[0], batches);
        }
        System.out.println(gse2batches);
        
        String srcFile = "output/ImmuneExposureGeneExpression_072021.txt";
        String targetFile = "output/ImmuneExposureGeneExpression_082521.txt";
        FileReader fr = new FileReader(srcFile);
        BufferedReader br = new BufferedReader(fr);
        PrintWriter pr = new PrintWriter(targetFile);
        String line = br.readLine();
        // Add a new column
        pr.println(line + "\tBatch Factor");
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split("\t");
            // GSE
            if (gse2batches.containsKey(tokens[20])) {
                String gsmId = tokens[18];
                int batch = findBatch(gsmId, gse2batches.get(tokens[20]));
                if (batch == -1)
                    throw new IllegalStateException("Cannot find batch for " + gsmId);
                pr.println(line + "\t" + tokens[20] + "_" + batch);
            }
            else {
                // Use GSE directly
                pr.println(line + "\t" + tokens[20]);
            }
        }
        pr.close();
        br.close();
        fr.close();
    }
    
    private int findBatch(String gsmId, List<Integer[]> batches) {
        int batch = -1;
        Integer id = new Integer(gsmId.substring(3));
        for (int i = 0; i < batches.size(); i++) {
            Integer[] range = batches.get(i);
            if (id >= range[0] && id <= range[1])
                return i;
        }
        return batch;
    }
    
    
    /**
     * Some GSEs don't have pre-processed expression data in the matrix format. We will exclude this.
     * @throws IOException
     */
    @Test
    public void filterGSEs() throws IOException {
        String[] escapedGSEs = {
//                "GSE41080" // https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE41080. 91 GSMs will be filtered out.
        		"GSE22121",
        		"GSE125921",
        		"GSE62110",
        		"GSE64361",
        		"GSE64514",
        		"GSE65440"
        };
        String input = "output/ImmuneExposureGeneExpression_071921.txt";
        String output = "output/ImmuneExposureGeneExpression_072021.txt";
        
        input = "output/ImmuneExposureGeneExpression_090921.txt";
        output = "output/ImmuneExposureGeneExpression_091421.txt";
        
        int colIndex = 20;
        filterRows(input, output, escapedGSEs, colIndex);
    }
    
    /**
     * Some GPL platforms cannot be handled in the current R script: including
     * RNA-seq and peptide arrays.
     * @throws IOException
     */
    @Test
    public void filterGPLs() throws IOException {
        String[] escapePlatforms = new String[] {
                "GPL16791", // Illumina HiSeq 2500 (Homo sapiens): RNA-seq
                "GPL16497" // Utz lab influenza peptide array: peptide array
        };
        String input = "output/ImmuneExposureGeneExpression_071221.txt";
        String output = "output/ImmuneExposureGeneExpression_071921.txt";
        int colIndex = 21;
        filterRows(input, output, escapePlatforms, colIndex);
    }

    protected void filterRows(String input, String output, String[] filterKeys, int colIndex) throws FileNotFoundException, IOException {
        FileReader fr = new FileReader(input);
        BufferedReader br = new BufferedReader(fr);
        PrintWriter pr = new PrintWriter(output);
        String line = br.readLine();
        pr.println(line);
        Map<String, Integer> platform2deleted = new HashMap<>();
        Set<String> filterSet = Stream.of(filterKeys).collect(Collectors.toSet());
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split("\t");
            String platform = tokens[colIndex];
            if (!filterSet.contains(platform)) {
                // Remove some unexpected attached spaces
                line = Stream.of(tokens).map(t -> t.trim()).collect(Collectors.joining("\t"));
                pr.println(line);
                continue;
            }
            platform2deleted.compute(platform, (k, v) -> {
               if (v == null)
                   v = 0;
               v ++;
               return v;
            });
        }
        pr.close();
        br.close();
        fr.close();
        System.out.println("Filtered samples according to keys:");
        platform2deleted.forEach((k, v) -> System.out.println(k + ": " + v));
    }
    
    @Test
    public void generateUniGeneMap() throws IOException {
        String dir = "/Volumes/ssd/datasets/NCBI/UniGene/";
        String src = dir + "Hs.data";
        String dest = dir + "UniGeneMapperWithAcc_071321.txt";
        String dest1 = dir + "UniGeneMapper_071321.txt";
        FileReader fr = new FileReader(src);
        BufferedReader br = new BufferedReader(fr);
        PrintWriter pr = new PrintWriter(dest);
        PrintWriter pr1 = new PrintWriter(dest1);
        pr.println("UniGene\tGene Symbol\tGene_ID\tAccession");
        pr1.println("UniGene\tGene Symbol\tGene_ID");
        String id = null, gene = null, geneId = null, line = null;
        while ((line = br.readLine()) != null) {
            if (line.startsWith("ID ")) {
                id = line.substring("ID".length()).trim();
            }
            else if (line.startsWith("GENE ")) {
                gene = line.substring("GENE".length()).trim();
            }
            else if (line.startsWith("GENE_ID ")) {
                geneId = line.substring("GENE_ID".length()).trim();
                pr1.println(id + "\t" + 
                            (gene == null ? "" : gene) + "\t" + 
                            (geneId == null ? "" : geneId));
            }
            else if (line.startsWith("SEQUENCE    ACC=")) {
                int index = line.indexOf(".");
                String acc = line.substring("SEQUENCE    ACC=".length(), index);
                if (id == null)
                    throw new IllegalStateException(acc + " is not linked to ID!");
                pr.println(id + "\t" + 
                           (gene == null ? "" : gene) + "\t" + 
                           (geneId == null ? "" : geneId) + "\t" + 
                           acc);
            }
            else if (line.startsWith("//")) {
                id = gene = geneId = null;
            }
        }
        pr.close();
        pr1.close();
        br.close();
        fr.close();
    }
    
    @Test
    public void checkGSEOverlaps() throws IOException {
        Map<String, Set<String>> gse2gsms = new HashMap<>();
        String fileName = "output/ImmuneExposureGeneExpression_061721.txt";
        Files.lines(Paths.get(fileName))
             .skip(1)
             .map(line -> line.split("\t"))
             .forEach(tokens -> {
                 String gsm = tokens[18];
                 String gse = tokens[20];
                 gse2gsms.compute(gse, (key, set) -> {
                     if (set == null)
                         set = new HashSet<>();
                     set.add(gsm);
                     return set;
                 });
             });
        System.out.println("Total GSE accessions: " + gse2gsms.size());
        
        // Perform some pairwise overlapping analysis
        List<String> gseList = gse2gsms.keySet().stream().sorted().collect(Collectors.toList());
        // The first gse should be replaced by the second gse since the second gse contains all GSM ids in the first gse
        // and it has the larger accession number
        Map<String, String> gse2gse = new HashMap<>();
        for (int i = 0; i < gseList.size() - 1; i++) {
            String gse1 = gseList.get(i);
            Set<String> gsms1 = gse2gsms.get(gse1);
            for (int j = i + 1; j < gseList.size(); j++) {
                String gse2 = gseList.get(j);
                Set<String> gsms2 = gse2gsms.get(gse2);
                Set<String> shared = new HashSet<>(gsms2);
                shared.retainAll(gsms1);
                if (shared.size() == 0)
                    continue;
                System.out.println(gse1 + "\t" + gse2 + "\t" +
                                   gsms1.size() + "\t" + gsms2.size() + "\t" +
                                   shared.size());
                checkSubsume(gse1, gse2, gsms1, gsms2, shared, gse2gse);
            }
        }
        System.out.println();
        gse2gse.keySet().stream().sorted().forEach(g -> System.out.println(g + " -> " + gse2gse.get(g)));
        
        // Based on the above mapping to update the output
        String outFileName = "output/ImmuneExposureGeneExpression_071221.txt";
        PrintWriter pr = new PrintWriter(outFileName);
        FileReader fr = new FileReader(fileName);
        BufferedReader br = new BufferedReader(fr);
        String line = br.readLine();
        pr.println(line); // Just the header
        Set<String> processedGSMs = new HashSet<>();
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split("\t");
            String gsm = tokens[18];
            if (processedGSMs.contains(gsm))
                continue;
            processedGSMs.add(gsm);
            String gse = tokens[20];
            String mapped = getMappedGSE(gse, gse2gse);
            if (!gse.equals(mapped)) {
                tokens[20] = mapped;
                line = String.join("\t", tokens);
            }
            pr.println(line);
        }
        br.close();
        fr.close();
        pr.close();
    }
    
    private String getMappedGSE(String gse, Map<String, String> gse2gse) {
        String rtn = gse;
        while (gse2gse.containsKey(rtn)) {
            rtn = gse2gse.get(rtn);
        }
        return rtn;
    }
    
    private void checkSubsume(String gse1,
                              String gse2,
                              Set<String> gsms1,
                              Set<String> gsms2,
                              Set<String> shared, 
                              Map<String, String> gse2gse) {
        if (shared.size() == 0)
            return;
        if (gsms1.size() != shared.size() && gsms2.size() != shared.size())
            return; // There are some extra ids in one of these two sets
        if (gsms1.size() > gsms2.size())
            gse2gse.put(gse2, gse1); // The second should be replaced by the first
        else if (gsms1.size() < gsms2.size())
            gse2gse.put(gse1, gse2); // The first should be replaced by the second
        else {
            // In this case, we will pick up the id that is smaller to avoid pruning.
            // The id that is larger should contain more GSMs in GEO.
            Integer id1 = new Integer(gse1.substring(3));
            Integer id2 = new Integer(gse2.substring(3));
            if (id1 < id2)
                gse2gse.put(gse2, gse1);
            else
                gse2gse.put(gse1, gse2); // There are only two cases
        }
    }
    
    @Test
    public void checkGSE2GPLMap() throws IOException {
        Map<String, Set<String>> gse2gpl = new HashMap<>();
        String fileName = "output/ImmuneExposureGeneExpression_061721.txt";
        Files.lines(Paths.get(fileName))
             .skip(1)
             .forEach(line -> {
                 String[] tokens = line.split("\t");
                 String gse = tokens[20];
                 String gpl = tokens[21];
                 gse2gpl.compute(gse, (key, set) -> {
                     if (set == null)
                         set = new HashSet<>();
                     set.add(gpl);
                     return set;
                 });
             });
        gse2gpl.forEach((gse, gpls) -> {
            if (gpls.size() > 1)
                System.out.println(gse + "\t" + String.join(", ", gpls));
        });
        // Output from the above results:
//        GSE22121    GPL9700, GPL10465
//        GSE13699    GPL6104, GPL6883
//        GSE29619    GPL570, GPL3921, GPL13158
//        GSE18323    GPL570, GPL571
//        GSE48024    GPL10558, GPL6947
    }
    
    /**
     * Do some filtering and then map GSM ids to GSE and GPL based on the GEOSqlLite database.
     * @throws Exception
     */
    @Test
    public void filterAndMapGSMSamples() throws Exception {
        Connection conn = getGEODbConnection();
        String gseSql = "SELECT gse, gsm FROM gse_gsm WHERE gsm = ?";
        String gplSql = "SELECT gpl, gsm FROM gsm WHERE gsm = ?";
        String gplTitleSql = "SELECT title FROM gpl WHERE gpl = ?";
        PreparedStatement gseStat = conn.prepareStatement(gseSql);
        PreparedStatement gplStat = conn.prepareStatement(gplSql);
        PreparedStatement gplTitleStat = conn.prepareStatement(gplTitleSql);
        
        String srcFileName = "output/ImmuneExposureGeneExpression_060921.txt";
        String targetFilename = "output/ImmuneExposureGeneExpression_061721.txt";
        List<String> lines = Files.readAllLines(Paths.get(srcFileName));
        System.out.println("Total lines: " + lines.size());
        PrintWriter pr = new PrintWriter(targetFilename);
        String header = lines.get(0);
        pr.println(header + "\tGSE Accession\tGPL Accession\tGPL Title" );
        Set<String> escapeStudies = Stream.of("SDY406", "SDY675", "SDY773").collect(Collectors.toSet());
        Set<String> allGSEs = new HashSet<>();
        Set<String> allGPLs = new HashSet<>();
        Set<String> allStudies = new HashSet<>();
        int totalMoreThanOneGSE = 0;
        Set<Integer> moreThanOneCounter = new HashSet<>();
        int totalOutputLines = 0;
        for (int i = 1; i < lines.size(); i++) {
            String line = lines.get(i);
            String[] tokens = line.split("\t");
            // Filter 1: Make sure "Exposure Material Id" starting with VO
            if (!tokens[8].startsWith("VO"))
                continue;
            // Filter 2: Ignore three studies, which are DNA-sequences of antibody genes
            // SDY460, SDY675, SDY773
            if (escapeStudies.contains(tokens[0]))
                continue;
            // Filter 3: Make sure expression ids starting with GSM
            if (!tokens[18].startsWith("GSM"))
                continue;
            allStudies.add(tokens[0]);
            // Need to do trim since there is a space after the id in some cases
            String gsm = tokens[18].trim(); 
            String gpl = queryGPL(gsm, gplStat);
            allGPLs.add(gpl);
            String gplTitle = queryGPLTitle(gpl, gplTitleStat);
            gseStat.setString(1, gsm);
            ResultSet result = gseStat.executeQuery();
            int counter = 0;
            while (result.next()) {
                // We may get more than one GSE
                String gse = result.getString("gse");
                allGSEs.add(gse);
                pr.println(line + "\t" + 
                           gse + "\t" + 
                           gpl + "\t" + 
                           gplTitle);
                counter ++;
                totalOutputLines ++;
            }
            result.close();
            if (counter > 1) {
                System.out.println(gsm + " has " + counter + " GSEs.");
                totalMoreThanOneGSE ++;
                moreThanOneCounter.add(counter);
            }
        }
        gseStat.close();
        gplStat.close();
        gplTitleStat.close();
        conn.close();
        pr.close();
        System.out.println("Total GSM ids having more than one GSE: " + totalMoreThanOneGSE);
        System.out.println("The counter is: " + moreThanOneCounter);
        System.out.println("Total output lines: " + totalOutputLines);
        System.out.println("Total GSEs: " + allGSEs.size());
        System.out.println("Total GPLs: " + allGPLs.size());
        System.out.println("Total Studies: " + allStudies.size());
    }
    
    private String queryGPLTitle(String gpl,
                                 PreparedStatement stat) throws Exception {
        stat.setString(1, gpl);
        ResultSet result = stat.executeQuery();
        String rtn = null;
        if (result.next()) {
            rtn = result.getString(1);
        }
        result.close();
        return rtn;
    }
    
    private String queryGPL(String gsm,
                            PreparedStatement stat) throws SQLException {
        stat.setString(1, gsm);
        ResultSet results = stat.executeQuery();
        int counter = 0;
        String gpl = null;
        while (results.next()) {
            gpl = results.getString(1);
            counter ++;
        }
        results.close();
        if (gpl == null)
            throw new IllegalStateException("Cannot find a gpl for " + gsm);
        if (counter > 1)
            throw new IllegalStateException("More than one gpl found for " + gsm);
        return gpl;
    }
    
    private Connection getGEODbConnection() throws SQLException {
        String fileName = "/Users/wug/GEOmetadb.sqlite";
        String url = "jdbc:sqlite:" + fileName;
        Connection conn = DriverManager.getConnection(url);
        return conn;
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
                               .append(expression.getRepositoryAccession().trim()).append("\t")
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
    
    private class VaccineRelatedStats {
    	String vaccine;
    	String vaccineId;
    	Set<String> studies = new HashSet<>();
    	Set<String> races = new HashSet<>();
    	double minAge = Double.MAX_VALUE;
    	double maxAge = 0.0d;
    	double minDays = Double.MAX_VALUE;
    	double maxDays = -minDays;
    	Set<String> cellTypes = new HashSet<>();
    	Set<String> cellSubTypes = new HashSet<>();
    	Set<String> biosamples = new HashSet<>();
    	Set<String> subjects = new HashSet<>();
    	Set<String> gses = new HashSet<>();
    	Set<String> gsms = new HashSet<>();
    	Set<String> gpls= new HashSet<>();
    	
    	public VaccineRelatedStats() {
    	}
    }

}
