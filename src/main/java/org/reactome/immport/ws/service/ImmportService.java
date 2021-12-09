package org.reactome.immport.ws.service;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class ImmportService {
	private final Logger logger = Logger.getLogger(ImmportService.class);
	
	@Autowired
	private ImmportServiceConfig config;
	
	private String biosampleMetadataString = "";
	
	public ImmportService() {}
	
	public String getBiosampleMetadata() {
		if(biosampleMetadataString == null || biosampleMetadataString == "") loadBiosampleMetadata();
		
		return biosampleMetadataString;
	}
	
	/**
	 * Loads in biosample_metadata.csv and stores in biosampleMetadataString for caching before sending
	 * to front end. Stores string were columns are \t separated and rows are \n separated
	 */
	private void loadBiosampleMetadata() {
	    try {
	        logger.info("Starting loading " + config.getBiosampleMetatdataFileLocation() + "...");
	        InputStream is = getClass().getClassLoader().getResourceAsStream(config.getBiosampleMetatdataFileLocation());
	        InputStreamReader isr = new InputStreamReader(is);
	        Iterable<CSVRecord> records = CSVFormat.EXCEL.parse(isr);
	        for(CSVRecord record : records) {
	            String row = "";
	            for(String cell : record) {
	                row += cell + "\t";
	            }
	            biosampleMetadataString += row.replaceAll("\"", "") + "\n";
	        }
	        isr.close();
	        is.close();
	        logger.info("Done loading.");
	    } catch (FileNotFoundException e) {
	        logger.error(e);
	        e.printStackTrace();
	    } catch (IOException e) {
	        logger.error(e);
	        e.printStackTrace();
	    }
	}
	
	public String loadTestDiffExpResults() throws IOException {
		InputStream is = getClass().getClassLoader().getResourceAsStream(config.getTestDiffGeneExpFileLocation());
		InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        List<String> lines = new ArrayList<>();
        String line = null;
        while ((line = br.readLine()) != null) {
        	lines.add(line);
        }
        br.close();
        isr.close();
        is.close();
        return String.join("\n", lines);
	}
}
