package org.reactome.immport.ws.service;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.apache.log4j.Logger;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

import org.reactome.immport.ws.service.ImmportServiceConfig;

@Component
public class ImmportService {
	private final Logger logger = Logger.getLogger(ImmportService.class);
	
	@Autowired
	private ImmportServiceConfig config;
	
	private String biosampleMetadataString = "";
	
	public ImmportService() {}
	
	public String getBiosampleMetadata() {
		if(biosampleMetadataString == null || biosampleMetadataString == "") loadBiosampleMetatdata();
		
		return biosampleMetadataString;
	}
	
	private void loadBiosampleMetatdata() {
		File file = new File(config.getBiosampleMetatdataFileLocation());
		if(!file.exists()) return;
				
		try (BufferedReader reader = new BufferedReader(new FileReader(config.getBiosampleMetatdataFileLocation()))) {
			String line;
			while((line = reader.readLine()) != null) {
				biosampleMetadataString += line.replaceAll("\"", "") + "\n";
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
