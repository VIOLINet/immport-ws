package org.reactome.immport.ws.service;

/**
 * 
 * @author brunsont
 *
 */
public class ImmportServiceConfig {

	private String biosampleMetatdataFileLocation;
	
	public ImmportServiceConfig() {
		
	}

	public String getBiosampleMetatdataFileLocation() {
		return biosampleMetatdataFileLocation;
	}

	public void setBiosampleMetatdataFileLocation(String biosampleMetatdataFileLocation) {
		this.biosampleMetatdataFileLocation = biosampleMetatdataFileLocation;
	}
	
}
