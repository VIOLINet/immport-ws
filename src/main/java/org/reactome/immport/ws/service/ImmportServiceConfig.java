package org.reactome.immport.ws.service;

/**
 * 
 * @author brunsont
 *
 */
public class ImmportServiceConfig {

	private String biosampleMetatdataFileLocation;
	private String testDiffGeneExpFileLocation;


	public ImmportServiceConfig() {
	}
	
	public String getTestDiffGeneExpFileLocation() {
		return testDiffGeneExpFileLocation;
	}

	public void setTestDiffGeneExpFileLocation(String testDiffGeneExpFileLocation) {
		this.testDiffGeneExpFileLocation = testDiffGeneExpFileLocation;
	}
	
	public String getBiosampleMetatdataFileLocation() {
		return biosampleMetatdataFileLocation;
	}

	public void setBiosampleMetatdataFileLocation(String biosampleMetatdataFileLocation) {
		this.biosampleMetatdataFileLocation = biosampleMetatdataFileLocation;
	}
	
}
