package org.reactome.immport.ws.model.queries;

public class BioSampleCollectionTime {

	private double time;
	private String units;
	
	public BioSampleCollectionTime(double time, String units) {
		this.time = time;
		this.units = units;
	}

	public double getTime() {
		return time;
	}

	public void setTime(double time) {
		this.time = time;
	}

	public String getUnits() {
		return units;
	}

	public void setUnits(String units) {
		this.units = units;
	}	
}
