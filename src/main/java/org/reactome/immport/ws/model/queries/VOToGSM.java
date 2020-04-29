package org.reactome.immport.ws.model.queries;

public class VOToGSM {

	String gsm;
	String exposureMaterialId;

	
	public VOToGSM(String gsm, String exposureMaterialId) {
		this.gsm = gsm;
		this.exposureMaterialId = exposureMaterialId;
	}

	public String getGsm() {
		return gsm;
	}

	public void setGsm(String gsm) {
		this.gsm = gsm;
	}

	public String getExposureMaterialId() {
		return exposureMaterialId;
	}

	public void setExposureMaterialId(String exposureMaterialId) {
		this.exposureMaterialId = exposureMaterialId;
	}
}
