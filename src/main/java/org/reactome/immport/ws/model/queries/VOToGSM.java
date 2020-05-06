package org.reactome.immport.ws.model.queries;

public class VOToGSM {

	String gsm;
	String exposureMaterialId;
	String gender;
	String time;

	
	public VOToGSM(String gsm, String exposureMaterialId, String gender, String time) {
		this.gsm = gsm;
		this.exposureMaterialId = exposureMaterialId;
		this.gender = gender;
		this.time = time;
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

	public String getGender() {
		return gender;
	}

	public void setGender(String gender) {
		this.gender = gender;
	}

	public String getTime() {
		return time;
	}

	public void setTime(String time) {
		this.time = time;
	}
}
