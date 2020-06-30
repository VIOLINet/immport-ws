package org.reactome.immport.ws.model.queries;

public class VOToGSM {

	String gsm;
	String exposureMaterialId;
	String gender;
	String time;
	String studyAccession;
	String expSampleRepositoryName;

	
	public VOToGSM(String gsm, String exposureMaterialId, String gender, String time, String studyAccession, String expSampleRepositoryName) {
		this.gsm = gsm;
		this.exposureMaterialId = exposureMaterialId;
		this.gender = gender;
		this.time = time;
		this.studyAccession = studyAccession;
		this.expSampleRepositoryName = expSampleRepositoryName;
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

	public String getStudyAccession() {
		return studyAccession;
	}

	public void setStudyAccession(String studyAccession) {
		this.studyAccession = studyAccession;
	}

	public String getExpSampleRepositoryName() {
		return expSampleRepositoryName;
	}

	public void setExpSampleRepositoryName(String expSampleRepositoryName) {
		this.expSampleRepositoryName = expSampleRepositoryName;
	}
}
