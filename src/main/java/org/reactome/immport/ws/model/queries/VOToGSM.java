package org.reactome.immport.ws.model.queries;

public class VOToGSM {

	String studyId;
	String gsm;
	String exposureMaterialId;
	String studyTimeCollected;
	String studyTimeCollectedUnit;
	String race;
	String gender;
	
	public VOToGSM(String studyId, String gsm, String exposureMaterialId, String studyTimeCollected, String studyTimeCollectedUnit, String race, String gender) {
		this.studyId = studyId;
		this.gsm = gsm;
		this.exposureMaterialId = exposureMaterialId;
		this.studyTimeCollected = studyTimeCollected;
		this.studyTimeCollectedUnit = studyTimeCollectedUnit;
		this.race = race;
		this.gender = gender;
	}

	public String getGsm() {
		return gsm;
	}

	public void setGsm(String gsm) {
		this.gsm = gsm;
	}

	public String getRace() {
		return race;
	}

	public void setRace(String race) {
		this.race = race;
	}

	public String getGender() {
		return gender;
	}

	public void setGender(String gender) {
		this.gender = gender;
	}

	public String getStudyTimeCollected() {
		return studyTimeCollected;
	}

	public void setStudyTimeCollected(String studyTimeCollected) {
		this.studyTimeCollected = studyTimeCollected;
	}

	public String getStudyTimeCollectedUnit() {
		return studyTimeCollectedUnit;
	}

	public void setStudyTimeCollectedUnit(String studyTimeCollectedUnit) {
		this.studyTimeCollectedUnit = studyTimeCollectedUnit;
	}

	public String getExposureMaterialId() {
		return exposureMaterialId;
	}

	public void setExposureMaterialId(String exposureMaterialId) {
		this.exposureMaterialId = exposureMaterialId;
	}

	public String getStudyId() {
		return studyId;
	}

	public void setStudyId(String studyId) {
		this.studyId = studyId;
	}
}
