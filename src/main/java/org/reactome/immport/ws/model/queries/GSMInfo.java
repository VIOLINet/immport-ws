package org.reactome.immport.ws.model.queries;

public class GSMInfo {

	private String studyAccesssion;
	private String studyMinAge;
	private String studyMaxAge;
	private String subjectId;
	private String subjectGender;
	private String subjectRace;
	private String voExposureId;
	private String repositoryAccession;
	private Float studyTimeCollected;
	private String studyTimeCollectedUnit;
	private String expSampleAccession;
	private boolean expToMultipleGSM;
	
	public GSMInfo(String studyAccesssion, String studyMinAge, String studyMaxAge, String subjectId, String subjectGender, String subjectRace, String voExposureId, String repositoryAccession, Float studyTimeCollected,
			String studyTimeCollectedUnit, String expSampleAccession, boolean expToMultipleGSM) {
		super();
		this.studyAccesssion = studyAccesssion;
		this.studyMinAge = studyMinAge;
		this.studyMaxAge = studyMaxAge;
		this.subjectId = subjectId;
		this.subjectGender = subjectGender;
		this.subjectRace = subjectRace;
		this.voExposureId = voExposureId;
		this.repositoryAccession = repositoryAccession;
		this.studyTimeCollected = studyTimeCollected;
		this.studyTimeCollectedUnit = studyTimeCollectedUnit;
		this.expSampleAccession = expSampleAccession;
		this.expToMultipleGSM = expToMultipleGSM;
	}
	public String getStudyAccesssion() {
		return studyAccesssion;
	}
	public void setStudyAccesssion(String studyAccesssion) {
		this.studyAccesssion = studyAccesssion;
	}
	public String getSubjectId() {
		return subjectId;
	}
	public void setSubjectId(String subjectId) {
		this.subjectId = subjectId;
	}
	public String getRepositoryAccession() {
		return repositoryAccession;
	}
	public void setRepositoryAccession(String repositoryAccession) {
		this.repositoryAccession = repositoryAccession;
	}
	public Float getStudyTimeCollected() {
		return studyTimeCollected;
	}
	public void setStudyTimeCollected(Float studyTimeCollected) {
		this.studyTimeCollected = studyTimeCollected;
	}
	public String getStudyTimeCollectedUnit() {
		return studyTimeCollectedUnit;
	}
	public void setStudyTimeCollectedUnit(String studyTimeCollectedUnit) {
		this.studyTimeCollectedUnit = studyTimeCollectedUnit;
	}
	public String getExpSampleAccession() {
		return expSampleAccession;
	}
	public void setExpSampleAccession(String expSampleAccession) {
		this.expSampleAccession = expSampleAccession;
	}
	public boolean isExpToMultipleGSM() {
		return expToMultipleGSM;
	}
	public void setExpToMultipleGSM(boolean expToMultipleGSM) {
		this.expToMultipleGSM = expToMultipleGSM;
	}
	public String getStudyMinAge() {
		return studyMinAge;
	}
	public void setStudyMinAge(String studyMinAge) {
		this.studyMinAge = studyMinAge;
	}
	public String getStudyMaxAge() {
		return studyMaxAge;
	}
	public void setStudyMaxAge(String studyMaxAge) {
		this.studyMaxAge = studyMaxAge;
	}
	public String getSubjectGender() {
		return subjectGender;
	}
	public void setSubjectGender(String subjectGender) {
		this.subjectGender = subjectGender;
	}
	public String getSubjectRace() {
		return subjectRace;
	}
	public void setSubjectRace(String subjectRace) {
		this.subjectRace = subjectRace;
	}
	public String getVoExposureId() {
		return voExposureId;
	}
	public void setVoExposureId(String voExposureId) {
		this.voExposureId = voExposureId;
	}
	/**
	 * Overrides toString to return a tab separated String of the object's data
	 */
	@Override
	public String toString() {
		return this.getStudyAccesssion() + "\t" +
			   this.getStudyMinAge() + "\t" + 
			   this.getStudyMaxAge() + "\t" +
			   this.getSubjectId() + "\t" +
			   this.getSubjectGender() + "\t" +
			   this.getSubjectRace() + "\t" +
			   this.getVoExposureId() + "\t" +
			   this.getRepositoryAccession() + "\t" +
			   this.getStudyTimeCollected() + "\t" +
			   this.getStudyTimeCollectedUnit() + "\t" +
			   this.getExpSampleAccession() + "\t" +
			   this.isExpToMultipleGSM();
	}
	
	
}
