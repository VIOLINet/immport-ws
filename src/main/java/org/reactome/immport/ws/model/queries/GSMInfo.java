package org.reactome.immport.ws.model.queries;

public class GSMInfo {

	private String studyAccesssion;
	private String subjectId;
	private String repositoryAccession;
	private Float studyTimeCollected;
	private String studyTimeCollectedUnit;
	private String expSampleAccession;
	private boolean expToMultipleGSM;
	
	public GSMInfo(String studyAccesssion, String subjectId, String repositoryAccession, Float studyTimeCollected,
			String studyTimeCollectedUnit, String expSampleAccession, boolean expToMultipleGSM) {
		super();
		this.studyAccesssion = studyAccesssion;
		this.subjectId = subjectId;
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
	
	/**
	 * Overries toString to return a tab separated String of the object's data
	 */
	@Override
	public String toString() {
		return this.getStudyAccesssion() + "\t" +
			   this.getSubjectId() + "\t" + 
			   this.getRepositoryAccession() + "\t" +
			   this.getStudyTimeCollected() + "\t" +
			   this.getStudyTimeCollectedUnit() + "\t" +
			   this.getExpSampleAccession() + "\t" +
			   this.isExpToMultipleGSM();
	}
	
	
}
