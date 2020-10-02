package org.reactome.immport.ws.model.wrapper;

import com.fasterxml.jackson.annotation.JsonProperty;

/**
 * 
 * @author brunsont
 *
 */
public class BiosampleObject {
	@JsonProperty("gsm")
	private String gsm;
	@JsonProperty("gpl")
	private String gpl;
	@JsonProperty("gse")
	private String gse;
	@JsonProperty("platform_ame")
	private String platformName;
	@JsonProperty("platform_desc")
	private String platformDesc;
	@JsonProperty("immport_study_accession")
	private String studyAccession;
	@JsonProperty("immport_arm_accession")
	private String armAccession;
	@JsonProperty("immport_subject_accession")
	private String subjectAccession;
	@JsonProperty("immport_biosample_accession")
	private String biosampleAccession;
	@JsonProperty("immport_immune_exposure_materia_id")
	private String immuneExposureMaterialId;
	@JsonProperty("vo_vaccine_name")
	private String vaccineName;
	@JsonProperty("vo_vaccine_type")
	private String vaccineType;
	@JsonProperty("immport_vaccination_time")
	private Double vaccinationTime;
	@JsonProperty("immport_vaccination_time_unit")
	private String vaccinationTimeUnit;
	@JsonProperty("race")
	private String race;
	@JsonProperty("pathogen")
	private String pathogen;
	@JsonProperty("study_design")
	private String studyDesign;
	@JsonProperty("day_0_def")
	private String day0Def;
	@JsonProperty("has_metadata")
	private String hasMetadata;
	@JsonProperty("gender")
	private String gender;
	
	public BiosampleObject() {
	}

	public BiosampleObject(String gsm, String gpl, String gse, String platformName, String platformDesc,
			String studyAccession, String armAccession, String subjectAccession, String biosampleAccession,
			String immuneExposureMaterialId, String vaccineName, String vaccineType, Double vaccinationTime,
			String vaccinationTimeUnit, String race, String pathogen, String studyDesign, String day0Def,
			String hasMetadata, String gender) {
		super();
		this.gsm = gsm;
		this.gpl = gpl;
		this.gse = gse;
		this.platformName = platformName;
		this.platformDesc = platformDesc;
		this.studyAccession = studyAccession;
		this.armAccession = armAccession;
		this.subjectAccession = subjectAccession;
		this.biosampleAccession = biosampleAccession;
		this.immuneExposureMaterialId = immuneExposureMaterialId;
		this.vaccineName = vaccineName;
		this.vaccineType = vaccineType;
		this.vaccinationTime = vaccinationTime;
		this.vaccinationTimeUnit = vaccinationTimeUnit;
		this.race = race;
		this.pathogen = pathogen;
		this.studyDesign = studyDesign;
		this.day0Def = day0Def;
		this.hasMetadata = hasMetadata;
		this.gender = gender;
	}

	public String getGsm() {
		return gsm;
	}

	public void setGsm(String gsm) {
		this.gsm = gsm;
	}

	public String getGpl() {
		return gpl;
	}

	public void setGpl(String gpl) {
		this.gpl = gpl;
	}

	public String getGse() {
		return gse;
	}

	public void setGse(String gse) {
		this.gse = gse;
	}

	public String getPlatformName() {
		return platformName;
	}

	public void setPlatformName(String platformName) {
		this.platformName = platformName;
	}

	public String getPlatformDesc() {
		return platformDesc;
	}

	public void setPlatformDesc(String platformDesc) {
		this.platformDesc = platformDesc;
	}

	public String getStudyAccession() {
		return studyAccession;
	}

	public void setStudyAccession(String studyAccession) {
		this.studyAccession = studyAccession;
	}

	public String getArmAccession() {
		return armAccession;
	}

	public void setArmAccession(String armAccession) {
		this.armAccession = armAccession;
	}

	public String getSubjectAccession() {
		return subjectAccession;
	}

	public void setSubjectAccession(String subjectAccession) {
		this.subjectAccession = subjectAccession;
	}

	public String getBiosampleAccession() {
		return biosampleAccession;
	}

	public void setBiosampleAccession(String biosampleAccession) {
		this.biosampleAccession = biosampleAccession;
	}

	public String getImmuneExposureMaterialId() {
		return immuneExposureMaterialId;
	}

	public void setImmuneExposureMaterialId(String immuneExposureMaterialId) {
		this.immuneExposureMaterialId = immuneExposureMaterialId;
	}

	public String getVaccineName() {
		return vaccineName;
	}

	public void setVaccineName(String vaccineName) {
		this.vaccineName = vaccineName;
	}

	public String getVaccineType() {
		return vaccineType;
	}

	public void setVaccineType(String vaccineType) {
		this.vaccineType = vaccineType;
	}

	public Double getVaccinationTime() {
		return vaccinationTime;
	}

	public void setVaccinationTime(Double vaccinationTime) {
		this.vaccinationTime = vaccinationTime;
	}

	public String getVaccinationTimeUnit() {
		return vaccinationTimeUnit;
	}

	public void setVaccinationTimeUnit(String vaccinationTimeUnit) {
		this.vaccinationTimeUnit = vaccinationTimeUnit;
	}

	public String getRace() {
		return race;
	}

	public void setRace(String race) {
		this.race = race;
	}

	public String getPathogen() {
		return pathogen;
	}

	public void setPathogen(String pathogen) {
		this.pathogen = pathogen;
	}

	public String getStudyDesign() {
		return studyDesign;
	}

	public void setStudyDesign(String studyDesign) {
		this.studyDesign = studyDesign;
	}

	public String getDay0Def() {
		return day0Def;
	}

	public void setDay0Def(String day0Def) {
		this.day0Def = day0Def;
	}

	public String getHasMetadata() {
		return hasMetadata;
	}

	public void setHasMetadata(String hasMetadata) {
		this.hasMetadata = hasMetadata;
	}

	public String getGender() {
		return gender;
	}

	public void setGender(String gender) {
		this.gender = gender;
	}	
}
