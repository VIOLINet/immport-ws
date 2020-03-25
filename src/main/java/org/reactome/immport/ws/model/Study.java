package org.reactome.immport.ws.model;

import java.util.Set;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.OneToMany;
import javax.persistence.Table;

import com.fasterxml.jackson.annotation.JsonBackReference;
import com.fasterxml.jackson.annotation.JsonIgnoreProperties;
import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.annotation.JsonInclude.Include;

@Entity
@Table(name = "study")
@JsonInclude(Include.NON_NULL)
public class Study {
    @Id
    @Column(name = "study_accession")
    private String accession;
    @Column(name = "brief_title")
    private String briefTitle;
    @Column(name = "age_unit")
    private String ageUnit;
    @Column(name = "brief_description")
    private String briefDescription;
    @Column(name = "condition_studied")
    private String conditionStudied;
    private String doi;
    @Column(name = "gender_included")
    private String genderIncluded;
    @Column(name = "intervention_agent")
    private String interventionAgent;
    @Column(name = "maximum_age")
    private String maximumAge;
    @Column(name = "minimum_age")
    private String minimumAge;
    @Column(name = "official_title")
    private String officialTitle;
    @Column(name = "target_enrollment")
    private Integer targetEnrollment;
    @OneToMany(fetch = FetchType.LAZY)
    @JsonIgnoreProperties({"hibernateLazyInitializer", "handler"})
    @JoinColumn(name = "study_accession")
    @JsonBackReference // Ignore these relationships during json serilization
    private Set<Experiment> experiments;
    @OneToMany(fetch = FetchType.LAZY)
    @JoinColumn(name = "study_accession")
    @JsonIgnoreProperties({"hibernateLazyInitializer", "handler"})
    @JsonBackReference
    private Set<ArmOrCohort> arms;
    
    public Study() {
    }

    public Set<ArmOrCohort> getArms() {
        return arms;
    }

    public void setArms(Set<ArmOrCohort> arms) {
        this.arms = arms;
    }

    public Set<Experiment> getExperiments() {
        return experiments;
    }

    public void setExperiments(Set<Experiment> experiments) {
        this.experiments = experiments;
    }

    public String getAgeUnit() {
        return ageUnit;
    }

    public void setAgeUnit(String ageUnit) {
        this.ageUnit = ageUnit;
    }

    public String getBriefDescription() {
        return briefDescription;
    }

    public void setBriefDescription(String briefDescription) {
        this.briefDescription = briefDescription;
    }

    public String getConditionStudied() {
        return conditionStudied;
    }

    public void setConditionStudied(String conditionStudied) {
        this.conditionStudied = conditionStudied;
    }

    public String getDoi() {
        return doi;
    }

    public void setDoi(String doi) {
        this.doi = doi;
    }

    public String getGenderIncluded() {
        return genderIncluded;
    }

    public void setGenderIncluded(String genderIncluded) {
        this.genderIncluded = genderIncluded;
    }

    public String getInterventionAgent() {
        return interventionAgent;
    }

    public void setInterventionAgent(String interventionAgent) {
        this.interventionAgent = interventionAgent;
    }

    public String getMaximumAge() {
        return maximumAge;
    }

    public void setMaximumAge(String maximumAge) {
        this.maximumAge = maximumAge;
    }

    public String getMinimumAge() {
        return minimumAge;
    }

    public void setMinimumAge(String minimumAge) {
        this.minimumAge = minimumAge;
    }

    public String getOfficialTitle() {
        return officialTitle;
    }

    public void setOfficialTitle(String officialTitle) {
        this.officialTitle = officialTitle;
    }

    public Integer getTargetEnrollment() {
        return targetEnrollment;
    }

    public void setTargetEnrollment(Integer targetEnrollment) {
        this.targetEnrollment = targetEnrollment;
    }

    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public String getBriefTitle() {
        return briefTitle;
    }

    public void setBriefTitle(String briefTitle) {
        this.briefTitle = briefTitle;
    }

}
