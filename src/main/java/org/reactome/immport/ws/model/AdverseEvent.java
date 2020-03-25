package org.reactome.immport.ws.model;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Table;

@Entity
@Table(name = "adverse_event")
public class AdverseEvent {
    @Id
    @Column(name = "adverse_event_accession")
    private String accession;
    private String causality;
    private String description;
    @Column(name = "end_study_day")
    private Float endStudyDay;
    @Column(name = "end_time")
    private String endTime;
    @Column(name = "name_preferred")
    private String namePreferred;
    @Column(name = "name_reported")
    private String nameReported;
    @Column(name = "outcome_preferred")
    private String outcomePreferred;
    @Column(name = "outcome_reported")
    private String outcomeReported;
    @Column(name = "relation_to_nonstudy_treatment")
    private String relationToNonStudyTreatment;
    @Column(name = "relation_to_study_treatment")
    private String relationToStudyTreatment;
    @Column(name = "severity_preferred")
    private String severityPreferred;
    @Column(name = "severity_reported")
    private String severityReported;
    @Column(name = "start_study_day")
    private Float startStudyDay;
    @Column(name = "start_time")
    private String startTime;
    @ManyToOne
    @JoinColumn(name = "study_accession")
    private Study study;
    @ManyToOne
    @JoinColumn(name = "subject_accession")
    private Subject subject;
    
    public AdverseEvent() {
    }

    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public String getCausality() {
        return causality;
    }

    public void setCausality(String causality) {
        this.causality = causality;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public Float getEndStudyDay() {
        return endStudyDay;
    }

    public void setEndStudyDay(Float endStudyDay) {
        this.endStudyDay = endStudyDay;
    }

    public String getEndTime() {
        return endTime;
    }

    public void setEndTime(String endTime) {
        this.endTime = endTime;
    }

    public String getNamePreferred() {
        return namePreferred;
    }

    public void setNamePreferred(String namePreferred) {
        this.namePreferred = namePreferred;
    }

    public String getNameReported() {
        return nameReported;
    }

    public void setNameReported(String nameReported) {
        this.nameReported = nameReported;
    }

    public String getOutcomePreferred() {
        return outcomePreferred;
    }

    public void setOutcomePreferred(String outcomePreferred) {
        this.outcomePreferred = outcomePreferred;
    }

    public String getOutcomeReported() {
        return outcomeReported;
    }

    public void setOutcomeReported(String outcomeReported) {
        this.outcomeReported = outcomeReported;
    }

    public String getRelationToNonStudyTreatment() {
        return relationToNonStudyTreatment;
    }

    public void setRelationToNonStudyTreatment(String relationToNonStudyTreatment) {
        this.relationToNonStudyTreatment = relationToNonStudyTreatment;
    }

    public String getRelationToStudyTreatment() {
        return relationToStudyTreatment;
    }

    public void setRelationToStudyTreatment(String relationToStudyTreatment) {
        this.relationToStudyTreatment = relationToStudyTreatment;
    }

    public String getSeverityPreferred() {
        return severityPreferred;
    }

    public void setSeverityPreferred(String severityPreferred) {
        this.severityPreferred = severityPreferred;
    }

    public String getSeverityReported() {
        return severityReported;
    }

    public void setSeverityReported(String severityReported) {
        this.severityReported = severityReported;
    }

    public Float getStartStudyDay() {
        return startStudyDay;
    }

    public void setStartStudyDay(Float startStudyDay) {
        this.startStudyDay = startStudyDay;
    }

    public String getStartTime() {
        return startTime;
    }

    public void setStartTime(String startTime) {
        this.startTime = startTime;
    }

    public Study getStudy() {
        return study;
    }

    public void setStudy(Study study) {
        this.study = study;
    }

    public Subject getSubject() {
        return subject;
    }

    public void setSubject(Subject subject) {
        this.subject = subject;
    }

}
