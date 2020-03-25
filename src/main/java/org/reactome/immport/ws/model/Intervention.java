package org.reactome.immport.ws.model;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Table;

@Entity
@Table(name = "intervention")
public class Intervention {
    @Id
    @Column(name = "intervention_accession")
    private String accession;
    @Column(name = "compound_name_reported")
    private String compoundNameReported;
    @Column(name = "compound_role")
    private String compoundRole;
    private Float duration;
    @Column(name = "duration_unit")
    private String durationUnit;
    @Column(name = "end_day")
    private Float endDay;
    @Column(name = "end_time")
    private String endTime;
    @Column(name = "start_day")
    private Float startDay;
    @Column(name = "start_time")
    private String startTime;
    private String status;
    @ManyToOne
    @JoinColumn(name = "study_accession")
    private Study study;
    @ManyToOne
    @JoinColumn(name = "subject_accession")
    private Subject subject;

    public Intervention() {
        
    }

    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public String getCompoundNameReported() {
        return compoundNameReported;
    }

    public void setCompoundNameReported(String compoundNameReported) {
        this.compoundNameReported = compoundNameReported;
    }

    public String getCompoundRole() {
        return compoundRole;
    }

    public void setCompoundRole(String compoundRole) {
        this.compoundRole = compoundRole;
    }

    public Float getDuration() {
        return duration;
    }

    public void setDuration(Float duration) {
        this.duration = duration;
    }

    public String getDurationUnit() {
        return durationUnit;
    }

    public void setDurationUnit(String durationUnit) {
        this.durationUnit = durationUnit;
    }

    public Float getEndDay() {
        return endDay;
    }

    public void setEndDay(Float endDay) {
        this.endDay = endDay;
    }

    public String getEndTime() {
        return endTime;
    }

    public void setEndTime(String endTime) {
        this.endTime = endTime;
    }

    public Float getStartDay() {
        return startDay;
    }

    public void setStartDay(Float startDay) {
        this.startDay = startDay;
    }

    public String getStartTime() {
        return startTime;
    }

    public void setStartTime(String startTime) {
        this.startTime = startTime;
    }

    public String getStatus() {
        return status;
    }

    public void setStatus(String status) {
        this.status = status;
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
