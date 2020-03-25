package org.reactome.immport.ws.model;

import java.util.Set;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.JoinTable;
import javax.persistence.ManyToMany;
import javax.persistence.ManyToOne;
import javax.persistence.Table;

@Entity
@Table(name = "biosample")
public class BioSample {
    @Id
    @Column(name = "biosample_accession")
    private String accession;
    private String description;
    private String name;
    //TODO: Need to discuss if planned_visit is used in our project
    @ManyToOne
    @JoinColumn(name = "study_accession")
    private Study study;
    @Column(name = "study_time_collected")
    private Float studyTimeCollected;
    @Column(name = "study_time_collected_unit")
    private String studyTimeCollectedUnit;
    @Column(name = "study_time_t0_event")
    private String studyTimeT0Event;
    @Column(name = "study_time_t0_event_specify")
    private String studyTimeT0EventSpecify;
    @ManyToOne
    @JoinColumn(name = "subject_accession")
    private Subject subject;
    private String subtype;
    private String type; // Foreign key directly used
    // The relationship between BioSample and ExpSample is a real many-to-many relationship!
    @ManyToMany
    @JoinTable(name = "expsample_2_biosample",
               joinColumns = {@JoinColumn(name = "biosample_accession")},
               inverseJoinColumns = {@JoinColumn(name = "expsample_accession")})
    private Set<ExpSample> expSamples;
    
    public BioSample() {
        
    }

    public Set<ExpSample> getExpSamples() {
        return expSamples;
    }

    public void setExpSamples(Set<ExpSample> expSamples) {
        this.expSamples = expSamples;
    }

    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public Study getStudy() {
        return study;
    }

    public void setStudy(Study study) {
        this.study = study;
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

    public String getStudyTimeT0Event() {
        return studyTimeT0Event;
    }

    public void setStudyTimeT0Event(String studyTimeT0Event) {
        this.studyTimeT0Event = studyTimeT0Event;
    }

    public String getStudyTimeT0EventSpecify() {
        return studyTimeT0EventSpecify;
    }

    public void setStudyTimeT0EventSpecify(String studyTimeT0EventSpecify) {
        this.studyTimeT0EventSpecify = studyTimeT0EventSpecify;
    }

    public Subject getSubject() {
        return subject;
    }

    public void setSubject(Subject subject) {
        this.subject = subject;
    }

    public String getSubtype() {
        return subtype;
    }

    public void setSubtype(String subtype) {
        this.subtype = subtype;
    }

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }
}
