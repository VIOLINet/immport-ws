package org.reactome.immport.ws.model;

import java.util.Set;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.OneToMany;
import javax.persistence.Table;

import com.fasterxml.jackson.annotation.JsonManagedReference;

@Entity
@Table(name = "arm_or_cohort")
public class ArmOrCohort {
    @Id
    @Column(name = "arm_accession")
    private String accession;
    private String description;
    private String name;
    private String type;
    @ManyToOne
    @JoinColumn(name = "study_accession")
    @JsonManagedReference
    private Study study;
    @OneToMany(mappedBy = "arm",
               cascade = CascadeType.ALL)
    private Set<ArmToSubject> armToSubjects;
    @OneToMany
    @JoinColumn(name = "arm_accession")
    private Set<ImmuneExposure> exposures;
    
    public ArmOrCohort() {
    }

    public Set<ImmuneExposure> getExposures() {
        return exposures;
    }

    public void setExposures(Set<ImmuneExposure> exposures) {
        this.exposures = exposures;
    }

    public Set<ArmToSubject> getArmToSubjects() {
        return armToSubjects;
    }

    public void setArmToSubjects(Set<ArmToSubject> armToSubjects) {
        this.armToSubjects = armToSubjects;
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

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }

    public Study getStudy() {
        return study;
    }

    public void setStudy(Study study) {
        this.study = study;
    }
    
}
