package org.reactome.immport.ws.model;

import java.util.Set;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.OneToMany;
import javax.persistence.Table;

import com.fasterxml.jackson.annotation.JsonManagedReference;

@Entity
@Table(name = "experiment")
public class Experiment {
    @Id
    @Column(name = "experiment_accession")
    private String accession;
    private String description;
    @Column(name = "measurement_technique")
    private String measureTechnique;
    private String name;
    // Link to study
    @ManyToOne
    @JoinColumn(name = "study_accession")
    @JsonManagedReference
    private Study study;
    @OneToMany
    @JoinColumn(name = "experiment_accession")
    private Set<ExpSample> expSamples;
    
    public Experiment() {
    }

    public Set<ExpSample> getExpSamples() {
        return expSamples;
    }

    public void setExpSamples(Set<ExpSample> expSamples) {
        this.expSamples = expSamples;
    }

    public Study getStudy() {
        return study;
    }

    public void setStudy(Study study) {
        this.study = study;
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

    public String getMeasureTechnique() {
        return measureTechnique;
    }

    public void setMeasureTechnique(String measureTechnique) {
        this.measureTechnique = measureTechnique;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

}
