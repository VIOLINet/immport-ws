package org.reactome.immport.ws.model;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Table;

@Entity
@Table(name = "immune_exposure")
public class ImmuneExposure {
    @Id
    @Column(name = "exposure_accession")
    private String accession;
    @ManyToOne
    @JoinColumn(name = "arm_accession")
    private ArmOrCohort armOrCohort;
    @Column(name = "disease_ontology_id")
    private String diseaseOntologyId;
    @Column(name = "disease_preferred")
    private String diseasePreferred;
    @Column(name = "disease_reported")
    private String diseaseReported;
    @Column(name = "disease_stage_preferred")
    private String diseaseStagePreferred;
    @Column(name = "disease_stage_reported")
    private String diseaseStageReported;
    @Column(name = "exposure_material_id")
    private String exposureMaterialId;
    @Column(name = "exposure_material_reported")
    private String exposureMaterialReported;
    @Column(name = "exposure_process_preferred")
    private String exposureProcessPreferred;
    @Column(name = "exposure_process_reported")
    private String exposureProcessReported;
    @ManyToOne
    @JoinColumn(name = "subject_accession")
    private Subject subject;
    
    public ImmuneExposure() {
        
    }

    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public ArmOrCohort getArmOrCohort() {
        return armOrCohort;
    }

    public void setArmOrCohort(ArmOrCohort armOrCohort) {
        this.armOrCohort = armOrCohort;
    }

    public String getDiseaseOntologyId() {
        return diseaseOntologyId;
    }

    public void setDiseaseOntologyId(String diseaseOntologyId) {
        this.diseaseOntologyId = diseaseOntologyId;
    }

    public String getDiseasePreferred() {
        return diseasePreferred;
    }

    public void setDiseasePreferred(String diseasePreferred) {
        this.diseasePreferred = diseasePreferred;
    }

    public String getDiseaseReported() {
        return diseaseReported;
    }

    public void setDiseaseReported(String diseaseReported) {
        this.diseaseReported = diseaseReported;
    }

    public String getDiseaseStagePreferred() {
        return diseaseStagePreferred;
    }

    public void setDiseaseStagePreferred(String diseaseStagePreferred) {
        this.diseaseStagePreferred = diseaseStagePreferred;
    }

    public String getDiseaseStageReported() {
        return diseaseStageReported;
    }

    public void setDiseaseStageReported(String diseaseStageReported) {
        this.diseaseStageReported = diseaseStageReported;
    }

    public String getExposureMaterialId() {
        return exposureMaterialId;
    }

    public void setExposureMaterialId(String exposureMaterialId) {
        this.exposureMaterialId = exposureMaterialId;
    }

    public String getExposureMaterialReported() {
        return exposureMaterialReported;
    }

    public void setExposureMaterialReported(String exposureMaterialReported) {
        this.exposureMaterialReported = exposureMaterialReported;
    }

    public String getExposureProcessPreferred() {
        return exposureProcessPreferred;
    }

    public void setExposureProcessPreferred(String exposureProcessPreferred) {
        this.exposureProcessPreferred = exposureProcessPreferred;
    }

    public String getExposureProcessReported() {
        return exposureProcessReported;
    }

    public void setExposureProcessReported(String exposureProcessReported) {
        this.exposureProcessReported = exposureProcessReported;
    }

    public Subject getSubject() {
        return subject;
    }

    public void setSubject(Subject subject) {
        this.subject = subject;
    }

}
