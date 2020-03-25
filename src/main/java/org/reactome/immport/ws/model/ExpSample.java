package org.reactome.immport.ws.model;

import java.util.Set;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.JoinTable;
import javax.persistence.ManyToMany;
import javax.persistence.ManyToOne;
import javax.persistence.OneToMany;
import javax.persistence.Table;

@Entity
@Table(name = "expsample")
public class ExpSample {
    @Id
    @Column(name = "expsample_accession")
    private String accession;
    private String description;
    private String name;
    @ManyToOne
    @JoinColumn(name = "result_schema", referencedColumnName = "name")
    private ResultSchema resultSchema;
    @ManyToOne
    @JoinColumn(name = "experiment_accession")
    private Experiment experiment;
    @OneToMany
    @JoinColumn(name = "expsample_accession")
    private Set<SampleGeneExpression> geneExpressions;
    @ManyToMany
    @JoinTable(name = "expsample_2_biosample",
               joinColumns = {@JoinColumn(name = "expsample_accession")},
               inverseJoinColumns = {@JoinColumn(name = "biosample_accession")})
    private Set<BioSample> bioSamples;
    
    public ExpSample() {
    }

    public Set<BioSample> getBioSamples() {
        return bioSamples;
    }

    public void setBioSamples(Set<BioSample> bioSamples) {
        this.bioSamples = bioSamples;
    }

    public Set<SampleGeneExpression> getGeneExpressions() {
        return geneExpressions;
    }

    public void setGeneExpressions(Set<SampleGeneExpression> geneExpressions) {
        this.geneExpressions = geneExpressions;
    }

    public Experiment getExperiment() {
        return experiment;
    }

    public void setExperiment(Experiment experiment) {
        this.experiment = experiment;
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

    public ResultSchema getResultSchema() {
        return resultSchema;
    }

    public void setResultSchema(ResultSchema resultSchema) {
        this.resultSchema = resultSchema;
    }
    
}
