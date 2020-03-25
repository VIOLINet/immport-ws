package org.reactome.immport.ws.model;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Table;

@Entity
@Table(name = "expsample_public_repository")
public class SampleGeneExpression {
    @Id
    @Column(name = "result_id")
    private Long id;
    @ManyToOne
    @JoinColumn(name = "expsample_accession")
    private ExpSample expSample;
    @Column(name = "repository_accession")
    private String repositoryAccession;
    @ManyToOne
    @JoinColumn(name = "repository_name")
    private PublicRepository repository;
    
    public SampleGeneExpression() {
        
    }

    public Long getId() {
        return id;
    }

    public void setId(Long id) {
        this.id = id;
    }

    public ExpSample getExpSample() {
        return expSample;
    }

    public void setExpSample(ExpSample expSample) {
        this.expSample = expSample;
    }

    public String getRepositoryAccession() {
        return repositoryAccession;
    }

    public void setRepositoryAccession(String repositoryAccession) {
        this.repositoryAccession = repositoryAccession;
    }

    public PublicRepository getRepository() {
        return repository;
    }

    public void setRepository(PublicRepository repository) {
        this.repository = repository;
    }
    
}
