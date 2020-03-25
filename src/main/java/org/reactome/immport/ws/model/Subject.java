package org.reactome.immport.ws.model;

import java.util.Set;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.OneToMany;
import javax.persistence.Table;

@Entity
@Table(name = "subject")
public class Subject {
    @Id
    @Column(name = "subject_accession")
    private String accession;
    @Column(name = "ancestral_population")
    // It is labeled as a foreign key to LK_ANCESTRAL_POPULATION. Apparently it is not enforced and used
    // as a free text in the current version of the database (March, 2020)
    private String ancestralPopulation;
    private String description;
    private String ethnicity; // Thought it is a foreign key, however, used it directly
    private String gender;
    private String race; // Foreign key used directly here
    @Column(name = "race_specify")
    private String specifiedRace;
    private String species; // Foreign key used directly here
    private String strain;
    @Column(name = "strain_characteristics")
    private String strainCharacteristics;
    @OneToMany
    @JoinColumn(name = "subject_accession")
    private Set<BioSample> biosamples;
    @OneToMany
    @JoinColumn(name = "subject_accession")
    private Set<AdverseEvent> adverseEvents;
    @OneToMany
    @JoinColumn(name = "subject_accession")
    private Set<Intervention> interventions;
    @OneToMany
    @JoinColumn(name = "subject_accession")
    private Set<ImmuneExposure> immuneExposures;
    
    public Subject() {
    }

    public Set<ImmuneExposure> getImmuneExposures() {
        return immuneExposures;
    }

    public void setImmuneExposures(Set<ImmuneExposure> immuneExposures) {
        this.immuneExposures = immuneExposures;
    }

    public Set<Intervention> getInterventions() {
        return interventions;
    }

    public void setInterventions(Set<Intervention> interventions) {
        this.interventions = interventions;
    }

    public Set<AdverseEvent> getAdverseEvents() {
        return adverseEvents;
    }

    public void setAdverseEvents(Set<AdverseEvent> adverseEvents) {
        this.adverseEvents = adverseEvents;
    }

    public Set<BioSample> getBiosamples() {
        return biosamples;
    }

    public void setBiosamples(Set<BioSample> biosamples) {
        this.biosamples = biosamples;
    }

    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public String getAncestralPopulation() {
        return ancestralPopulation;
    }

    public void setAncestralPopulation(String ancestralPopulation) {
        this.ancestralPopulation = ancestralPopulation;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getEthnicity() {
        return ethnicity;
    }

    public void setEthnicity(String ethnicity) {
        this.ethnicity = ethnicity;
    }

    public String getGender() {
        return gender;
    }

    public void setGender(String gender) {
        this.gender = gender;
    }

    public String getRace() {
        return race;
    }

    public void setRace(String race) {
        this.race = race;
    }

    public String getSpecifiedRace() {
        return specifiedRace;
    }

    public void setSpecifiedRace(String specifiedRace) {
        this.specifiedRace = specifiedRace;
    }

    public String getSpecies() {
        return species;
    }

    public void setSpecies(String species) {
        this.species = species;
    }

    public String getStrain() {
        return strain;
    }

    public void setStrain(String strain) {
        this.strain = strain;
    }

    public String getStrainCharacteristics() {
        return strainCharacteristics;
    }

    public void setStrainCharacteristics(String strainCharacteristics) {
        this.strainCharacteristics = strainCharacteristics;
    }

}
