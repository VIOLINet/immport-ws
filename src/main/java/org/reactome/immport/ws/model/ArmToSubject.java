package org.reactome.immport.ws.model;

import java.util.Objects;

import javax.persistence.Column;
import javax.persistence.EmbeddedId;
import javax.persistence.Entity;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.MapsId;
import javax.persistence.Table;

/**
 * This class is created so that we can find the ages for subjects. The implementation of this class is based on 
 * this web page: https://vladmihalcea.com/the-best-way-to-map-a-many-to-many-association-with-extra-columns-when-using-jpa-and-hibernate/.
 * @author wug
 *
 */
@Entity
@Table(name = "arm_2_subject")
public class ArmToSubject {
    
    @EmbeddedId
    private ArmToSubjectId id;
    @ManyToOne
    @MapsId("armAccession")
    @JoinColumn(name = "arm_accession") // This is important to point out what columns are used together with MapsId.
    private ArmOrCohort arm;
    @ManyToOne
    @MapsId("subjectAccession")
    @JoinColumn(name = "subject_accession")
    private Subject subject;
    @Column(name = "age_event")
    private String ageEvent;
    @Column(name = "age_event_specify")
    private String ageEventSpecify;
    @Column(name = "age_unit")
    private String ageUnit;
    @Column(name = "max_subject_age")
    private Float maxSubjectAge;
    @Column(name = "min_subject_age")
    private Float minSubjectAge;
    @Column(name = "subject_phenotype")
    private String subjectPhenotype;
    
    private ArmToSubject() {
    }
    
    public ArmToSubject(ArmOrCohort arm, Subject subject) {
        this.arm = arm;
        this.subject = subject;
        this.id = new ArmToSubjectId(arm.getAccession(), subject.getAccession());
    }

    public String getAgeEventSpecify() {
        return ageEventSpecify;
    }

    public void setAgeEventSpecify(String ageEventSpecify) {
        this.ageEventSpecify = ageEventSpecify;
    }

    public String getAgeUnit() {
        return ageUnit;
    }

    public void setAgeUnit(String ageUnit) {
        this.ageUnit = ageUnit;
    }

    public Float getMaxSubjectAge() {
        return maxSubjectAge;
    }

    public void setMaxSubjectAge(Float maxSubjectAge) {
        this.maxSubjectAge = maxSubjectAge;
    }

    public Float getMinSubjectAge() {
        return minSubjectAge;
    }

    public void setMinSubjectAge(Float minSubjectAge) {
        this.minSubjectAge = minSubjectAge;
    }

    public String getSubjectPhenotype() {
        return subjectPhenotype;
    }

    public void setSubjectPhenotype(String subjectPhenotype) {
        this.subjectPhenotype = subjectPhenotype;
    }

    public ArmToSubjectId getId() {
        return id;
    }

    public void setId(ArmToSubjectId id) {
        this.id = id;
    }

    public ArmOrCohort getArm() {
        return arm;
    }

    public void setArm(ArmOrCohort arm) {
        this.arm = arm;
    }

    public Subject getSubject() {
        return subject;
    }

    public void setSubject(Subject subject) {
        this.subject = subject;
    }

    public String getAgeEvent() {
        return ageEvent;
    }

    public void setAgeEvent(String ageEvent) {
        this.ageEvent = ageEvent;
    }
    
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || o.getClass() != getClass())
            return false;
        ArmToSubject that = (ArmToSubject) o;
        return Objects.equals(arm, that.arm) && Objects.equals(subject, that.subject);
    }
    
    @Override
    public int hashCode() {
        return Objects.hash(arm, subject);
    }

}
