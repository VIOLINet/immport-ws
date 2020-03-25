package org.reactome.immport.ws.model;

import java.io.Serializable;
import java.util.Objects;

import javax.persistence.Column;
import javax.persistence.Embeddable;

/**
 * The implementation of this class is based on:
 * https://vladmihalcea.com/the-best-way-to-map-a-many-to-many-association-with-extra-columns-when-using-jpa-and-hibernate/.
 * @author wug
 *
 */
@Embeddable
public class ArmToSubjectId implements Serializable {

    @Column(name = "arm_accession")
    private String armAccession;
    @Column(name = "subject_accession")
    private String subjectAccession;
    
    private ArmToSubjectId() {
    }
    
    public ArmToSubjectId(String armAccession, String subjectAccession) {
        this.armAccession = armAccession;
        this.subjectAccession = subjectAccession;
    }

    public String getArmAccession() {
        return armAccession;
    }

    public void setArmAccession(String armAccession) {
        this.armAccession = armAccession;
    }

    public String getSubjectAccession() {
        return subjectAccession;
    }

    public void setSubjectAccession(String subjectAccession) {
        this.subjectAccession = subjectAccession;
    }
    
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass())
            return false;
        ArmToSubjectId other = (ArmToSubjectId) o;
        return Objects.equals(armAccession, other.armAccession) &&
               Objects.equals(subjectAccession, other.subjectAccession);
    }
    
    @Override
    public int hashCode() {
        return Objects.hash(armAccession, subjectAccession);
    }
    
    @Override
    public String toString() {
        return armAccession + "+" + subjectAccession;
    }
    
}
