package org.reactome.immport.ws.model.analysis;

import com.fasterxml.jackson.annotation.JsonProperty;

/**
 * 
 * @author brunsont
 *
 */
public class BiosampleAnalysis {
	
	@JsonProperty("gene_name")
	private String geneName;
	@JsonProperty("logFC")
	private Double logFC;
	@JsonProperty("AveExpr")
	private Double averageExpression;
	@JsonProperty("t")
	private Double t;
	@JsonProperty("P.Value")
	private Double pVal;
	@JsonProperty("adj.P.Val")
	private Double adjustedPVal;
	@JsonProperty("B")
	private Double B;
	
	public BiosampleAnalysis() {
		super();
	}

	public BiosampleAnalysis(String geneName, Double logFC, Double averageExpression, Double t, Double pVal, Double adjustedPVal,
			Double b) {
		super();
		this.geneName = geneName;
		this.logFC = logFC;
		this.averageExpression = averageExpression;
		this.t = t;
		this.pVal = pVal;
		this.adjustedPVal = adjustedPVal;
		this.B = b;
	}

	public String getGeneName() {
		return geneName;
	}

	public void setGeneName(String geneName) {
		this.geneName = geneName;
	}

	public Double getLogFC() {
		return logFC;
	}

	public void setLogFC(Double logFC) {
		this.logFC = logFC;
	}

	public Double getAverageExpression() {
		return averageExpression;
	}

	public void setAverageExpression(Double averageExpression) {
		this.averageExpression = averageExpression;
	}

	public Double getT() {
		return t;
	}

	public void setT(Double t) {
		this.t = t;
	}

	public Double getpVal() {
		return pVal;
	}

	public void setpVal(Double pVal) {
		this.pVal = pVal;
	}

	public Double getAdjustedPVal() {
		return adjustedPVal;
	}

	public void setAdjustedPVal(Double adjustedPVal) {
		this.adjustedPVal = adjustedPVal;
	}

	public Double getB() {
		return B;
	}

	public void setB(Double b) {
		B = b;
	}
}
