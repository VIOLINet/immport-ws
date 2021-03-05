package org.reactome.immport.ws.model.analysis;

import com.fasterxml.jackson.annotation.JsonProperty;

/**
 * 
 * @author brunsont
 *
 */
public class BiosampleAnalysis {
	
	private String geneName;
	private Double logFC;
	private Double averageExpression;
	private Double t;
	private Double pVal;
	private Double adjustedPVal;
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

	@JsonProperty("gene_name")
	public String getGeneName() {
		return geneName;
	}

	@JsonProperty("gene_name")
	public void setGeneName(String geneName) {
		this.geneName = geneName;
	}

	@JsonProperty("logFC")
	public Double getLogFC() {
		return logFC;
	}

	@JsonProperty("logFC")
	public void setLogFC(Double logFC) {
		this.logFC = logFC;
	}

	@JsonProperty("AveExpr")
	public Double getAverageExpression() {
		return averageExpression;
	}

	@JsonProperty("AveExpr")
	public void setAverageExpression(Double averageExpression) {
		this.averageExpression = averageExpression;
	}

	@JsonProperty("t")
	public Double getT() {
		return t;
	}

	@JsonProperty("t")
	public void setT(Double t) {
		this.t = t;
	}

	@JsonProperty("pValue")
	public Double getpVal() {
		return pVal;
	}

	@JsonProperty("P.Value")
	public void setpVal(Double pVal) {
		this.pVal = pVal;
	}

	@JsonProperty("adjPValue")
	public Double getAdjustedPVal() {
		return adjustedPVal;
	}

	@JsonProperty("adj.P.Val")
	public void setAdjustedPVal(Double adjustedPVal) {
		this.adjustedPVal = adjustedPVal;
	}

	@JsonProperty("B")
	public Double getB() {
		return B;
	}

	@JsonProperty("B")
	public void setB(Double b) {
		B = b;
	}
}
