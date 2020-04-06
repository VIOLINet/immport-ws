package org.reactome.immport.ws.model.queries;

public class VOToGSM {

	String gsm;
	String race;
	
	public VOToGSM(String gsm, String race) {
		this.gsm = gsm;
		this.race = race;
	}

	public String getGsm() {
		return gsm;
	}

	public void setGsm(String gsm) {
		this.gsm = gsm;
	}

	public String getRace() {
		return race;
	}

	public void setRace(String race) {
		this.race = race;
	}
	
	
}
