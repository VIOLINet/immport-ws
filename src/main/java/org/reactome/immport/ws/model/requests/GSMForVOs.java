package org.reactome.immport.ws.model.requests;

import java.util.Collection;
import java.util.Set;

public class GSMForVOs {

	private Collection<String> voIds;
	private Collection<String> genderList;
	private Set<String> times;
	
	public GSMForVOs() {}

	public Collection<String> getVoIds() {
		return voIds;
	}

	public void setVoIds(Collection<String> voIds) {
		this.voIds = voIds;
	}

	public Collection<String> getGenderList() {
		return genderList;
	}

	public void setGenderList(Collection<String> genderList) {
		this.genderList = genderList;
	}

	public Set<String> getTimes() {
		return times;
	}

	public void setTimes(Set<String> times) {
		this.times = times;
	}
	
}
