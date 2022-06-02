package biological_units;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class Protein {
//	String name;
	String ensID;
	String cardID;
	HashMap<String, Double> reprogrammingValues = new HashMap<>();
	HashMap<String, Double> cancerValues = new HashMap<>();
	List<String> complexIDs = new ArrayList<>();
	
	public Protein(String ensID){
		this.ensID=ensID;
	}
//	public String getName() {
//		return name;
//	}
//	public void setName(String name) {
//		this.name = name;
//	}
	public String getEnsID() {
		return ensID;
	}
	public void setEnsID(String ensID) {
		this.ensID = ensID;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((ensID == null) ? 0 : ensID.hashCode());
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Protein other = (Protein) obj;
		if (ensID == null) {
			if (other.ensID != null)
				return false;
		} else if (!ensID.equals(other.ensID))
			return false;
		return true;
	}
	public String getCardID() {
		return cardID;
	}
	public void setCardID(String cardID) {
		this.cardID = cardID;
	}
	public HashMap<String, Double> getReprogrammingValues() {
		return reprogrammingValues;
	}
	public void setReprogrammingValues(String condition, Double values) {
		reprogrammingValues.put(condition, values);
	}
	public HashMap<String, Double> getCancerValues() {
		return cancerValues;
	}
	public void setCancerValues(String condition, Double values) {
		cancerValues.put(condition, values);
	}
	public List<String> getComplexIDs() {
		return complexIDs;
	}
	public void setComplexID(String complexID) {
		complexIDs.add(complexID);
	}
	public boolean isInCancerDataset(){
		if(cancerValues.isEmpty())
			return false;
		else
			return true;
	}
	public boolean isInReprogDataset(){
		if(reprogrammingValues.isEmpty())
			return false;
		else
			return true;
	}
}
