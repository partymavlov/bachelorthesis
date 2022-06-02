package biological_units;

import java.util.ArrayList;
import java.util.List;

public class Gene {
	String ensmblID;
	String uniProtID;
	List<String> proteinNames = new ArrayList<String>();

	public Gene(String ensmblID) {
		this.ensmblID = ensmblID;
	}

	public String getEnsmblID() {
		return ensmblID;
	}

	public String getUniProtID() {
		return uniProtID;
	}

	public List<String> getProteins() {
		return proteinNames;
	}

	public void setUniProtID(String id) {
		this.uniProtID = id;
	}

	public void setProteins(List<String> prots) {
		this.proteinNames = prots;
	}

	public void addProtein(String prot) {
		this.proteinNames.add(prot);
	}
}
