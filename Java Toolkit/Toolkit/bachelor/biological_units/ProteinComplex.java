package biological_units;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ProteinComplex {
	String id;
	String name;
	int size;
	List<String> members = new ArrayList<>();

	public ProteinComplex(String id, String name, int size, List<String> members) {
		this.id = id;
		this.name = name;
		this.size = size;
		this.members = members;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getSize() {
		return size;
	}

	public void setSize(int size) {
		this.size = size;
	}

	public List<String> getMembers() {
		return members;
	}

	public void setMembers(List<String> members) {
		this.members = members;
	}
}
