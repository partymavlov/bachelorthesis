package parsers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import javax.swing.CellEditor;

import biological_units.Gene;

public class TableConverter {
	String reference;
	String cancFile;
	String reprogFile;
	static HashMap<String, String> gene2ProteinNames = new HashMap<>();
	static List<Gene> genes = new ArrayList<Gene>();
	static List<String> uniIDs = new ArrayList<>();
	static Set<String> check = new HashSet<>();
	static int cellLinesCounter;

	public TableConverter() {
	}

	public void parseReference(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
		while (br.readLine() != null) {
			line = br.readLine();
			if (!line.split("\t")[0].startsWith("ENSG"))
				continue;
			String geneID = line.split("\t")[0];
			String uniID = line.split("\t")[2];
			String protName = line.split("\t")[5];
			gene2ProteinNames.put(protName, geneID);
			Gene g = new Gene(geneID);
			g.setUniProtID(uniID);
			genes.add(g);
			uniIDs.add(uniID);
		}
		br.close();
	}

	public void parseProcessedTarget(String tar) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(tar)));
		String line;
		while ((line = br.readLine()) != null) {
			// System.out.println(line);
			if (!check.contains(line)) {
				check.add(line);
				if (uniIDs.contains(line))
					cellLinesCounter++;
			}
		}
		br.close();
	}
	public void parseUnprocessedTarget(String tar) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(tar)));
		String line;
		while ((line = br.readLine()) != null) {
			// System.out.println(line);
			if (!check.contains(line)) {
				check.add(line);
				if (uniIDs.contains(line))
					cellLinesCounter++;
			}
		}
		br.close();
	}
	public static void main(String[] args) throws IOException {
		TableConverter tc = new TableConverter();
		tc.parseReference(args[0]);
		tc.parseProcessedTarget(args[1]);
		// for (Entry<String, String> entry : gene2ProteinNames.entrySet())
		// System.out.println(entry.getKey() + " " + entry.getValue());
		// for (Gene g : genes)
		// System.out.println(g.getEnsmblID() + " " + g.getUniProtID());
//		for (String s : uniIDs)
//			System.out.println(s);
		System.out.println(cellLinesCounter);
	}
}
