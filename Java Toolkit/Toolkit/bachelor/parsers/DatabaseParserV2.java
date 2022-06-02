package parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;

import org.apache.commons.lang3.StringUtils;

import java.util.Scanner;
import java.util.Set;

import biological_units.Protein;
import biological_units.ProteinComplex;
import utils.Utils;

public class DatabaseParserV2 {

	static HashMap<String, ProteinComplex> complexes = new HashMap<>();
	static HashMap<String, List<String>> scores = new HashMap<String, List<String>>();
	static HashMap<String, Protein> proteins = new HashMap<>();
	static HashMap<String, String> card2Ens = new HashMap<>();
	static List<String> genes4Mapping = new ArrayList<String>();
	static int countProts = 0;

	public DatabaseParserV2() {
	}

	public void parseComplexes(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
		br.readLine();
		while ((line = br.readLine()) != null) {
			String id = line.split("\t")[0];
			int size = Integer.parseInt(line.split("\t")[1]);
			String name = line.split("\t")[3];
			String[] memb = line.split("\t")[8].split(" ");
			List<String> members = new ArrayList<>();
			for (int i = 0; i < memb.length; i++) {
				Protein prot = new Protein(memb[i]);
				members.add(memb[i]);
				if (!proteins.containsKey(memb[i]))
					proteins.put(memb[i], prot);
				proteins.get(memb[i]).setComplexID(id);
			}
			complexes.put(id, new ProteinComplex(id, name, size, members));
		}
		br.close();
	}

	public void compareToExternal(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
		br.readLine();
		while ((line = br.readLine()) != null) {
			// System.out.println(line);
			if (!genes4Mapping.contains(line))
				genes4Mapping.add(line);
		}
		br.close();
		for (Entry<String, ProteinComplex> entry : complexes.entrySet()) {
			scores.put(entry.getKey(), new ArrayList<>());
			for (String s : entry.getValue().getMembers()) {
				if (genes4Mapping.contains(s)) {
					countProts++;
					scores.get(entry.getKey()).add(s);
				}
			}
		}

	}

	public void parseCardIDs(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
//		br.readLine();
		while ((line = br.readLine()) != null) {
			String card = line.split("\t")[0];
			String ens = line.split("\t")[1];
			// System.out.println(line);
			if (proteins.containsKey(ens)) {
				proteins.get(ens).setCardID(card);
				card2Ens.put(card, ens);
			}
		}
		br.close();
	}
	
	public void parseProteomicsReprog(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
//		br.readLine();
		while ((line = br.readLine()) != null) {
			String card = line.split("\t")[2];
			if(!card2Ens.containsKey(card))
				continue;
			String s1 = line.split("\t")[4];
			String s2 = line.split("\t")[4];
			String s3 = line.split("\t")[4];
			String s4 = line.split("\t")[4];
			String s5 = line.split("\t")[4];
			if(s1.contains("NaN") || s2.contains("NaN") || s3.contains("NaN") || s4.contains("NaN") || s5.contains("NaN"))
				continue;
			Double d0d3 = Double.parseDouble(line.split("\t")[4]);
			Double d3d6 = Double.parseDouble(line.split("\t")[5]);
			Double d6d9 = Double.parseDouble(line.split("\t")[6]);
			Double d9d12 = Double.parseDouble(line.split("\t")[7]);
			Double d12ips = Double.parseDouble(line.split("\t")[8]);
			proteins.get(card2Ens.get(card)).setReprogrammingValues("D0D3", d0d3);
			proteins.get(card2Ens.get(card)).setReprogrammingValues("D3D6", d3d6);
			proteins.get(card2Ens.get(card)).setReprogrammingValues("D6D9", d6d9);
			proteins.get(card2Ens.get(card)).setReprogrammingValues("D9D12", d9d12);
			proteins.get(card2Ens.get(card)).setReprogrammingValues("D12iPS", d12ips);
		}
		br.close();
	}
	public void parseProteomicsCancer(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
//		br.readLine();
		while ((line = br.readLine()) != null) {
			String card = line.split("\t")[2];
			if(!card2Ens.containsKey(card))
				continue;
			Double d0d3 = Double.parseDouble(line.split("\t")[4]);
			Double d3d6 = Double.parseDouble(line.split("\t")[5]);
			Double d6d9 = Double.parseDouble(line.split("\t")[6]);
			Double d9d12 = Double.parseDouble(line.split("\t")[7]);
			Double d12ips = Double.parseDouble(line.split("\t")[8]);
			proteins.get(card2Ens.get(card)).setReprogrammingValues("D0D3", d0d3);
			proteins.get(card2Ens.get(card)).setReprogrammingValues("D3D6", d3d6);
			proteins.get(card2Ens.get(card)).setReprogrammingValues("D6D9", d6d9);
			proteins.get(card2Ens.get(card)).setReprogrammingValues("D9D12", d9d12);
			proteins.get(card2Ens.get(card)).setReprogrammingValues("D12iPS", d12ips);
		}
		br.close();
	}

	public static void main(String[] args) throws IOException {
		int counter = 0;
		int counter1 = 0;
		DatabaseParserV2 dp = new DatabaseParserV2();
		dp.parseComplexes(args[0]);
		Scanner scanner = new Scanner(System.in);
		System.out.println("Choose dataset from: 11 Cancer Cell Lines dataset(1); Reprogramming dataset(2)...");
		String dataset = scanner.next();
		if (dataset.equals("1"))
			dp.compareToExternal(args[1]);
		else if (dataset.equals("2"))
			dp.compareToExternal(args[2]);
		for (Entry<String, ProteinComplex> entry : complexes.entrySet()) {
			// System.out.println(entry.getValue().getSize() + "-" +
			// scores.get(entry.getKey()).size());
			// if (entry.getValue().getSize() !=
			// entry.getValue().getMembers().size())
			// counter++;
			// if (entry.getValue().getSize() ==
			// entry.getValue().getMembers().size())
			// System.out.println(entry.getValue().getId() + " " +
			// entry.getValue().getSize() + " "
			// + entry.getValue().getName() + " " +
			// entry.getValue().getMembers());
			// if (entry.getValue().getSize() !=
			// scores.get(entry.getKey()).size())
			// counter++;
			if (scores.get(entry.getKey()).size() > 4) {
				counter++;
				counter1 += scores.get(entry.getKey()).size();
			}
			// if (entry.getValue().getSize() ==
			// scores.get(entry.getKey()).size())
			// System.out.println("yeah");
		}
		System.out.println("Protein complexes found in chosen dataset: " + counter);
		System.out.println("Proteins found in chosen dataset: " + counter1);
		System.out.println(genes4Mapping.size());

	}
}
