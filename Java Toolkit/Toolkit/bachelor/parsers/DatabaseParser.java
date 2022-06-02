package parsers;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.Set;

import biological_units.ProteinComplex;
import utils.Utils;

public class DatabaseParser {
	static HashMap<String, ProteinComplex> complexes = new HashMap<>();
	static HashMap<String, List<String>> scores = new HashMap<String, List<String>>();
	static List<String> genes4Mapping = new ArrayList<String>();
	static int countProts = 0;

	public DatabaseParser() {
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
			for (int i = 0; i < memb.length; i++)
				members.add(memb[i]);
			complexes.put(id, new ProteinComplex(id, name, size, Utils.convertEnsToUniGeneIDs(members)));
		}
		br.close();
	}

	public void compareToExternal(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
		br.readLine();
		while ((line = br.readLine()) != null) {
			// System.out.println(line);
			genes4Mapping.add(line);
		}
		br.close();
		for (Entry<String, ProteinComplex> entry : complexes.entrySet()) {
			scores.put(entry.getKey(), new ArrayList<>());
			for (String s : entry.getValue().getMembers()) {
				if (genes4Mapping.contains(s)) {
					countProts++;
					scores.get(entry.getKey()).add("1");
				}
			}
		}

	}
	
	public void compareToExternal1(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
		br.readLine();
		while ((line = br.readLine()) != null) {
			// System.out.println(line);
			genes4Mapping.add(line);
		}
		br.close();
		for (Entry<String, ProteinComplex> entry : complexes.entrySet()) {
			scores.put(entry.getKey(), new ArrayList<>());
			for (String s : entry.getValue().getMembers()) {
				if (genes4Mapping.contains(s)) {
					countProts++;
					scores.get(entry.getKey()).add("1");
				}
			}
		}

	}

	public void compareToExternalSpec(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
		br.readLine();
		while ((line = br.readLine()) != null) {
			// System.out.println(line.split("\t")[3]);
			String c01 = line.split("\t")[10];
			String c02 = line.split("\t")[11];
			String c03 = line.split("\t")[12];
			String c11 = line.split("\t")[13];
			String c12 = line.split("\t")[14];
			String c13 = line.split("\t")[15];
			String c21 = line.split("\t")[16];
			String c22 = line.split("\t")[17];
			String c23 = line.split("\t")[18];
			String c31 = line.split("\t")[19];
			String c32 = line.split("\t")[20];
			String c33 = line.split("\t")[21];
			String c41 = line.split("\t")[22];
			String c42 = line.split("\t")[23];
			String c43 = line.split("\t")[24];
			String c51 = line.split("\t")[25];
			String c52 = line.split("\t")[26];
			String c53 = line.split("\t")[27];
			String c61 = line.split("\t")[28];
			String c62 = line.split("\t")[29];
			String c63 = line.split("\t")[30];
			String c71 = line.split("\t")[31];
			String c72 = line.split("\t")[32];
			String c73 = line.split("\t")[33];
			String c81 = line.split("\t")[34];
			String c82 = line.split("\t")[35];
			String c83 = line.split("\t")[36];
			String c91 = line.split("\t")[37];
			String c92 = line.split("\t")[38];
			String c93 = line.split("\t")[39];
			String c101 = line.split("\t")[40];
			String c102 = line.split("\t")[41];
			String c103 = line.split("\t")[42];
			if (!line.split("\t")[3].isEmpty() && line.split("\t")[3].length() > 5)
				if ((!c01.contains("NaN") && !c02.contains("NaN")) || (!c01.contains("NaN") && !c03.contains("NaN"))
						|| (!c03.contains("NaN") && !c02.contains("NaN")))
					if ((!c11.contains("NaN") && !c12.contains("NaN")) || (!c11.contains("NaN") && !c13.contains("NaN"))
							|| (!c13.contains("NaN") && !c12.contains("NaN")))
						if ((!c21.contains("NaN") && !c22.contains("NaN"))
								|| (!c21.contains("NaN") && !c23.contains("NaN"))
								|| (!c23.contains("NaN") && !c22.contains("NaN")))
							if ((!c31.contains("NaN") && !c32.contains("NaN"))
									|| (!c31.contains("NaN") && !c33.contains("NaN"))
									|| (!c33.contains("NaN") && !c32.contains("NaN")))
								if ((!c41.contains("NaN") && !c42.contains("NaN"))
										|| (!c41.contains("NaN") && !c43.contains("NaN"))
										|| (!c43.contains("NaN") && !c42.contains("NaN")))
									if ((!c51.contains("NaN") && !c52.contains("NaN"))
											|| (!c51.contains("NaN") && !c53.contains("NaN"))
											|| (!c53.contains("NaN") && !c52.contains("NaN")))
										if ((!c61.contains("NaN") && !c62.contains("NaN"))
												|| (!c61.contains("NaN") && !c63.contains("NaN"))
												|| (!c63.contains("NaN") && !c62.contains("NaN")))
											if ((!c71.contains("NaN") && !c72.contains("NaN"))
													|| (!c71.contains("NaN") && !c73.contains("NaN"))
													|| (!c73.contains("NaN") && !c72.contains("NaN")))
												if ((!c81.contains("NaN") && !c82.contains("NaN"))
														|| (!c81.contains("NaN") && !c83.contains("NaN"))
														|| (!c83.contains("NaN") && !c82.contains("NaN")))
													if ((!c91.contains("NaN") && !c92.contains("NaN"))
															|| (!c91.contains("NaN") && !c93.contains("NaN"))
															|| (!c93.contains("NaN") && !c92.contains("NaN")))
														if ((!c101.contains("NaN") && !c102.contains("NaN"))
																|| (!c101.contains("NaN") && !c103.contains("NaN"))
																|| (!c103.contains("NaN") && !c102.contains("NaN")))
															genes4Mapping.add(line.split("\t")[3]);

		}
		br.close();
		for (Entry<String, ProteinComplex> entry : complexes.entrySet()) {
			scores.put(entry.getKey(), new ArrayList<>());
			for (String s : entry.getValue().getMembers()) {
				for (String r : genes4Mapping) {
					if (r.contains(s)) {
						countProts++;
						scores.get(entry.getKey()).add("1");
					}
				}
			}
		}

	}

	public void compareToExternalSpec2(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
		br.readLine();
		while ((line = br.readLine()) != null) {
			genes4Mapping.add(line);
			

		}
		br.close();
		for (Entry<String, ProteinComplex> entry : complexes.entrySet()) {
			scores.put(entry.getKey(), new ArrayList<>());
			for (String s : entry.getValue().getMembers()) {
				for (String r : genes4Mapping) {
					if (r.contains(s)) {
						countProts++;
						scores.get(entry.getKey()).add("1");
					}
				}
			}
		}

	}

	public static void main(String[] args) throws IOException {
		int counter = 0;
		int counter1 = 0;
		DatabaseParser dp = new DatabaseParser();
		dp.parseComplexes(args[0]);
		Scanner scanner = new Scanner(System.in);
		System.out.println("Choose dataset from: 11 Cancer Cell Lines dataset(1); Reprogramming dataset(2)...");
		String dataset = scanner.next();
		if(dataset.equals("1")){
		System.out.print("Enter desired method (1/2)...");
		String method = scanner.next();

		if (dataset.equals("1") && method.equals("1"))
			dp.compareToExternal(args[1]);
		else if (dataset.equals("1") && method.equals("2"))
			dp.compareToExternalSpec(args[2]);
		else
			System.out.println("Please enter only available options");
		}
		else if (dataset.equals("2"))
			dp.compareToExternalSpec2(args[3]);
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
		System.out.println("Protein complexes found in chosen dataset: " + (279 - counter));
		System.out.println("Proteins found in chosen dataset: " + counter1);
		System.out.println(genes4Mapping.size());

	}
}
