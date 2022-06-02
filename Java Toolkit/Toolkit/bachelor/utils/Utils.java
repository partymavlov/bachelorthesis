package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import org.apache.commons.lang3.StringUtils;

import biological_units.ProteinComplex;

public class Utils {
	static int counter0;
	static int counter1;
	static int counter2;
	static int counter3;
	static int counter4;
	static int counter5;
	static int counter6;
	static int counter7;
	static int counter8;
	static int counter9;
	static int counter10;

	public Utils() {

	}

	public static List<String> convertEnsToUniGeneIDs(List<String> ensGenes) throws IOException {
		List<String> output = new ArrayList<>();
		HashMap<String, String> idPairs = new HashMap<>();
		BufferedReader br = new BufferedReader(new FileReader(
				new File("C:/Users/Martin/Desktop/Lernen/6. Semester/Bachelor/Uniprot_Genes_reviewed+unreviewed.tab")));
		String line;
		while ((line = br.readLine()) != null) {
			if (!line.startsWith("ENS"))
				continue;
			String ens = line.split("\t")[0];
			String uni = line.split("\t")[2];
			idPairs.put(ens, uni);
		}
		br.close();
		for (String s : ensGenes) {
			if (idPairs.containsKey(s))
				output.add(idPairs.get(s));
			else {
				output.add(s);
				System.out.println("oops");
			}
		}
		return output;
	}

	public static List<String> convertUniToEnsGeneIDs(List<String> uniGenes) throws IOException {
		List<String> output = new ArrayList<>();
		HashMap<String, String> idPairs = new HashMap<>();
		BufferedReader br = new BufferedReader(new FileReader(
				new File("C:/Users/Martin/Desktop/Lernen/6. Semester/Bachelor/Uniprot_Genes_reviewed+unreviewed.tab")));
		String line;
		while ((line = br.readLine()) != null) {
			if (!line.startsWith("ENS"))
				continue;
			String ens = line.split("\t")[0];
			String uni = line.split("\t")[2];
			idPairs.put(ens, uni);
		}
		br.close();
		for (String s : uniGenes) {
			if (idPairs.containsKey(s))
				output.add(idPairs.get(s));
			else {
				output.add(s);
				System.out.println("oops");
			}
		}
		return output;
	}

	public void filterNaNs1(String ref, String out) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String line;
		br.readLine();
		while ((line = br.readLine()) != null) {
			// System.out.println(line.split("\t")[3]);
			String ids = null;
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
															ids = line.split("\t")[2];
			if (ids != null)
				for (String s : ids.split(";"))
					bw.write(s + "\r\n");

		}
		bw.close();
		br.close();

	}

	public double median(List<Double> numbers) {
		if (numbers.size() == 3) {
			return Math.max(Math.min(numbers.get(0),numbers.get(1)), Math.min(Math.max(numbers.get(0),numbers.get(1)),numbers.get(2)));
		} else if (numbers.size() == 2)
			return (numbers.get(0) + numbers.get(1)) / 2;
		else
			return 0.0;
	}

	public void filterNaNs2(String ref, String out) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		// BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String line;
		br.readLine();
		while ((line = br.readLine()) != null) {
			// System.out.println(line.split("\t")[3]);
			String ids = null;
			String c01 = line.split("\t")[16];
			String c02 = line.split("\t")[17];
			String c11 = line.split("\t")[18];
			String c12 = line.split("\t")[19];
			String c21 = line.split("\t")[20];
			String c22 = line.split("\t")[21];
			String c31 = line.split("\t")[22];
			String c32 = line.split("\t")[23];
			String c41 = line.split("\t")[24];
			String c42 = line.split("\t")[25];
			// if ((!c01.contains("NaN") || !c02.contains("NaN")))
			// if ((!c11.contains("NaN") || !c12.contains("NaN")))
			// if ((!c21.contains("NaN") || !c22.contains("NaN")))
			// if ((!c31.contains("NaN") || !c32.contains("NaN")))
			// if ((!c41.contains("NaN") || !c42.contains("NaN")))
			// if(!c01.contains("NaN") && !c02.contains("NaN") &&
			// !c11.contains("NaN") && !c12.contains("NaN") &&
			// !c21.contains("NaN") && !c22.contains("NaN") &&
			// !c31.contains("NaN") && !c32.contains("NaN") &&
			// !c41.contains("NaN") && !c42.contains("NaN"))
			// if(StringUtils.countMatches(line, "NaN") == 0 &&
			// line.split("\t")[2].length()>1)
			// bw.write(line.split("\t")[2] + "\r\n");
			if (StringUtils.countMatches(line, "NaN") == 0)
				counter0++;
			if (StringUtils.countMatches(line, "NaN") == 1)
				counter1++;
			if (StringUtils.countMatches(line, "NaN") == 2)
				counter2++;
			if (StringUtils.countMatches(line, "NaN") == 3)
				counter3++;
			if (StringUtils.countMatches(line, "NaN") == 4)
				counter4++;
			if (StringUtils.countMatches(line, "NaN") == 5)
				counter5++;
			if (StringUtils.countMatches(line, "NaN") == 6)
				counter6++;
			if (StringUtils.countMatches(line, "NaN") == 7)
				counter7++;
			if (StringUtils.countMatches(line, "NaN") == 8)
				counter8++;
			if (StringUtils.countMatches(line, "NaN") == 9)
				counter9++;
			if (StringUtils.countMatches(line, "NaN") == 10)
				counter10++;
		}
		// bw.close();
		br.close();

	}

	public static void main(String[] args) throws IOException {
		Utils u = new Utils();
		u.filterNaNs2(args[0], args[1]);
		int total = counter0 + counter1 + counter2 + counter3 + counter4 + counter5 + counter6 + counter7 + counter8
				+ counter9 + counter10;
		System.out.println("0 zeroes: " + counter0);
		System.out.println("1 zeroes: " + counter1);
		System.out.println("2 zeroes: " + counter2);
		System.out.println("3 zeroes: " + counter3);
		System.out.println("4 zeroes: " + counter4);
		System.out.println("5 zeroes: " + counter5);
		System.out.println("6 zeroes: " + counter6);
		System.out.println("7 zeroes: " + counter7);
		System.out.println("8 zeroes: " + counter8);
		System.out.println("9 zeroes: " + counter9);
		System.out.println("10 zeroes: " + counter10);
		System.out.println("Total: " + total);
	}
}
