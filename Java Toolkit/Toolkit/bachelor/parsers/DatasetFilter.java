package parsers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import biological_units.ProteinComplex;

public class DatasetFilter {
	List<String> ensIDs = new ArrayList<>();
	List<String> cardIDs = new ArrayList<>();

	public DatasetFilter() {

	}

	public void parseRef(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
		br.readLine();
		while ((line = br.readLine()) != null) {
			String[] memb = line.split("\t")[8].split(" ");
			for (int i = 0; i < memb.length; i++)
				ensIDs.add(memb[i]);
		}
		br.close();
	}

	public void parseMap(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
		while ((line = br.readLine()) != null) {
			String ens = line.split("\t")[1];
			String card = line.split("\t")[0];
			if (ensIDs.contains(ens))
				cardIDs.add(card);
		}
		br.close();
	}

	public void writeTable(String ref, String out) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String line;
		while ((line = br.readLine()) != null) {
			String card = line.split("\t")[2];
			for (String s : cardIDs) {
				if (card.contains(s)) {
					bw.write(line + "\r\n");
					break;
				}
			}
		}
		br.close();
	}
	public static void main(String[] args) throws IOException {
	DatasetFilter df = new DatasetFilter();
	df.parseRef(args[0]);
	df.parseMap(args[1]);
	df.writeTable(args[2], args[3]);
	}
}
