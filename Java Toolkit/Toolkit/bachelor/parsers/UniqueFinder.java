package parsers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

import biological_units.Gene;

public class UniqueFinder {
	public static Set<String> IDs = new HashSet<>();

	public UniqueFinder() {
	}

	public void parseReference(String ref) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		String line;
		while ((line = br.readLine()) != null) {
			// String[] ids = line.split(" ");
			// System.out.println(ids.length);
			// for(int i = 0; i < ids.length; i++)
			if (!IDs.contains(line))
				IDs.add(line);
		}
		br.close();
	}

	public void parseReference2(String ref, String out) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
		BufferedWriter bw = new BufferedWriter(new FileWriter(out));
		String line;
		while ((line = br.readLine()) != null) {
			if (!IDs.contains(line)) {
				IDs.add(line);
				bw.write(line + "\r\n");
			}
		}
		bw.close();
		br.close();
	}

	public static void main(String[] args) throws IOException {
		UniqueFinder uf = new UniqueFinder();
		try {
			Scanner scanner = new Scanner(System.in);
			System.out.println("Enter file path...");
			String dataset = scanner.next();
			uf.parseReference(dataset);
		} catch (FileNotFoundException e) {
			uf.parseReference2(args[0], args[1]);
		}
		System.out.println(IDs.size());
	}
}
