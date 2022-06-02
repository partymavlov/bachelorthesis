package parsers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import biological_units.Gene;

public class TableOrganizer {
	public TableOrganizer() {

	}

	public void parseFile1(String ref, String output) {
		try {
			try {
				BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
				BufferedWriter bw = new BufferedWriter(new FileWriter(output));

				String line;
				while ((line = br.readLine()) != null) {
					if (!line.isEmpty() && line != null) {
						String[] ids = line.split(";");
						for (int i = 0; i < ids.length; i++){
							if(ids[i].contains("/"))
								ids[i] = ids[i].split("-")[0];
							bw.write(ids[i] + "\r\n");
					
						}}
				}
				bw.close();
				br.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} catch (

		IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	public void parseFile2(String ref, String output) {
		try {
			try {
				BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
				BufferedWriter bw = new BufferedWriter(new FileWriter(output));

				String line;
				while ((line = br.readLine()) != null) {
					if (!line.isEmpty() && line != null) {
						String[] ids = line.split(";");
						for (int i = 0; i < ids.length; i++){
							if(ids[i].contains("-"))
								ids[i] = ids[i].split("-")[0];
							bw.write(ids[i] + "\r\n");
					
						}}
				}
				bw.close();
				br.close();
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} catch (

		IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public static void main(String[] args) throws IOException {
		TableOrganizer to = new TableOrganizer();
		to.parseFile1(args[0], args[1]);
	}
}
