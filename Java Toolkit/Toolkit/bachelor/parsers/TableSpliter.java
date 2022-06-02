package parsers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class TableSpliter {
	public TableSpliter() {
	}

	public void parseFile(String ref, String output) {
		try {
			try {
				BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
				BufferedWriter bw = new BufferedWriter(new FileWriter(output));

				String line;
				while ((line = br.readLine()) != null) {
					String[] genes = line.split(" ");
					for (int i = 0; i < genes.length; i++)
						bw.write(genes[i] + "\r\n");

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
	public static void main(String[] args) {
		TableSpliter ts = new TableSpliter();
		ts.parseFile(args[0], args[1]);
	}
}
