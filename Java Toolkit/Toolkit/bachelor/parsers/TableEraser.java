package parsers;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class TableEraser {
	public TableEraser(){
		
	}
	public void parseFile(String ref, String output) {
		try {
			try {
				BufferedReader br = new BufferedReader(new FileReader(new File(ref)));
				BufferedWriter bw = new BufferedWriter(new FileWriter(output));

				String line;
				while ((line = br.readLine()) != null) {
					if(!line.startsWith("IP"))
						continue;
					
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

}
