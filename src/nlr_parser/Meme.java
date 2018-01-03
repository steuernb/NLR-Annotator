package nlr_parser;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class Meme {
	
	
	
	public void writeMeme(File directory)throws IOException{
		
		
		
		InputStream inputStream = getClass().getResourceAsStream("meme.xml");
		
		
		
		
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(directory, "meme.xml")));
		BufferedReader in = new BufferedReader(new InputStreamReader(inputStream));
		for(String inputline = in.readLine(); inputline != null; inputline = in.readLine()){
			out.write(inputline);
			out.newLine();
			
		}
		
		out.close();
				
		inputStream.close();		
				
	}
	
	
	
	
	
	
	
	
	
	
}	