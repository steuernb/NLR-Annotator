package nlr_parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class MastXMLReader {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			MastXMLReader reader = new MastXMLReader(new File("/Users/steuernb/Documents/projects/MastParser/mast.xml"));
			MastEntry entry = reader.readEntry();
			while(entry != null){
				System.out.println(entry.getName() + "\t" + entry.getHits().size());
				entry = reader.readEntry();
			}
			reader.close();
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	
	
	
	
	File mastXMLFile;
	BufferedReader in;
	private String inputline;
	
	
	
	
	
	public MastXMLReader(File inputFile)throws IOException{
		this.mastXMLFile = inputFile;
		this.in = new BufferedReader(new FileReader(inputFile));
		
		this.inputline = in.readLine();
		
		/*
		 * parse out the used algorithm
		 */
		
		String def = "";
		while(inputline!= null && !inputline.contains("<sequence ")){
			def = def + inputline.trim();
			
			inputline = in.readLine();
			
		}
		try{
			def = def + inputline.split("<sequence ",-1)[0];
		}catch (NullPointerException e ) {}  //do nothing in case there was no nibbler in the sequence file
		
		
		
	}
	
	
	
	
	public File getFile(){
		return this.mastXMLFile;
	}
	
	
	
	
	
	
	
	
	
	public MastEntry readEntry()throws IOException{
		
		
		
		
		/*
		 * Initialize the Mast entry
		 */
		
		String entryString = getEntryString();
		
		
		if(entryString == null){
			return null;
		}
		
		
		MastEntry mastEntry= new MastEntry(entryString);
		
		
		
		
		
		return mastEntry;
	}
	
	
	
	private String getEntryString()throws IOException{
		StringBuilder builder = new StringBuilder();
		
		
		
		/*
		 * move to start of iteration
		 */
		while(inputline !=null && !inputline.contains("<sequence")){
			inputline = in.readLine();
		}
		if(inputline == null){
			return null;
		}
		
		/*
		 * read the entire iteration
		 */
		
		
		
		builder.append("<sequence" +  inputline.split("<sequence", -1)[1]);
		inputline = in.readLine();
		
		while (inputline != null && !inputline.contains("</sequence>")){
			builder.append(inputline.trim());
			inputline = in.readLine();
		}
		
		
		
		builder.append(inputline.split("</sequence>", -1)[0].trim());		
		
		
		return builder.toString();
	}
	
	
	
	
	
	public void close()throws IOException{
		in.close();
	}
	
	
	
	
	
	
	
	
	
	
	
	
}
