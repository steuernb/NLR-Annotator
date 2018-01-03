package support;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

/**
 * 
 * 
 * @author steuernb
 *
 */
public class FastaReader implements BioSequenceReader {


	BufferedReader in;
	String lastLine;
	
	
	public FastaReader(File inputFile)throws IOException{
		
		this.initialize(inputFile);
	}
	
	public FastaReader(BufferedReader in)throws IOException{
		
		this.in = in;
		
		this.lastLine = in.readLine();
		while(lastLine != null && !lastLine.trim().startsWith(">")){
			this.lastLine=in.readLine();
		}
				
	}
	
	
	public BioSequence readEntry()throws IOException{
		if (lastLine == null){
			return null;
		}
		
		
		String identifier = lastLine.split("\\s")[0].substring(1);
		String description = "";
		if( lastLine.split("\\s").length > 1){
			description = lastLine.split("\\s", 2)[1];
		}
		
		
		StringBuilder builder = new StringBuilder();
		lastLine = in.readLine();
		while(lastLine != null && !lastLine.trim().startsWith(">")){
			builder.append(lastLine.trim());
			lastLine = in.readLine();
		}
		
		BioSequence seq = new BioSequence(identifier, builder);
		if( description.length() > 0){
			seq.setDescription(description);
		}
		
		return seq;
	}

	
	
	public void close()throws IOException{
		this.in.close();
	}
	
	
	
	private void initialize(File inputFile)throws IOException{
		
		//check if the fastq is gzipped
			
		FileInputStream fis = new FileInputStream(inputFile);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();

		
		if(gzip){
			this.in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile))));
		}else{
			this.in = new BufferedReader((new FileReader(inputFile)));
		}
		
		
		
		this.lastLine = in.readLine();
		while(lastLine != null && !lastLine.trim().startsWith(">")){
			this.lastLine=in.readLine();
		}
				
	}
	
	
	
	
	
}
