package support;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class FastqReader implements BioSequenceReader{

	BufferedReader in;
	String lastLine;
	int qualityOffset;
		
		
	public FastqReader(File inputFile)throws IOException{
		
		this.initialize(inputFile);
		this.qualityOffset = 33;
		
	}
		
	public FastqReader(File inputFile, int qualityOffset)throws IOException{
		this.initialize(inputFile);
		this.qualityOffset = qualityOffset;
	}	

	
	
	
	/**
	 * 
	 * Reads the next entry in the file. It is assumed that each sequence consists of four lines and the third line starts with a +
	 * 
	 * @return
	 * 			The Biosequence object.
	 * @throws IOException
	 */
	public BioSequence readEntry()throws IOException{
		if (lastLine == null){
			return null;
		}
		
		// read identifier and description (if exists)
		String identifier = lastLine.split("\\s")[0].substring(1);
		String description = "";
		if( lastLine.split("\\s").length > 1){
			description = lastLine.split("\\s", 2)[1];
		}
		
		
		StringBuilder builder = new StringBuilder();
		
		//read sequence
		lastLine = in.readLine();
		if( lastLine != null){
			builder.append(lastLine);
		}else{
			return null;
		}
		
		//check if third line starts with a +
		lastLine = in.readLine();
		if(lastLine == null || !lastLine.trim().startsWith("+")){
			return null;
		}
		
		
		BioSequence seq = new BioSequence(identifier, builder);
		if( description.length() > 0){
			seq.setDescription(description);
		}
		
		
		//set the quality
		lastLine = in.readLine();
		if( lastLine != null){
			seq.setQuality(lastLine, qualityOffset);
		}else{
			return null;
		}
		
		lastLine = in.readLine();
		while(lastLine != null && !lastLine.trim().startsWith("@")){
			lastLine = in.readLine();
		}
		
		return seq;
	}

	
	/**
	 * Close the FastqReader
	 * 
	 * @throws IOException
	 */
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
				
				while(lastLine != null && !lastLine.trim().startsWith("@")){
					this.lastLine=in.readLine();
				}
				
	}
	
	
	
	
	



	
	
	
	
	
	
}
