package support;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

/**
 *   FastaReader - Read a fasta file sequence by sequence
 *   
 *   Copyright (C) 2021  John Innes Centre
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 * 
 *  @author Burkhard Steuernagel
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
