package annotator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import support.BioSequence;


/**
 *   NLR-Annotator - find signatures of NLR-type resistance genes in genomic sequence
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
public class SequenceChopper {

	 int fragmentLength; //the length of each sub-sequence
	 int overlap; //the overlap between subsequences
	
	
	
	BufferedReader bufferedReader;
	
	StringBuffer currentSequence;
	String currentID;
	int baseCharacter;
	int offset;
	
	
	
	
	
	
	/**
	 * 
	 * Initialize the SequenceChopper with a fasta file as input. The chopper will create fragments of given sequences in the 
	 * fasta file according to given length and required overlap. 
	 * 
	 * 
	 * @param inputFastaFile
	 * @param fragmentLength
	 * @param overlap
	 * @throws IOException
	 */
	public SequenceChopper(File inputFastaFile, int fragmentLength, int overlap)throws IOException {
		this.fragmentLength = fragmentLength;
		this.overlap = overlap;
		
		initialize(inputFastaFile);
	}
	
	private void initialize(File inputFastaFile)throws IOException{
		
		
		
		FileInputStream fis = new FileInputStream(inputFastaFile);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		if(gzip){
			this.bufferedReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFastaFile))));
		}else{
			this.bufferedReader = new BufferedReader((new FileReader(inputFastaFile)));
		}		
		
		this.baseCharacter=bufferedReader.read();
		while((char)this.baseCharacter != '>' ){
			this.baseCharacter=bufferedReader.read();
		}
		
		this.readIdentifier();
		this.currentSequence = new StringBuffer();
		offset = 0;
		
	}
	
	
	
	/**
	 * 
	 * Get the next sequence fragment from the chopper. 
	 * The sequence will have the same identifier as the sequence it was taken from with an underscore '_' and the 0-based offset of the fragment.
	 * 
	 * @return
	 * 		The sequence fragment in form of a BioSequence
	 * @throws IOException
	 */
	public BioSequence nextSequence() throws IOException{
		
		if(this.baseCharacter <0) {//we are at end of file and don't want to read more
			return null;
		}
		
		
		while (this.baseCharacter >= 0) {
			char c = (char)this.baseCharacter;
			
			if(c == '>'){//hit next sequence
				BioSequence sequence = new BioSequence(this.currentID+"_"+offset, this.currentSequence.toString().toUpperCase());
				offset = 0;
				currentSequence = new StringBuffer();
				
				this.readIdentifier();
				return sequence;
			}else
			if( this.currentSequence.length() == fragmentLength) { //chop
				
				BioSequence sequence = new BioSequence(this.currentID+"_"+offset, this.currentSequence.toString().toUpperCase());
				offset = offset+fragmentLength - overlap;
				currentSequence = new StringBuffer(currentSequence.substring(fragmentLength - overlap));
				
				return sequence;
				
				
			}else {//all good just append next base
				if( c != '\n') {
					this.currentSequence.append(c);
				}
				
			}
			
			
			this.baseCharacter  = bufferedReader.read();
		}//end of file
		
		BioSequence sequence = new BioSequence(this.currentID+"_"+offset, this.currentSequence.toString().toUpperCase());
		
		
		return sequence;
		
		
	}
	
	
	
	
	private void readIdentifier() throws IOException{
		String header = "";
		while( (char) this.baseCharacter != '\n') {
			header = header + (char)this.baseCharacter;
			this.baseCharacter=bufferedReader.read();
		}
		this.currentID =  header.split("\\s")[0].substring(1);
	}
	
	
	public void close() throws IOException{
		this.bufferedReader.close();
	}
	
}
