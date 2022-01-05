package annotator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import support.BioSequence;

public class SequenceChopper {

	 int fragmentLength; //the length of each sub-sequence
	 int overlap; //the overlap between subsequences
	
	
	
	BufferedReader bufferedReader;
	
	StringBuffer currentSequence;
	String currentID;
	int baseCharacter;
	int offset;
	
	
	
	
	public static void main(String[] args) {
		
		try {
			
		 
			File file = new File("/Users/steuernb/Documents/projects/Annotator2/test_data/Rgenes.CDS.fasta");
			SequenceChopper chopper = new SequenceChopper(file, 100,2);
			for(BioSequence seq = chopper.nextSequence(); seq != null; seq = chopper.nextSequence()) {
				System.out.println(chopper.baseCharacter);
				System.out.println(seq.getFastaString());
				
			}
			chopper.close();
		}catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
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
