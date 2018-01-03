package nlr_annotator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

import support.BioSequence;
import support.CLI;
import support.CLIParseException;


public class ChopSequence {

	
	
	public static void main(String[] args) {
		
		
		
		String helpString =		"Usage:\n-i <inputFasta>\t\t\tInput file in fasta format. Mandatory parameter\n"+
				
				"-o <output.fasta>\t\tOutput file in fasta format. Mandatory parameter\n"+
				"-l <int>\t\t\tlength of subsequences. Default 20000\n"+
				"-p <int>\t\t\toverlap between subsequences. Default 5000\n";

		
		
		
		
		
		
		CLI cli = new CLI();
		
		cli.parseOptions(args);
		
		int sequenceLength = 20000;
		int overlap = 5000;
		
		
		try {
			if( cli.hasOption("i")){
				File in = new File(cli.getArg("i"));
				
				if( cli.hasOption("o")){
					File out = new File(cli.getArg("o"));
					
					
					try{
						if(cli.hasOption("l")){
							sequenceLength = Integer.parseInt(cli.getArg("l"));
						}
					}
					catch(NumberFormatException e){
						System.err.println("argument for option l should be of type int. Proceeding with default value for sequence length: 10000");
						
					}
					
					try {
						if( cli.hasOption("p")){
							overlap = Integer.parseInt(cli.getArg("p"));
						}
					}	
					catch(NumberFormatException e){
						System.err.println("argument for option p should be of type int. Proceeding with default value for overlap: 3000");
					}
					
					chop(in,out, sequenceLength, overlap);
					
					
					
					
				}else{
					
					throw new CLIParseException("Option o not specified. o is mandatory");
				}
				
				
			}else{
				throw new CLIParseException("Option i not specified. i is mandatory");
			}
			
			
			
			
		} catch (IOException e) {
			e.printStackTrace();
			System.err.println(helpString);
		}
		catch(CLIParseException e){
			e.printStackTrace();
			System.err.println(helpString);
		}
		
		
		
		

	}

	
	
	
	public static void chop(File inputFile, File outputFile, int sequenceLength, int overlap)throws IOException{
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		BufferedReader in;
		
		FileInputStream fis = new FileInputStream(inputFile);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		if(gzip){
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile))));
		}else{
			in = new BufferedReader((new FileReader(inputFile)));
		}		
		
		
		
		int baseCharacter=in.read();
		while((char)baseCharacter != '>' ){
			baseCharacter=in.read();
		}
		baseCharacter=in.read();
		
		StringBuffer sequence = new StringBuffer("");
		boolean inHeader = true;
		long offset = 0;
		String identifier = "";
		String header = "";
		
		while(baseCharacter >=0){
			char c = (char)baseCharacter;
			if(inHeader){
				if(c == '\n'){
					inHeader = false;
					identifier = header.split("\\s")[0];
				}else{
					header = header + c;
				}
			}else{
				if(c == '>'){
					BioSequence seq = new BioSequence(identifier + "_" + offset, sequence.toString().toUpperCase());
					out.write(seq.getFastaString(100));
					inHeader = true;
					header = "";
					offset = 0;
					sequence = new StringBuffer("");
					
				}else{
					if(c != '\n'){
						sequence.append(c);
					}
					if(sequence.length() == sequenceLength){
						BioSequence seq = new BioSequence(identifier + "_" + offset, sequence.toString().toUpperCase());
						out.write(seq.getFastaString(100));
						
						sequence = new StringBuffer(sequence.substring(sequenceLength - overlap));
						offset = offset + sequenceLength - overlap;
					}
					
				}
				
			}
			
			
			
			baseCharacter = in.read();
		}
		
		BioSequence seq = new BioSequence(identifier + "_" + offset, sequence.toString().toUpperCase());
		out.write(seq.getFastaString(100));
		//System.out.println(seq.getIdentifier() + "\t" + seq.getLength());
		
		
		in.close();
		out.close();
		
	}
	
	
	
	
	

	
	
}
