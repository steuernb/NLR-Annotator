package support;

import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;



/**
 * The basic represenation of a sequence. This can be DNA or amino acid. Additionally to the sequence itself the per-base quality can also be stored. 
 * 
 * @author steuernb
 *
 */
public class BioSequence implements Comparable<BioSequence> {

	
	String identifier;
	StringBuilder sequence;
	String description;
	int[] quality;
	Hashtable<String,String> geneticCode;
	
	/**
	 * 
	 * Most simple invocation of a Biosequence. Just the sequence itself and an identifier.
	 * 
	 * @param identifier
	 * 			A String that denotes the name of the sequence. Ideally it should not contain any space or special characters. New line characters in the identifier will not be compatible with some methods.
	 * @param sequence
	 * 			The sequence itself.
	 */
	public BioSequence(String identifier, StringBuilder sequence) {
		
		this.identifier = identifier;
		this.sequence = new StringBuilder(sequence);
		
	}

	/**
	 * 
	 * 
	 * 
	 * @param identifier
	 * @param description
	 * @param sequence
	 */
	public BioSequence(String identifier,String description,  StringBuilder sequence) {
		
		this.identifier = identifier;
		this.sequence = new StringBuilder(sequence);
		this.description = description;
		
	}
	
	public BioSequence(String identifier, String sequence) {
		
		this.identifier = identifier;
		this.sequence = new StringBuilder(sequence);
		
	}
	
	
	   public int compareTo(BioSequence arg0) {
				BioSequence o2 =  arg0;
		    	int compared= 0;
		    	if(this.getLength()< o2.getLength()){
		    		compared = 1;
		    	}else
		    	if(this.getLength()> o2.getLength()){
		    		compared = -1;
		    	}else
		    	if(this.getLength()== o2.getLength()){
		    		compared = 0;
		    	}
		    	
		      return compared;		
		   }
	
	
	public String getIdentifier(){
		return this.identifier;
	}
	
	public String getDescription(){
		if(this.description!= null){
			return this.description;
		}else{
			return "";
		}
	}
	
	public String getSequence(){
		return this.sequence.toString();
	}
	
	public void appendSequence(String sequence){
		this.sequence.append(sequence);
	}
	
	public void setIdentifier(String identifier){
		this.identifier = identifier;
	}
	
	public void setDescription(String description){
		this.description = description;
	}
	
	public void setSequence(String sequence){
		this.sequence = new StringBuilder(sequence);
	}
	
	public void setSequence(StringBuilder sequence){
		this.sequence = sequence;
	}
	
	public void setQuality(int[] quality){
		this.quality = quality;
	}
	
	public void setQuality(String quality, int offset){
		this.quality = new int[quality.length()];
		char[] c = quality.toCharArray();
		for( int i = 0; i< this.quality.length; i++){
			this.quality[i] = c[i] - offset;
		}
	}
	
	public int[] getQuality (){
		return this.quality;
	}
	
	
	
	public String getQuality(int offset){
		String s = "";
		for( int i = 0; i< quality.length; i++){
			char c = Character.toChars( (quality[i] + offset)   )[0];
			s = s + c;
		}
		return s;
	}

	public String getFastaString(){
		String header = ">" + identifier;
		if( this.description != null && this.description.length()>0){
			header = header + " " + this.description;
		}
		
		return header + "\n" + this.sequence + "\n";
		
	}
	
	public String getFastaString(int linelength){
		String header = ">" + identifier;
		if( this.description != null && this.description.length()>0){
			header = header + " " + this.description;
		}
		
		
		 String formattedSequence = new String(this.sequence.substring(0, Math.min(this.sequence.length(),linelength) ));
	     String restSequence = new String(this.sequence.substring(Math.min(this.sequence.length(),linelength)));
	     int length = this.sequence.length() - linelength;

	    
	     
	     while (linelength <= length) {
	    	 formattedSequence = formattedSequence + "\n" + restSequence.substring(0, linelength);
	    	 restSequence = restSequence.substring(linelength, length);
	            length = length - linelength;
	        }
	        formattedSequence = formattedSequence + "\n" + restSequence;

	        if (!formattedSequence.endsWith("\n")){
	        	formattedSequence = formattedSequence+"\n";
	        }
		return header + "\n" + formattedSequence;
		
	}
	
	public int getLength(){
		return this.sequence.length();
	}
	
	
	/**
	 * 
	 * returns the fastq representation of the biosequence. The offset for the quality is Sanger (+33).
	 * 
	 * @return
	 * 		The fastq string of the sequence, ending with a new line. If quality had not been initialized, return value is null.
	 * 
	 */
	public String getFastqString(){
		if(this.lengthIsSynchronized()){
			return "@" + identifier + "\n" + sequence.toString() + "\n" + this.getQuality(33) + "\n";
		}else{
			return null;
		}
	}
	
	
	
	
	
	
	
	/**
	 * 
	 * Checks if the length of the sequence is consistent with the number of quality entries.
	 * 
	 * @return
	 * 		true if length is consistent, false if not or quality is not initialized.
	 * 
	 */
	public boolean lengthIsSynchronized(){
		
		if( this.quality == null){
			return false;
		}
		
		if( sequence.length() == quality.length){
			return true;
		}else{
			return false;
		}
	}
	
	
	/**
	 * calculates the number of bases in the sequence that are not an a,t,g or c. method is case insensitive.
	 * 
	 * @return
	 * 	number of non standard bases.
	 */
	  public int getNumberOfNonStandardBases(){
	    	int numNonStandardBases=0;
	    	Pattern p = Pattern.compile("[^ATGCatgc]+");
	    	Matcher m = p.matcher(this.getSequence());
	    	
	    	while(m.find()){
	    		numNonStandardBases = numNonStandardBases + m.group(0).length();
	    	}
	    	
	    	return numNonStandardBases;
	    }
	
	
	
	
	
	 /**
     * checks if a sequence is likely to be a nucleotide sequence. If more than half of the first 500 characters are A,T,G,C or N, than it is considered to be nucleotides
     * 
     * @return
     */
    public boolean isDNA(){
    	String s = this.sequence.substring(0, Math.min(500, this.sequence.length()));
    	Pattern p = Pattern.compile("([ATGCNatgcn]+)");
    	Matcher m = p.matcher(s);
    	int l = 0;
    	while(m.find()){
    		l = l + m.group(0).length();
    	}
    	if(l > s.length()/2){
    		return true;
    	}else{
    		return false;
    	}
    	
    }
	
	
    /**
     * translate an nucleotide sequence into an amino acid sequence. 
     * This results in 6 different reading frames: three on the forward strand and three on the reverse complementary strand.
     * 
     * 
     * @return An Array of sex FastaSequences
     */
    public BioSequence[] translate2Protein(){
    	BioSequence[] sixFrameTranslation = new BioSequence[6];
    	
    	StringBuilder seq = new StringBuilder(this.getSequence());
    	
    	int nucleotideLength = seq.length();
    	
    	StringBuilder seq1 = new StringBuilder(nucleotideLength/3);
    	StringBuilder seq2 = new StringBuilder(nucleotideLength/3);
    	StringBuilder seq3 = new StringBuilder(nucleotideLength/3);
    	StringBuilder seq4 = new StringBuilder(nucleotideLength/3);
    	StringBuilder seq5 = new StringBuilder(nucleotideLength/3);
    	StringBuilder seq6 = new StringBuilder(nucleotideLength/3);
    	
    	int terminus = seq.length()%3;
    	int end=0;
    	
    	if(terminus ==0){end=-2;}
    	if(terminus ==1){end=-3;}
    	if(terminus ==2){end=-4;}   	    	
    	for (int i = 0; i< seq.length()+end; i=i+3){
    		seq1.append(this.translateTriplet(seq.substring(i, i+3)));
    	}
    	sixFrameTranslation[0]= new BioSequence(this.getIdentifier()+"_frame+0", seq1);
    	//System.out.println("translated first");
    	if(terminus ==0){end=-4;}
    	if(terminus ==1){end=-2;}
    	if(terminus ==2){end=-3;}
    	for (int i = 1; i< seq.length()+end; i=i+3){
    		seq2.append(this.translateTriplet(seq.substring(i, i+3)));
    	}
    	sixFrameTranslation[1]= new BioSequence(this.getIdentifier()+"_frame+1", seq2);
    	
    	if(terminus ==0){end=-3;}
    	if(terminus ==1){end=-4;}
    	if(terminus ==2){end=-2;}
    	for (int i = 2; i< seq.length()+end; i=i+3){
    		seq3.append(this.translateTriplet(seq.substring(i, i+3)));
    	}
    	sixFrameTranslation[2]= new BioSequence(this.getIdentifier()+"_frame+2", seq3);
    	
    	
    	seq = new StringBuilder(this.getReverseComplementarySequence());
    	
    	if(terminus ==0){end=-2;}
    	if(terminus ==1){end=-3;}
    	if(terminus ==2){end=-4;}   	    	
    	for (int i = 0; i< seq.length()+end; i=i+3){
    		seq4.append(this.translateTriplet(seq.substring(i, i+3)));
    	}
    	sixFrameTranslation[3]= new BioSequence(this.getIdentifier()+"_frame-0", seq4);
    	
    	if(terminus ==0){end=-4;}
    	if(terminus ==1){end=-2;}
    	if(terminus ==2){end=-3;}
    	for (int i = 1; i< seq.length()+end; i=i+3){
    		seq5.append(this.translateTriplet(seq.substring(i, i+3)));
    	}
    	sixFrameTranslation[4]= new BioSequence(this.getIdentifier()+"_frame-1", seq5);
    	
    	if(terminus ==0){end=-3;}
    	if(terminus ==1){end=-4;}
    	if(terminus ==2){end=-2;}
    	for (int i = 2; i< seq.length()+end; i=i+3){
    		seq6.append(this.translateTriplet(seq.substring(i, i+3)));
    	}
    	sixFrameTranslation[5]= new BioSequence(this.getIdentifier()+"_frame-2", seq6);
    	
    	
    	
    	
    	return sixFrameTranslation;
    }

    
    public char translateTriplet(String s)throws StringIndexOutOfBoundsException{
    	
    	char aa = 'X';
    	if(s.length()!= 3){
    		throw new StringIndexOutOfBoundsException ("Length of the inputString has to be 3. Otherwise it is not a triplet and cannot be translated.");
    	}
    	
    	if(this.geneticCode == null){
    		this.loadGeneticCodeTable();
    	}
    	s = s.toUpperCase();
    	if( this.geneticCode.containsKey(s)){
    		return this.geneticCode.get(s).charAt(0);
    	}
    	
    	
    	return aa;
    }
    
    
  
    
  
    
   

    public String getReverseComplementarySequence(){
    	StringBuilder seq = new StringBuilder(getTranscribeSequence());
    	
    	return seq.reverse().toString();
    }
    
    /**
     * 
     * @return
     */
    public String getTranscribeSequence() {
        StringBuilder seq = new StringBuilder(getSequence());
        StringBuilder transcribe = new StringBuilder(seq.length());

        for (int i = 0; i < seq.length(); i++) {
            char t = seq.charAt(i);
            switch (t) {
            case 'A':
                transcribe.append('T');
                break;
            case 'T':
                transcribe.append('A');
                break;
            case 'G':
                transcribe.append('C');
                break;
            case 'C':
                transcribe.append('G');
                break;
            case 'a':
                transcribe.append('t');
                break;
            case 't':
                transcribe.append('a');
                break;
            case 'g':
                transcribe.append('c');
                break;
            case 'c':
                transcribe.append('g');
                break;    
            default:
            	transcribe.append(t);
                break;
            }
        }
        return transcribe.toString();
    }
    
    
    public void loadGeneticCodeTable(){
    	this.geneticCode = new Hashtable<String,String>();
    	
    	
    	this.geneticCode.put("AAA", "K");
    	this.geneticCode.put("AAT", "N");
    	this.geneticCode.put("AAG", "K");
    	this.geneticCode.put("AAC", "N");
    	this.geneticCode.put("ATA", "I");
    	this.geneticCode.put("ATT", "I");
    	this.geneticCode.put("ATG", "M");
    	this.geneticCode.put("ATC", "I");
    	this.geneticCode.put("AGA", "R");
    	this.geneticCode.put("AGT", "S");
    	this.geneticCode.put("AGG", "R");
    	this.geneticCode.put("AGC", "S");
    	this.geneticCode.put("ACA", "T");
    	this.geneticCode.put("ACT", "T");
    	this.geneticCode.put("ACG", "T");
    	this.geneticCode.put("ACC", "T");
    	this.geneticCode.put("TAA", "*");
    	this.geneticCode.put("TAT", "Y");
    	this.geneticCode.put("TAG", "*");
    	this.geneticCode.put("TAC", "Y");
    	this.geneticCode.put("TTA", "L");
    	this.geneticCode.put("TTT", "F");
    	this.geneticCode.put("TTG", "L");
    	this.geneticCode.put("TTC", "F");
    	this.geneticCode.put("TGA", "*");
    	this.geneticCode.put("TGT", "C");
    	this.geneticCode.put("TGG", "W");
    	this.geneticCode.put("TGC", "C");
    	this.geneticCode.put("TCA", "S");
    	this.geneticCode.put("TCT", "S");
    	this.geneticCode.put("TCG", "S");
    	this.geneticCode.put("TCC", "S");
    	this.geneticCode.put("GAA", "E");
    	this.geneticCode.put("GAT", "D");
    	this.geneticCode.put("GAG", "E");
    	this.geneticCode.put("GAC", "D");
    	this.geneticCode.put("GTA", "V");
    	this.geneticCode.put("GTT", "V");
    	this.geneticCode.put("GTG", "V");
    	this.geneticCode.put("GTC", "V");
    	this.geneticCode.put("GGA", "G");
    	this.geneticCode.put("GGT", "G");
    	this.geneticCode.put("GGG", "G");
    	this.geneticCode.put("GGC", "G");
    	this.geneticCode.put("GCA", "A");
    	this.geneticCode.put("GCT", "A");
    	this.geneticCode.put("GCG", "A");
    	this.geneticCode.put("GCC", "A");
    	this.geneticCode.put("CAA", "Q");
    	this.geneticCode.put("CAT", "H");
    	this.geneticCode.put("CAG", "Q");
    	this.geneticCode.put("CAC", "H");
    	this.geneticCode.put("CTA", "L");
    	this.geneticCode.put("CTT", "L");
    	this.geneticCode.put("CTG", "L");
    	this.geneticCode.put("CTC", "L");
    	this.geneticCode.put("CGA", "R");
    	this.geneticCode.put("CGT", "R");
    	this.geneticCode.put("CGG", "R");
    	this.geneticCode.put("CGC", "R");
    	this.geneticCode.put("CCA", "P");
    	this.geneticCode.put("CCT", "P");
    	this.geneticCode.put("CCG", "P");
    	this.geneticCode.put("CCC", "P");
    	
    	
    	
    	
    	
    	
    }

    public String toString(){
    	return this.getFastaString();
    }
    
    
  
}
