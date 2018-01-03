package nlr_parser;


import java.io.File;
import java.io.IOException;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import support.BioSequence;
import support.FastaReader;

public class MastFile {

		
	
	File mastFile;
	private Hashtable<String, MastEntry> entries;
	Hashtable<String, MastMotifHitList > nblrrs;
	
	double pvalue_threshold;
	
	
	
	
	
	public MastFile(File mastFile)throws IOException{
		this.mastFile         = mastFile;
		this.entries          = new Hashtable<String, MastEntry>();
		this.pvalue_threshold = 1e-5;      //default value
		this.readEntries();
	}
	
	public MastFile(File mastFile,  double pvalue_threshold)throws IOException{
		
		this.pvalue_threshold = pvalue_threshold;
		this.mastFile         = mastFile;
		this.entries          = new Hashtable<String, MastEntry>();
		this.readEntries();
	}
	
	private void readEntries()throws IOException{
		MastXMLReader reader = new MastXMLReader(this.mastFile);
		MastEntry entry = reader.readEntry();
		while(entry != null){
			entries.put(entry.getName(), entry);
			entry = reader.readEntry();
		}
		reader.close();
	}
	
	
	public void addAASequences(File sequenceFile)throws IOException{
		
		FastaReader fastaReader = new FastaReader(sequenceFile);
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			if( nblrrs.containsKey(seq.getIdentifier() )){
				nblrrs.get(seq.getIdentifier() ).addAASequence(seq.getSequence());
			}
			
		}
		fastaReader.close();
		
	}
	
	
	public void addDnaSequences(File sequenceFile)throws IOException{
		
		FastaReader fastaReader = new FastaReader(sequenceFile);
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			
			if( nblrrs.containsKey(seq.getIdentifier() + "+" )){
				nblrrs.get(seq.getIdentifier() + "+").addDnaSequence(seq.getSequence());
			}
			
			if(nblrrs.containsKey(seq.getIdentifier() +"-")){
				nblrrs.get(seq.getIdentifier() + "-").addDnaSequence(seq.getReverseComplementarySequence());
			}
			
		}
		fastaReader.close();
		
	}
	
	
	public Hashtable<String,MastMotifHitList> getNibblers(){
		return this.nblrrs;
	}
	
	
	/**
	 *
	 * 
	 * @throws IOException
	 */
	public void mergeDNAEntries()throws IOException{
		nblrrs = new Hashtable<String,MastMotifHitList>();
		if(entries.size()==0){
			this.readEntries();
		}
		
		for(Enumeration<String> myenum = entries.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			String name = key.split("_frame")[0];
			String strand = key .split("_frame")[1].substring(0,1);
			
			Vector<MastMotifHit> v = entries.get(key).getHits();
			MastMotifHitList list = new MastMotifHitList();
			list.setSequenceLength(entries.get(key).getLength());
			
			if(nblrrs.containsKey(name+strand)){
				list= nblrrs.get(name+strand);
			}
			for(Enumeration<MastMotifHit> myenum2 = v.elements(); myenum2.hasMoreElements();){
				MastMotifHit mmh = myenum2.nextElement();
				if(mmh.getPvalue()<=this.pvalue_threshold){
					list.addMotif(mmh);
				}
			}
			if(list.getSize()>0){
				nblrrs.put(name+strand, list);
			}
			
		}
		
	}
	
	
		
	public void prepareProteinEntries()throws IOException {
		nblrrs = new Hashtable<String,MastMotifHitList>();
		if(entries.size()==0){
			this.readEntries();
		}
		
		for(Enumeration<String> myenum = entries.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			String name = key;
			
			
			Vector<MastMotifHit> v = entries.get(key).getHits();
			MastMotifHitList list = new MastMotifHitList();
			list.setSequenceLength(entries.get(key).getLength());
			
			if(nblrrs.containsKey(name)){
				list= nblrrs.get(name);
			}
			for(Enumeration<MastMotifHit> myenum2 = v.elements(); myenum2.hasMoreElements();){
				MastMotifHit mmh = myenum2.nextElement();
				if(mmh.getPvalue()<=this.pvalue_threshold){
					list.addMotif(mmh);
				}
			}
			if(list.getSize()>0){
				nblrrs.put(name, list);
			}
			
		}
	}

	
	
	
	
	
	/*
	
	public void writeOutputGFF(File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		out.write("##gff-version 2");
		out.newLine();
		out.write("##source-version MastParser " + NLRParser.getVersion());
		out.newLine();
		
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm");
		Date date = new Date();
		out.write("##date "+dateFormat.format(date)); 
		out.newLine();
		out.write("##Type DNA");
		out.newLine();
		out.write("#seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute");
		out.newLine();
		
		
		for(Enumeration<String> myenum1 = this.nblrrs.keys(); myenum1.hasMoreElements();){
			String sequenceName = myenum1.nextElement();
			MastMotifHitList list = this.nblrrs.get(sequenceName);
		
			for( Enumeration<MastMotifHit> myenum2 = list.getMotifs().elements(); myenum2.hasMoreElements();){
				MastMotifHit hit = myenum2.nextElement();
				out.write(sequenceName +"\tMastParser\tMastMotif\t"  );
				
				int frame = Integer.parseInt(hit.getSequence_name().split("_frame")[1].substring(1,2) );
				
				boolean isForwardStrand = true;
				if(hit.getSequence_name().split("_frame")[1].substring(0, 1).equalsIgnoreCase("-")){
					isForwardStrand= false;
				}
				int aa_start = hit.getPos();
				int aa_end = hit.getPos() + hit.getMatch().length();
				int nucl_start = -1;
				int nucl_end = -1;
				try{
					nucl_start = MastFile.aminoAcidPositionToNucleotidePosition(aa_start, frame,isForwardStrand, list.getSequenceLength(), true);
					nucl_end   = MastFile.aminoAcidPositionToNucleotidePosition(aa_end, frame,isForwardStrand, list.getSequenceLength(), false);
					if( !isForwardStrand){
						nucl_start = MastFile.aminoAcidPositionToNucleotidePosition(aa_end, frame,isForwardStrand, list.getSequenceLength(), false);
						nucl_end   = MastFile.aminoAcidPositionToNucleotidePosition(aa_start, frame,isForwardStrand, list.getSequenceLength(), true);
					}
				}catch (OutOfFrameException e ){
					System.out.println(hit.getSequence_name() +" does not encode a proper frame");
				}
				String strand = "+";
				if(!isForwardStrand){
					strand = "-";
				}
				
				
				out.write(nucl_start + "\t" + nucl_end + "\t" +  hit.getPvalue() + "\t" +strand +"\t" + (Math.abs(frame) ) + "\t" + "name "+hit.getMotif() );
				out.newLine();
				
				
				
			}
			
		}	
		
		
		
		
		
		
		
		out.close();
		
		
		
		
	}
	*/
	/**
	 * 
	 * translates an aminoacid position back to the corresponding nucleotide position. It can be defined if the output position should be the first or the last base of the triplet
	 * 
	 * 
	 * 
	 * @param aaPos
	 * 			amino acid position
	 * @param frame
	 * 			the frame of the sequence (i.e. -3, -2, -1, 1, 2 or 3)
	 * @param sequenceLength
	 * 			the length of the amino acid sequence. In the reverse frames the position is dependent on the length.
	 * @param start
	 * 			true will give the position of the first base in the triplet, false the third.
	 * @return
	 * 			the corresponding nucleotide position.
	 * @throws OutOfFrameException
	 * 			thrown if the denoted reading frame is bullshit. 
	 */
	public static int aminoAcidPositionToNucleotidePosition(int aaPos, int frame,boolean isForwardStrand, int sequenceLength, boolean start)throws OutOfFrameException{
		
		
		
		if(  frame < -2 || frame >2){
			throw new OutOfFrameException("You gave this method a frame of " + frame + ". We are talking about reading frames on nucleotide sequences here.");
		}
	
		
		if( start){
		
			if(isForwardStrand){
				return ((aaPos-1) * 3 + frame + 1);
			
			}else{
			
				return (sequenceLength - ( (aaPos-1) *3) - (Math.abs(frame) )  );
			}
			
			
		}else{
			
			if(frame >0){
				return ((aaPos-1) * 3 + frame) +3;
			}
			
			if(frame <0){
				return (sequenceLength - ( (aaPos-1) *3) - (Math.abs(frame) ) -2  );
			}
			
			
		}
		return -1;
	}
	

	
	
	
	
}
