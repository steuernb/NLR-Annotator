package nlr_parser;


import java.io.File;
import java.io.IOException;
import java.util.Enumeration;
import java.util.Hashtable;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import support.BioSequence;
import support.FastaReader;

public class MastFile {

		
	
	File mastFile;
	
	Hashtable<String, MastMotifHitList > nblrrs;
	
	double pvalue_threshold;
	
	
	
	
	
	public MastFile(File mastFile)throws IOException, ParserConfigurationException, SAXException{
		this.mastFile         = mastFile;
	
		this.pvalue_threshold = 1e-5;      //default value
		nblrrs = new Hashtable<String,MastMotifHitList>();
		readEntries(pvalue_threshold);
		
	}
	
	public MastFile(File mastFile,  double pvalue_threshold)throws IOException, ParserConfigurationException, SAXException{
		
		this.pvalue_threshold = pvalue_threshold;
		this.mastFile         = mastFile;
		
		nblrrs = new Hashtable<String,MastMotifHitList>();
		readEntries(pvalue_threshold);
		
		
		
	}
	
	private void readEntries(double pvalue)throws IOException, ParserConfigurationException, SAXException{
		Hashtable<String, MastMotifHitList> mastMotifHitLists = this.parseXML();
		
		for(Enumeration<String> myenum = mastMotifHitLists.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			MastMotifHitList list = mastMotifHitLists.get(key);
			MastMotifHitList fw = new MastMotifHitList();
			MastMotifHitList rv = new MastMotifHitList();
			
			for(Enumeration<MastMotifHit> myenum2 = list.getMotifs().elements(); myenum2.hasMoreElements();){
				MastMotifHit hit = myenum2.nextElement();
				if( hit.getPvalue() > pvalue){
					continue;
				}
				if( hit.forwardStrand){
					fw.addMotif(hit);
				}else{
					rv.addMotif(hit);
				}
				
				if( fw.getSize()>0){
					fw.sort();
					nblrrs.put(key+"+", fw);
				}
				if( rv.getSize()>0){
					rv.sort();
					nblrrs.put(key+"-", rv);
				}
			}
			
			
		}
		
	}
	
	
	
	/**
	 * called by readEntries() and reads in the xml file produced by MAST. New versions 9.1 and 11.2 are supported.
	 * 
	 * @return
	 * @throws IOException
	 * @throws SAXException
	 * @throws ParserConfigurationException
	 */
	private Hashtable<String, MastMotifHitList> parseXML()throws IOException, SAXException, ParserConfigurationException{
		
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		Document dom = db.parse(this.mastFile);
		Hashtable<String,MastMotifHitList> motifLists = new Hashtable<String,MastMotifHitList>();
		
		Element rootElement = dom.getDocumentElement();
		System.out.println(rootElement.getNodeName());
		
		NodeList nodeList = rootElement.getElementsByTagName("sequences");
		for( int i = 0; i< nodeList.getLength(); i++){
			Element elementSequences = (Element) nodeList.item(i);
			
			
			NodeList nodeListSequences = elementSequences.getElementsByTagName("sequence");
			for( int j = 0; j< nodeListSequences.getLength(); j++){
				Element elementSequence = (Element)nodeListSequences.item(j);
				String sequence_id = elementSequence.getAttribute("name");
				int length = Integer.parseInt(elementSequence.getAttribute("length"));
				
				String sequenceName = sequence_id.split("_frame")[0];
				MastMotifHitList list = new MastMotifHitList();
				if(!motifLists.containsKey(sequenceName)){
					motifLists.put(sequenceName, list);
				}
				
				
				char[] sequence = new char[length];
				for( int index = 0; index< sequence.length; index++){
					sequence[index]= '-';
				}
				
				
				
				
				NodeList nodeListSeg = elementSequence.getElementsByTagName("seg");
				for( int k = 0; k< nodeListSeg.getLength(); k++){
					Element elementSeg = (Element)nodeListSeg.item(k);
					int start = Integer.parseInt(elementSeg.getAttribute("start"));
					char[] s = elementSeg.getAttribute("data").toCharArray();
					for(int index = 0; index < s.length; index ++){
						sequence[start-1+index] = s[index];
					}
					
					NodeList nodeListHit = elementSeg.getElementsByTagName("hit");
					for( int l=0; l< nodeListHit.getLength(); l++){
						Element elementHit = (Element)nodeListHit.item(l);
						
						int pos = Integer.parseInt(elementHit.getAttribute("pos"));
						int motif = 0;
						if( elementHit.hasAttribute("motif")){
							motif = Integer.parseInt(elementHit.getAttribute("motif").split("_")[1]);
						}
						if( elementHit.hasAttribute("idx")){
							motif = Integer.parseInt(elementHit.getAttribute("idx")) + 1;
						}
						
						double pvalue = Double.parseDouble(elementHit.getAttribute("pvalue"));
						String match = elementHit.getAttribute("match");
						
						boolean forwardStrand = true;
						int frame = 0;
								
						try{
							String string= sequence_id.split("_frame")[1];
							if(string.substring(0, 1).equalsIgnoreCase("-")){
								forwardStrand = false;
							}
							frame = Integer.parseInt(string.substring(1,2));
						}catch(ArrayIndexOutOfBoundsException e){}
						
						MastMotifHit hit = new MastMotifHit(sequenceName, pos,  motif, pvalue, match, forwardStrand, frame);
						motifLists.get(sequenceName).addMotif(hit);
					}
					
					
				}
				
				
				
				
				
				motifLists.get(sequenceName).addAASequence(sequence.toString());
				
				
				
			}
			
			
			
		}
		
		return motifLists;
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
