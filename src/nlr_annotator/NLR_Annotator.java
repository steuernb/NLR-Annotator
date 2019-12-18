package nlr_annotator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;
import java.util.zip.GZIPInputStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import nlr_parser.MastMotifHit;
import nlr_parser.MastMotifHitList;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import support.BioSequence;
import support.BlastHSP;
import support.BlastHit;
import support.BlastIPKReader;
import support.BlastIteration;
import support.BlastReader;
import support.CLI;
import support.CLIParseException;
import support.FastaReader;

public class NLR_Annotator {

	
	
	public static final double version = 0.7;
	
	Vector<NLR_MotifList> nlrs;
	Hashtable<String,Vector<NLR_Motif>> motifs;
	
	final int lrrMotifRank = 14;
	
	public NLR_Annotator(File choppedNlrParserXML)throws ParserConfigurationException, IOException, SAXException{
		
		nlrs = new Vector<NLR_MotifList>();
		
		Hashtable<String, MastMotifHitList> motifLists = readXML(choppedNlrParserXML);
		this. motifs = consolidateLists(motifLists);
		this.findNLRs(500, 2500, 50000);
		
		
		
	}
	
public NLR_Annotator(File choppedNlrParserXML, int distanceWithinMotifCombination, int distanceForElongating, int distanceBetweenMotifCombinations)throws ParserConfigurationException, IOException, SAXException{
		
		nlrs = new Vector<NLR_MotifList>();
		
		Hashtable<String, MastMotifHitList> motifLists = readXML(choppedNlrParserXML);
		this. motifs = consolidateLists(motifLists);
		this.findNLRs(distanceWithinMotifCombination, distanceForElongating, distanceBetweenMotifCombinations);
		
		
		
	}
	
	
	public double getVersion(){
		return this.version;
	}
	
	public int getNumberOfNLRs(){
		return nlrs.size();
	}
	
	
	public Vector<NLR_MotifList> getNLRs(){
		return nlrs;
	}
	
	public Hashtable<String,Vector<NLR_Motif>> getMotifs(){
		return this.motifs;
	}
	
	public void findNLRs(int distanceWithinMotifCombination,int distanceForElongating, int distanceBetweenMotifCombinations){
		
		for(Enumeration<String> myenum1 = this.motifs.keys(); myenum1.hasMoreElements();){
			String contigName = myenum1.nextElement();
		//	String contigName = "scaffold20093";   //debug
			Vector<NLR_Motif> v = motifs.get(contigName);
			
			Vector<NLR_Motif> forwardMotifs = new Vector<NLR_Motif>();
			Vector<NLR_Motif> reverseMotifs = new Vector<NLR_Motif>();
			
			for(Enumeration<NLR_Motif> myenum2 = v.elements(); myenum2.hasMoreElements();){
				NLR_Motif motif = myenum2.nextElement();
				if(!motif.isForwardStrand()){
					reverseMotifs.add(motif);
					
				}else{
					forwardMotifs.add(motif);
				}
			}
			Collections.sort(forwardMotifs);
			Collections.sort(reverseMotifs);
			Collections.reverse(reverseMotifs);
			
			
			
			
			
			
			Vector<Vector<NLR_Motif>> forwardNibblers = motifs2NLRs( getCombinationList(forwardMotifs, distanceWithinMotifCombination), distanceForElongating, distanceBetweenMotifCombinations);
			Vector<Vector<NLR_Motif>> reverseNibblers = motifs2NLRs( getCombinationList(reverseMotifs, distanceWithinMotifCombination), distanceForElongating, distanceBetweenMotifCombinations);
			
			int counta = 1;
			
			for(Enumeration<Vector<NLR_Motif>> myenum2 = forwardNibblers.elements(); myenum2.hasMoreElements();){
				
				NLR_MotifList list = new NLR_MotifList(contigName, contigName +"_nlr_"+counta, myenum2.nextElement());
				list.removeRedundantMotifs();
				list = findAttachedMotifs(list, forwardMotifs, distanceForElongating);
				list = findIntermediateMotifs(list, forwardMotifs);
				list = removeDoubleMotifs(list);
				this.nlrs.add(list);
				counta++;
			}
			for(Enumeration<Vector<NLR_Motif>> myenum2 = reverseNibblers.elements(); myenum2.hasMoreElements();){
				NLR_MotifList list = new NLR_MotifList(contigName, contigName +"_nlr_"+counta, myenum2.nextElement());
				list.removeRedundantMotifs();
				list = findAttachedMotifs(list, reverseMotifs, distanceForElongating);
				list = findIntermediateMotifs(list, reverseMotifs);
				list = removeDoubleMotifs(list);
				this.nlrs.add(list);
				counta++;
			}
		
		}
		
	}
	
	
	
	private NLR_MotifList removeDoubleMotifs(NLR_MotifList list){
		
		NLR_MotifList newlist = new NLR_MotifList(list.getSequenceName(), list.getNlrID());
		
		int lastMotifRank = -1;
		for( Enumeration<NLR_Motif > myenum = list.getMotifs().elements(); myenum.hasMoreElements();){
			NLR_Motif motif = myenum.nextElement();
			if( this.getMotifRank(motif) > lastMotifRank || (this.getMotifRank(motif) == lastMotifRank && lastMotifRank == this.lrrMotifRank) ){
				newlist.addMotif(motif);
				
			}else{
				//System.out.println(list.getSequenceName() + "\t" + motif.getMotifNumber() + "\t" + list.getMotifs().firstElement().getDnaStart());
			}
			lastMotifRank = this.getMotifRank(motif);
		}
		
		
		return newlist;
	}
	
	private NLR_MotifList findIntermediateMotifs(NLR_MotifList list, Vector<NLR_Motif> wholeList){
		
		int index = 0;
		int end = list.getMotifs().size()-1;
		
		while( index < end  ){
			NLR_Motif m1 = list.getMotifs().get(index);
			NLR_Motif m2 = list.getMotifs().get(index + 1);
			
			for( int i = 0; i< wholeList.size(); i++){
				NLR_Motif motif = wholeList.get(i);
				
				if(motif.getDnaStart() > Math.min(m1.getDnaStart(), m2.getDnaStart()) && motif.getDnaStart() < Math.max(m1.getDnaStart(), m2.getDnaStart()) ){
					if( this.getMotifRank(m1.getMotifNumber()) < this.getMotifRank(motif.getMotifNumber())   && this.getMotifRank(m2.getMotifNumber()) > this.getMotifRank(motif.getMotifNumber())){
						list.addMotif(motif);
						end++;
						index--;
						break;
					}
					
				}
				
			}
			
			index++;
			
		}
		
		
		
		
		return list;
	}
	
	
	private NLR_MotifList findAttachedMotifs(NLR_MotifList list, Vector<NLR_Motif> wholeList, int distanceForElongating){
		
		
		NLR_Motif firstMotif = list.getMotifs().firstElement();
		NLR_Motif lastMotif = list.getMotifs().lastElement();
		
		int firstIndex = -1;
		int lastIndex = -1;
		for( int i = 0; i< wholeList.size(); i++){
			NLR_Motif motif = wholeList.get(i);
			if(motif.getMotifNumber() == firstMotif.getMotifNumber() && motif.getDnaStart() == firstMotif.getDnaStart()){
				firstIndex = i;
			}
			if(motif.getMotifNumber() == lastMotif.getMotifNumber() && motif.getDnaStart() == lastMotif.getDnaStart()){
				lastIndex = i;
				break;
			}
		}
		
		
		for( int i =firstIndex-1; i >=0; i--){
			NLR_Motif motif = wholeList.get(i);
			if(  Math.abs(firstMotif.getDnaStart() - motif.getDnaStart()) <= distanceForElongating  ){
				if(this.getMotifRank(motif) < this.getMotifRank(firstMotif)){
					list.addMotif(motif);
					firstMotif = motif;
				}
			}else{
				break;
			}
		}
		
		for( int i = lastIndex+1; i< wholeList.size(); i++){
			NLR_Motif motif = wholeList.get(i);
			if(  Math.abs(lastMotif.getDnaStart() - motif.getDnaStart()) <= distanceForElongating  ){
				if( this.getMotifRank(motif) > this.getMotifRank(lastMotif) ||    (this.getMotifRank(motif) == this.getMotifRank(lastMotif) && this.getMotifRank(motif)==lrrMotifRank  ) ){
					list.addMotif(motif);
					lastMotif = motif;
				}
			}
			
		}
		
		
		
		return list;
	}
	
	
	
	private Vector<Vector<NLR_Motif>> motifs2NLRs(Vector<NLR_Motif[]> motifCombinations, int distanceForElongating, int distanceBetweenMotifCombinations){
		Vector<Vector<NLR_Motif>> nlrs = new Vector<Vector<NLR_Motif>>();
		
		int lastRank = -2;
		long lastPosition = -1;
		
		Vector<NLR_Motif> putative_nlr = new Vector<NLR_Motif>();
		for(Enumeration<NLR_Motif[]> myenum = motifCombinations.elements(); myenum.hasMoreElements();){
			NLR_Motif[] a = myenum.nextElement();
			/*
			System.out.print(lastRank+":"+lastPosition + "\t" + getMotifRank(a[0])+":"+ a[0].getDnaStart() + "\t" + a[0].getMotifNumber());
			for( int i = 1; i< a.length; i++){
				System.out.print(","+a[i].getMotifNumber());
			}
			System.out.println();
			*/
			
			
			if(  (getMotifRank(a[0].getMotifNumber()) <= lastRank && getMotifRank(a[0].getMotifNumber())!=lrrMotifRank) || Math.abs(a[0].getDnaStart()-lastPosition) > distanceBetweenMotifCombinations  ){
						
				
				
				
				if( hasNBARC(putative_nlr)){
					nlrs.add(putative_nlr); 
				}
				
				putative_nlr = new Vector<NLR_Motif>();
			}
			
			lastRank = getMotifRank(a[0].getMotifNumber());
			lastPosition = a[a.length-1].getDnaStart();
			
			if( isLrrMotifCombination(a)){  //if we are at LRRs, the rank can be even. Thus we reduce the ran if we are at LRRs
				lastRank--;
			}
			for( int i = 0; i< a.length; i++){
				putative_nlr.add(a[i]);
			}
		}
		
		
		
		
		if( putative_nlr.size()>0   && hasNBARC(putative_nlr)){
			nlrs.add(putative_nlr);
		}
		
		return nlrs;
	}
	
	

	
	
	/*
	private String getMotifCombinationString(NLR_Motif[] a){
		String s = a[0].getMotifNumber() +"";
		for( int i = 1; i< a.length; i++){
			s = s + "," + a[i].getMotifNumber();
		}
		return s;
	}
	*/
	
	private boolean isLrrMotifCombination(NLR_Motif[] a){
		for( int i = 0; i< a.length; i++){
			if( getMotifRank(a[i].getMotifNumber()) != lrrMotifRank){
				return false;
			}
		}
		return true;
	}
	private int getMotifRank(NLR_Motif motif){
		return getMotifRank(motif.getMotifNumber());
	}
	
	private boolean isNbarcMotif(NLR_Motif motif){
		int i = motif.getMotifNumber();
		if( i == 1){return true;}
		if( i == 6){return true;}
		if( i == 4){return true;}
		if( i == 5){return true;}
		if( i == 10){return true;}
		if( i == 3){return true;}
		if( i == 12){return true;}
		if( i == 2){return true;}
		
		
		return false;
	}
	
	private boolean isLRRMotif(NLR_Motif motif){
		int i = motif.getMotifNumber();
		if( i == 9){return true;}
		if( i == 11){return true;}
		if( i == 19){return true;}
		return false;
	}
	
	private int getMotifRank(int motifNumber){
		if(motifNumber == 1){ return 4;}
		if(motifNumber == 2){ return 11;}
		if(motifNumber == 3){ return 9;}
		if(motifNumber == 4){ return 6;}
		if(motifNumber == 5){ return 7;}
		if(motifNumber == 6){ return 5;}
		if(motifNumber == 7){ return 13;}
		if(motifNumber == 8){ return 12;}
		if(motifNumber == 9){ return lrrMotifRank;}
		if(motifNumber == 10){ return 8;}
		if(motifNumber == 11){ return lrrMotifRank;}
		if(motifNumber == 12){ return 10;}
		if(motifNumber == 13){ return 3;}
		if(motifNumber == 14){ return 3;}
		if(motifNumber == 15){ return 2;}
		if(motifNumber == 16){ return 2;}
		if(motifNumber == 17){ return 1;}
		if(motifNumber == 18){ return 1;}
		if(motifNumber == 19){ return lrrMotifRank;}
		if(motifNumber == 20){ return 15;}
		return -1;
	}
	
	
	private boolean containsCCMotif(String motifCombinationString){
		if( motifCombinationString.contains("17") || motifCombinationString.contains("16") || motifCombinationString.contains(",6") || motifCombinationString.contains(",2") ){
			return true;
		}
		return false;
	}
	
	private boolean containsTIRMotif(String motifCombinationString){
		if( motifCombinationString.contains("18") || motifCombinationString.contains("15") || motifCombinationString.contains("13")){
			return true;
		}
		return false;
	}
	
	private Vector<NLR_Motif[]> getCombinationList(Vector<NLR_Motif> mylist, int distanceWithinMotifCombination){
		Vector<NLR_Motif[]> v = new Vector<NLR_Motif[]>();
		for( int i = 0; i< mylist.size(); i++){
			boolean foundOne = false;
			if( i+2< mylist.size() && Math.abs(mylist.get(i).getDnaStart()-mylist.get(i+1).getDnaStart())< distanceWithinMotifCombination  && Math.abs(mylist.get(i+1).getDnaStart()-mylist.get(i+2).getDnaStart())<distanceWithinMotifCombination){
				
				String s = mylist.get(i).getMotifNumber() + "," + mylist.get(i+1).getMotifNumber() + "," + mylist.get(i+2).getMotifNumber();
				
				
				int d1 = getMotifRank(mylist.get(i+1).getMotifNumber()) - getMotifRank(mylist.get(i).getMotifNumber());
				int d2 = getMotifRank(mylist.get(i+2).getMotifNumber()) - getMotifRank(mylist.get(i+1).getMotifNumber());
				
				if(     ( d1 > 0 && d1 < 3   )
				      &&( d2 > 0 && d2 < 3  || (d2 ==0 && getMotifRank(mylist.get(i+1)) ==lrrMotifRank ) )
					   //if CC then no TIR motif OR no CC
					&&	 (  (containsCCMotif(s) && !containsTIRMotif(s)) || !containsCCMotif(s)       )	
				){
					
					NLR_Motif[] a = {mylist.get(i), mylist.get(i+1), mylist.get(i+2)};
					foundOne = true;
					v.add(a);
				}
			}
			
			
			if(!foundOne && i+1<mylist.size() && Math.abs(mylist.get(i).getDnaStart() - mylist.get(i+1).getDnaStart()) <  distanceWithinMotifCombination  ){
					String s = mylist.get(i).getMotifNumber() + "," + mylist.get(i+1).getMotifNumber();
					
					if( 	s.equalsIgnoreCase("18,15")
						||	s.equalsIgnoreCase("15,13")
						||	s.equalsIgnoreCase("13,1")
						||	s.equalsIgnoreCase("17,16")
						||	s.equalsIgnoreCase("16,1")
						||	s.equalsIgnoreCase("17,1")
						||	s.equalsIgnoreCase("1,6")
						
					){
						
						NLR_Motif[] a = {mylist.get(i), mylist.get(i+1)};
						v.add(a);
						foundOne = true;
					}
				}
			
			
		}
		
		
		
		
		return v;
	
	}

	
	
	private boolean hasNBARC(Vector<NLR_Motif> v){
		String s = "";
		for(Enumeration<NLR_Motif> myenum= v.elements(); myenum.hasMoreElements();){
			s= s + "," + myenum.nextElement().getMotifNumber();
		}
		
		if(s.contains(",1,6,4")){return true;}
		if(s.contains(",6,4,5")){return true;}
		if(s.contains(",4,5,10")){return true;}
		if(s.contains(",5,10,3")){return true;}
		if(s.contains(",10,3,12")){return true;}
		if(s.contains(",3,12,2")){return true;}
		if(s.contains(",1,4,5")){return true;}
		if(s.contains(",12,2,8")){return true;}
		if(s.contains(",2,8,7")){return true;}
		return false;
		
	}
	
	
	
	private Hashtable<String, MastMotifHitList> readXML(File xmlFile)throws ParserConfigurationException, IOException, SAXException{
		Hashtable<String, MastMotifHitList> motifLists = new Hashtable<String, MastMotifHitList>();
		
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		Document dom = db.parse(xmlFile);
		
		
		Element rootElement = dom.getDocumentElement();
		NodeList nodeList = rootElement.getElementsByTagName("MastMotifHitList");
		for( int i = 0; i< nodeList.getLength(); i++){
			Element mastMotifHitList = (Element) nodeList.item(i);
			MastMotifHitList list = new MastMotifHitList(mastMotifHitList);
			motifLists.put(mastMotifHitList.getAttribute("name"), list);
		}
		
		return motifLists;
	}
	
	
	private Hashtable<String,Vector<NLR_Motif>> consolidateLists(Hashtable<String, MastMotifHitList> motifLists){
		Hashtable<String,Vector<NLR_Motif>> h = new Hashtable<String,Vector<NLR_Motif>>();
		
		for(Enumeration<String> myenum1 = motifLists.keys(); myenum1.hasMoreElements();){
			String key = myenum1.nextElement();
			
			String[] split = key.split("_");
			String sequenceID = split[0];
			for( int i =1; i<split.length-1;i++){
				sequenceID = sequenceID + "_" + split[i];
 			}
			
			long offset = Long.parseLong(split[split.length-1].substring(0,split[split.length-1].length()-1));
			boolean forwardStrand = true;
			if( split[split.length-1].substring(split[split.length-1].length()-1).equalsIgnoreCase("-")){
				forwardStrand = false;
			}
			
			
			
			MastMotifHitList list = motifLists.get(key);
			String sequence = list.getDNASequence();
			int length = list.getDNASequence().length();
			
			
			//System.out.println(key +"\t" + sequenceID + "\t" + offset + "\t" + forwardStrand);
			
			for(Enumeration<MastMotifHit> myenum2 = list.getMotifs().elements(); myenum2.hasMoreElements();){
				MastMotifHit hit = myenum2.nextElement();
				
				
				int nucl_start = hit.getDNAPos();
				int nucl_end = hit.getDNAPos() + hit.getMatch().length()*3;
				
				
				if(!forwardStrand){
					nucl_start = length - nucl_end;
					nucl_end = length - hit.getDNAPos();
				}
				
				int frame = hit.getFrame();
				long frameshift = offset % 3;
				if(frameshift == 1){
					if(frame == 0){frame =1;}
					if(frame == 1){frame =2;}
					if(frame == 2){frame =0;}
				}
				if(frameshift == 2){
					if(frame == 0){frame =2;}
					if(frame == 1){frame =0;}
					if(frame == 2){frame =1;}
				}
				NLR_Motif motif = new NLR_Motif(sequenceID, hit.getMotif(), nucl_start+offset,nucl_end+offset, hit.getPvalue(), forwardStrand, frame);
				
					System.out.println(sequenceID);
				
					motif.setSequence(sequence.substring(nucl_start, nucl_end)); 
					if(motif.getStrand().equalsIgnoreCase("-")){
						motif.setSequence(new BioSequence("", sequence).getReverseComplementarySequence().substring(nucl_start, nucl_end));
					}
	
				
				Vector<NLR_Motif> v = new Vector<NLR_Motif>();
				if( h.containsKey(sequenceID)){
					v = h.get(sequenceID);
				}
				v.add(motif);
				Collections.sort(v);
				h.put(sequenceID, v);
				
			}
			
			
		}
		
		Hashtable<String, Vector<NLR_Motif>> h1 = new Hashtable<String, Vector<NLR_Motif>> ();
		for(Enumeration<String> myenum1 = h.keys(); myenum1.hasMoreElements();){
			String key = myenum1.nextElement();
			Vector<NLR_Motif> v = h.get(key);
			Vector<NLR_Motif> v1 = new Vector<NLR_Motif> ();
			v1.add(v.get(0));
			int lastMotif = v.get(0).getMotifNumber();
			long lastPosition = v.get(0).getDnaStart();
			for(Enumeration<NLR_Motif> myenum2 = v.elements(); myenum2.hasMoreElements();){
				NLR_Motif motif = myenum2.nextElement();
				if(motif.getMotifNumber() == lastMotif && motif.getDnaStart() == lastPosition ){
					continue;
				}else{
					v1.add(motif);
					lastMotif = motif.getMotifNumber();
					lastPosition = motif.getDnaStart();
				}
			}
			
			
			h1.put(key, v1);
		}
		
		
		
		
		return h1;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/* ***************************************************************************
	 *                   Output Methods
	 *****************************************************************************/
	
	
	
	public void writeNlrGFF(File outputFile, boolean completeOnly)throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		out.write("##gff-version 2");
		out.newLine();
		out.write("##source-version NLR-Annotator 2" );
		out.newLine();
		
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm");
		Date date = new Date();
		out.write("##date "+dateFormat.format(date)); 
		out.newLine();
		out.write("##Type DNA");
		out.newLine();
		out.write("#seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute");
		out.newLine();
		
		for(Enumeration <NLR_MotifList> myenum = this.nlrs.elements(); myenum.hasMoreElements();){
			
			NLR_MotifList list = myenum.nextElement();
			
			if( !completeOnly || list.isCompleteNLR()){
				out.write(list.getSequenceName() +"\tNLR_Annotator\tNBSLRR\t" );
				
				String strand = "+";
				long start = Math.min( Math.min(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.min( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
				long end = Math.max( Math.max(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.max( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
				
				start = start + 1; //GFF is 1-based!
				
				
				if( !list.getMotifs().get(0).isForwardStrand()){
					strand = "-";
				}
				out.write(start + "\t" + end + "\t.\t"+strand+"\t.\tname="+ list.getNlrID());
				out.newLine();
			}
			
			
		}
		
		
		
		out.close();
		
		
	}
	
	
	public void writeNlrBED(File outputFile)throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		out.write("#track name=\"NLR_Loci\"");
		out.newLine();
		out.write("#itemRgb=\"On\"");
		out.newLine();
		
		for(Enumeration <NLR_MotifList> myenum = this.nlrs.elements(); myenum.hasMoreElements();){
			
			NLR_MotifList list = myenum.nextElement();
			
			
			
			String strand = "+";
			
			//positions in BED are supposed to be zero-based. 
			long start = Math.min( Math.min(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.min( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
			long end = Math.max( Math.max(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.max( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
			
			
			
			if( !list.getMotifs().get(0).isForwardStrand()){
				strand = "-";
			}
			
			
			String blockCount = list.getMotifs().size()+"";
			String blockStarts = "";
			String blockSizes = "";
			
			
			
			for(Enumeration<NLR_Motif> myenum2 = list.getMotifs().elements(); myenum2.hasMoreElements();){
				NLR_Motif motif = myenum2.nextElement();
				
				if( list.getMotifs().firstElement().isForwardStrand()){
					blockStarts = blockStarts + "," + ( motif.getDnaStart()  - start);
					blockSizes = blockSizes + "," + (motif.getDnaEnd()-motif.getDnaStart());
				}else{
					blockStarts =  "," + ( motif.getDnaStart()  - start) + blockStarts ;
					blockSizes = "," + (motif.getDnaEnd()-motif.getDnaStart()) + blockSizes;
				}
				
			}
			blockStarts = blockStarts.substring(1);
			blockSizes = blockSizes.substring(1);
			
			
			
			
			String color = "0,255,0";//green
			if(!list.isCompleteNLR()){
				color = "255,128,0";//orange
			}
			
			if( list.hasStopCodon()){
				color = "255,0,0";
			}
			
			out.write(list.getSequenceName() +"\t" + start +"\t" + end +"\t"+list.getNlrID()+"\t0\t"+ strand +"\t" + start +"\t" + end +"\t" + color + "\t" + blockCount + "\t" + blockSizes+"\t" +blockStarts );
			out.newLine();
		
			
		
		}
		
		
		
		out.close();
		
		
	}
	
	public void writeMotifGFF(Hashtable<String,Vector<NLR_Motif>> motifs, File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("##gff-version 2");
		out.newLine();
		out.write("##source-version NLR-Annotator " + this.getVersion() );
		out.newLine();
		
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm");
		Date date = new Date();
		out.write("##date "+dateFormat.format(date)); 
		out.newLine();
		out.write("##Type DNA");
		out.newLine();
		out.write("#seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute");
		out.newLine();
		
		for(Enumeration<String> myenum1 = motifs.keys(); myenum1.hasMoreElements();){
			String key = myenum1.nextElement();
			Vector<NLR_Motif> v = motifs.get(key);
			for(Enumeration<NLR_Motif> myenum2 = v.elements();myenum2.hasMoreElements();){
				NLR_Motif motif = myenum2.nextElement();
				out.write(motif.getGFFString());
				out.newLine();
				
				
			}
			
		}
		
		
		out.close();
	}
	
	public void writeMotifBED(Hashtable<String,Vector<NLR_Motif>> motifs, File outputFile, boolean annotateSTOP)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("#track name=\"NLR_Motifs\"");
		out.newLine();
		out.write("#itemRgb=\"On\"");
		out.newLine();
		
		for(Enumeration<String> myenum1 = motifs.keys(); myenum1.hasMoreElements();){
			String key = myenum1.nextElement();
			Vector<NLR_Motif> v = motifs.get(key);
			for(Enumeration<NLR_Motif> myenum2 = v.elements();myenum2.hasMoreElements();){
				NLR_Motif motif = myenum2.nextElement();
				out.write(motif.getBEDString());
				out.newLine();
				if( annotateSTOP && motif.getAASequence().contains("*")){
					out.write(motif.getBlackMotif());
					out.newLine();
				}
				
			}
			
		}
		
		
		out.close();
	}
	
	
	public void writeReportTxt( File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		for(Enumeration <NLR_MotifList> myenum = this.nlrs.elements(); myenum.hasMoreElements();){
			NLR_MotifList list = myenum.nextElement();
			long start = Math.min( Math.min(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.min( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
			long end = Math.max( Math.max(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.max( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
		   String strand = "+";
			if( !list.getMotifs().get(0).isForwardStrand()){
				strand = "-";
			}
			String complete = "complete";
			if( !list.isCompleteNLR()){
				complete = "partial";
			}
			
			if( list.hasStopCodon()){
				complete = complete +" (pseudogene)";
			}
			
			out.write(list.getSequenceName()+"\t"+list.getNlrID() +"\t"+complete+"\t" + start + "\t" + end +"\t" + strand + "\t" + list.getMotifList());
			out.newLine();
			
			
		}	
		out.close();
		
	}
	
	
	
	public void writeNbarcMotifAlignment(File outputFasta)throws IOException{
		/*
		1	6	4	5	10	3	12	2	
		PIWGMGGVGKTTLARAVYNDP
		HFDCRAWVCVSQQYDMKKVLRDIIQQVGG
		YLVVLDDVWDTDQWD
		NGSRIIITTRNKHVANYMCT
		LSHEESWQLFHQHAF
		CGGLPLAIKVWGGMLAGKQKT
		IMPVLRLSYHHLPYH
		LKPCFLYCAIFPEDYMIDKNKLIWLWMAE
	*/
		
		String motif1 = "PIWGMGGVGKTTLARAVYNDP";
		String motif6 = "HFDCRAWVCVSQQYDMKKVLRDIIQQVGG";
		String motif4 = "YLVVLDDVWDTDQWD";
		String motif5 = "NGSRIIITTRNKHVANYMCT";
		String motif10 = "LSHEESWQLFHQHAF";
		String motif3 = "CGGLPLAIKVWGGMLAGKQKT";
		String motif12 = "IMPVLRLSYHHLPYH";
		String motif2 = "LKPCFLYCAIFPEDYMIDKNKLIWLWMAE";
		String[] motifReference = {motif1, motif6, motif4, motif5, motif10, motif3, motif12, motif2};
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFasta));
		
		out.write("#mega");
		out.newLine();
		out.write("!Title annotator5_nbarcmotifs alignment;");
		out.newLine();
		out.write("!Format DataType=protein;");
		out.newLine();
		out.newLine();
		
		out.write("#reference_sequence");
		out.newLine();
		out.write(motifReference[0]);
		for( int i = 1; i< motifReference.length;i++){
			out.write( motifReference[i]);
		}
		
		out.newLine();
		out.newLine();
		
		for(Enumeration<NLR_MotifList> myenum1 = this.nlrs.elements(); myenum1.hasMoreElements();){
			
			
			
			
			NLR_MotifList list = myenum1.nextElement();
			/*
			if( !list.getNlrID().equalsIgnoreCase("scaffold18474_nlr_1")){
				continue;
			}
			*/
			//System.out.println(list.getMotifList());
			if(!list.isCompleteNLR()){
				continue;
			}
			
			
			String sequence = "";
			
			int index =0;
			NLR_Motif motif = list.getMotifs().get(0);
			while(motif.getMotifNumber()!=1){
				index++;
				motif = list.getMotifs().get(index);
			}
			
			sequence = motif.getAASequence() ;
			
			index++;
			
			
			int[] motifOrder = {1,	6,	4,	5,	10,	3,	12,	2	};
			for( int j = 1; j<motifOrder.length; j++){
				motif = list.getMotifs().get(index);
				if( motif.getMotifNumber() != motifOrder[j]){
					String s = "";
					for( int i= 0; i<motifReference[j].length(); i++){
						s = s +"-";
					}
					sequence = sequence +s ;
					
				}else{
					sequence = sequence + motif.getAASequence();
					if (motifReference[j].length() != motif.getAASequence().length()){
						System.out.println( list.getNlrID() + "\t" + motif.getAASequence() );
					}
					index ++;
					
					
				}
				
			}
			
			
			
			
			
			
			
			
			if( !sequence.contains("*")){
				out.write("#" + list.getNlrID());
				out.newLine();
				out.write(sequence);
				out.newLine();
				out.newLine();
			}
			
			
			
			
			
			
			
			
			
		}
		
		
		out.close();
		
		
		
		
	}
	
	
	public void writeNbarcMultipleAlignmentFasta(File outputFasta)throws IOException{
		/*
		1	6	4	5	10	3	12	2	
		PIWGMGGVGKTTLARAVYNDP
		HFDCRAWVCVSQQYDMKKVLRDIIQQVGG
		YLVVLDDVWDTDQWD
		NGSRIIITTRNKHVANYMCT
		LSHEESWQLFHQHAF
		CGGLPLAIKVWGGMLAGKQKT
		IMPVLRLSYHHLPYH
		LKPCFLYCAIFPEDYMIDKNKLIWLWMAE
	*/
		
		String motif1 = "PIWGMGGVGKTTLARAVYNDP";
		String motif6 = "HFDCRAWVCVSQQYDMKKVLRDIIQQVGG";
		String motif4 = "YLVVLDDVWDTDQWD";
		String motif5 = "NGSRIIITTRNKHVANYMCT";
		String motif10 = "LSHEESWQLFHQHAF";
		String motif3 = "CGGLPLAIKVWGGMLAGKQKT";
		String motif12 = "IMPVLRLSYHHLPYH";
		String motif2 = "LKPCFLYCAIFPEDYMIDKNKLIWLWMAE";
		String[] motifReference = {motif1, motif6, motif4, motif5, motif10, motif3, motif12, motif2};
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFasta));
		
		out.write(">NP_001021202.1");
		out.newLine();
        out.write("FLHGRAGSGKSVIASQALSKS-----------------------------TLFVFDDVVQEETIRLRLRCLVTTRDVEISNAASQ"
        		+ "--------------------------------------------------------------------------------");
        out.newLine();
			
		out.write(">reference_sequence");
		out.newLine();
		out.write(motifReference[0]);
		for( int i = 1; i< motifReference.length;i++){
			out.write( motifReference[i]);
		}
		
		out.newLine();
		
		
		for(Enumeration<NLR_MotifList> myenum1 = this.nlrs.elements(); myenum1.hasMoreElements();){
			
			
			
			
			NLR_MotifList list = myenum1.nextElement();
			
			
			if(!list.isCompleteNLR()){
				continue;
			}
			
			
			String sequence = "";
			
			int index =0;
			NLR_Motif motif = list.getMotifs().get(0);
			while(motif.getMotifNumber()!=1){
				index++;
				motif = list.getMotifs().get(index);
			}
			
			sequence = motif.getAASequence() ;
			
			index++;
			
			
			int[] motifOrder = {1,	6,	4,	5,	10,	3,	12,	2	};
			for( int j = 1; j<motifOrder.length; j++){
				motif = list.getMotifs().get(index);
				if( motif.getMotifNumber() != motifOrder[j]){
					String s = "";
					for( int i= 0; i<motifReference[j].length(); i++){
						s = s +"-";
					}
					sequence = sequence +s ;
					
				}else{
					sequence = sequence + motif.getAASequence();
					if (motifReference[j].length() != motif.getAASequence().length()){
						System.out.println( list.getNlrID() + "\t" + motif.getAASequence() );
					}
					index ++;
					
					
				}
				
			}
			
			
			
			
			
			
			
			
		
			/*   initially I would not use loci with a stop codon. But they are still interesting....
			if( !sequence.contains("*") && !sequence.contains("X")){ //exclude pseudogenes and genes where the sequence has a gap
				
				out.write(new BioSequence( list.getNlrID() , sequence ).getFastaString()   );
			}
			*/
			
			//replace nucleotide gaps and stop codons by gaps in alignment
			sequence.replaceAll("\\*", "_");
			sequence.replaceAll("X", "_");
			out.write(new BioSequence( list.getNlrID() , sequence ).getFastaString()   );
			
			
			
			
			
			
			
		}
		
		
		out.close();
		
		
		
		
	}
	
	public void writeNbarcMotifs(File outputFasta)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFasta));
		
		for(Enumeration<NLR_MotifList> myenum1 = this.nlrs.elements(); myenum1.hasMoreElements();){
						
			NLR_MotifList list = myenum1.nextElement();
			if(!list.isCompleteNLR()){
				continue;
			}
			
			String motifList = "";
			String sequence = "";
			
			for(Enumeration<NLR_Motif> myenum2 = list.getMotifs().elements(); myenum2.hasMoreElements();){
				NLR_Motif motif = myenum2.nextElement();
				if(this.isNbarcMotif(motif)){
					motifList = motifList + ","+motif.getMotifNumber();
					sequence = sequence + motif.getAASequence();
				}
			}
			BioSequence seq = new BioSequence(list.getNlrID(), sequence);
			seq.setDescription("sequence="+list.getSequenceName()+" start="+list.getMotifs().firstElement().getDnaStart()+" motifs="+motifList.substring(1));
			if(!seq.getSequence().contains("*")){
				out.write(seq.getFastaString());
			}
			
			
		}
		
		
		out.close();
		
	}
	
	
	
	
	
	/* ********************************************************************************************
	 * *********** Some diagnostic methods ********************************************************
	 */
	
	
	
	
	public Hashtable<String,Vector<long[]>> getMotifRegions( int motifExtension, boolean completeOnly){
		Hashtable<String,Vector<long[]>> regions = new Hashtable<String,Vector<long[]>>();
		
		for(Enumeration<NLR_MotifList> myenum1 = this.nlrs.elements(); myenum1.hasMoreElements();){
			NLR_MotifList nlr = myenum1.nextElement();
			
			if(!completeOnly && !nlr.isCompleteNLR()){
				continue;
			}
			
			String scaffold = nlr.getSequenceName();
			
			for(Enumeration<NLR_Motif> myenum2 = nlr.getMotifs().elements(); myenum2.hasMoreElements();){
				NLR_Motif motif = myenum2.nextElement();
				long start =Math.max(0, Math.min(motif.getDnaStart(), motif.getDnaEnd())-motifExtension) ;
				long end   = Math.max(motif.getDnaStart(), motif.getDnaEnd())+motifExtension;
				
				
				
				long[] a = {start, end};
				Vector<long[]> v = new Vector<long[]>();
				if(regions.containsKey(scaffold)){
					v = regions.get(scaffold);
				}
				v = this.addRegion(v, a);
				
				
				//System.out.println(v.size());
				regions.put(scaffold, v);
			}
			
		
			
		}
		
		
		return regions;
		
	}
	
	public Hashtable<String,Vector<long[]>> checkCoverage(File mpileupFile, int minCoverage)throws IOException{
		
		
		System.out.print("Extractiong regions for file " + mpileupFile.getName() + ": ");
		long covered = 0;
		Hashtable<String,Vector<long[]>> h = new Hashtable<String,Vector<long[]>>();
		
		String currentScaffold = "";
		long startInterval = -1;
		long lastPosition = -1;
		Vector<long[]> v = new Vector<long[]>();
		BufferedReader in = new BufferedReader(new FileReader(mpileupFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\t");
			String scaffold = split[0];
			if(!scaffold.equalsIgnoreCase(currentScaffold)){
				if(startInterval != -1){
					long[] a = {startInterval, lastPosition};
					covered = covered + 1 + lastPosition - startInterval;
					v.add(a);
					h.put(currentScaffold, v);
					startInterval = -1;
					v = new Vector<long[]>();
					currentScaffold = scaffold;
				}
			}
			
			
			
			long position = Long.parseLong(split[1]);
			int cov = Integer.parseInt(split[3]);
			if( cov>=minCoverage){
				if( startInterval == -1){
					startInterval = position;
				}else{
					if(position != lastPosition + 1){
						long[] a = {startInterval, lastPosition};
						covered = covered + 1 + lastPosition - startInterval;
						v.add(a);
						startInterval = position;
					}
				}
				
				lastPosition = position;
			}
			
		}

		if( startInterval != -1){
			long[] a = {startInterval, lastPosition};
			covered = covered + 1 + lastPosition - startInterval;
			v.add(a);
		}
		if( v.size() > 0){
			h.put(currentScaffold, v);
		}
		
		in.close();
		System.out.println(covered + " bp covered;");
		
		return h;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public Hashtable<String,Vector<long[]>> getBaitCoverage(File blastFile, int minHSPLength, double minIdentity)throws IOException{
		
		Hashtable<String,Vector<long[]>> h = new Hashtable<String,Vector<long[]>>();
		
		BlastReader reader = new BlastIPKReader(blastFile);
		for(BlastIteration it = reader.readIteration(); it != null; it = reader.readIteration()){
			for(Enumeration<BlastHit> myenum1 = it.getHits().elements(); myenum1.hasMoreElements();){
				BlastHit hit = myenum1.nextElement();
				for(Enumeration<BlastHSP > myenum2 = hit.getHSPs().elements(); myenum2.hasMoreElements();){
					BlastHSP hsp = myenum2.nextElement();
					
					if( hsp.getHspLength() >= minHSPLength && hsp.getIdentityPercentage() >= minIdentity){
						long start = Math.min( hsp.getHitStart(), hsp.getHitEnd());
						long end = Math.max( hsp.getHitStart(), hsp.getHitEnd());
						
						Vector<long[]> v = new Vector<long[]>();
						if( h.containsKey(hsp.getHitName())){
							v = h.get(hsp.getHitName());
						}
						long[] a = {start, end};
						
						v= this.addRegion(v, a);
						h.put(hsp.getHitName(), v);
					}
					
				}
				
			}
			
			
		}
		
		long positionsCoveredByBaits = 0;
		
		for(Enumeration<Vector<long[]> > myenum1 = h.elements(); myenum1.hasMoreElements(); ){
			for(Enumeration<long[]> myenum2 = myenum1.nextElement().elements(); myenum2.hasMoreElements();){
				long[] a = myenum2.nextElement();
				positionsCoveredByBaits = positionsCoveredByBaits + 1 + a[1] - a[0];
			}
		}
		
		System.out.println(blastFile.getName() + " has " + positionsCoveredByBaits + " positions covered by baits");
		
		return h;
	}
	
	private static Vector<long[]> addRegion(Vector<long[]> v, long[] a ){
		Collections.sort(v , new Comparator<long[]>(){public int compare(long[] a, long[] b){if(a[0] <b[0]){	return -1;}	if(a[0] >b[0]){	return 1;}return 0;}} );
		Vector<long[] > w = new Vector<long[]>();
		
		
		boolean added = false;
		for(Enumeration<long[]> myenum = v.elements(); myenum.hasMoreElements();){
			long[] b = myenum.nextElement();
			
			if( (b[0] <= a[0] && b[1] >= a[0])  ||   (b[0] <= a[1] && b[1]>=a[1] )   ||  (a[0]< b[0] && a[1]>b[1])  ){
				added = true;
				long[] c = {  Math.min(a[0], b[0])  , Math.max(a[1], b[1])};
				
				boolean added2 = false;
				while(myenum.hasMoreElements()){
					long[] d = myenum.nextElement(); //check if the new region bridges the two old.
					if( d[0] <= c[1] ){
						c[0] = Math.min(c[0], d[0]);
						c[1] = Math.max(c[1], d[1]);
					}else{
						w.add(c);
						w.add(d);
						added2 = true;
						break;
					}
				}
				if(!added2){
					w.add(c);
				}
				
			}else{
				w.add(b);
			}
			
			
		}
		
		if(!added){
			w.add(a);
		}
		return w;
		
	}
	
	
	
	
	
	public Hashtable<String,Vector<long[]>> compareRegions( Hashtable<String,Vector<long[]>>  refRegions, Hashtable<String,Vector<long[]>> testRegions){
		System.out.print("Compare regions");
		Hashtable<String,Vector<long[]>> h = new Hashtable<String,Vector<long[]>>();
		
		long refPositions = 0;
		long coveredPositions = 0;
		
		
		for(Enumeration<String> myenum1 = refRegions.keys(); myenum1.hasMoreElements();){
			String scaffold = myenum1.nextElement();
			Vector<long[]> v = refRegions.get(scaffold);
			Vector<long[]> x = new Vector<long[]>();
			for(Enumeration<long[]> myenum2=v.elements(); myenum2.hasMoreElements();){
				long[] a = myenum2.nextElement();
				refPositions = refPositions + a[1] - a[0] + 1;
				Vector<long[]> w = testRegions.get(scaffold);
				if( w!= null){
					for(Enumeration<long[]> myenum3 = w.elements(); myenum3.hasMoreElements();){
						long[] b = myenum3.nextElement();
						if( Math.max(b[1],a[1]) - Math.min(a[0],b[0]) < (a[1]-a[0]) + (b[1]-b[0]) ){//they overlap
							coveredPositions = coveredPositions + Math.min(a[1],b[1])+1 - Math.max(a[0], b[0]);
							long[] c = {Math.max(a[0], b[0]), Math.min(a[1],b[1])};
							x.add(c);
						}
					}
				}
				
				
			}
			if(x.size()>0){
				h.put(scaffold, x);
			}
			
			
		}
		
		System.out.println(": " + coveredPositions +" out of "+refPositions+"(" + ((double)coveredPositions *100/refPositions)+"%)." );
		
		return h;
	}
	
	
	
	
	public void writeNLRLoci(File genomeSequence, File outputFile, int flankingSequence)throws IOException{
		Hashtable<String, Vector<NLR_MotifList>> h = new Hashtable<String, Vector<NLR_MotifList>>();
		
		for(Enumeration<NLR_MotifList> myenum = nlrs.elements(); myenum.hasMoreElements();){
			NLR_MotifList list = myenum.nextElement();
			String scaffold = list.getSequenceName();
			if( !h.containsKey(scaffold)){
				h.put(scaffold,new Vector<NLR_MotifList>());
			}
			h.get(scaffold).add(list);
		}
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		BufferedReader in;
		
		FileInputStream fis = new FileInputStream(genomeSequence);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		if(gzip){
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(genomeSequence))));
		}else{
			in = new BufferedReader((new FileReader(genomeSequence)));
		}		
		
		
		
		int baseCharacter=in.read();
		String identifier= "";
		StringBuffer sequence = new StringBuffer(); 
		boolean inSequence = false;
		while( baseCharacter > 0){
			char c = (char)baseCharacter;
			if( c == '>'){
				
				if(inSequence){
					
					Vector<NLR_MotifList> v = h.get(identifier);
					for( Enumeration<NLR_MotifList> myenum = v.elements(); myenum.hasMoreElements();){
						NLR_MotifList list = myenum.nextElement();
						
						
						int p1 = (int) list.getMotifs().firstElement().getDnaStart();
						int p2 = (int) list.getMotifs().firstElement().getDnaEnd();
						int p3 = (int) list.getMotifs().lastElement().getDnaStart();
						int p4 = (int) list.getMotifs().lastElement().getDnaEnd();
						
						int start = Math.max(0 , Math.min(Math.min(p1, p2), Math.min(p3, p4)) - flankingSequence);
						int end = Math.min(sequence.length(), Math.max(Math.max(p1, p2), Math.max(p3, p4)) + flankingSequence);
						
						BioSequence seq = new BioSequence(list.getNlrID(), sequence.substring(start, end));
						String strand = "strand:+";
						if(!list.getMotifs().firstElement().isForwardStrand()){
							strand = "strand:-";
							seq.setSequence(seq.getReverseComplementarySequence());
						}
						seq.setDescription(identifier + " " + start + "-" + end + " " +strand + " " + list.getMotifList());
						out.write(seq.getFastaString(100));
						
					}
					
				}
				
				String header = in.readLine();
				identifier = header.split("\\s")[0];
				sequence = new StringBuffer();
				if(h.containsKey(identifier)){
					inSequence = true;
				}else{
					inSequence = false;
				}
				baseCharacter = in.read();
			}
			if(c != '\n' && inSequence && c != '>'){
				sequence = sequence.append(c);
			}
			
			baseCharacter = in.read();
		}
		in.close();
		
		if(inSequence){
			Vector<NLR_MotifList> v = h.get(identifier);
			for( Enumeration<NLR_MotifList> myenum = v.elements(); myenum.hasMoreElements();){
				NLR_MotifList list = myenum.nextElement();
				
				
				int p1 = (int) list.getMotifs().firstElement().getDnaStart();
				int p2 = (int) list.getMotifs().firstElement().getDnaEnd();
				int p3 = (int) list.getMotifs().lastElement().getDnaStart();
				int p4 = (int) list.getMotifs().lastElement().getDnaEnd();
				
				int start = Math.min(Math.min(p1, p2), Math.min(p3, p4)) - flankingSequence;
				int end = Math.max(Math.max(p1, p2), Math.max(p3, p4)) + flankingSequence;
				//System.out.println(list.getNlrID() +"\t"+sequence.length());
				BioSequence seq = new BioSequence("","");
				try{
					 seq = new BioSequence(list.getNlrID(), sequence.substring(start, end));
				}catch(StringIndexOutOfBoundsException e){
					 seq = new BioSequence(list.getNlrID(), sequence.substring(Math.max(0,start), Math.min(end, sequence.length())));
				}
				String strand = "strand:+";
				if(!list.getMotifs().firstElement().isForwardStrand()){
					strand = "strand:-";
					seq.setSequence(seq.getReverseComplementarySequence());
				}
				seq.setDescription(identifier + " " + start + "-" + end + " " +strand + " " + list.getMotifList());
				out.write(seq.getFastaString(100));
				
			}
		}
		
		
		
		
		out.close();
		
		
	}
	
	
	
	public void writeNLRLoci(File genomeSequence, File outputFile, int motifFlanking, int lrrMotifFlanking, int flankingRegion , boolean predictGeneEnds)throws IOException{
		Hashtable<String, Vector<NLR_MotifList>> h = new Hashtable<String, Vector<NLR_MotifList>>();
		
		for(Enumeration<NLR_MotifList> myenum = nlrs.elements(); myenum.hasMoreElements();){
			NLR_MotifList list = myenum.nextElement();
			String scaffold = list.getSequenceName();
			if( !h.containsKey(scaffold)){
				h.put(scaffold,new Vector<NLR_MotifList>());
			}
			h.get(scaffold).add(list);
		}
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		BufferedReader in;
		
		FileInputStream fis = new FileInputStream(genomeSequence);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		if(gzip){
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(genomeSequence))));
		}else{
			in = new BufferedReader((new FileReader(genomeSequence)));
		}		
		
		
		
		int baseCharacter=in.read();
		String identifier= "";
		StringBuffer sequence = new StringBuffer(); 
		boolean inSequence = false;
		while( baseCharacter > 0){
			char c = (char)baseCharacter;
			if( c == '>'){
				
				if(inSequence){
					
					Vector<NLR_MotifList> v = h.get(identifier);
					writeSequences(out, sequence, v,motifFlanking,  lrrMotifFlanking,  flankingRegion);
				}
				
				String header = in.readLine();
				identifier = header.split("\\s")[0];
				sequence = new StringBuffer();
				if(h.containsKey(identifier)){
					inSequence = true;
				}else{
					inSequence = false;
				}
				baseCharacter = in.read();
			}
			if(c != '\n' && inSequence && c != '>'){
				sequence = sequence.append(c);
			}
			
			baseCharacter = in.read();
		}
		in.close();
		
		if(inSequence){
			Vector<NLR_MotifList> v = h.get(identifier);
			writeSequences(out, sequence, v,motifFlanking,  lrrMotifFlanking,  flankingRegion);
		}
		
		
		
		
		out.close();
		
	}
	
	
	
	private void writeSequences(BufferedWriter out, StringBuffer sequence, Vector<NLR_MotifList> v, int motifFlanking, int lrrMotifFlanking, int outsideFlanking)throws IOException{
		for(Enumeration<NLR_MotifList> myenum1 = v.elements(); myenum1.hasMoreElements();){
			NLR_MotifList list = myenum1.nextElement();
			
			int p1 = (int) list.getMotifs().firstElement().getDnaStart();
			int p2 = (int) list.getMotifs().firstElement().getDnaEnd();
			int p3 = (int) list.getMotifs().lastElement().getDnaStart();
			int p4 = (int) list.getMotifs().lastElement().getDnaEnd();
			
			int start = Math.max(0 , Math.min(Math.min(p1, p2), Math.min(p3, p4)) - outsideFlanking);
			int end = Math.min(sequence.length(), Math.max(Math.max(p1, p2), Math.max(p3, p4)) + outsideFlanking);
		
			int offset = start;
			
			
			int[] a = new int[end-start];
			for( int i = 0; i< outsideFlanking && i < a.length; i++){
				a[i]++;
			}
			for( int i = a[a.length-1]; i> a.length-outsideFlanking && i>=0; i--){
				a[i]++;
			}
			for(Enumeration<NLR_Motif> myenum2= list.getMotifs().elements(); myenum2.hasMoreElements();){
				NLR_Motif motif = myenum2.nextElement();
				int flanking = motifFlanking;
				if(isLRRMotif(motif)){
					flanking = lrrMotifFlanking;
				}
				for(int i = Math.max(0,  (int)motif.getDnaStart() - offset - flanking) ; i< Math.min(a.length, (int)motif.getDnaEnd() - offset + flanking); i++ ){
					a[i]++;
				}
			}
			int exonNumber = 0;
			int begin = 0;
			boolean writing = true;
			for( int i = 0; i< a.length; i++){
				if( a[i] >0){
					if(!writing){
						writing = true;
						begin = i;
					}
				}else{
					if(writing){
						writing = false;
						exonNumber++;
						BioSequence seq = new BioSequence(list.getNlrID() + "_exon" +exonNumber, sequence.substring(offset+begin, Math.min(offset+i, sequence.length())));
						seq.setDescription("scaffold="+list.getSequenceName() + " start="+(offset+begin)+ " end="+(offset+i) + " direction=forward");
						if(!list.getMotifs().firstElement().isForwardStrand()){
							seq.setDescription("scaffold="+list.getSequenceName() +  " start="+(offset+begin)+ " end="+(offset+i) + " direction=reverse");
							seq.setSequence(seq.getReverseComplementarySequence());
						}
						out.write(seq.getFastaString(100));
						//System.out.println(seq.getIdentifier() + "\t" + seq.getLength() + "\t" + seq.getDescription());
					}
				}
				
			}
			if(writing){
				exonNumber++;
				BioSequence seq = new BioSequence(list.getNlrID() + "_exon" +exonNumber, sequence.substring(offset+begin, Math.min(offset+a.length, sequence.length())));
				seq.setDescription("scaffold="+list.getSequenceName() + " start="+(offset+begin)+ " end="+(offset+a.length) + " direction=forward");
				if(!list.getMotifs().firstElement().isForwardStrand()){
					seq.setDescription("scaffold="+list.getSequenceName() + " start="+(offset+begin)+ " end="+(offset+a.length) + " direction=reverse");
					seq.setSequence(seq.getReverseComplementarySequence());
				}
				out.write(seq.getFastaString(100));
			}
			
		}
		
		
	}
		
	

	
	public void writeBaitCoverage(File blastFile, File outputFile)throws IOException{
		Hashtable<String, Vector<long[]>> baitRegions = this.getBaitCoverage(blastFile, 100, 80);
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		for(Enumeration<NLR_MotifList> myenum = this.nlrs.elements(); myenum.hasMoreElements();){
			NLR_MotifList list = myenum.nextElement();
			Vector<long[]> regions = new Vector<long[]>();
			for(Enumeration<NLR_Motif> myenum2 = list.getMotifs().elements(); myenum2.hasMoreElements();){
				NLR_Motif motif = myenum2.nextElement();
				long[] a = {Math.min(motif.getDnaStart(), motif.getDnaEnd()), Math.max(motif.getDnaStart(), motif.getDnaEnd())};
				regions = addRegion(regions,a);
			}
			long numPositions = 0;
			long coveredPositions = 0;
			
			Vector<long[]>w = baitRegions.get(list.getSequenceName());
			
			for(Enumeration<long[]> myenum1 = regions.elements(); myenum1.hasMoreElements();){
				long[] a = myenum1.nextElement();
				numPositions = numPositions + a[1] - a[0] +1;
				if(w != null){
					for(Enumeration<long[]> myenum2 = w.elements(); myenum2.hasMoreElements();){
						long[] b = myenum2.nextElement();
						if( Math.max(b[1],a[1]) - Math.min(a[0],b[0]) < (a[1]-a[0]) + (b[1]-b[0]) ){//they overlap
							coveredPositions = coveredPositions + Math.min(a[1],b[1])+1 - Math.max(a[0], b[0]);
						}
					}
				}
			}
			String complete = "complete";
			if(!list.isCompleteNLR()){
				complete = "partial";
			}
			out.write(list.getSequenceName() + "\t" + list.getNlrID() + "\t"+complete+ "\t" + numPositions + "\t" + coveredPositions + "\t" + (coveredPositions *100/numPositions));
			out.newLine();
			
		}
		
		
		
		out.close();
		
		
		
	}
	
	
	public Hashtable<String, Vector<Long>> generateMotifDistanceMatrix()throws IOException{
		Hashtable<String, Vector<Long>> h = new Hashtable<String, Vector<Long>>();
		
		for(Enumeration<NLR_MotifList> myenum = this.nlrs.elements(); myenum.hasMoreElements();){
			
			NLR_MotifList list = myenum.nextElement();
			
			if(!list.isCompleteNLR()){
				continue;
			}
			
			int lastElement = -1;
			long lastPosition = -1;
			
			for(Enumeration<NLR_Motif> myenum2 = list.motifs.elements(); myenum2.hasMoreElements();){
				NLR_Motif motif = myenum2.nextElement();
				if( lastElement !=-1){
					long distance = Math.abs(motif.getDnaStart()-lastPosition);
					String a = lastElement+"_"+ motif.getMotifNumber();
					if( distance > 5000){
						System.out.println(list.getSequenceName() + "\t"+a + "\t" + list.getNlrID() + "\t" + distance );
					}
					
					Vector<Long> v = new Vector<Long>();
					if( h.containsKey(a)){
						v = h.get(a);
					}
					v.add(new Long(distance));
					h.put(a, v);
				}
				lastElement = motif.getMotifNumber();
				lastPosition = motif.getDnaStart();
				
			}
			
			
		}
		
		
		
		for(Enumeration<String> myenum = h.keys(); myenum.hasMoreElements();){
			String a = myenum.nextElement();
			
			if( h.get(a).size()>100){
				BufferedWriter out = new BufferedWriter(new FileWriter("/Users/steuernb/Documents/projects/ChineseSpring/IWGSC_assembly/intronTests/distance"+a));
	
				
				
				long min = 100000;
				long max = 0;
				long sum = 0;
				for( Enumeration<Long> myenum2 = h.get(a).elements(); myenum2.hasMoreElements();){
					long l = myenum2.nextElement();
					
						out.write(l+"");
						out.newLine();
						
					if( l<min){
						min = l;
					}
					if(l> max){
						max = l;
					}
					sum = sum + l;
				}
				out.close();
			}
			//	System.out.println(a+ "\t" + h.get(a).size() + "\t" + min + "\t" + max + "\t" + (sum/h.get(a).size()) + "\t" + Math.round(sd(h.get(a))));
			
			
		}
		
		System.out.println(h.size());
		
		
		return h;
	}
	
	
	
	
	
	
	public void writeGffForExtractedSequence(File inputFasta, File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("##gff-version 2");out.newLine();
		out.write("##source-version NLR-Annotator "+this.getVersion());out.newLine();
		out.write("##date 2016-07-05 15:44");out.newLine();
		out.write("##Type DNA");out.newLine();
		out.write("#seqname	source	feature	start	end	score	strand	frame	attribute");out.newLine();
		
		FastaReader fastaReader = new FastaReader(inputFasta);
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			//System.out.println(seq.getIdentifier() + "\t" +seq.getDescription());
			String[] split = seq.getDescription().trim().split("\\s");
			//scaffold18474	NLR_Annotator_old1	NBSLRR	881278	884500	.	-	.	name=scaffold18474_nlr_1
			//scaffold=scaffold179 start=5199 end=6131 direction=forward
			String scaffold = split[0].split("=")[1];
			int start = Math.min(Integer.parseInt(split[1].split("=")[1]), Integer.parseInt(split[2].split("=")[1]) );
			int end  =  Math.max(Integer.parseInt(split[1].split("=")[1]), Integer.parseInt(split[2].split("=")[1]) );
			String direction = "+";
			if( !split[3].split("=")[1].equalsIgnoreCase("forward")){
				direction = "-";
			}
			out.write(scaffold + "\tNLR_Annotator\tNBSLRR\t" + (start+1) + "\t" + end + "\t.\t" + direction + "\t.\tname="+seq.getIdentifier()  );
			out.newLine();
		}
		fastaReader.close();
		out.close();
		
		
	}
	
	
	
	public void writeNbarcMotifSequences(File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
	
		for(Enumeration<NLR_MotifList> myenum = this.nlrs.elements(); myenum.hasMoreElements();){
			NLR_MotifList list = myenum.nextElement();
			String motifs = "";
			String motifNumbers = "";
			
			for(Enumeration<NLR_Motif> myenum2 = list.getMotifs().elements(); myenum2.hasMoreElements();){
				NLR_Motif motif = myenum2.nextElement();
				if(this.isNbarcMotif(motif)){
					motifNumbers = motifNumbers + ","+motif.getMotifNumber();
					motifs = motifs + motif.getSequence();
					
					
				}
			}
			
			
		}
		
		
		out.close();
		
	}
	
	
	public static void main(String[] args) {
		try {
			File inputFile = new File("/Volumes/group-scratch/Matthew-Hartley/steuernb/Eucalyptus/nlr_annotation/Ccitriodorav2_intermediate.xml");
			File genomeSequence = new File("/Volumes/group-scratch/Matthew-Hartley/steuernb/Eucalyptus/Ccitriodora.prelim.v2.1.annot/CCV_018_PMs_scaffs.idFix.fa");
			File outputFile = new File("/Volumes/group-scratch/Matthew-Hartley/steuernb/Eucalyptus/pipeline/step1/cc.nlr.fasta");
			NLR_Annotator annotator = new NLR_Annotator(inputFile, 500, 2500, 10000);
			annotator.writeReportTxt(new File("/Volumes/group-scratch/Matthew-Hartley/steuernb/Eucalyptus/pipeline/step1/cc.nlr.txt"));
			annotator.writeNLRLoci(genomeSequence, outputFile, 250, 250, 250, false);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void mainCLI(String[] args){
		CLI cli = new CLI();
		
		
		
		String helpString =		"NLR-Annotator " + version + "\n"+
								
								"-i <inputFasta>\t\t\tInput file in xml format. This is what comes out from NLR-Parser -c of the chopped inputSequence\n"+
								"-o <output.txt>\t\t\toutput file in tabular format\n"+
								"-g <output.gff>\t\t\toutput file in gff format\n"+
								"-b <output.bed>\t\t\toutput file in bed format\n"+
								"-m <output.motifs.bed>\t\t\toutput file of the motifs in bed format\n"+
								"-a <output.nbarkMotifAlignment.fasta>\t\toutput file of the nb-arc motifs as multiple alignment.\n"+
								"-f <genome.fasta> <output.nlr.fasta> <(int) flanking>\twrite fasta of nlr loci.\n"+
								"-distanceWithinMotifCombination <int>\t(default:500)\n"+
								"-distanceForElongating <int>\t\t(default:2500)\n" + 
								"-distanceBetweenMotifCombinations <int>\t(default:50000)";
			
			
		try{
			
			
			cli.parseOptions(args);
			if(!cli.hasOption("i") ){
				throw new CLIParseException("Missing options. -i is mandatory");
			}
			
			File inputFile = new File( cli.getArg("i") );
			
			
			

			int distanceWithinMotifCombination=500;
			int distanceForElongating=2500;
			int distanceBetweenMotifCombinations=10000;
			
			if(cli.hasOption("distanceWithinMotifCombination")){
				distanceWithinMotifCombination = Integer.parseInt(cli.getArg("distanceWithinMotifCombination"));
			}
			if(cli.hasOption("distanceForElongating")){
				distanceForElongating = Integer.parseInt(cli.getArg("distanceForElongating"));
			}
			if(cli.hasOption("distanceBetweenMotifCombinations")){
				distanceBetweenMotifCombinations = Integer.parseInt(cli.getArg("distanceBetweenMotifCombinations"));
			}
			
			
			
			NLR_Annotator annotator = new NLR_Annotator(inputFile, distanceWithinMotifCombination, distanceForElongating, distanceBetweenMotifCombinations);
			
			if( cli.hasOption("o")){
				annotator.writeReportTxt(new File(cli.getArg("o")));
			}
			
			if( cli.hasOption("g") ){
				annotator.writeNlrGFF(new File(cli.getArg("g")), false);
			}
			
			if(cli.hasOption("b")){
				annotator.writeNlrBED(new File(cli.getArg("b")));
			}
			
			if(cli.hasOption("m")){
				annotator.writeMotifBED(annotator.motifs, new File(cli.getArg("m")), false);
			}
			
			if( cli.hasOption("a")){
				annotator.writeNbarcMultipleAlignmentFasta(new File(cli.getArg("a")));
			}
			
			if( cli.hasOption("f")){
				Vector<String> s = cli.getArgs("f");
				annotator.writeNLRLoci(new File(s.get(0)), new File(s.get(1)), Integer.parseInt(s.get(2)));
			}
			
			
		}catch(CLIParseException e){
			e.printStackTrace();
			
			System.err.println(helpString);
			
			
			
		}
		catch(IOException e){
			e.printStackTrace();
			System.err.println(helpString);
		}
		catch(ParserConfigurationException e){
			e.printStackTrace();
		}
		
		catch(SAXException e){
			e.printStackTrace();
		}
		
		
		
	}
	
	
	
}
