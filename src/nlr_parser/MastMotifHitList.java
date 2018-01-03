package nlr_parser;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import support.BioSequence;

public class MastMotifHitList {
	Vector<MastMotifHit> motivlist;
	int aAsequenceLength;
	boolean sorted;
	
	
	String dnaSequence;
	String aaSequence;
	
	
	
	
	public MastMotifHitList(Vector<MastMotifHit> v){
		this.motivlist = v;
		this.sort();
		
	}
	public MastMotifHitList(){
		this.motivlist = new Vector<MastMotifHit>();
		sorted = true;
	}
	
	public MastMotifHitList(Element xmlElement){
		
			Vector<MastMotifHit> hits = new Vector<MastMotifHit>();
			
			Element mastMotifHits = (Element) xmlElement.getElementsByTagName("MastMotifHits").item(0);
			
			NodeList hitNodeList = mastMotifHits.getElementsByTagName("MastMotifHit");
			for( int j = 0; j< hitNodeList.getLength(); j++){
				Element hitElement = (Element)hitNodeList.item(j);
				MastMotifHit hit = new MastMotifHit(hitElement);
				hits.add(hit);
									
			}
			
			this.motivlist = hits;
			
		
			Element sequenceElement = (Element) xmlElement.getElementsByTagName("Sequence").item(0);
			if(sequenceElement != null){
				if( sequenceElement.hasAttribute("DNA")){
					this.dnaSequence = (sequenceElement.getAttribute("DNA"));
				}
				if( sequenceElement.hasAttribute("AminoAcid")){
					this.aaSequence = (sequenceElement.getAttribute("AminoAcid"));
				}
			}
	}
	
	
	
	public void addMotif(MastMotifHit hit){
		motivlist.add(hit);
		sorted = false;
	}
	
	public void sort(){
		Collections.sort(motivlist);
		sorted = true;
	}
	

	
	
	
	
	
	public void setSequenceLength(int length){
		this.aAsequenceLength = length;
	}
	public int getAASequenceLength(){
		return this.aAsequenceLength;
	}
	
	public void addAASequence(String aaSequence){
		this.aaSequence = aaSequence;
	}
	public void addDnaSequence(String dnaSequence){
		this.dnaSequence = dnaSequence;
	}
	public String getDNASequence(){
		return this.dnaSequence;
	}
	
	public static MastMotifHitList mergeLists( MastMotifHitList list1 , MastMotifHitList list2){
		MastMotifHitList list = new MastMotifHitList();
		
		for(Enumeration<MastMotifHit> myenum = list1.getMotifs().elements(); myenum.hasMoreElements();){
			list.addMotif(myenum.nextElement());
		}
		
		for(Enumeration<MastMotifHit> myenum = list2.getMotifs().elements(); myenum.hasMoreElements();){
			list.addMotif(myenum.nextElement());
		}
		
		list.sort();
		
		return list;
	}
	
	
	
	
	
	public  boolean isNibbler(){
		
		
		if(!sorted){
			this.sort();
		}
				
		
		Vector<int[]> motivCombinations = new Vector<int[]>();
		int[] b = {17, 16}; motivCombinations.add(b);
		b = new int[3]; b[0] = 1; b[1] = 6; b[2] = 4; motivCombinations.add(b);
		b = new int[3]; b[0] = 6; b[1] = 4; b[2] = 5; motivCombinations.add(b);
		b = new int[3]; b[0] = 4; b[1] = 5; b[2] = 10;motivCombinations.add(b);
		b = new int[3]; b[0] = 5; b[1] = 10; b[2] = 3;motivCombinations.add(b);
		b = new int[3]; b[0] = 10; b[1] = 3; b[2] = 12; motivCombinations.add(b);
		b = new int[3]; b[0] = 3; b[1] = 12; b[2] = 2; motivCombinations.add(b);
		b = new int[3]; b[0] = 12; b[1] = 2; b[2] = 8; motivCombinations.add(b);
		b = new int[3]; b[0] = 2; b[1] = 8; b[2] = 7; motivCombinations.add(b);
		b = new int[3]; b[0] = 8; b[1] = 7; b[2] = 9; motivCombinations.add(b);
		b = new int[3]; b[0] = 7; b[1] = 9; b[2] = 11;motivCombinations.add(b);
		b = new int[2]; b[0] = 9; b[1] = 11;  motivCombinations.add(b);
		b = new int[2]; b[0] = 18; b[1] = 15; motivCombinations.add(b);
		b = new int[2]; b[0] = 15; b[1] = 13; motivCombinations.add(b);
		b = new int[2]; b[0] = 13; b[1] = 1; motivCombinations.add(b);
		b = new int[3]; b[0] = 1; b[1] = 4; b[2]=5 ; motivCombinations.add(b);
		
		
		
		
		
		int[] a = new int[motivlist.size()];
		for( int i = 0; i< a.length; i++){
			a[i] = Integer.parseInt(motivlist.get(i).getMotif().split("_")[1]);
		}
		
		
		for( Enumeration<int[]> myenum = motivCombinations.elements(); myenum.hasMoreElements();){
			b = myenum.nextElement();
			for (int i = 0; i<= a.length-b.length; i++ ){
				if( a[i] == b[0]){
					boolean found = true;
					
					for(int j = 1; j<b.length; j++){
						if(b[j] != a[i+j]){
							found = false;
						}
					}
					
					if(found){
						return true;
					}
					
				}
				
			}
			
			
		}	
		return false;
	}
	
	
	/**
	 * This checks if the NBLRR is complete based on the existence of certain motifs. Note that this method not specifically checks if it is a NBLRR. 
	
	 * @return
	 */
	public  boolean isComplete(){
		if(!sorted){
			sort();
		}
		
		int[] a = new int[motivlist.size()];
		for( int i = 0; i< a.length; i++){
			a[i] = Integer.parseInt(motivlist.get(i).getMotif().split("_")[1]);
		}
		
		boolean hasNTerminal = false;
		boolean hasNBARC = false;
		//boolean hasLRR = false;
		
		for(int i = 0; i< a.length; i++){
			int m = a[i];
			if(m ==16 || m == 18){
				hasNTerminal = true;
			}
			if(m == 3 && hasNTerminal){
				hasNBARC=true;
			}
			if(m == 11 && hasNBARC){
				//hasLRR = true;
				return true;
			}
			
		}
		return false;
		
	}
	
	public  int[] getMotifList(){
		
		
		int[] a = new int[motivlist.size()];
		for( int i = 0; i< a.length; i++){
			a[i] = Integer.parseInt(motivlist.get(i).getMotif().split("_")[1]);
		}
		return a;
	}	
	
	public  String getNBLRR_Class(){
		if(!sorted){
			sort();
		}
		
		int[] a = new int[motivlist.size()];
		for( int i = 0; i< a.length; i++){
			a[i] = Integer.parseInt(motivlist.get(i).getMotif().split("_")[1]);
		}
		
		
		
		
		for( int i = 0; i< a.length; i++){
			
			
			try{
				if(a[i] == 2 || (a[i]==1 && a[i+1] ==6  ) ||  (a[i] == 17 && a[i+1] ==16) ){
					return "CNL";
				}
			}catch (ArrayIndexOutOfBoundsException e){}	
			
		}
		
		
		for( int i = 0; i< a.length; i++){
			try{
				if(a[i] ==18 && a[i+1] == 15 && a[i+2] ==13){
					return "TNL";
				}
			}catch (ArrayIndexOutOfBoundsException e){}
			
				
			try{
				if( a[i] ==1 && a[i+1]!=6   && (a[i+1] ==4   || a[i+2] ==4)    ){
					return "TNL";
				}
			}catch (ArrayIndexOutOfBoundsException e){}	
			
		}
		
		
		
		
		return "N/A";
		
	}
	
	/**
	 * Get the frame of the first motif in the list. 
	 * 
	 * @return
	 */
	public int getFirstMotifFrame(){
		if(!sorted){
			sort();
		}
		return this.motivlist.get(0).getFrame();
	}
	
	/**
	 * Get the frame of hte last motif in the list.
	 * 
	 * @return
	 */
	public int getLastMotifFrame(){
		if(!sorted){
			sort();
		}
		return this.motivlist.get(this.motivlist.size()-1).getFrame();
	}
	
	/**
	 * gives the name of the sequence and the position the first motif is found. Sequence-name and position are separated with " : "
	 * 
	 * 
	 * 		A vector of MastMotifHits
	 * @return
	 */
	public  int getAAStart(){
		if(!sorted){
			sort();
		}
		
		MastMotifHit hit = this.motivlist.get(0);
		return  hit.getPos();
	}
	
	public  int getDNAStart(){
		if(!sorted){
			sort();
		}
		
		MastMotifHit hit = this.motivlist.get(0);
		return hit.getDNAPos();
	}
	
	
	
	public  int getAAEnd(){
		if(!sorted){
			sort();
		}
		MastMotifHit hit = motivlist.get(motivlist.size()-1);
		return  (hit.getPos() + hit.getMatch().length());
	}
	
	
	public  int getDNAEnd(){
		if(!sorted){
			sort();
		}
		MastMotifHit hit = motivlist.get(motivlist.size()-1);
		return  (hit.getDNAPos() + (hit.getMatch().length())*3  );
	}
	
	
	/**
	 * gives the sequence starting from 10 aa before the Motif 1 and 200 aa after beginning of Motif 1. 
	 * 
	 * 
	 * @return
	 */
	public  String getNBARC_Sequence(){
		if(!sorted){
			sort();
		}
		for( Enumeration<MastMotifHit> myenum = motivlist.elements(); myenum.hasMoreElements();){
			MastMotifHit hit = myenum.nextElement();
			if(this.aaSequence == null && this.dnaSequence == null){
				return "N/A";
			}
			if(Integer.parseInt(hit.getMotif().split("_")[1]) == 1 ){
				if(this.aaSequence == null){
					this.aaSequence = new BioSequence("", this.dnaSequence).translate2Protein()[hit.frame].getSequence();
					
				}
				int start = Math.max(0, hit.getPos()-11);
				int end = Math.min(hit.getPos() + 200, this.aaSequence.length()-1);
				
				
				return this.aaSequence.substring(start, end);
			}
		}
		return "N/A";
	}

	public  boolean hasPotentiallyTwoNBLRRs(){
		if(!sorted){
			sort();
		}
		
		
		
		boolean foundLRR = false;
		
		int[] a = getMotifList();
		
		for( int i = 0; i< a.length; i++){
			
			if(a[i] == 9   ||  a[i] == 11 || a[i] ==19){
				foundLRR=true;
			}
			if(a[i] ==1 && foundLRR){
				return true;
			}
		}
		
		
		return false;
	}
	
	public String getMotifListString(){
		
		if(!sorted){
			sort();
		}
		
		int[] a = getMotifList();
		String s= "";
		for( int i = 0; i< a.length; i++){
			s = s + "," + a[i];
		}
		return s.substring(1);
		
	}
	
	public int getSize(){
		return this.motivlist.size();
	}
	public Vector<MastMotifHit> getMotifs(){
		return this.motivlist;
	}
	
	
	
	
	public int getForwardStrandPosition(MastMotifHit hit){
		if( this.dnaSequence == null){
			return -1;
		}
		
		if(!hit.forwardStrand){
			return this.dnaSequence.length()-hit.getDNAPos();
		}
		
		return hit.getDNAPos();
	}
	
	
	
	
	
	public String getMotifSequence(MastMotifHit hit){
		
		String aaseq = this.aaSequence;
		if( this.aaSequence == null){
			if(this.dnaSequence==null){
				return null;
			}
			
			BioSequence[] a = new BioSequence("", this.dnaSequence).translate2Protein();
			
			
			if(hit.isForwardStrand()){
				
				aaseq = a[hit.getFrame()].getSequence();
			}else{
				
				aaseq = a[hit.getFrame() + 3].getSequence();
				
			}
			
		}
		
		return aaseq.substring(hit.getPos(), hit.getPos() + hit.getMatch().length());
		
		
		
	}
	
	public static Hashtable<String, MastMotifHitList> importFromXML(File inputFile)throws IOException,ParserConfigurationException, SAXException {
		Hashtable<String, MastMotifHitList> h = new Hashtable<String, MastMotifHitList>();
		
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		Document dom = db.parse(inputFile);
		
		
		Element rootElement = dom.getDocumentElement();
		
		
		NodeList nodeList = rootElement.getElementsByTagName("MastMotifHitList");
		
		for( int i = 0; i< nodeList.getLength(); i++){
			Element listElement = (Element)nodeList.item(i);
			String name = listElement.getAttribute("name");
			
			MastMotifHitList list = new MastMotifHitList(listElement);
			
			h.put(name, list);
		}
		
		
		
		
		return h;
	}
	
	/**
	 * 
	 * 
	 * @param key
	 * 			usually MastMotifHitLists are stored in a hashtable this is the key of that hashtable, usually corresponding to the sequencename plus the frame.
	 * @param dom
	 * 			The Document object used for creating the XML.
	 * @return
	 */
	public Element getXMLElement(String key, Document dom){
		Element listElement = dom.createElement("MastMotifHitList");
		listElement.setAttribute("name", key);
		
		if( this.aaSequence!= null || this.dnaSequence != null){
			Element sequenceElement = dom.createElement("Sequence");
			if( this.aaSequence!= null){
				sequenceElement.setAttribute("AminoAcid", this.aaSequence);
			}
			if( this.dnaSequence!= null){
				sequenceElement.setAttribute("DNA", this.dnaSequence);
			}
			listElement.appendChild(sequenceElement);
		}
		
		if(!this.sorted){
			this.sort();
		}
		
		Element mastMotifHits = dom.createElement("MastMotifHits");
		listElement.appendChild(mastMotifHits);
		
		for(Enumeration<MastMotifHit> myenum2 = this.getMotifs().elements(); myenum2.hasMoreElements(); ){
			MastMotifHit hit = myenum2.nextElement();
			Element hitElement = hit.getXMLElement(dom);
			mastMotifHits.appendChild(hitElement);
		}
		return listElement;
	}
	
	
	public static void exportToXML(Hashtable<String, MastMotifHitList> h, File outputFile)throws IOException, ParserConfigurationException, TransformerConfigurationException, TransformerException{
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		Document dom = db.newDocument();
		
		Element rootElement = dom.createElement("NLRParser");
		dom.appendChild(rootElement);
		
		for(Enumeration<String> myenum1 = h.keys(); myenum1.hasMoreElements();){
			String key = myenum1.nextElement();
			MastMotifHitList list = h.get(key);
			
			Element listElement = list.getXMLElement(key, dom);
			
			rootElement.appendChild(listElement);
		}
		
		
		
		
		DOMSource source = new DOMSource(rootElement) ;
		StreamResult result = new StreamResult(outputFile);
		Transformer transformer = TransformerFactory.newInstance().newTransformer();
		if (dom.getDoctype() != null) {
		    String systemValue = (new File (dom.getDoctype().getSystemId())).getName();
		    transformer.setOutputProperty(OutputKeys.DOCTYPE_SYSTEM, systemValue);
		}
		transformer.setOutputProperty(OutputKeys.INDENT, "yes"); //without this there is no newline after elements.
		
		transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");  //not sure what happens here but this makes the nice whitespace hierarchy in the output txt.
		transformer.transform(source, result);
		
	}
	
	
	
	
	
	
	
	
	
}
