package nlr_annotator;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

import support.BioSequence;

public class NLR_Motif implements Comparable<NLR_Motif>{

	String sequenceName;
	int motifNumber;
	long dnaStart;
	long dnaEnd;
	double pvalue;
	boolean forwardStrand;
	int frame;
	
	String sequence;
	
	
	
	public NLR_Motif(String sequenceName, int motifNumber, long dnaStart,
			long dnaEnd, double pvalue, boolean forwardStrand, int frame) {
		super();
		this.sequenceName = sequenceName;
		this.motifNumber = motifNumber;
		this.dnaStart = dnaStart;
		this.dnaEnd = dnaEnd;
		this.pvalue = pvalue;
		this.forwardStrand = forwardStrand;
		this.frame = frame;
	}



	public NLR_Motif(Element xmlElement){
		
		this.sequenceName = xmlElement.getAttribute("sequence_name");
		this.motifNumber = Integer.parseInt(xmlElement.getAttribute("motif_number"));
		this.dnaStart  = Integer.parseInt(xmlElement.getAttribute("start"));
		this.dnaEnd = Integer.parseInt(xmlElement.getAttribute("end"));
		this.pvalue = Double.parseDouble(xmlElement.getAttribute("pvalue"));
		this.frame = Integer.parseInt(xmlElement.getAttribute("frame"));
		
		this.forwardStrand = true;
		if(xmlElement.getAttribute("strand").equalsIgnoreCase("-")){
			forwardStrand = false;
		}
		
	}
	
	

	
	public int compareTo(NLR_Motif motif2){
		if( this.isForwardStrand() && !motif2.isForwardStrand()){
			return -1;
		}
		if( !this.isForwardStrand() && motif2.isForwardStrand()){
			return 1;
		}
		if( this.dnaStart < motif2.getDnaStart()){
			return -1;
		}else if( this.dnaStart > motif2.getDnaStart()){
			return 1;
		}else{
			return 0;
		}
	}





	public boolean isNbarcMotif(){
		int i = this.getMotifNumber();
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











	public String getSequenceName() {
		return sequenceName;
	}


	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}


	public int getMotifNumber() {
		return motifNumber;
	}



	public void setMotifNumber(int motifNumber) {
		this.motifNumber = motifNumber;
	}


	public long getDnaStart() {
		return dnaStart;
	}



	public void setDnaStart(long dnaStart) {
		this.dnaStart = dnaStart;
	}


	public long getDnaEnd() {
		return dnaEnd;
	}



	public void setDnaEnd(long dnaEnd) {
		this.dnaEnd = dnaEnd;
	}



	public double getPvalue() {
		return pvalue;
	}



	public void setPvalue(double pvalue) {
		this.pvalue = pvalue;
	}


	public boolean isForwardStrand() {
		return forwardStrand;
	}

	public String getStrand(){
		if(this.isForwardStrand()){
			return "+";
		}else{
			return "-";
		}
	}

	public void setForwardFrame(boolean forwardStrand) {
		this.forwardStrand = forwardStrand;
	}



	public int getFrame() {
		return frame;
	}



	public void setFrame(int frame) {
		this.frame = frame;
	}


	public String getRGBMotifColor(){
		int[] a = getRGBMotifColor(this.getMotifNumber());
		
		return  a[0] +"," +a[1] +","+a[2];
	}
	

	public Element getXMLElement(Document dom){
		Element hitElement = dom.createElement("NLR_Motif");
		
		hitElement.setAttribute("sequence_name", this.getSequenceName());
		hitElement.setAttribute("start", this.getDnaStart()+"");
		hitElement.setAttribute("end", this.getDnaEnd()+"");
		String strand = "+";
		if(!this.isForwardStrand()){
			strand = "-";
		}
		hitElement.setAttribute("strand", strand);
		hitElement.setAttribute("frame", this.getFrame()+"");
		hitElement.setAttribute("motif_number", this.getMotifNumber() + "");
		hitElement.setAttribute("pvalue", this.getPvalue()+"");
		
		return hitElement;
		
		
	}



	public String getBEDString(){
		
		return sequenceName +"\t" + this.getDnaStart() +"\t" + this.getDnaEnd() +"\tmotif_"+this.getMotifNumber()+"\t0\t"+ this.getStrand() +"\t" + this.getDnaStart() +"\t" + this.getDnaEnd() +"\t" + this.getRGBMotifColor();
		
	}
	public String getBlackMotif(){
		return sequenceName +"\t" + this.getDnaStart() +"\t" + this.getDnaEnd() +"\tSTOP\t0\t"+ this.getStrand() +"\t" + this.getDnaStart() +"\t" + this.getDnaEnd() +"\t0,0,0" ;
	}


	public String getGFFString(){
		//eqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute
		return sequenceName + "\tNLR_Annotator\tMast_Motif\t" + (this.getDnaStart()+1) +"\t" + (this.getDnaEnd()) +"\t" + this.getPvalue() + "\t" +this.getStrand() + "\t" + this.getFrame() +"\tname=" + this.getMotifNumber(); 
		
	}


	public String getSequence(){
		return this.sequence;
	}

	public void setSequence(String sequence){
		this.sequence = sequence;
	}
	
	public String getAASequence(){
		String seq = new BioSequence("", this.getSequence()).translate2Protein()[0].getSequence();
		if(this.getStrand()=="-"){
			seq = new BioSequence("", this.getSequence()).translate2Protein()[3].getSequence();
		
		
		}
		return seq;
	}



	/**
	 * Deliver the color code for motifs as we know it from the mast html output. 
	 * 
	 * @param motifNumber
	 * @return
	 */
	public static int[] getRGBMotifColor(int motifNumber){
		if( motifNumber ==1 ){ int[] a = {0,255,255}; return a;}
		if( motifNumber ==2 ){ int[] a = {0,0,255}; return a;}
		if( motifNumber ==3 ){ int[] a = {255,0,0}; return a;}
		if( motifNumber ==4 ){ int[] a = {255,0,255}; return a;}
		if( motifNumber ==5 ){ int[] a = {255,255,0}; return a;}
		if( motifNumber ==6 ){ int[] a = {0,255,0}; return a;}
		if( motifNumber ==7 ){ int[] a = {0,128,128}; return a;}
		if( motifNumber ==8 ){ int[] a = {68,68,68}; return a;}
		if( motifNumber ==9 ){ int[] a = {0,128,0}; return a;}
		if( motifNumber ==10 ){ int[] a = {192,192,192}; return a;}
		if( motifNumber ==11 ){ int[] a = {128,0,128}; return a;}
		if( motifNumber ==12 ){ int[] a = {128,128,0}; return a;}
		if( motifNumber ==13 ){ int[] a = {0,0,128}; return a;}
		if( motifNumber ==14 ){ int[] a = {128,0,0}; return a;}
		if( motifNumber ==15 ){ int[] a = {255,255,255}; return a;}
		if( motifNumber ==16 ){ int[] a = {0,255,255}; return a;}
		if( motifNumber ==17 ){ int[] a = {0,0,255}; return a;}
		if( motifNumber ==18 ){ int[] a = {255,0,0}; return a;}
		if( motifNumber ==19 ){ int[] a = {255,0,255}; return a;}
		if( motifNumber ==20 ){ int[] a = {255,255,0}; return a;}
		
		return new int[3];
	}
	
	
	/**
	 * return a ranking of motifs. I.e. TIR/CC before NB-ARC before linker before LRR. The return value of 0 means it should not be ranked.
	 * 
	 * 
	 */
	public static int getRanking(long motifNumber){
		
		if( motifNumber ==1 ) {return 2; }
		if( motifNumber ==2 ) {return 2; }
		if( motifNumber ==3 ) {return 2; }
		if( motifNumber ==4 ) {return 2; }
		if( motifNumber ==5 ) {return 2; }
		if( motifNumber ==6 ) {return 2; }
		if( motifNumber ==7 ) {return 3; }
		if( motifNumber ==8 ) {return 3; }
		if( motifNumber ==9 ) {return 4; }
		if( motifNumber ==10 ){return 2; }
		if( motifNumber ==11 ){return 4; }
		if( motifNumber ==12 ){return 2; }
		if( motifNumber ==13 ){return 1; }
		if( motifNumber ==14 ){return 0; }
		if( motifNumber ==15 ){return 1; }
		if( motifNumber ==16 ){return 1; }
		if( motifNumber ==17 ){return 1; }
		if( motifNumber ==18 ){return 1; }
		if( motifNumber ==19 ){return 4; }
		if( motifNumber ==20 ){return 0; }

		return 0;
		
		
	}
	
	
	
	
	
	
}
