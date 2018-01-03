package nlr_parser;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.w3c.dom.Document;
import org.w3c.dom.Element;

public class MastMotifHit implements Comparable<MastMotifHit>{

	String sequence_name;
	int pos;
	int gap;
	String motif;
	double pvalue;
	String match;
	
	boolean forwardStrand;
	int frame;
	
	public MastMotifHit(String sequence_name, String hitString){
		Pattern p = Pattern.compile("<hit\\s+pos=\"(\\d+)\"\\s+gap=\"(\\d+)\"\\s+motif=\"(\\w+)\"\\s+pvalue=\"([\\w-\\.]+)\"\\s+match=\"([\\s+]+)\"/>");
		Matcher m = p.matcher(hitString);
		m.find();
		this.pos=Integer.parseInt(m.group(1));
		this.gap=Integer.parseInt(m.group(2));
		this.motif = m.group(3);
		this.pvalue = Double.parseDouble(m.group(4));
		this.match = m.group(5);
		this.sequence_name = sequence_name;
		
		String s= this.sequence_name.split("_frame")[1];
		forwardStrand = true;
		if(s.substring(0, 1).equalsIgnoreCase("-")){
			forwardStrand = false;
		}
		this.frame = Integer.parseInt(s.substring(1,2));
		
			
	}

	
	
	public MastMotifHit(String sequence_name, int pos, int gap, String motif,double pvalue, String match, boolean forwardStrand, int frame) {
		
		this.sequence_name = sequence_name;
		this.pos = pos;
		this.gap = gap;
		this.motif = motif;
		this.pvalue = pvalue;
		this.match = match;
		this.forwardStrand = forwardStrand;
		this.frame = frame;
		
	}


	public MastMotifHit(Element xmlElement){
		boolean forwardStrand=true;
		if(xmlElement.getAttribute("strand").equalsIgnoreCase("-")){
			forwardStrand = false;
		}
		
		this.sequence_name = xmlElement.getAttribute("sequence_name"); 
		this.pos = Integer.parseInt(xmlElement.getAttribute("pos"));
		this.gap = Integer.parseInt(xmlElement.getAttribute("gap"));
		this.motif = xmlElement.getAttribute("motif");
		this.pvalue = Double.parseDouble(xmlElement.getAttribute("pvalue"));
		this.match = xmlElement.getAttribute("match");
		this.forwardStrand = forwardStrand;
		this.frame = Integer.parseInt(xmlElement.getAttribute("frame") );
									
		//System.out.println(this.sequence_name + "\t" + this.pos);
	}
	

	public String toString(){
		return "<hit pos=\""+this.pos+"\" gap=\""+this.gap+"\" motif=\""+this.motif+"\" pvalue=\""+this.pvalue+"\" match=\""+this.match+"\"/>";
	}
	
	public String getSequence_name() {
		return sequence_name;
	}

	public int getPos() {
		return pos;
	}

	public int getDNAPos(){
		return (this.pos -1)*3 + this.frame;
	}
	
	public int getGap() {
		return gap;
	}

	public String getMotif() {
		return motif;
	}

	public double getPvalue() {
		return pvalue;
	}

	public String getMatch() {
		return match;
	}
	
	public boolean isForwardStrand(){
		return this.forwardStrand;
	}
	
	public String getStrand(){
		String strand = "forward";
		if(!this.isForwardStrand()){
			strand = "reverse";
		}
		return strand;
	}
	
	public int compareTo(MastMotifHit m2){
		if( this.getPos() < m2.getPos()){
			return -1;
		}
		else if( this.getPos() > m2.getPos()){
			return 1;
		}else{
			return 0;
		}
		
		
	}
	
	public int getFrame(){
		return this.frame;
	}
	
	public void setFrame(int frame){
		this.frame = frame;
	}
	
	public Element getXMLElement(Document dom){
		Element hitElement = dom.createElement("MastMotifHit");
		
		hitElement.setAttribute("sequence_name", this.getSequence_name());
		hitElement.setAttribute("pos", this.getPos()+"");
		hitElement.setAttribute("gap", this.getGap()+"");
		hitElement.setAttribute("motif", this.getMotif());
		hitElement.setAttribute("pvalue", this.getPvalue()+"");
		hitElement.setAttribute("match", this.getMatch());
		String strand = "+";
		if(!this.isForwardStrand()){
			strand = "-";
		}
		hitElement.setAttribute("strand", strand);
		hitElement.setAttribute("frame", this.getFrame()+"");
		
		return hitElement;
		
		
	}
	
	
	/**
	 * Deliver the color code for motifs as we know it from the mast html output. 
	 * 
	 * @param motifNumber
	 * @return
	 */
	public  int[] getRGBMotifColor(){
		
		int motifNumber = Integer.parseInt(this.motif.split("_")[1]);
		
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
	
	
}
