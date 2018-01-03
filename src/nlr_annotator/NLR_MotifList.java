package nlr_annotator;

import java.util.Collections;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Vector;

public class NLR_MotifList {

	String sequenceName;
	String id;
	Vector<NLR_Motif> motifs;
	
	
	public NLR_MotifList(String sequenceName, String id){
		this.sequenceName = sequenceName;
		this.id = id;
		this.motifs = new Vector<NLR_Motif>();
	}
	
	public NLR_MotifList(String sequenceName, String id, Vector<NLR_Motif> v){
		this.sequenceName = sequenceName;
		this.id = id;
		this.motifs = v;
	}
	
	
	public void addMotif(NLR_Motif motif){
		this.motifs.add(motif);
		Collections.sort(motifs);
		if( !motifs.firstElement().isForwardStrand()){
			Collections.reverse(motifs);
		}
	}
	
	
	
	public Vector<NLR_Motif> getMotifs(){
		return this.motifs;
	}
	
	public String getSequenceName(){
		return sequenceName;
	}
	
	public String getNlrID(){
		return  id;
	}
	
	/**
	 * gives a 2 dimensional array with [number of motifs][2]. first row is motif number, second row is DNA start.
	 * 
	 * @return
	 */
	public long[][] getMotifListArary(){
		long[][] a = new long[motifs.size()][2];
		for( int i = 0; i< a.length; i++){
			a[i][0] = motifs.get(i).getMotifNumber();
			a[i][1] = motifs.get(i).getDnaStart();
		}
		return a;
	}
	
	public String getMotifList(){
		String s = "";
		for(Enumeration<NLR_Motif> myenum = this.motifs.elements(); myenum.hasMoreElements();){
			s = s + "," + myenum.nextElement().getMotifNumber();
		}
		return s.substring(1);
	}
	
	public String getBED(){
		String s = "";
		for(Enumeration<NLR_Motif> myenum = this.motifs.elements(); myenum.hasMoreElements();){
			s = s +myenum.nextElement().getBEDString()+"\n";
		}
		return s;
	}
	
	
	public void removeFirstMotif(){
		this.motifs.remove(0);
	}
	public void removeLastMotif(){
		this.motifs.remove(motifs.size()-1);
	}
	
	
	public  boolean isCompleteNLR(){
		
		String s= this.getMotifList();
		
		if( s.contains(",1,") || s.startsWith("1,")){
			if( s.contains(",9") || s.contains(",11") || s.contains(",19")){
				return true;
			}
			
		}
		
		
		return false;
	}
	
	
	public void removeRedundantMotifs(){
		HashSet<String> h = new HashSet<String>();
		Vector<NLR_Motif> v = new Vector<NLR_Motif>();
		for(Enumeration<NLR_Motif> myenum = this.motifs.elements(); myenum.hasMoreElements();){
			NLR_Motif motif = myenum.nextElement();
			String s = motif.getMotifNumber() + "\t" + motif.getDnaStart() + "\t" +motif.getDnaEnd();
			if(!h.contains(s)){
				v.add(motif);
			}
			h.add(s);
		}
		
		this.motifs = v;
	}
	
	
	public void removeInconsistentTirs(){
		if(this.getMotifList().contains(",6,") || this.getMotifList().contains(",2,")){
			while(motifs.firstElement().getMotifNumber() == 18 || motifs.firstElement().getMotifNumber() == 15 || motifs.firstElement().getMotifNumber() == 13){
				this.motifs.remove(0);
			}
			
			
		}
		
	}
	
	
	public int getLength(){
		return this.motifs.size();
	}
	
	
	
	public boolean hasStopCodon(){
		
		for(Enumeration<NLR_Motif > myenum = this.getMotifs().elements(); myenum.hasMoreElements();){
			NLR_Motif motif = myenum.nextElement();
			if(motif.getAASequence().contains("*")){
				return true;
			}
		}
		
		return false;
	}
}
