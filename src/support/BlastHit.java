package support;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

public class BlastHit implements Comparable<BlastHit> {
	String queryID;
	String hitID;
	
	int hitRank;
	
	Vector<BlastHSP> hsps;
	
	
	public BlastHit(String queryID, String hitID) {
		super();
		this.queryID = queryID;
		this.hitID = hitID;
		this.hsps = new Vector<BlastHSP>();
	}

	/*
	public BlastHit(String xml, String queryID, String queryDescription, int queryLength){
			
		
		Pattern p = Pattern.compile("<Hit_id>(.*)</Hit_id>");
		Matcher m = p.matcher(xml);
		if(m.find()){
			this.hitID = m.group();
		}
		p= Pattern.compile("<Hit_len>[0-9]+</Hit_len>");
		m=p.matcher(xml);
		if(m.find()){
			this.hitLength = Integer.parseInt(m.group());
		}
		this.queryID = queryID;
		this.queryDescription = queryDescription;
		this.queryLength = queryLength;
		this.hitDescription = ("N/A");
		
		p= Pattern.compile("<Hit_num>[0-9]+</Hit_num>");
		m=p.matcher(xml);
		if(m.find()){
			this.hitNum = Integer.parseInt(m.group());
		}
		
		p = Pattern.compile("<HSP>.*<HSP>");
		m = p.matcher(xml);
		while(m.find()){
			BlastHSP hsp = new BlastHSP(m.group());
			this.addHSP(hsp);
		}
	}
	
	
	*/
	
	
	
	
	public boolean addHSP(BlastHSP hsp){
		if(hsp == null){
			return false;
		}
		//System.out.print("\tadding HSP " + hsp.getQueryName() + " " + hsp.getHitName()  +" " +hsp.getHspRank());
		if(hsp.getQueryName().equalsIgnoreCase(this.queryID) && hsp.getHitName().equalsIgnoreCase(this.hitID)){
			
			
			if(hsps.size() + 1 == hsp.getHspRank() || hsp.getHspRank() == 0){
				hsps.add(hsp);
			}else{
				Vector<BlastHSP> v = new Vector<BlastHSP>();
				boolean added = false;
				for(Enumeration<BlastHSP> myenum =  hsps.elements(); myenum.hasMoreElements();){
					BlastHSP myHSP = myenum.nextElement();
					if(added|| myHSP.getBitScore() >= hsp.getBitScore()){
						v.add(myHSP);
					}else{
						v.add(hsp);
						added = true;
						v.add(myHSP);
						
					}
					
				}
			}
			
		}else{
			return false;
		}
		
		
		
		return true;
	}
	
	
	
	/**
	 * Adds a BlastHSP to this BlastHit if queryID, hitID are matching and the hit does not already contain this hsp
	 * @param b 
	 * 		the BlastHSP that will be added
	 * @return
	 * 		boolean: true if the hps was added, false if not.
	 */
	public boolean addHSP_old(BlastHSP b){
		//boolean addable = false;
		if( b == null){
			return false;
		}
		
		
		//do not add if query and hit do not match
		if       (  !b.getQueryName().equalsIgnoreCase(this.queryID) || 
			        !b.getHitName().equalsIgnoreCase(this.hitID)     
			       
			        ){
			return false;
		}
		
		
		//do not add if this hsp is already in the vector
		for(Enumeration<BlastHSP> myenum = hsps.elements(); myenum.hasMoreElements();){
			BlastHSP h = myenum.nextElement();
			
			if(b.getHitStart() == h.getHitStart() &&
			   b.getHitEnd() == h.getHitEnd() &&
			   b.getQueryStart() == h.getQueryStart ()&&
			   b.getQueryEnd() == h.getQueryEnd()){
				return false;
			}
		}
		
		
		
		//add the new hsp according to its
		Vector<BlastHSP> newHSPs = new Vector<BlastHSP>();
		boolean added = false;
		for(Enumeration<BlastHSP> myenum = this.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP hsp = myenum.nextElement();
			if(!added){
				if(hsp.getRawScore()<= b.getRawScore()){
					newHSPs.add(hsp);
				}else{
					newHSPs.add(b);
					newHSPs.add(hsp);
					added=true;
				}
				
			}else{
				newHSPs.add(hsp);
			}
			
		}
		
		
		if(!added){
			hsps.add(b);
		}
		
		this.hsps = newHSPs;
		
		return true;
	}
	
	
	public void setHitRank(int hitRank){
		this.hitRank = hitRank;
	}
	
	
	/**
	 * returns the rank of this hit. If this information is not available the return value is -1
	 * @return
	 */
	public int getHitRank(){
		if( hitRank >0 ){
			return hitRank;
		}else{
			return -1;
		}
	}
	
	public int getPercentHitCoverage(){
		int[] hitCoverage = this.getHitCoverage();
		
		int j = 0;
		for(int i = 0; i <  hitCoverage.length; i++){
			if(hitCoverage[i]>0){
				j++;
			}
		}
		
		return j*100/hitCoverage.length;
		
	}
	
	
	
	public int getPercentQueryCoverage(){
		int[] queryCoverage = this.getQueryCoverage();
		int j = 0;
		for(int i = 0; i <  queryCoverage.length; i++){
			if(queryCoverage[i]>0){j++;}
		}
		
		return j*100/queryCoverage.length;
		
	}
	
	public int getPercentLargestDistinctQueryCoverage(){
		int longest = 0;
		for(Enumeration<BlastHSP> myenum = hsps.elements(); myenum.hasMoreElements() ;){
			BlastHSP hsp= myenum.nextElement();
			int length = hsp.getHspLength();
			if(length > longest){
				longest = length;
			}
		}
		return longest*100/this.getQueryLength();
	}
	
	
	public int getQueryStart(){
		int queryStart = this.getQueryLength()+10;
		for(Enumeration<BlastHSP> myenum = this.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP b = myenum.nextElement();
			if(b.getQueryStart()< queryStart){
				queryStart = b.getQueryStart();
			}
		}
		return queryStart;
	}
	
	public int getHitStart(){
		int hitStart = this.getHitLength()+10;
		for(Enumeration<BlastHSP> myenum = this.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP b = myenum.nextElement();
			if(b.getHitStart()< hitStart){
				hitStart = b.getHitStart();
			}
		}
		return hitStart;
	}
	
	
	public int getQueryEnd(){
		int queryEnd = 0;
		for(Enumeration<BlastHSP> myenum = this.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP b = myenum.nextElement();
			if(b.getQueryEnd()> queryEnd){
				queryEnd = b.getQueryEnd();
			}
		}
		return queryEnd;
	}
	
	public int getHitEnd(){
		int hitEnd = 0;
		for(Enumeration<BlastHSP> myenum = this.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP b = myenum.nextElement();
			if(b.getHitEnd()> hitEnd){
				hitEnd = b.getHitEnd();
			}
		}
		return hitEnd;
	}
	/*
	public void setHitNum(int hitNum){
		this.hitNum = hitNum;
	}
	
	public int getHitNum(){
		if(hitNum <1){
			return -1;
		}
		return hitNum;
	}
	*/
	  
	 
	public String getQueryDescription()throws NullPointerException{
		
		return hsps.get(0).getQueryDescription();
	}
	
	public String getHitDescription()throws NullPointerException{
		
		return hsps.get(0).getHitDescription();
	}
	
	public String getQueryID() {
		return queryID;
	}

	

	public String getHitID() {
		return hitID;
	}

	public Vector<BlastHSP> getHSPs(){
		return hsps;
	}

	public Hashtable<Integer,BlastHSP> getHspsByRank() {
		Hashtable<Integer, BlastHSP> h = new Hashtable<Integer, BlastHSP>();
		int i = 1;
		for(Enumeration<BlastHSP> myenum = hsps.elements(); myenum.hasMoreElements();){
			h.put(i, myenum.nextElement());
			i++;
		}
		return h;
	}
	
	
	public double getLargestPercentIdentity(){
		double d = 0.0;
		for(Enumeration<BlastHSP> myenum = hsps.elements(); myenum.hasMoreElements();){
			BlastHSP h=myenum.nextElement();
			if(h.getIdentityPercentage()>d){
				d=h.getIdentityPercentage();
			}
			
		}
		return d;
	}
	
	/**
	 * returns a set of HSPs that have identity percentage larger or equal than given minimum identityPercentage. 
	 * The keys in the table are ints starting from 1 to size of hashtable, and preserver the original ranking.
	 * 
	 * @param minIdentity
	 * @return
	 */
	public Hashtable<Integer, BlastHSP> getHSPsByRank(double minIdentityPercentage){
		Hashtable<Integer, BlastHSP> h = new Hashtable<Integer,BlastHSP>();
		int i = 1;
		for(Enumeration<BlastHSP> myenum = hsps.elements(); myenum.hasMoreElements();){
			BlastHSP hsp = myenum.nextElement();
			if(hsp.getIdentityPercentage() >= minIdentityPercentage){
				h.put(i, myenum.nextElement());
				i++;
			}	
		}
		
		return h;
	}
	
	
	public BlastHSP getHSP(int rank){
		return hsps.get(rank-1);
	}

	
	/*
	public BlastHSP getRandomBestHSP(){
		
		if(hsps.size() ==1 || hsps.get(1).getNumericEvalue()< hsps.get(2).getNumericEvalue()){
			return hsps.get(1);
		}
		Vector<BlastHSP> v = new Vector<BlastHSP>();
		v.add(hsps.get(1));
		v.add(hsps.get(2));
		double evalue = hsps.get(1).getNumericEvalue();
		for( int i =3; i<= hsps.size(); i++){
			if(hsps.get(i).getNumericEvalue()<= evalue){
				v.add(hsps.get(i));
			}
		}
		
		int index = new Random().nextInt(v.size());
		
		return v.get(index);
		
	}
	*/
	
	
	public int getQueryLength() {
		return hsps.get(0).getQueryLength();
	}


	public int getHitLength() {
		return hsps.get(0).getHitLength();
	}

	public int[] getHitCoverage(){
		return getHitCoverage(0);
	}
	public int[] getQueryCoverage(){
		return getQueryCoverage(0);
	}
	
	
	
	public int[] getExtendedHitCoverage(){
		int[] hitCoverage = new int[this.getHitLength()];
		for(Enumeration<BlastHSP> myenum = this.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP b = myenum.nextElement();
			
			hitCoverage = mergeCoverages(hitCoverage, b.getHitCoverage());
			
		}
		return hitCoverage;
		
	}
	
	public int[] getHitCoverage(double minPercentIdentity){
		int[] hitCoverage = new int[this.getHitLength()];
		
		for(Enumeration<BlastHSP> myenum = this.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP b = myenum.nextElement();
			if(b.getIdentityPercentage()>= minPercentIdentity){
				hitCoverage = mergeCoverages(hitCoverage, b.getHitCoverage());
			}
		}
		return hitCoverage;
	}
	
	public int[] getQueryCoverage(double minPercentIdentity){
		int[] queryCoverage = new int[this.getQueryLength()];
		
		for(Enumeration<BlastHSP> myenum = this.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP b = myenum.nextElement();
			
			if(b.getIdentityPercentage()>= minPercentIdentity){
				
				queryCoverage = mergeCoverages(queryCoverage, b.getQueryCoverage());
			}
		}
		return queryCoverage;
	}
	
	
	
	/**
	 * Returns the evalue of the blast hit. This is the smallest evalue of all HSPs
	 * @return
	 */
	public double getEvalue(){
		double evalue = 10000;
		for(Enumeration<BlastHSP> myenum = this.hsps.elements();myenum.hasMoreElements();){
			double e = myenum.nextElement().getNumericEvalue();
			if( e < evalue){
				evalue = e;
			}
			
		}
		return evalue;
	}
	
	
	public double getBitScore(){
		double score = 1;
		for(Enumeration<BlastHSP> myenum = this.hsps.elements();myenum.hasMoreElements();){
			double e = myenum.nextElement().getBitScore();
			if( e > score){
				score = e;
			}
			
		}
		return score;
	}
	
	public double getRawScore(){
		double score = 1;
		for(Enumeration<BlastHSP> myenum = this.hsps.elements();myenum.hasMoreElements();){
			double e = myenum.nextElement().getRawScore();
			if( e > score){
				score = e;
			}
			
		}
		return score;
		
	}
	
	public static int[] mergeCoverages(int[] cov1, int[] cov2)throws ArrayIndexOutOfBoundsException{
		int[] newcov = new int[cov1.length];
		if(cov1.length != cov2.length){
			throw new ArrayIndexOutOfBoundsException (" you cannot merge two array with no equal length here");
		}
		
		for(int i = 0; i < newcov.length; i++){
			newcov[i] = cov1[i]+cov2[i];
		}
		
		return newcov;
	}
	
	
	public static double getPercentCoverage(int[] cov){
		int j = 0;
		for(int i = 0; i <  cov.length; i++){
			if(cov[i]>0){j++;}
		}
		double pcov = (double) j * 100.0 /cov.length;
		pcov = (Math.round(pcov*100)) / 100; 
		
		
		return pcov;
		
	
	
		
	}
	

	/**
	 * @deprecated
	 * compare two Blast hits with the same query and return weather first or second is better. If the hits cover different parts of the query it is also returned
	 * @param hit1
	 * @param hit2
	 * @return
	 *      an int describing the result:
	 *      		1:	First BlastHit is better
	 *      		2:	Second BlastHit is better
	 *      		0:	Both hits cover different parts of the query
	 *      	   -1:	Query of the two BlastHits is not equal
	 */
	public static int compareTwoHitResults(BlastHit hit1, BlastHit hit2){
		if( !hit1.getQueryID().equalsIgnoreCase(hit2.getQueryID())){
			return -1;
		}
		
		int[] q1 = hit1.getQueryCoverage();
		int[] q2 = hit2.getQueryCoverage();
		int numOverlapping = 0;
		for( int i = 0; i< q1.length;i++){
			if(q1[i] >0 && q2[i] > 0){
				numOverlapping++;
			}
		}
		if(numOverlapping <10){
			return 0;
		}
		
		if( hit1.getPercentQueryCoverage() >= hit2.getPercentQueryCoverage()){
			return 1;
		}else{
			return 2;
		}
		
		
		
		
	}
	
	
	
	/**
	 * 
	 * return the set of BlastHits representing the query in the best way. 
	 * It takes into account the hightest level of query coverage (% query coverage), 
	 * the highest identity (Num Identical of all HSPs * 100 / Sum HSP lengths) 
	 * and the lowest redundance (Percent of bases of the query sequence represented by only one HSP)
	 * 
	 * @param v
	 * 		input Vector of all BlastHits
	 * @return
	 * 		output Vector of the best combinatino of BlastHits
	 */
	public static Vector<BlastHit> getBestHitCombination(Vector<BlastHit> v){
		
		Hashtable<Vector<Integer>, Integer> result = new Hashtable<Vector<Integer>, Integer>();
		
		
		for(int i = 0; i< v.size(); i++){
			for ( int j = v.size()-1; j>=i;j--){
				int[] q = new int[v.get(0).getQueryCoverage().length];
				int numIdentical = 0;
				int sumLengths = 0;
				Vector<Integer> vector = new Vector<Integer>();
				String s = "";
				for( int k = i; k<=j; k++){
					s = s + " " + v.get(k).getHitID();
					vector.add(k);
					q = BlastHit.mergeCoverages(q, v.get(k).getQueryCoverage());
					for( Enumeration<BlastHSP> myenum = v.get(k).getHSPs().elements(); myenum.hasMoreElements();){
						BlastHSP hsp = myenum.nextElement();
						numIdentical = numIdentical + hsp.getNumIdentical();
						sumLengths =  sumLengths + hsp.getHspLength();
					}
					
				}
				int sum =0;
				int num =0;
				for( int m = 0; m< q.length;m++){
					sum = sum  + q[m];
					if( q[m] ==1){
						num++;
					}
				}
				
				double value = BlastHit.getPercentCoverage(q) + (numIdentical *100.0 /sumLengths) + (num*100.0/v.get(0).getQueryLength());
				result.put(vector, (int)value);
				//System.out.println(((int)BlastHit.getPercentCoverage(q))+  "%\t" + (numIdentical *100 /sumLengths)+"%\t"+(num*100/v.get(0).getQueryLength())+"%\t"+s);
			}
			
			
		}
		
		Vector<Integer> vector = new Vector <Integer>();
		double best = 0;
		for(Enumeration<Vector<Integer>> myenum = result.keys(); myenum.hasMoreElements();){
			Vector<Integer> key = myenum.nextElement();
			double d = result.get(key);
			if( d > best){
				best = d;
				vector = key;
			}
		}
		
		Vector<BlastHit > resultVector = new Vector<BlastHit>();
		for(Enumeration<Integer> myenum = vector.elements(); myenum.hasMoreElements();){
			int i = myenum.nextElement();
			resultVector.add(v.get(i));
		}
		
		return resultVector;
	}
	
	
	
	
	
	
	
	public static double getCombinedQueryCoverage(Vector <BlastHit> v){
		int[] cov = new int[v.get(0).getQueryCoverage().length];
		for(Enumeration<BlastHit> myenum = v.elements(); myenum.hasMoreElements();){
			BlastHit b = myenum.nextElement();
			cov = mergeCoverages(cov, b.getQueryCoverage());
		}
		return BlastHit.getPercentCoverage(cov);
	}
	
	public static double getCombinedIdentity(Vector <BlastHit> v){
		int numIdentical = 0;
		int sumLengths = 0;
		for(Enumeration<BlastHit> myenum = v.elements(); myenum.hasMoreElements();){
			BlastHit b = myenum.nextElement();
			for(Enumeration<BlastHSP> myenum2 = b.getHSPs().elements(); myenum2.hasMoreElements();){
				BlastHSP hsp = myenum2.nextElement();
				numIdentical = numIdentical + hsp.getNumIdentical();
				sumLengths = sumLengths + hsp.getHspLength();
			}
		}
		return ( numIdentical * 100.0 / sumLengths  );
	}
	
	public static double getCombinedIdentity(BlastHit b){
		int numIdentical = 0;
		int sumLengths = 0;
		
			
			for(Enumeration<BlastHSP> myenum2 = b.getHSPs().elements(); myenum2.hasMoreElements();){
				BlastHSP hsp = myenum2.nextElement();
				numIdentical = numIdentical + hsp.getNumIdentical();
				sumLengths = sumLengths + hsp.getHspLength();
			}
		
		return ( numIdentical * 100.0 / sumLengths  );
	}
	
	
	public static double getNonRedundantlyCoveredBases(Vector <BlastHit> v){
		int[] cov = new int[v.get(0).getQueryCoverage().length];
		for(Enumeration<BlastHit> myenum = v.elements(); myenum.hasMoreElements();){
			BlastHit b = myenum.nextElement();
			cov = mergeCoverages(cov, b.getQueryCoverage());
		}
		int num=0;
		
		for( int i = 0; i< cov.length; i++){
			if( cov[i]==1){
				num++;
			}
		}
		return ( num * 100.0 / cov.length);
	}
	
	
	
	
	public String getGalaxyTSV(){
		String s = "";
		for( int i = 1; i<= this.getHSPs().size(); i++){
			s = s + this.getHspsByRank().get(i).getGalaxyTSV() + "\n";
		}
		return s;
	}
	public String getIPKParserTSV(){
		StringBuilder s = new StringBuilder();
		for(Enumeration<BlastHSP> myenum = this.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP hsp = myenum.nextElement();
			s.append(hsp.getIPKParserTSV() + "\n");
			
		}
		
		return s.toString();
	}
	public static String getHeaderForIPKParserTSV(){
		return "query_name\tquery_description\tquery_length\tnum_hits\tdatabase_name\talgorithm\thit_name\thit_description\thit_length\tnum_hsps\tevalue\t"+
		           "score\tnum_identical\tnum_positives\tnum_gaps\tbit_score\tidentity_percentage\thsp_rank\tquery_string\talignment_string\thit_string\t"+
				   "hsp_length\thsp_hit_length\thsp_query_length\thit_start\tquery_start\thit_end\tquery_end\thit_frame\tquery_frame\thit_strand\tquery_strand\n";	
	}  
	
	
	
	public void discardLowIdentityHSPs(double minPercentIdentity){
		Vector<BlastHSP> v = new Vector<BlastHSP>();
		int rankCounta=0;
		for(Enumeration<BlastHSP> myenum = this.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP hsp = myenum.nextElement();
			if(hsp.getIdentityPercentage()>=minPercentIdentity){
				rankCounta++;
				hsp.setHspRank(rankCounta);
				v.add(hsp);
			}
		}
		this.hsps = v;
	}
	
	
	
	
	/**
	 * 
	 * @param arg0
	 * @return
	 */
	  public int compareTo(BlastHit arg0) {
			BlastHit o2 =  arg0;
			
			if( this.getQueryID().compareTo(o2.getQueryID()) != 0){
				return this.getQueryID().compareTo(o2.getQueryID());
			}
			if(this.getBitScore() > o2.getBitScore()){
				return 1;
			}else if(this.getBitScore()<o2.getBitScore()){
				return -1;
			}else{
				return 0;
			}
	    	
	     
	   }
	   
	
	
	public void removeHSP(double identityThreshold){
		Vector<BlastHSP> v = new Vector<BlastHSP>();
		for(Enumeration<BlastHSP> myenum = this.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP hsp = myenum.nextElement();
			if(hsp.getIdentityPercentage()>=identityThreshold){
				v.add(hsp);
			}
		}
		this.hsps = v;
	}
	
	
}
