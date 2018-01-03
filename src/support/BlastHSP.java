package support;

import java.util.Enumeration;
import java.util.Vector;

public class BlastHSP implements Comparable<BlastHSP> {

	
	
	String queryName;
	String queryDescription;
	String hitName;
	String hitDescription;
	
	
	double rawScore;
	double bitScore;
	
	int queryLength;
	int hitLength;
	int hspLength;
	
	String evalue;
	
	int queryStart;
	int hitStart; 
	int queryEnd;
	int hitEnd;
	
	int hspRank;
	
	
	String algorithm;
	String database_name; 
	
	
	
	
	int hitStrand;
	int queryStrand;
	int hitFrame;
	int queryFrame;
	int numIdentical;
	int numMismatches;
	int numGaps;
	int numGapOpenings;
	int numPositive;
	
	String queryString;
	String hitString;
	String alignmentString;

	
	String allHitIDs;
	int numHits;
	int numHSPs;
	
	
	
	public BlastHSP(String queryName, String hitName, String queryDescription, String hitDescription, 
			        int hspLength, int queryLength, int hitLength, 
			        int numIdentical, int numMisMatches, int numGaps, int numGapOpenings, int numPositives,
			        int queryStart, int queryEnd, int hitStart, int hitEnd,
			        double rawScore, double bitScore, String evalue, String queryString, String hitString,
			        int queryFrame, int hitFrame, int queryStrand, int hitStrand){
		
		this.queryName = queryName;
		this.hitName = hitName;
		this.queryDescription = queryDescription;
		this.hitDescription = hitDescription;
		this.hspLength = hspLength;
		this.queryLength = queryLength;
		this.hitLength = hitLength;
		this.numIdentical = numIdentical;
		this.numMismatches = numMisMatches;
		this.numGaps = numGaps;
		this.queryStart = queryStart;
		this.queryEnd = queryEnd;
		this.hitStart = hitStart;
		this.hitEnd = hitEnd;
		this.rawScore = rawScore;
		this.bitScore = bitScore;
		this.evalue = evalue;
		this.queryString = queryString;
		this.hitString = hitString;
		this.numGapOpenings = numGapOpenings;
		this.allHitIDs = "N/A";
		this.numPositive = numPositives;
		this.queryFrame = queryFrame;
		this.hitFrame = hitFrame;
		this.queryStrand = queryStrand;
		this.hitStrand = hitStrand;
	}
	

	public String getIPKParserTSV(){
		return  this.getQueryName() + "\t" +
			    this.getQueryDescription()+ "\t" +
			    this.getQueryLength()+ "\t" +	
			    this.getNumHits()+ "\t" +
			    this.getDatabaseName()+ "\t" +
			    this.getAlgorithm()+ "\t" +
			    this.getHitName()  + "\t" +
			    this.getHitDescription()+ "\t" +
				this.getHitLength()+ "\t" +
			    this.getNumHSPs()+ "\t" +
				this.getEvalue()+ "\t" +
			    this.getRawScore()+ "\t" +
				this.getNumIdentical()+ "\t" +
			    this.getNumPositives()+ "\t" +
				this.getNumberOfGaps()+ "\t" +
			    this.getBitScore()+ "\t" +
				this.getIdentityPercentage()+ "\t" +
			    this.getHspRank()+ "\t" +
				this.getQueryString()+ "\t" +
			    this.getAlignmentString()+ "\t" +
				this.getHitString()+ "\t" +
				this.getHspLength()+ "\t" +
				this.getHSPHitLength()+ "\t" +
				this.getHSPQueryLength()+ "\t" +
				this.getHitStart()+ "\t" +
				this.getQueryStart()+ "\t" +
				this.getHitEnd()+ "\t" +
				this.getQueryEnd()+ "\t" +
				this.getHitFrame()+ "\t" +
				this.getQueryFrame()+ "\t" +
				this.getHitStrand()+ "\t" +
				this.getQueryStrand();
	}
	
	public int compareTo(BlastHSP otherHSP){
		if(this.getBitScore()> otherHSP.getBitScore()){
			return 1;
		}else if(this.getBitScore()< otherHSP.getBitScore()){
			return -1;
		}else{return 0;
		
		}
	}
	
	
	public String getHeaderForIPKParserTSV(){
		return "query_name\tquery_description\tquery_length\tnum_hits\tdatabase_name\talgorithm\thit_name\thit_description\thit_length\tnum_hsps\tevalue\t"+
	           "score\tnum_identical\tnum_positives\tnum_gaps\tbit_score\tidentity_percentage\thsp_rank\tquery_string\talignment_string\thit_string\t"+
			   "hsp_length\thsp_hit_length\thsp_query_length\thit_start\tquery_start\thit_end\tquery_end\thit_frame\tquery_frame\thit_strand\tquery_strand";	
	}
	
	public int getNumHits(){
		return this.numHits;
	}

	public int getNumHSPs(){
		return this.numHSPs;
	}
	public void setNumHits(int numHits){
		this.numHits = numHits;
	}
	
	public void setNumHits(BlastIteration it){
		this.numHits = it.getHits().size();
	}
	
	
	public void setNumHSPs(int num){
		this.numHSPs = num;
	}
	
	public void setNumHSPs(BlastHit hit){
		this.numHSPs = hit.getHSPs().size();
	}
	
		
	public int getHSPHitLength(){
		return Math.max(this.getHitEnd() ,this.getHitStart() ) - Math.min(this.getHitEnd() ,this.getHitStart() );
	}
	public int getHSPQueryLength(){
		return this.getQueryEnd() - this.getQueryStart();
	}
	
	public int getHitFrame() {
		return hitFrame;
	}


	public int getQueryFrame() {
		return queryFrame;
	}


	public double getPositivePercentage(){
		return (double)this.numPositive *100.0/(double)this.hspLength;
	}
	
	public int getNumPositives(){
		return this.numPositive;
	}

	public String getAllHitIDs(){
		if(allHitIDs == null){
			return "N/A";
		}
		return allHitIDs;
	}
	
	public int getNumberOfMismatches(){
		return this.numMismatches;
	}
	public int getNumberOfGapOpenings(){
		return this.numGapOpenings;
	}
	
	public String getEvalue() {
		if(evalue.startsWith("e")){
			evalue = "1"+evalue;
		}
		return evalue;
	}
	
	/*
	public String getAlignment() {
		return alignment;
	}
	*/
	
	public double getNumericEvalue(){
		return Double.parseDouble(this.getEvalue());
	}
	
	
	public double getBitScore(){
		return this.bitScore;
	}
	
	public double getRawScore(){
		return this.rawScore;
	}
	

	public int getQueryEnd() {
		return queryEnd;
	}

	

	public String getQueryName() {
		return queryName;
	}

	
	public String getQueryDescription() {
		return queryDescription;
	}

	

	public String getHitName() {
		return hitName;
	}

	

	public String getHitDescription() {
		return hitDescription;
	}


	public int getQueryLength() {
		return queryLength;
	}

	public int getHitLength() {
		return hitLength;
	}

	public int getHspLength() {
		return hspLength;
	}

	public int getQueryStart() {
		return queryStart;
	}

	public int getHitStart() {
		return hitStart;
	}


	public int getHitEnd() {
		return hitEnd;
	}


	public int getHspRank() {
		return hspRank;
	}

	public void setHspRank(int hspRank) {
		this.hspRank = hspRank;
	}

	
	public int getHspLength_IgnoreNs(){
		int l = 0;
		char[] c1 = this.getQueryString().toCharArray();
		char[] c2 = this.getHitString().toCharArray();
		for( int i = 0; i< c1.length;i++){
			if(c1[i]!='N' && c2[i]!='N'){
				l++;
			}
		}
		return l;
	}
	
	/*
	public String getIdentityString() {
		return identityString;
	}

	public void setIdentityString(String identityString) {
		this.identityString = identityString;
	}
    */
	

	

	public double getIdentityPercentage() {
		return (double) this.numIdentical * 100.0 / (double)this.hspLength;
	}


	public int getNumberOfGaps() {
		return numGaps;
	}

	

	public int getHitStrand() {
		if(hitStrand == 0){
			if(hitFrame<0){
				hitStrand = -1;
			}
			if(hitFrame > 0){
				hitStrand = 1;
			}
		}
		
		return hitStrand;
	}


	public int getQueryStrand() {
		
		if(queryStrand == 0){
			if(queryFrame<0){
				queryStrand = -1;
			}
			if(queryFrame > 0){
				queryStrand = 1;
			}
		}
		
		
		return queryStrand;
	}



	public String getQueryString() {
		return queryString;
	}

	

	public String getHitString() {
		return hitString;
	}
	
	
	public String getAlignmentString(){
		if(alignmentString!= null){
			return this.alignmentString;
		}
		char[] q = this.getQueryString().toCharArray();
		char[] h = this.getHitString().toCharArray();
		
		if( q.length != h.length){
			return null;
		}
		
		String r = "";
		for( int i = 0; i< q.length; i++){
			if(q[i] == h[i]){
				r = r + "|";
			}else{
				r = r + " ";
			}
		}
		return r;
	}
	
	public void setAlignmentString(String alignmentString){
		this.alignmentString = alignmentString;
	}
	
	
	public int getNumIdentical(){
		return this.numIdentical;
	}
	
	public int[] getQueryCoverage(){
		int[] cov = new int[this.getQueryLength()];
		
		int start = this.getQueryStart();
		int end = this.getQueryEnd();
		if(start>end){
			end = start;
			start = this.getQueryEnd();
		}
		
		for( int i = this.getQueryStart()-1; i< this.getQueryEnd();i++){
			try{
				cov[i]++;
			}catch(ArrayIndexOutOfBoundsException e){throw new ArrayIndexOutOfBoundsException (e.getMessage() + " in  Query "+this.getQueryName() + "; Hit "+this.getHitName()+"; HSP "+this.getHspRank() );}	
		}
		return cov;
	}
	public int[] getHitCoverage(){
		int[] cov = new int[this.getHitLength()];
		
		int start = this.getHitStart();
		int end = this.getHitEnd();
		
		if( start > end){
			end = start;
			start = this.getHitEnd();
		}
		
		for( int i = start-1; i< end;i++){
			cov[i]++;
		}
		return cov;
	}
	
	public int[] getExtendedHitCoverage(){
		int[] cov = new int[this.getHitLength()];
		for( int i = this.getHitStart()-1; i< this.getHitEnd();i++){
			cov[i]++;
		}
		return cov;
	}
	
	
	
	public double getPercentQueryCoverage(){
		int[] cov = this.getQueryCoverage();
		int j = 0;
		for(int i = 0; i <  cov.length; i++){
			if(cov[i]>0){j++;}
		}
		double pcov = (double) j * 100.0 /cov.length;
		pcov = (Math.round(pcov*100)) / 100; 
		
		
		return pcov;
		
	}
	
	
	
	public Vector<String[]> getMisMatches(){
		Vector<String[]> v = new Vector<String[]>();
		
		char[] querySeq = this.getQueryString().toCharArray();
		char[] homolSeq = this.getAlignmentString().toCharArray();
		char[] hitSeq   = this.getHitString().toCharArray();
		
		int numQgaps = 0;
		int numHgaps = 0;
		
		for(int i = 0; i< homolSeq.length; i++){
			if(homolSeq[i] != '|'){
				int qPos = this.getQueryStart()+i-numQgaps;
				int hPos = this.getHitStart()+i-numQgaps;
				String q = (querySeq[i]+"").toUpperCase();
				String h = (hitSeq[i]+"").toUpperCase();
				
				if( this.getHitStrand()==-1){
					hPos = this.getHitEnd()-i+numHgaps;
				}
				String[] s = {this.getQueryName(),this.getHitName(),qPos+"", hPos+"", q,h};
				v.add(s);
			}
			if(querySeq[i]=='-'){
				numQgaps++;
			}
			if(hitSeq[i]=='-'){
				numHgaps++;
			}
			
		}
		
		
		
		
		return v;
	}
	
	
	public void setQueryDescription(String queryDescription){
		this.queryDescription = queryDescription;
	}
	public void setHitDescription(String hitDescription){
		this.hitDescription = hitDescription;
	}
	
	public void setQueryName(String queryName){
		this.queryName = queryName;
	}
	
	public void setAllHitIDs(String hitIDs){
		this.allHitIDs = hitIDs;
	}
	public void setAllHitIDs(BlastIteration it){
		StringBuilder builder = new StringBuilder();
		for(Enumeration<BlastHit> myenum = it.getHits().elements();myenum.hasMoreElements();){
			builder.append(";"+myenum.nextElement().getHitID());
			
		}
		this.allHitIDs = builder.toString().substring(1);
	}
	
	
	
	
	
	
	
	
	
	
	public static BlastHSP setHSPFromGalaxy(String galaxyLine){
		String[] split = galaxyLine.split("\t");
		BlastHSP hsp = new BlastHSP(
				split[0], split[1], "", "", 
		        Integer.parseInt(split[3]), Integer.parseInt(split[22]), Integer.parseInt(split[23]), 
		        Integer.parseInt(split[14]), Integer.parseInt(split[4]), Integer.parseInt(split[16]), Integer.parseInt(split[5]), Integer.parseInt(split[15]),
		        		Integer.parseInt(split[6]), Integer.parseInt(split[7]), Integer.parseInt(split[8]), Integer.parseInt(split[9]),
		        Double.parseDouble(split[13]), Double.parseDouble(split[11]), split[10], split[20], split[21],
		        Integer.parseInt(split[18]), Integer.parseInt(split[19]), 0, 0
				);
		
		
		return hsp;
	}
	
	public static BlastHSP setHSPFromIPKParser(String input)throws NumberFormatException{
		String[] split = input.split("\t");
			try {
			BlastHSP hsp = new BlastHSP(
					split[0], split[6], split[1], split[7], 
					Integer.parseInt(split[21]), Integer.parseInt(split[2]), Integer.parseInt(split[8]), 
			        Integer.parseInt(split[12]), -1, Integer.parseInt(split[14]), -1, Integer.parseInt(split[13]),
			        Integer.parseInt(split[25]),Integer.parseInt(split[27]), Integer.parseInt(split[24]), Integer.parseInt(split[26]),
			        Double.parseDouble(split[15]), Double.parseDouble(split[11]), split[10], split[18], split[20],
			        Integer.parseInt(split[29]), Integer.parseInt(split[28]), Integer.parseInt(split[31]), Integer.parseInt(split[30])
					);
			
			hsp.setHspRank(Integer.parseInt(split[17]));
			hsp.setNumHits(Integer.parseInt(split[3]));
			hsp.setNumHSPs(Integer.parseInt(split[9]));
			hsp.setAlgorithm(split[5]);
			hsp.setAlignmentString(split[19]);
			hsp.setDatabaseName(split[4]);
			return hsp;
		} catch (Exception e) {
			throw new NumberFormatException("Parsing of this line failed:\n" + input);
		}
		
		
		
	}
	
	
	public String getGalaxyTSV(){
		return this.getQueryName() + "\t" +
			   this.getHitName()  + "\t" +
			   this.getIdentityPercentage() + "\t" +
			   this.getHspLength() + "\t" +
			   this.getNumberOfMismatches() + "\t" +
			   this.getNumberOfGapOpenings() + "\t" +
			   this.getQueryStart() + "\t" +
			   this.getQueryEnd() + "\t" +
			   this.getHitStart() + "\t" +
			   this.getHitEnd() + "\t" +
			   this.getEvalue() + "\t" +
			   this.getBitScore() + "\t" +
			   this.getAllHitIDs() + "\t" +
			   this.getRawScore() + "\t" +
			   this.getNumIdentical() + "\t" +
			   this.getNumPositives() + "\t" +
			   this.getNumberOfGaps() + "\t" +
			   this.getPositivePercentage() + "\t" +
			   this.getQueryFrame() + "\t" +
			   this.getHitFrame() + "\t" +
			   this.getQueryString() + "\t" +
			   this.getHitString() + "\t" +
			   this.getQueryLength() + "\t" +
			   this.getHitLength() + "\t" +
			   this.getQueryStrand() + "\t" +
			   this.getHitStrand() ;
	}
	
	
	public void setAlgorithm(String algorithm){
		this.algorithm = algorithm;
	}
	public String getAlgorithm(){
		if(this.algorithm != null){
			return this.algorithm;
		}
		return "N/A";
	}
	
	public void setDatabaseName(String db_name){
		this.database_name = db_name;
	}
	public String getDatabaseName(){
		if(this.database_name != null){
			return this.database_name;
		}
		return "N/A";
	}
	
	public String getGffQuery(){
		String strand = "+";
		if(this.getHitStrand() <0){
			strand = "-";
		}
		
		String frame = ".";
		if( this.getQueryFrame() != 0){
			frame = (Math.abs(this.getQueryFrame() -1)) + "";
		}
		return this.getQueryName() + "\t" + this.algorithm + "\t" +  "blast_hsp" + "\t" + this.getQueryStart() + "\t" + this.getQueryEnd() + "\t" +this.getIdentityPercentage() + "\t" + strand + "\t" + frame + "\t" + "name " + this.getHitName();
	}
	public String getGffHit(){
		String strand = "+";
		
		int start = this.getHitStart();
		int end = this.getHitEnd();
		if(this.getHitStrand() <0){
			strand = "-";
			start = this.getHitEnd();
			end = this.getHitStart();
		}
		
		String frame = ".";
		if( this.getHitFrame() != 0){
			frame = (Math.abs(this.getQueryFrame() -1)) + "";
		}
		return this.getHitName() + "\t" + this.algorithm + "\t" +  "blast_hsp" + "\t" + start + "\t" + end + "\t" +this.getIdentityPercentage() + "\t" + strand + "\t" + frame + "\t" + "name " + this.getQueryName();
	}
	
	
	
	
	
	
}
