package support;

import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

public class BlastIteration {

	
	Vector<BlastHit> hits;
	String query;
	int queryLength;
	String algorithm;
	
	public BlastIteration(String query){
		this.query = query;
		hits = new Vector<BlastHit>();
		
	}
	
	public boolean addHit(BlastHit hit){
		if (hit.getQueryID().equalsIgnoreCase(query)){
			hits.add(hit);
			return true;
		}else{
			return false;
		}
		
	}

	
	public String getQueryID(){
		return query;
	}
	
	public int getQueryLength(){
		if(queryLength != 0){
			return queryLength;
		}else{
			return hits.get(0).getQueryLength();
		}
		
	}
	
	public void setQueryLength(int queryLength){
		this.queryLength = queryLength;
	}
	
	public Vector<BlastHit> getHits(){
		return hits;
	}
	
	
	public void removeHSP(double identityThreshold){
		Vector<BlastHit> newHits= new Vector<BlastHit>();
		for(Enumeration<BlastHit> myenum1 = this.getHits().elements(); myenum1.hasMoreElements();){
			BlastHit b = myenum1.nextElement();
			b.removeHSP(identityThreshold);
			if(b.getHSPs().size()>0){
				newHits.add(b);
			}
		}
		this.hits=newHits;
	}
	
	
	/**
	 * returns a vector of BlastHits from this iteration where the best HSP is longer or equal then a minimum required alignment length
	 * 
	 * @param minAlignmentLength
	 * @return
	 */
	public Vector<BlastHit> getHits(int minAlignmentLength){
		Vector<BlastHit> v = new Vector<BlastHit>();
		for(Enumeration<BlastHit> myenum = this.hits.elements(); myenum.hasMoreElements();){
			BlastHit hit = myenum.nextElement();
			if(hit.getHSP(1).getHspLength() >= minAlignmentLength){
				v.add(hit);
			}
		}
		return v;
	}
	
	public Vector<BlastHit> getHits(double minIdentityPercentage){
		Vector<BlastHit> v = new Vector<BlastHit>();
		for(Enumeration<BlastHit> myenum = this.hits.elements(); myenum.hasMoreElements();){
			BlastHit hit = myenum.nextElement();
			if(hit.getHSP(1).getIdentityPercentage() >= minIdentityPercentage){
				v.add(hit);
			}
		}
		return v;
	}
	
	public Vector<BlastHit> getHits(int minAlignmentLength, HashSet<String> hitFilter){
		Vector<BlastHit> v = new Vector<BlastHit>();
		for(Enumeration<BlastHit> myenum = this.hits.elements(); myenum.hasMoreElements();){
			BlastHit hit = myenum.nextElement();
			if(hit.getHSP(1).getHspLength() >= minAlignmentLength && hitFilter.contains(hit.getHitID())){
				v.add(hit);
			}
		}
		return v;
	}
	
	public int getNumberOfHits(){
		return hits.size();
	}
	
	
	public String getIPKParserTSV(){
		StringBuilder s = new StringBuilder();
		for( int i = 0; i< this.getHits().size(); i++){
			s.append( this.getHits().get(i).getIPKParserTSV() );
		}
		return s.toString();
	}
	/*
	public String getHeaderForIPKParserTSV(){
		return "query_name\tquery_description\tquery_length\tnum_hits\tdatabase_name\talgorithm\thit_name\thit_description\thit_length\tnum_hsps\tevalue\t"+
		           "score\tnum_identical\tnum_positives\tnum_gaps\tbit_score\tidentity_percentage\thsp_rank\tquery_string\talignment_string\thit_string\t"+
				   "hsp_length\thsp_hit_length\thsp_query_length\thit_start\tquery_start\thit_end\tquery_end\thit_frame\tquery_frame\thit_strand\tquery_strand";	
	}
	*/
	public static String getHeaderForIPKParserTSV(){
		return "query_name\tquery_description\tquery_length\tnum_hits\tdatabase_name\talgorithm\thit_name\thit_description\thit_length\tnum_hsps\tevalue\t"+
		           "score\tnum_identical\tnum_positives\tnum_gaps\tbit_score\tidentity_percentage\thsp_rank\tquery_string\talignment_string\thit_string\t"+
				   "hsp_length\thsp_hit_length\thsp_query_length\thit_start\tquery_start\thit_end\tquery_end\thit_frame\tquery_frame\thit_strand\tquery_strand";	
	}
	
	public String getQueryDescription(){
		return hits.get(0).getQueryDescription();
	}
	
	public int[] getQueryCoverage(){
		int[] cov = new int[ this.getQueryLength()];
			
		for(Enumeration<BlastHit> myenum1 = hits.elements(); myenum1.hasMoreElements();){
			BlastHit hit = myenum1.nextElement();
			cov = BlastHit.mergeCoverages(cov,hit.getQueryCoverage());
		}
		return cov;
	}
	
	public int[] getQueryCoverage(double identityPercentage){
		int[] cov = new int[ this.getQueryLength()];
			
		for(Enumeration<BlastHit> myenum1 = hits.elements(); myenum1.hasMoreElements();){
			BlastHit hit = myenum1.nextElement();
			cov = BlastHit.mergeCoverages(cov,hit.getQueryCoverage(identityPercentage));
		}
		return cov;
	}
	
	
	public double getPercentQueryCoverage(){
		
		int counta = 0;
		int[] cov = this.getQueryCoverage();
		for( int i = 0; i< cov.length; i++){
			if(cov[i]>0){
				counta++;
			}
		}
		return (Math.round((double) counta*10000.0/cov.length))/100;
		
	}
	
	
	
	
	
	
	
	
	public String getAlgorithm(){
		if(algorithm == null){
			return "N/A";
		}
		return this.algorithm;
	}
	public void setAlgorithm(String algorithm){
		this.algorithm = algorithm;
	}
	
	
	
	public void discardLowIdentityHSPs(double minPercentIdentity){
		
		Vector<BlastHit> v = new Vector<BlastHit>();
		int hitRankCounta = 0;
		for(Enumeration<BlastHit> myenum1 = this.getHits().elements(); myenum1.hasMoreElements();){
			BlastHit hit = myenum1.nextElement();
			hit.discardLowIdentityHSPs(minPercentIdentity);
			if(hit.getHSPs().size()>0){
				hitRankCounta++;
				hit.setHitRank(hitRankCounta);
				v.add(hit);
			}
		}
		this.hits = v;
		
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
	public Vector<BlastHit> getBestHitCombination(double minPercentIdentity){
		
		Vector<BlastHit> v = this.getHits(minPercentIdentity);
		
		Hashtable<Vector<Integer>, Integer> result = new Hashtable<Vector<Integer>, Integer>();
		
		
		for(int i = 0; i< v.size(); i++){
			//System.out.println("index i: " + i + " " + v.size());
			for ( int j = v.size()-1; j>=i;j--){
				//System.out.println("index j: " + j + " " + v.size());
				int[] q = new int[v.get(0).getQueryCoverage().length];
				int numIdentical = 0;
				int sumLengths = 0;
				Vector<Integer> vector = new Vector<Integer>();
				String s = "";
				for( int k = i; k<=j; k++){
					//System.out.println("index k: " + k + " " + v.size());
					s = s + " " + v.get(k).getHitID();
					vector.add(k);
					//System.out.println("num hsps " + v.get(k).getHSPs().size());
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
	public Vector<BlastHit> getBestHitCombination(int minAlignmentLength, HashSet<String> hitFilter){
		
		Vector<BlastHit> v = this.getHits(minAlignmentLength, hitFilter);
		
		Hashtable<Vector<Integer>, Integer> result = new Hashtable<Vector<Integer>, Integer>();
		
		
		for(int i = 0; i< v.size(); i++){
			//System.out.println("index i: " + i + " " + v.size());
			for ( int j = v.size()-1; j>=i;j--){
				//System.out.println("index j: " + j + " " + v.size());
				int[] q = new int[v.get(0).getQueryCoverage().length];
				int numIdentical = 0;
				int sumLengths = 0;
				Vector<Integer> vector = new Vector<Integer>();
				String s = "";
				for( int k = i; k<=j; k++){
					//System.out.println("index k: " + k + " " + v.size());
					s = s + " " + v.get(k).getHitID();
					vector.add(k);
					//System.out.println("num hsps " + v.get(k).getHSPs().size());
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
	
	
	
	
	public String getGff(){
		String s = "";
		for(Enumeration<BlastHit> myenum1= this.getHits().elements(); myenum1.hasMoreElements();){
			BlastHit hit = myenum1.nextElement();
			for(Enumeration<BlastHSP> myenum2 = hit.getHSPs().elements(); myenum2.hasMoreElements();){
				BlastHSP hsp = myenum2.nextElement();
				s = s+ hsp.getGffQuery() + "\n";
			}
		}
		return s;
	}
	
	public String getGff(int numHits){
		//System.out.println(this.getQueryID());
		String s = "";
		int count = 0;
		for(Enumeration<BlastHit> myenum1= this.getHits().elements(); myenum1.hasMoreElements();){
			BlastHit hit = myenum1.nextElement();
			
			for(Enumeration<BlastHSP> myenum2 = hit.getHSPs().elements(); myenum2.hasMoreElements();){
				BlastHSP hsp = myenum2.nextElement();
				s = s+ hsp.getGffQuery() + "\n";
				
			}
			count++;
			//System.out.println(count);
			if(count>=numHits){
				break;
			}
		}
		return s;
	}
	
	
	
}
