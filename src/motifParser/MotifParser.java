package motifParser;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import core.Motif;
import core.MotifList;
import core.SignatureDefinition;
import support.BioSequence;

/**
 *   NLR-Annotator - find signatures of NLR-type resistance genes in genomic sequence
 *   
 *   Copyright (C) 2021  John Innes Centre
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 * 
 * 
 *  @author Jitender Cheema
 *  @author Burkhard Steuernagel
 * 
 */
public class MotifParser {

	private int k;
	 
    private Map<Integer, Double> non_overlapping_pvalues;

	
	MotifDefinition motifDefinition;
	BioSequence aminoAcidSequence;
	
	public MotifParser(MotifDefinition motifDefinition, BioSequence aminoAcidSequence)throws IOException {
		
	    this.motifDefinition = motifDefinition; 
	    this.aminoAcidSequence = aminoAcidSequence;
	}
	
	

	
	public MotifList findMotifs( ) throws IOException{
		
		
		
		
		MotifList motifList = new MotifList(aminoAcidSequence.getIdentifier());
			
		double thresh = 0.0001D;
	        

	    if( aminoAcidSequence.getLength() < this.motifDefinition.getMaximumMotifLength()) {
	    	
	       return motifList;
	       
	    }
	        

	       
	    int sequenceLength = aminoAcidSequence.getLength();
	    String protein =  aminoAcidSequence.getSequence();
	                
	            	
	   
	    for (String m : motifDefinition.getMotifNames()) {

	    	int width = motifDefinition.getLength(m);
	    	double lowest_pval = 2.0;
	        for (int i = 0; i < sequenceLength - width + 1; i++) {
	        	double p = pval_at(protein, m, i);
	        	if (p < lowest_pval) { 
	        		lowest_pval = p; 
	        	} // best p
	        }
	    }

	    
	    Map<Integer, Integer> hits = get_motifs_stats(protein, sequenceLength, thresh);
	    	
	    	
    	for(int q = 0; q < sequenceLength; q++) {
    		
    		
    		if ( hits.get(q) > 0) {
    			
    			int motif_number =  hits.get(q);
    			
	            double non_overlap_p = non_overlapping_pvalues.get(q);
	            //System.out.println( "MotifParser: motif_number " +motif_number + "; q="+q+"; non_overlap_p="+non_overlap_p );           
	            if( non_overlap_p < 1E-5) {
                        	  
	            	//String motif_namer = String.format("motif_%d", motif_number);
	            	String motif_namer = this.motifDefinition.getMotifNames().get(motif_number -1);
	            	//if(aminoAcidSequence.getIdentifier().equalsIgnoreCase("AT1G51480.1")) {
	            	//	System.out.println(aminoAcidSequence.getIdentifier() + "\t" + motif_namer + "\t" + aminoAcidSequence.getLength() + "\t" + q + "\t" + motifDefinition.getLength(motif_namer));
	            	//}
	            	
	            	Motif motif = new Motif(	motif_namer, 
	            								aminoAcidSequence.getIdentifier(), 
                    		   					q+1, 
	                        		   			aminoAcidSequence.getSequence().substring(q,(q+motifDefinition.getLength(motif_namer))),
                    		   					non_overlap_p
	                        		   					);
	                        	
	            	motifList.addMotif(motif);
	            }
    		}
    	}
	
    	//System.out.println("MotifList " + motifList.getMotifListName() + ": " + motifList.getNumberOfMotifs());
    	
	   return motifList;
	   
		
		
	}
	
	
	
	
	
	public MotifDefinition getMotifDefinition() {
		return this.motifDefinition;
	}
	

	
	private Map<Integer, Integer> get_motifs_stats(String protein, int  sequenceLength, double thresh){

        Map<Integer, Integer> hits = new HashMap<Integer, Integer>();
        for (int i = 0; i < sequenceLength; i++) {
            hits.put(i, 0);
        }
        Map<Integer, Double> pvalues = new HashMap<Integer, Double>();
        for (int i = 0; i < sequenceLength; i++) {
            pvalues.put(i, 1.0d);
        }
        int imotif;
        int ws = 0;
        int maxws= 0;
       
        boolean ok_to_mark = false;
        double prod;
        int m;
        
        
        for(imotif = 0; imotif < this.motifDefinition.getMotifNames().size(); imotif++){
            String mot = this.motifDefinition.getMotifNames().get(imotif);
          
       
        	ws = motifDefinition.getLength(mot);
            maxws = Math.max(ws, maxws);
       

            for (int j = 0; j < sequenceLength - ws + 1; j++) {
                double pvalue = pval_at(protein, mot, j);
                if (pvalue < thresh) {
                    ok_to_mark = true;
                    prod = 1.0;

                    int first = Math.max(0, j - maxws + 1);
                    int last  = (int) Math.min(j + ws, sequenceLength);

                    // get product of p-values overlap right
                    for(int  k=j; k < last; k++){
                        if (ok_to_mark) {
                            m = Math.abs(hits.get(k)) - 1;
                            if (hits.get(k) != 0) {
                                prod *= pvalues.get(k);
                                if (pvalue >= prod) {
                                    ok_to_mark = false;
                                }
                            }
                        }
                    }
                    // get product left
                    for(int k = first; k < j; k++){
                        if (ok_to_mark) {
                            m = Math.abs(hits.get(k)) - 1;
                            if ( m >= 0 &&  motifDefinition.getLength(this.motifDefinition.getMotifNames().get(m)) > j - k ){
                               prod *= pvalues.get(k);
                               if (pvalue >= prod){
                                   ok_to_mark = false;
                               }
                            }

                        }
                    }
                    //  Marking
                    if (ok_to_mark){
                    	
                        hits.put(j, imotif +1 );
                        pvalues.put(j, pvalue );

                        // remove from right
                        for(k=j+1; k < last; k++){
                            hits.put(k, 0);
                            pvalues.put(k, 1d);
                        }

                        // remove from left
                        for(k=first; k < j; k++){
                            m = Math.abs(hits.get(k)) - 1;
                            if ( m >= 0 &&  motifDefinition.getLength(this.motifDefinition.getMotifNames().get(m)) > j - k ) {
                                hits.put(k, 0);
                                pvalues.put(k, 1d);
                            }
                        }
                    } // ok_to_mark



                } // best p
            } // imotif for
        }

        non_overlapping_pvalues = pvalues;
        return hits;
    } 

	
	
	

	    private double pval_at(String sequence, String motif, int pos){

	        int width = motifDefinition.getLength(motif);
	        String kmer  = sequence.substring(pos, pos+width);
	        double kmer_score = 0.0d;
	        for (int i = 0; i < kmer.length(); i++) {
	            char base= kmer.charAt(i);
	            
	           
	            //double score = pos2score.get(String.format("%s@%d@%s", motif, i, base));
	            double score = motifDefinition.getScore(motif, i, base);
	            kmer_score += score;
	        }
	        int kmer_score_int = (int) kmer_score;
	        //  get position p-value of score
	        double pvalue  =  motifDefinition.getCDF(motif, kmer_score_int);
	        //double pvalue  = cdf_store.get(String.format("%s@%d", motif, kmer_score_int));
	        return pvalue;
	    }
	
	
}
