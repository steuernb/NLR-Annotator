package core;

import java.util.Collections;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import annotator.AnnotatorSignatureDefinition;




/**
 * 
 * 
 * 
 *
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
 *  @author Burkhard Steuernagel
 * 
 *
 *
 *
 *
 *
 * this data structure holds a set of motifs. Typically, this holds all motifs of one sequence. 
 * In NLR-Annotator, it is also used to define a set of motifs that resembles an NLR locus.
 *
 *
 *
 *
 *
 *
 */
public class MotifList implements Comparable<MotifList> {

	
	
	String motifListName; // this is the place holder for the name of an NLR locus. 
	
	
	
	Vector<Motif> motifs; // the set of motifs in this list
	boolean sorted; //keep track of whether this set of motifs is sorted by position in the sequence
	
	
	
	/**
	 * 
	 * Creates this object with a given list of motifs.
	 * 
	 * @param motifListName
	 * @param motifs
	 */
	public MotifList(String motifListName, Vector<Motif> motifs) {
		this.motifListName = motifListName;
		this.motifs = motifs;
		sorted = false;
		this.sort();
		
	}
	
	
	/**
	 * 
	 * Creates this object without motifs, in order to add them subsequently
	 * 
	 * @param motifListName
	 */
	public MotifList(String motifListName) {
		this.motifListName = motifListName;
		this.motifs = new Vector<Motif>();
		sorted = false;
	}
	
	public Vector<Motif> getMotifs(){
		return this.motifs;
	}
	
	public void addMotif(Motif motif) {
		this.motifs.add(motif);
		sorted = false;
		sort();
		
	}
	
	public boolean addMotifs(MotifList motifList) {
		return addMotifs(motifList.getMotifs());
	}
	
	public boolean addMotifs(Vector<Motif> motifList) {
		
		
		this.motifs.addAll(motifList);
		
		sorted = false;
		sort();
		return true;
	}
	
	public Motif firstMotif() {
		
		if(!sorted) {
			sort();
		}
		return this.motifs.firstElement();
	}
	public Motif lastMotif() {
		
		if(!sorted) {
			sort();
		}
		return this.motifs.lastElement();
	}
	
	
	
	
	
	public String getMotifString() {
		
		this.sort();
		
		
		
		String s = "";
		for(Iterator<Motif> iterator = this.motifs.iterator(); iterator.hasNext();) {
			Motif motif = iterator.next();
			String motif_id = motif.getID();
			s = s + "," + motif_id;
		}
		
		if( s.length()>0) {
			return s.substring(1);
		}else {
			return null;
		}
		
	}
	
	public String getMotifListString() {
		if(this.motifs.size() == 0) {
			return " no motifs found ";
		}
		this.sort();
		
		String s = "";
		
		for(Iterator<Motif> iterator = this.motifs.iterator(); iterator.hasNext();) {
			Motif motif = iterator.next();
			s = s +","+motif.getID();
		}	
		
		return s.substring(1);
	}
	
	
	
	/**
	 * 
	 * @deprecated
	 * 
	 * @return
	 */
	public String getReportString(){
		
		this.sort();
		
		
		String reportString = "";
		String sequence_id = null;
		
		String motifList = "";
		
		for(Iterator<Motif> iterator = this.motifs.iterator(); iterator.hasNext();) {
			Motif motif = iterator.next();
			
			
			if(sequence_id ==null) {
				sequence_id = motif.getProteinSequenceID();
			}else {
				if(!sequence_id.equalsIgnoreCase(motif.getProteinSequenceID())) {
					return "error - inconsistent sequence IDs for motifs";
				}
			}
			
			motifList = motifList + "," + motif.getID().substring(6);
			
		}
		if(motifList.length() > 0) {
			motifList = motifList.substring(1);
		}
		reportString = sequence_id  + "\t" + motifList;
		
		return reportString;
	}
	
	
	
	
	
	
	
	
	public boolean has_NLR_signature(SignatureDefinition signatureDefinition) {
		
		
		if(signatureDefinition.hasSignature(this)) {
			return true;
		}else {
			return false;
		}
		
	}
	
	public void sort() {
		if(!sorted) {
			Collections.sort(motifs);
			sorted = true;
		}
		
	}
	
	public boolean isSorted() {
		return this.sorted;
	}
	
	
	/**
	 * 
	 * set the DNA parameters for each motif in this list
	 * 
	 * 
	 * @param dna_sequence_id
	 * @param offset
	 * @param fragmentLength
	 * @param frame
	 * @param forwardStrand
	 */
	public void setDNA(String dna_sequence_id, long offset, int fragmentLength, int frame, boolean forwardStrand) {
		
		for(Iterator<Motif> iterator = this.motifs.iterator(); iterator.hasNext();) {
			
			Motif motif = iterator.next();
			motif.setDNA(dna_sequence_id, offset, fragmentLength, frame, forwardStrand);
		
		}
		
		
	}
	
	
	public String getMotifListName() {
		return this.motifListName;
	}
	
	
	public void setMotifListName(String motifListName) {
		this.motifListName = motifListName;
	}
	
	
	
	public String getSequenceName() {
		return this.motifs.firstElement().getDnaSequenceID();
	}
	
	
	
	public boolean isForwardStrand() {
		
		return this.motifs.firstElement().forwardStrand;
	}
	
	
	
	
	/**
	 * 
	 * Check if an NLR is complete with regards to the canonical denifition of an NLR (i.e. NB-ARC domain and LRRs)
	 * 
	 * This concept might need to be revised. Complete/Partial is misleading since there are functional genes that are classified as "partial" with this definition.
	 * 
	 * 
	 * @param annotatorSignatureDefinition
	 * @return
	 */
	public  boolean isCompleteNLR(AnnotatorSignatureDefinition annotatorSignatureDefinition){
		
		//check if lowest rank NBARC motif is present. That would be the p-loop
		
		boolean hasPloop = false;
		boolean hasLRR = false;
		
		
		for(Iterator<Motif> iterator = motifs.iterator(); iterator.hasNext();) {
			Motif motif = iterator.next();
			
			if( annotatorSignatureDefinition.isLRR(motif)) {
				hasLRR = true;
			}
			if( annotatorSignatureDefinition.isPloop(motif)) {
				hasPloop = true;
			}
			
		}
		
		
		
		
		
		return hasPloop && hasLRR ;
	}
	
	
	public boolean hasStopCodon(){
		
		for(Iterator<Motif > iterator = this.getMotifs().iterator(); iterator.hasNext();){
			Motif motif = iterator.next();
			//System.out.println(motif.getID() +"\t" +motif.getDnaStart() + "\t" + motif.getProteinSequence() + "\t" + motif.score);
			if(motif.getProteinSequence().contains("*")){
				return true;
			}
		}
		
		return false;
	}
	
	public int getNumberOfMotifs() {
		return this.motifs.size();
	}
	
	
	
	
	public void removeRedundantMotifs(){
		HashSet<String> h = new HashSet<String>();
		Vector<Motif> v = new Vector<Motif>();
		for(Enumeration<Motif> myenum = this.motifs.elements(); myenum.hasMoreElements();){
			Motif motif = myenum.nextElement();
			String s = motif.getID() + "\t" + motif.getDnaStart() + "\t" +motif.getDnaEnd();
			if(!h.contains(s)){
				v.add(motif);
			}
			h.add(s);
		}
		
		this.motifs = v;
	}
	
	
	
	
	/**
	 * 
	 * Check if this MotifList overlaps with another Vector of Motifs.
	 * This method is called when NLR seeds are being merged: NLR_Annotator.findNLRs_substep2_mergeSeeds()
	 * 
	 * 
	 * @param motifList
	 * @return
	 */
	public boolean overlapsWith(Vector<Motif> motifList) {
		for(Iterator<Motif> iterator1 = this.motifs.iterator(); iterator1.hasNext();) {
			Motif motif1 = iterator1.next();
			for(Iterator<Motif> iterator2 = motifList.iterator(); iterator2.hasNext();) {
				Motif motif2 = iterator2.next();
				if( motif1.equals(motif2)) {
					return true;
				}
			}
			
			
		}
		
		return false;
	}
	
	
	/**
	 * Check if this MotifList overlaps with another Vector of Motifs.
	 * If at least one motif is shared, it checks consistency if merged
	
	 * 
	 * 
	 * @param motifList
	 * @return
	 */
	public boolean overlapsWith(MotifList motifList, AnnotatorSignatureDefinition annotatorSignatureDefinition) {
		HashSet<Motif> test = new HashSet<Motif>();
		for(Iterator<Motif> iterator1 = this.motifs.iterator(); iterator1.hasNext();) {
			Motif motif1 = iterator1.next();
			test.add(motif1);
			for(Iterator<Motif> iterator2 = motifList.getMotifs().iterator(); iterator2.hasNext();) {
				Motif motif2 = iterator2.next();
				test.add(motif2);
				if( motif1.equals(motif2)) {
					Vector<Motif > v= new Vector<Motif>();
					v.addAll(test);
					Collections.sort(v);
					int lastRank = annotatorSignatureDefinition.getMotifRank(v.firstElement());
					for( int i = 1; i< v.size(); i++) {
						int rank = annotatorSignatureDefinition.getMotifRank(v.get(i));
						if( rank > lastRank || annotatorSignatureDefinition.isLRR(v.get(i))) {
							lastRank = rank;
						}else {
							return false;
						}
					}
					
					
					
					return true;
				}
			}
			
			
		}
		
		return false;
	}
	
	/**
	 * 
	 * check if two motif combinations can be merged. This will involve the correct distance and the right order according to annotator signature definition
	 * 
	 * 
	 * 
	 * @param motifList
	 * 			The other motif list
	 * @param distanceBetweenMotifCombinations
	 * 			The shortest distance between DNA starts of any pair of motifs from the two lists.
	 * @param annotatorSignatureDefinition
	 * 			The AnnotatorSignatureDefinition that contains all iformation about motif combinations.
	 * @return
	 */
	public boolean canBeMergedWith(MotifList motifList, int distanceBetweenMotifCombinations, AnnotatorSignatureDefinition annotatorSignatureDefinition) {
		
		sort();
		motifList.sort();
		
		
		//check strand;
		if( this.isForwardStrand() && !motifList.isForwardStrand()) {
			return false;
		}
		if(!this.isForwardStrand() && motifList.isForwardStrand()){
			return false;
		}
		
		//check distance
		long distance1 = Math.abs(this.firstMotif().getDnaStart() - motifList.firstMotif().getDnaStart());
		long distance2 = Math.abs(this.firstMotif().getDnaStart() - motifList.lastMotif().getDnaStart());
		long distance3 = Math.abs(this.lastMotif().getDnaStart() - motifList.firstMotif().getDnaStart());
		long distance4 = Math.abs(this.lastMotif().getDnaStart() - motifList.lastMotif().getDnaStart());
		
		
		
		if( distance1 > distanceBetweenMotifCombinations
			&& distance2 > distanceBetweenMotifCombinations
			&& distance3 > distanceBetweenMotifCombinations
			&& distance4 > distanceBetweenMotifCombinations
			
		) {
			
			return false;
			
			
			
		}
		
		//check consistency	
			
		HashSet<Motif> h = new HashSet<Motif>();
		h.addAll(this.getMotifs());
		h.addAll(motifList.getMotifs());
		Vector<Motif> v = new Vector<Motif>();
		v.addAll(h);
		Collections.sort(v);
		int lastRank = annotatorSignatureDefinition.getMotifRank(v.firstElement());
		for( int i = 1; i< v.size(); i++) {
			int rank = annotatorSignatureDefinition.getMotifRank(v.get(i));
			if( rank > lastRank || annotatorSignatureDefinition.isLRR(v.get(i))) {
				lastRank = rank;
			}else {
				
				return false;
			}
		}
		
		
		/*
		if(this.isForwardStrand() && motifList.isForwardStrand()) { //on the forward strand
			
			
			if(this.lastMotif().getDnaStart() <= motifList.firstMotif().getDnaStart()) { // other motif list is behind
				if( motifList.firstMotif().getDnaStart() - this.lastMotif().getDnaStart() <= distanceBetweenMotifCombinations) { // and in good distance
					
					if(annotatorSignatureDefinition.getMotifRank( this.lastMotif() ) <  annotatorSignatureDefinition.getMotifRank(motifList.firstMotif())) { //even the motif ranks match
						return true;
					}
					
					
					
				}
				
				
			}else if( this.firstMotif().getDnaStart() >= motifList.lastMotif().getDnaStart() ) { //other motif list is before
				if(this.firstMotif().getDnaStart() - motifList.lastMotif().getDnaStart() <= distanceBetweenMotifCombinations ) { //and in good distance
				
					if(annotatorSignatureDefinition.getMotifRank( motifList.lastMotif())  < annotatorSignatureDefinition.getMotifRank(this.firstMotif()) ) { //even the motif ranks match
						return true;
					}
				}
			}
			
			
			
		}
		
		
		
		
		
		if(!this.isForwardStrand() && !motifList.isForwardStrand()) { //on the reverse strand
			
			
			if(this.firstMotif().getDnaStart() <= motifList.lastMotif().getDnaStart()) { // other motif list is behind (keep in mind, motifs are sorted by position on chromosome.)
				if( motifList.lastMotif().getDnaStart() - this.firstMotif().getDnaStart() <= distanceBetweenMotifCombinations) { // and in good distance
					
					if(annotatorSignatureDefinition.getMotifRank( this.firstMotif() ) >  annotatorSignatureDefinition.getMotifRank(motifList.lastMotif())) { //even the motif ranks match
						return true;
					}
					
					
					
				}
				
				
			}else if( this.lastMotif().getDnaStart() >= motifList.firstMotif().getDnaStart() ) { //other motif list is before
				if(this.lastMotif().getDnaStart() - motifList.firstMotif().getDnaStart() <= distanceBetweenMotifCombinations ) { //and in good distance
				
					if(annotatorSignatureDefinition.getMotifRank( motifList.firstMotif())  > annotatorSignatureDefinition.getMotifRank(this.lastMotif()) ) { //even the motif ranks match
						return true;
					}
				}
			}
			
			
			
			
		}
		*/
		
		return true;
	}
	
	
	
	/**
	 * 
	 * get the distance between two seeds. Also check if they fit together.
	 * If the distance is beyond threshold for distance or if they don't fit by motif rank return -1. Otherwise, return the distance.
	 * 
	 * 
	 * @param list
	 * 			The other MotifList that is compared to this one.
	 * @param distanceBetweenMotifCombinations
	 * 			Threshold for distance of the two motif lists
	 * @param annotatorSignatureDefinition
	 * 			The set of rules for NLR-Annotator. Here, this is used to determine the rank of a motif.
	 * @return
	 * 		Either -1 if the motifLists don't fit or the distance.
	 */
	public long getDistance(MotifList list,int distanceBetweenMotifCombinations, AnnotatorSignatureDefinition annotatorSignatureDefinition) {
		long distance = -1;
		if(this.isForwardStrand() && list.isForwardStrand()) { //forward strand
			if(this.lastMotif().getDnaStart() < list.firstMotif().getDnaStart()) { //other list behind
				
				distance = list.firstMotif().getDnaStart() - this.lastMotif().getDnaStart();
				if( distance <= distanceBetweenMotifCombinations) {
					
					if( annotatorSignatureDefinition.getMotifRank(this.lastMotif()) < annotatorSignatureDefinition.getMotifRank(list.firstMotif())) {
						return distance;
					}else {
						return -1;
					}
					
					
				}
				
				
			}else if( this.firstMotif().getDnaStart() > list.lastMotif().getDnaStart()) { //other list before
				
				distance = this.firstMotif().getDnaStart() - list.lastMotif().getDnaStart();
				if( distance <= distanceBetweenMotifCombinations) {
					
					if( annotatorSignatureDefinition.getMotifRank(this.firstMotif()) > annotatorSignatureDefinition.getMotifRank(list.lastMotif())) {
						return distance;
					}
					
					
				}
				
			}else {//overlapping? should not happen
				System.out.println("Debug: checking distance between seeds. This case should not happen. Position: " + list.firstMotif().getDnaStart());
				return 0;
			}
			
			
		}else if(!this.isForwardStrand() && !list.isForwardStrand()) { //reverse strand
			
			if( list.firstMotif().getDnaStart() < this.lastMotif().getDnaStart()) { 
				distance = this.lastMotif().getDnaStart() - list.firstMotif().getDnaStart();
				if(distance <= distanceBetweenMotifCombinations) {
					if(annotatorSignatureDefinition.getMotifRank(list.firstMotif()) > annotatorSignatureDefinition.getMotifRank(this.lastMotif())) {
						return distance;
					}
					
				}
				
			}else if ( this.firstMotif().getDnaStart() < list.lastMotif().getDnaStart()) {
				distance = list.lastMotif().getDnaStart() - this.firstMotif().getDnaStart();
				if( distance <= distanceBetweenMotifCombinations) {
					if(annotatorSignatureDefinition.getMotifRank(this.firstMotif()) > annotatorSignatureDefinition.getMotifRank(list.lastMotif())) {
						return distance;
					}
				}
			}else {
				//System.out.println("Debug: checking distance between seeds (rv strand). This case should not happen. Position: " + list.firstMotif().getDnaStart());
				return 0;
			}
		}
		
		
		return -1;
	}
	
	
	public long getDistance(Vector<Motif> motifs,int distanceBetweenMotifCombinations, AnnotatorSignatureDefinition annotatorSignatureDefinition) {
		MotifList list = new MotifList("temp", motifs);
		return getDistance(list, distanceBetweenMotifCombinations, annotatorSignatureDefinition);
		
		
	}
	
	
	/**
	 * 
	 * Check if any of the motifs of this list are part of a given HashSet.
	 * 
	 * @param motifs
	 * @return
	 */
	public boolean containsUsedMotif(HashSet<Motif> usedmotifs) {
		
		for(Iterator<Motif> iterator = motifs.iterator(); iterator.hasNext();) {
			Motif motif = iterator.next();
			if(usedmotifs.contains(motif)) {
			
				return true;
			}
		}
		return false;
		
	}
	
	
	
	
	@Override
	public int compareTo(MotifList list) {
		
		if( this.isForwardStrand() && !list.isForwardStrand()) {
			return -1;
		}
		if( !this.isForwardStrand() && list.isForwardStrand()) {
			return 1;
		}
		
		if(this.isForwardStrand()) { //both forward
			if( this.firstMotif().getDnaStart() < list.firstMotif().getDnaStart()) {
				return -1;
			}else
			if( this.firstMotif().getDnaStart()  > list.firstMotif().getDnaStart()) {
				return 1;
			}else { //first motif same start position. get longer one first.
				if(this.getMotifs().size() > list.getMotifs().size() ) {
					return -1;
				}else
				if( this.getMotifs().size() < list.getMotifs().size()) {
					return 1;
				}else {
					return 0;
				}
				
			}
			
		}else { //both reverse
			
			
			if( this.firstMotif().getDnaStart() < list.firstMotif().getDnaStart()) {
				return 1;
			}else
			if( this.firstMotif().getDnaStart()  > list.firstMotif().getDnaStart()) {
				return -1;
			}else { //first motif same start position. get longer one first.
				if(this.getMotifs().size() > list.getMotifs().size() ) {
					return -1;
				}else
				if( this.getMotifs().size() < list.getMotifs().size()) {
					return 1;
				}else {
					return 0;
				}
				
			}
		}
		
	}
	
}
