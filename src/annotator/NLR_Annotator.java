package annotator;

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
import java.util.Date;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;
import java.util.concurrent.ExecutionException;
import java.util.zip.GZIPInputStream;

import core.Motif;
import core.MotifList;
import core.SignatureDefinition;
import motifParser.MotifDefinition;
import motifParser.MotifParser;
import motifParser.MotifParserParallelExecutor;
import support.BioSequence;
import support.commandLineInterface.CLI;
import support.commandLineInterface.CLIParseException;




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
 *  @author Burkhard Steuernagel
 * 
 */
public class NLR_Annotator {

	
	public static final double version = 2.0;
	
	
	AnnotatorSignatureDefinition annotatorSignatureDefinition;
	SignatureDefinition signatureDefinition;
	MotifDefinition motifDefinition;
	
	
	Vector<MotifList> nlrs;
	HashMap<String,Vector<Motif>> motifs;
	
	
	
	public NLR_Annotator(	File inputFasta, 
							File motFile, 
							File storeFile, 
							int fragmentLength, 
							int overlap, 
							int numSequencesPerThread, 
							int numThreads,
							int distanceWithinMotifCombination,
							int distanceForElongating, 
							int distanceBetweenMotifCombinations)throws IOException, InterruptedException, ExecutionException {
		
		long startTime = System.nanoTime();
		
		
		nlrs = new Vector<MotifList>();
		
		this.annotatorSignatureDefinition = new AnnotatorSignatureDefinition();
		this.signatureDefinition = new SignatureDefinition();
		this.motifDefinition = new MotifDefinition(motFile, storeFile);
		
		MotifParserParallelExecutor executor = new MotifParserParallelExecutor(motifDefinition, signatureDefinition, inputFasta, numThreads, numSequencesPerThread, fragmentLength, overlap);
		this.motifs = executor.executeDNA();
		executor.close();
		
		long timestamp2 = System.nanoTime();
		System.out.println("Ran MotifParser with " +numThreads+" threads in " +((timestamp2  - startTime)/1000000000) +" seconds" );
		
		this.removeRedundantMotifs();
		
		
		this.findNLRs(distanceWithinMotifCombination, distanceForElongating, distanceBetweenMotifCombinations);
		
		
		
	}
	
	public NLR_Annotator(File importFile, int fragmentLength, 
						int overlap,
		
						int distanceWithinMotifCombination,
						int distanceForElongating, 
						int distanceBetweenMotifCombinations)throws IOException{
		
		long startTime = System.nanoTime();
		
		
		this.annotatorSignatureDefinition = new AnnotatorSignatureDefinition();
		
		this.importMotifParserResult(importFile);
		
		
		long timestamp1 = System.nanoTime();
		System.out.println("pre-calculated motifs loaded. " + (((double)timestamp1  - (double)startTime)/1000000000));
		
		this.findNLRs(distanceWithinMotifCombination, distanceForElongating, distanceBetweenMotifCombinations);
		
		
	}
	
	/**
	@deprecated
	*/
	public NLR_Annotator(File inputFasta, 
						File motFile, 
						File storeFile, 
						int fragmentLength, 
						int overlap,
		
						int distanceWithinMotifCombination,
						int distanceForElongating, 
						int distanceBetweenMotifCombinations
					)throws IOException {
		
		
		long startTime = System.nanoTime();
		
		
		this.annotatorSignatureDefinition = new AnnotatorSignatureDefinition();
		this.signatureDefinition = new SignatureDefinition();
		this.motifDefinition = new MotifDefinition(motFile, storeFile);
		
		motifs = new HashMap<String, Vector<Motif>>();
		
		long timestamp1 = System.nanoTime();
		System.out.println("initialized. " + (((double)timestamp1  - (double)startTime)/1000000000));
		
		
		this.runMotifParser(inputFasta, fragmentLength, overlap);
		
		long timestamp2 = System.nanoTime();
		System.out.println("Ran MotifParser in " +((timestamp2  - timestamp1)/1000000000) +" seconds" );
		
		this.removeRedundantMotifs();
		
		
		this.findNLRs(distanceWithinMotifCombination, distanceForElongating, distanceBetweenMotifCombinations);
		
		
		
		
		
	}
	
	
	
	
	
	public static double getVersion(){
		return version;
	}
	
	public int getNumberOfNLRs(){
		return nlrs.size();
	}
	
	
	public Vector<MotifList> getNLRs(){
		return nlrs;
	}
	
	public HashMap<String,Vector<Motif>> getMotifs(){
		return this.motifs;
	}
	
	
	
	private void runMotifParser(File inputFastaFile,int fragmentLength, int overlap) throws IOException{
		SequenceChopper chopper = new SequenceChopper(inputFastaFile,  fragmentLength,  overlap);
		
		for(BioSequence seq = chopper.nextSequence(); seq != null; seq = chopper.nextSequence() ) {
		
			BioSequence[] proteinSeqs = seq.translate2Protein();
			
			
			Vector<MotifList> mylist = new Vector<MotifList>(); //store all motifs
			boolean foundInForwardStrand = false;
			boolean foundInReverseStrand = false;
			
			for (int i = 0; i< proteinSeqs.length; i++) {
				BioSequence proteinSequence = proteinSeqs[i];
				
				MotifParser mp = new MotifParser(motifDefinition, proteinSequence);
				MotifList list = mp.findMotifs();
				mylist.add(list);
				
				if( list.has_NLR_signature(signatureDefinition)) {
					
					
					
					String[] split = list.getMotifListName().split("_");
					int l = split.length-1;
					String strand = split[l].substring(5, 6);
				
					if( strand.equalsIgnoreCase("+")) {
						foundInForwardStrand = true;
					}else {
						foundInReverseStrand = true;
					}	
					
				}
				
			}	
			
			if(foundInReverseStrand || foundInForwardStrand) {
				for(Iterator<MotifList> iterator2 = mylist.iterator(); iterator2.hasNext();) {
					
					MotifList list = iterator2.next();
					String[] split = list.getMotifListName().split("_");
					int l = split.length-1;
					String strand = split[l].substring(5, 6);
				
					boolean forwardStrand = true;
					if(strand.equalsIgnoreCase("-")) {
						forwardStrand = false;
					}
					int frame = Integer.parseInt(split[l].substring(6));
					
					int offset = Integer.parseInt(split[l-1]);
					
					String id = split[0];
					for( int j = 1; j< l-1; j++) { //in case id has an underscore in name
						id = id + "_" + split[j];
					}
					
					if( (forwardStrand && foundInForwardStrand)   || (!forwardStrand && foundInReverseStrand) ){
						if( ! motifs.containsKey(id)) {
							motifs.put(id,new Vector<Motif>());
						}
						for(Iterator<Motif> iterator = list.getMotifs().iterator(); iterator.hasNext();) {
							Motif motif = iterator.next();
							motif.setDNA(id, offset, fragmentLength, frame, forwardStrand);
							motifs.get(id).add(motif);
						}
					}
					
					
					
					
					
				}
			}	
				
					
					
				
				
			
			
			
		}
		//this.printMotifs();
	}
	
	
	/**
	 * 
	 * debug method. Prints all motifs in a motifs
	 * 
	 
	 */
	public void printMotifs() {
		
		for(Iterator<String> iterator = this.motifs.keySet().iterator(); iterator.hasNext();){
			String chr = iterator.next();
			Vector<Motif> v = this.motifs.get(chr);
			
			for(Iterator<Motif> iterator2 = v.iterator(); iterator2.hasNext();) {
				Motif motif = iterator2.next();
				System.out.println(chr + "\t" + motif.getMotifString());
			}
			
			
		}
		
		
	}
	
	/**
	 * 
	 * debug method. Prints all motifs in a Vector<Motif>
	 * 
	 * @param motifs
	 */
	public static void printMotifs(Vector<Motif> motifs) {
		
		for(Iterator<Motif> iterator = motifs.iterator(); iterator.hasNext();){
			Motif motif = iterator.next();
			System.out.println( motif.getMotifString());
			
			
			
		}
		
		
	}
	
	
	
	
	
	
	
	
	
	
	private void findNLRs(int distanceWithinMotifCombination,int distanceForElongating, int distanceBetweenMotifCombinations) {
		
		for(Iterator<String> iterator1 = this.motifs.keySet().iterator(); iterator1.hasNext();){
			
			String dna_seq_id = iterator1.next();
			Vector<Motif> v = this.motifs.get(dna_seq_id);
			
			
			
		
			
			
			Collections.sort(v);
			Vector<MotifList> motifSeeds = findNLRs_substep1_findSeeds(v, distanceWithinMotifCombination);
			
			try {
				writeBED(motifSeeds, new File("/Users/steuernb/Documents/projects/Annotator2/test_data/seeds1.bed"));
			} catch (Exception e) {}
			
			Vector<MotifList> mergedSeeds = findNLRs_subsep2a_mergeSeeds(motifSeeds, distanceBetweenMotifCombinations);
			
			try {
				writeBED(mergedSeeds, new File("/Users/steuernb/Documents/projects/Annotator2/test_data/seeds2.bed"));
			} catch (Exception e) {}
			
			
			Vector<MotifList> preNlrs = findNLRs_substep3_elongate(mergedSeeds, v, distanceForElongating);
			try {
				writeBED(preNlrs, new File("/Users/steuernb/Documents/projects/Annotator2/test_data/seeds3.bed"));
			} catch (Exception e) {}
			
			
			
			/*
			
			Vector<Motif> forwardMotifs = new Vector<Motif>();
			Vector<Motif> reverseMotifs = new Vector<Motif>();
			
			for(Enumeration<Motif> myenum2 = v.elements(); myenum2.hasMoreElements();){
				Motif motif = myenum2.nextElement();
				if(!motif.isForwardStrand()){
					reverseMotifs.add(motif);
					
				}else{
					forwardMotifs.add(motif);
				}
			}
			
			
			
		
			
			
			
			
			
			
			
			Collections.sort(forwardMotifs);
			Collections.sort(reverseMotifs);
			//Collections.reverse(reverseMotifs);
			
			
			Vector<Vector<Motif>> forwardMotifSeeds = findNLRs_substep1_findSeeds(forwardMotifs, distanceWithinMotifCombination);
			
			Vector<MotifList> forwardMergedSeeds = findNLRs_substep2_mergeSeeds(forwardMotifSeeds, distanceBetweenMotifCombinations);
			
			try {
				writeBED(forwardMergedSeeds, new File("/Users/steuernb/Documents/projects/Annotator2/test_data/debug1.bed"));
			}catch (IOException e) {
				e.printStackTrace();
			}
			
			Vector<MotifList> forwardPreNlrs = findNLRs_substep3_elongate(forwardMergedSeeds, forwardMotifs, distanceForElongating);
			try {
				writeBED(forwardPreNlrs, new File("/Users/steuernb/Documents/projects/Annotator2/test_data/debug2.bed"));
			}catch (IOException e) {
				e.printStackTrace();
			}
			
			Vector<Vector<Motif>> reverseMotifSeeds = findNLRs_substep1_findSeeds(reverseMotifs, distanceWithinMotifCombination);
			Vector<MotifList> reverseMergedSeeds = findNLRs_substep2_mergeSeeds(reverseMotifSeeds, distanceBetweenMotifCombinations);
			Vector<MotifList> reversePreNlrs = findNLRs_substep3_elongate(reverseMergedSeeds, reverseMotifs, distanceForElongating);
			
			
			if(this.nlrs == null) {
				nlrs = new Vector<MotifList>();
			}
			
			Vector<MotifList> allPreNlrs = new Vector<MotifList>();
			allPreNlrs.addAll(forwardPreNlrs);
			allPreNlrs.addAll(reversePreNlrs);
			*/
			
			
			
			if(this.nlrs == null) {
				nlrs = new Vector<MotifList>();
			}
			
			
			Collections.sort(preNlrs);
			
			
			int suffix = 1;
			for( Iterator<MotifList> iterator = preNlrs.iterator(); iterator.hasNext();) {
				MotifList motifList = iterator.next();
				motifList.sort();
				motifList.setMotifListName(dna_seq_id + "_nlr" + suffix );
				suffix++;
				nlrs.add(motifList);
			}
			
			
		}
	}
	
	/**
	 * 
	 * find all combinations of motifs that would be an NBARC motif seed; 
	 * The list of combinations are defined in AnnotatorSignatureDefinition.nbarcMotifIDCombinations 
	 * 
	 * 
	 * 
	 * @param motifs
	 * @param distanceWithinMotifCombination
	 * @return
	 */
	private Vector<MotifList> findNLRs_substep1_findSeeds(Vector<Motif> motifs ,int distanceWithinMotifCombination ){
		
		int counta = 0;
		
		Vector<MotifList> motifSeeds = new Vector<MotifList>();
		
		for(int i = 0; i< motifs.size(); i++) {
			Motif motif = motifs.get(i);
			
			
			
			String s = motif.getID();
			Vector<Motif> potentialMotifSeed = new Vector<Motif>();
			potentialMotifSeed.add(motif);
			
			for( int j = i+1; j<motifs.size(); j++) {
				Motif motif2 = motifs.get(j);
				
				if( Math.abs(motif.getDnaStart() - motif2.getDnaStart()) <= distanceWithinMotifCombination) {
					s = s + "," + motif2.getID();
					potentialMotifSeed.add(motif2);
					motif = motif2;
					
					
					
					if( this.annotatorSignatureDefinition.isSeed(s)) {
						counta++;
						
						MotifList list = new MotifList(potentialMotifSeed.firstElement().getDnaSequenceID()+"_nlr"+counta, potentialMotifSeed);
						motifSeeds.add(list);
						
						
						
						break;
					}
				}else {
					break;
				}
			}
			
			
		}
		
		Collections.sort(motifSeeds);
		return motifSeeds;
	}
	
	
	
	private Vector<MotifList> findNLRs_subsep2a_mergeSeeds(Vector<MotifList> motifSeeds, int distanceBetweenMotifCombinations){
		Vector<MotifList> returnLists = new Vector<MotifList>();
		
		
		
		while(motifSeeds.size()>0) {
			MotifList motifList = motifSeeds.firstElement();
			motifSeeds.remove(0);
			
			while(motifSeeds.size() >0) {
				MotifList list = motifSeeds.firstElement();
				if( motifList.canBeMergedWith(list, distanceBetweenMotifCombinations, annotatorSignatureDefinition)) {
					motifList.addMotifs(list);
					motifSeeds.remove(0);
					
				}else {
					break;
				}
			}
			
			returnLists.add(motifList);
			
			
			
			
		}
		
		
		return returnLists;
	}
	
	
	
	private Vector<MotifList> findNLRs_substep3_elongate(Vector<MotifList> mergedMotifSeeds, Vector<Motif> motifs, int distanceForElongating){
		Vector<MotifList> pre_nlrs = new Vector<MotifList>();
		
		HashSet<Motif> usedMotifs = new HashSet<Motif>(); //keep in mind motifs that have been used. don't add them again.
		
		
		for(Iterator<MotifList> iterator = mergedMotifSeeds.iterator(); iterator.hasNext();) {
			MotifList motifList =iterator.next();
			
			if( motifList.containsUsedMotif(usedMotifs)) {
				continue;
			}
			
			
			// go to first motif
			int index = 0;
			while(index < motifs.size()) { //go to first motif in list
				Motif motif = motifs.get(index);
				if( motif.equals(motifList.firstMotif())) {
					break;
				}
				index++;
			}
			
			// find matching motifs before
			int firstMotifIndex = index;
			index--;
			while( index >=0 ) {
				Motif motif = motifs.get(index);
				if( Math.abs(motifList.firstMotif().getDnaStart() - motif.getDnaStart()) <=distanceForElongating ) { //in distance
					if( annotatorSignatureDefinition.getMotifRank(motif) < annotatorSignatureDefinition.getMotifRank(motifList.firstMotif())) {
						
						if(!usedMotifs.contains(motif) &&( !motif.hasStop() || annotatorSignatureDefinition.isNBARC(motif))) {
							motifList.addMotif(motif);
							
						}
						
					}
					
				}else {
					break;
				}
				index--;
			}
			
			//go to last motif
			index = firstMotifIndex;
			while(index < motifs.size()) {
				Motif motif = motifs.get(index);
				if( motif.equals(motifList.lastMotif())) {
					break;
				}
				index++;
			}
			
			// find matching motifs after
			index++;
			while( index < motifs.size()) {
				Motif motif = motifs.get(index);
				
				if( Math.abs(motifList.lastMotif().getDnaStart() - motif.getDnaStart()) <=distanceForElongating) {//in distance
					int r1 = annotatorSignatureDefinition.getMotifRank(motif);
					int r2 = annotatorSignatureDefinition.getMotifRank(motifList.lastMotif());
					if( r1 > r2 || (r1==r2 && annotatorSignatureDefinition.isLRR(motif))) {
						if(!usedMotifs.contains(motif) &&( !motif.hasStop() || annotatorSignatureDefinition.isNBARC(motif))) {
							motifList.addMotif(motif);
							
						}
					}
				}else {
					break;
				}
				
				
				
				index++;
			}
			usedMotifs.addAll(motifList.getMotifs());
			pre_nlrs.add(motifList);
		}
		
		
		/*
		for( Iterator<MotifList> iterator1 = mergedMotifSeeds.iterator(); iterator1.hasNext();) {
			MotifList motifList = iterator1.next();
			
			
			//find range to next seeds
			
			long seed_start = Math.min(motifList.firstMotif().getDnaStart(), motifList.lastMotif().getDnaStart());
			long seed_end = Math.max(motifList.firstMotif().getDnaStart(), motifList.lastMotif().getDnaStart());
			
			long r1 = 0;
			long r2 = 1000000000;
			
			for( Iterator<MotifList> iterator2 = mergedMotifSeeds.iterator(); iterator2.hasNext();) {
				MotifList list = iterator2.next();
				long i1 = Math.max(list.firstMotif().getDnaStart(), list.lastMotif().getDnaStart());
				long i2 = Math.min(list.firstMotif().getDnaStart(), list.lastMotif().getDnaStart());
				
				if( i1 < seed_start && i1 >r1) {
					r1 = i1;
				}else if( i2 > seed_end && i2 < r2) {
					r2 = i2;;
				}
			}
			
			
			
			for( Iterator<Motif> iterator2 = motifs.iterator(); iterator2.hasNext();) {
				Motif motif = iterator2.next();
				
				if( motif.getDnaStart() < r1 ) {
					continue;
				}
				if( motif.getDnaStart() > r2) {
					continue;
				}
				
				int motifRank = annotatorSignatureDefinition.getMotifRank(motif);
				if(motif.hasStop()) { //if it has a stop codon and is none of the "essential" motifs, then ignore
					String category = annotatorSignatureDefinition.getCategory(motif);
					if( 	!category.equalsIgnoreCase("NBARC") &&
							!category.equalsIgnoreCase("LINKER") &&
							!category.equalsIgnoreCase("TIR") &&
							!category.equalsIgnoreCase("CC")) {
						continue;
					}
				}
				
				//System.out.print("\t" + motif.getID() + "\t" + motif.getDnaStart() + "\t" + motif.getStrand());
				
				
				if( motifList.isForwardStrand()) {
					if( !motif.isForwardStrand()) {
						continue;
					}
					
					if( motifList.firstMotif().getDnaStart() >  motif.getDnaStart()) { //before seed
				//		System.out.print("\tbefore seed");
						if( motifList.firstMotif().getDnaStart() - motif.getDnaStart() <= distanceForElongating) { //within range
						//	System.out.print("\twithin range");
							if(annotatorSignatureDefinition.getMotifRank(motifList.firstMotif() ) > motifRank) {
						//		System.out.print("\tADDED!");
								motifList.addMotif(motif);
							}
						}
						
					}else if (motifList.lastMotif().getDnaEnd() >= motif.getDnaStart()) { //within seed
						motifList.sort();
						
						for( int i  = 0; i< motifList.getMotifs().size()-1; i++) {
							Motif motif1 = motifList.getMotifs().get(i);
							Motif motif2 = motifList.getMotifs().get(i+1);
							
							if( motif1.getDnaStart() <= motif.getDnaStart() && motif.getDnaStart() <= motif2.getDnaStart()) {
								if( annotatorSignatureDefinition.getMotifRank(motif1) <= annotatorSignatureDefinition.getMotifRank(motif) &&
									annotatorSignatureDefinition.getMotifRank(motif) <= annotatorSignatureDefinition.getMotifRank(motif2)	) {
									motifList.addMotif(motif);
							//		System.out.print("\tADDED within!");
							//		System.out.println();
									break;
								}
							}
							
						}
						
						
					}else { //behind seed
						if( motif.getDnaStart() - motifList.lastMotif().getDnaStart() <= distanceForElongating) { //in distance
							if( annotatorSignatureDefinition.getMotifRank(motifList.lastMotif())  <= motifRank) { //rank fits. note that LRRs can be repeated. so this is a <= here.
								motifList.addMotif(motif);
							}
						}
					}
					
					
					
					
				}else {//reverse strand! According to sorting we will start with larger DNA positions
					if(motif.isForwardStrand()) {
						continue;
					}
					
					if( motif.getDnaStart() > motifList.firstMotif().getDnaStart() ) { //is before
				//		System.out.print("\tisbefore");
						if( motif.getDnaStart() - motifList.firstMotif().getDnaStart() <= distanceForElongating ) { //within distance
				//			System.out.println("\tis in range");
							if( motifRank < annotatorSignatureDefinition.getMotifRank(motifList.firstMotif())) {
								motifList.addMotif(motif);
							}
						}
						
					}else if ( motif.getDnaStart() >= motifList.lastMotif().getDnaStart()) { // motif is within seed
					//	System.out.print("\tis within");
						for( int i  = 0; i< motifList.getMotifs().size()-1; i++) {
							Motif motif1 = motifList.getMotifs().get(i);
							Motif motif2 = motifList.getMotifs().get(i+1);
							
							if( motif1.getDnaStart() >= motif.getDnaStart() && motif.getDnaStart() >= motif2.getDnaStart()) {
								if( annotatorSignatureDefinition.getMotifRank(motif1) <= annotatorSignatureDefinition.getMotifRank(motif) &&
									annotatorSignatureDefinition.getMotifRank(motif) <= annotatorSignatureDefinition.getMotifRank(motif2)	) {
									motifList.addMotif(motif);
									break;
								}
							}
							
						}
						
						
					}else { //motif is after seed
				//		System.out.print("\tis after");
						if( motifList.lastMotif().getDnaStart() - motif.getDnaStart() <= distanceForElongating) { //within distance
							if( annotatorSignatureDefinition.getMotifRank(motifList.lastMotif()) <= motifRank ) { // using <=; LRRs can have same rank
								motifList.addMotif(motif);
							}
							
						}else { //out of distance. we can stop going through the motif list
					//		System.out.println();
							break;
						}
						
						
					}
					
					
					
				}
				
				
				
				
				
				
				
				
			//	System.out.println();
				
			}
			
			motifList.removeRedundantMotifs();
			pre_nlrs.add(motifList);
			
		}
		
		try {
			writeBED(pre_nlrs, new File("/Users/steuernb/Documents/projects/Annotator2/test_data/debug/seeds4.bed"));
		}catch (IOException e) {}
		
		*/
		return pre_nlrs;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
		
	
	

	

	
	



	
	/**
	 * 
	 * Go through the motifs HashMap and remove duplicate motifs. Those will come because we find motifs in overlapping sequence fragments.  
	 * 
	 * 
	 */
	private void removeRedundantMotifs() {
		
		for(Iterator<String > iterator= this.motifs.keySet().iterator(); iterator.hasNext();) {
			
			String seq_id = iterator.next();
			Vector<Motif> mymotifs = this.motifs.get(seq_id);
			Collections.sort(mymotifs);
			
			
			int index = 0;
			
			while (index < mymotifs.size()){
				
				Motif motif1 = mymotifs.get(index);
				
				if( mymotifs.size()>index+1){
					Motif motif2 = mymotifs.get(index+1);
					
					if( motif1.getID().equalsIgnoreCase(motif2.getID()) 
							&& motif1.getDnaStart() == motif2.getDnaStart()) {
						mymotifs.remove(index+1);
						
					}
					
					
				}
				index++;
			}
			this.motifs.put(seq_id, mymotifs);
			
		}
		
		
		
		
		
		
		
	}
	
	

	/*
	
	private HashMap<String,Vector<Motif>> consolidateLists(HashMap<String, MotifList> motifLists){
		HashMap<String,Vector<Motif>> h = new HashMap<String,Vector<Motif>>();
		
		for(Iterator<String> iterator1 = motifLists.keySet().iterator(); iterator1.hasNext();){
			String key = iterator1.next();
			
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
			
			
			
			MotifList list = motifLists.get(key);
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
	
	
	
	
	
	
	
	*/
	
	
	
	
	
	
	
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
		
		for(Iterator <MotifList> iterator = this.nlrs.iterator(); iterator.hasNext();){
			
			MotifList list = iterator.next();
			
			if( !completeOnly || list.isCompleteNLR(this.annotatorSignatureDefinition)){
				out.write(list.getSequenceName() +"\tNLR_Annotator\tNBSLRR\t" );
				
				String strand = "+";
				long start = Math.min( Math.min(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.min( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
				long end = Math.max( Math.max(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.max( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
				
				start = start + 1; //GFF is 1-based!
				
				String domains = this.annotatorSignatureDefinition.getDomainString(list);
				
				if( !list.getMotifs().get(0).isForwardStrand()){
					strand = "-";
				}
				out.write(start + "\t" + end + "\t.\t"+strand+"\t.\tname="+ list.getMotifListName()+";nlrClass="+domains);
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
		
		for(Iterator <MotifList> iterator = this.nlrs.iterator(); iterator.hasNext();){
			
			MotifList list = iterator.next();
			
			Vector<Motif> v= new Vector<Motif>();
			v.addAll(list.getMotifs());
			
			String strand = "+";
			
			//positions in BED are supposed to be zero-based. 
			long start = Math.min( Math.min(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.min( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
			long end = Math.max( Math.max(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.max( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
			
			
			
			if( !list.getMotifs().get(0).isForwardStrand()){
				strand = "-";
				
			}
			
			
			String blockCount = v.size()+"";
			String blockStarts = "";
			String blockSizes = "";
			
			
			
			for(Iterator<Motif> iterator2 = v.iterator(); iterator2.hasNext();){
				Motif motif = iterator2.next();
				
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
			if(!list.isCompleteNLR(this.annotatorSignatureDefinition)){
				color = "255,128,0";//orange
			}
			
			if( list.hasStopCodon()){
				color = "255,0,0";
			}
			
			out.write(list.getSequenceName() +"\t" + start +"\t" + end +"\t"+list.getMotifListName()+"\t0\t"+ strand +"\t" + start +"\t" + end +"\t" + color + "\t" + blockCount + "\t" + blockSizes+"\t" +blockStarts );
			out.newLine();
		
			
		
		}
		
		
		
		out.close();
		
		
	}
	
	
	
	
	public static void writeBED(Vector<MotifList> nlrs, File outputFile)throws IOException {
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		
		out.write("#track name=\"Debug NLR-Annotator\"");
		out.newLine();
		out.write("#itemRgb=\"On\"");
		out.newLine();
		
		for(Iterator <MotifList> iterator = nlrs.iterator(); iterator.hasNext();){
			
			MotifList list = iterator.next();
			
			Vector<Motif> v= new Vector<Motif>();
			v.addAll(list.getMotifs());
			
			String strand = "+";
			
			//positions in BED are supposed to be zero-based. 
			long start = Math.min( Math.min(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.min( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
			long end = Math.max( Math.max(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.max( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
			
			
			
			if( !list.isForwardStrand()){
				strand = "-";
				
			}
			
			
			String blockCount = v.size()+"";
			String blockStarts = "";
			String blockSizes = "";
			
			
			
			for(Iterator<Motif> iterator2 = v.iterator(); iterator2.hasNext();){
				Motif motif = iterator2.next();
				
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
			
			
			
			
			String color = "0,0,255";//blue
			
			
			out.write(list.getSequenceName() +"\t" + start +"\t" + end +"\t"+list.getMotifListName()+"\t0\t"+ strand +"\t" + start +"\t" + end +"\t" + color + "\t" + blockCount + "\t" + blockSizes+"\t" +blockStarts );
			out.newLine();
		
			
		
		}
		
		
		
		out.close();
		
		
	}
	
	
	
	public void writeMotifGFF(Hashtable<String,Vector<Motif>> motifs, File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("##gff-version 2");
		out.newLine();
		out.write("##source-version NLR-Annotator " + getVersion() );
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
			Vector<Motif> v = motifs.get(key);
			for(Enumeration<Motif> myenum2 = v.elements();myenum2.hasMoreElements();){
				Motif motif = myenum2.nextElement();
				out.write(motif.getGFFString());
				out.newLine();
				
				
			}
			
		}
		
		
		out.close();
	}
	
	
	
	
	
	/**
	 * 
	 * Write a BED file that contains all motifs
	 * 
	 * @param outputFile
	 * 				The output BED file
	 * @param annotateSTOP
	 * 				if true, motifs that contain a stop codon are colored black
	 * @throws IOException
	 */
	public void writeMotifBED( File outputFile, boolean annotateSTOP)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("#track name=\"NLR_Motifs\"");
		out.newLine();
		out.write("#itemRgb=\"On\"");
		out.newLine();
		
		for(Iterator<String> iterator1 = motifs.keySet().iterator(); iterator1.hasNext();){
			String key = iterator1.next();
			Vector<Motif> v = motifs.get(key);
			for(Iterator<Motif> iterator2 = v.iterator();iterator2.hasNext();){
				Motif motif = iterator2.next();
				out.write(motif.getBEDString(this.annotatorSignatureDefinition));
				out.newLine();
				if( annotateSTOP && motif.getProteinSequence().contains("*")){
					out.write(motif.getBlackMotif());
					out.newLine();
				}
				
			}
			
		}
		
		
		out.close();
	}
	
	
	public void writeReportTxt( File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		for(Iterator <MotifList> iterator = this.nlrs.iterator(); iterator.hasNext();){
			
			MotifList list = iterator.next();
			long start = Math.min( Math.min(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.min( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
			long end = Math.max( Math.max(list.getMotifs().firstElement().getDnaStart(), list.getMotifs().firstElement().getDnaEnd()), Math.max( list.getMotifs().lastElement().getDnaStart(),list.getMotifs().lastElement().getDnaEnd())  ) ;
		   String strand = "+";
			if( !list.getMotifs().get(0).isForwardStrand()){
				strand = "-";
			}
			String domains = this.annotatorSignatureDefinition.getDomainString(list);
			//System.out.println(domains);
			
			
			out.write(list.getSequenceName()+"\t"+list.getMotifListName() +"\t"+domains+"\t" + start + "\t" + end +"\t" + strand + "\t" + list.getMotifListString());
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
		
		for(Iterator <MotifList> iterator = this.nlrs.iterator(); iterator.hasNext();){
			
			MotifList list = iterator.next();
			
			if(!list.isCompleteNLR(this.annotatorSignatureDefinition)){
				continue;
			}
			
			
			String sequence = "";
			
			int index =0;
			Motif motif = list.getMotifs().firstElement();
			while(!this.annotatorSignatureDefinition.isPloop(motif)){
				index++;
				motif = list.getMotifs().get(index);
			}
			
			sequence = motif.getProteinSequence() ;
			
			index++;
			
			
			String[] motifOrder = this.annotatorSignatureDefinition.getNbarcMotifOrder();
			for( int j = 1; j<motifOrder.length; j++){
				motif = list.getMotifs().get(index);
				if( !motif.getID().equalsIgnoreCase(motifOrder[j] )){
					String s = "";
					for( int i= 0; i<motifReference[j].length(); i++){
						s = s +"-";
					}
					sequence = sequence +s ;
					
				}else{
					sequence = sequence + motif.getProteinSequence();
					if (motifReference[j].length() != motif.getProteinSequence().length()){
						System.out.println( list.getMotifListName() + "\t" + motif.getProteinSequence() );
					}
					index ++;
					
					
				}
				
			}
			
			
			
			
			
			
			
			
			if( !sequence.contains("*")){
				out.write("#" + list.getMotifListName());
				out.newLine();
				out.write(sequence);
				out.newLine();
				out.newLine();
			}
			
			
			
			
			
			
			
			
			
		}
		
		
		out.close();
		
		
		
		
	}
	
	
	public void writeNbarcMultipleAlignmentFasta(File outputFasta,  boolean includeCed4)throws IOException{
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
		
		
		String[] motifOrder = this.annotatorSignatureDefinition.getNbarcMotifOrder();
		
		String[] motifReference = new String[motifOrder.length];
		for( int i = 0;i< motifOrder.length; i++) {
			motifReference[i] = this.annotatorSignatureDefinition.getDefaultSequence(motifOrder[i]);
		}
		
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFasta));
		
		if(includeCed4) {
			out.write(">NP_001021202.1");
			out.newLine();
	        out.write("FLHGRAGSGKSVIASQALSKS-----------------------------TLFVFDDVVQEETIRLRLRCLVTTRDVEISNAASQ"
	        		+ "--------------------------------------------------------------------------------");
	        out.newLine();
			
		}
		
       
		
		
		for(Iterator <MotifList> iterator = this.nlrs.iterator(); iterator.hasNext();){
			
			MotifList list = iterator.next();
		
			
			if(!list.isCompleteNLR(this.annotatorSignatureDefinition)){
				continue;
			}
			
			
			String sequence = "";
			
			int index =0;
			Motif motif = list.getMotifs().get(0);
			while(this.annotatorSignatureDefinition.isPloop(motif)){
				index++;
				motif = list.getMotifs().get(index);
			}
			
			sequence = motif.getProteinSequence() ;
			
			index++;
			
			
			
			for( int j = 1; j<motifOrder.length; j++){
				motif = list.getMotifs().get(index);
				if( !motif.getID().equalsIgnoreCase(motifOrder[j])){
					String s = "";
					for( int i= 0; i<motifReference[j].length(); i++){
						s = s +"-";
					}
					sequence = sequence +s ;
					
				}else{
					sequence = sequence + motif.getProteinSequence();
					if (motifReference[j].length() != motif.getProteinSequence().length()){
						System.out.println( list.getMotifListName() + "\t" + motif.getProteinSequence() );
					}
					index ++;
					
					
				}
				
			}
			
			
			
			
			
			
			
			
		
			
			
			//replace nucleotide gaps and stop codons by gaps in alignment
			sequence.replaceAll("\\*", "_");
			sequence.replaceAll("X", "_");
			out.write(new BioSequence( list.getMotifListName() , sequence ).getFastaString()   );
			
			
			
		}
		
		
		out.close();
		
		
		
		
	}
	
	
	
	
	
	
	
	public void writeNLRLoci(File genomeSequence, File outputFile, int flankingSequence)throws IOException{
		Hashtable<String, Vector<MotifList>> h = new Hashtable<String, Vector<MotifList>>();
		
		for(Enumeration<MotifList> myenum = nlrs.elements(); myenum.hasMoreElements();){
			MotifList list = myenum.nextElement();
			String scaffold = list.getSequenceName();
			if( !h.containsKey(scaffold)){
				h.put(scaffold,new Vector<MotifList>());
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
					
					Vector<MotifList> v = h.get(identifier);
					for( Enumeration<MotifList> myenum = v.elements(); myenum.hasMoreElements();){
						MotifList list = myenum.nextElement();
						
						
						int p1 = (int) list.getMotifs().firstElement().getDnaStart();
						int p2 = (int) list.getMotifs().firstElement().getDnaEnd();
						int p3 = (int) list.getMotifs().lastElement().getDnaStart();
						int p4 = (int) list.getMotifs().lastElement().getDnaEnd();
						
						int start = Math.max(0 , Math.min(Math.min(p1, p2), Math.min(p3, p4)) - flankingSequence);
						int end = Math.min(sequence.length(), Math.max(Math.max(p1, p2), Math.max(p3, p4)) + flankingSequence);
						
						BioSequence seq = new BioSequence(list.getMotifListName() , sequence.substring(start, end));
						String strand = "strand:+";
						if(!list.getMotifs().firstElement().isForwardStrand()){
							strand = "strand:-";
							seq.setSequence(seq.getReverseComplementarySequence());
						}
						seq.setDescription(identifier + " " + start + "-" + end + " " +strand + " " + list.getMotifListString());
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
			Vector<MotifList> v = h.get(identifier);
			for( Enumeration<MotifList> myenum = v.elements(); myenum.hasMoreElements();){
				MotifList list = myenum.nextElement();
				
				
				int p1 = (int) list.getMotifs().firstElement().getDnaStart();
				int p2 = (int) list.getMotifs().firstElement().getDnaEnd();
				int p3 = (int) list.getMotifs().lastElement().getDnaStart();
				int p4 = (int) list.getMotifs().lastElement().getDnaEnd();
				
				int start = Math.min(Math.min(p1, p2), Math.min(p3, p4)) - flankingSequence;
				int end = Math.max(Math.max(p1, p2), Math.max(p3, p4)) + flankingSequence;
				//System.out.println(list.getNlrID() +"\t"+sequence.length());
				BioSequence seq = new BioSequence("","");
				try{
					 seq = new BioSequence(list.getMotifListName(), sequence.substring(start, end));
				}catch(StringIndexOutOfBoundsException e){
					 seq = new BioSequence(list.getMotifListName(), sequence.substring(Math.max(0,start), Math.min(end, sequence.length())));
				}
				String strand = "strand:+";
				if(!list.getMotifs().firstElement().isForwardStrand()){
					strand = "strand:-";
					seq.setSequence(seq.getReverseComplementarySequence());
				}
				seq.setDescription(identifier + " " + start + "-" + end + " " +strand + " " + list.getMotifListString());
				out.write(seq.getFastaString(100));
				
			}
		}
		
		
		
		
		out.close();
		
		
	}
	
	
	
	
	
	
	
	public void exportMotifParserResult(File outputFile)throws IOException{
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		for(Iterator<String> iterator = motifs.keySet().iterator(); iterator.hasNext();) {
			String key = iterator.next();
			Vector<Motif> v = motifs.get(key);
			
			for(Iterator<Motif> iterator2 = v.iterator(); iterator2.hasNext();) {
				Motif motif = iterator2.next();
				out.write(motif.getExportString());
				out.newLine();
			}
			
			
		}
			
		out.close();
	}
	
	public void importMotifParserResult(File inputFile)throws IOException{
		
		if(this.motifs == null) {
			this.motifs = new HashMap<String, Vector<Motif>>();
		}
		
		BufferedReader in = new BufferedReader(new FileReader(inputFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			Motif motif = new Motif(inputline);
			String chr = motif.getDnaSequenceID();
			if(!this.motifs.containsKey(chr)) {
				this.motifs.put(chr,new Vector<Motif>());
			}
			this.motifs.get(chr).add(motif);
			
		}

		in.close();
		
		
		
	}
	
	
		
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
		public static void cli(String[] args)throws ExecutionException, InterruptedException{
		CLI cli = new CLI();
		
		
		
		String helpString =		"NLR-Annotator " + version + "\n"+
								
								
								"-i <inputFasta>\t\t\tInput file in fasta format. \n"+
								"-x <mot.txt>\t\t\tConfiguration file mot.txt. \n"+
								"-y <store.txt>\t\t\tConfiguration file store.txt. \n"+
								"-t <numthreads>\t\t\tNumber of parallel threads. Default is 1. \n"+
								"-n <numSeqsPerThead>\t\t\tNumber of Sequences processed in each thread. Default is 1000. \n"+
								"-c <output.tsv>\t\t\texport pre-calculated motifs into file. These can be imported again. . \n"+
								
								
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
			
			NLR_Annotator annotator;
			cli.parseOptions(args);
			
			int fragmentLength = 20000;
			int overlap = 5000;

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
		
			
			if( cli.hasOption("c") && !cli.hasOption("i") ) {
				File importFile = new File(cli.getArg("c"));
				annotator = new NLR_Annotator(importFile, fragmentLength, overlap, distanceWithinMotifCombination, distanceForElongating, distanceBetweenMotifCombinations);
				
				
				
			}else {
				
				if(!cli.hasOption("x")) {
					throw new CLIParseException("Configuration file mot.txt is missing. Use -x to specify");
				}
				
				File motFile = new File(cli.getArg("x"));
				
				if(!cli.hasOption("y")) {
					throw new CLIParseException("Configuration file store.txt is missing. Use -y to specify");
				}
				
				File storeFile = new File(cli.getArg("y"));

				if(!cli.hasOption("i")) {
					throw new CLIParseException("No input file specified. Use -i to specify the fasta file you want to analyze");
				}
				
				File inputFastaFile = new File(cli.getArg("i"));
				
				int numThreads = 1;
				int numSequencesPerThread = 1000;
				
				if( cli.hasOption("t")) {
					
				
					if(cli.hasOption("n")) {
						numSequencesPerThread = Integer.parseInt(cli.getArg("n"));
					}
					
					 numThreads = Integer.parseInt(cli.getArg("t"));
				
					
				}
				
					annotator = new NLR_Annotator(inputFastaFile,  motFile,  storeFile,  fragmentLength,  overlap,  numSequencesPerThread,  numThreads, distanceWithinMotifCombination, distanceForElongating, distanceBetweenMotifCombinations);
				
				
			}
			
			
		
			
			
			if( cli.hasOption("c")){
				annotator.exportMotifParserResult(new File(cli.getArg("c")));
			}
			
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
				annotator.writeMotifBED( new File(cli.getArg("m")), false);
			}
			
			if( cli.hasOption("a")){
				annotator.writeNbarcMultipleAlignmentFasta(new File(cli.getArg("a")),  true);
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
		
		
		
	}

	
	
	public static void main(String[] args) {
		try {
			
			//run711N16();
			//runTAIRchr1();
			//runTAIRchr1_Import();
			//runTAIRchr1_parallel();
			cli(args);
			//runTestFastaOut();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
	}
	
	
	
	
	
	
	
	
	public static void run711N16() throws IOException{
		File motFile = new File("/Users/steuernb/workspace/NLR-Annotator2/src/motifParser/mot.txt");
		File storeFile = new File("/Users/steuernb/workspace/NLR-Annotator2/src/motifParser/store.txt");

		File inputFastaFile = new File("/Users/steuernb/Documents/projects/Annotator2/test_data/711N16/BAC711N16.fasta");
		
		int fragmentLength = 20000;
		int overlap = 5000;
		int distanceWithinMotifCombination=500;
		int distanceForElongating=2500;
		int distanceBetweenMotifCombinations=10000;
		
		NLR_Annotator annotator = new NLR_Annotator(inputFastaFile, motFile, storeFile, fragmentLength, overlap, distanceWithinMotifCombination, distanceForElongating, distanceBetweenMotifCombinations);
		annotator.exportMotifParserResult(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/BAC711N16.MotifParser.tsv"));
		annotator.writeNlrBED(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/BAC711N16.nlr.bed"));
		annotator.writeMotifBED(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/BAC711N16.motifs.bed"), true);
		annotator.writeReportTxt(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/BAC711N16.nlr.txt"));
	}
	public static void runTAIRchr1_Import()throws IOException{
		File importFile = new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.MotifParser.tsv");
		
		int fragmentLength = 20000;
		int overlap = 5000;
		int distanceWithinMotifCombination=500;
		int distanceForElongating=2500;
		int distanceBetweenMotifCombinations=10000;
		
		NLR_Annotator annotator = new NLR_Annotator(importFile, fragmentLength, overlap, distanceWithinMotifCombination, distanceForElongating, distanceBetweenMotifCombinations);
		
		annotator.writeNlrBED(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.nlr.bed"));
		annotator.writeMotifBED(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.motifs.bed"), true);
		annotator.writeReportTxt(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.nlr.txt"));
	}
	public static void runTAIRchr1()throws IOException{
		File motFile = new File("/Users/steuernb/workspace/NLR-Annotator2/src/motifsearch/mot.txt");
		File storeFile = new File("/Users/steuernb/workspace/NLR-Annotator2/src/motifsearch/store.txt");

		File inputFastaFile = new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.fas");
		
		int fragmentLength = 20000;
		int overlap = 5000;
		int distanceWithinMotifCombination=500;
		int distanceForElongating=2500;
		int distanceBetweenMotifCombinations=10000;
		
		NLR_Annotator annotator = new NLR_Annotator(inputFastaFile, motFile, storeFile, fragmentLength, overlap, distanceWithinMotifCombination, distanceForElongating, distanceBetweenMotifCombinations);
		annotator.exportMotifParserResult(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.MotifParser.tsv"));
		annotator.writeNlrBED(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.nlr.bed"));
		annotator.writeMotifBED(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.motifs.bed"), true);
		annotator.writeReportTxt(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.nlr.txt"));
	}
	public static void runTAIRchr1_parallel()throws IOException, ExecutionException,InterruptedException{
		File motFile = new File("/Users/steuernb/workspace/NLR-Annotator2/src/motifsearch/mot.txt");
		File storeFile = new File("/Users/steuernb/workspace/NLR-Annotator2/src/motifsearch/store.txt");

		File inputFastaFile = new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.fas");
		
		int fragmentLength = 20000;
		int overlap = 5000;
		int distanceWithinMotifCombination=500;
		int distanceForElongating=2500;
		int distanceBetweenMotifCombinations=10000;
		int numThreads = 2;
		int numSequencesPerThread = 100;
		
		NLR_Annotator annotator = new NLR_Annotator( inputFastaFile,  motFile,  storeFile,  fragmentLength,  overlap,  numSequencesPerThread,  numThreads, distanceWithinMotifCombination, distanceForElongating, distanceBetweenMotifCombinations);
		annotator.exportMotifParserResult(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.MotifParser_p.tsv"));
		annotator.writeNlrBED(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.nlr_p.bed"));
		annotator.writeMotifBED(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.motifs_p.bed"), true);
		annotator.writeReportTxt(new File("/Users/steuernb/Documents/projects/Annotator2/test_data/TAIR10_chr1.nlr_p.txt"));
	}
	
	
	
	
}
