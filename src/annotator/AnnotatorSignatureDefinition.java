package annotator;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import core.Motif;
import core.MotifList;


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
public class AnnotatorSignatureDefinition {

	
	
	HashMap<String, Integer> motifRanks;
	
	Vector<String[]> motifIDCombinations;
	HashSet<String> motifIDCombinationHash;
	
	HashMap<String, String> motifCategories; //keys are motif_ids, values are a category, such as "NBARC", "LRR"
	
	HashMap<String, String> defaultMotifSequences;
	
	HashMap<String, int[]> motifRgbColors;
	
	String ploopMotif;
	
	
	public AnnotatorSignatureDefinition() {
		
		this.loadDefaultMotifRanks();
		this.loadDefaultMotifIDCombinations();
		this.loadDefaultMotifCategories();
		this.loadDefaultDefaultMotifSequences();
		this.loadDefaultMotifRgbColors();
		
		ploopMotif = "motif_1";
		
	}
	
	
	
	public String getMotifRgbColor(String motif_id) {
		int[] rgb = this.motifRgbColors.get(motif_id);
		return rgb[0] + "," + rgb[1] + "," + rgb[2];
	}
	
	
	
	public int getMotifRank(String motif_id) {
		
		return this.motifRanks.get(motif_id);
		
	}
	
	
	public int getMotifRank(Motif motif) {
		
		return this.motifRanks.get(motif.getID());
		
	}
	
	public boolean isSeed(String s) {
		if(this.motifIDCombinationHash.contains(s)) {
			return true;
		}else {
			return false;
		}
	}
	
	
	
	
	
	
	
	
	public boolean isConsistent(Vector<Motif> motifList) {
		Collections.sort(motifList);
		
		int lastRank = getMotifRank(motifList.firstElement());
		for( int i = 1; i< motifList.size(); i++) {
			Motif motif = motifList.get(i);
			int rank = getMotifRank(motif);
			if(rank >=lastRank) {
				lastRank = rank;
			}else {
				return false;
			}
		}
		
		return true;
	}
	
	
	/**
	 * 
	 * check if a motif List contains an NB-ARC domain. If it contains at least 3 NB-ARC-associated motifs, this will return true;
	 * 
	 * @param motifs
	 * @return
	 */
	public boolean hasNBARC(Vector<Motif> motifList) {
		
		int nbarcCounta = 0;
		
		for(Iterator<Motif> iterator= motifList.iterator(); iterator.hasNext();){
			Motif motif = iterator.next();
			if( this.getCategory(motif).equalsIgnoreCase("NBARC")) {
				nbarcCounta++;
			}
			if(nbarcCounta >=3) {
				return true;
			}
		}
		
		return false;
	}
	
	
	
	
	public String getDomainString(MotifList motifList) {
		String s = "";
		String currentCategory = "";
		for(Iterator<Motif> iterator = motifList.getMotifs().iterator(); iterator.hasNext();) {
			Motif motif = iterator.next();
			String category = this.getCategory(motif);
			
			if(!category.equalsIgnoreCase(currentCategory) && !category.equalsIgnoreCase("NA") && !category.equalsIgnoreCase("LINKER")) {
				s = s + "-" + category;
			}
			currentCategory = category;
		}
		return s.substring(1);
		
	}
	
	
	public String getCategory(String motif_id) {
		return this.motifCategories.get(motif_id);
	}
	
	public String getCategory(Motif motif) {
		return this.motifCategories.get(motif.getID());
	}
	
	
	public boolean isLRR(Motif motif) {
		if(this.motifCategories.get(motif.getID()).equalsIgnoreCase("LRR") ) {
			return true;
		}else {
			return false;
		}
	}
	
	
	public boolean isPloop(Motif motif) {
		
		if( ploopMotif == null) {
			this.guessPloopMotif();
		}	
		
		if(this.ploopMotif.equalsIgnoreCase(motif.getID())) {
			return true;
		}else {
			return false;
		}
	
		
	}
	
	public boolean isNBARC(Motif motif) {
		if(this.motifCategories.get(motif.getID()).equalsIgnoreCase("NBARC")) {
			return true;
		}else {
			return false;
		}
	}
	
	
	/**
	 * 
	 * Get the consensus amino acid sequence for a motif. This might be interesting for pylogenetics. Currently not in use.
	 * 
	 * 
	 * @param motif_id
	 * 		ID of the motif we want to get the sequence for.
	 * @return
	 * 		The consensus amino acid sequence for this motif
	 */
	public String getDefaultSequence(String motif_id) {
		
		return this.defaultMotifSequences.get(motif_id);
		
	}
	
	
	/**
	 * 
	 * Assume that the NBARC motif with lowest rank is the ploop. 
	 * 
	 * 
	 * 
	 */
	private void guessPloopMotif() {
		
		String motifID = null;
		int rank = 1000000;
		
		for(Iterator<String> iterator = this.motifRanks.keySet().iterator(); iterator.hasNext();) {
			
			String current_motifID = iterator.next();
			int currentRank = this.motifRanks.get(current_motifID);
			
			if(this.motifCategories.get(current_motifID).equalsIgnoreCase("NBARC") && currentRank < rank) {
				rank = currentRank;
				motifID = current_motifID;
						
			}
			
		}
		
		this.ploopMotif = motifID;
		
	}
	
	
	
	
	
	
	private void loadDefaultMotifIDCombinations() {
		
		this.motifIDCombinations = new Vector<String[]>();
		
		this.motifIDCombinations.add(new String[] {"motif_1", "motif_6", "motif_4"});
		this.motifIDCombinations.add(new String[] {"motif_6", "motif_4", "motif_5"});
		this.motifIDCombinations.add(new String[] {"motif_4", "motif_5", "motif_10"});
		this.motifIDCombinations.add(new String[] {"motif_5", "motif_10", "motif_3"});
		this.motifIDCombinations.add(new String[] {"motif_10", "motif_3", "motif_12"});
		this.motifIDCombinations.add(new String[] {"motif_3", "motif_12", "motif_2"});
		this.motifIDCombinations.add(new String[] {"motif_1", "motif_4", "motif_5"});
		this.motifIDCombinations.add(new String[] {"motif_12", "motif_2", "motif_8"});
		this.motifIDCombinations.add(new String[] {"motif_2", "motif_8", "motif_7"});
		this.motifIDCombinations.add(new String[] {"motif_18", "motif_15", "motif_13"});
		
		/*
		this.motifIDCombinations.add(new String[] {"motif_17", "motif_16"});
		this.motifIDCombinations.add(new String[] {"motif_16", "motif_1"});
		this.motifIDCombinations.add(new String[] {"motif_1", "motif_6"});
		this.motifIDCombinations.add(new String[] {"motif_18", "motif_15"});
		this.motifIDCombinations.add(new String[] {"motif_15", "motif_13"});
		this.motifIDCombinations.add(new String[] {"motif_13", "motif_1"});
		*/
		this.motifIDCombinations.add(new String[] {"motif_1", "motif_6"});
		this.ploopMotif="motif_1";
		
		
		
		this.motifIDCombinationHash = new HashSet<String>();
		for(Iterator<String[]> iterator = this.motifIDCombinations.iterator(); iterator.hasNext();) {
			String[] a = iterator.next();
			String s = a[0];
			for( int i = 1; i< a.length; i++) {
				s = s + "," +a[i];
			}
			this.motifIDCombinationHash.add(s);
		}
		
		
	}
	
	private void loadDefaultMotifRanks(){
		
		if( this.motifRanks == null) {
			this.motifRanks = new HashMap<String, Integer>();
		}
		
		this.motifRanks.put("motif_1",4);
		this.motifRanks.put("motif_2",11);
		this.motifRanks.put("motif_3",9);
		this.motifRanks.put("motif_4",6);
		this.motifRanks.put("motif_5",7);
		this.motifRanks.put("motif_6",5);
		this.motifRanks.put("motif_7",13);
		this.motifRanks.put("motif_8",12);
		this.motifRanks.put("motif_9",14);
		this.motifRanks.put("motif_10",8);
		this.motifRanks.put("motif_11",14);
		this.motifRanks.put("motif_12",10);
		this.motifRanks.put("motif_13",3);
		this.motifRanks.put("motif_14",3);
		this.motifRanks.put("motif_15",2);
		this.motifRanks.put("motif_16",2);
		this.motifRanks.put("motif_17",1);
		this.motifRanks.put("motif_18",1);
		this.motifRanks.put("motif_19",14);
		this.motifRanks.put("motif_20",15);
	
	}
	
	
	private void loadDefaultMotifCategories() {
		
		if(this.motifCategories == null) {
			this.motifCategories = new HashMap<String, String>();
		}
		
		this.motifCategories.put("motif_1", "NBARC");
		this.motifCategories.put("motif_2", "NBARC");
		this.motifCategories.put("motif_3", "NBARC");
		this.motifCategories.put("motif_4", "NBARC");
		this.motifCategories.put("motif_5", "NBARC");
		this.motifCategories.put("motif_6", "NBARC");
		this.motifCategories.put("motif_7", "LINKER");
		this.motifCategories.put("motif_8", "LINKER");
		this.motifCategories.put("motif_9", "LRR");
		this.motifCategories.put("motif_10", "NBARC");
		this.motifCategories.put("motif_11", "LRR");
		this.motifCategories.put("motif_12", "NBARC");
		this.motifCategories.put("motif_13", "TIR");
		this.motifCategories.put("motif_14", "NA");
		this.motifCategories.put("motif_15", "TIR");
		this.motifCategories.put("motif_16", "CC");
		this.motifCategories.put("motif_17", "CC");
		this.motifCategories.put("motif_18", "TIR");
		this.motifCategories.put("motif_19", "LRR");
		this.motifCategories.put("motif_20", "NA");
		
		
	}
	
	
	public void loadDefaultDefaultMotifSequences() {
		
		if(this.defaultMotifSequences == null) {
			this.defaultMotifSequences = new HashMap<String, String>();
		}
		
		this.defaultMotifSequences.put("motif_1", "PIWGMGGVGKTTLARAVYNDP");
		this.defaultMotifSequences.put("motif_6", "HFDCRAWVCVSQQYDMKKVLRDIIQQVGG");
		this.defaultMotifSequences.put("motif_4", "YLVVLDDVWDTDQWD");
		this.defaultMotifSequences.put("motif_5", "NGSRIIITTRNKHVANYMCT");
		this.defaultMotifSequences.put("motif_10", "LSHEESWQLFHQHAF");
		this.defaultMotifSequences.put("motif_3", "CGGLPLAIKVWGGMLAGKQKT");
		this.defaultMotifSequences.put("motif_12", "IMPVLRLSYHHLPYH");
		this.defaultMotifSequences.put("motif_2", "LKPCFLYCAIFPEDYMIDKNKLIWLWMAE");


	}
	
	
	
	public void loadDefaultMotifRgbColors() {
		
		if(this.motifRgbColors == null) {
			this.motifRgbColors = new HashMap<String, int[]>();
		}
		
		this.motifRgbColors.put("motif_1"  ,new int[]{0,255,255});
		this.motifRgbColors.put("motif_2"  ,new int[]{0,0,255}); 
		this.motifRgbColors.put("motif_3"  ,new int[]{255,0,0}); 
		this.motifRgbColors.put("motif_4"  ,new int[]{255,0,255});
		this.motifRgbColors.put("motif_5"  ,new int[]{255,255,0});
		this.motifRgbColors.put("motif_6"  ,new int[]{0,255,0}); 
		this.motifRgbColors.put("motif_7"  ,new int[]{0,128,128});
		this.motifRgbColors.put("motif_8"  ,new int[]{68,68,68}); 
		this.motifRgbColors.put("motif_9"  ,new int[]{0,128,0}); 
		this.motifRgbColors.put("motif_10" ,new int[]{192,192,192});
		this.motifRgbColors.put("motif_11" ,new int[]{128,0,128});
		this.motifRgbColors.put("motif_12" ,new int[]{128,128,0});
		this.motifRgbColors.put("motif_13" ,new int[]{0,0,128}); 
		this.motifRgbColors.put("motif_14" ,new int[]{128,0,0}); 
		this.motifRgbColors.put("motif_15" ,new int[]{255,255,255});
		this.motifRgbColors.put("motif_16" ,new int[]{0,255,255});
		this.motifRgbColors.put("motif_17" ,new int[]{0,0,255}); 
		this.motifRgbColors.put("motif_18" ,new int[]{255,0,0}); 
		this.motifRgbColors.put("motif_19" ,new int[]{255,0,255});
		this.motifRgbColors.put("motif_20" ,new int[]{255,255,0});

		
		
	}
	
	
	/**
	 * 
	 * check if a rank is the same as an LRR motif.
	 * 
	 * @param inputMotifRank
	 * @return
	 */
	public boolean isLrrMotifRank(int inputMotifRank) {
		for(Iterator<String> iterator = this.motifRanks.keySet().iterator(); iterator.hasNext();) {
			String key = iterator.next();
			int rank = this.motifRanks.get(key);
			if(rank == inputMotifRank) {
				if(this.motifCategories.get(key).equalsIgnoreCase("LRR")) {
					return true;
				}else {
					return false;
				}
			}
		}
		
		return false;
	}
	
	
	public String[] getNbarcMotifOrder() {
		
		Vector<String> v = new Vector<String>();
		for(Iterator<String> iterator = this.motifCategories.keySet().iterator(); iterator.hasNext();) {
			String motif_id = iterator.next();
			if(this.motifCategories.get(motif_id).equalsIgnoreCase("NBARC")) {
				v.add(motif_id);
			}
		}
		
		Collections.sort(v, new Comparator<String>() {
												public int compare(String s1, String s2) {
													
													if(motifRanks.get(s1) < motifRanks.get(s2)) {
														return -1;
													}
													if(motifRanks.get(s1) > motifRanks.get(s2)) {
														return 1;
													}
													
													return 0;}} );
		
		
		String[] nbarcMotifOrder = new String[v.size()];
		for( int i = 0; i< v.size(); i++) {
			nbarcMotifOrder[i] = v.get(i);
		}
		
		return nbarcMotifOrder;
	}
	
	
	
}
