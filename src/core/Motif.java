package core;

import annotator.AnnotatorSignatureDefinition;


/**
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
 *
 *
 *
 * Motif holds a motif that has been found in a sequence. The id is the name of the motif, e.g. motif_1 from the meme.xml.
 * 
 * This data structure can also hold information about DNA sequence if a motif is found in translated DNA. This object will be initialized with
 * the protein information but can be extended with DNA information by the setDNA() method. 
 *
 *
 *
 *
 *
 *
 */
public class Motif implements Comparable<Motif> {

	
	String id; // name of the motif, e.g. motif_1 from the meme.xml
	String protein_sequence_id; // identifier of the protein sequence
	int position; //start position in the protein sequence
	
	String protein_sequence; //the amino acid sequence of the motif
	double pvalue; // a p-value
	double score; // a score
	
	
	String dna_sequence_id; // in case this motif is from a translated DNA fragment, this is the identifier of the DNA sequence
	long dnaStart; // leftmost position in the DNA seqeuence
	long dnaEnd; // rightmost position in the DNA sequence
	boolean forwardStrand; // true if it is the forward strand, false if reverse
	int frame; // the frame. this is 0, 1 or 2.
	boolean dnaParametersSet; // if the DNA parameters are set, this should be set to true. 
	
	
	
	
	
	/**
	 * 
	 *  Constructor for a Motif object
	 * 
	 * 
	 * 
	 * @param id
	 * 			motif name
	 * @param protein_sequence_id
	 * 			identifier of the protein sequence
	 * @param position
	 * 			leftmost position in the protein sequence
	 * @param protein_sequence
	 * 			amino acid sequence of the motif (the version found in this sequence)
	 * @param score
	 * 			a score
	 */
	public Motif(String id, String protein_sequence_id, int position, String protein_sequence, double pvalue) {
		super();
		this.id = id;
		this.protein_sequence_id = protein_sequence_id;
		this.position = position;
		this.protein_sequence = protein_sequence;
		this.pvalue = pvalue;
		dnaParametersSet = false;
	}
	
	
	
	public Motif(String importString) {
		super();
		
		String[] split = importString.split("\t");
		this.id = split[0];
		this.protein_sequence_id = split[1];
		this.position = Integer.parseInt(split[2]);
		this.protein_sequence = split[3];
		this.pvalue = Double.parseDouble(split[4]);
		this.score = Double.parseDouble(split[5]);
		
		if( split.length>6 && split[6].length()>0) {
			this.dna_sequence_id = split[6];
			this.dnaStart = Integer.parseInt(split[7]);
			this.dnaEnd = Integer.parseInt(split[8]);
			this.forwardStrand = true;
			if(split[9].equalsIgnoreCase("-")) {
				forwardStrand = false;
			}
			this.frame = Integer.parseInt(split[10]);
			this.dnaParametersSet = true;
		}else {
			this.dnaParametersSet = false;
		}
		
			
	}
	
	
	public boolean isForwardStrand() {
		return this.forwardStrand;
	}
	
	public long getDnaStart() {
		return this.dnaStart;
	}
	
	public long getDnaEnd() {
		return this.dnaEnd;
	}
	
	public void setPvalue(double pvalue) {
		this.pvalue = pvalue;
	}
	
	@Override
	public boolean equals(Object motifObject) {
		
		if( motifObject.getClass() != this.getClass()) {
			return false;
		}
		
		Motif motif = (Motif) motifObject;
		if (this.id.equalsIgnoreCase(motif.getID()) &&  this.dnaStart == motif.getDnaStart()) {
			return true;
		}else {
			return false;
		}
	}
	
	
	@Override
	public int hashCode() {
		
		String s = dna_sequence_id + dnaStart+id;
		
		return s.hashCode();
		
	}
	
	
	public String getStrand(){
		if(this.isForwardStrand()){
			return "+";
		}else{
			return "-";
		}
	}
	
	
	
	public void setDNA(String dna_sequence_id, long offset, int fragmentLength, int frame, boolean forwardStrand) {
		
		this.dna_sequence_id = dna_sequence_id;
		
		this.forwardStrand = forwardStrand;
		this.frame = frame;
		
		if(forwardStrand) {
			dnaStart = (this.position -1)*3 + frame + offset;
			dnaEnd = (this.position + protein_sequence.length() -1) * 3 + frame + offset;
		}else {
			
			dnaStart = offset + fragmentLength - (( position + protein_sequence.length()-1) *3 + frame);
			dnaEnd = offset + fragmentLength - (( position -1) *3 + frame);
		}
		
		
		dnaParametersSet = true;
	}
	
	
	
	
	public String getID() {
		return this.id;
	}
	
	
	
	
	
	
	public int getStart() {
		return this.position;
		
	}
	public String getProteinSequenceID() {
		return this.protein_sequence_id;
	}
	
	public String getDnaSequenceID() {
		return this.dna_sequence_id;
	}
	
	
	public String getProteinSequence() {
		return this.protein_sequence;
	}
	
	
	
	public boolean meetsThreshold(double pvalueThreshold) {
		if( pvalue <= pvalueThreshold) {
			return true;
		}else {
			return false;
		}
	}
	
	public int compareTo(Motif motif2) {
		
		if(dnaParametersSet) {
			if( this.dna_sequence_id.equalsIgnoreCase(motif2.dna_sequence_id)) {
				
				if( this.forwardStrand && !motif2.forwardStrand) {
					return -1;
				}else
				if( !this.forwardStrand && motif2.forwardStrand) {
					return 1;
				}else		
				if( this.forwardStrand) { //both forward strand
					if(this.dnaStart < motif2.dnaStart) {
						return -1;
					}else if(this.dnaStart > motif2.dnaStart) {
						return 1;
					}else {
						return 0;
					}
					
				}else { //both reverse strand
					if(this.dnaStart > motif2.dnaStart) {
						return -1;
					}else if(this.dnaStart < motif2.dnaStart) {
						return 1;
					}else {
						return 0;
					}
				}
				
				
				
			}else {
				return this.protein_sequence_id.compareTo(motif2.protein_sequence_id);
			}
		}else {
			if( this.protein_sequence_id.equalsIgnoreCase(motif2.protein_sequence_id)) {
				
				if(this.position < motif2.position) {
					return -1;
				}else if(this.position > motif2.position) {
					return 1;
				}else {
					return 0;
				}
				
				
			}else {
				return this.protein_sequence_id.compareTo(motif2.protein_sequence_id);
			}
		}
		
		
		
	}
	
	
	public boolean hasStop() {
		if(this.getProteinSequence().contains("*")) {
			return true;
		}else {
			return false;
		}
	}
	
	
	
	public String getMotifString() {
		String strand = "+";
		if(!this.isForwardStrand()) {
			strand = "-";
		}
		
		String motifString = this.protein_sequence_id +"\t"+ id + "\t" +  this.dnaStart + "\t" + this.dnaEnd + "\t" + strand + "\t" + frame;
		return motifString;
	}
	
	public String getBEDString(AnnotatorSignatureDefinition annotatorSignatureDefinition){
		
		return this.dna_sequence_id +"\t" + this.getDnaStart() +"\t" + this.getDnaEnd() +"\t"+this.getID()+"\t"+this.pvalue+"\t"+ this.getStrand() +"\t" + this.getDnaStart() +"\t" + this.getDnaEnd() +"\t" + annotatorSignatureDefinition.getMotifRgbColor(id);
		
	}
	public String getBlackMotif(){
		return this.dna_sequence_id +"\t" + this.getDnaStart() +"\t" + this.getDnaEnd() +"\tSTOP\t0\t"+ this.getStrand() +"\t" + this.getDnaStart() +"\t" + this.getDnaEnd() +"\t0,0,0" ;
	}
	
	
	public String getGFFString(){
		//eqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute
		return this.dna_sequence_id + "\tNLR_Annotator\tMast_Motif\t" + (this.getDnaStart()+1) +"\t" + (this.getDnaEnd()) +"\t" + this.score + "\t" +this.getStrand() + "\t" + this.frame +"\tname=" + this.getID(); 
		
	}

	public String getExportString() {
		String s = id + "\t" + protein_sequence_id + "\t" + position + "\t" + protein_sequence + "\t" + pvalue + "\t" + score;
		
		if(dnaParametersSet) {
			s = s + "\t" + dna_sequence_id + "\t" + dnaStart + "\t" + dnaEnd + "\t" + this.getStrand() + "\t" + frame;
		}else {
			s = s + "\t\t\t\t\t";
		}
		return s;
	}
}
