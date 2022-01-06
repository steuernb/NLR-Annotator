package core;

import java.util.Iterator;
import java.util.Vector;



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
public class SignatureDefinition {

	
	Vector<String[]> signatures;
	
	
	public SignatureDefinition(Vector<String[]> signatures) {
		this.signatures = signatures;
	}
	public SignatureDefinition() {
		this.signatures = new Vector<String[]>();
		this.loadDefaultSignature();
	}
	
	
	public boolean hasSignature(MotifList motifList) {
		
		for(Iterator<String[]> iterator = this.signatures.iterator(); iterator.hasNext();) {
			String[] signature = iterator.next();
			
			
			motifList.sort();
			
			Vector<Motif> inputList = motifList.motifs;
			
			
			for( int i = 0; i< inputList.size(); i++) {
				
				
				boolean hasSignature = true;
				for( int j = 0; j<signature.length; j++) {
					
					try {
						if( !signature[j].equalsIgnoreCase( inputList.get(i+j).getID()) ){
							hasSignature = false;
							break;
						}
					}catch(ArrayIndexOutOfBoundsException e) {
						hasSignature = false;
						break;
					}
					
					
					
				}
				if(hasSignature) {
					return true;
				}
			}
		}
		
		
		return false;
	}
	
	
	public void loadDefaultSignature() {
		Vector<String[]> motivCombinations = new Vector<String[]>();
		String[] b = {"motif_17", "motif_16"}; motivCombinations.add(b);
		b = new String[2]; b[0] = "motif_1"; b[1] = "motif_6";  motivCombinations.add(b);
		b = new String[3]; b[0] = "motif_1"; b[1] = "motif_6"; b[2] = "motif_4"; motivCombinations.add(b);
		b = new String[3]; b[0] = "motif_6"; b[1] = "motif_4"; b[2] = "motif_5"; motivCombinations.add(b);
		b = new String[3]; b[0] = "motif_4"; b[1] = "motif_5"; b[2] = "motif_10";motivCombinations.add(b);
		b = new String[3]; b[0] = "motif_5"; b[1] = "motif_10"; b[2] = "motif_3";motivCombinations.add(b);
		b = new String[3]; b[0] = "motif_10"; b[1] = "motif_3"; b[2] = "motif_12"; motivCombinations.add(b);
		b = new String[3]; b[0] = "motif_3"; b[1] = "motif_12"; b[2] = "motif_2"; motivCombinations.add(b);
		b = new String[3]; b[0] = "motif_12"; b[1] = "motif_2"; b[2] = "motif_8"; motivCombinations.add(b);
		b = new String[3]; b[0] = "motif_2"; b[1] = "motif_8"; b[2] = "motif_7"; motivCombinations.add(b);
		b = new String[3]; b[0] = "motif_8"; b[1] = "motif_7"; b[2] = "motif_9"; motivCombinations.add(b);
		b = new String[3]; b[0] = "motif_7"; b[1] = "motif_9"; b[2] = "motif_11";motivCombinations.add(b);
		b = new String[2]; b[0] = "motif_9"; b[1] = "motif_11";  motivCombinations.add(b);
		b = new String[2]; b[0] = "motif_11"; b[1] = "motif_9";  motivCombinations.add(b);
		b = new String[2]; b[0] = "motif_18"; b[1] = "motif_15"; motivCombinations.add(b);
		b = new String[2]; b[0] = "motif_15"; b[1] = "motif_13"; motivCombinations.add(b);
		b = new String[2]; b[0] = "motif_13"; b[1] = "motif_1"; motivCombinations.add(b);
		b = new String[3]; b[0] = "motif_1"; b[1] = "motif_4"; b[2]="motif_5" ; motivCombinations.add(b);
		
		this.signatures = motivCombinations;

		
	}
	
}
