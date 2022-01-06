package motifParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
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
 *  @author Jitender Cheema
 * 
 */
public class MotifDefinition {

	HashMap<String, int[][]> pwm;
	HashMap<String, double[]> cdf;
	
	HashMap<String, Integer> motifLengths;
	Vector<String> motifNames;
	
	
	public MotifDefinition(File pwmFile, File cdfFile) throws IOException{
		
		this.loadMotifLengths(pwmFile);
		this.loadPWM(pwmFile);
		this.loadCDF(cdfFile);
		
		
		
		
		
		
	}
	
	
	
	public int getLength(String motif) {
		
		return motifLengths.get(motif);
	}
	
	
	
	
	/**
	 * 
	 * 
	 * get the maximum length a motif can have. This value is important for the MotifParser, which cannot handle sequences shorter than max motif length.
	 * 
	 * 
	 * @return
	 * 		maximum motif length
	 */
	public int getMaximumMotifLength() {
		int length = 0;
		
		for(Iterator<String> iterator = motifLengths.keySet().iterator(); iterator.hasNext();) {
			String motifName = iterator.next();
			if( motifLengths.get(motifName) > length) {
				length = motifLengths.get(motifName);
			}
		}
		return length;
	}
	
	
	/**
	 * 
	 * 
	 * 
	 * 
	 * @param motif
	 * 			name of the motif
	 * @param position
	 * 			the position within the motif. This is 0-based
	 * @param amino_acid
	 * 		 	the amino acid
	 * @return
	 * 			the score for this amino acid at this position in this motif.
	 */
	public int getScore(String motif, int position, char amino_acid) {
		//System.out.println(motif + "\t" + amino_acid + "\t" +  position + "\t" );
		if( (int) amino_acid - 65 <0) {
			return 0;
		}else {
			return this.pwm.get(motif)[amino_acid-65][position];
		}
		
	}
	
	
	/**
	 * 
	 * 
	 * @param motif
	 * @param i
	 * @return
	 */
	public double getCDF(String motif, int i) {
		
		return this.cdf.get(motif)[i];
	
	}
	
	private void loadMotifLengths(File inputFile) throws IOException{
		this.motifNames = new Vector<String>();
		motifLengths = new HashMap<String, Integer>();
		
		BufferedReader in = new BufferedReader(new FileReader(inputFile));

		
		
		
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\\s+")[0].split("@");
			String motif = split[0];
			int pos = Integer.parseInt(split[1]) +1; //the mot.txt file is 0-based
			
			if(!motifLengths.containsKey(motif)) {
				this.motifNames.add(motif);
				motifLengths.put(motif, pos);
				
				
			}else if ( pos > motifLengths.get(motif)) {
				motifLengths.put(motif, pos);
			}
				
					
					
			
				
				
			
			
			
		}

		in.close();
		/*
		for( int i = 1; i<= 20; i++) {
			this.motifNames.add("motif_"+i);
		}
		*/
	}
	
	private void loadPWM(File inputFile)throws IOException {
		
		pwm = new HashMap<String, int[][]>();
		
		
		BufferedReader in = new BufferedReader(new FileReader(inputFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			String[] split = inputline.split("\\s+");
			//motif_18@32@S 19
			
			String[] split2 = split[0].split("@");
			
			String motif = split2[0];
			int position = Integer.parseInt(split2[1]);
			char aa = split2[2].toCharArray()[0];
			int score = Integer.parseInt(split[1]);
			
			if(!pwm.containsKey(motif)) {
				pwm.put(motif, new int[26][motifLengths.get(motif)]);
			}
			
			pwm.get(motif)[aa-65][position] = score; 
			
		}

		in.close();
		
		
		
	}
	
	
	private void loadCDF(File inputFile)throws IOException{
		HashMap<String, Integer> maxLengths = new HashMap<String, Integer>();
		
		BufferedReader in = new BufferedReader(new FileReader(inputFile));

		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			//motif_1@11 0.9795331304536605
			String[] split1 = inputline.split("\\s+");
			String[] split2 = split1[0].split("@");
			String motif = split2[0];
			int pos = Integer.parseInt(split2[1]);
			//double d = Double.parseDouble(split1[1]);
		
			if(!maxLengths.containsKey(motif) || pos > maxLengths.get(motif)) {
				
					maxLengths.put(motif,pos);
				
			}
		}

		
		in.close();
		
		
		in = new BufferedReader(new FileReader(inputFile));

		this.cdf = new HashMap<String, double[]>();
		
		
		for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
			
			String[] split1 = inputline.split("\\s+");
			String[] split2 = split1[0].split("@");
			String motif = split2[0];
			int pos = Integer.parseInt(split2[1]);
			double d = Double.parseDouble(split1[1]);
		
			if(!cdf.containsKey(motif)) {
				cdf.put(motif, new double[maxLengths.get(motif)+1]);
			}
			cdf.get(motif)[pos] = d;
		}

		in.close();
		
		
		
	}
	
	public Vector<String> getMotifNames(){
		return this.motifNames;
	}
	
}
