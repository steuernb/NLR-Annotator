package support;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class BlastIPKReader implements BlastReader {

	
	

	File inputFile;
	BufferedReader in;
	BlastHSP lastHSP;
	
	
	public BlastIPKReader(File inputFile)throws IOException{
		this.initialize(inputFile);
		
		
		String inputline = in.readLine();
		while(inputline.startsWith("#") || inputline.startsWith("query_name\tquery_description\tquery_length")){
			inputline = in.readLine();
		}
		//System.out.println(inputline);
		lastHSP =BlastHSP.setHSPFromIPKParser(inputline);
	}
	
	
	/**
	 * 
	 * \check whether the file is gzipped or not and initialize the buffered reader
	 * 
	 * @param file
	 * @throws IOException
	 */
	private void initialize(File file)throws IOException{
		
		this.inputFile = file;
		
		//check if the fastq is gzipped
		FileInputStream fis = new FileInputStream(file);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		
		
			
		if(gzip){
			this.in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
		}else{
			this.in = new BufferedReader((new FileReader(file)));
		}
	}
	
	
	
	
	
	
	
	
	public BlastIteration readIteration()throws IOException{
		
		
		if(lastHSP== null){
			return null;
		}
		
		int hitRankCounter =0;
		
		BlastIteration blastIteration = new BlastIteration(lastHSP.getQueryName());
		blastIteration.setQueryLength(lastHSP.getQueryLength());
		
		BlastHit blastHit = new BlastHit(lastHSP.getQueryName(), lastHSP.getHitName());
		hitRankCounter++;
		blastHit.setHitRank(hitRankCounter);
		
		
		blastHit.addHSP(lastHSP);
		
		
		//this bit is probably not necessary but just to make sure comments are skipped and headers from concatenated files are skipped as well.
		String inputline = in.readLine();
		while( inputline != null && (inputline.startsWith("#") || inputline.startsWith("query_name\tquery_description\tquery_length"))){
			inputline = in.readLine();
		}
		if(inputline == null){
			lastHSP = null;
			blastIteration.addHit(blastHit);
			return blastIteration;
		}
		BlastHSP hsp = BlastHSP.setHSPFromIPKParser(inputline);
		lastHSP = hsp;
		
		while( hsp.getQueryName().equalsIgnoreCase(blastIteration.getQueryID())){
			if( hsp.getHitName().equalsIgnoreCase(blastHit.getHitID())){
				blastHit.addHSP(hsp);
			}else{
				blastIteration.addHit(blastHit);
				blastHit = new BlastHit(hsp.getQueryName(), hsp.getHitName());
				hitRankCounter++;
				blastHit.setHitRank(hitRankCounter);
				blastHit.addHSP(BlastHSP.setHSPFromIPKParser(inputline));
			}
			
			
			inputline = in.readLine();
			while( inputline != null && (inputline.startsWith("#") || inputline.startsWith("query_name\tquery_description\tquery_length"))){
				inputline = in.readLine();
			}
			if(inputline == null){
				lastHSP = null;
				blastIteration.addHit(blastHit);
				return blastIteration;
			}
			hsp = BlastHSP.setHSPFromIPKParser(inputline);
			lastHSP = hsp;
		}
		blastIteration.addHit(blastHit);
		
		
		
		
		return blastIteration;
	}
	
	
	/**
	 * @deprecated
	 * @return
	 * @throws IOException
	 */
	public BlastIteration readIteration_old()throws IOException{
		
				
		if(lastHSP == null){
			return null;
		}
				
		
		
		BlastIteration iteration = new BlastIteration(lastHSP.getQueryName());
		BlastHit b = new BlastHit(lastHSP.getQueryName(), lastHSP.getHitName());
		b.addHSP(lastHSP);
		
		String inputline = in.readLine();
		while(inputline != null &&( inputline.startsWith("#") || inputline.startsWith("query_name\tquery_description\tquery_length"))){
			inputline = in.readLine();
		}
		
		
		if(inputline != null){
			lastHSP = BlastHSP.setHSPFromIPKParser(inputline);
		}else{
			lastHSP= null;
			iteration.addHit(b);
			return iteration;
		}
		
		while(lastHSP.getQueryName().equalsIgnoreCase(iteration.getQueryID())){
			
			if(b.getHitID().equalsIgnoreCase(lastHSP.getQueryName())){
				b.addHSP(lastHSP);
			}else{
				iteration.addHit(b);
				b=new BlastHit(lastHSP.getQueryName(), lastHSP.getHitName());
				b.addHSP(lastHSP);
			}
		
			inputline = in.readLine();
		
			while(inputline != null &&( inputline.startsWith("#") || inputline.startsWith("query_name\tquery_description\tquery_length"))){
				inputline = in.readLine();
			}
			if(inputline != null){
				lastHSP = BlastHSP.setHSPFromIPKParser(inputline);
			}else{
				lastHSP = null;
				iteration.addHit(b);
				return iteration;
			}
		}
		
		return iteration;
	}
	
	
	

	
	
	
	
	
	
	public void close()throws IOException{
		in.close();
		}
	
	
	
	
	
	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		

	}

}
