package support;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class BlastXMLReader implements BlastReader {

	File blastFile;
	BufferedReader in;
	String algorithm;
	String database_name;
	String inputline;
	
	
	
	
	
	public BlastXMLReader(File inputFile)throws IOException{
		this.blastFile = inputFile;
		
		
		
		FileInputStream fis = new FileInputStream(inputFile);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		
		
		
		
		if(gzip){
			this.in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFile))));
		}else{
			this.in = new BufferedReader(new FileReader(inputFile));
		}
		
		
		
		
		this.inputline = in.readLine();
		
		/*
		 * parse out the used algorithm
		 */
		
		String def = "";
		while(inputline!= null && !inputline.contains("<BlastOutput_iterations>")){
			def = def + inputline.trim();
			inputline = in.readLine();
		}
		
		def = def + inputline.split("<BlastOutput_iterations>",-1)[0];
		
		this.algorithm     = def.split("<BlastOutput_program>")[1].split("</BlastOutput_program>")[0];
		this.database_name = def.split("<BlastOutput_db>")[1].split("</BlastOutput_db>")[0];
		
		
	}
	
	
	public BlastXMLReader(InputStreamReader inputStreamReader)throws IOException{
		this.blastFile = null;
		this.in = new BufferedReader((inputStreamReader));
		
		this.inputline = in.readLine();
		
		/*
		 * parse out the used algorithm
		 */
		
		String def = "";
		while(inputline!= null && !inputline.contains("<BlastOutput_iterations>")){
			def = def + inputline.trim();
			inputline = in.readLine();
		}
		
		def = def + inputline.split("<BlastOutput_iterations>",-1)[0];
		
		this.algorithm     = def.split("<BlastOutput_program>")[1].split("</BlastOutput_program>")[0];
		this.database_name = def.split("<BlastOutput_db>")[1].split("</BlastOutput_db>")[0];
		
	}
	
	public File getFile(){
		return this.blastFile;
	}
	
	public BlastIteration readIteration()throws IOException{
		return readIteration(false);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public BlastIteration readIteration(boolean returnEmpty)throws IOException{
		
		
		
		
		/*
		 * Initialize the blast iteration
		 */
		
		String iteration = getIterationString();
		
		if( !returnEmpty){
			while(iteration != null && iteration.contains("<Iteration_message>No hits found</Iteration_message>") ){
				iteration = getIterationString();
				if( iteration == null){
					return null;
				}
			}
			
			
		}
		if(iteration == null){
			return null;
		}
		
		
		String query = iteration.split("<Iteration_query-def>")[1].split("</Iteration_query-def>")[0];
		
		String query_id = query.split("\\s")[0];
		
		String query_description = "";
		try{ query_description = query.substring(query_id.length()+1);}catch (StringIndexOutOfBoundsException e ){}  //catch if there is no description
		int query_length = Integer.parseInt(iteration.split("<Iteration_query-len>")[1].split("</Iteration_query-len>")[0]);
		
		
		BlastIteration blastIteration = new BlastIteration(query_id);
		blastIteration.setAlgorithm(this.algorithm);
		blastIteration.setQueryLength(query_length);
		
		
		String[] hits = iteration.split("<Hit>");
		/*
		 * read the hits and add them to the iteration
		 */
		for( int i = 1; i< hits.length; i++){
			
			String hit = hits[i].split("<Hit_def>")[1].split("</Hit_def>")[0];
			String hit_id = hit.split("\\s")[0];
			String hit_description = "";
			try{ hit_description = hit.substring(hit_id.length() +1);}catch (StringIndexOutOfBoundsException e ){}  //catch if there is no description
			int hit_length = Integer.parseInt(hits[i].split("<Hit_len>")[1].split("</Hit_len>")[0]);
			int hit_rank = Integer.parseInt(hits[i].split("<Hit_num>")[1].split("</Hit_num>")[0]);
			
			BlastHit blastHit = new BlastHit(query_id, hit_id);
			blastHit.setHitRank(hit_rank);
			
			/*
			 * read the HSPs for Hit
			 */
			String[] hsps = hits[i].split("<Hsp>");
			for( int j = 1; j<hsps.length; j++){
				int hsp_rank = Integer.parseInt(hsps[j].split("<Hsp_num>")[1].split("</Hsp_num>")[0]);
				double bitScore = Double.parseDouble(hsps[j].split("<Hsp_bit-score>")[1].split("</Hsp_bit-score>")[0]);
				double rawScore = Double.parseDouble(hsps[j].split("<Hsp_score>")[1].split("</Hsp_score>")[0]);
				String evalue = hsps[j].split("<Hsp_evalue>")[1].split("</Hsp_evalue>")[0];
				int query_start = Integer.parseInt(hsps[j].split("<Hsp_query-from>")[1].split("</Hsp_query-from>")[0]);
				int query_end = Integer.parseInt(hsps[j].split("<Hsp_query-to>")[1].split("</Hsp_query-to>")[0]);
				int hit_start = Integer.parseInt(hsps[j].split("<Hsp_hit-from>")[1].split("</Hsp_hit-from>")[0]);
				int hit_end = Integer.parseInt(hsps[j].split("<Hsp_hit-to>")[1].split("</Hsp_hit-to>")[0]);
				
				int query_frame = Integer.parseInt(hsps[j].split("<Hsp_query-frame>")[1].split("</Hsp_query-frame>")[0]);
				int hit_frame = Integer.parseInt(hsps[j].split("<Hsp_hit-frame>")[1].split("</Hsp_hit-frame>")[0]);
				int numIdentical = Integer.parseInt(hsps[j].split("<Hsp_identity>")[1].split("</Hsp_identity>")[0]);
				int numPositive = Integer.parseInt(hsps[j].split("<Hsp_positive>")[1].split("</Hsp_positive>")[0]);
				int numGaps = Integer.parseInt(hsps[j].split("<Hsp_gaps>")[1].split("</Hsp_gaps>")[0]);
				int alignment_length = Integer.parseInt(hsps[j].split("<Hsp_align-len>")[1].split("</Hsp_align-len")[0]);
				String query_string = hsps[j].split("<Hsp_qseq>")[1].split("</Hsp_qseq>")[0];
				String hit_string = hsps[j].split("<Hsp_hseq>")[1].split("</Hsp_hseq>")[0];
				String alignment_string = hsps[j].split("<Hsp_midline>")[1].split("</Hsp_midline>")[0];
				
				BlastHSP blastHSP = new BlastHSP( query_id,  hit_id,  query_description,  hit_description, 
				         alignment_length,  query_length,  hit_length, 
				         numIdentical,  -1,  numGaps,  -1,  numPositive,
				         query_start,  query_end,  hit_start,  hit_end,
				         rawScore,  bitScore,  evalue,  query_string,  hit_string,
				         query_frame,  hit_frame,  0,  0);
				
				blastHSP.setHspRank(hsp_rank);
				blastHSP.setNumHSPs(hsps.length-1);
				blastHSP.setNumHits(hits.length-1);
				blastHSP.setAlignmentString(alignment_string);
				blastHSP.setAlgorithm(this.algorithm);
				blastHSP.setDatabaseName(this.database_name);
				blastHit.addHSP(blastHSP);
				
			}
			
			
			
			
			
			blastIteration.addHit(blastHit);
		}
			
			
				
		
		
		
		
		
		
		
		
		return blastIteration;
	}
	
	
	
	private String getIterationString()throws IOException{
		StringBuilder builder = new StringBuilder();
		
		/*
		 * move to start of iteration
		 */
		while(inputline !=null && !inputline.contains("<Iteration>")){
			inputline = in.readLine();
		}
		if(inputline == null){
			return null;
		}
		
		/*
		 * read the entire iteration
		 */
		builder.append(inputline.split("<Iteration>", -1)[1]);
		inputline = in.readLine();
		
		while (inputline != null && !inputline.contains("</Iteration>")){
			builder.append(inputline.trim());
			inputline = in.readLine();
		}
		
		
		
		builder.append(inputline.split("</Iteration>", -1)[0].trim());		
		
		
		return builder.toString();
	}
	
	
	
	
	
	public void close()throws IOException{
		in.close();
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			

			
			BlastXMLReader reader = new BlastXMLReader(new File("/Users/steuernb/Documents/projects/test/test.tblastn.xml"));
			BlastIteration it = reader.readIteration(true);
		
			while(it != null){
				System.out.println(it.getQueryID() + " " + it.getHits().size() + " " + it.getQueryDescription());
				it = reader.readIteration(true);
			
				
			}
			
			reader.close();
			
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		

	}

}
