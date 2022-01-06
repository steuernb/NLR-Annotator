package motifParser;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import annotator.SequenceChopper;
import core.Motif;
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
 *  @author Burkhard Steuernagel
 * 
 */
public class MotifParserParallelExecutor {

	int numThreads;
	int numSequencesPerCall;
	
	MotifDefinition motifDefinition;
	SignatureDefinition signatureDefinition;
	File inputFastaFile;
	
	File tempdir;
	
	
	int fragmentLength;
	int overlap;
	
	
	public MotifParserParallelExecutor(MotifDefinition motifDefinition,
									   SignatureDefinition signatureDefinition, 
									   File inputFastaFile, 
									   int numThreads, 
									   int numSequencesPerCall,
									   int fragmentLength,
									   int overlap
									   ) {
		
		this.numSequencesPerCall = numSequencesPerCall;
		this.numThreads = numThreads;
		this.motifDefinition = motifDefinition;
		this.signatureDefinition = signatureDefinition;
		this.inputFastaFile = inputFastaFile;
		this.fragmentLength = fragmentLength;
		this.overlap = overlap;
		
		this.tempdir = new File( "temp_mast1" );
		while(tempdir.exists()){
			tempdir = new File("temp_mast" + ((int) Math.round( Math.random() *100000000  )));
		}
		tempdir.mkdir();
		
		
		System.err.println("DEBUG: Temporary directory is "+tempdir.getAbsolutePath() );
		
		
	}
	
	
	public HashMap<String,Vector<Motif>> executeDNA()throws IOException, InterruptedException, ExecutionException{
		ExecutorService executor = Executors.newFixedThreadPool(this.numThreads);
		List<Future<File>> list = new ArrayList<Future<File>>();
		
		int fileCounta = 1;
		int seqCounta = 0;
		File tempSeq = new File(tempdir, "tempSeq"+fileCounta+".fasta");
		File returnSeq =  new File(tempdir, "tempSeq"+fileCounta+".tsv");
		BufferedWriter out = new BufferedWriter(new FileWriter( tempSeq ));
		
		SequenceChopper chopper = new SequenceChopper(inputFastaFile, fragmentLength, overlap);
		for(BioSequence seq = chopper.nextSequence(); seq != null; seq = chopper.nextSequence() )	{
			
			seqCounta++;
			out.write(seq.getFastaString(100));
						
			if(seqCounta>=this.numSequencesPerCall){
				seqCounta = 0;
				out.close();
				MotifParserThread mt = new MotifParserThread( motifDefinition, signatureDefinition, tempSeq,returnSeq );
				Callable<File> worker = mt;
				Future <File> submit = executor.submit(worker);
				list.add(submit);
				fileCounta++;
				tempSeq = new File(tempdir, "tempSeq"+fileCounta+".fasta");
				returnSeq = new File(tempdir, "tempSeq"+fileCounta+".tsv");
				out = new BufferedWriter(new FileWriter( tempSeq ));
			}
		}	
		
		chopper.close();
		out.close();
		
		//last call with remaining sequences
		MotifParserThread mt = new MotifParserThread( motifDefinition, signatureDefinition ,  tempSeq, returnSeq);
		Callable<File> worker = mt;
		Future <File> submit = executor.submit(worker);
		list.add(submit);
		
		
		HashMap<String,Vector<Motif>> h = new HashMap<String,Vector<Motif>>();
		
		for(Future<File> future : list){
			File file = future.get();
			
			BufferedReader in = new BufferedReader(new FileReader(file));

			for (String inputline = in.readLine(); inputline != null; inputline = in.readLine()) {
				Motif motif = new Motif(inputline);
				if(!h.containsKey(motif.getDnaSequenceID())) {
					h.put(motif.getDnaSequenceID(), new Vector<Motif>());
				}
				h.get(motif.getDnaSequenceID()).add(motif);
				
			}

			in.close();
			
				
			
		}
		
				
		executor.shutdown();
			
		return h;
	}

	
	public boolean close(){
		
		System.err.println("DEBUG: Temporary directory "+tempdir.getAbsolutePath() +" is now being deleted.");
		
		return recoursiveDelete(tempdir);
	}
	

	private boolean recoursiveDelete(File dir){
		File[] files = dir.listFiles();
		for( int i = 0; i< files.length; i++){
			File file = files[i];
			if(file.isDirectory()){
				recoursiveDelete(file);
			}else{
				file.delete();
			}
		}
		
		return dir.delete();
	}
	
	/*
	
	public Vector< MotifList> execute( ) throws IOException, ExecutionException, InterruptedException{
		
		
		Vector<MotifList> resultMotifList = new Vector<MotifList>();
		
		ExecutorService executor = Executors.newFixedThreadPool(this.numThreads);
		List<Future<Vector<MotifList>>> list = new ArrayList<Future<Vector<MotifList>>>();
		
		Vector<BioSequence> seqs = new Vector<BioSequence>();
		
		FastaReader fastaReader = new FastaReader(inputFastaFile);
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			
			
			seqs.add(seq);
			
			if(seqs.size() >= numSequencesPerCall) {
				Vector<BioSequence> copy_of_seqs = new Vector<BioSequence>();
				copy_of_seqs.addAll(seqs);
				
				//MotifParserThread mt = new MotifParserThread(motifParser, signatureDefinition, copy_of_seqs, minSequenceLength);
				Callable<Vector<MotifList>> worker = new MotifParserThread( motifDefinition, signatureDefinition, copy_of_seqs);
				Future <Vector<MotifList>> submit = executor.submit(worker);
				list.add(submit);
				
				seqs.clear();
				
			}
			
			
			
			
		}
		fastaReader.close();
		
		
		if( seqs.size() >0) {
			Vector<BioSequence> copy_of_seqs = new Vector<BioSequence>();
			copy_of_seqs.addAll(seqs);
			
			Callable<Vector<MotifList>> worker = new MotifParserThread( motifDefinition, signatureDefinition, copy_of_seqs);
			Future <Vector<MotifList>> submit = executor.submit(worker);
			list.add(submit);
			seqs.clear();
		}
		
		
		for(Future<Vector<MotifList>> future : list){
			
			
			
			
			Vector<MotifList> mymotifs = future.get();
			
			
			
			
			//System.out.println("call: " +  mymotifs.size() + " motif lists");
			for( Iterator<MotifList> iterator = mymotifs.iterator(); iterator.hasNext();) {
				MotifList motifList = iterator.next();
				if( motifList.has_NLR_signature(signatureDefinition)) {
					resultMotifList.add(motifList);
					System.out.println(motifList.getMotifListName() + "\t" + motifList.getMotifString() );
				}
				
			}
				
						
			
		}
		
				
		executor.shutdown();
		
		
		
		
		return resultMotifList;
		
	}
	
	*/
	
	
}
