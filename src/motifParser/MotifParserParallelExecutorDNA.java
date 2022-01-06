package motifParser;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import annotator.SequenceChopper;
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
 *  @author Burkhard Steuernagel
 * 
 */
public class MotifParserParallelExecutorDNA {

	
	int numThreads;
	int numSequencesPerCall;
	
	int fragmentLength;
	int overlap;
	
	MotifDefinition motifDefinition;
	SignatureDefinition signatureDefinition;
	File inputFastaFile;
	
	
	public MotifParserParallelExecutorDNA(MotifDefinition motifDefinition,
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
		
	}
	
	

/*

	public Vector< MotifList> execute( ) throws IOException, ExecutionException, InterruptedException{
		
		
		Vector<MotifList> resultMotifList = new Vector<MotifList>();
		
		ExecutorService executor = Executors.newFixedThreadPool(this.numThreads);
		List<Future<Vector<MotifList>>> list = new ArrayList<Future<Vector<MotifList>>>();
		
		Vector<BioSequence> seqs = new Vector<BioSequence>();
		
		SequenceChopper chopper = new SequenceChopper(inputFastaFile, this.fragmentLength, this.overlap);
		
		
		for (BioSequence seq = chopper.nextSequence()  ; seq != null; seq = chopper.nextSequence()) {
			
			BioSequence[] translatedSeqs = seq.translate2Protein();
			
			for( int i = 0; i< translatedSeqs.length; i++) {
				seqs.add(translatedSeqs[i]);
			}
			
			
			
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
		chopper.close();
		
		
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
					//System.out.println(motifList.getMotifListName() + "\t" + motifList.getMotifString() );
				}
				
			}
				
						
			
		}
		
				
		executor.shutdown();
		
		
		
		
		return resultMotifList;
		
	}
	
	
	*/
	
	
}
