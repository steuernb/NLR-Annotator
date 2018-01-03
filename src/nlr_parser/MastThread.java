package nlr_parser;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.Callable;



public class MastThread implements Callable<File>{
	
	File sequenceFile;
	File tempDir ;
	
	File mastExe;
	File memeXML;

	String splitPattern;
	double pvalue_threshold;
	
	public MastThread( File sequenceFile, File tempDir , File mastExe, File memeXML){
		this.mastExe = mastExe;
		this.memeXML = memeXML;
		
		
		this.sequenceFile = sequenceFile;
		this.tempDir = tempDir;
		
		
	}
	
	
	public File call()throws InterruptedException, IOException{
		
		//System.out.println("Running Mast on " + sequenceFile.getName());
		
		int dirCounta = 1;
		File outputDir = new File(this.tempDir, "output_"+sequenceFile.getName() + "_" + dirCounta);
		while(outputDir.exists()){
			dirCounta++;
			outputDir = new File(this.tempDir, "output_"+sequenceFile.getName() + "_" + dirCounta);
		}
		
		
		
		Process process = new ProcessBuilder(this.mastExe.getAbsolutePath(), this.memeXML.getAbsolutePath(), this.sequenceFile.getAbsolutePath(), "-o", outputDir.getAbsolutePath() , "-ev", "1000").start();
		//Process process = new ProcessBuilder(this.mastExe.getAbsolutePath(), this.memeXML.getAbsolutePath(), this.sequenceFile.getAbsolutePath(), "-o", outputDir.getAbsolutePath() ).start();
		process.waitFor();
		
		
		
		File outputFile = new File(outputDir, "mast.xml");
		while(!outputFile.exists()){
			Thread.sleep(100);
		}
		
		
		
		
		
		
		return outputFile;
			
	
		
	}
	
	
	
	
	
	
	
}
