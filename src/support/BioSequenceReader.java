package support;

import java.io.IOException;


/**
 * 
 * 
 * @author Burkhard Steuernagel
 *
 */
public  interface BioSequenceReader {

	
	
	
	public  BioSequence readEntry() throws IOException;
	public void close()throws IOException;
	
	
	
}
