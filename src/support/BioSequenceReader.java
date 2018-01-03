package support;

import java.io.IOException;

public  interface BioSequenceReader {

	
	
	
	public  BioSequence readEntry() throws IOException;
	public void close()throws IOException;
	
	
	
}
