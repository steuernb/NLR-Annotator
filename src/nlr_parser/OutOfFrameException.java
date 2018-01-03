package nlr_parser;

public class OutOfFrameException extends Exception {
	static final long serialVersionUID = 42L;
	
	public OutOfFrameException(){
		super();
		
	}

    public OutOfFrameException (String message) {
        super (message);
        
    }

    public OutOfFrameException (Throwable cause) {
        super (cause);
    }

    public OutOfFrameException (String message, Throwable cause) {
        super (message, cause);
    }
	
	
	
}
