package support.commandLineInterface;

public class CLIParseException extends Exception {
	static final long serialVersionUID = 42L;
	
	public CLIParseException(){
		super();
		
	}

    public CLIParseException (String message) {
        super (message);
        
    }

    public CLIParseException (Throwable cause) {
        super (cause);
    }

    public CLIParseException (String message, Throwable cause) {
        super (message, cause);
    }
}
