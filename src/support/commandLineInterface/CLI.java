package support.commandLineInterface;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

public class CLI {

	Hashtable<String,Vector<String>> options;
	
	
	
	public CLI(){
		options = new Hashtable<String,Vector<String>>();
	}
	
	
	public boolean hasOption(String option){
		if(options.containsKey(option)){
			return true;
		}
		return false;
	}
	
	
	
	public String getArg(String option){
		return options.get(option).get(0);
	}
	
	public Vector<String> getArgs(String option){
		return options.get(option);
	}
	
	public void parseOptions(String[] args){
				
		boolean inOption = false;
		String optionValue = "";
		Vector<String> arguments = new Vector<String>();
		
		for( int i = 0; i< args.length; i++){
			if(args[i].startsWith("-")){
				if(inOption){
					options.put(optionValue, arguments);
				}
				inOption = true;
				optionValue = args[i];
				while(optionValue.startsWith("-")){
					optionValue = optionValue.substring(1);
				}
				
				arguments = new Vector<String>();
			}else if(!args[i].trim().equalsIgnoreCase("")){
				arguments.add(args[i]);
			}
		}
		if(inOption){
			options.put(optionValue, arguments);
		}
	}
	
	
	public String report(){
		String s = "";
		
		for(Enumeration<String> myenum1 = this.options.keys(); myenum1.hasMoreElements();){
			String key = myenum1.nextElement();
			Vector<String> v = options.get(key);
			s = s + "-"+key;
			for(Enumeration<String> myenum2 = v.elements(); myenum2.hasMoreElements();){
				s = s + " " + myenum2.nextElement();
			}
			s = s + "\n";
		}
		
		
		return s;
	}
	
	
}
