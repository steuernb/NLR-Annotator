package support.commandLineInterface;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;


/**
 *   CLI -  a light weight cmd line interface
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
