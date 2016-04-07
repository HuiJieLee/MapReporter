/**
 * 
 */
package reporter;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.*;

/**
 * Store substitution types subset as ArrayList of ArrayList.
 * @author Hui-Jie Lee
 *
 */
public class Main {
	private MappingParser parse;
	
	
	/**
	 * Constructor
	 */
	public Main(MappingParser parse) {
		this.setParse(parse);
	}

	/** 
	 * Main class
	 * @param args[0]: start index
	 * @param args[1]: end index
 	 * @param args[2]: number of iterations
	 * @param args[3]: prefix of filename
	 * @param args[4]: outgroup filename, currently will produce incorrect tree topology if an ourgroup file is not given.
	 * @param args[5]: gtr parameter file
	 */
	public static void main(String args[]) {
		int start = Integer.parseInt(args[0]);
		int end = Integer.parseInt(args[1]);
		int C = Integer.parseInt(args[2]);
		String name = args[3];
		String outgroup = args[4];
		String gtr_param = args[5];		
		if (args.length == 6) {	
			MappingParser parse = new MappingParser(start, end, C, name, outgroup);	
			
			Main main = new Main(parse);
			// call GTR to calculate log P(M^{(c)}, X|mu_GTR) for each mapping.
			SufficientStatistics gtrSuff = new SufficientStatistics(true, parse.getBranchNum(), C, parse.getPropStates(), parse.getNumberOfChanges());
			Object[] param = main.gtrParameter(gtr_param);
			double[] pi = (double[]) param[0];
			double[] R = (double[]) param[1];		
			GTR gtr = new GTR(pi, R, gtrSuff, parse.getRootStateCount(), parse.getBranchLengthFirst());

			//output sufficient stat files and gtr weights
			printSufficientStatistics(parse, start, end, C, gtr);

		}  else {	//error
			System.out.println("Argument error.");
		}
			
		
	}
	
	/**
	 * This method prints sufficient statistics to files.
	 * It creates four output files. 
	 * @param parse
	 */
	private static void printSufficientStatistics(MappingParser parse, int start, int end, int C, GTR gtr) {
		int branchNum = parse.getBranchNum();
		double[][][] prop = parse.getPropStates();		// double[branchNum][numTypeStates][C]; numTypeStates = 6
		int[][][] count = parse.getNumberOfChanges();	//int[branchNum][numTypeChanges][C]; numTypeChanges = 18
		int[][] rootCount = parse.getRootStateCount();	//int[4][C]; //G,C,T,A
		double[] w_gtr = gtr.getLogL();
		
		File output;
		output = new File("PropState"+start+"_"+end+".txt");
		PrintStream print = null;
		try {
			print = new PrintStream(output);
			for (int j = 0; j < branchNum; j++) {
				for (int l = 0; l < 6; l++) {
					//print.println("Branch = " + j + " ; State = "+l);
					for (int c = 0; c < C; c++) {
						print.print(prop[j][l][c]);
						print.print(" ");
					}
					print.println();
				}
			}
		} catch (FileNotFoundException e) {
			System.out.println("Problem creating PropState file!");
		} finally {
  	        if (print != null) print.close();
  	    }
		
		output = new File("NumChange"+start+"_"+end+".txt");
		try {
			print = new PrintStream(output);
			for (int j = 0; j < branchNum; j++) {
				for (int k = 0; k < 18; k++) {
					//print.println("Branch = " + j + " ; Type = "+k);
					for (int c = 0; c < C; c++) {
						print.print(count[j][k][c]);
						print.print(" ");
					}
					print.println();
				}
			}
		} catch (FileNotFoundException e) {
			System.out.println("Problem creating NumChange file!");
		} finally {
  	        if (print != null) print.close();
  	    }
		
		output = new File("RootState"+start+"_"+end+".txt");
		try {
			print = new PrintStream(output);
			for (int l = 0; l < 4; l++) {
				//print.println("Root State = " + l);
				for (int c = 0; c < C; c++) {
					print.print(rootCount[l][c]);
					print.print(" ");
				}
				print.println();
			}
		} catch (FileNotFoundException e) {
			System.out.println("Problem creating RootState file!");
		} finally {
  	        if (print != null) print.close();
  	    }
		
		output = new File("GTRweight"+start+"_"+end+".txt");
		try {
			print = new PrintStream(output);
			for (int c = 0; c < C; c++) {
				print.print(w_gtr[c]);
				print.print(" ");
			}
		} catch (FileNotFoundException e) {
			System.out.println("Problem creating GTRweight file!");
		} finally {
  	        if (print != null) print.close();
  	    }
	}
	
	
	
	
	/**
	 * This function parse the file that stores substitution types grouping info.
	 * @param filename
	 * @return an arraylist of arraylist storing grouping info. each list contains a list of elements in that subset.
	 */
	public ArrayList<ArrayList<Integer>> grouping(String filename){
	    ArrayList<ArrayList<Integer>> outer = new ArrayList<ArrayList<Integer>>();
	    ArrayList<Integer> inner = new ArrayList<Integer>();   
	    int N_G;
		    
	    InputStream inStream = this.getClass().getResourceAsStream(new File("../" + filename).getPath().toString());
    	BufferedReader r = new BufferedReader(new InputStreamReader(inStream));
    	String line;
		try {
			line = r.readLine();
			//first line contains the number of subsets
             N_G = Integer.parseInt(line);             
             for (int i = 0; i < N_G; i++) {
             	line = r.readLine(); // each line contains numbers in that subset
             	String[] elements = line.split(" ");
             	for (int j = 0; j < elements.length; j++) {
             		int element = Integer.parseInt(elements[j]);
             		inner.add(element-1);	//index = element - 1, i.e. type = 0, ..., 8 instead of 1,..., 9
             	}
             	outer.add(inner);
             	inner = new ArrayList<Integer>();
             }
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
            try {
                inStream.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
	
	    return outer;
	}
	
	/**
	 * This file parse the input file that stores GTR parameter values
	 * First line is pi pi_G, pi_C, pi_T, pi_A
	 * Second line is R R_GC/R_CG, R_GT/R_TG, R_GA/R_AG, R_CT/R_TC, R_CA/R_AC, R_TA/R_AT
	 * @param filename
	 * @return
	 */
	public Object[] gtrParameter(String filename) {
		double[] pi = new double[4];
		double[] R = new double[6];
		
		InputStream inStream = this.getClass().getResourceAsStream(new File("../" + filename).getPath().toString());
    	BufferedReader r = new BufferedReader(new InputStreamReader(inStream));
    	String line;
		try {
			line = r.readLine();
			String[] elements = line.split(" ");
			if(elements.length != 4) {
				System.out.println("Incorrect format in GTR parameter values! (pi)");
			}
			for (int i = 0; i < 4; i++) {
				double element = Double.parseDouble(elements[i]);
				pi[i] = element;
             }
			line = r.readLine();
			elements = line.split(" ");
			if(elements.length != 6) {
				System.out.println("Incorrect format in GTR parameter values! (R)");
			}
			for (int i = 0; i < 6; i++){
				double element = Double.parseDouble(elements[i]);
				R[i] = element;
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} finally {
            try {
                inStream.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
		
		return new Object[]{pi, R};
	}

	public MappingParser getParse() {
		return parse;
	}

	public void setParse(MappingParser parse) {
		this.parse = parse;
	}

}
