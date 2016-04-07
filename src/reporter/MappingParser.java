package reporter;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StringReader;

/**
 * This class reads control files and multiple .map files, store # of each type of
 * nucleotide substitution by summarizing multiple sites (from .map files). 
 * There are N sites in an alignment. 
 * There are C trees / steps in a .map files (single site).
 * @author Hui-Jie
 *
 */
public class MappingParser {
	
	/** Start index */
	private int start;
	/** End index*/
	private int end;
	/** Number of iterations / trees in a file*/
	private int C;
	/** Filename prefix */
	private String name;
    /** Number of branches */
	private int branchNum;
    /** Store proportion of states*/
	private double[][][] propStates;
	/** Store time of states */
	private double[][][] timeStates;
	/** Store branch lengths */
	private double[][] br;
	/** Store number of changes */
	private int[][][] numberOfChanges;
	/** Store the tree structure to produce newick format for Multidivtime output file */
	private Tree tree[];
	/** Store the filename of the outgroup */
	private String outgroupFileName;
	/** The number of types of changes. 12 for SingleSiteParser, 4 for TripletParser */
	private int numTypeChanges;
	/** The number of types of states. 4 for SingleSiteParser, 3 for TripletParser */
	private int numTypeStates;
	/** Store the root state count for each mapping */
	private int[][] rootStateCount;
	
	/**
	 * Constructor
	 * @param start
	 * @param end
	 * @param C
	 * @param name
	 * @param outgroupFileName
	 */
	public MappingParser(int start, int end, int C, String name, String outgroupFileName) {
		this.start = start;
		this.end = end;
		this.C = C;
		this.name = name;
		this.outgroupFileName = outgroupFileName;
		this.numTypeChanges = 18;
		this.numTypeStates = 6;
		this.tree = new Tree[3];
		this.rootStateCount = new int[4][C]; //G,C,T,A
		
		System.out.println(new File("").getPath().toString());
		
		String inputMap = name+"_"+start+".map";

		try {
			// returns the ClassLoader object associated with this Class
	        //ClassLoader cLoader = this.getClass().getClassLoader();
	        // input stream
			//InputStream inStream = cLoader.getResourceAsStream("./parse_phylobayes/"+inputMap);
			//InputStream inStream = MappingParser.class.getClassLoader().getResourceAsStream(inputMap);
			InputStream inStream = this.getClass().getResourceAsStream(new File("../" + inputMap).getPath().toString());
			BufferedReader r = new BufferedReader(new InputStreamReader(inStream));
			//need to figure out the number of branches first so that 
			//i can declare the size of the array to store info
			String line = r.readLine();
			StringReader str = new StringReader(line);
            TreeParser tp = new TreeParser(str, outgroupFileName);
            //store the tree structure here
            //note that this tree shares the same node/branch numbering and ancestral with all other trees
            this.tree[0] = tp.tokenize();
            this.branchNum = tree[0].getNumBranches();
            r.close();
            inStream.close();
		} catch (IOException e) {
            e.printStackTrace();
        } 
		
		try {
			setUp();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Get number of branches
     * @return branchNum
	 */
	public int getBranchNum() {
		return branchNum;
	}
	
	/**
	 * Get number of changes and time in states for each mapping
	 * 
	 * @throws IOException 
	 *             
	 */
	public void setUp() throws IOException {
		numberOfChanges = new int[branchNum][numTypeChanges][C];
		//changesInGroups = new int[numGroup][branchNum][C]; 
		timeStates = new double[branchNum][numTypeStates][C];
		br = new double[branchNum][C];
		propStates = new double[branchNum][numTypeStates][C];
		
		String inputMap = null;
		//File f = null;
		//move window from site 1 (index 0) to site N-2 (index N-3)
		for (int i = start; i < (end-1); i++) {		
			System.out.println("i="+i);
			String[] line = new String[3];
			StringReader[] str = new StringReader[3];
			BufferedReader[] r = new BufferedReader[3];
			InputStream[] inStream = new InputStream[3];
			for (int j = i; j < (i+3); j++) { //construct 3 bufferedreader
				inputMap = name+"_"+j+".map";
				//f = new File(inputMap);
				//r[j-i] = new BufferedReader(new FileReader(f));
				//inStream[j-i] = MappingParser.class.getClassLoader().getResourceAsStream(inputMap);
				inStream[j-i] = this.getClass().getResourceAsStream(new File("../" + inputMap).getPath().toString());					
				r[j-i] = new BufferedReader(new InputStreamReader(inStream[j-i]));
			}
				
			for (int j = 0; j < C; j++) { //read C trees for each file 
				try {
					for (int k = 0; k < 3; k++) { //repeat for each file (total 3 files at a time)
						//read first tree
						line[k] = r[k].readLine();
						str[k] = new StringReader(line[k]);
		                TreeParser tp = new TreeParser(str[k], outgroupFileName);
		                tree[k] = tp.tokenize();
		                //read second tree, discard this tree
                        //REMOVE THIS LINE IF PHYLOBAYES HAS BEEN CHANGED TO INCLUDE ONLY ONE MAPPING PER MCMC ITERATION
		                line[k] = r[k].readLine();
		                //read "" and discard it.
		                line[k] = r[k].readLine();
					}
					SiteParser parse = new TripletParser(tree);
						
					for (int l = 0; l < branchNum; l++) {
						//br = new double[branchNum][C];
		               	br[l][j] = parse.getBranchLengths()[l]; 
		               	for (int k = 0; k < numTypeStates; k++) {
		               		//numberOfChanges = new int[4][trees[0].getNumBranches()];
		               		numberOfChanges[l][k][j] += parse.getNumberOfChanges()[l][k];
		               		//timeOfStates = new double[3][trees[0].getNumBranches()];
		               		timeStates[l][k][j] += parse.getTimeOfStates()[l][k];
		               		//propStates = new double[3][trees[0].getNumBranches()];
		               		propStates[l][k][j] += parse.getPropStates()[l][k];
		                }//end k
		                for (int k = numTypeStates; k < numTypeChanges; k++) {
		                	numberOfChanges[l][k][j]+=parse.getNumberOfChanges()[l][k];
		                }//end k

		             }//end l	
					
					//increase root state count
					rootStateCount[parse.getRootState()][j]++;
						
				} catch (IOException e) {
					  // TODO Auto-generated catch block
					  e.printStackTrace();
				} 
					
			}//end j
				
			//close bufferedreader
			for (int k = 0; k < 3; k++) {
				r[k].close();
				inStream[k].close();
			}
				
		} //end i
			
	}
	
	
    /**
     * Return the number of changes for each substitution type (scenario 1) on branches for each MCMC iteration.
     * @return numberOfChanges
     */
	public int[][][] getNumberOfChanges() {
		return numberOfChanges;
	}
	
	
    /**
	 * Get proportion of time the second position is a CpG / non-CpG site for each branch for each MCMC iteration.
	 * @return propStates
	 */
	public double[][][] getPropStates() {
		return propStates;
	}
	
    /**
	 * Get branch lengths for each branch for each MCMC iteration.
	 * @return br
	 */
	public double[][] getBranchLengths() {
		return br;
	}
	
	/**
	 * Get branch lengths for the first iteration. 
	 * Note, in the case of fixed branch lengths, the branch lengths should be the same for different iterations.
	 * @return branch lengths
	 */
	public double[] getBranchLengthFirst() {
		double[] branchLength = new double[branchNum];
		for (int j = 0; j < branchNum; j++) {
			branchLength[j] = br[j][0];
		}
		return branchLength;
	}


    /**
     * Return the total time in each states on each branch for each MCMC iteration
     * @return timeStates
     */
	public double[][][] getTimeStates() {
		return timeStates;
	}
	
	/**
	 * Return the tree 
	 * @return tree
	 */
	public Tree getTree(){
		return tree[1];
	}	
	
	public int[][] getRootStateCount() {
		return rootStateCount;
	}
	
	
}

