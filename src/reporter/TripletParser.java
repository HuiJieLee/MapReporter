package reporter;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * Parse three consecutive sites /triplet at a time.
 * There are 6 states: non-CpG G site, non-CpG C site, non-CpG T site, non-CpG A site, CpG G site, CpG C site
 * 18 kinds of changes: 
 * non-CpG site: G -> C (Index 0), C -> G (Index 1), G -> T (Index 2), C -> A (Index 3), T -> A (Index 4), A -> T (Index 5)
 *   			 T -> G (Index 6), A -> C (Index 7), G -> A (Index 8), C -> T (Index 9), A -> G (Index 10), T -> C (Index 11)
 * CpG site:     G -> C (Index 12), C -> G (Index 13), G -> T (Index 14), C -> A (Index 15), G -> A (Index 16), C -> T (Index 17)
 * 
 * @author Hui-Jie Lee
 *
 */
public class TripletParser implements SiteParser {
	
    /** Store three trees for three sites **/
	private Tree[] trees;
    
    /** Store number of changes  
     *  dim = (# of branches) x 18
     **/
	private int[][] numberOfChanges;
    
    /** Store the total time spent in each state on branches
     *  dim = (# of branches) x 6
     **/
	private double[][] timeOfStates;
    
    /** Store the branch lengths for each branch
     *  dim = # of branches
     **/
	private double[] branchLengths;
    
    /** Store the proportion of time in each state on branches
     *  dim = (# of branches) x 6
     **/
	private double[][] propStates;
    
	/** Store path state as an array of Arraylist of String, each String contains 3 characters (triplet). 
	 *  There are branchNum of such Arraylist for each branch.
     **/
	private ArrayList<String>[] pathState;
    
	/** Store path time as an array of arraylist of double.
	 * There are branchNum of such arraylist for each branch 
     **/
	private ArrayList<Double>[] pathTime;
	
	/** Get state G/C/T/A at the root */
	private int rootState;
    
	/**
	 * Constructor
	 * @param trees: 3 trees for 3 sites
	 */
	public TripletParser(Tree[] trees) {
		this.trees = new Tree[3];
		for (int i = 0; i < 3; i++) {
			this.trees[i] = trees[i];
		}
		
		numberOfChanges = new int[trees[0].getNumBranches()][18];
		timeOfStates = new double[trees[0].getNumBranches()][6];
		branchLengths = new double[trees[0].getNumBranches()];
		propStates = new double[trees[0].getNumBranches()][6];
	
		pathState = (ArrayList<String>[])new ArrayList[trees[0].getNumBranches()];
		pathTime = (ArrayList<Double>[])new ArrayList[trees[0].getNumBranches()]; 
		setRootState();

		for(int i = 0; i < trees[0].getNumBranches(); i++) {
			pathState[i] = new ArrayList<String>();
			pathTime[i] = new ArrayList<Double>();
		}
		
		constructWholePath();
	
		setBranchLengths(computeBranchLength());
		setNumberOfChanges(computeNumberOfChanges());
		setTimeOfStates(computeTimeOfStates());
		setPropStates(computePropState());
	}
    
    
	/** Return the number of branches in the tree
	 * @return number of branches
	 */
	public int numberOfBranches() {
		return branchLengths.length;
	}

	
	/**
	 * Wrapper class for constructPath(). Construct path for all branches.
	 */
	public void constructWholePath() {
		for (int i = 0; i < numberOfBranches(); i++) {
			constructPath(i);
		}
	}
	
	/**
	 * Construct the path state and path time by combining three sites/trees for a given branch.
	 * path state: triplet starting from the parent node to current node
	 * path time: time interval between substitution events.
	 *            The last time will be branch length of site 2 minus the time at the last substitution event.
	 * @param branchIndex index of the given branch (the index of the node that ends the branch)
	 */
	public void constructPath(int branchIndex) {
		ArrayList<String>[] state = (ArrayList<String>[])new ArrayList[3];
		ArrayList<Double>[] time = (ArrayList<Double>[])new ArrayList[3];
		HashMap<Double, Integer> map = new HashMap<Double, Integer>();
		
		//create a HashMap: key = time when substitutions occurred, value = site where substitutions occurred
		//do not add the last time in the HashMap b/c they will be the same for all sites.
		//i.e. last time = branch length
		for (int i = 0; i < 3; i ++) {
			state[i] = trees[i].getNodeByNodeNum(branchIndex).getPathState();
			time[i] = cumulativeSum(trees[i].getNodeByNodeNum(branchIndex).getPathTime());
			for (int j = 0; j < time[i].size()-1; j ++) {
				map.put(time[i].get(j), i);
			}
		}
		
		//starting state (triplet)
		int[] index = new int[3]; //track the index of states for three sites
		pathState[branchIndex].add(state[0].get(index[0])+state[1].get(index[1])+state[2].get(index[2]));
		
		//sort the HashMap by key => create a TreeMap (ordered) from HashMap
        Map<Double, Integer> treemap = new TreeMap<Double, Integer>(map);
        
        //iterate the TreeMap
        Set set = treemap.entrySet();
        Iterator iterator = set.iterator();
        double previous = 0.0;
        while(iterator.hasNext()) {
             Map.Entry me = (Map.Entry)iterator.next();
             //System.out.print(me.getKey() + ": ");
             //System.out.println(me.getValue());
             pathTime[branchIndex].add(((Double) me.getKey() - previous)); //store time interval length but not cumulative time
             previous = (Double) me.getKey();
             int position = (Integer) me.getValue();
             index[position] ++; //increment the index of pathState[position]
             char[] lastState = pathState[branchIndex].get(pathState[branchIndex].size()-1).toCharArray();
             lastState[position] = state[position].get(index[position]).charAt(0); //replace character at position
             pathState[branchIndex].add(String.valueOf(lastState));             
        }
        //add last time interval (branch length - last event time)
        pathTime[branchIndex].add(trees[1].getBranchLength()[branchIndex]-previous);
        //add last state. i.e. state at current node
        pathState[branchIndex].add(state[0].get(state[0].size()-1)+state[1].get(state[1].size()-1)+state[2].get(state[2].size()-1));        
		
	}
	
	/**
	 * Return pathState for a given branch
	 * @param branchIndex
	 * @return
	 */
	public ArrayList<String> getPathStateAtBranch(int branchIndex) {
		//constructPath(branchIndex);
		return pathState[branchIndex];
	}
	
	/**
	 * Return pathTime for a given branch
	 * @param branchIndex
	 * @return
	 */
	public ArrayList<Double> getPathTimeAtBranch(int branchIndex) {
		//constructPath(branchIndex);
		return pathTime[branchIndex];
	} 
	
	/**
	 * Return an Arraylist which contains the cumulative sum of the given Arraylist.
	 * e.g. given {1, 2, 3, 4, 5} would return {1, 3, 6, 10, 15}
	 * @param time
	 * @return cumulative sum
	 */
	public ArrayList<Double> cumulativeSum(ArrayList<Double> time) {
		ArrayList<Double> list = new ArrayList<Double>();
		int length = time.size();
		double sum = 0;
		for (int i = 0; i < length; i ++) {
			sum += time.get(i);
			list.add(sum);
		}
		return list;
	}
	
	
	/**
	 * Assign number of changes for each type of changes on each branch.
	 * There are 18 types of changes.
	 * @return numberOfChanges
	 */
	private int[][] computeNumberOfChanges() {
		int[][] changes = new int[numberOfBranches()][18];
		for (int j = 0; j < numberOfBranches(); j++) {
			changes[j] = computeNumberOfChangesAtBranch(j);
		}		
		return changes;
	}
	
	/**
	 * Count the number of given type of changes on a given branch
	 * non-CpG site: G -> C (Index 0), C -> G (Index 1), G -> T (Index 2), C -> A (Index 3), T -> A (Index 4), A -> T (Index 5)
	 * 				 T -> G (Index 6), A -> C (Index 7), G -> A (Index 8), C -> T (Index 9), A -> G (Index 10), T -> C (Index 11)
	 * CpG site:     G -> C (Index 12), C -> G (Index 13), G -> T (Index 14), C -> A (Index 15), G -> A (Index 16), C -> T (Index 17)
	 * 
	 * @param branchIndex
	 * @return number of changes for each type on a given branch
	 */
	private int[] computeNumberOfChangesAtBranch(int branchIndex) {
		int[] count = new int[18];
		
		for (int i = 0; i < pathState[branchIndex].size() - 1; i++) {
			if(!pathState[branchIndex].get(i).contains("CG")) { //non-CpG site
				if(pathState[branchIndex].get(i).charAt(1) == 'G') {
					if(pathState[branchIndex].get(i+1).charAt(1) == 'C') {
						count[0] ++;
					} else if(pathState[branchIndex].get(i+1).charAt(1) == 'T') {
						count[2] ++;
					} else if(pathState[branchIndex].get(i+1).charAt(1) == 'A') {
						count[8] ++;
					}
				} else if(pathState[branchIndex].get(i).charAt(1) == 'C') {
					if(pathState[branchIndex].get(i+1).charAt(1) == 'G') {
						count[1] ++;
					} else if(pathState[branchIndex].get(i+1).charAt(1) == 'T') {
						count[9] ++;
					} else if(pathState[branchIndex].get(i+1).charAt(1) == 'A') {
						count[3] ++;
					}
				} else if(pathState[branchIndex].get(i).charAt(1) == 'T') {
					if(pathState[branchIndex].get(i+1).charAt(1) == 'G') {
						count[6] ++;
					} else if(pathState[branchIndex].get(i+1).charAt(1) == 'C') {
						count[11] ++;
					} else if(pathState[branchIndex].get(i+1).charAt(1) == 'A') {
						count[4] ++;
					}
				} else if(pathState[branchIndex].get(i).charAt(1) == 'A') {
					if(pathState[branchIndex].get(i+1).charAt(1) == 'G') {
						count[10] ++;
					} else if(pathState[branchIndex].get(i+1).charAt(1) == 'C') {
						count[7] ++;
					} else if(pathState[branchIndex].get(i+1).charAt(1) == 'T') {
						count[5] ++;
					}
				}
			} else { //CpG site
				if(pathState[branchIndex].get(i).charAt(1) == 'G') {
					if(pathState[branchIndex].get(i+1).charAt(1) == 'C') {
						count[12] ++;
					} else if(pathState[branchIndex].get(i+1).charAt(1) == 'T') {
						count[14] ++;
					} else if(pathState[branchIndex].get(i+1).charAt(1) == 'A') {
						count[16] ++;
					}
				} else if(pathState[branchIndex].get(i).charAt(1) == 'C') {
					if(pathState[branchIndex].get(i+1).charAt(1) == 'G') {
						count[13] ++;
					} else if(pathState[branchIndex].get(i+1).charAt(1) == 'T') {
						count[17] ++;
					} else if(pathState[branchIndex].get(i+1).charAt(1) == 'A') {
						count[15] ++;
					}
				}
				
			}
			
		}
		
		return count;
	}
	
	/**
	 * Assign times of each state on each branch.
	 * There are 6 states: non-CpG G, non-CpG C, non-CpG T, non-CpG A, CpG G, CpG C
	 * @return timeOfStates
	 */
	private double[][] computeTimeOfStates() {
		double[][] times = new double[numberOfBranches()][6];
		for (int j = 0; j < numberOfBranches(); j++) {
			times[j] = computeTimeOfStatesAtBranch(j);
		}				
		return times;
	}
	
	/**
	 * Compute time of each state at the second position on a given branch
	 * @param branchIndex
	 * @return
	 */
	private double[] computeTimeOfStatesAtBranch(int branchIndex) {
		double[] timeState = new double[6];
		
		for (int i = 0; i < pathState[branchIndex].size()-1; i++) {
			if(!pathState[branchIndex].get(i).contains("CG")) { //non-CpG
				if(pathState[branchIndex].get(i).charAt(1) == 'G') {
					timeState[0] += pathTime[branchIndex].get(i);
				} else if(pathState[branchIndex].get(i).charAt(1) == 'C') {
					timeState[1] += pathTime[branchIndex].get(i);
				} else if(pathState[branchIndex].get(i).charAt(1) == 'T') {
					timeState[2] += pathTime[branchIndex].get(i);
				} else if(pathState[branchIndex].get(i).charAt(1) == 'A') {
					timeState[3] += pathTime[branchIndex].get(i);
				} 
			} else {	//CpG
				if(pathState[branchIndex].get(i).charAt(1) == 'G') {
					timeState[4] += pathTime[branchIndex].get(i);
				} else if(pathState[branchIndex].get(i).charAt(1) == 'C') {
					timeState[5] += pathTime[branchIndex].get(i);
				}
			}
		}
		
		return timeState;
	}
	
	/**
	 * Assign branch lengths for each branch. 
	 * Branch lengths are supposed to be the same for all three sites, but in case of rounding error, I pick site 2.
	 * @return branchLengths
	 */
	private double[] computeBranchLength() {
		double[] br = new double[numberOfBranches()];
		for (int i = 0; i < trees[1].getNumBranches(); i++) {
			br[i] = trees[1].getBranchLength()[i];
		}
		return br;
	}
	
	/**
	 * Compute the proportion of time that the second position is a non-CpG G, non-CpG C, non-CpG T, non-CpG A, CpG G, CpG C site.
	 * @return prop
	 */
	private double[][] computePropState() {
		double[][] prop = new double[numberOfBranches()][6];
		for (int i = 0; i < numberOfBranches(); i++) {
			for (int j = 0; j < 6; j++) {
				prop[i][j] = timeOfStates[i][j]/branchLengths[i];
			}
		}		
		return prop;
	}
	
	/**
	 * Get proportion of time the second position is a CpG / non-CpG site for each branch. 
	 * dim = (# of branches) x 6
	 * @return propStates
	 */
	public double[][] getPropStates() {
		return propStates;
	}
	
    /**
     * Set proportion of states
     * @param propStates to be set
     */
	public void setPropStates(double[][] propStates) {
		this.propStates = propStates;
	}

    /**
     * Return the number of changes for each substitution type (scenario 1) on branches
     * @return numberOfChanges
     */
	public int[][] getNumberOfChanges() {
		return numberOfChanges;
	}
	
	
    /** 
     * Set the number of changes for each substitution type (scenario 1) on branches
     * @param numberOfChanges
     */
	private void setNumberOfChanges(int[][] numberOfChanges) {
		this.numberOfChanges = numberOfChanges;
	}

    /**
     * Return the total time in each states on each branch
     * @return timeOfStates
     */
	public double[][] getTimeOfStates() {
		return timeOfStates;
	}
	
    /**
     * Set the total time in each states on each branch
     * @param timeOfStates
     */
	private void setTimeOfStates(double[][] timeOfStates) {
		this.timeOfStates = timeOfStates;
	}

    /**
     * Return the branch length on each branch
     * @return branchLengths
     */
	public double[] getBranchLengths() {
		return branchLengths;
	}
	
    /**
     * Set the branch length on each branch
     * @param branchLengths
     */
	private void setBranchLengths(double[] branchLengths) {
		this.branchLengths = branchLengths;
	}
	
	/**
	 * Set rootStateCount by the state of the second position
	 */
	private void setRootState() {
		String s = trees[1].getRoot().getState();
		if(s.equals("G")) {
			this.rootState = 0;
		} else if (s.equals("C")) {
			this.rootState = 1;
		} else if (s.equals("T")) {
			this.rootState = 2;
		} else if (s.equals("A")) {
			this.rootState = 3;
		}
		
	}
	
	public int getRootState() {
		return rootState;
	}

}

