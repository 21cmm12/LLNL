package main;
import java.util.Scanner;
public class NZ {

	public static int NUM_TRIALS = 100000; //How many simulation trials you want to run (can be changed by user)
	public static int L = 30; // Linear dimension (can be changed by user)
	static int numNodes = 0;

	public static final int NUM_NEIGHBORS = 8; //4 for square lattice, change this for different lattice shapes
	public static final int DIM_LATTICE = 2; //Change this for different dimensional latties
	public static int numVert = L * L; //Change this for different lattice shapes

	int[] ptr = new int[numVert]; // ptr[i] is the root vertex of the cluster containing i
	int[] sizeClus = new int[numVert]; // sizeClus[i] is the size of the cluster with root i (nonsense at non-root points)
	int[][] disp = new int[numVert][DIM_LATTICE]; // The displacement to the root vertex
	int[][] c0 = new int[numVert][DIM_LATTICE]; // The wrapping vector, c0, of the root of a given cluster (nonsense at non-root points)
	int numEdges = NUM_NEIGHBORS * numVert / 2;
	int[][] edgeList = new int[numEdges][2 + DIM_LATTICE]; //edgeList[i] is a particular edge (v,w) followed by its x and y displacement from v to w
	int[] order = new int[numEdges]; // Occupation order



	public NZ(double[] vArr, double[] qArr) {
		Node head = new Node();
		boundaries();
		for (int i = 0; i < NUM_TRIALS; i++) {
			permutation();
			percolate(head);
		}
		for(int i = 0; i < vArr.length; i++) {
			summation(head, vArr[i], qArr[i]);
		}
	}

	public void boundaries() { // Change this function only (and NUM_NEIGHBORS, DIM_LATTICE, numVert) for different lattice shapes... right now it is square matching
		for (int i = 0; i < numVert; i++) {
			edgeList[4 * i][0] = i;
			edgeList[4 * i][1] = (i + 1) % numVert;
			edgeList[4 * i][2] = 1;
			edgeList[4 * i][3] = 0; //Adds the edge between this vertex and the vertex to the right
			edgeList[4 * i + 1][0] = i;
			edgeList[4 * i + 1][1] = (i + L) % numVert;
			edgeList[4 * i + 1][2] = 0;
			edgeList[4 * i + 1][3] = 1; //Adds the edge between this vertex and the vertex below
			edgeList[4 * i + 2][0] = i;
			edgeList[4 * i + 2][1] = (i + L + 1) % numVert;
			edgeList[4 * i + 2][2] = 1;
			edgeList[4 * i + 2][3] = 1; //Adds the edge between this vertex and the vertex below and to the right
			edgeList[4 * i + 3][0] = i;
			edgeList[4 * i + 3][1] = (i + L - 1) % numVert;
			edgeList[4 * i + 3][2] = -1;
			edgeList[4 * i + 3][3] = 1; //Adds the edge between this vertex and the vertex below and to the left
			//if ((i + 1) % L == 0)
			//	edgeList[2 * i][1] = i - L + 1; //This edge case was in NZ code but I don't think it is needed
		}
	}


	public void permutation() { //All this does is makes a random order of 1-numEdges.
		for (int i = 0; i < numEdges; i++)
			order[i] = i;
		for (int i = 0; i < numEdges; i++) {
			int j = (int) (i + (numEdges - i) * Math.random());
			int temp = order[i];
			order[i] = order[j];
			order[j] = temp;
		}
	}

	public int[] findroot(int i) {//The return value will be {xptr, yptr, root}
		if (ptr[i] == i) {
			int[] ans = new int[DIM_LATTICE + 1];
			ans[DIM_LATTICE] = i;
			return ans;
		}
		int[] ans = findroot(ptr[i]);//Also path compresses on the way down
		for(int j = 0; j < DIM_LATTICE; j++) {//Including path compressing pointers
			ans[j] += disp[i][j];
			disp[i][j] = ans[j];
		}
		ptr[i] = ans[DIM_LATTICE];
		return ans;
	}

	public void percolate(Node head) {
		//AFunction aFunction = (x, y) -> {
		//	return true;
		//};// Replace this with the 0D or 2D function, or any other A function
		Node current = head;
		int s1, s2 = 0;
		int r1, r2;
		int maxClusterDim = 0;
		//int big = 0;
		for (int i = 0; i < numVert; i++) {
			ptr[i] = i;
			for (int j = 0; j < DIM_LATTICE; j++) {
				disp[i][j] = 0;
				c0[i][j] = 0;
			}
		}
		for (int i = 0; i < numEdges; i++) {
			int edgeNum = order[i];
			//boolean horizontal = ((edgeNum % 2) == 0); //Taken out to make this work with new lattice shapes
			s1 = edgeList[edgeNum][0];
			s2 = edgeList[edgeNum][1];
			r1 = findroot(s1)[DIM_LATTICE];
			r2 = findroot(s2)[DIM_LATTICE];
			boolean combinedTwo = (r1 != r2);
			if (combinedTwo) { //Union the two clusters
				boolean oldRootWrapped = false;
				if (sizeClus[r1] > sizeClus[r2]) {//Then r1 will be the root of the new cluster
					ptr[r2] = r1;
					sizeClus[r1] += sizeClus[r2];
					//Now we update displacement pointers from s2 to s1, the new root
					for (int j = 0; j < DIM_LATTICE; j++) {
						disp[r2][j] = disp[s1][j] - disp[s2][j] - edgeList[edgeNum][2 + j];
						if (c0[r1][j] != 0) {
							oldRootWrapped = true;
						}
					}

					if (oldRootWrapped) {//If the cluster being joined wrapped, the new root gets the wrapping
						for (int j = 0; j < DIM_LATTICE; j++) {
							c0[r1][j] = c0[r2][j];
						}
					}
				} else {
					ptr[r1] = r2;
					sizeClus[r2] += sizeClus[r1];
					for (int j = 0; j < DIM_LATTICE; j++) {
						disp[r1][j] = disp[s2][j] - disp[s1][j] + edgeList[edgeNum][2 + j];
						if (c0[r2][j] != 0) {
							oldRootWrapped = true;
						}
					}
					if (oldRootWrapped) {
						for (int j = 0; j < DIM_LATTICE; j++) {
							c0[r2][j] = c0[r1][j];
						}
					}
				}
				//if (sizeClus[r1] > big)
				//	big = sizeClus[r1]; //Can check for largest cluster size
			} else {
				if(maxClusterDim == 0) {
					for (int j = 0; j < DIM_LATTICE; j++) {
						int dispCoord = disp[s2][j] - disp[s1][j] + edgeList[edgeNum][2 + j];
						if (dispCoord != 0) {
							maxClusterDim = 1;
						}
					}
					if (maxClusterDim == 1) {//If we just wrapped for the first time, store the c0
						for (int j = 0; j < DIM_LATTICE; j++) {//Better that this is a seperate loop so we aren't remembering dispCoords when this doesn't trigger.
							c0[r1][j] = disp[s2][j] - disp[s1][j] + edgeList[edgeNum][2 + j];
						}
					}
				} else if (maxClusterDim == 1) {//This part (P(2D)) only works for 2-dimensional planar lattice, and probably not even then :(
					int crossProd = c0[r1][0] * (disp[s2][1] - disp[s1][1] + edgeList[edgeNum][2 + 1]) - c0[r1][1] * (disp[s2][0] - disp[s1][0] + edgeList[edgeNum][2 + 0]);
					if (crossProd != 0) {
						maxClusterDim = 2;
					}
				}
			}
			current = operate(current, combinedTwo, (maxClusterDim == 0)); //aFunction.A(x, y)); //Change this to what we are testing for
		}
	}

	public Node operate(Node current, boolean combinedTwo, boolean aHappened) {
		if (combinedTwo) {// Now we move to (n+1,C-1)
			if (current.rightAndDown != null) {
				current = current.rightAndDown;
			} else {
				current = new Node(current, Pos.UP);
			}
		} else {
			if (current.right != null) {
				current = current.right;
			} else {
				current = new Node(current, Pos.LEFT);
			}
		}
		current.happen(aHappened);
		return current;
	}

	public void summation(Node head, double v, double q) {
		double sumOfB = 0.0;
		double sumOfBA = 0.0;
		Node leftEdge = head;
		double p = v / (1 + v);
		int nMax = (int) (numVert * p * (1 - (1 - p) * L * L / numVert * Math.log(q)));
		// get to nMax, cMax
		int n = 0;
		for (n = 0; n < nMax; n++) {
			if (leftEdge.right != null)
				leftEdge = leftEdge.right;
			else
				leftEdge = leftEdge.rightAndDown;
		}

		leftEdge.b = 1.0;
		Node start = leftEdge;

		while (leftEdge.right != null || leftEdge.rightAndDown != null) {
			Node traversing = leftEdge;
			while (traversing.hasDown()) {
				traversing.down().b = traversing.b / q;
				traversing = traversing.down();
			}
			n++;
			if (leftEdge.right != null) {
				leftEdge = leftEdge.right;
				leftEdge.b = leftEdge.left.b * (numVert - n + 1) / n * v;
			} else {
				leftEdge = leftEdge.rightAndDown;
				leftEdge.b = leftEdge.leftAndUp.b * (numVert - n + 1) / n * v / q;
			}
			
		}

		leftEdge = start;
		n = nMax;

		while (leftEdge.hasDown())
			leftEdge = leftEdge.down();

		while (leftEdge.left != null || leftEdge.leftAndUp != null) {
			Node traversing = leftEdge;
			while (traversing.hasUp()) {
				traversing.up().b = traversing.b * q;
				traversing = traversing.up();
			}
			n--;
			if (leftEdge.left != null) {
				leftEdge = leftEdge.left;
				leftEdge.b = leftEdge.right.b * (n + 1) / (numVert - n) / v;
			} else {
				leftEdge = leftEdge.leftAndUp;
				leftEdge.b = leftEdge.rightAndDown.b * (n + 1) / (numVert - n) / v * q;
			}
			
		}

		leftEdge = head;
		while (leftEdge != null) {
			while (leftEdge != null && leftEdge.rightAndDown == null) {
				// Process leftEdge
				sumOfB += leftEdge.b * leftEdge.numVisits;
				sumOfBA += leftEdge.b * leftEdge.numTimesA;
				leftEdge = leftEdge.right;
			}
			Node traversing = leftEdge;
			while (traversing != null) {
				// Process traversing
				sumOfB += traversing.b * traversing.numVisits;
				sumOfBA += traversing.b * traversing.numTimesA;
				traversing = traversing.right;
			}
			if(leftEdge != null)
			leftEdge = leftEdge.rightAndDown;
		}
		System.out.println("P(A) is about " + (sumOfBA / sumOfB) + " for p = " + (v / (1 + v)) + " and q = " + q);
	}
	
	





	public static void main(String[] args) {
		Scanner s = new Scanner(System.in);
		System.out.println("Input L Value: ");
		String strL = s.nextLine();
		if (strL == "") {
			System.out.println("Assuming L = " + L);
		} else {
			try {
				L = Integer.parseInt(strL);
			} catch (NumberFormatException e) {
				System.err.println("Error: L must be a positive integer");
				System.exit(1);
			}
			if (L < 1) {
				System.err.println("Error: L must be a positive integer");
				System.exit(1);
			}
		}
		numVert = L * L;
		System.out.println("Input the number of trials: ");
		String strNumTri = s.nextLine();
		if (strNumTri == "") {
			System.out.println("Assuming NUM_TRIALS = " + NUM_TRIALS);
		} else {
			try {
				NUM_TRIALS = Integer.parseInt(strNumTri);
			} catch (NumberFormatException e) {
				System.err.println("Error: NUM_TRIALS must be a positive integer");
				System.exit(1);
			}
			if (NUM_TRIALS < 1) {
				System.err.println("Error: NUM_TRIALS must be a positive integer");
				System.exit(1);
			}
		}
		System.out.println("Input the number of (p, q) pairs: ");
		int NUM_PAIRS = 1;
		String strNumPar = s.nextLine();
		if (strNumPar == "") {
			System.out.println("Assuming NUM_PAIRS = " + NUM_PAIRS);
		} else {
			try {
				NUM_PAIRS = Integer.parseInt(strNumPar);
			} catch (NumberFormatException e) {
				System.err.println("Error: NUM_PAIRS must be a positive integer");
				System.exit(1);
			}
			if (NUM_PAIRS < 1) {
				System.err.println("Error: NUM_PAIRS must be a positive integer");
				System.exit(1);
			}
		}
		double[] vArr = new double[NUM_PAIRS];
		double[] qArr = new double[NUM_PAIRS];
		System.out.println("Input each (p, q) pair on a line separated by a space: ");
		for (int i = 0; i < NUM_PAIRS; i++) {
			String str = s.nextLine();
			if (str == "") {
				System.out.println("Assuming p = 0.5 and q = 1.0");
				vArr[i] = 1.0;
				qArr[i] = 1.0;
			} else {
				String[] strs = str.split(" ");
				if (strs.length < 2) {
					System.err.println("Error: input two arguments!");
					System.exit(1);
				}
				double p = 1, q = 1;
				try {
					p = Double.parseDouble(strs[0]);
				} catch (NumberFormatException e) {
					System.err.println("Error: Argument " + strs[0] + " must be a double, p.");
					System.exit(1);
				}
				if (p < 0.0 || p > 1.0) {
					System.err.println("Error: p must be between 0 and 1");
					System.exit(1);
				}
				vArr[i] = p / (1 - p);
				try {
					q = Double.parseDouble(strs[1]);
				} catch (NumberFormatException e) {
					System.err.println("Error: Argument " + strs[1] + " must be a double, q.");
					System.exit(1);
				}
				if (q <= 0.0) {
					System.err.println("Error: q must be positive");
					System.exit(1);
				}
				qArr[i] = q;
			}
		}
		System.out.println("Simulating...");
		long startTime = System.currentTimeMillis();
		new NZ(vArr, qArr);
		long endTime = System.currentTimeMillis();
		System.out.println("Nodes created: " + numNodes);
		System.out.println("Time taken: " + (endTime - startTime) + " ms");
	}

}
