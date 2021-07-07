package main;

enum Pos {
	LEFT, UP
};

public class Node {
	int numTimesA;
	int numVisits;
	public Node leftAndUp;
	public Node left;
	public Node right;
	public Node rightAndDown;
	double b = 0.0;

	private Node(Node leftNode, Node upNode) {
		this();
		left = leftNode;
		left.right = this;
		leftAndUp = upNode;
		leftAndUp.rightAndDown = this;
		fill();
	}

	public Node(Node parent, Pos position) {
		this();
		switch (position) {
		case LEFT:
			left = parent;
			left.right = this;
			break;
		case UP:
			leftAndUp = parent;
			leftAndUp.rightAndDown = this;
			break;
		}
		fill();
	}

	public Node() { NZ.numNodes++;
	} // ONLY CALL TO GENERATE THE HEAD

	private void fill() {
		if (left != null && left.rightAndDown != null) {
			if (left.rightAndDown.right != null) {
				rightAndDown = left.rightAndDown.right;
				rightAndDown.leftAndUp = this;
			} else new Node(left.rightAndDown, this);
		}
		if (leftAndUp != null && leftAndUp.right != null) {
			if (leftAndUp.right.rightAndDown != null) {
				right = leftAndUp.right.rightAndDown;
				right.left = this;
			} else new Node(this, leftAndUp.right);
		}
	}

	public void happen(boolean aHappened) {
		if (aHappened)
			numTimesA++;
		numVisits++;
	}

	public double getFraction() {
		return ((double) numTimesA) / numVisits;
	}

	public Node down() {
		if (rightAndDown == null)
			return null;
		return rightAndDown.left;
	}

	public boolean hasDown() {
		return down() != null;
	}

	public Node up() {
		if (leftAndUp == null)
			return null;
		return leftAndUp.right;
	}

	public boolean hasUp() {
		return up() != null;
	}
}
