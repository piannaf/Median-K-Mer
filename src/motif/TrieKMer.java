package motif;

import java.lang.Comparable;

import motif.MedianKMer.Distance;


/**
 * A class to build and maintain a tree representing all k-mers found in a DNA
 * sequence set
 */
public class TrieKMer {

	private final TrieNode root; // the root of the tree structure
	private final Alphabet alpha; // the alphabet from which all k-mers are
									// constructed
	
	/**
	 * Method for adding a k-mer to the tree structure. Adds a count of one to
	 * EACH applicable node (including non-leaf nodes). Adds nodes to the tree
	 * if the k-mer has not been added before. You should change the code of
	 * this method. You may change the signature if you need to.
	 * 
	 * @param kmer
	 *            the k-mer to add
	 *            
	 * @param depth
	 *			  the maximum depth
	 */
	public void putKMer(TrieNode parent, KMer kmer) {
		int[] symbols = kmer.getKMer();
		for (int i = 0; i < symbols.length; i++) {
			TrieNode current = parent.children[symbols[i]];
			
			if (current == null) {
				current = new TrieNode(this.alpha);
				current.symbol = symbols[i];
				parent.children[symbols[i]] = current;
			}
			
			current.count++;
				
			parent = current;
		}
	}
	
	/**
	 * Helper method for constructing a TrieKMer object
	 * 
	 * @param DNA
	 * 			The DNA sequence whose k-mers will be added to the trie
	 * 
	 * @param K
	 * 			The length of each k-mer to add to the tree
	 */
	public void buildTrie(DNASequence DNA, int K) {
		/* 
		 * The number of possible k-mers for the given length is
		 * |DNA| - K
		 */
		for (int i = 0; i <= DNA.getLength() - K; i++) {
			/* Get next kmer */
			KMer kmer = new KMer(DNA, i, K);
			
			/* put kmer in trie */
			putKMer(this.root, kmer);
		}
	}

	public int[] sortedChildren(KMer prefix) {
		/* The path to follow */
		int[] path = prefix.getKMer();
		
		/* The current node of the trie */
		TrieNode current = this.root;
		
		/* The depth which has been traveled */
		int depth = 0;
		
		/*
		 * Continue traversing the path until either we have traversed the
		 * whole path or we have reached a leaf node
		 */
		while(depth < prefix.getLevel() && current.children[path[depth]] != null) {
			/* go to the correct child and increase the depth*/
			current = current.children[path[depth++]];
		}
		
		if (depth > prefix.getLevel()) {
			/* Followed the path but didn't find prefix */
			return null;
		} else {
			/* Found prefix ending at current node */
			TrieNode[] sorted = current.sortedChildren();
			int[] result = new int[sorted.length];
			
			for (int i = 0; i < sorted.length; i++) {
				result[i] = sorted[i].symbol;
			}
			return result;
		}
	}
	
	
	/**
	 * Method for querying the tree structure. Finds a node representing the
	 * k-mer and returns the count. Should return 0 if the k-mer is not found.
	 * Note that the k-mer may be shorter than the depth of the tree. In that
	 * case it returns the count of the node representing the final symbol of
	 * the k-mer.
	 * 
	 * You should change the code of this method. You may change the signature
	 * if you need to.
	 * 
	 * @param kmer
	 *            the query k-mer
	 * @return the number of times this exact k-mer has been observed
	 */
	public int getCount(KMer kmer) {
		/* The path t0 follow */
		int[] path = kmer.getKMer();
		
		/* The current node of the trie*/
		TrieNode current = this.root;
		
		/* 
		 * The count.
		 * If there is nothing in the tree then it's not possible to find the 
		 * kmer so default to 0
		 */
		int count = 0;
		
		/* The depth which has been traveled */
		int depth = 0;
		
		/*
		 * Continue traversing the path until either we have traversed the
		 * whole path or we have reached a leaf node
		 */
		while(depth < kmer.getLevel() && current.children[path[depth]] != null) {
			/* go to the correct child and increase the depth*/
			current = current.children[path[depth++]];
		}
		
		if (depth < kmer.getLevel()) {
			/* Followed the path but didn't find kmer */
			return 0;
		} else {
			/* Found kmer ending at current node */
			count = current.count;
		}
		
		return count;
	}

	/**
	 * Extracts all sub-sequences of the specified length (depth) and then
	 * constructs a trie representing all of them. Note that the counts in the
	 * trie should correspond to the number of times a letter has occurred in
	 * each position of a full sub-sequence, e.g. if the sequence is "GGAGCACA"
	 * and k=2, then the resulting trie will be ((A:2(C:2) C:2(A:2) G:2(A:1
	 * G:1)) It should make use of putKMer.
	 * 
	 * @param seqs
	 *            the sequences
	 * @param depth
	 *            the depth of the tree and the maximum length of the counted
	 *            k-mers
	 * @see #putKMer(KMer)
	 */
	public TrieKMer(DNASequence[] seqs, int depth) {
		if (seqs.length < 1) // no sequences are provided
			this.alpha = null; // hence, no alphabet can be identified
		else {
			// check and set the alphabet
			this.alpha = seqs[0].getAlphabet();
			for (int i = 0; i < seqs.length; i++) {
				if (!alpha.equals(seqs[i].getAlphabet()))
					throw new RuntimeException("Alphabets are not equal");
			}
		}
		
		/* The root of the tree begins life as a leaf */
		root = new TrieNode(this.alpha);
		
		/* Build up the tree sequence by sequence */
		for (DNASequence seq : seqs) {
			buildTrie(seq, depth);
		}
	}
}

/**
 * A class for representing node instances (each representing a symbol that is
 * part of a k-mer).
 */
class TrieNode {
	public int symbol;	// The letter of k-mer stored in this node
	public int count;	// The frequency of the prefix stored in this node
	public TrieNode[] children;
	public Alphabet alpha;
	
	public TrieNode(Alphabet alpha) {
		symbol = -1;	// No letter is stored here
		count = 0;
		children = new TrieNode[alpha.getSize()];
		this.alpha = alpha;
	}
	
	public boolean isLeaf() {
		for(TrieNode child : children) {
			if (child != null) {
				return false;
			}
		}
		return true;
	}
	
	public String toString() {
		return "" + count;
	}
	
	public TrieNode[] sortedChildren() {
		TrieNode[] sorted = new TrieNode[this.children.length];
		for (int i = 0; i < sorted.length; i++) {
			if (this.children[i] == null) {
				sorted[i] = new TrieNode(this.alpha);
			} else {
				sorted[i] = this.children[i];
			}
		}
		
		for (int i = 1; i < sorted.length; i++){
		  int j = i;
		  TrieNode current = sorted[i];
		  while ((j > 0) && (sorted[j-1].count > current.count)){
		    sorted[j] = sorted[j-1];
		    j--;
		  }
		  sorted[j] = current;
	    }
		
		return sorted;
	}
}
