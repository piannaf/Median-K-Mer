package motif;

/**
 * A class for representing a string of k symbols.
 */
public class KMer {

	private final Alphabet alpha; // the alphabet from which symbols are drawn
	private int[] kmer; // the string (by index to alphabet)

	/**
	 * Constructs an empty k-mer
	 * 
	 * @param alpha
	 *            the alphabet from which valid symbols are drawn
	 * @param k
	 *            the length of the k-mer
	 */
	public KMer(Alphabet alpha, int k) {
		this.alpha = alpha;
		this.kmer = new int[k];
		for (int i = 0; i < k; i++)
			this.kmer[i] = -1; // mark each element with "no-symbol"
	}
	
	/**
	 * Constructs a k-mer from another k-mer starting at the beginning to the
	 * specified length. If k is larger than the length of the source, only
	 * up to the length of the source will be copied.
	 * 
	 * @param kmer
	 * 			the source k-mer
	 * 
	 * @param k
	 * 			the length of the new k-mer
	 */
	public KMer(KMer kmer, int k) {
		/* Use the same alphabet */
		this.alpha = kmer.alpha;
		
		/* get the minimum length */
		k = k > kmer.getK() ? kmer.getK() : k;
		
		/* construct the new kmer */
		this.kmer = new int[k];
		for (int i = 0; i < k; i++) {
			this.kmer[i] = kmer.kmer[i];
		}
	}

	/**
	 * Constructs a k-mer from a sub-sequence
	 * 
	 * @param seq
	 *            the complete sequence
	 * @param pos
	 *            the start position of the sub-sequence
	 * @param k
	 *            the length (width) of the sub-sequence
	 */
	public KMer(DNASequence seq, int pos, int k) {
		this.alpha = seq.getAlphabet();
		this.kmer = new int[k];
		for (int i = 0; i < k; i++)
			this.kmer[i] = seq.getSymbolIndex(pos + i);
	}

	/**
	 * Constructs a fully specified k-mer
	 * 
	 * @param alpha
	 *            the alphabet from which valid symbols are drawn
	 * @param kmer
	 *            an array of indices to symbols
	 */
	public KMer(Alphabet alpha, int[] kmer) {
		this.alpha = alpha;
		this.kmer = kmer;
	}

	/**
	 * Compares two k-mers (this and other)
	 * 
	 * @param other
	 *            the other k-mer (to compare with)
	 * @return true if the two k-mers are identical in regard to string AND
	 *         alphabet, false otherwise
	 */
	public boolean equals(KMer other) {
		if (other.alpha.equals(this.alpha)) {
			for (int i = 0; i < kmer.length; i++)
				if (other.kmer[i] != this.kmer[i])
					return false;
			return true;
		}
		return false;
	}

	/**
	 * Retrieves the string of symbol indices
	 * 
	 * @return the string of symbol indices
	 */
	public int[] getKMer() {
		return kmer;
	}

	/**
	 * Checks if the k-mer is completely specified, i.e. has no empty positions.
	 * 
	 * @return true if all symbols are specified (that is, no "no-symbols": -1),
	 *         false otherwise
	 */
	public boolean isComplete() {
		if (kmer[kmer.length - 1] < 0)
			return false;
		return true;
	}

	/**
	 * Retrieves the number of symbols that are specified in a possibly partial
	 * k-mer (until the first -1 element is observed). If k-mer is partial, this
	 * number will be less than k.
	 * 
	 * @return the number of symbols in the k-mer
	 * @see #getK() to retrieve the number of positions in the k-mer including
	 *      unspecified (-1) elements
	 */
	public int getLevel() {
		for (int i = 0; i < kmer.length; i++)
			if (kmer[i] < 0)
				return i;
		return kmer.length;
	}

	/**
	 * @return the length of the k-mer including unspecified elements.
	 * @see #getLevel() to retrieve the number of specified elements in a
	 *      partial k-mer.
	 */
	public int getK() {
		return kmer.length;
	}

	/**
	 * Make an array-based clone of the string (including empty elements)
	 * 
	 * @return a clone of the k-mer
	 */
	private int[] cloneKMer() {
		int[] nkmer = new int[kmer.length];
		for (int i = 0; i < kmer.length; i++)
			nkmer[i] = kmer[i];
		return nkmer;
	}

	/**
	 * Printable string representation of the k-mer instance
	 */
	public String toString() {
		StringBuffer sbuf = new StringBuffer();
		sbuf.append(alpha.toChar(kmer));
		return sbuf.toString();
	}

	/**
	 * Expands the current, partial k-mer into all k-mers that can be
	 * constructed by adding a single symbol from the alphabet. The order of
	 * k-mers follows the order defined for the alphabet.
	 * 
	 * @return an array of k-mers extending the current by one symbol
	 */
	public KMer[] getExtensions() {
		if (!this.isComplete()) {
			KMer[] children = new KMer[alpha.getSize()];
			int level = getLevel();
			for (int i = 0; i < alpha.getSize(); i++) {
				int[] nkmer = cloneKMer();
				nkmer[level] = i;
				children[i] = new KMer(alpha, nkmer);
			}
			return children;
		}
		return null;
	}

	/**
	 * Creates all (overlapping) k-mers that are found in a DNA sequence.
	 * 
	 * @param seq
	 *            the sequence
	 * @param k
	 *            the length of k-mers
	 * @return an array of k-mers
	 */
	public static KMer[] factory(DNASequence seq, int k) {
		KMer[] kmers = new KMer[seq.getLength() - k + 1];
		for (int i = 0; i < seq.getLength() - k + 1; i++)
			kmers[i] = new KMer(seq, i, k);
		return kmers;
	}
	
	public KMer pushSymbol(int symbol) {
		int[] clone = this.cloneKMer();
		//System.out.println(this.isComplete());
		if (!this.isComplete()) {
			clone[this.getLevel()] = symbol;
			//System.out.println(symbol);
		}
		return new KMer(this.alpha, clone);
	}

}
