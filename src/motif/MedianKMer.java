package motif;

import java.io.IOException;
import java.util.Arrays;

/**
 * Class for searching for a "median" k-mer in DNA sequence data.
 */
public class MedianKMer {

	private final DNASequence[] seqs; // the sequence data that are searched
	private final Alphabet alpha; // the alphabet that each of the sequences uses

	/**
	 * Constructs an instance of the class, prepares for search and checks for
	 * problems. Do not modify the "signature" of this constructor. However, you
	 * may change the code it contains.
	 * 
	 * @param seqs
	 *            the sequence data
	 */
	public MedianKMer(DNASequence[] seqs) {
		this.seqs = seqs;
		Alphabet tmp = null;
		for (DNASequence seq : seqs)
			if (tmp == null)
				tmp = seq.getAlphabet();
			else if (tmp != seq.getAlphabet())
				throw new RuntimeException(
						"Sequences are using different alphabets");
		this.alpha = tmp;
	}

	/**
	 * Start search for median k-mer of a specified length You may NOT modify
	 * the "signature" of this method. However, you may change the code it
	 * contains.
	 * 
	 * @param k
	 *            the length of the sought k-mer
	 * @return the minimum distance and the median k-mer that rendered that
	 *         distance
	 */
	public Distance findMedianKMer(int k) {
		TrieKMer trie = new TrieKMer(seqs, k);
		
		//TESTS
		//KMer p = new KMer(alpha, new int[] {3, 1, 2, 2, 3}); // 'TCGGT'
		//KMer s = new KMer(alpha, new int[] {1}); // 'C'
		//KMer p = new KMer(alpha, new int[] {0, 2}); // 'AG'
		//KMer s = new KMer(alpha, new int[] {2}); // 'G'
		//
		//int[][] ds = getExtDists(p, s);
		//System.out.println(Arrays.deepToString(ds) + " -> " + getMedianDistance(ds));
		//END TESTS
		
		return findMedianKMer(null, k, null); // Naive Branch and Bound
		//return findMedianKMer(null, null, k, null, trie); // Trie Branch and Bound
	}
	
	/**
	 * Sample as below but using a Trie
	 */
	public Distance findMedianKMer(KMer prefix, int[][] p_dist, int k, Distance old_best, TrieKMer trie) {
		if (prefix == null)
			prefix = new KMer(alpha, k); // empty k-mer

		int distance;
		if (p_dist == null) {
			distance = 0;
			/* Find distance for unknown prefix */
			for (DNASequence seq : this.seqs) {
				distance += getDistance(seq, prefix);
			}
		} else {
			distance = getMedianDistance(p_dist);
		}
		
		Distance best_dist = new Distance(k * this.seqs.length + 1, prefix);

		if (prefix.isComplete()) {
			return new Distance(distance, prefix);
		}
		
		if( old_best!=null) {
			if (distance >= old_best.actual) {
				return old_best;
			}	
		}
		
		//KMer[] extensions = prefix.getExtensions();
		/* Sort extensions based on frequency in trie */
		KMer[] extensions = sortedExtensions(prefix, trie);
		
		for (KMer kmer : extensions) {	
			if( old_best!=null) {
				if (best_dist.actual > old_best.actual)
					best_dist = old_best;
			}
			
			
			Distance temp_dist = findMedianKMer(kmer, null, k, best_dist, trie);
			if (temp_dist.actual < best_dist.actual) {
				best_dist = temp_dist;
			}
		}

		return best_dist;
	}
	
	/** 
	 * Get the resulting distances from extending a prefix with a suffix 
	 */
	public int[][] getExtDists(KMer prefix, KMer suffix) {
		/* Initialize everything */
		int prefix_size = prefix.getLevel();
		int[][] prefix_dists = new int[seqs.length][];
		int[][] suffix_dists = new int[seqs.length][];
		for (int i = 0; i < seqs.length; i++) {
			prefix_dists[i] = getDistances(seqs[i], prefix);
			suffix_dists[i] = getDistances(seqs[i], suffix);
		}
		
		/* How many distance arrays to compare */
		int outer_len = prefix_dists.length;
		int inner_len = prefix_dists[0].length;
		
		/* the result */
		int[][] result = new int[outer_len][inner_len - 1];
		
		for (int i = 0; i < outer_len; i++) {
			//System.out.println(Arrays.toString(prefix_dists[i]));
			//System.out.println(Arrays.toString(suffix_dists[i]));
			
			/* Loop over suffix distance array */
			for (int j = prefix_size; j < inner_len + prefix_size - 1; j++) {
				//System.out.print(suffix_dists[i][j] + "+" + prefix_dists[i][j-prefix_size] + " ");
				result[i][j-prefix_size] = suffix_dists[i][j] + prefix_dists[i][j-prefix_size];
			}
			//System.out.println();
		}
		return result;
	}
	
	public KMer[] sortedExtensions(KMer prefix, TrieKMer trie) {
		KMer[] sorted = prefix.getExtensions().clone();
		
		for (int i = 1; i < sorted.length; i++){
			  int j = i;
			  KMer current = sorted[i];
			  while ((j > 0) && 
					  (trie.getCount(sorted[j-1]) < trie.getCount(current))){
			    sorted[j] = sorted[j-1];
			    j--;
			  }
			  sorted[j] = current;
		    }
			
			return sorted;
	}

	/**
	 * Search for median k-mer of a specified length, from the given prefix. You
	 * may modify the "signature" of this method and change the code it
	 * contains.
	 * 
	 * @param prefix
	 *            the prefix from which the search is started, e.g. if prefix is
	 *            "AC", we explore "ACA", "ACC", "ACG" and "ACT" (and beyond
	 *            through recursive calls).
	 * @param k
	 *            the length of the sought k-mer, e.g. if prefix is "AC" and k
	 *            is 2, the search ends, if k is 3 we explore longer words.
	 * @return the minimum distance and the median k-mer that rendered that
	 *         distance
	 */
	public Distance findMedianKMer(KMer prefix, int k, Distance old_best) // <== you may change this definition
	{
		if (prefix == null)
			prefix = new KMer(alpha, k); // empty k-mer

		if (prefix.isComplete()) {
			int distance = 0;
			for (DNASequence seq : this.seqs) {
				distance += getDistance(seq, prefix);
			}
			
			return new Distance(distance, prefix);
		}
		Distance best_dist = new Distance(k * this.seqs.length + 1, prefix);
		
		if( old_best!=null) {
			int distance = 0;	// total distance over all DNA sequences
			
			for (DNASequence sequence : this.seqs)
				distance += getDistance(sequence, prefix);
			
			if (distance >= old_best.actual) {
				return old_best;
			}	
		}
		Distance temp_dist;
		for (KMer kmer : prefix.getExtensions()) {	
			/* Keep old_best throughout traversal */
			if( old_best!=null) {
				if (best_dist.actual >= old_best.actual)
					best_dist = old_best;
			}
			
			temp_dist = findMedianKMer(kmer, k, best_dist);
			if (temp_dist.actual <= best_dist.actual) {
				best_dist = temp_dist;
			}
		}

		return best_dist;
	}

	/**
	 * Determines the number of mismatched letters when the specified word is
	 * aligned to each position of the sequence. That is, it returns an array of
	 * the same size as the sequence (minus the positions that would extend over
	 * the end when the word is aligned), with each element containing the
	 * minimum "Hamming distance" between the sequence and the word at that
	 * position.
	 * 
	 * @param seq
	 *            the sequence that is searched
	 * @param word
	 *            the word that is aligned to the sequence
	 * @return an array with all the Hamming distances
	 */
	public static int[] getDistances(DNASequence seq, KMer word) {
		int N = seq.getLength();
		int K = word.getLevel();
		int[] distances = new int[N - K + 1];
		
		for (int i = 0; i < (N - K + 1); i++) {
			int count = 0;
			
			for (int j = 0; j < K; j++) {
				if (seq.getSymbolChar(i + j) != word.toString().charAt(j)) {
					count++;
				}
			}
			distances[i] = count;
		}
		return distances;
	}

	public static int getMedianDistance(int[][] distances) {
		int result = 0;
		
		for (int[] distance : distances) {
			result += distance[getMinPosition(distance)];
		}
		
		return result;
	}
	
	/**
	 * Determines the number of mismatched letters when the specified word is
	 * aligned optimally to the sequence. That is, the minimum
	 * "Hamming distance."
	 * 
	 * @param seq
	 *            the sequence that is searched
	 * @param word
	 *            the word that is aligned to the sequence
	 * @return the Hamming distance between the aligned word and the sequence
	 *         (extending over the length of the word)
	 */
	public static int getDistance(DNASequence seq, KMer word) {
		int N = seq.getLength();
		int K = word.getLevel();
		int mismatches = K;
		
		for (int i = 0; i < (N - K + 1); i++) {
			int count = 0;
			
			for (int j = 0; j < K; j++) {
				
				if (seq.getSymbolIndex(i + j) != word.getKMer()[j]) {
					if (mismatches < count) break;
					count++;
				}
			}
			mismatches = count < mismatches ? count : mismatches;
		}
		return mismatches;
	}

	/**
	 * Helper method that finds the position in a distance array that has the
	 * smallest distance. Do not modify the "signature" of this constructor.
	 * There is no need to modify the code it contains.
	 * 
	 * @param distances
	 * @return the first position/index with the smallest value
	 * @see #printReport(KMer)
	 */
	public static int getMinPosition(int[] distances) {
		int pos = 0;
		if (distances == null)
			throw new RuntimeException("Invalid distance array");
		if (distances.length < 1)
			throw new RuntimeException("Invalid distance array");
		int min = distances[pos];
		for (int i = 1; i < distances.length; i++) {
			if (distances[i] < min) {
				min = distances[i];
				pos = i;
			}
		}
		return pos;
	}

	/**
	 * Retrieves the alphabet that applies to this class.
	 * 
	 * @return the alphabet
	 */
	public Alphabet getAlphabet() {
		return alpha;
	}

	/**
	 * A method that prints some stats relating to a specific k-mer when used
	 * against the sequence data set provided to the constructor.
	 * 
	 * @param kmer
	 */
	public void printReport(KMer kmer) {
		System.out.println("REPORT for " + kmer);
		int[][] counts = new int[kmer.getK()][alpha.getSize()];
		int total = 0;
		for (int i = 0; i < seqs.length; i++) {
			int[] d = getDistances(seqs[i], kmer);
			int pos = getMinPosition(d);
			total += d[pos];
			KMer found = new KMer(seqs[i], pos, kmer.getK());
			System.out.println(found + "\t" + d[pos] + "\t@ " + pos + "\tin "
					+ seqs[i]);
			for (int j = 0; j < found.getK(); j++)
				counts[j][found.getKMer()[j]]++;
		}
		System.out.println("Distance: " + total);
		System.out.println("Counts (can use as TomTom input): ");
		for (int i = 0; i < alpha.getSize(); i++) {
			for (int j = 0; j < counts.length; j++)
				System.out.print(String.format("%5d", counts[j][i]));
			System.out.println();
		}
	}

	/**
	 * A command line application that accepts a number of parameters. 
	 * -f <filename> 
	 * -k <length-of-k-mer> 
	 * -q <k-mer> 
	 * See usage message for more information.
	 * 
	 * @param args
	 */
	public static void main(String[] args) {

		int k = 10; // default length of k-mer
		String file = null; // FASTA file
		String kstr = null; // user-specified k-mer
		DNASequence[] seqs = null; // loaded sequences

		// parse the parameters
		for (int i = 0; i < args.length; i++) {
			if (args[i].charAt(0) == '-') // option
			{
				switch (args[i].charAt(1)) {
				case 'k':
					if (i + 1 < args.length)
						k = Integer.parseInt(args[++i]);
					break;
				case 'f':
					if (i + 1 < args.length)
						file = args[++i];
					break;
				case 'q':
					if (args[i].length() > 2)
						if (args[i].charAt(2) == 's') {
						}
					if (i + 1 < args.length)
						kstr = args[++i];
					break;
				default:
					System.err.println("Unknown option \"-" + args[i].charAt(1)
							+ "\"");
				}
			}
		}

		if (file == null && kstr == null) {
			System.err
					.println("Usage: MedianKMer -f <filename> { -q <query-k-mer> | -k <length> }");
			System.err.println("where <filename> is a FASTA file");
			System.err
					.println("-q will report on matches with the specified <query-k-mer> (-qs produces a short report)");
			System.err
					.println("-k will search for the best median k-mer where k=<length>");
			System.exit(1);
		}

		if (file != null) {
			try {
				seqs = DNASequence.readFile(new Alphabet(), file); // read a
				// FASTA
				// file with
				// sequences
			} catch (IOException e) {
				System.err.println(e.getMessage());
				System.exit(2);
			}
		}

		System.out.println("Started at "
				+ new java.util.Date(System.currentTimeMillis()));
		MedianKMer ms = null;

		if (seqs != null)
			ms = new MedianKMer(seqs);
		else {
			System.err
					.println("Usage: MedianKMer -f <filename> { -q[s] <query-k-mer> | -k <length> }");
			System.exit(3);
		}

		if (kstr != null) // query
		{
			ms.printReport(new KMer(ms.getAlphabet(), ms.getAlphabet().toIndex(
					kstr.toCharArray())));
		} else // search
		{
			Distance dist = ms.findMedianKMer(k); // start searching
			System.out.println(dist); // print result
		}

		System.out.println("Ended at "
				+ new java.util.Date(System.currentTimeMillis()));
	}

	/**
	 * Holder of score and the path leading to those scores. You may modify the
	 * code for this but keep the original constructor signature.
	 */
	class Distance {

		final int actual; // the actual distance of this k-mer
		final KMer path; // the k-mer to which the distance applies

		/**
		 * Constructs an instance that combines a distance and the applicable
		 * k-mer.
		 * 
		 * @param actual
		 *            the actual distance
		 * @param kmer
		 *            the k-mer
		 */
		public Distance(int actual, KMer kmer) {
			this.actual = actual;
			this.path = kmer;
		}

		public String toString() {
			StringBuffer sbuf = new StringBuffer();
			sbuf.append(path.toString() + ":" + actual);
			return sbuf.toString();
		}
	}
}