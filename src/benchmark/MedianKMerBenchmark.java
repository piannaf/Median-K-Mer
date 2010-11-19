package benchmark;

import motif.Alphabet;
import motif.DNASequence;
import motif.MedianKMer;

import java.io.*;
import java.util.*;

public class MedianKMerBenchmark {
	final static char[] alphabet = { 'A', 'C', 'G', 'T' };
	
	public static void main(String[] args) throws IOException {
		// warm-up
		runTest(3,1);
		runTest(4,2);
		//runTest(5,7);
		
		FileWriter fw = new FileWriter("temp.csv", false);
		BufferedWriter bw = new BufferedWriter(fw);
		
		for (int n = 4; n <= 25; n++) {
			for (int k = 5; k <= 12; k++) {
				double time = runTest(n, k);
				String line = n + "," + k + "," + time;
				System.out.println(line);
				bw.write(line);
				bw.newLine();
				bw.flush();
				
				/* Don't care to do more if more than 30 seconds */
				if (time > 30000) break;
			}
		}
	}
	
	/**
	 * Runs a single benchmark
	 * @param n
	 * 			the number of sequences
	 * @param k
	 * 			the length of the k-mer
	 * @return average time to complete
	 */
	public static double runTest(int n, int k) {
		long average_time = 0;
		int iterations = 300/(n*k);	//Assume largest test will be 25*12
		
		Alphabet alpha = new Alphabet();
		DNASequence[] seqs = new DNASequence[n];
		
		for (int i = 0; i < n; i++) {
			seqs[i] = new DNASequence(alpha, "s" + i, randomSeq());
		}

		MedianKMer m = new MedianKMer(seqs);
		
		for (int i = 0; i < iterations; i++) {
			//System.out.println(i + ".\tRunning findMedianKmer for n=" + n + ", k=" + k);
			long start = System.currentTimeMillis();
			m.findMedianKMer(k);
			average_time += System.currentTimeMillis() - start;
		}
		
		
		return average_time/(double)iterations;
	}
	
	static Random generator = new Random(System.currentTimeMillis());
	
	private static char[] randomSeq() {
		char[] seq = new char[100];
		
		
		for(int i = 0; i < 100; i++) {
			seq[i] = alphabet[generator.nextInt(alphabet.length)];
		}
		//System.out.println("random: " + String.valueOf(seq));
		return seq;
	}
}
