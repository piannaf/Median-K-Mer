/**
 * 
 */
package motif;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Test;

/**
 */
public class TrieKMerTest {

	@Test
	public void testTrieKMer0() {
		Alphabet alpha = new Alphabet();
		DNASequence[] seqs = new DNASequence[] { new DNASequence(alpha, "s1",
				new char[] { 'G', 'G', 'A', 'C', 'A', 'C', 'A' }) };
		TrieKMer tkm = new TrieKMer(seqs, 2);
		// we test with a 3-mer 'AC'
		assertEquals(2, tkm.getCount(new KMer(alpha, new int[] { 0, 1 })));
	}

	@Test
	public void testTrieKMer1() {
		Alphabet alpha = new Alphabet();
		DNASequence[] seqs = new DNASequence[] {
				new DNASequence(alpha, "s1", 
						new char[] { 'T', 'C', 'G', 'G', 'A', 'C' }),
				new DNASequence(alpha, "s2", 
						new char[] { 'A', 'G', 'G', 'T', 'T', 'G' }),
				new DNASequence(alpha, "s3", 
						new char[] { 'T', 'A', 'A', 'G', 'G', 'C' }) };
		TrieKMer tkm = new TrieKMer(seqs, 4);
		// we test with a 1-mer 'A' (note that since the A in s1 does not appear
		// as the first letter in a 4-mer and is thus not counted
		assertEquals(3, tkm.getCount(new KMer(alpha, new int[] { 0 })));
	}

	@Test
	public void testTrieKMer2() {
		Alphabet alpha = new Alphabet();
		DNASequence[] seqs = new DNASequence[] { new DNASequence(alpha, "s1",
				new char[] { 'T', 'C', 'G', 'G', 'A', 'C', 'A', 'C', 'A' }) };
		TrieKMer tkm = new TrieKMer(seqs, 3);
		// we test with a 3-mer 'ACA'
		assertEquals(2, tkm.getCount(new KMer(alpha, new int[] { 0, 1, 0 })));
	}

	@Test
	public void testTrieKMer3() {
		Alphabet alpha = new Alphabet();
		try {
			DNASequence[] seqs = DNASequence.readFile(alpha,
					"data/malT_5.fasta");
			TrieKMer tkm = new TrieKMer(seqs, 3);
			// we test with a 3-mer 'TTT'
			assertEquals(9, tkm
					.getCount(new KMer(alpha, new int[] { 3, 3, 3 })));
		} catch (IOException e) {
			fail("Sequence file not found: " + e.getMessage());
		}
	}

}
