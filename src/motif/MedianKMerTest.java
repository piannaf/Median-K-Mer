/**
 * 
 */
package motif;

import java.io.IOException;

import org.junit.Test;
import junit.framework.TestCase;

/**
 */
public class MedianKMerTest extends TestCase {

	public MedianKMerTest(String name) {
		super(name);
	}

	public void testGetDistance1() {
		Alphabet alpha = new Alphabet();
		KMer km = new KMer(alpha, new int[] { 1, 2, 0 }); // "CGA"
		assertEquals(1, MedianKMer.getDistance(new DNASequence(alpha, "s1",
				new char[] { 'T', 'C', 'G', 'G', 'A', 'C' }), km));
	}

	public void testGetDistance2() {
		Alphabet alpha = new Alphabet();
		KMer km = new KMer(alpha, new int[] { 1, 2, 0, 0 }); // "CGAA"
		assertEquals(4, MedianKMer.getDistance(new DNASequence(alpha, "s1",
				new char[] { 'A', 'A', 'T', 'T', 'C', 'C' }), km));
	}

	public void testGetDistance3() {
		Alphabet alpha = new Alphabet();
		KMer km = new KMer(alpha, new int[] { 0, 0 }); // "AA"
		assertEquals(0, MedianKMer.getDistance(new DNASequence(alpha, "s1",
				new char[] { 'T', 'A', 'T', 'G', 'A', 'A' }), km));
	}

	public void testGetDistance4() {
		Alphabet alpha = new Alphabet();
		KMer km = new KMer(alpha, new int[] { 0, 0 }); // "AA"
		int[] d = MedianKMer.getDistances(new DNASequence(alpha, "s1",
				new char[] { 'T', 'A', 'T', 'G', 'A', 'A' }), km);
		int pos = MedianKMer.getMinPosition(d);
		assertEquals(4, pos);
	}

	public void testFindMedianKMer1() {
		Alphabet alpha = new Alphabet();
		DNASequence[] seqs = new DNASequence[] {
				new DNASequence(alpha, "s1", 
						new char[] { 'A', 'C', 'G', 'G', 'A', 'C' }),
				new DNASequence(alpha, "s2", 
						new char[] { 'A', 'G', 'G', 'A', 'C', 'G' }),
				new DNASequence(alpha, "s3", 
						new char[] { 'A', 'A', 'C', 'G', 'C', 'C' }) };

		MedianKMer m = new MedianKMer(seqs);
		MedianKMer.Distance d = m.findMedianKMer(3);
		assertEquals("ACG:0", d.toString());
	}

	public void testFindMedianKMer2() {
		Alphabet alpha = new Alphabet();
		DNASequence[] seqs = new DNASequence[] {
				new DNASequence(alpha, "s1", 
						new char[] { 'T', 'C', 'G', 'G', 'A', 'C' }),
				new DNASequence(alpha, "s2", 
						new char[] { 'A', 'G', 'G', 'T', 'T', 'G' }),
				new DNASequence(alpha, "s3", 
						new char[] { 'T', 'A', 'A', 'G', 'G', 'C' }) };

		MedianKMer m = new MedianKMer(seqs);
		MedianKMer.Distance d = m.findMedianKMer(3);
		assertEquals("AGG:1", d.toString());
	}

	public void testFindMedianKMer3() {
		Alphabet alpha = new Alphabet();
		DNASequence[] seqs = new DNASequence[] {
				new DNASequence(alpha, "s1", 
						new char[] { 'T', 'C', 'G', 'G', 'A', 'C' }),
				new DNASequence(alpha, "s2", 
						new char[] { 'A', 'C', 'G', 'T', 'T', 'G' }),
				new DNASequence(alpha, "s3", 
						new char[] { 'T', 'A', 'A', 'G', 'T', 'C' }) };

		MedianKMer m = new MedianKMer(seqs);
		MedianKMer.Distance d = m.findMedianKMer(6);
		assertEquals("TCGGTC:6", d.toString());
	}
	
	public void testFindMedianKMer4() {
		Alphabet alpha = new Alphabet();
		try {
			DNASequence[] seqs = DNASequence.readFile(alpha,
					"data/malT_5.fasta");
			MedianKMer m = new MedianKMer(seqs);
			MedianKMer.Distance d = m.findMedianKMer(5);
			assertEquals(1, d.actual);
		} catch (IOException e) {
			fail("Sequence file not found: " + e.getMessage());
		}
	}

	// The following test takes about 2 seconds on my 3yo macbook pro when
	// improvements in Problem 4 have been completed.
	@Test(timeout = 3000)
	public void testFindMedianKMer5() {
		Alphabet alpha = new Alphabet();
		try {
			DNASequence[] seqs = DNASequence.readFile(alpha,
					"data/malT_5.fasta");
			MedianKMer m = new MedianKMer(seqs);
			MedianKMer.Distance d = m.findMedianKMer(9);
			assertEquals(14, d.actual);
		} catch (IOException e) {
			fail("Sequence file not found: " + e.getMessage());
		}
	}
}
