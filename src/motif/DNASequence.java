package motif;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

/**
 * A class for representing biological sequence data
 */
public class DNASequence {

	private final Alphabet alpha;
	private final String name; // name of sequence
	private final int[] seq; // the sequence of symbols as represented by index

	/**
	 * Constructs a DNA sequence instance.
	 * 
	 * @param alpha
	 *            the alphabet from which valid symbols are drawn
	 * @param name
	 *            the name of the sequence
	 * @param seq
	 *            the character string representing the symbols of the sequence
	 *            (must be instances of alphabet).
	 * @see Alphabet
	 */
	public DNASequence(Alphabet alpha, String name, char[] seq) {
		this.alpha = alpha;
		this.name = name;
		// convert and check that the sequence is valid
		this.seq = alpha.toIndex(seq);
	}

	/**
	 * Retrieves the index of the symbol found at the specified position.
	 * 
	 * @param position
	 *            position of symbol 0..n-1 where n is the length of the
	 *            sequence
	 * @return the index of the symbol
	 * @throws DNASequenceRuntimeException
	 *             if an invalid position is given
	 */
	public int getSymbolIndex(int position) {
		if (position >= 0 && position < seq.length)
			return seq[position];
		else
			throw new DNASequenceRuntimeException(this,
					"Attempt to retrieve invalid index " + position + " in \"" + name + "\"");
	}

	/**
	 * Retrieves the indices of all the symbols in the sequence 0..n-1 where n
	 * is the length of the sequence
	 * 
	 * @return the indices
	 */
	public int[] getSymbolIndices() {
		return seq;
	}

	/**
	 * Retrieves the character representation of all the symbols in the sequence
	 * 0..n-1 where n is the length of the sequence
	 * 
	 * @return the printable characters
	 */
	public char[] getSymbolChars() {
		return alpha.toChar(getSymbolIndices());
	}

	/**
	 * Retrieve the printable character of the symbol found at the specified
	 * position.
	 * 
	 * @param position
	 *            position of symbol 0..n-1 where n is the length of the
	 *            sequence
	 * @return the character of the symbol
	 * @throws DNASequenceRuntimeException
	 *             if an invalid position is given
	 */
	public char getSymbolChar(int position) {
		int symindex = getSymbolIndex(position);
		if (symindex >= 0 && symindex < alpha.getSize())
			return alpha.getSymbol(symindex);
		else
			throw new DNASequenceRuntimeException("Invalid symbol index \""
					+ symindex + "\" at position " + position);
	}

	/**
	 * Retrieve the length of the DNA sequence.
	 * 
	 * @return the length (number of symbols)
	 */
	public int getLength() {
		if (seq != null)
			return seq.length;
		else
			throw new DNASequenceRuntimeException(this,
					"Invalid symbol sequence in \"" + name + "\"");
	}

	/**
	 * Retrieves the alphabet that is used for this sequence
	 * 
	 * @return the alphabet
	 */
	public Alphabet getAlphabet() {
		return alpha;
	}

	/**
	 * Printable representation of sequence
	 */
	public String toString() {
		return name + " (" + seq.length + ")";
	}

	/**
	 * Reads DNA sequences from a file on the FASTA standard format
	 * 
	 * @param filename
	 *            the name of the file
	 * @return an array of instance of {@link #DNASequence}
	 * @throws IOException
	 *             if the file operation fails
	 */
	public static DNASequence[] readFile(Alphabet alpha, String filename)
			throws IOException {
		List<DNASequence> seqs = new ArrayList<DNASequence>();
		BufferedReader br = new BufferedReader(new FileReader(
				new File(filename)));
		// buffer variables to hold recently read data
		String name = null;
		StringBuffer buf = null;

		int row = 0;
		String line = br.readLine();
		while (line != null) {
			row++;
			line = line.trim(); // remove any spaces, tabs etc at the ends
			if (line.startsWith(">")) {
				if (buf != null) // there is data in the buffer, we need to
				// store it before processing the new entry
				{
					try {
						seqs.add(new DNASequence(alpha, name, buf.toString()
								.toCharArray()));
					} catch (DNASequenceRuntimeException e) {
						System.err.println("Ignored " + name + ": "
								+ e.getMessage());
					}
					buf = null;
					name = null;
				}
				try {
					StringTokenizer stok = new StringTokenizer(line, " \t");
					name = stok.nextToken().substring(1);
				} catch (NoSuchElementException e) {
					throw new RuntimeException("Invalid format in file "
							+ filename + " at row " + row);
				}
				buf = new StringBuffer();
			} else {
				if (buf != null) {
					buf.append(line);
				}
			}
			line = br.readLine();
		}
		if (buf != null) // there is data in the buffer, we need to store it
						 // before processing the new entry
		{
			try {
				seqs.add(new DNASequence(alpha, name, buf.toString()
						.toCharArray()));
			} catch (DNASequenceRuntimeException e) {
				System.err.println("Ignored " + name + ": " + e.getMessage());
			}
			buf = null;
			name = null;
		}
		DNASequence[] all = new DNASequence[seqs.size()];
		seqs.toArray(all);
		return all;
	}

	/**
	 * Example application that simply loads a FASTA file and prints out the
	 * sequences in it.
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		for (int i = 0; i < args.length; i++) {
			try {
				DNASequence[] seqs = DNASequence.readFile(new Alphabet(),
						args[i]);
				System.out.println("Read " + seqs.length + " sequences from "
						+ args[i]);
				for (DNASequence s : seqs)
					System.out.println(s);
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
	}
}

/**
 * A class that holds information about sequence-related runtime exceptions
 */
class DNASequenceRuntimeException extends RuntimeException {
	private static final long serialVersionUID = 1L;
	public DNASequence s = null;

	public DNASequenceRuntimeException(DNASequence s, String msg) {
		super(msg);
		this.s = s;
	}

	public DNASequenceRuntimeException(String msg) {
		super(msg);
	}
}
