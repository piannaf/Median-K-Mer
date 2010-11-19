package motif;

/**
 * A class for representing valid symbols (to be used for constructing
 * sequences, k-mers etc)
 */
public class Alphabet {

	private char[] SYMBOLS = { 'A', 'C', 'G', 'T' }; // default characters to

	// specify the a
	// sequence, note index
	// of array is used

	/**
	 * Constructs the default DNA alphabet consisting of symbols A, C, G and T
	 */
	public Alphabet() {
	}

	/**
	 * Constructs a custom-made alphabet
	 * 
	 * @param symbols
	 *            the symbols that the alphabet is composed of (the order is
	 *            important)
	 */
	public Alphabet(char[] symbols) {
		this.SYMBOLS = symbols;
	}

	/**
	 * Checks if the current and the specified alphabets are identical
	 * 
	 * @param other
	 *            the other alphabet
	 * @return true if the alphabets are identical, false otherwise
	 */
	public boolean equals(Alphabet other) {
		if (SYMBOLS.length != other.getSize())
			return false;
		for (int i = 0; i < SYMBOLS.length; i++)
			if (SYMBOLS[i] != other.getSymbol(i))
				return false;
		return true;
	}

	/**
	 * Retrieves the full list of symbols (order is important)
	 * 
	 * @return an array of all symbols in the alphabet
	 */
	public char[] getSymbols() {
		return SYMBOLS;
	}

	/**
	 * Retrieves the character used for symbol index
	 * 
	 * @param index
	 *            the symbol index
	 * @return the printable character
	 */
	public char getSymbol(int index) {
		if (index >= 0 && index < SYMBOLS.length)
			return SYMBOLS[index];
		else
			return '-';
	}

	/**
	 * Retrieves the size of the alphabet (the number of unique symbols)
	 * 
	 * @return the alphabet size
	 */
	public int getSize() {
		return SYMBOLS.length;
	}

	/**
	 * Utility method for translating a sequence represented by characters, to a
	 * sequence represented by indices
	 * 
	 * @param seq
	 *            the character sequence
	 * @return the index sequence
	 * @throws DNASequenceRuntimeException
	 *             if an invalid character is observed
	 */
	public int[] toIndex(char[] seq) {
		int[] indices = new int[seq.length];
		// do some error checking
		for (int i = 0; i < seq.length; i++) {
			// assume that the char is NOT valid
			int index = -1;
			for (int j = 0; j < SYMBOLS.length; j++) {
				if (seq[i] == SYMBOLS[j]) {
					index = j;
					break;
				}
			}
			if (index == -1)
				throw new DNASequenceRuntimeException("Invalid symbol \""
						+ seq[i] + "\" at position " + i);
			indices[i] = index;
		}
		return indices;
	}

	/**
	 * Utility method for translating a sequence represented by indices, to a
	 * sequence represented by characters
	 * 
	 * @param seq
	 *            the index sequence
	 * @return the character sequence
	 * @throws DNASequenceRuntimeException
	 *             if an invalid index is observed
	 */
	public char[] toChar(int[] seq) {
		char[] arr = new char[seq.length];
		for (int i = 0; i < seq.length; i++) {
			if (seq[i] >= SYMBOLS.length)
				throw new DNASequenceRuntimeException("Invalid index \""
						+ seq[i] + "\" at position " + i);
			if (seq[i] < 0)
				arr[i] = '-';
			else
				arr[i] = SYMBOLS[seq[i]];
		}
		return arr;
	}

}

/**
 * A class that holds information about alphabet-related exceptions
 */
class AlphabetRuntimeException extends RuntimeException {
	private static final long serialVersionUID = 1L;

	public AlphabetRuntimeException(String msg) {
		super(msg);
	}
}
