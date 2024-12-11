package bioobjects;

import constantes.AminoAcid;

import java.util.ArrayList;

/**
 * This class represents a protein with its sequence of amino acids, its name and a description
 * @author 	BENOIST Emile, TESSIER Dominique
 *
 */
public final class Protein {

	private String real_sequence;
	private int[] sequence;
	private String name;
	private String protein_description;
	
	/*
	 * Creation d'une proteine à partir d'une "sequence" d'AAs (lettre en majuscule). L'argument "reference" contient l'ensemble 
	 * des informations sur la proteine qui sont presentes dans le fichier et qui ne sont pas la sequence.
	 */
	/**
	 * Constructor
	 * @param sequence: the amino acids sequence
	 * @param name: its name
	 * @param reference: its reference
	 */
	public Protein(String sequence, String name ,String reference) {
		
		this.name = name;
		protein_description = reference;
		this.real_sequence = sequence;
		
		/*
		 * An temporary integer version of the amino acids sequence
		 */
		ArrayList<Integer> sequence_temp = new ArrayList<Integer>(sequence.length());
		
		char c;
		for (int i = sequence.length() - 1; i >= 0; --i)
		{
			c = sequence.charAt(i);
			if (AminoAcid.exists(c))
			{
				sequence_temp.add(AminoAcid.getNumber(c));
			}
		}
		
		/*
		 * The temporary sequence is copied in the non temporary one with a fixed length
		 */
		this.sequence = new int[sequence_temp.size()];
		for (int i = 0; i < this.sequence.length; ++i)
		{
			this.sequence[i] = sequence_temp.get(i);
		}
	}
	
	public String getName() {return this.name;}
	
	/**
	 * Getter of the protein length
	 * @return protein sequence length
	 */
	public int getLength() {return this.sequence.length;}
	
	/**
	 * Return the reference of the integer version of the protein sequence
	 * @return protein sequence with integers
	 */
	public int[] getSequence() {return this.sequence;}
	
	/**
	 * Return  the string version of the protein sequence
	 * @return protein sequence with characters
	 */
	public final String getRealSequence() {return this.real_sequence;}
	
	/**
	 * Return the integer version of the subsequence starting (resp. ending) at position 'start' (resp. 'end'), both include
	 * @param start: the beginning of the subsequence
	 * @param end: the end of the subsequence
	 * @return the integer subsequence
	 */
	public int[] getSubSequence(int start, int end)
	{
		start = Math.max(0, start);
		end = Math.min(this.sequence.length - 1 ,end);
		int[] sub_sequence = new int[end - start + 1];
		for (int i = start; i <= end; ++i)
		{
			sub_sequence[i - start] = this.sequence[i];
		}
		return sub_sequence;
	}
	
	/**
	 * Return the amino acid at the position in parameter
	 * @param position: the position in the sequence of the returned amino acid
	 * @return the amino acid at the given position
	 */
	public AminoAcid getAminoAcid(int position) {return AminoAcid.get(this.sequence[position]);}
	
	/**
	 * Return the integer version of the amino acid at the position in parameter
	 * @param position: the position in the sequence of the returned amino acid
	 * @return the integer version of the amino acid at the given position
	 */
	public final int getAminoAcidNumber(int position) {return this.sequence[position];}	
	
	/**
	 * Return the letter corresponding to the amino acid at the position in parameter
	 * @param position: the position in the sequence of the amino acid
	 * @return the letter corresponding to the amino acid at the given position
	 */
	public final String getLetter(int position) {return Character.toString(AminoAcid.getLetter(this.sequence[position]));}
	
	/**
	 * 
	 * Return the string version of the subsequence starting (resp. ending) at position 'start' (resp. 'end'), both include
	 * @param start: the beginning of the subsequence
	 * @param end: the end of the subsequence
	 * @return the string subsequence
	 */
	public final String subSequence(int start ,int end)
	{
		String sub_sequence = "";
		for (int i = start; i <= end; ++i)
		{
			sub_sequence = sub_sequence.concat(Character.toString(AminoAcid.getLetter(this.sequence[i])));
		}
		return sub_sequence;
	}
	
	/**
	 * Print the protein sequence with letters
	 */
	public void print()
	{
		for (int i = 0; i < this.sequence.length; ++i)
		{
			System.out.println(AminoAcid.getLetter(this.sequence[i]));
		}
		System.out.println();
	}
	
	/**
	 * Print the protein description
	 */
	public void printProteinDescription()
	{
		 System.out.println(this.protein_description);
	}
}
