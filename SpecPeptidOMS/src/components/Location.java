package components;

import constantes.AminoAcid;

public class Location {
	
	private int[] protein_sequence;
	private int beginning;
	private int length;

	/*
	public Location(int[] protein_sequence, int beginning, int end)
	{
		this.protein_sequence = protein_sequence;
		this.beginning = Math.max(0, beginning);
		this.length = Math.min(protein_sequence.length - 1 ,end) - beginning + 1;
	}
	*/
	
	public void overwrite(int[] protein_sequence, int beginning, int end)
	{
		this.protein_sequence = protein_sequence;
		this.beginning = Math.max(0, beginning);
		this.length = Math.min(protein_sequence.length - 1 ,end) - this.beginning + 1;
	}
	
	public int getLength() {return this.length;}
	
	public int get(int position) {
		return this.protein_sequence[position + this.beginning];
	}
	
	@Override
	public String toString()
	{
		String ans = "";
		for (int i = 0; i < this.length; ++i)
		{
			ans = ans.concat(Character.toString(AminoAcid.getLetter(this.get(i))));
		}
		return ans;
	}
}
