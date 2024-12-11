package components;

import java.text.DecimalFormat;

/*
 * Contient toutes les informations liées à la présence d'un AA dans un spectre. "right_column" et "left_column" correspondent respectivement aux
 * colonnes formant le pic de droite et le pic de gauche de l'AA. "condition" est un attribut utilise dans la version 3peaks (voir la methode "runOneRealignmentThreePeaksVersion"
 * de "SpecGlobThread"). Il s'agit d'un entier sur 32 bits dont les 24 premiers bits sont utilises comme des booleens. Les 22 derniers bits correspondent aux
 * 22 AAs possible. En appelant 'X', l'AA correspondant a un objet "AAPosition", tout AA 'Y' precedant 'X' dans le spectre voit sont bit correspondant dans "condition" mit a 1.
 * Les autres sont mis a 0. Si plusieurs AAs precedent 'X', c'est qu'ils sont tous situe juste avant lui (ils ne forment pas une grande sequence mais plusieurs sequences de taille 2).
 * Les 2 premiers bits sont utilises lorsque 'X' n'est precede d'aucun AA. Soit dans le cas ou 'X' est situe au debut du spectre, auquel cas, le deuxieme bit est mit a 1 (les autres a 0),
 * soit parce qu'il n'est simplement pas precede par un autre AA, auquel cas le premier bit est mit a 1 (les autres a 0).
 * L'objectif de cela est de facilement pouvoir determiner si un AA d'une proteine doit etre considere s'il est present dans le spectre, en fonction de l'AA precedent dans la proteine
 * et de cette "condition". L'acide amine est concidere dans 2 cas :
 *  - le bit correspondant a l'AA precedent dans la protein est a 1 (impossible au debut de la proteine)
 *  - l'un des 2 premier bit est a 1
 *  Pour verifier ces conditions, il suffit de faire un 'et bit a bit' entre "condition" et un entier 32 bits avec les 2 premier bits a 1 et un 1 sur le bit correspondant
 *  a l'AA precedent dans la proteine (si il y en a un). Si le resultat de cette operation est 0, alors il ne faut pas conciderer l'AA, sinon, il doit l'etre.
 *  Les fait de distinguer le cas ou l'AA se trouve au debut de spectre et le cas ou il ne l'est pas et qu'il n'est pas precede par un AA dans le spectre est utile pour autre chose.
 */
public final class AAPosition {
	
	private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.00");

	private int right_column;
	private int left_column;
	private int left_peak;
	private int right_peak;
	private double right_mass;
	private double left_mass;
	private int condition;
	
	public AAPosition() {}

	public AAPosition(int right_column ,int left_column, int right_peak, int left_peak, double right_mass, double left_mass)
	{
		this.right_column = right_column;
		this.left_column = left_column;
		this.right_peak = right_peak;
		this.left_peak = left_peak;
		this.right_mass = right_mass;
		this.left_mass = left_mass;
	}
	
	public void overwrite(int right_column ,int left_column, int right_peak, int left_peak, double right_mass, double left_mass)
	{
		this.right_column = right_column;
		this.left_column = left_column;
		this.right_peak = right_peak;
		this.left_peak = left_peak;
		this.right_mass = right_mass;
		this.left_mass = left_mass;
	}
	
	public int getRightPeak() {return this.right_peak;}
	
	public int getLeftPeak() {return this.left_peak;}
	
	public int getRightColumn() {return this.right_column;}
	
	public int getLeftColumn() {return this.left_column;}
	
	public int getCondition() {return this.condition;}
	
	public void setCondition(int condition) 
	{
		if (condition == 0)
		{
			this.condition = 1;
		}	
		else
		{
			this.condition = condition;
		}
	}	
	
	/*
	 * Les colonnes attribues a l'AAPosition en premier lieu ont le meme numero que leur pic correspondant. Etant donne que certaines colonnes sont supprimees
	 * juste apres, il est necessaire de reattribuer un numero aux colonnes de maniere a ce que les numeros se suivent. Cette methode met a jour les numeros
	 * des colonnes de l'AAPosition.
	 */
	public void updateColumn(int[] column_correspondance)
	{
		this.right_column = column_correspondance[this.right_column] - 1;
		if (column_correspondance[this.left_column] != 0)
		{
			this.left_column = column_correspondance[this.left_column] - 2;
		}
		else
		{
			--this.left_column;
			while (column_correspondance[this.left_column] == 0)
			{
				--this.left_column;
			}
			this.left_column = column_correspondance[this.left_column] - 1;
		}	
	}
	
	@Override
	public String toString()
	{
		if ((this.condition & 1) == 1)
		{
			return "(" + DECIMAL_FORMAT.format(this.left_mass) + " ," + DECIMAL_FORMAT.format(this.right_mass) + ")";// " + Integer.toBinaryString(this.condition);
		}
		else
		{
			return "(" + DECIMAL_FORMAT.format(this.left_mass) + " ," + DECIMAL_FORMAT.format(this.right_mass) + ")*";// + Integer.toBinaryString(this.condition);
		}
	}
}
