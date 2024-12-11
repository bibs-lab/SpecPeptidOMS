package constantes;

import java.util.HashMap;

/*
 * @author Emile Benoist / Dominique Tessier
 */

/*
 * Chaque AA possede une masse totale et une masse de modification definie dans les parametres.
 * La masse total inclut la masse de base plus la masse de modification. Les 22 AAs sont stockes
 * dans un tableau statique et sont donc lies a un numero (leur position dans ce tableau).
 */
public final class AminoAcid
{
	private static double NTER_MASS = 1.007825032241 + Parameters.NTERMODIF();
	private static double CTER_MASS = 1.007825032241 + 15.99491461956 + Parameters.CTERMODIF();
	
	/*
	 * Le tableau statique contenant les 22 AAs
	 */
	private static final AminoAcid[] list = {
			new AminoAcid(57.021463721083 ,Parameters.GMODIF()), //G
			new AminoAcid(71.037113785565 ,Parameters.AMODIF()), //A
			new AminoAcid(87.032028405125 ,Parameters.SMODIF()), //S
			new AminoAcid(97.052763850047 ,Parameters.PMODIF()), //P
			new AminoAcid(99.068413914529 ,Parameters.VMODIF()), //V
			new AminoAcid(101.047678469607 ,Parameters.TMODIF()), //T
			new AminoAcid(103.009184785565 ,Parameters.CMODIF()), //C
			new AminoAcid(113.084063979011 ,Parameters.IMODIF()), //I
			new AminoAcid(113.084063979011 ,Parameters.LMODIF()), //L
			new AminoAcid(114.042927442166 ,Parameters.NMODIF()), //N
			new AminoAcid(115.026943024685 ,Parameters.DMODIF()), //D
			new AminoAcid(128.058577506648 ,Parameters.QMODIF()), //Q
			new AminoAcid(128.094963016052 ,Parameters.KMODIF()), //K
			new AminoAcid(129.042593089167 ,Parameters.EMODIF()), //E
			new AminoAcid(131.040484914529 ,Parameters.MMODIF()), //M
			new AminoAcid(137.058911859647 ,Parameters.HMODIF()), //H
			new AminoAcid(147.068413914529 ,Parameters.FMODIF()), //F
			new AminoAcid(156.101111025652 ,Parameters.RMODIF()), //R
			new AminoAcid(163.063328534089 ,Parameters.YMODIF()), //Y
			new AminoAcid(186.07931295157 ,Parameters.WMODIF()), //W
			new AminoAcid(168.964198469607 ,Parameters.UMODIF()), //U
			new AminoAcid(255.158291550141 ,Parameters.OMODIF()) //O
	};
	
	/*
	 * Une table de hashage statique pour retrouver rapidement le numero associe a un AA (via une lettre)
	 */
	@SuppressWarnings("serial")
	private static final HashMap<Character, Integer> numbers = new HashMap<Character, Integer>() {{
		put('G' ,0);
		put('A' ,1);
		put('S' ,2);
		put('P' ,3);
		put('V' ,4);
		put('T' ,5);
		put('C' ,6);
		put('I' ,7);
		put('L' ,8);
		put('N' ,9);
		put('D' ,10);
		put('Q' ,11);
		put('K' ,12);
		put('E' ,13);
		put('M' ,14);
		put('H' ,15);
		put('F' ,16);
		put('R' ,17);
		put('Y' ,18);
		put('W' ,19);
		put('U' ,20);
		put('O' ,21);
	}};
	
	/*
	 * Un tableau statique qui permet de retrouver rapidement la lettre d'un AA en fonction de sont numero
	 */
	private static final char[] letter = {'G' ,'A' ,'S' ,'P' ,'V' ,'T' ,'C' ,'I' ,'L' ,'N' ,'D' ,'Q' ,'K' ,'E' ,'M' ,'H' ,'F' ,'R' ,'Y' ,'W' ,'U' ,'O'};
	
	private static final double smallest_mass = 57.021463721083;
	
	private double mass;
	
	public AminoAcid(double mass ,double modification)
	{
		this.mass = mass + modification;
	}
	
	public double getMass() {return this.mass;}
	
	/*
	 * Retourne le nombre total d'AA
	 */
	public static int getCount() {return AminoAcid.list.length;}
	
	/*
	 * Retourne l'objet AA correspondant au numero en parametre
	 */
	public static AminoAcid get(int number) {return AminoAcid.list[number];}
	
	/*
	 * Retourne l'objet AA correspondant a la lettre en parametre
	 */
	public static AminoAcid get(char letter) {return AminoAcid.list[AminoAcid.numbers.get(letter)];}
	
	/*
	 * Retourne le numero de l'AA correspondant a la lettre en parametre
	 */
	public static int getNumber(char letter) {return AminoAcid.numbers.get(letter);}
	
	/*
	 * Retourne la lettre de l'AA correspondant au numero en parametre
	 */
	public static char getLetter(int number) {return AminoAcid.letter[number];}
	
	//public static double getHMass() {return 1.007825032241;}
	
	public static double getHMass() {return 1.007825032241;}
	
	public static double getHplusMass() {return 1.007276466879;}
	
	public static double getNTermMass() {return NTER_MASS;}
	
	public static double getCTermMass() {return CTER_MASS;}
	
	public static double SMALLESTMASS() {return smallest_mass;}
	
	public static boolean exists(char letter) {return numbers.containsKey(letter);}
	
	/*
	public static double getSequenceMass(String sequence)
	{
		double mass = 0.0;
		for (int i = 0; i < sequence.length(); ++i)
		{
			mass += AminoAcids.valueOf(Character.toString(sequence.charAt(i))).getMass();
		}
		return mass;
	}
	
	public static double getSequenceMass(AminoAcids[] sequence)
	{
		double mass = 0.0;
		for (AminoAcids aa : sequence)
		{
			mass += aa.getMass();
		}
		return mass;
	}
	*/
	
	public static void reset()
	{
		NTER_MASS = 1.007825032241 + Parameters.NTERMODIF();
		CTER_MASS = 1.007825032241 + 15.99491461956 + Parameters.CTERMODIF();
		list[0] = new AminoAcid(57.021463721083 ,Parameters.GMODIF());
		list[1] = new AminoAcid(71.037113785565 ,Parameters.AMODIF());
		list[2] = new AminoAcid(87.032028405125 ,Parameters.SMODIF());
		list[3] = new AminoAcid(97.052763850047 ,Parameters.PMODIF());
		list[4] = new AminoAcid(99.068413914529 ,Parameters.VMODIF());
		list[5] = new AminoAcid(101.047678469607 ,Parameters.TMODIF());
		list[6] = new AminoAcid(103.009184785565 ,Parameters.CMODIF());
		list[7] = new AminoAcid(113.084063979011 ,Parameters.IMODIF());
		list[8] = new AminoAcid(113.084063979011 ,Parameters.LMODIF());
		list[9] = new AminoAcid(114.042927442166 ,Parameters.NMODIF());
		list[10] = new AminoAcid(115.026943024685 ,Parameters.DMODIF());
		list[11] = new AminoAcid(128.058577506648 ,Parameters.QMODIF());
		list[12] = new AminoAcid(128.094963016052 ,Parameters.KMODIF());
		list[13] = new AminoAcid(129.042593089167 ,Parameters.EMODIF());
		list[14] = new AminoAcid(131.040484914529 ,Parameters.MMODIF());
		list[15] = new AminoAcid(137.058911859647 ,Parameters.HMODIF());
		list[16] = new AminoAcid(147.068413914529 ,Parameters.FMODIF());
		list[17] = new AminoAcid(156.101111025652 ,Parameters.RMODIF());
		list[18] = new AminoAcid(163.063328534089 ,Parameters.YMODIF());
		list[19] = new AminoAcid(186.07931295157 ,Parameters.WMODIF());
		list[20] = new AminoAcid(168.964198469607 ,Parameters.UMODIF());
		list[21] = new AminoAcid(255.158291550141 ,Parameters.OMODIF());		
	}
}