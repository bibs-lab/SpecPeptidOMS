package components;

import constantes.Parameters;

/*
 * Cette objet contient une liste de scenario, formant ainsi, une interpretation complete du spectre.
 * Lorsque un alignement entre une proteine est un spectre est effectue, une matrice de scenario est mise a jour.
 * En passant a la proteine suivante (ou au spectre suivant), il est necessaire d'effectuer une sauvegarde des bonnes interpretations
 * (liste de scenarios consecutif) de la matrice pour que ceux la ne soit pas ecrase par l'alignement suivant.
 * L'objet "ScenarioBackup" correspond a la sauvegarde de l'une de ces interpretation. Lorsque qu'une interpretation interessante est trouvee,
 * l'une des methodes "update" est appelee pour mettre a jour un "pointeur" vers le scenario de la matrice qui termine dans l'interpretation.
 * Lorsque que l'alignement est termine, une copie de l'interpretation pointee est effectuee ici pour la sauvegarder. 
 */
public final class LocationBackup {
	
	private static int MIN_SCENARIO_SCORE = Parameters.MIN_SCENARIO_SCORE();
	
	private int start_line;
	private int end_line;
	private int score;
	private int tree_id;
	private int protein_id;
	private int spectrum_id;
	
	/*
	 * Un "ScenarioBackup" est initialise avec une liste de scenarios vides. Le nombre de ces scenarios est fixe au nombre de colonne que contient
	 * la matrice de score car une interpretation ne peut pas contenir plus de scenarios qu'il n'y a de colonne (les cas 'not found' ne sont pas compte
	 * comme des scenarios).
	 */
	public LocationBackup()
	{	
		this.score = MIN_SCENARIO_SCORE;
		this.protein_id = -1;
		this.spectrum_id = -1;
	}
	
	public void update(int start_line ,int end_line ,int score ,int tree_id)
	{
		this.start_line = start_line;
		this.end_line = end_line;
		this.score = score;
		this.tree_id = tree_id;
	}
	
	public void update(int start_line ,int end_line ,int score ,int tree_id ,int protein_number, int spectrum_id)
	{
		this.start_line = start_line;
		this.end_line = end_line;
		this.score = score;
		this.tree_id = tree_id;
		this.protein_id = protein_number;
		this.spectrum_id = spectrum_id;
	}
	
	public int getStartLine() {return this.start_line;}
	
	public int getEndLine() {return this.end_line;}
		
	public int getScore() {return this.score;}
	
	public int getTreeID() {return this.tree_id;}
	
	public int getProteinId() {return this.protein_id;}
	
	public int getSpectrumId() {return this.spectrum_id;}
	
	@Override
	public String toString()
	{
		return String.format("[%d ,%d], score : %d ,protein %d, spectrum %d ,tree %d", this.start_line ,this.end_line ,this.score ,this.protein_id, this.spectrum_id ,this.tree_id);
	}
	
	public static void reset()
	{
		MIN_SCENARIO_SCORE = Parameters.MIN_SCENARIO_SCORE();
	}
}
