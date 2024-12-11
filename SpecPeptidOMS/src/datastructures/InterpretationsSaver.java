package datastructures;

import components.Alignment;
import components.LocationBackup;
import constantes.Parameters;
import loaders.Proteins;
import bioobjects.NativeSpectrum;

import java.io.PrintWriter;
import java.text.DecimalFormat;

/*
 * l'objet "ScenarioSaver" sert a stocker les 'x' meilleurs interpretations d'un spectre en dynamique, c'est à dire que ces 'x' meilleures peuvent changer en direct ('x' est a définir dans les parametres).
 * Un "ScenarioSaver" est principalement compose d'un tas binaire permettant d'effectuer en temps logarithmique en 'x', la modification dynamique des meilleures interpretations.
 * Chaque interpretation est stockee sous la forme d'un "ScenarioBackup"
 */
public final class InterpretationsSaver {
	
	private static int MIN_SCENARIO_SCORE = Parameters.MIN_SCENARIO_SCORE();
	private static double LOCATIONS_SCORE_THRESHOLD = Parameters.LOCATIONS_SCORE_THRESHOLD();
	private static int NB_LOCATIONS_SAVED = Parameters.NB_LOCATIONS_SAVED();
	private static int NB_INTERPRETATIONS_SAVED = Parameters.NB_INTERPRETATIONS_SAVED();
	
	private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.00");

	private LocationBackup[] locations_heap; // le tas binaire
	private LocationBackup current_heap_element;
	private int min_score_location; // le score de la moins bonne des interpretation stockee dans le tas
	private int max_score_location;
	private double locations_score_threshold;
	
	private int[] trees_position; // une liste permettant de connaitre la position d'une interpretation dans le tas en fonction de son numero d'arbre
	private int current_protein_number; // permet de retenir le numero de la protein sur laquelle est effectue l'alignement courant
	private NativeSpectrum current_spectrum; // permet de retenir le numero du spectre sur lequel est effectue l'alignement courant 
	
	private String[] simple_alignment_hits_modified;
	private String[] simple_alignment_sequences;
	private int[] simple_alignment_scores;
	private int[] simple_alignment_nb_common_peaks;
	
	private String[] non_aligned_mass_hits_modified;
	private String[] non_aligned_mass_sequences;
	private int[] non_aligned_mass_scores;
	private int[] non_aligned_mass_nb_common_peaks;
	private double[] non_aligned_masses;
	
	private String[] hit_modified_and_real_sequence;
	
	public InterpretationsSaver()
	{
		this.locations_heap = new LocationBackup[NB_LOCATIONS_SAVED]; // la taille du tas correspond au nombre maximum de localisations potentielles a retenir par spectre
		// le tas est initialise avec des "ScenarioBackup" par default
		for (int i = 0; i < this.locations_heap.length; ++i)
		{
			this.locations_heap[i] = new LocationBackup();
		}
		this.min_score_location = MIN_SCENARIO_SCORE;
		this.max_score_location = MIN_SCENARIO_SCORE;
		simple_alignment_hits_modified = new String[NB_INTERPRETATIONS_SAVED];
		simple_alignment_sequences = new String[NB_INTERPRETATIONS_SAVED];
		simple_alignment_scores = new int[NB_INTERPRETATIONS_SAVED];
		simple_alignment_nb_common_peaks = new int[NB_INTERPRETATIONS_SAVED];
		
		non_aligned_mass_hits_modified = new String[NB_INTERPRETATIONS_SAVED];
		non_aligned_mass_sequences = new String[NB_INTERPRETATIONS_SAVED];
		non_aligned_mass_scores = new int[NB_INTERPRETATIONS_SAVED];
		non_aligned_mass_nb_common_peaks = new int[NB_INTERPRETATIONS_SAVED];
		non_aligned_masses = new double[NB_INTERPRETATIONS_SAVED];
	}
	
	/*
	 * methode appelee lorsque l'on change de spectre (on pointe alors vers l'objet "trees_position" du nouveau spectre)
	 */
	public void updateSpectrum(int[] trees_position ,NativeSpectrum spectrum)
	{
		this.trees_position = trees_position;
		this.current_spectrum = spectrum;
	}
	
	/*
	 * methode appelee lorsque l'on change de proteine 
	 */
	public void updateProtein(int protein_number)
	{
		this.current_protein_number = protein_number;
	}
	
	/*
	 * cette methode tente d'ajouter une interpretation dans le tas. C'est possible si et seulement si son score est superieur au moins bon score des interpretations retenus jusqu'a maintenant.
	 * Si c'est le cas, on verifie qu'une autre interpretation 'x' du meme arbre n'est pas deja presente dans le tas. Si ce n'est pas le cas, alors la nouvelle interpretation prend simplement la place de la moins bonne.
	 * Si c'est le cas et si la nouvelle interpretation a un meilleur score que 'x', alors elle remplace 'x'
	 */
	public int tryAddLocation(int start_line ,int end_line ,int score ,int tree_id)
	{
		if (this.trees_position[tree_id] != -1 && this.locations_heap[this.trees_position[tree_id]].getProteinId() == this.current_protein_number)
		{
			if (this.locations_heap[this.trees_position[tree_id]].getScore() < score)
			{
				this.locations_heap[this.trees_position[tree_id]].update(start_line ,end_line ,score ,tree_id);
				this.max_score_location = Math.max(this.max_score_location, score);
				this.repositionLocation(this.trees_position[tree_id]);
			}
		}
		else
		{
			if (this.locations_heap[0].getProteinId() == this.current_protein_number)
			{
				this.trees_position[this.locations_heap[0].getTreeID()] = -1;
			}
			this.locations_heap[0].update(start_line, end_line, score, tree_id ,this.current_protein_number ,this.current_spectrum.getID());
			this.max_score_location = Math.max(this.max_score_location, score);
			this.repositionLocation(0);
		}
		return this.min_score_location;
	}
	
	/*
	 * methode interne liee au traitement du tas
	 */
	private void repositionLocation(int position)
	{
		//LocationBackup heap_element = this.locations_heap[position];
		this.current_heap_element = this.locations_heap[position];
		int left_son_position = (position << 1) + 1;
		int score_min_son;
		boolean direction_left = true;
		if (left_son_position < this.locations_heap.length)
		{
			score_min_son = this.locations_heap[left_son_position].getScore();
			if (left_son_position + 1 < this.locations_heap.length && score_min_son > this.locations_heap[left_son_position + 1].getScore())
			{
				score_min_son = this.locations_heap[left_son_position + 1].getScore();
				direction_left = false;
			}
		}
		else
		{
			score_min_son = this.current_heap_element.getScore();
		}
		while (this.current_heap_element.getScore() > score_min_son)
		{
			if (direction_left)
			{
				this.locations_heap[position] = this.locations_heap[left_son_position];
				if (this.locations_heap[position].getProteinId() == this.current_protein_number)
				{				
					this.trees_position[this.locations_heap[position].getTreeID()] = position;
				}
				position = left_son_position;
			}
			else
			{
				this.locations_heap[position] = this.locations_heap[left_son_position + 1];
				if (this.locations_heap[position].getProteinId() == this.current_protein_number)
				{
					this.trees_position[this.locations_heap[position].getTreeID()] = position;
				}
				position = left_son_position + 1;
			}
			left_son_position = (position << 1) + 1;
			if (left_son_position < this.locations_heap.length)
			{
				score_min_son = this.locations_heap[left_son_position].getScore();
				direction_left = true;
				if (left_son_position + 1 < this.locations_heap.length && score_min_son > this.locations_heap[left_son_position + 1].getScore())
				{
					score_min_son = this.locations_heap[left_son_position + 1].getScore();
					direction_left = false;
				}
			}
			else
			{
				score_min_son = this.current_heap_element.getScore();
			}
		}
		this.min_score_location = this.locations_heap[0].getScore();
		this.locations_heap[position] = this.current_heap_element;
		this.trees_position[this.current_heap_element.getTreeID()] = position;
	}
	
	public LocationBackup[] getLocations()
	{
		this.locations_score_threshold = this.max_score_location * LOCATIONS_SCORE_THRESHOLD;
		this.max_score_location = MIN_SCENARIO_SCORE;
		return this.locations_heap;
	}
	
	public double getLocationsScoreThreshold() {return this.locations_score_threshold;}
	
	public void saveSimpleAlignment(int interpretation_number, Alignment alignment)
	{
		if (alignment.getScore() > MIN_SCENARIO_SCORE)
		{
			//String[] res = alignment.computeHitModifiedAndRealSequence(spectrum);
			this.hit_modified_and_real_sequence = alignment.computeHitModifiedAndRealSequence();
			this.simple_alignment_hits_modified[interpretation_number] = this.hit_modified_and_real_sequence[0];
			this.simple_alignment_sequences[interpretation_number] = this.hit_modified_and_real_sequence[1];
			this.simple_alignment_scores[interpretation_number] = alignment.getScore();
			//this.simple_alignment_nb_common_peaks[interpretation_number] = alignment.nbPeaksInCommon(spectrum);
			alignment.computeNbPeaksInCommon();
			this.simple_alignment_nb_common_peaks[interpretation_number] = alignment.getNBPeaksInCommon();
		}
		else
		{
			this.simple_alignment_hits_modified[interpretation_number] = null;
		}
	}
	
	public void saveNonAlignedMass(int interpretation_number, Alignment alignment)
	{
		if (alignment.getScore() > MIN_SCENARIO_SCORE)
		{
			//String[] res = alignment.computeHitModifiedAndRealSequence(spectrum);
			this.hit_modified_and_real_sequence = alignment.computeHitModifiedAndRealSequence(); 
			this.non_aligned_mass_hits_modified[interpretation_number] = this.hit_modified_and_real_sequence[0];
			this.non_aligned_mass_sequences[interpretation_number] = this.hit_modified_and_real_sequence[1];
			this.non_aligned_mass_scores[interpretation_number] = alignment.getScore();
			//this.non_aligned_mass_nb_common_peaks[interpretation_number] = alignment.nbPeaksInCommon(spectrum);
			alignment.computeNbPeaksInCommon();
			this.non_aligned_mass_nb_common_peaks[interpretation_number] = alignment.getNBPeaksInCommon();
			this.non_aligned_masses[interpretation_number] = alignment.getTransformedSpectrum().getNonAlignedMass();
		}
		else
		{
			this.non_aligned_mass_hits_modified[interpretation_number] = null;
		}
	}
	
	/*
	 * methode appelee lorsque les resultats courant doivent etre sauvegardes dans le fichier de resultat
	 */
	
	public void save(PrintWriter output)
	{
		for (int interpretation_number = 0; interpretation_number < NB_INTERPRETATIONS_SAVED; ++interpretation_number)
		{
			if (this.simple_alignment_hits_modified[interpretation_number] != null)
			{
				output.write(String.format("%s;%d;%d", this.current_spectrum.getTitle(), this.current_spectrum.getScan(), this.current_spectrum.getID()));
				output.write(String.format(";%s;%s;%s;%d;%d", this.simple_alignment_sequences[interpretation_number], Proteins.foundPeptide(this.simple_alignment_sequences[interpretation_number]), this.simple_alignment_hits_modified[interpretation_number], this.simple_alignment_scores[interpretation_number], this.simple_alignment_nb_common_peaks[interpretation_number]));
				//if (this.non_aligned_mass_nb_common_peaks[interpretation_number] >= this.simple_alignment_nb_common_peaks[interpretation_number] + 3)
				if (this.non_aligned_mass_hits_modified[interpretation_number] != null)
				{
					output.write(String.format(";%s;%s;%s;%s;%d;%d\n", DECIMAL_FORMAT.format(this.non_aligned_masses[interpretation_number]), this.non_aligned_mass_sequences[interpretation_number], Proteins.foundPeptide(this.non_aligned_mass_sequences[interpretation_number]), this.non_aligned_mass_hits_modified[interpretation_number], this.non_aligned_mass_scores[interpretation_number], this.non_aligned_mass_nb_common_peaks[interpretation_number]));
				}
				else
				{
					output.write(";;;;;;;\n");
				}
			}
		}
	}
	
	public static void reset()
	{
		MIN_SCENARIO_SCORE = Parameters.MIN_SCENARIO_SCORE();
		LOCATIONS_SCORE_THRESHOLD = Parameters.LOCATIONS_SCORE_THRESHOLD();
		NB_LOCATIONS_SAVED = Parameters.NB_LOCATIONS_SAVED();
		NB_INTERPRETATIONS_SAVED = Parameters.NB_INTERPRETATIONS_SAVED();
	}
}
