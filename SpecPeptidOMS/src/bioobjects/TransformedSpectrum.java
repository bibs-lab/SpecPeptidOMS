package bioobjects;

import java.util.ArrayList;
import java.util.Arrays;

import components.AAPosition;
import components.ComplementaryPeak;
import constantes.AminoAcid;
import constantes.Parameters;
import loaders.Spectra;

/**
 * This class represent a spectrum after its pre-treatment.
 * Some of the attributes of a TransformedSpectrum can be stored in static structure because they are treated sequentially.
 * @author 	BENOIST Emile, TESSIER Dominique
 *
 */
public class TransformedSpectrum {
	
	private static double ACCURACY = Parameters.ACCURACY();
	private static int REAL_TIME_SAVE = Parameters.REAL_TIME_SAVE();
	
	private int[] peak_mask;
	private double non_aligned_mass;
	
	private int id;
	private double precursor_mass;
	private int middle_peak;
	private AAPosition[][] amino_acids_positions;
	private int[] amino_acids_positions_number; 
	private boolean[] double_column;
	private double[] column_mass;
	private double last_mass;
	private int nb_complementary_peaks;
	private int column_count;
	
	
	//private static int[][] peak_masks;
	private static ComplementaryPeak[][] all_peaks;
	private static boolean[][] checked;
	private static ArrayList<AAPosition>[] amino_acids_positions_temp;
	private static int[][] ending;
	private static boolean[][] selected_peaks;
	private static int[][] column_correspondance;
	
	@SuppressWarnings("unchecked")
	/**
	 * This method initialize all the static attributes by giving then a fixed size that can be known after loading the spectrum file.
	 */
	public static void initialize()
	{
		int max_native_peak_count = Spectra.getMaxNativePeakCount();
		//TransformedSpectrum.peak_masks = new int[REAL_TIME_SAVE][max_native_peak_count];
		TransformedSpectrum.all_peaks = new ComplementaryPeak[REAL_TIME_SAVE][];
		for (int i = 0; i < TransformedSpectrum.all_peaks.length; ++i)
		{
			TransformedSpectrum.all_peaks[i] = new ComplementaryPeak[(max_native_peak_count << 1) + 2];
			for (int j = 0; j < TransformedSpectrum.all_peaks[i].length; ++j)
			{
				TransformedSpectrum.all_peaks[i][j] = new ComplementaryPeak();
			}
		}
		TransformedSpectrum.checked = new boolean[REAL_TIME_SAVE][max_native_peak_count];
		TransformedSpectrum.amino_acids_positions_temp = new ArrayList[AminoAcid.getCount()];
		for (int i = 0; i < TransformedSpectrum.amino_acids_positions_temp.length; ++i)
		{
			TransformedSpectrum.amino_acids_positions_temp[i] = new ArrayList<AAPosition>();
		}
		TransformedSpectrum.ending = new int[REAL_TIME_SAVE][(Spectra.getMaxNativePeakCount() << 1) + 2];
		TransformedSpectrum.selected_peaks = new boolean[REAL_TIME_SAVE][(Spectra.getMaxNativePeakCount() << 1) + 2];
		TransformedSpectrum.column_correspondance = new int[REAL_TIME_SAVE][(Spectra.getMaxNativePeakCount() << 1) + 2];
	}
	
	/**
	 * The first constructor that build the first transformed version of a native spectrum whose the attributes are in parameter.
	 * It correspond to the pre-treatment of a spectrum
	 * @param peaks: the peaks list of the native spectrum
	 * @param precursor_mass_charged: the precursor mass/charge ratio
	 * @param id: the native spectrum id
	 * @param charge: the charge of the precursor
	 */
	public TransformedSpectrum(double[] peaks, double precursor_mass_charged, int id, int charge)
	{
		this.non_aligned_mass = 0.0;
		this.id = id;
		this.precursor_mass = (precursor_mass_charged * charge) - (charge * AminoAcid.getHplusMass());
		this.last_mass = this.complementary(AminoAcid.getHplusMass());
		double first_mass = AminoAcid.getCTermMass() + AminoAcid.getNTermMass() + AminoAcid.getHplusMass();
		ComplementaryPeak peak_temp;
		
		Arrays.fill(TransformedSpectrum.checked[0], 0, peaks.length, false);
		int i;
		int j = peaks.length - 1;
		TransformedSpectrum.all_peaks[0][0].overwrite(first_mass, false ,0);
		TransformedSpectrum.all_peaks[0][1].overwrite(this.last_mass, false ,peaks.length + 1);
		this.nb_complementary_peaks = 2;
		int current_add_count;
		double native_mass;
		double complementary_mass;
		for (i = 0; i < peaks.length; ++i)
		{
			if (!TransformedSpectrum.checked[0][i])
			{
				TransformedSpectrum.checked[0][i] = true;
				native_mass = peaks[i];
				complementary_mass = this.complementary(native_mass);
				if (native_mass < complementary_mass)
				{
					this.middle_peak = i + 1;
				}
				while (j > 0 && Math.abs(peaks[j] - complementary_mass) > Math.abs(peaks[j - 1] - complementary_mass))
				{
					--j;
				}
				//0 -> les 2 / 1 -> le natif / -1 -> le complémentaire				
				if (Math.abs(peaks[j] - complementary_mass) > ACCURACY)
				{
					TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks].overwrite(complementary_mass ,false ,i + 1);
					TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks + 1].overwrite(native_mass ,false ,i + 1);
					current_add_count = 2;
				}
				//else if (checked[j])
				else if (TransformedSpectrum.checked[0][j])
				{
					if (j == i)
					{
						TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks].overwrite(native_mass ,false ,i + 1);
						current_add_count = 1;
					}
					else
					{
						TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks].overwrite(complementary_mass ,false ,i + 1);
						TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks + 1].overwrite(native_mass ,false ,i + 1);
						current_add_count = 2;
					}
				}
				else
				{
					TransformedSpectrum.checked[0][j] = true;
					TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks].overwrite(complementary_mass ,true ,i + 1);
					TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks + 1].overwrite(native_mass ,true ,i + 1);
					current_add_count = 2;
				}
				if (current_add_count == 1)
				{
					if (TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks].get_mass() <= first_mass || this.last_mass <= TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks].get_mass())
					{
						current_add_count = 0;
					}
				}
				else
				{
					if (TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks].get_mass() <= first_mass || this.last_mass <= TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks].get_mass())
					{
						peak_temp = TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks];
						TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks] = TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks + 1];
						TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks + 1] = peak_temp;
						--current_add_count;
					}
					if (TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks + current_add_count - 1].get_mass() <= first_mass || this.last_mass <= TransformedSpectrum.all_peaks[0][this.nb_complementary_peaks + current_add_count - 1].get_mass())
					{
						--current_add_count;
					}
				}
				this.nb_complementary_peaks += current_add_count;
			}
		}
		Arrays.sort(TransformedSpectrum.all_peaks[0] ,0 ,this.nb_complementary_peaks);
		Arrays.fill(TransformedSpectrum.ending[0], 0, this.nb_complementary_peaks, 0);
		Arrays.fill(TransformedSpectrum.selected_peaks[0], 0, this.nb_complementary_peaks, false);
		int nb_aa_founded = 0;
		this.amino_acids_positions_number = new int[AminoAcid.getCount()];
		int k;
		for (k = 0; k < TransformedSpectrum.amino_acids_positions_temp.length; ++k)
		{
			TransformedSpectrum.amino_acids_positions_temp[k].clear();
			
		}
		this.column_count = 1;
		double difference;
		for (j = 1; j < this.nb_complementary_peaks; ++j) // Cette partie peut sans doute être améliorée !!
		{
			for (i = 0; i < j; ++i)
			{
				for (k = 0; k < AminoAcid.getCount(); ++k)
				{
					difference = (TransformedSpectrum.all_peaks[0][j].get_mass() - TransformedSpectrum.all_peaks[0][i].get_mass()) - AminoAcid.get(k).getMass();
					if (Math.abs(difference) <= Parameters.ACCURACY() && (this.amino_acids_positions_number[k] == 0 || TransformedSpectrum.amino_acids_positions_temp[k].get(this.amino_acids_positions_number[k] - 1).getRightColumn() < j))
					{
						++nb_aa_founded;
						TransformedSpectrum.ending[0][j] |= 1 << (k + 2);
						TransformedSpectrum.amino_acids_positions_temp[k].add(new AAPosition(j, i, TransformedSpectrum.all_peaks[0][j].getNativePeak(), TransformedSpectrum.all_peaks[0][i].getNativePeak(), TransformedSpectrum.all_peaks[0][j].get_mass(), TransformedSpectrum.all_peaks[0][i].get_mass()));
						++this.amino_acids_positions_number[k];
						if(!TransformedSpectrum.selected_peaks[0][j])
						{
							TransformedSpectrum.selected_peaks[0][j] = true;
							++this.column_count;
						}
					}
				}
			}
		}
		

		if (nb_aa_founded >= Parameters.NB_MIN_AA_FOUND()) // Le spectre est cree si et seulement si suffisemment d'AAs ont ete trouve dedans
		{
			// Le contenu de "amino_acids_positions_temp" est copie dans la matrice "this.amino_acids_positions" pour plus d'efficacite
			this.amino_acids_positions = new AAPosition[AminoAcid.getCount()][];
			for (k = 0; k < AminoAcid.getCount(); ++k)
			{
				this.amino_acids_positions[k] = new AAPosition[this.amino_acids_positions_number[k]];
				for (i = 0; i < this.amino_acids_positions_number[k]; ++i)
				{
					this.amino_acids_positions[k][i] = TransformedSpectrum.amino_acids_positions_temp[k].get(i);
					// L'appel a la methode "SetCondition" est lie a la version 3peaks (voir la documentation de la classe AAPosition)
					if (this.amino_acids_positions[k][i].getLeftColumn() == 0)
					{
						this.amino_acids_positions[k][i].setCondition(2);
					}
					else
					{
						this.amino_acids_positions[k][i].setCondition(TransformedSpectrum.ending[0][this.amino_acids_positions[k][i].getLeftColumn()]);
					}
				}
			}
			
			// Seul les pics utiles (a droite d'un AA dans le spectre) generent une colonne dans la matrice de score
			Arrays.fill(TransformedSpectrum.column_correspondance[0], 0);
			TransformedSpectrum.column_correspondance[0][0] = 1;
			this.double_column = new boolean[this.column_count];
			// "this.column_mass" sert pour le calcul des shifts
			this.column_mass = new double[this.column_count];
			this.column_mass[0] = TransformedSpectrum.all_peaks[0][0].get_mass();
			int column_index = 1;
			i = 0;
			while (column_index <= this.column_count - 1)
			{
				if (TransformedSpectrum.selected_peaks[0][i])
				{
					if (TransformedSpectrum.all_peaks[0][i].isDouble())
					{
						this.double_column[column_index] = true;
					}
					this.column_mass[column_index] = TransformedSpectrum.all_peaks[0][i].get_mass();
					TransformedSpectrum.column_correspondance[0][i] = column_index + 1;
					++column_index;
				}
				++i;
			}
			
			for (k = 0; k < AminoAcid.getCount(); ++k)
			{
				for (i = 0; i < this.amino_acids_positions_number[k]; ++i)
				{
					this.amino_acids_positions[k][i].updateColumn(TransformedSpectrum.column_correspondance[0]);
				}
			}
		}
	}
	
	/**
	 * The second constructor that only initializes the size of certain structures. The spectrum obtained will then be used to contain
	 * a TransformedSpectrum with some modification (possibly ignored peaks and possible non-aligned mass).
	 */
	public TransformedSpectrum()
	{
		this.peak_mask = new int[Spectra.getMaxNativePeakCount()];
		this.amino_acids_positions = new AAPosition[AminoAcid.getCount()][];
		for (int i = 0; i < this.amino_acids_positions.length; ++i)
		{
			this.amino_acids_positions[i] = new AAPosition[(Spectra.getMaxNativePeakCount() << 1) + 1];
			for (int j = 0; j < this.amino_acids_positions[i].length; ++j)
			{
				this.amino_acids_positions[i][j] = new AAPosition();
			}
		}
		this.amino_acids_positions_number = new int[AminoAcid.getCount()];
		this.double_column = new boolean[(Spectra.getMaxNativePeakCount() << 1) + 2];
		this.column_mass = new double[(Spectra.getMaxNativePeakCount() << 1) + 2];
	}
	
	
	/**
	 * This method overwrite a TransformedSpectrum with a new one. Each parameters corresponds to an attribute of the new spectrum.
	 * The 'id' must correspond to the id of the native spectrum. The parameters 'peak_mask' is used when an alignment rely on a peak and its complementary.
	 * More precisely, if a series of peak i(1), i(2), ..., i(n) and their complementary are all used in the same alignment, the algorithm generate all possible
	 * TransformedSpectrum with only i(1) or only its complementary, same thing for i(2), ..., i(n).
	 * A 'peak_mask' associate an integer to each peak : 0 if the peak and its complementary are needed, 1 if only the original peak is needed, -1 if only the complementary is needed.
	 * The parameter 'non_aligned_mass' correspond to a supposed error in the precursor mass. The new TransformedSpectrum is built taken into account this supposed error. 
	 * @param peaks: the list of peaks of the native spectrum
	 * @param precursor_mass_charged: the mass of the charged precursor
	 * @param id: the id of the native spectrum
	 * @param charge: the charge of the precursor
	 * @param peak_mask
	 * @param non_aligned_mass
	 */
	public void overwrite(double[] peaks, double precursor_mass_charged, int id, int charge, int[] peak_mask, double non_aligned_mass)
	{
		System.arraycopy(peak_mask, 0, this.peak_mask, 0, peak_mask.length);
		this.non_aligned_mass = non_aligned_mass;
		
		this.id = id;
		this.precursor_mass = ((precursor_mass_charged) * charge) - (charge * AminoAcid.getHplusMass()) - non_aligned_mass;
		this.last_mass = this.complementary(AminoAcid.getHplusMass());
		double first_mass = AminoAcid.getCTermMass() + AminoAcid.getNTermMass() + AminoAcid.getHplusMass();
		ComplementaryPeak peak_temp;
		
		int current_modulo_id = id % REAL_TIME_SAVE;
		Arrays.fill(TransformedSpectrum.checked[current_modulo_id], false);
		int i;
		int j = peaks.length - 1;
		while (j >= 0 && peaks[j] >= this.last_mass)
		{
			TransformedSpectrum.checked[current_modulo_id][j] = true;
			--j;
		}
		
		TransformedSpectrum.all_peaks[current_modulo_id][0].overwrite(AminoAcid.getCTermMass() + AminoAcid.getNTermMass() + AminoAcid.getHplusMass(), false ,0);
		TransformedSpectrum.all_peaks[current_modulo_id][1].overwrite(this.last_mass, false ,peaks.length + 1);
		this.nb_complementary_peaks = 2;
		int current_add_count;
		double native_mass;
		double complementary_mass;
		for (i = 0; i < peaks.length; ++i)
		{
			if (!TransformedSpectrum.checked[current_modulo_id][i])
			{
				current_add_count = 0;
				TransformedSpectrum.checked[current_modulo_id][i] = true;
				native_mass = peaks[i];
				complementary_mass = this.complementary(native_mass);
				if (native_mass < complementary_mass)
				{
					this.middle_peak = i + 1;
				}
				while (j > 0 && Math.abs(peaks[j] - complementary_mass) > Math.abs(peaks[j - 1] - complementary_mass))
				{
					--j;
				}
				//0 -> les 2 / 1 -> le natif / -1 -> le complémentaire				
				if (Math.abs(peaks[j] - complementary_mass) > ACCURACY)
				{
					if (this.peak_mask[i] <= 0)
					{
						TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks].overwrite(complementary_mass ,false ,i + 1);
						++current_add_count;
					}
					if (this.peak_mask[i] >= 0)
					{
						TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks + current_add_count].overwrite(native_mass ,false ,i + 1);
						++current_add_count;
					}
				}
				else if (TransformedSpectrum.checked[current_modulo_id][j])
				{
					if (j == i)
					{
						TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks].overwrite(native_mass ,false ,i + 1);
						++current_add_count;
					}
					else
					{
						if (this.peak_mask[i] <= 0)
						{
							TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks].overwrite(complementary_mass ,false ,i + 1);
							++this.nb_complementary_peaks;
						}
						if (this.peak_mask[i] >= 0)
						{
							TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks + current_add_count].overwrite(native_mass ,false ,i + 1);
							++current_add_count;
						}
					}
				}
				else
				{
					TransformedSpectrum.checked[current_modulo_id][j] = true;
					if (this.peak_mask[i] <= 0)
					{
						TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks].overwrite(complementary_mass ,true ,i + 1);
						++current_add_count;
					}
					if (this.peak_mask[i] >= 0)
					{
						TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks + current_add_count].overwrite(native_mass ,true ,i + 1);
						++current_add_count;
					}
				}
				if (current_add_count == 1)
				{
					if (TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks].get_mass() <= first_mass || this.last_mass <= TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks].get_mass())
					{
						current_add_count = 0;
					}
				}
				else
				{
					if (TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks].get_mass() <= first_mass || this.last_mass <= TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks].get_mass())
					{
						peak_temp = TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks];
						TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks] = TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks + 1];
						TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks + 1] = peak_temp;
						--current_add_count;
					}
					if (TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks + current_add_count - 1].get_mass() <= first_mass || this.last_mass <= TransformedSpectrum.all_peaks[current_modulo_id][this.nb_complementary_peaks + current_add_count - 1].get_mass())
					{
						--current_add_count;
					}
				}
				this.nb_complementary_peaks += current_add_count;
			}
		}
		Arrays.sort(TransformedSpectrum.all_peaks[current_modulo_id] ,0 ,this.nb_complementary_peaks);
		
		Arrays.fill(TransformedSpectrum.ending[current_modulo_id], 0, this.nb_complementary_peaks, 0);
		Arrays.fill(TransformedSpectrum.selected_peaks[current_modulo_id], 0, this.nb_complementary_peaks, false);
		Arrays.fill(this.amino_acids_positions_number, 0);
		this.column_count = 1;
		int k;
		double difference;
		for (j = 1; j < this.nb_complementary_peaks; ++j) // Cette partie peut sans doute être améliorée !!
		{
			for (i = 0; i < j; ++i)
			{
				for (k = 0; k < AminoAcid.getCount(); ++k)
				{
					difference = (TransformedSpectrum.all_peaks[current_modulo_id][j].get_mass() - TransformedSpectrum.all_peaks[current_modulo_id][i].get_mass()) - AminoAcid.get(k).getMass();
					if (Math.abs(difference) <= Parameters.ACCURACY() && (this.amino_acids_positions_number[k] == 0 || this.amino_acids_positions[k][this.amino_acids_positions_number[k] - 1].getRightColumn() < j))
					{
						TransformedSpectrum.ending[current_modulo_id][j] |= 1 << (k + 2);
						this.amino_acids_positions[k][this.amino_acids_positions_number[k]].overwrite(j, i, TransformedSpectrum.all_peaks[current_modulo_id][j].getNativePeak(), TransformedSpectrum.all_peaks[current_modulo_id][i].getNativePeak(), TransformedSpectrum.all_peaks[0][j].get_mass(), TransformedSpectrum.all_peaks[0][i].get_mass());
						++this.amino_acids_positions_number[k];
						if(!TransformedSpectrum.selected_peaks[current_modulo_id][j])
						{
							TransformedSpectrum.selected_peaks[current_modulo_id][j] = true;
							++this.column_count;
						}
					}
				}
			}
		}
		// Le contenu de "amino_acids_positions_temp" est copie dans la matrice "this.amino_acids_positions" pour plus d'efficacite
		for (k = 0; k < AminoAcid.getCount(); ++k)
		{
			for (i = 0; i < this.amino_acids_positions_number[k]; ++i)
			{
				// L'appel a la methode "SetCondition" est lie a la version 3peaks (voir la documentation de la classe AAPosition)
				if (this.amino_acids_positions[k][i].getLeftColumn() == 0)
				{
					this.amino_acids_positions[k][i].setCondition(2);
				}
				else
				{
					this.amino_acids_positions[k][i].setCondition(TransformedSpectrum.ending[current_modulo_id][this.amino_acids_positions[k][i].getLeftColumn()]);
				}
			}
		}
		
		// Seul les pics utiles (a droite d'un AA dans le spectre) generent une colonne dans la matrice de score
		Arrays.fill(TransformedSpectrum.column_correspondance[current_modulo_id], 0);
		TransformedSpectrum.column_correspondance[current_modulo_id][0] = 1;
		Arrays.fill(this.double_column, 0, this.column_count, false);
		// "this.column_mass" sert pour le calcul des shifts
		Arrays.fill(this.column_mass, 0, this.column_count, 0.0);
		this.column_mass[0] = TransformedSpectrum.all_peaks[current_modulo_id][0].get_mass();
		int column_index = 1;
		i = 0;
		while (column_index <= this.column_count - 1)
		{
			if (TransformedSpectrum.selected_peaks[current_modulo_id][i])
			{
				if (TransformedSpectrum.all_peaks[current_modulo_id][i].isDouble())
				{
					this.double_column[column_index] = true;
				}
				this.column_mass[column_index] = TransformedSpectrum.all_peaks[current_modulo_id][i].get_mass();
				TransformedSpectrum.column_correspondance[current_modulo_id][i] = column_index + 1;
				++column_index;
			}
			++i;
		}
		
		for (k = 0; k < AminoAcid.getCount(); ++k)
		{
			for (i = 0; i < this.amino_acids_positions_number[k]; ++i)
			{
				this.amino_acids_positions[k][i].updateColumn(TransformedSpectrum.column_correspondance[current_modulo_id]);
			}
		}
	}
	
	/**
	 * The getter of the id of the spectrum
	 * @return spectrum id
	 */
	public int getID() {return this.id;}

	/**
	 * Return if the spectrum respects or not the user parameters (minimum number of peaks/amino acids). If it's not the case, it return false, otherwise true
	 * @return true if the spectrum respects the user defined parameters, otherwise false
	 */
	public boolean isUsefull() {return this.double_column != null;}
	
	/**
	 * Return the number of column needed in the matrix for the alignment of the spectrum
	 * @return the number of column
	 */
	public int getColumnCount() {return this.column_count;}
	
	// Beaucoup appelée !
	/**
	 * Return if the column in parameter correspond to a double peak or not
	 * @param column: the number of the column
	 * @return true if the column corresponds to a double peak, otherwise false
	 */
	public boolean isDoubleColumn(int column) {return this.double_column[column];}
	
	/**
	 * Return the complete list of amino acid positions in the spectrum. This is a list of list whose the first dimension correspond to the number of an amino acid 
	 * and the second dimension correspond the the different position of an amino acid in the spectrum. Each position is characterized by an AAPosition object and 
	 * the list of position for a given amino acid is sorted from the left to the right in the spectrum.
	 * @return The list of amino acid positions in the spectrum
	 */
	public AAPosition[][] getAminoAcidsPositions() {return this.amino_acids_positions;}
	
	/**
	 * Return the list that associate each amino acid to the number of time that it is present in the spectrum
	 * @return the list
	 */
	public int[] getAminoAcidsPositionsNumber() {return this.amino_acids_positions_number;}
	
	// Beaucoup appelée !
	/**
	 * Compute and return the mass shift between the 'left_column' and the left peak of the amino acid 'aa' given that the right peak of 'aa' correspond to the 'right_column'.
	 * @param left_column: beginning of the shift
	 * @param right_column: end of the shift plus 'aa'
	 * @param aa: the amino acid considered in the computation of the shift
	 * @return the shift value
	 */
	public double getShift(int left_column ,int right_column ,AminoAcid aa)
	{
		return this.column_mass[right_column] - aa.getMass() - this.column_mass[left_column];
	}
	
	/**
	 * Return the mass shift between the 'last_column' and the hypothetical last y-ion peak 
	 * @param last_column: the beginning of the shift
	 * @return the shift value
	 */
	public double getLastShift(int last_column)
	{
		double last_shift = this.last_mass - this.column_mass[last_column];
		if (last_shift > Parameters.ACCURACY()) {
			return last_shift;
		}
		else
		{
			return 0.0;
		}
	}
	
	/**
	 * The middle peak correspond to the leftmost peak for which the complementary peak is to its left. It is useful when a given alignment used at least one peak and its complementary (see algorithm.SpecGlobThread.cleanPeaksUsedSeveralTimes).
	 * @return the number of the middle peak
	 */
	public int getMiddlePeak() {return this.middle_peak;}
	
	/**
	 * Return the peak_mask list
	 * @param peak_mask list
	 */
	public int[] getPeakMask() {return this.peak_mask;}
	
	/**
	 * Setter of the attribute peak_mask
	 * @param peak_mask
	 */
	public void setPeakMask(int[] peak_mask)
	{
		if (this.peak_mask == null)
		{
			this.peak_mask = peak_mask;
		}
	} 
	
	/**
	 * Return the hypothetical error in the precursor mass.
	 * @return non_aligned_mass value
	 */
	public double getNonAlignedMass() {return this.non_aligned_mass;}
	
	/**
	 * Return the mass of the column in parameter
	 * @param column: the number of the column
	 * @return the mass value
	 */
	public double getColumnMass(int column) {return this.column_mass[column];}
	
	/**
	 * Compute the number of peaks (complementary included) that shared the mass in parameter
	 * @param mass: the researched mass
	 * @return the number of peak with the mass in parameter
	 */
	public int massCount(double mass)
	{
		double[] peaks = Spectra.getNatives()[this.id].getPeaks();
		int left = - 1;
		int right = peaks.length;
		int counter = 0;
		boolean found = false;
		while (!found && left + 1 < right)
		{
			if (mass > peaks[(left + right) / 2] + ACCURACY)
			{
				left = (left + right) / 2;
			}
			else if (mass < peaks[(left + right) / 2] - ACCURACY)
			{
				right = (left + right) / 2;
			}
			else
			{
				counter = 1;
				found = true;
			}
		}
		left = - 1;
		right = peaks.length;
		found = false;
		mass = this.complementary(mass);
		while (!found && left + 1 < right)
		{
			if (mass > peaks[(left + right) / 2] + ACCURACY)
			{
				left = (left + right) / 2;
			}
			else if (mass < peaks[(left + right) / 2] - ACCURACY)
			{
				right = (left + right) / 2;
			}
			else
			{
				++counter;
				found = true;
			}
		}
		return counter;
	}
	
	/**
	 * A method that compute the mass of the complementary peak of a peak whose the mass is in parameter
	 * @param mass: the mass of the peak of which we wish to obtain the complementary peak
	 * @return the mass of the complementary peak of the peak in parameter
	 */
	private double complementary(double mass)
	{
		return this.precursor_mass - mass + 2 * AminoAcid.getHplusMass();
	}
	
	public static void reset()
	{
		ACCURACY = Parameters.ACCURACY();
		REAL_TIME_SAVE = Parameters.REAL_TIME_SAVE();
	}
}
