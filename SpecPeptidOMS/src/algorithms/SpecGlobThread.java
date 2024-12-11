package algorithms;

import java.util.ArrayList;
import java.util.Arrays;

import bioobjects.NativeSpectrum;
import bioobjects.Protein;
import bioobjects.TransformedSpectrum;
import components.AAPosition;
import components.Alignment;
import components.LocationBackup;
import components.InterpretationBackup;
import components.Location;
import components.Scenario;
import constantes.AminoAcid;
import constantes.Parameters;
import datastructures.InterpretationsSaver;
import loaders.Proteins;
import loaders.Spectra;


/**
 * This class contains the methods performing the alignment between a spectrum and all the proteins.
 * This alignment is divided in several steps :<br/> 
 *  - A first alignment with a first score system on all the proteins returning the best locals sub-sequences spectrum matches (LssSM) that can match the spectrum. The number of returned locations is a user parameter.<br/> 
 *  - A second alignment with a second score system on all the best LssSM found in the previous step. The best of them are saved in the results file. The number of saved alignment is a user parameter.<br/>
 *  - For each alignment saved in the previous step, we try to improve it by supposing that there exists a non-aligned mass in it. If it is possible to improve it, the line in the results file corresponding to the alignment is completed by this new alignment.<br/> 
 *  When a spectrum is processed, the class SpecGlobExecutor the alignment method with the next spectrum until each spectrum is treated.<br/>
 *  This class implements the Runnable Interface so several spectra can be treated in parallel. The number of thread is a user parameter.
 * @author 	BENOIST Emile, TESSIER Dominique
 * 
 * @version 1.0
 */
public final class SpecGlobThread implements Runnable
{	
	
	private static int NB_LOCATIONS_SAVED = Parameters.NB_LOCATIONS_SAVED();
	private static int NB_INTERPRETATIONS_SAVED = Parameters.NB_INTERPRETATIONS_SAVED();
	
	private static double ACCURACY = Parameters.ACCURACY();
	private static int MIN_SCENARIO_SCORE = Parameters.MIN_SCENARIO_SCORE();
	private static int MAX_REALIGNMENT_SIZE = Parameters.MAX_REALIGNMENT_SIZE();
	private static int MAX_REALIGNMENT_SIZE_FIRST_COLUMN = Parameters.MAX_REALIGNMENT_SIZE_FIRST_COLUMN();
	
	/* Score system used to compute scores in the dynamic programming matrix D */
	private static int CERTAINLY_FOUND_MAIN = Parameters.CERTAINLY_FOUND_MAIN();
	private static int FOUND_MAIN = Parameters.FOUND_MAIN();
	private static int CERTAINLY_FOUND_WITH_SHIFT_MAIN = Parameters.CERTAINY_FOUND_WITH_SHIFT_MAIN();
	private static int FOUND_WITH_SHIFT_MAIN = Parameters.FOUND_WITH_SHIFT_MAIN();
	private static int NOT_FOUND_MAIN = Parameters.NOT_FOUND_MAIN();
	private static int INITIALISATION_SCORE_MAIN = Parameters.INITIALISATION_SCORE_MAIN();	
	private static int CERTAINLY_FOUND_POST = Parameters.CERTAINLY_FOUND_POST();
	private static int FOUND_POST = Parameters.FOUND_POST();
	private static int CERTAINLY_FOUND_WITH_SHIFT_POST = Parameters.CERTAINY_FOUND_WITH_SHIFT_POST();
	private static int FOUND_WITH_SHIFT_POST = Parameters.FOUND_WITH_SHIFT_POST();
	private static int NOT_FOUND_POST = Parameters.NOT_FOUND_POST();
	private static int INITIALISATION_SCORE_POST = Parameters.INITIALISATION_SCORE_POST();
	
	/* During the post-processing step, some processing may be on or off */
	private static boolean SHUT_DOWN_3_PEAKS_VERSION = Parameters.SHUT_DOWN_3_PEAKS_VERSION();
	private static boolean SHUT_DOWN_NON_ALIGNED_MASS = Parameters.SHUT_DOWN_NON_ALIGNED_MASS();
	private static boolean SHUT_DOWN_PEAKS_CLEANING = Parameters.SHUT_DOWN_PEAKS_CLEANING();
	
	private static double SMALLEST_MASS = AminoAcid.SMALLESTMASS();
	private static int SURPLUS = Parameters.SURPLUS();

	private static int number_of_proteins = Proteins.getNumberOfProteins();
	private static Protein[] proteins = Proteins.get();
	
	private NativeSpectrum current_native_spectrum;
	private TransformedSpectrum current_transformed_spectrum;
	private Alignment current_alignment;
	private Location current_location;
	private int current_position_in_bucket;
	private int[] current_peak_mask;
	private double current_non_aligned_mass;
	
	private TransformedSpectrum first_transformed_spectrum;
	private TransformedSpectrum transformed_spectrum_temp;
	
	private InterpretationsSaver interpretations_saver;
	private InterpretationBackup best_interpretation_backup;
	private ArrayList<ArrayList<Scenario>> scenarios_matrice;
	private int[] trees_position;
	private int[] interest_cells; // i : row_number / i+1 : score / i+2 : tree_id / i+3 : beginning
	private double[] current_masses;
	private int[] scenario_buffer;
	
	private Location[] locations;
	private TransformedSpectrum[] transformed_spectra_bucket;
	private Alignment[] alignments_bucket;
	
	private int[] counter_used_peaks;
	private int[] used_peaks;
	private ArrayList<Integer>[] peaks_groups;

	/**
	 *  A SpecGlobThread object is designed to perform the alignment of several spectra on all the proteins and an alignment is based on structures such as tabular or ArrayList whose directly depends on the spectrum size.
	 * To limited the used of the garbage collector, each of these structures are created only one time and will be used for the alignment of all spectra.
	 * These structures are initialized here and there length are generally computed from the largest spectrum. 
	 * @param number: number of the thread calling this method (not used for now).
	 */
	@SuppressWarnings("unchecked")
	
	public SpecGlobThread(int number)
	{
		this.current_peak_mask = new int[Spectra.getMaxNativePeakCount()];
		this.best_interpretation_backup = new InterpretationBackup();
		this.scenarios_matrice = new ArrayList<ArrayList<Scenario>>(Spectra.getMaxColumnCount() + 1);
		for (int row_number = 0; row_number < Spectra.getMaxColumnCount() + 1; ++row_number)
		{
			this.scenarios_matrice.add(new ArrayList<Scenario>(Spectra.getMaxColumnCount() + 1));
			for (int column_number = 0; column_number < Spectra.getMaxColumnCount() + 1; ++column_number)
			{
				this.scenarios_matrice.get(row_number).add(new Scenario());
			}
		}
		this.trees_position = new int[(Proteins.getMaxLength() + 1)*(Spectra.getMaxColumnCount() + 1)];
		this.interest_cells = new int[Spectra.getMaxColumnCount() << 2];
		this.current_masses = new double[Proteins.getMaxLength()];
		this.scenario_buffer = new int[Spectra.getMaxColumnCount() << 2];

		this.counter_used_peaks = new int[Spectra.getMaxNativePeakCount() + 2];
		this.used_peaks = new int[(Spectra.getMaxColumnCount() - 1) << 1];
		this.peaks_groups = new ArrayList[Spectra.getMaxNativePeakCount() << 1];
		for (int i = 0; i < this.peaks_groups.length; ++i)
		{
			this.peaks_groups[i] = new ArrayList<Integer>();
		}
		
		this.locations = new Location[NB_LOCATIONS_SAVED];
		for (int i = 0; i < NB_LOCATIONS_SAVED; ++i)
		{
			this.locations[i] = new Location();
		}
		//this.bests_locations_spectra = new TransformedSpectrum[NB_LOCATIONS_SAVED];
		this.transformed_spectra_bucket = new TransformedSpectrum[Math.max(NB_LOCATIONS_SAVED + 2, NB_INTERPRETATIONS_SAVED + 4)];
		for (int i = 0; i < this.transformed_spectra_bucket.length; ++i)
		{
			this.transformed_spectra_bucket[i] = new TransformedSpectrum();
		}
		this.alignments_bucket = new Alignment[this.transformed_spectra_bucket.length];
		for (int i = 0; i < this.alignments_bucket.length; ++i)
		{
			this.alignments_bucket[i] = new Alignment();
		}
	}
	
	/**
	 * This method is called to move to the next spectrum.
	 * The role of a locations saver is to store the best alignments of a given spectrum, so, the LssSM saver is also changed at this point.
	 * @param native_spectrum the next spectrum that will be treated, in it native version
	 * @param transformed_spectrum the first transformed version of the native spectrum in parameter
	 * @param locations_saver the locations saver used to store the best alignment(s) of this next spectrum
	 */
	public void updateSpectrum(NativeSpectrum native_spectrum, TransformedSpectrum transformed_spectrum ,InterpretationsSaver locations_saver)
	{
		this.current_native_spectrum = native_spectrum;
		if (native_spectrum != null)
		{
			this.first_transformed_spectrum = transformed_spectrum;
			this.current_transformed_spectrum = transformed_spectrum;
			//this.current_peak_mask = this.first_transformed_spectrum.getPeakMask();
			this.interpretations_saver = locations_saver;
			locations_saver.updateSpectrum(this.trees_position, native_spectrum);
		}
	}

	/*
	private void afficherBucket(int end)
	{
		int counter = 0;
		System.out.print("[");
		for (TransformedSpectrum spectre : this.transformed_spectra_bucket)
		{
			if (spectre.getPeakMask().hashCode() == 368921242)
			{
				counter++;
			}
			System.out.print(spectre.getPeakMask().hashCode());
			System.out.print("    ");
		}
		System.out.println("]");
		if (counter > 1)
		{
			System.out.println("############____________###############_____________");
			System.out.println("############____________###############_____________");
			System.out.println("############____________###############_____________");
			System.out.println("############____________###############_____________");
		}
	}
	*/
	
	/**
	 * Each step of the treatment of the current spectrum are performed here (see the description of the class). 
	 */
	@Override
	public void run()
	{
		try
		{
			SpecGlobExecutor.getNext(this);
			LocationBackup[] locations_backup;
			LocationBackup location_backup;
			double locations_score_threshold;
			int[] first_peak_mask = new int[Spectra.getMaxNativePeakCount()];
			ArrayList<Double> possible_non_aligned_masses;
			while (this.current_native_spectrum != null)
			{
				System.out.println("__");
				/*
				 * The preliminary treatment correspond to the alignment of the current spectrum with all the proteins.
				 * After it, the interpretations saver contains the best locations found.
				 * Three methods have been tested. We retained the VERSION_PRELIMINARY_TREATMENT() = 1 for publication
				 */
				if (Parameters.VERSION_PRELIMINARY_TREATMENT() == 0)
				{
					this.preliminaryTreatmentv0();
				}
				else if (Parameters.VERSION_PRELIMINARY_TREATMENT() == 1)
				{
					this.preliminaryTreatmentv1();
				}
				else if (Parameters.VERSION_PRELIMINARY_TREATMENT() == 2)
				{
					this.preliminaryTreatmentv2();
				}
				
				locations_backup = this.interpretations_saver.getLocations();
				locations_score_threshold = this.interpretations_saver.getLocationsScoreThreshold();
				
				/*
				 * For each location, a second alignment is performed with a second system of scores.
				 */
				for (int current_location_number = 0; current_location_number < locations_backup.length; ++current_location_number)
				{
					this.current_position_in_bucket = current_location_number;
					location_backup = locations_backup[current_location_number];
					
					/*
					 * The verification that the location correspond to the current spectrum and not to a precedent. 
					 */
					if (location_backup.getSpectrumId() == this.current_native_spectrum.getID() && location_backup.getScore() >= locations_score_threshold)
					{
						this.current_transformed_spectrum = this.first_transformed_spectrum;
						this.current_peak_mask = first_peak_mask;
						
						this.current_location = this.locations[current_location_number];
						this.current_location.overwrite(proteins[location_backup.getProteinId()].getSequence(), location_backup.getStartLine() - 1, location_backup.getEndLine() - 1);
						this.expandMatriceHeight(this.current_location.getLength() + 1);
						
						/*
						 * The final treatment correspond to the alignment on the current location with the second system of score.
						 * The result is store in this.best_interpretation_backup and in this.scenarios_matrice
						 */
						this.finalTreatment();

						/*
						 * The verification that the alignment has a high enough score.
						 */
						if (this.best_interpretation_backup.getScore() > MIN_SCENARIO_SCORE)
						{
							this.current_alignment = this.alignments_bucket[current_location_number];
							this.current_alignment.overwrite(this.best_interpretation_backup, this.current_transformed_spectrum, this.scenarios_matrice, this.current_location);
							
							if (!SHUT_DOWN_PEAKS_CLEANING)
							{
								/*
								 * A correction can be apply to the alignment if it rely several times on the same peak in the spectrum.
								 */
								this.cleanPeaksUsedSeveralTimes();
							}
							
							/*
							 * If the alignment has a high enough score even after the previous correction, it is store with the other alignments, sorted in descending order of scores.
							 */
							if (this.alignments_bucket[this.current_position_in_bucket].getScore() > MIN_SCENARIO_SCORE)
							{
								this.current_transformed_spectrum.setPeakMask(first_peak_mask);
							}
						}
						else
						{
							this.alignments_bucket[this.current_position_in_bucket].clear();
						}
					}
					else
					{
						this.alignments_bucket[this.current_position_in_bucket].clear();
					}
				}
				
				/*
				 * A save of the best alignment is performed in this.interpretations_saver.
				 */
				Arrays.sort(this.alignments_bucket, 0, NB_LOCATIONS_SAVED);
				for (int interpretation_number = 0; interpretation_number < NB_INTERPRETATIONS_SAVED; ++interpretation_number)
				{
					this.interpretations_saver.saveSimpleAlignment(interpretation_number, this.alignments_bucket[interpretation_number]);
				}
				
				if (!SHUT_DOWN_NON_ALIGNED_MASS)
				{
					this.current_position_in_bucket = NB_INTERPRETATIONS_SAVED + 1;
					for (int interpretation_number = 0; interpretation_number < NB_INTERPRETATIONS_SAVED; ++interpretation_number)
					{
						this.current_alignment = this.alignments_bucket[interpretation_number];
						if (this.current_alignment.getScore() > MIN_SCENARIO_SCORE)
						{
							this.current_transformed_spectrum = this.current_alignment.getTransformedSpectrum();
							this.current_location = this.current_alignment.getLocation();
							this.alignments_bucket[NB_INTERPRETATIONS_SAVED].clear();
							
							/*
							 * We generate all the masses that can possibly correspond to a non-aligned mass and for each of them, a new transformed spectrum is constructed.
							 * A new alignment is then performed on them and the current location and the best one is keep.
							 */
							possible_non_aligned_masses = this.current_alignment.computePossibleNonAlignedMasses();
							for (int mass_number = 0; mass_number < possible_non_aligned_masses.size(); ++mass_number)
							{
								this.current_non_aligned_mass = possible_non_aligned_masses.get(mass_number);
								this.current_peak_mask = this.current_transformed_spectrum.getPeakMask();
								this.current_transformed_spectrum = this.transformed_spectra_bucket[NB_INTERPRETATIONS_SAVED + 1];
								this.current_native_spectrum.transform(this.transformed_spectra_bucket[NB_INTERPRETATIONS_SAVED + 1], this.current_peak_mask, this.current_non_aligned_mass);
								this.expandMatriceWidth(this.current_transformed_spectrum.getColumnCount() + 1);
								
								/*
								 * We use the second system of scores to perform the alignment.
								 */
								this.finalTreatment();
								
								/*
								 * The verification that the alignment has a high enough score.
								 */
								if (this.best_interpretation_backup.getScore() > MIN_SCENARIO_SCORE)
								{
									this.current_alignment = this.alignments_bucket[NB_INTERPRETATIONS_SAVED + 1];
									this.current_alignment.overwrite(this.best_interpretation_backup, this.current_transformed_spectrum, this.scenarios_matrice, this.current_location);
									
									if (!SHUT_DOWN_PEAKS_CLEANING)
									{
										/*
										 * A correction can be apply to the alignment if it rely several times on the same peak in the spectrum.
										 */
										this.cleanPeaksUsedSeveralTimes();
									}
									
									/*
									 * Only the best alignment is kept. The best is the one with the best score and with the most peaks in common in case of a tie.
									 */
									if (this.current_alignment.compareTo(this.alignments_bucket[NB_INTERPRETATIONS_SAVED]) == -1)
									{
										this.alignments_bucket[NB_INTERPRETATIONS_SAVED + 1] = this.alignments_bucket[NB_INTERPRETATIONS_SAVED];
										this.alignments_bucket[NB_INTERPRETATIONS_SAVED] = this.current_alignment;
										this.transformed_spectra_bucket[NB_INTERPRETATIONS_SAVED + 1] = this.transformed_spectra_bucket[NB_INTERPRETATIONS_SAVED];
										this.transformed_spectra_bucket[NB_INTERPRETATIONS_SAVED] = this.current_transformed_spectrum;
									}
								}
							}
							
							/*
							 * If the new alignment is better than the original, the new one is also saved in this.interpretation_saver. 
							 */
							this.alignments_bucket[NB_INTERPRETATIONS_SAVED].computeNbPeaksInCommon();
							if (this.alignments_bucket[interpretation_number].getScore() > this.alignments_bucket[NB_INTERPRETATIONS_SAVED].getScore() || this.alignments_bucket[interpretation_number].getNBPeaksInCommon() >= this.alignments_bucket[NB_INTERPRETATIONS_SAVED].getNBPeaksInCommon())
							{
								this.alignments_bucket[NB_INTERPRETATIONS_SAVED].clear();
							}
							this.interpretations_saver.saveNonAlignedMass(interpretation_number, this.alignments_bucket[NB_INTERPRETATIONS_SAVED]);
							
							this.current_non_aligned_mass = 0.0;
						}
					}
				}
				SpecGlobExecutor.getNext(this);
			}
		}
		catch (Exception e)
		{
			SpecGlobExecutor.Interrupt(e);
		}
	}
	
	/**
	 * (Version classique)
	 * The method performing the first alignment of the current spectrum and all the proteins.
	 * At the end of the call of this method, the interpretation_saver contains the best(s) LssSM for the current spectrum.
	 */
	private void preliminaryTreatmentv0()
	{
		Protein protein;
		int tree_id;
		int[] sequence;
		int sequence_length;
		int current_condition ,position_condition;
		int row_number;
		int aa_number;
		AminoAcid aa;
		int current_temp_scenario;
		int last_max_shift_score;
		int last_max_column;
		int last_column;
		AAPosition[][] all_aa_positions = this.current_transformed_spectrum.getAminoAcidsPositions();
		int[] all_aa_positions_number = this.current_transformed_spectrum.getAminoAcidsPositionsNumber();
		AAPosition[] aa_positions;
		int aa_positions_length;
		int aa_position_number;
		AAPosition position;
		int current_left_column ,current_right_column;
		int current_score_shift, current_score_found;
		int first_realignment_length = 0;// ,current_first_realignement;
		double current_shift = 0;
		int best_column;
		int best_score;
		int current_score;
		int column_count = this.current_transformed_spectrum.getColumnCount();
		int min_max_score = MIN_SCENARIO_SCORE;	
		
		for (int protein_number = 0; protein_number < number_of_proteins; ++protein_number)
		{
			if (false)// && row_number > 255 && row_number < 280)
			{
				System.out.print("masses   :");
				for (int j = 0; j < current_transformed_spectrum.getColumnCount(); ++j)
				{
					if (j%10 == 0 && j > 0)
					{
						System.out.print("   ");
					}
					SpecGlobThread.printCase((int)this.current_transformed_spectrum.getColumnMass(j));
				}
				System.out.println();
				System.out.print("nb       :");
				for (int j = 0; j < current_transformed_spectrum.getColumnCount(); ++j)
				{
					if (j%10 == 0 && j > 0)
					{
						System.out.print("   ");
					}
					SpecGlobThread.printCase(j);
				}
				System.out.println();
				System.out.println();
			}
			this.interpretations_saver.updateProtein(protein_number);
			protein = proteins[protein_number];
			
			tree_id = 0;
			
			/*
			 * Initialization of the first row
			 */
			this.interest_cells[1] = 0;
			for (int j = 4; j < (column_count<<2); j += 4)
			{
				this.interest_cells[j] = 0;
				this.interest_cells[j+1] = INITIALISATION_SCORE_MAIN;
			}
			
			sequence = protein.getSequence();
			sequence_length = sequence.length;
			current_condition = 3;
			
			for (row_number = 0; row_number < sequence_length;)
			{
				this.interest_cells[0] = row_number;
				aa_number = sequence[row_number];
				aa = AminoAcid.get(aa_number);
				++row_number;
				current_temp_scenario = 0;
				last_max_shift_score = -1;
				last_max_column = -1;
				last_column = -1;
				aa_positions = all_aa_positions[aa_number];
				aa_positions_length = all_aa_positions_number[aa_number];
				
				/*
				 * For each position of the current amino acid in the current spectrum.
				 */
				for (aa_position_number = 0; aa_position_number < aa_positions_length; ++aa_position_number)
				{
					position = aa_positions[aa_position_number];
					position_condition = position.getCondition();
					//System.out.println(position.toString());
					
					/*
					 * The criterion of the 3 peaks version must be respected.
					 */
					if ((position_condition & current_condition) != 0 || SHUT_DOWN_3_PEAKS_VERSION)
					{
						//System.out.println("traité");
						
						first_realignment_length = 0;
						current_left_column = position.getLeftColumn();
						current_right_column = position.getRightColumn();
						current_shift = 0.0;
						
						/*
						 * The score of a shift and a found depend on the fact that the corresponding peak in the spectrum is double or not.
						 */
						if (this.current_transformed_spectrum.isDoubleColumn(current_right_column))
						{
							current_score_shift = CERTAINLY_FOUND_WITH_SHIFT_MAIN;
							current_score_found = CERTAINLY_FOUND_MAIN;
						}
						else
						{
							current_score_shift = FOUND_WITH_SHIFT_MAIN;
							current_score_found = FOUND_MAIN;
						}
						
						
						/*
						 * First possibility :
						 * not found case
						 */
						best_column = current_right_column;
						int current_right_column4 = current_right_column<<2;
						best_score = this.interest_cells[current_right_column4 + 1] + NOT_FOUND_MAIN*(row_number - this.interest_cells[current_right_column4]);
						
						/*
						 * Second possibility :
						 * found case
						 * 
						 * The condition check that the left peak of the current AAPosition has a column
						 */
						if ((position_condition & 1) == 0)
						{   int current_left_columnInTab = (current_left_column + 1) << 2;
							if (this.interest_cells[current_left_columnInTab] < row_number - 1)
							{
								/*
								 * If the found case follow a shift in the protein. 
								 */
								current_score = this.interest_cells[current_left_columnInTab + 1] + NOT_FOUND_MAIN*((row_number - 1) - this.interest_cells[current_left_columnInTab]) + current_score_shift;
							}
							else
							{
								current_score = this.interest_cells[current_left_columnInTab + 1] + current_score_found;
							}
							/*
							 * Saving the new score if it's the new best one
							 */
							if (current_score >= best_score)
							{
								best_score = current_score;
								best_column = current_left_column + 1;
							}
						}
						
						/*
						 * Computation of the best position for a shift
						 */
						for (int j = last_column + 1; j <= current_left_column; ++j)
						{
							int j4 = j<<2;
							current_score = this.interest_cells[j4 + 1] + NOT_FOUND_MAIN*((row_number - 1) - this.interest_cells[j4]);
							if (current_score > last_max_shift_score)
							{
								last_max_shift_score = current_score;
								last_max_column = j;
							}
						}
						last_column = current_left_column;
						
						/*
						 * Third possibility :
						 * a realignment from the best shift or a shift if the realignment is impossible.
						 * Note that, if the best position for a shift is at column 0, the realignment mechanism is different.
						 * In the normal case, the algorithm rely on the key-cell just above it but in the column 0, it has to try on all possible cells above it.
						 */
						if (current_left_column >= 0)
						{
							if (MAX_REALIGNMENT_SIZE >= 2)
							{
								current_shift = this.current_transformed_spectrum.getShift(last_max_column ,current_right_column ,aa);
								if (last_max_column == 0)
								{
									if (current_score_found > best_score)
									{
										while (first_realignment_length < MAX_REALIGNMENT_SIZE && row_number - first_realignment_length >= 2 && current_shift - ACCURACY > this.current_masses[first_realignment_length])
										{
											++first_realignment_length;
										}
										if (first_realignment_length < MAX_REALIGNMENT_SIZE && row_number - first_realignment_length >= 2 && Math.abs(current_shift - this.current_masses[first_realignment_length]) <= ACCURACY)
										{
											best_score = current_score_found;
											best_column = 0;
											++first_realignment_length;
										}
										else if (current_score_shift > best_score)
										{
											best_score = current_score_shift;
											best_column = 0;
											first_realignment_length = 0;
										}
									}
								}
								else
								{
									first_realignment_length = row_number - this.interest_cells[last_max_column<<2] - 2;
									if (first_realignment_length >= 0 && first_realignment_length < MAX_REALIGNMENT_SIZE && Math.abs(this.current_masses[first_realignment_length] - current_shift) <= ACCURACY)
									{
										current_score = this.interest_cells[(last_max_column<<2) + 1] + current_score_found;
									}
									else
									{
										current_score = last_max_shift_score + current_score_shift;//this.caseScore(row_number - 1, last_max_column) + current_score_shift;
									}
									if (current_score > best_score)
									{
										best_score = current_score;
										best_column = last_max_column;
									}
								}									
							}
							else if (last_max_shift_score + current_score_shift > best_score)
							{
								best_score = last_max_shift_score + current_score_shift;
								best_column = last_max_column;
							}
						}
						
						/*
						 * If the best of the 3 choices is not the first one (rely on the cell just above), we compute a new key-cell and add it in this.scenario_buffer.
						 */
						if (best_column != current_right_column)
						{
							this.scenario_buffer[current_temp_scenario] = current_right_column;
							this.scenario_buffer[current_temp_scenario + 1] = best_score;
							if (best_column == 0)
							{
								this.trees_position[tree_id] = -1;
								this.scenario_buffer[current_temp_scenario + 2] = tree_id;
								this.scenario_buffer[current_temp_scenario + 3] = row_number - first_realignment_length - SURPLUS - (int)(current_shift/SMALLEST_MASS);
								++tree_id;
							}
							else
							{
								this.scenario_buffer[current_temp_scenario + 2] = this.interest_cells[(best_column<<2) + 2];
								this.scenario_buffer[current_temp_scenario + 3] = this.interest_cells[(best_column<<2) + 3];
							}
							current_temp_scenario += 4;
						}
					}
				}
				current_condition = (1 << (aa_number + 2)) + 3;
				
				/*
				 * When all the position of the current row are computed, the key-cells in this.scenario_buffer are stored in this.interest_cells.
				 */
				for (int i = 0; i < current_temp_scenario; i += 4)
				{
					best_column = this.scenario_buffer[i];
					this.interest_cells[best_column<<2] = row_number;
					this.interest_cells[(best_column<<2) + 1] = this.scenario_buffer[i + 1];
					this.interest_cells[(best_column<<2) + 2] = this.scenario_buffer[i + 2];
					this.interest_cells[(best_column<<2) + 3] = this.scenario_buffer[i + 3];

					/*
					 * If a key-cell contains a good enough score, a new LssSM is stored in this.interpretations_saver.
					 */
					if (this.scenario_buffer[i + 1] > min_max_score)
					{
						min_max_score = this.interpretations_saver.tryAddLocation(this.scenario_buffer[i + 3] ,row_number + SURPLUS + (int)(current_transformed_spectrum.getLastShift(best_column)/SMALLEST_MASS) ,this.scenario_buffer[i + 1] ,this.scenario_buffer[i + 2]);
					}
				}
				if (MAX_REALIGNMENT_SIZE >= 2)
				{
					/*
					 * Updating current_masses for the next row
					 */	
					for (int i = MAX_REALIGNMENT_SIZE - 1; i >= 1; --i)
					{
						this.current_masses[i] = this.current_masses[i - 1] + aa.getMass();
					}
					this.current_masses[0] = aa.getMass();
				}
				
				if (false)// && row_number > 255 && row_number < 280)
				{
					this.interest_cells[0] = row_number;
					SpecGlobThread.printCase(row_number);
					System.out.print(" " + protein.getLetter(row_number - 1) + " :");
					//System.out.print(row_number + " " + protein.getLetter(row_number - 1) + " :");
					for (int j = 0; j < current_transformed_spectrum.getColumnCount(); ++j)
					{
						if (j%10 == 0 && j > 0)
						{
							System.out.print("   ");
						}
						SpecGlobThread.printCase(this.interest_cells[(j<<2) + 1] + NOT_FOUND_MAIN*(row_number - this.interest_cells[j<<2]));
					}
					for (AAPosition aa_position : all_aa_positions[aa_number])
					{
						if (aa_position != null)
						{
							System.out.print("  " + aa_position.toString());
						}
					}
					System.out.println();
					if (row_number%10 == 0)
					{
						System.out.println();
					}
				}		
			}
		}
	}

	/*
	 * A shift with column 1 is automatically tested at first for sake of speed
	 */
	
	private void preliminaryTreatmentv1()
	{
		Protein protein;
		int tree_id;
		int[] sequence;
		int sequence_length;
		int current_condition ,position_condition;
		int row_number;
		int aa_number;
		AminoAcid aa;
		int current_temp_scenario;
		int last_max_shift_score;
		int last_max_column;
		int last_column;
		AAPosition[][] all_aa_positions = this.current_transformed_spectrum.getAminoAcidsPositions();
		int[] all_aa_positions_number = this.current_transformed_spectrum.getAminoAcidsPositionsNumber();
		AAPosition[] aa_positions;
		int aa_positions_length;
		int aa_position_number;
		AAPosition position;
		int current_left_column ,current_right_column;
		int current_score_shift, current_score_found;
		int realignment_length = 0;
		double current_shift = 0;
		int best_column;
		int best_score;
		int current_score;
		int column_count = this.current_transformed_spectrum.getColumnCount();
		int min_max_score = MIN_SCENARIO_SCORE;	

		for (int protein_number = 0; protein_number < number_of_proteins; ++protein_number)
		{
			if (false)
			{
				System.out.print("masses   :");
				for (int j = 0; j < current_transformed_spectrum.getColumnCount(); ++j)
				{
					if (j%10 == 0 && j > 0)
					{
						System.out.print("   ");
					}
					SpecGlobThread.printCase((int)this.current_transformed_spectrum.getColumnMass(j));
				}
				System.out.println();
				System.out.print("nb       :");
				for (int j = 0; j < current_transformed_spectrum.getColumnCount(); ++j)
				{
					if (j%10 == 0 && j > 0)
					{
						System.out.print("   ");
					}
					SpecGlobThread.printCase(j);
				}
				System.out.println();
				System.out.println();
			}
			this.interpretations_saver.updateProtein(protein_number);
			protein = proteins[protein_number];
			
			tree_id = 0;
			this.interest_cells[1] = 0;
			
			/*
			 * Initialization of the first row
			 */
			for (int j = 4; j < (column_count<<2); j += 4)
			{
				this.interest_cells[j] = 0;
				this.interest_cells[j+1] = INITIALISATION_SCORE_MAIN;
			}
			sequence = protein.getSequence();
			sequence_length = sequence.length;
			current_condition = 3;
			
			/*
			 *  For each amino acid in the current protein sequence
			 */
			for (row_number = 0; row_number < sequence_length;)
			{
				this.interest_cells[0] = row_number;
				aa_number = sequence[row_number];
				aa = AminoAcid.get(aa_number);
				++row_number;
				current_temp_scenario = 0;
				last_max_shift_score = -1;
				last_max_column = -1;
				last_column = -1;
				aa_positions = all_aa_positions[aa_number];
				aa_positions_length = all_aa_positions_number[aa_number];
				
				/*
				 * For each position where the current amino acid is found in the current spectrum
				 */
				for (aa_position_number = 0; aa_position_number < aa_positions_length; ++aa_position_number)
				{
					position = aa_positions[aa_position_number];
					position_condition = position.getCondition();
					/*
					 * This condition correspond to the restriction do to the 3 peaks version
					 */
					if ((position_condition & current_condition) != 0 || SHUT_DOWN_3_PEAKS_VERSION)
					{
						realignment_length = 0;
						current_left_column = position.getLeftColumn();
						current_right_column = position.getRightColumn();
						current_shift = 0.0;
						
						/*
						 * The score for a shift or a found depend on the peak of the current right column. This peak can be double or not
						 */
						if (this.current_transformed_spectrum.isDoubleColumn(current_right_column))
						{
							current_score_shift = CERTAINLY_FOUND_WITH_SHIFT_MAIN;
							current_score_found = CERTAINLY_FOUND_MAIN;
						}
						else
						{
							current_score_shift = FOUND_WITH_SHIFT_MAIN;
							current_score_found = FOUND_MAIN;
						}
						
						/*
						 * First possibility :
						 * not found case
						 */
						best_column = current_right_column;
						best_score = this.interest_cells[(current_right_column<<2) + 1] + NOT_FOUND_MAIN*(row_number - this.interest_cells[current_right_column<<2]);
						
						/*
						 * Second possibility :
						 * found case
						 * 
						 * The condition check that the left peak of the current AAPosition has a column
						 */
						if ((position_condition & 1) == 0)
						{
							/*
							 * If it relies on a cell corresponding to a not found case, the score for the found case corresponds to the score of a shift
							 */
							if (this.interest_cells[(current_left_column + 1)<<2] < row_number - 1)
							{
								current_score = this.interest_cells[((current_left_column + 1)<<2) + 1] + NOT_FOUND_MAIN*((row_number - 1) - this.interest_cells[(current_left_column + 1)<<2]) + current_score_shift;
							}
							else
							{
								current_score = this.interest_cells[((current_left_column + 1)<<2) + 1] + current_score_found;
							}
							/*
							 * Saving the new score if it's the new best one
							 */
							if (current_score >= best_score)
							{
								best_score = current_score;
								best_column = current_left_column + 1;
							}
						}
						
						/*
						 * Found the new best position for a shift
						 */
						for (int j = last_column + 1; j <= current_left_column; ++j)
						{
							current_score = this.interest_cells[(j<<2) + 1] + NOT_FOUND_MAIN*((row_number - 1) - this.interest_cells[j<<2]);
							if (current_score > last_max_shift_score)
							{
								last_max_shift_score = current_score;
								last_max_column = j;
							}
						}
						last_column = current_left_column;
						
						/*
						 * If a shift is possible, i.e. the current left column is greater than -1
						 */
						if (current_left_column >= 0)
						{
							/*
							 * Trying a realignment from the first column if it is allowed
							 */
							if (MAX_REALIGNMENT_SIZE_FIRST_COLUMN >= 2)
							{
								current_shift = this.current_transformed_spectrum.getShift(0 ,current_right_column ,aa);
								if (current_score_found > best_score)
								{
									while (realignment_length < MAX_REALIGNMENT_SIZE_FIRST_COLUMN && row_number - realignment_length >= 2 && current_shift - ACCURACY > this.current_masses[realignment_length])
									{
										++realignment_length;
									}
									/*
									 * If the realignment works, saving the new score if it's the new best one
									 */
									if (realignment_length < MAX_REALIGNMENT_SIZE_FIRST_COLUMN && row_number - realignment_length >= 2 && Math.abs(current_shift - this.current_masses[realignment_length]) <= ACCURACY)
									{
										best_score = current_score_found;
										best_column = 0;
										++realignment_length;
									}
								}
							}
							/*
							 * Trying a realignment from the best shift position if it is allowed, ...
							 */
							if (MAX_REALIGNMENT_SIZE >= 2)
							{
								/*
								 * If the best shift position is the first column, there is no need to try a realignment because it has already been tried just before
								 */
								if (last_max_column != 0)
								{
									current_shift = this.current_transformed_spectrum.getShift(last_max_column ,current_right_column ,aa);
									realignment_length = row_number - this.interest_cells[last_max_column<<2] - 2;
									if (realignment_length >= 0 && realignment_length < MAX_REALIGNMENT_SIZE && Math.abs(this.current_masses[realignment_length] - current_shift) <= ACCURACY)
									{
										/*
										 * The realignment works
										 */
										current_score = this.interest_cells[(last_max_column<<2) + 1] + current_score_found;
									}
									else
									{
										/*
										 * The realignment doesn't works
										 */
										current_score = last_max_shift_score + current_score_shift;
									}
									/*
									 * Saving the new score if it's the new best one
									 */
									if (current_score > best_score)
									{
										best_score = current_score;
										best_column = last_max_column;
									}
								}
								/*
								 * Trying a shift from the first column
								 */
								else if (current_score_shift > best_score)
								{
									best_score = current_score_shift;
									best_column = 0;
									realignment_length = 0;
								}
							}
							/*
							 * ... if not, just try a shift from the best shift position
							 */
							else if (last_max_shift_score + current_score_shift > best_score)
							{
								best_score = last_max_shift_score + current_score_shift;
								best_column = last_max_column;
							}
						}
						
						/*
						 * If the best case is not a not found case, it is stored in the differents structures
						 */
						if (best_column != current_right_column)
						{
							this.scenario_buffer[current_temp_scenario] = current_right_column;
							this.scenario_buffer[current_temp_scenario + 1] = best_score;
							if (best_column == 0)
							{
								this.trees_position[tree_id] = -1;
								this.scenario_buffer[current_temp_scenario + 2] = tree_id;
								this.scenario_buffer[current_temp_scenario + 3] = row_number - realignment_length - SURPLUS - (int)(current_shift/SMALLEST_MASS);
								++tree_id;
							}
							else
							{
								this.scenario_buffer[current_temp_scenario + 2] = this.interest_cells[(best_column<<2) + 2];
								this.scenario_buffer[current_temp_scenario + 3] = this.interest_cells[(best_column<<2) + 3];
							}
							current_temp_scenario += 4;
						}
					}
				}
				
				current_condition = (1 << (aa_number + 2)) + 3;
				/*
				 * Saving all the scenarios for the current row
				 */
				for (int i = 0; i < current_temp_scenario; i += 4)
				{
					best_column = this.scenario_buffer[i];
					this.interest_cells[best_column<<2] = row_number;
					this.interest_cells[(best_column<<2) + 1] = this.scenario_buffer[i + 1];
					this.interest_cells[(best_column<<2) + 2] = this.scenario_buffer[i + 2];
					this.interest_cells[(best_column<<2) + 3] = this.scenario_buffer[i + 3];

					/*
					 * Saving the alignment if it is qualitative
					 */
					if (this.scenario_buffer[i + 1] > min_max_score)
					{
						min_max_score = this.interpretations_saver.tryAddLocation(this.scenario_buffer[i + 3] ,row_number + SURPLUS + (int)(current_transformed_spectrum.getLastShift(best_column)/SMALLEST_MASS) ,this.scenario_buffer[i + 1] ,this.scenario_buffer[i + 2]);
					}
				}
				if (MAX_REALIGNMENT_SIZE >= 2)
				{
					/*
					 * Updating current_masses for the next row
					 */
					for (int i = MAX_REALIGNMENT_SIZE - 1; i >= 1; --i)
					{
						this.current_masses[i] = this.current_masses[i - 1] + aa.getMass();
					}
					this.current_masses[0] = aa.getMass();
				}
				
				if (false)
				{
					this.interest_cells[0] = row_number;
					SpecGlobThread.printCase(row_number);
					System.out.print(" " + protein.getLetter(row_number - 1) + " :");
					//System.out.print(row_number + " " + protein.getLetter(row_number - 1) + " :");
					for (int j = 0; j < current_transformed_spectrum.getColumnCount(); ++j)
					{
						if (j%10 == 0 && j > 0)
						{
							System.out.print("   ");
						}
						SpecGlobThread.printCase(this.interest_cells[(j<<2) + 1] + NOT_FOUND_MAIN*(row_number - this.interest_cells[j<<2]));
					}
					for (AAPosition aa_position : all_aa_positions[aa_number])
					{
						if (aa_position != null)
						{
							System.out.print("  " + aa_position.toString());
						}
					}
					System.out.println();
					if (row_number%10 == 0)
					{
						System.out.println();
					}
				}
			}
		}
	}
	
	/*
	 * All the possible shifts are tested with the possibility of gaps on each shift - this method takes much more time, but gives slightly better results. We estimated that the gain is not worth considering the additional time required to process spectra.
	 */
	private void preliminaryTreatmentv2()
	{
		Protein protein;
		int tree_id;
		int[] sequence;
		int sequence_length;
		int current_condition ,position_condition;
		int row_number;
		int aa_number;
		AminoAcid aa;
		int current_temp_scenario;
		AAPosition[][] all_aa_positions = this.current_transformed_spectrum.getAminoAcidsPositions();
		int[] all_aa_positions_number = this.current_transformed_spectrum.getAminoAcidsPositionsNumber();
		AAPosition[] aa_positions;
		int aa_positions_length;
		int aa_position_number;
		AAPosition position;
		int current_left_column ,current_right_column;
		int current_score_shift, current_score_found;
		int realignment_length = 0;// ,current_first_realignement;
		double current_shift = 0;
		int best_column;
		int best_score;
		int current_score;
		int column_count = this.current_transformed_spectrum.getColumnCount();
		int min_max_score = MIN_SCENARIO_SCORE;	
		
		for (int protein_number = 0; protein_number < number_of_proteins; ++protein_number)
		{
			if (false)// && row_number > 255 && row_number < 280)
			{
				System.out.print("masses   :");
				for (int j = 0; j < current_transformed_spectrum.getColumnCount(); ++j)
				{
					if (j%10 == 0 && j > 0)
					{
						System.out.print("   ");
					}
					SpecGlobThread.printCase((int)this.current_transformed_spectrum.getColumnMass(j));
				}
				System.out.println();
				System.out.print("nb       :");
				for (int j = 0; j < current_transformed_spectrum.getColumnCount(); ++j)
				{
					if (j%10 == 0 && j > 0)
					{
						System.out.print("   ");
					}
					SpecGlobThread.printCase(j);
				}
				System.out.println();
				System.out.println();
			}
			this.interpretations_saver.updateProtein(protein_number);
			protein = proteins[protein_number];
			
			tree_id = 0;
			this.interest_cells[1] = 0;
			
			for (int j = 4; j < (column_count<<2); j += 4)
			{
				this.interest_cells[j] = 0;
				this.interest_cells[j+1] = INITIALISATION_SCORE_MAIN;
			}
			sequence = protein.getSequence();
			sequence_length = sequence.length;
			current_condition = 3;
			for (row_number = 0; row_number < sequence_length;)
			{
				this.interest_cells[0] = row_number;
				aa_number = sequence[row_number];
				aa = AminoAcid.get(aa_number);
				++row_number;
				current_temp_scenario = 0;
				aa_positions = all_aa_positions[aa_number];
				aa_positions_length = all_aa_positions_number[aa_number];
				for (aa_position_number = 0; aa_position_number < aa_positions_length; ++aa_position_number)
				{
					position = aa_positions[aa_position_number];
					position_condition = position.getCondition();
					if ((position_condition & current_condition) != 0 || SHUT_DOWN_3_PEAKS_VERSION)
					{
						realignment_length = 0;
						current_left_column = position.getLeftColumn();
						current_right_column = position.getRightColumn();
						current_shift = 0.0;
						
						if (this.current_transformed_spectrum.isDoubleColumn(current_right_column))
						{
							current_score_shift = CERTAINLY_FOUND_WITH_SHIFT_MAIN;
							current_score_found = CERTAINLY_FOUND_MAIN;
						}
						else
						{
							current_score_shift = FOUND_WITH_SHIFT_MAIN;
							current_score_found = FOUND_MAIN;
						}
						best_column = current_right_column;
						best_score = this.interest_cells[(current_right_column<<2) + 1] + NOT_FOUND_MAIN*(row_number - this.interest_cells[current_right_column<<2]);
						
						if ((position_condition & 1) == 0)
						{
							if (this.interest_cells[(current_left_column + 1)<<2] < row_number - 1)
							{
								// -8 / -6
								current_score = this.interest_cells[((current_left_column + 1)<<2) + 1] + NOT_FOUND_MAIN*((row_number - 1) - this.interest_cells[(current_left_column + 1)<<2]) + current_score_shift;
							}
							else
							{
								// +7 / +10
								current_score = this.interest_cells[((current_left_column + 1)<<2) + 1] + current_score_found;
							}
							if (current_score >= best_score)
							{
								best_score = current_score;
								best_column = current_left_column + 1;
							}
						}
						if (current_left_column >= 0) // si "current_left_column == -1" alors cela signifie que l'AA est situe au tout debut du spectre, il n'y a donc pas de shift a faire, donc pas de realignement non plus
						{
							for (int j = 1; j <= current_left_column; ++j)
							{
								current_shift = this.current_transformed_spectrum.getShift(j ,current_right_column ,aa);
								realignment_length = row_number - this.interest_cells[j << 2] - 2; // calcul de l'ecart entre la ligne precedente et la ligne de la derniere case d'interet de la colonne (taille du realignement - 1)
								// premier cas : realignement possible
								if (realignment_length >= 0 && Math.abs(current_masses[realignment_length] - current_shift) <= ACCURACY)
								{
									current_score = this.interest_cells[(j << 2) + 1] + current_score_found;
								}
								// deuxieme cas : realignement impossible
								else
								{
									current_score = this.interest_cells[(j << 2) + 1] + NOT_FOUND_POST*((row_number - 1) - this.interest_cells[j << 2]) + current_score_shift;
								}
								
								// dans tous les cas, on compare ce que l'on a obtenu avec le meilleur scenario trouve jusqu'a maintenant
								if (current_score > best_score)
								{
									best_score = current_score;
									best_column = j;
								}
							}
							
							if (current_score_found > best_score) // etant donne que le score des cases de la premiere colonne est de 0, le score engendre par un realignement depuis l'une de ces cases engendrera un score de 7 (ou 10).
								  // il n'est donc pas necessaire de tenter un realignement si l'on a deja trouve un scenario faisant mieux que 7 (ou 10)
							{
								// on remonte le long de la colonne 0
								current_shift = this.current_transformed_spectrum.getShift(0 ,current_right_column ,aa);
								realignment_length = 0;
								while (realignment_length < MAX_REALIGNMENT_SIZE && row_number - realignment_length >= 2 && current_shift - ACCURACY > this.current_masses[realignment_length])
								{
									++realignment_length;
								}
								// premier cas : le realignement est possible, alors il devient le nouveau meilleur scenario
								if (realignment_length < MAX_REALIGNMENT_SIZE && row_number - realignment_length >= 2 && Math.abs(current_shift - this.current_masses[realignment_length]) <= ACCURACY)
								{
									++realignment_length;
									best_score = current_score_found;
									best_column = 0;
								}
								// deuxieme cas : realignement impossible, alors on regarde si un simple shift peut etre le nouveau meilleur scenario
								else if (current_score_shift > best_score)
								{
									realignment_length = 0;
									best_score = current_score_shift;
									best_column = 0;
								}
							}
						}
						
						if (best_column != current_right_column)
						{
							this.scenario_buffer[current_temp_scenario] = current_right_column;
							this.scenario_buffer[current_temp_scenario + 1] = best_score;
							if (best_column == 0)
							{
								this.trees_position[tree_id] = -1;
								this.scenario_buffer[current_temp_scenario + 2] = tree_id;
								this.scenario_buffer[current_temp_scenario + 3] = row_number - realignment_length - SURPLUS - (int)(current_shift/SMALLEST_MASS);
								++tree_id;
							}
							else
							{
								this.scenario_buffer[current_temp_scenario + 2] = this.interest_cells[(best_column<<2) + 2];
								this.scenario_buffer[current_temp_scenario + 3] = this.interest_cells[(best_column<<2) + 3];
							}
							current_temp_scenario += 4;
						}
					}
				}
				current_condition = (1 << (aa_number + 2)) + 3;
				for (int i = 0; i < current_temp_scenario; i += 4)
				{
					best_column = this.scenario_buffer[i];
					this.interest_cells[best_column<<2] = row_number;
					this.interest_cells[(best_column<<2) + 1] = this.scenario_buffer[i + 1];
					this.interest_cells[(best_column<<2) + 2] = this.scenario_buffer[i + 2];
					this.interest_cells[(best_column<<2) + 3] = this.scenario_buffer[i + 3];

					// save the interpretation if it is qualitative
					if (this.scenario_buffer[i + 1] > min_max_score)
					{
						min_max_score = this.interpretations_saver.tryAddLocation(this.scenario_buffer[i + 3] ,row_number + SURPLUS + (int)(current_transformed_spectrum.getLastShift(best_column)/SMALLEST_MASS) ,this.scenario_buffer[i + 1] ,this.scenario_buffer[i + 2]);
					}
				}
				if (MAX_REALIGNMENT_SIZE >= 2)
				{
					
					for (int i = MAX_REALIGNMENT_SIZE - 1; i >= 1; --i)
					{
						this.current_masses[i] = this.current_masses[i - 1] + aa.getMass();
					}
					this.current_masses[0] = aa.getMass();
				}
				
				if (false)// && row_number > 255 && row_number < 280)
				{
					this.interest_cells[0] = row_number;
					SpecGlobThread.printCase(row_number);
					System.out.print(" " + protein.getLetter(row_number - 1) + " :");
					//System.out.print(row_number + " " + protein.getLetter(row_number - 1) + " :");
					for (int j = 0; j < current_transformed_spectrum.getColumnCount(); ++j)
					{
						if (j%10 == 0 && j > 0)
						{
							System.out.print("   ");
						}
						SpecGlobThread.printCase(this.interest_cells[(j<<2) + 1] + NOT_FOUND_MAIN*(row_number - this.interest_cells[j<<2]));
					}
					for (AAPosition aa_position : all_aa_positions[aa_number])
					{
						if (aa_position != null)
						{
							System.out.print("  " + aa_position.toString());
						}
					}
					System.out.println();
					if (row_number%10 == 0)
					{
						System.out.println();
					}
				}
			}
		}
	}
	
	/**
	 * The method performing the second alignment of the current spectrum and all the bests LssSM.
	 * At the end of the call of this method, the best_interpretation_backup contains the best alignment for the current LssSM.
	 */
	private void finalTreatment()
	{
		int sequence_length;
		int position_condition;
		int row_number;
		int aa_number;
		AminoAcid aa;
		int current_temp_scenario;
		AAPosition[][] all_aa_positions = this.current_transformed_spectrum.getAminoAcidsPositions();
		int[] all_aa_positions_number = this.current_transformed_spectrum.getAminoAcidsPositionsNumber();
		AAPosition[] aa_positions;
		int aa_positions_length;
		int aa_position_number;
		AAPosition position;
		int current_left_column ,current_right_column;
		int current_score_shift, current_score_found;
		int realignment_length;
		double current_shift;
		int best_column;
		int best_score;
		int current_score;		
		int column_count = this.current_transformed_spectrum.getColumnCount();
		this.interest_cells[1] = 0;
		
		this.best_interpretation_backup.reinit();
		
		/*
		 * Initialization of the first row
		 */
		for (int j = 2; j < (column_count << 1); j+=2)
		{
			this.interest_cells[j] = 0;
			this.interest_cells[j+1] = INITIALISATION_SCORE_POST;
		}
		sequence_length = this.current_location.getLength();
		if (false)
		{
			System.out.print("masses   :");
			for (int j = 0; j < current_transformed_spectrum.getColumnCount(); ++j)
			{
				if (j%10 == 0 && j > 0)
				{
					System.out.print("   ");
				}
				SpecGlobThread.printCase((int)this.current_transformed_spectrum.getColumnMass(j));
			}
			System.out.println();
			System.out.print("nb       :");
			for (int j = 0; j < current_transformed_spectrum.getColumnCount(); ++j)
			{
				if (j%10 == 0 && j > 0)
				{
					System.out.print("   ");
				}
				SpecGlobThread.printCase(j);
			}
			System.out.println();
			System.out.println();
		}
		
		/*
		 *  For each amino acid in the current protein sequence
		 */
		for (row_number = 0; row_number < sequence_length;)
		{
			this.interest_cells[0] = row_number;
			aa_number = this.current_location.get(row_number);
			aa = AminoAcid.get(aa_number);
			++row_number;
			current_temp_scenario = 0;
			aa_positions = all_aa_positions[aa_number];
			aa_positions_length = all_aa_positions_number[aa_number];
			
			/*
			 * For each position of the current amino acid in the current spectrum.
			 */
			for (aa_position_number = 0; aa_position_number < aa_positions_length; ++aa_position_number)
			{
				position = aa_positions[aa_position_number];
				position_condition = position.getCondition();
				current_left_column = position.getLeftColumn();
				current_right_column = position.getRightColumn();
				
				/*
				 * The score of a shift and a found depend on the fact that the corresponding peak in the spectrum is double or not.
				 */
				if (this.current_transformed_spectrum.isDoubleColumn(current_right_column))
				{
					current_score_shift = CERTAINLY_FOUND_WITH_SHIFT_POST;
					current_score_found = CERTAINLY_FOUND_POST;
				}
				else
				{
					current_score_shift = FOUND_WITH_SHIFT_POST;
					current_score_found = FOUND_POST;
				}
				
				realignment_length = 0;
				
				/*
				 * First possibility :
				 * not found case
				 */
				best_column = current_right_column;
				best_score = this.interest_cells[(current_right_column << 1) + 1] + NOT_FOUND_POST*(row_number - this.interest_cells[current_right_column << 1]);
				
				/*
				 * Second possibility :
				 * found case
				 * 
				 * The condition check that the left peak of the current AAPosition has a column
				 */
				if ((position_condition & 1) == 0)
				{
					/*
					 * If it rely on a cell corresponding to a not found case, the score for the found case correspond to the score of a shift
					 */
					if (this.interest_cells[(current_left_column + 1)<<1] < row_number - 1)
					{
						/*
						 * If the found case follow a shift in the protein. 
						 */
						current_score = this.interest_cells[((current_left_column + 1) << 1) + 1] + NOT_FOUND_POST*((row_number - 1) - this.interest_cells[(current_left_column + 1) << 1]) + current_score_shift;
					}
					else
					{
						current_score = this.interest_cells[((current_left_column + 1) << 1) + 1] + current_score_found;
					}
					/*
					 * Saving the new score if it's the new best one
					 */
					if (current_score >= best_score)
					{
						best_score = current_score;
						best_column = current_left_column + 1;
					}
				}
				
				/*
				 * If a shift is possible, i.e. the current left column is greater than -1
				 */
				if (current_left_column >= 0)
				{
					/*
					 * Trying a realignment for each possible shift, except for the first column
					 */
					for (int j = 1; j <= current_left_column; ++j)
					{
						current_shift = this.current_transformed_spectrum.getShift(j ,current_right_column ,aa);
						realignment_length = row_number - this.interest_cells[j << 1] - 2;
						if (realignment_length >= 0 && Math.abs(current_masses[realignment_length] - current_shift) <= ACCURACY)
						{
							/*
							 * The realignment works
							 */
							current_score = this.interest_cells[(j << 1) + 1] + current_score_found;
						}
						else
						{
							/*
							 * The realignment doesn't works
							 */
							current_score = this.interest_cells[(j << 1) + 1] + NOT_FOUND_POST*((row_number - 1) - this.interest_cells[j << 1]) + current_score_shift;
						}
						
						/*
						 * Saving the new score if it's the new best one
						 */
						if (current_score > best_score)
						{
							best_score = current_score;
							best_column = j;
						}
					}
					
					/*
					 * Trying the realignment for the column 0
					 */
					if (current_score_found > best_score)
					{
						current_shift = this.current_transformed_spectrum.getShift(0 ,current_right_column ,aa);
						realignment_length = 0;
						while (row_number - realignment_length >= 2 && current_shift - ACCURACY > this.current_masses[realignment_length])
						{
							++realignment_length;
						}
						/*
						 * If the realignment works, saving the new score if it's the new best one
						 */
						if (row_number - realignment_length >= 2 && Math.abs(current_shift - this.current_masses[realignment_length]) <= ACCURACY)
						{
							++realignment_length;
							best_score = current_score_found;
							best_column = 0;
						}
						/*
						 * If the realignment doesn't works, trying a shift and saving the corresponding score if its the new best one
						 */
						else if (current_score_shift > best_score)
						{
							realignment_length = 0;
							best_score = current_score_shift;
							best_column = 0;
						}
					}
				}
				
				/*
				 * If the best case is not a not found case, it is stored in the different structures
				 */
				if (best_column != current_right_column)
				{						
					this.scenario_buffer[current_temp_scenario] = current_right_column;
					this.scenario_buffer[current_temp_scenario + 1] = best_score;
					if (best_column == 0)
					{
						this.scenarios_matrice.get(row_number).get(current_right_column).update(row_number - 1 - realignment_length, 0, position);
					}
					else
					{
						this.scenarios_matrice.get(row_number).get(current_right_column).update(this.interest_cells[best_column << 1], best_column, position);
					}
					current_temp_scenario += 2;
				}
	        }
			
			/*
			 * Saving all the scenarios for the current row
			 */
			for (int i = 0; i < current_temp_scenario; i += 2)
			{
				best_column = this.scenario_buffer[i];
				this.interest_cells[best_column << 1] = row_number;
				this.interest_cells[(best_column << 1) + 1] = this.scenario_buffer[i + 1];
				
				/*
				 * Saving the alignment if it the new best one
				 */
				this.best_interpretation_backup.tryImprove(this.interest_cells[(best_column << 1) + 1], row_number, best_column);
			}
			
			if (MAX_REALIGNMENT_SIZE >= 2)
			{
				/*
				 * Updating current_masses for the next row
				 */
				for (int i = row_number - 1; i >= 1; --i)
				{
					this.current_masses[i] = this.current_masses[i - 1] + aa.getMass();
				}
				this.current_masses[0] = aa.getMass();
			}
			
			if (false)
			{
				this.interest_cells[0] = row_number;
				//System.out.print(row_number + " " + AminoAcid.getLetter(this.current_sequence[row_number - 1]) + " :");
				//System.out.print(row_number + " " + AminoAcid.getLetter(this.current_location.get(row_number)) + " :");
				SpecGlobThread.printCase(row_number);
				System.out.print(" " + AminoAcid.getLetter(this.current_location.get(row_number - 1)) + " :");
				for (int j = 0; j < this.current_transformed_spectrum.getColumnCount(); ++j)
				{
					if (j%10 == 0 && j > 0)
					{
						System.out.print("   ");
					}
					SpecGlobThread.printCase(this.interest_cells[(j << 1) + 1] + NOT_FOUND_POST*(row_number - this.interest_cells[j<<1]));
				}
				/*
				for (AAPosition aa_position : all_aa_positions[aa_number])
				{
					if (aa_position != null)
					{
						System.out.print("  " + aa_position.toString());
					}
				}
				*/
				for (int j = 0; j < aa_positions_length; ++j)
				{
					System.out.print("  " + all_aa_positions[aa_number][j].toString());
				}
				System.out.println();
				if (row_number%10 == 0)
				{
					System.out.println();
				}
			}
		}
	}

	/*
	 * Method used to print the score of a cell with a fixed quantity of characters.
	 */
	@SuppressWarnings("unused")
	private static void printCase(int element)
	{
		if (Math.abs(element) < 10)
		{
			if (element < 0)
			{
				System.out.print("   " + element + " ");
			}
			else
			{
				System.out.print("    " + element + " ");
			}
		}
		else if (Math.abs(element) < 100)
		{
			if (element < 0)
			{
				System.out.print("  " + element + " ");
			}
			else
			{
				System.out.print("   " + element + " ");
			}
		}
		else if (Math.abs(element) < 1000)
		{
			if (element < 0)
			{
				System.out.print(" " + element + " ");
			}
			else
			{
				System.out.print("  " + element + " ");
			}
		}
		else if (Math.abs(element) < 10000)
		{
			if (element < 0)
			{
				System.out.print(element + " ");
			}
			else
			{
				System.out.print(" " + element + " ");
			}
		}
	}

	/**
	 * This method is used to extends the height of the scenarios_matrix to the height in parameter
	 * @param height: the expected height for the matrix 
	 */
	private void expandMatriceHeight(int height)
	{
		int width = this.scenarios_matrice.get(0).size();
		ArrayList<Scenario> line;
		while (height > this.scenarios_matrice.size())
		{
			this.scenarios_matrice.add(new ArrayList<Scenario>(width));
			line = this.scenarios_matrice.get(this.scenarios_matrice.size() - 1);
			for (int column_number = 0; column_number < width; ++column_number)
			{
				line.add(new Scenario());
			}
			
		}
	}
	
	/**
	 * This method is used to extends the width of the scenarios_matrix to the width in parameter
	 * @param width: the expected width for the matrix  
	 */
	private void expandMatriceWidth(int width)
	{
		int delta = width - this.scenarios_matrice.get(0).size();
		if (delta > 0)
		{
			ArrayList<Scenario> expansion;
			for (int row_number = 0; row_number < this.scenarios_matrice.size(); ++row_number)
			{
				expansion = new ArrayList<Scenario>(delta);
				for (int i = 0; i < delta; ++i)
				{
					expansion.add(new Scenario());
				}
				this.scenarios_matrice.get(row_number).addAll(expansion);
			}
		}
	}
	
	/**
	 * This method is used to check if the current alignment has at least one peak used several times. If it's the case, new alignments are computed until no peak is used several times
	 */
	private void cleanPeaksUsedSeveralTimes()
	{
		int current_group_number;
		AAPosition aaPosition;
		
		do
		{
			/*
			 * Store the peaks used several times in the current alignment and divided them in several groups.
			 * A group contains peaks that are consecutive in the current alignment
			 */
			boolean has_redundant_peaks = false;
			Arrays.fill(this.counter_used_peaks, 0);
			int peak_index = 0;
			Scenario current_scenario;
			int scenario_number = 0;
			do
			{
				++scenario_number;
				current_scenario = this.current_alignment.getScenarios().get(scenario_number);
				aaPosition = current_scenario.getPosition();
				++this.counter_used_peaks[aaPosition.getRightPeak()];
				if (this.counter_used_peaks[aaPosition.getRightPeak()] == 2)
				{
					has_redundant_peaks = true;
				}
				this.used_peaks[peak_index] = aaPosition.getRightPeak();
				++peak_index;
				//if ((aaPosition.getCondition() & 1) == 1) // si la version 3 pics n'est plus appliquée, alors cette condition n'est plus valide ! A remplacer par ça (current_scenario.getPreviousColumnNumber() - 1 == current_scenario.getPosition().getLeftColumn()) si jamais
				if (current_scenario.getPreviousColumnNumber() != aaPosition.getLeftColumn() + 1)
				{
					++this.counter_used_peaks[aaPosition.getLeftPeak()];
					if (this.counter_used_peaks[aaPosition.getLeftPeak()] == 2)
					{
						has_redundant_peaks = true;
					}
					this.used_peaks[peak_index] = aaPosition.getLeftPeak();
					++peak_index;
				}
			}while (current_scenario.getPreviousColumnNumber() > 0);
			
			current_group_number = 0;
			if (has_redundant_peaks)
			{			
				int middle_peak = this.current_transformed_spectrum.getMiddlePeak();
				int left_index = 0, right_index = peak_index - 1;
				while (left_index < right_index)
				{
					while (this.counter_used_peaks[this.used_peaks[left_index]] == 1)
					{
						++left_index;
					}
					if (left_index < right_index)
					{
						while (this.used_peaks[left_index] != this.used_peaks[right_index])
						{
							--right_index;
						}
						do
						{
							if (this.used_peaks[left_index] <= middle_peak)
							{
								this.peaks_groups[current_group_number << 1].add(this.used_peaks[left_index]);
							}
							else
							{
								this.peaks_groups[(current_group_number << 1) + 1].add(this.used_peaks[left_index]);
							}
							++left_index;
							--right_index;
						}
						while (left_index < right_index && this.counter_used_peaks[this.used_peaks[left_index]] == 2 && this.used_peaks[left_index] == this.used_peaks[right_index]);
						++current_group_number;
					}
				}
			}
			
			/*
			 * If there is at least one peak used several times, computing a new alignment on a new version of the current spectrum
			 */
			if (current_group_number > 0)
			{
				this.alignments_bucket[this.current_position_in_bucket].clear();
				this.alignments_bucket[this.current_position_in_bucket + 1].clear();
				this.tryAllRedundantGroups(0, current_group_number);
				this.current_alignment = this.alignments_bucket[this.current_position_in_bucket + 1];
				if (this.current_alignment.getScore() > MIN_SCENARIO_SCORE)
				{
					this.alignments_bucket[this.current_position_in_bucket + 1] = this.alignments_bucket[this.current_position_in_bucket];
					this.alignments_bucket[this.current_position_in_bucket] = this.current_alignment;
					this.current_transformed_spectrum = this.transformed_spectra_bucket[this.current_position_in_bucket + 1];
					this.transformed_spectra_bucket[this.current_position_in_bucket + 1] = this.transformed_spectra_bucket[this.current_position_in_bucket];
					this.transformed_spectra_bucket[this.current_position_in_bucket] = this.current_transformed_spectrum;
					this.current_peak_mask = this.current_transformed_spectrum.getPeakMask();
				}
				else
				{
					current_group_number = 0;
				}
			}
		}
		while (current_group_number > 0);
	}
	
	/*
	 * A recursive method that try all possible combinations of peaks used several times. Each combination lead to a new spectrum with some removed peaks and a new alignment.
	 * At the end, the best alignment is store in the current_alignment variable
	 */
	private void tryAllRedundantGroups(int group_number, int number_of_groups)
	{
		if (group_number == number_of_groups)
		{
			/*
			 * For a given combination, a new spectrum is built and a new alignment is performed on it and the current LssSM
			 */
			this.current_transformed_spectrum = this.transformed_spectra_bucket[this.current_position_in_bucket + 2];
			this.current_native_spectrum.transform(this.current_transformed_spectrum, this.current_peak_mask, this.current_non_aligned_mass);
			this.finalTreatment();
			if (this.best_interpretation_backup.getScore() >= this.alignments_bucket[current_position_in_bucket + 1].getScore() && this.best_interpretation_backup.getScore() > MIN_SCENARIO_SCORE)
			{
				this.current_alignment = this.alignments_bucket[current_position_in_bucket + 2];
				this.current_alignment.overwrite(this.best_interpretation_backup, this.current_transformed_spectrum, this.scenarios_matrice, this.current_location);
				if (this.current_alignment.compareTo(this.alignments_bucket[current_position_in_bucket + 1]) == -1)
				{
					/*
					 * If this new alignment has a good enough score and enough peaks in common with the new spectrum, we store this alignment
					 */
					this.alignments_bucket[current_position_in_bucket + 2] = this.alignments_bucket[current_position_in_bucket + 1];
					this.alignments_bucket[current_position_in_bucket + 1] = this.current_alignment;
					this.transformed_spectra_bucket[current_position_in_bucket + 2] = this.transformed_spectra_bucket[current_position_in_bucket + 1];
					this.transformed_spectra_bucket[current_position_in_bucket + 1] = this.current_transformed_spectrum;
				}
			}
		}
		else	
		{
			ArrayList<Integer> group1 = this.peaks_groups[group_number << 1];
			ArrayList<Integer> group2 = this.peaks_groups[(group_number << 1) + 1];
			for (int i = 0; i < group1.size(); ++i)
			{
				this.current_peak_mask[group1.get(i) - 1] = 1;
			}
			for (int i = 0; i < group2.size(); ++i)
			{
				this.current_peak_mask[group2.get(i) - 1] = -1;
			}
			this.tryAllRedundantGroups(group_number + 1, number_of_groups);
			for (int i = 0; i < group1.size(); ++i)
			{
				this.current_peak_mask[group1.get(i) - 1] = -1;
			}
			for (int i = 0; i < group2.size(); ++i)
			{
				this.current_peak_mask[group2.get(i) - 1] = 1;
			}
			this.tryAllRedundantGroups(group_number + 1, number_of_groups);
			for (int i = 0; i < group1.size(); ++i)
			{
				this.current_peak_mask[group1.get(i) - 1] = 0;
			}
			for (int i = 0; i < group2.size(); ++i)
			{
				this.current_peak_mask[group2.get(i) - 1] = 0;
			}
			this.peaks_groups[group_number << 1].clear();
			this.peaks_groups[(group_number << 1) + 1].clear();
		}
	}
	
	public static void reset()
	{
		NB_LOCATIONS_SAVED = Parameters.NB_LOCATIONS_SAVED();
		NB_INTERPRETATIONS_SAVED = Parameters.NB_INTERPRETATIONS_SAVED();
		
		ACCURACY = Parameters.ACCURACY();
		MIN_SCENARIO_SCORE = Parameters.MIN_SCENARIO_SCORE();
		MAX_REALIGNMENT_SIZE = Parameters.MAX_REALIGNMENT_SIZE();
		MAX_REALIGNMENT_SIZE_FIRST_COLUMN = Parameters.MAX_REALIGNMENT_SIZE_FIRST_COLUMN();
		
		CERTAINLY_FOUND_MAIN = Parameters.CERTAINLY_FOUND_MAIN();
		FOUND_MAIN = Parameters.FOUND_MAIN();
		CERTAINLY_FOUND_WITH_SHIFT_MAIN = Parameters.CERTAINY_FOUND_WITH_SHIFT_MAIN();
		FOUND_WITH_SHIFT_MAIN = Parameters.FOUND_WITH_SHIFT_MAIN();
		NOT_FOUND_MAIN = Parameters.NOT_FOUND_MAIN();
		INITIALISATION_SCORE_MAIN = Parameters.INITIALISATION_SCORE_MAIN();
		
		CERTAINLY_FOUND_POST = Parameters.CERTAINLY_FOUND_POST();
		FOUND_POST = Parameters.FOUND_POST();
		CERTAINLY_FOUND_WITH_SHIFT_POST = Parameters.CERTAINY_FOUND_WITH_SHIFT_POST();
		FOUND_WITH_SHIFT_POST = Parameters.FOUND_WITH_SHIFT_POST();
		NOT_FOUND_POST = Parameters.NOT_FOUND_POST();
		INITIALISATION_SCORE_POST = Parameters.INITIALISATION_SCORE_POST();
		
		SHUT_DOWN_3_PEAKS_VERSION = Parameters.SHUT_DOWN_3_PEAKS_VERSION();
		SHUT_DOWN_NON_ALIGNED_MASS = Parameters.SHUT_DOWN_NON_ALIGNED_MASS();
		SHUT_DOWN_PEAKS_CLEANING = Parameters.SHUT_DOWN_PEAKS_CLEANING();
		
		SMALLEST_MASS = AminoAcid.SMALLESTMASS();
		SURPLUS = Parameters.SURPLUS();

		number_of_proteins = Proteins.getNumberOfProteins();
		proteins = Proteins.get();
	}
}