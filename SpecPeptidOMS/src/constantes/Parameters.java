package constantes;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.regex.Pattern;

public final class Parameters {
	
	static private String file_name = "execution_parameters.ini"; // le nom du fichier contenant les parametres. Il doit se trouver directement dans le dossier du projet
	
	// Datasets paths
	static private String proteins_folder_path; // chemin d'acces au dossier contenant la banque de proteines
	static private String proteins_file; // nom du fichier contenant la banque de proteines
	static private String spectra_folder_path; // chemin d'acces au dossier contenant la banque de spectre
	static private String spectra_file; // nom du fichier contenant la banque de spectre
	static private String results_folder_path; // chemin d'acces au dossier dans lequel les resultats seront enregistres
	static private String results_file = null; // nom du fichier dans lequel les resultats seront enregistres
	
	// Spectra selection
	static private int nb_peaks_min = 10; // le nombre minimum de pics que doit contenir un spectre pour etre considere
	static private int nb_min_aa_found = 10; // le nombre minimum d'AA devant etre present dans un spectre pour que celui-ci soit considere
	
	// Used filter
	static private String used_filter = "MostIntense"; // nom du filtre a appliquer sur les spectres (ne peut etre que "MostIntense" pour l'instant)
	
	// MostIntense Filter // les parametres lies au filtre "MostIntense" applique sur les spectres
	static private int nb_selected_peaks = 60; // le nombre de pics selectionnes au maximum par spectre
	
	// General settings
	static private double accuracy = 0.02; // la precision concernant la masse des pics
	static private int nb_threads = 1; // le nombre de threads utilises
	static private int max_realignment_size = 4; // la taille maximum d'un realignement
	static private int max_realignment_size_first_column = 4;
	static private int nb_peaks_max = 500; // le taille réservé en mémoire pour chaque spectre brut (nombre de pics). Un spectre plus gros demandera une réallacation mémoire
	static private int surplus = 5;
	
	// Debug mode
	
	static private boolean shut_down_3_peaks_version = false;
	static private boolean shut_down_peaks_cleaning = false;
	static private boolean shut_down_non_aligned_mass = false;
	static private int version_preliminary_treatment = 0;
	
	// Results
	static private int real_time_save = 500; // le nombre de spectre traites entre chaque sauvegarde
	static private int nb_locations_saved = 1; // le nombre de localisations potentielles retenues par spectre
	static private int nb_interpretations_saved = 1; // le nombre d'interpretations retenus par spectre 
	static private int min_scenario_score = 10; // le score minimum que doit atteindre un scenario pour etre considere
	static private double locations_score_threshold = 0.5;
	
	// Scores main traitement
	static private int certainly_found_main = 10;
	static private int found_main = 7;
	static private int certainly_found_with_shift_main = -6;
	static private int found_with_shift_main = -8;
	static private int not_found_main = -4;
	static private int initialization_score_main;
	
	// Scores post traitement
	static private int certainly_found_post = 10;
	static private int found_post = 7;
	static private int certainly_found_with_shift_post = -6;
	static private int found_with_shift_post = -8;
	static private int not_found_post = -4;
	static private int initialization_score_post;
	
	// Modifications
	static private double GModif = 0.0;
	static private double AModif = 0.0;
	static private double SModif = 0.0;
	static private double PModif = 0.0;
	static private double VModif = 0.0;
	static private double TModif = 0.0;
	static private double CModif = 0.0;
	static private double IModif = 0.0;
	static private double LModif = 0.0;
	static private double NModif = 0.0;
	static private double DModif = 0.0;
	static private double QModif = 0.0;
	static private double KModif = 0.0;
	static private double EModif = 0.0;
	static private double MModif = 0.0;
	static private double HModif = 0.0;
	static private double FModif = 0.0;
	static private double RModif = 0.0;
	static private double YModif = 0.0;
	static private double WModif = 0.0;
	static private double UModif = 0.0;
	static private double OModif = 0.0;
	static private double NTERModif = 0.0;
	static private double CTERModif = 0.0;
	
	// Lecture du fichier de parametre
	static public void load() throws FileNotFoundException ,IOException
	{
		InputStream inputFile;
		InputStreamReader inputStream;
		BufferedReader inputBuffer;
		try
		{
			System.out.println("**********Parameter File: "+getParamFile());
			inputFile = new FileInputStream(Parameters.getParamFile());
			inputStream = new InputStreamReader(inputFile);
			inputBuffer = new BufferedReader(inputStream);
			try
			{
				String middle_pattern = "\\s*=\\s*";
				String end_pattern = "(\\s*|\\s+%.*)";
				
				String path_pattern = ".*";
				String int_pattern = "(-)?\\d+";
				String double_pattern = "\\d+(\\.\\d+)?";
				String boolean_pattern = "(true|false)";
				String filter_pattern = "mostIntense";

				String line;
				while ((line = inputBuffer.readLine()) != null)
				{
					// Datasets paths
					if (Pattern.matches("proteinsFolderPath" + middle_pattern + path_pattern + end_pattern, line))
					{
						Parameters.proteins_folder_path = line.split(middle_pattern)[1].split("\\s+")[0];
					}
					else if (Pattern.matches("proteinsFile" + middle_pattern + path_pattern + end_pattern, line))
					{
						Parameters.proteins_file = line.split(middle_pattern)[1].split("\\s+")[0];
					}
					else if (Pattern.matches("spectraFolderPath" + middle_pattern + path_pattern + end_pattern, line))
					{
						Parameters.spectra_folder_path = line.split(middle_pattern)[1].split("\\s+")[0];
					}
					else if (Pattern.matches("spectraFile" + middle_pattern + path_pattern + end_pattern, line))
					{
						Parameters.spectra_file = line.split(middle_pattern)[1].split("\\s+")[0];
					}
					else if (Pattern.matches("resultsFolderPath" + middle_pattern + path_pattern + end_pattern, line))
					{
						Parameters.results_folder_path = line.split(middle_pattern)[1].split("\\s+")[0];
					}
					else if (Pattern.matches("resultsFile" + middle_pattern + path_pattern + end_pattern, line))
					{
						Parameters.results_file = line.split(middle_pattern)[1].split("\\s+")[0];
					}
					
					// Spectra selection
					else if (Pattern.matches("nbPeaksMin" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.nb_peaks_min = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("nbMinAAFound" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.nb_min_aa_found = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					
					// Used filter
					else if (Pattern.matches("usedFilter" + middle_pattern + filter_pattern + end_pattern, line))
					{
						Parameters.used_filter = line.split(middle_pattern)[1].split("\\s+")[0];
					}
					else if (Pattern.matches("nbSelectedPeaks" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.nb_selected_peaks = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					
					// General settings						
					else if (Pattern.matches("accuracy" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.accuracy = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("nbThreads" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.nb_threads = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("tolPeakMissing" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.max_realignment_size = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("tolPeakMissingFirstCol " + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.max_realignment_size_first_column = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("nbPeaksMax" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.nb_peaks_max = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("surplus" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.surplus = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					
					// Debug mode
					else if (Pattern.matches("shutDown3PeaksVersion" + middle_pattern + boolean_pattern + end_pattern, line))
					{
						Parameters.shut_down_3_peaks_version = Boolean.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("shutDownPeaksCleaning" + middle_pattern + boolean_pattern + end_pattern, line))
					{
						Parameters.shut_down_peaks_cleaning = Boolean.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("shutDownNonAlignedMass" + middle_pattern + boolean_pattern + end_pattern, line))
					{
						Parameters.shut_down_non_aligned_mass = Boolean.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("versionPreliminaryTreatment" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.version_preliminary_treatment = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					
					// Results
					else if (Pattern.matches("nbResultsAtOnce" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.real_time_save = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("nbLssSMSaved" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.nb_locations_saved = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("nbResultsReturned" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.nb_interpretations_saved = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("minScenarioScore" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.min_scenario_score = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("filterLssSMOnScore" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.locations_score_threshold = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					
					// Score main traitement
					else if (Pattern.matches("certainlyFoundMain" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.certainly_found_main = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("foundMain" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.found_main = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("certainlyFoundWithShiftMain" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.certainly_found_with_shift_main = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("foundWithShiftMain" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.found_with_shift_main = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("notFoundMain" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.not_found_main = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					
					// Score post traitement
					else if (Pattern.matches("certainlyFoundPost" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.certainly_found_post = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("foundPost" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.found_post = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("certainlyFoundWithShiftPost" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.certainly_found_with_shift_post = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("foundWithShiftPost" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.found_with_shift_post = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("notFoundPost" + middle_pattern + int_pattern + end_pattern, line))
					{
						Parameters.not_found_post = Integer.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					
					// Modifications
					else if (Pattern.matches("GModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.GModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("AModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.AModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("SModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.SModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("PModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.PModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("VModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.VModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("TModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.TModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("CModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.CModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("IModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.IModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("LModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.LModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("NModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.NModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("DModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.DModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("QModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.QModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("KModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.KModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("EModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.EModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("MModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.MModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("HModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.HModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("FModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.FModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("RModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.RModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("YModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.YModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("WModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.WModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("UModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.UModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("OModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.OModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("NTERModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.NTERModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
					else if (Pattern.matches("CTERModif" + middle_pattern + double_pattern + end_pattern, line))
					{
						Parameters.CTERModif = Double.valueOf(line.split(middle_pattern)[1].split("\\s+")[0]);
					}
				}
				Parameters.Initialize();
			}
			finally
			{
				inputBuffer.close();
			}
		}
		catch (FileNotFoundException e)
		{
			throw new FileNotFoundException("The parameters file is not found, please check the presence of this file in the project folder ('" + Parameters.file_name + "').");
		}
	}
	
	/*
	 * Methode appelee une fois les parametres charges pour effectuer certains traitements
	 */
	static private void Initialize()
	{
		Parameters.initialization_score_main = Parameters.found_with_shift_main - Parameters.certainly_found_main - 1;
		Parameters.initialization_score_post = Parameters.found_with_shift_post - Parameters.certainly_found_post - 1;
		if (results_file == null)
		{
			Parameters.results_file = Parameters.spectra_file.split("\\.")[0] + ".csv";
		}
		Parameters.nb_threads = Math.min(Parameters.nb_threads ,Parameters.real_time_save);
		Parameters.min_scenario_score = Math.max(Parameters.min_scenario_score, 0);
	}
	
	public static String PROTEINS_FOLDER_PATH() {return Parameters.proteins_folder_path;}
	
	public static String PROTEINS_FILE() {return Parameters.proteins_file;}
	
	public static String SPECTRA_FOLDER_PATH() {return Parameters.spectra_folder_path;}
	
	public static String SPECTRA_FILE() {return Parameters.spectra_file;}
	
	public static String RESULTS_FOLDER_PATH() {return Parameters.results_folder_path;}
	
	public static String RESULTS_FILE() {return Parameters.results_file;}
	
	public static int NB_PEAKS_MIN() {return Parameters.nb_peaks_min;}
	
	public static int NB_MIN_AA_FOUND() {return Parameters.nb_min_aa_found;}
	
	public static double ACCURACY() {return Parameters.accuracy;}
	
	public static int NB_THREADS() {return Parameters.nb_threads;}
	
	public static int MAX_REALIGNMENT_SIZE() {return Parameters.max_realignment_size;}
	
	public static int MAX_REALIGNMENT_SIZE_FIRST_COLUMN() {return Parameters.max_realignment_size_first_column;}
	
	public static int NB_PEAKS_MAX() {return Parameters.nb_peaks_max;}
	
	public static int SURPLUS() {return Parameters.surplus;}
	
	public static boolean SHUT_DOWN_3_PEAKS_VERSION() {return Parameters.shut_down_3_peaks_version;}
	
	public static boolean SHUT_DOWN_PEAKS_CLEANING() {return Parameters.shut_down_peaks_cleaning;}
	
	public static boolean SHUT_DOWN_NON_ALIGNED_MASS() {return Parameters.shut_down_non_aligned_mass;}
	
	public static int VERSION_PRELIMINARY_TREATMENT() {return Parameters.version_preliminary_treatment;}
	
	public static String USED_FILTER() {return Parameters.used_filter;}
	
	public static int NB_SELECTED_PEAKS() {return Parameters.nb_selected_peaks;}
	
	public static int REAL_TIME_SAVE() {return Parameters.real_time_save;}
	
	public static int NB_LOCATIONS_SAVED() {return Parameters.nb_locations_saved;}

	public static int NB_INTERPRETATIONS_SAVED() {return Parameters.nb_interpretations_saved;}
	
	public static int MIN_SCENARIO_SCORE() {return Parameters.min_scenario_score;}
	
	public static double LOCATIONS_SCORE_THRESHOLD() {return Parameters.locations_score_threshold;}
	
	public static int CERTAINLY_FOUND_MAIN() {return Parameters.certainly_found_main;}

	public static int FOUND_MAIN() {return Parameters.found_main;}
	
	public static int CERTAINY_FOUND_WITH_SHIFT_MAIN() {return Parameters.certainly_found_with_shift_main;}
	
	public static int FOUND_WITH_SHIFT_MAIN() {return Parameters.found_with_shift_main;}
	
	public static int NOT_FOUND_MAIN() {return Parameters.not_found_main;}
	
	public static int INITIALISATION_SCORE_MAIN() {return Parameters.initialization_score_main;}
	
	public static int CERTAINLY_FOUND_POST() {return Parameters.certainly_found_post;}

	public static int FOUND_POST() {return Parameters.found_post;}
	
	public static int CERTAINY_FOUND_WITH_SHIFT_POST() {return Parameters.certainly_found_with_shift_post;}
	
	public static int FOUND_WITH_SHIFT_POST() {return Parameters.found_with_shift_post;}
	
	public static int NOT_FOUND_POST() {return Parameters.not_found_post;}
	
	public static int INITIALISATION_SCORE_POST() {return Parameters.initialization_score_post;}
	
	public static double GMODIF() {return Parameters.GModif;}
	
	public static String getParamFile() {
		return file_name;
	}
	
	public static double AMODIF() {return Parameters.AModif;}
	
	public static double SMODIF() {return Parameters.SModif;}	
	
	public static double PMODIF() {return Parameters.PModif;}
	
	public static double VMODIF() {return Parameters.VModif;}
	
	public static double TMODIF() {return Parameters.TModif;}
	
	public static double CMODIF() {return Parameters.CModif;}
	
	public static double IMODIF() {return Parameters.IModif;}
	
	public static double LMODIF() {return Parameters.LModif;}
	
	public static double NMODIF() {return Parameters.NModif;}
	
	public static double DMODIF() {return Parameters.DModif;}
	
	public static double QMODIF() {return Parameters.QModif;}
	
	public static double KMODIF() {return Parameters.KModif;}
	
	public static double EMODIF() {return Parameters.EModif;}
	
	public static double MMODIF() {return Parameters.MModif;}
	
	public static double HMODIF() {return Parameters.HModif;}
	
	public static double FMODIF() {return Parameters.FModif;}
	
	public static double RMODIF() {return Parameters.RModif;}
	
	public static double YMODIF() {return Parameters.YModif;}
	
	public static double WMODIF() {return Parameters.WModif;}
	
	public static double UMODIF() {return Parameters.UModif;}
	
	public static double OMODIF() {return Parameters.OModif;}
	
	public static double NTERMODIF() {return Parameters.NTERModif;}
	
	public static double CTERMODIF() {return Parameters.CTERModif;}
	
	// Parameters can be modified by the user through GUI
	
	public static void setAModif(double aModif) {
		AModif = aModif;
	}
	
	public static void setAccuracy(double accuracy) {
		Parameters.accuracy = accuracy;
	}
	
	public static void setCModif(double cModif) {
		CModif = cModif;
	}
	
	public static void setDModif(double dModif) {
		DModif = dModif;
	}
	
	public static void setEModif(double eModif) {
		EModif = eModif;
	}

	public static void setFModif(double fModif) {
		FModif = fModif;
	}
	
	public static void setGModif(double gModif) {
		GModif = gModif;
	}

	public static void setHModif(double hModif) {
		HModif = hModif;
	}

	public static void setIModif(double iModif) {
		IModif = iModif;
	}
	
	public static void setKModif(double kModif) {
		KModif = kModif;
	}

	public static void setLModif(double lModif) {
		LModif = lModif;
	}	

	public static void setNModif(double nModif) {
		NModif = nModif;
	}
	
	public static void setMin_scenario_score(int min_scenario_score) {
		Parameters.min_scenario_score = min_scenario_score;
	}
	
	public static void setMModif(double mModif) {
		MModif = mModif;
	}
	public static void setNb_selected_peaks(int nb_selected_peaks) {
		Parameters.nb_selected_peaks = nb_selected_peaks;
	}
	
	public static void setNb_threads(int nb_threads) {
		Parameters.nb_threads = nb_threads;
	}
	

	public static void setOModif(double oModif) {
		OModif = oModif;
	}
	
	public static void setParamFile(String name) {
		file_name=name;
	}
	
	public static void setPModif(double pModif) {
		PModif = pModif;
	}

	public static void setProteins_file(String proteins_file) {
		Parameters.proteins_file = proteins_file;
	}

	public static void setProteins_folder_path(String proteins_folder_path) {
		Parameters.proteins_folder_path = proteins_folder_path;
	}

	public static void setQModif(double qModif) {
		QModif = qModif;
	}
	
	public static void setRModif(double rModif) {
		RModif = rModif;
	}
	
	public static void setResults_folder_path(String results_folder_path) {
		Parameters.results_folder_path = results_folder_path;
	}

	public static void setResults_file(String results_file) {
		Parameters.results_file = results_file;
	}

	public static void setSModif(double sModif) {
		SModif = sModif;
	}

	public static void setSpectra_file(String spectra_file) {
		Parameters.spectra_file = spectra_file;
	}

	public static void setTModif(double tModif) {
		TModif = tModif;
	}
	
	public static void setUModif(double uModif) {
		UModif = uModif;
	}

	public static void setVModif(double vModif) {
		VModif = vModif;
	}


	public static void setWModif(double wModif) {
		WModif = wModif;
	}

	public static void setYModif(double yModif) {
		YModif = yModif;
	}

	public static void setSpectra_folder_path(String spectra_folder_path) {
		Parameters.spectra_folder_path = spectra_folder_path;
	}
}
