package algorithms;

import constantes.Parameters;

import comparators.SimplePeakIntensityComparator;
import components.SimplePeak;

import java.util.ArrayList;
import java.util.Arrays;


/**
 * This class contains the different filtering functions to be applied to the peaks of a spectrum.
 * For now, only one method is implemented
 * @author 	BENOIST Emile, TESSIER Dominique
 */
public final class PeakFilters {
	
	private static SimplePeakIntensityComparator simple_peak_intensity_comparator = new SimplePeakIntensityComparator();
	private static double ACCURACY = Parameters.ACCURACY();

	/**
	 * Apply the filtering function specified in parameter file, on the list of peaks in parameter
	 * @param peaks: the list of peaks to filter
	 */
	public static double[] apply(ArrayList<SimplePeak> peaks)
	{
		if (Parameters.USED_FILTER().equals("mostIntense"))
		{
			peaks = PeakFilters.mostIntense(peaks);
		}
		else
		{
			// générer une erreur ici !!
			return null;
		}
		int i = 1;
		while (i < peaks.size())
		{
			if (Math.abs(peaks.get(i).getMass() - peaks.get(i - 1).getMass()) < ACCURACY)
			{
				peaks.remove(i);
			}
			else
			{
				++i;
			}
		}
		double[] masses = new double[Math.min(Parameters.NB_SELECTED_PEAKS(), peaks.size())];
		int indexMostIntense = peaks.size() - 1;
		
		for (i = 0; i < masses.length; ++i)
		{			
			masses[i] = peaks.get(indexMostIntense - i).getMass();
		}
		Arrays.sort(masses);
		return masses;
	}

	/**
	 * Filtering function that keep the most intense peaks. The number of selected peaks is specified in the parameters file.
	 *@param peaks: the list of peaks to filter
	 */
	private static ArrayList<SimplePeak> mostIntense(ArrayList<SimplePeak> peaks)
	{
		peaks.sort(simple_peak_intensity_comparator);
		return peaks;
	}
	
	public static void reset()
	{
		simple_peak_intensity_comparator = new SimplePeakIntensityComparator();
		ACCURACY = Parameters.ACCURACY();
	}
}
