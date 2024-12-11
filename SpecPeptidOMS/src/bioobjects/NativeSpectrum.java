package bioobjects;

import java.util.ArrayList;

import algorithms.PeakFilters;
import components.SimplePeak;

/**
 * This class represents a spectrum in its native form (its list of peaks is still filtered).
 * @author 	BENOIST Emile, TESSIER Dominique
 */
public class NativeSpectrum {
	
	private final String file;
	private final String title;
	private final int scan;
	private int id;
	private final double[] peaks;
	private final int charge;
	private final double precursor_mass_charged;

	/**
	 * @param file: the path to the file containing the spectrum
	 * @param title: the title of the spectrum
	 * @param scan: the scan of the spectrum
	 * @param peaks: the list of peaks
	 * @param precursor_mass_charged: the mass of the precursor charged
	 * @param charge: the charge of the precursor
	 */
	public NativeSpectrum (String file ,String title, int scan ,ArrayList<SimplePeak> peaks ,double precursor_mass_charged ,int charge)
	{
		this.file = file;
		this.title = title;
		this.scan = scan;
		this.id = 0;
		this.peaks = PeakFilters.apply(peaks);
		this.precursor_mass_charged = precursor_mass_charged;
		this.charge = charge;
	}
	
	/**
	 * A method used to fix the intern identifier of the spectrum
	 * @param id: the desired identifier  
	 */
	public void setID(int id) {this.id = id;}
	
	/**
	 * Getter of the file
	 * @return the file path
	 */
	public String getFile() {return this.file;}
	
	/**
	 * Getter of the scan
	 * @return the scan number
	 */
	public int getScan() {return this.scan;}
	
	/**
	 * Getter of the intern identifier
	 * @return the intern identifier
	 */
	public int getID() {return this.id;}
	
	/**
	 * Getter of the title
	 * @return the title
	 */
	public String getTitle() {return this.title;}
	
	/**
	 * Getter of the peaks list
	 * @return the peaks list
	 */
	public double[] getPeaks() {return this.peaks;}
	
	/**
	 * Return the peaks list size
	 * @return size of the peaks list
	 */
	public int getNumberOfNativePeaks() {return this.peaks.length;}
	
	/**
	 * Transforms the native spectrum into a transformed spectrum with its first attributes
	 * @return the transformed spectrum
	 */
	public TransformedSpectrum firstTransformation()
	{
		return new TransformedSpectrum(this.peaks, this.precursor_mass_charged, this.id, this.charge);
	}
	
	/**
	 * Transforms the native spectrum into a transformed spectrum depending on its first attributes and the parameters of this method
	 * @param spectrum: the parent transformed spectrum
	 * @param peak_mask: the peak mask to be apply on the peak list
	 * @param non_aligned_mass: the supposed error on the precursor mass
	 * @return the transformed spectrum
	 */
	public void transform(TransformedSpectrum spectrum, int[] peak_mask, double non_aligned_mass)
	{
		spectrum.overwrite(this.peaks, this.precursor_mass_charged, this.id, this.charge, peak_mask, non_aligned_mass);
	}
}
