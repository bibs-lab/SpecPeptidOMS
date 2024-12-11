package loaders;

import constantes.*;
import main.Main;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import GUI.SpecProtGUI;
import bioobjects.NativeSpectrum;
import bioobjects.TransformedSpectrum;
import components.SimplePeak;

import uk.ac.ebi.pride.tools.jmzreader.JMzReaderException;
import uk.ac.ebi.pride.tools.mzdata_wrapper.MzMlWrapper;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;

//import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;
//import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile.MzXMLScanIterator;
//import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLParsingException;
//import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLSpectrum;
//import uk.ac.ebi.pride.tools.mzxml_parser.mzxml.model.ParentFile;
//import uk.ac.ebi.pride.tools.mzxml_parser.mzxml.model.Scan;

/*
 * cette classe sert a charger un fichier de spectre au format 'mgf'. Les spectres sont pretraite puis stocke ici
 */
public final class Spectra {
	
	private static int NB_PEAKS_MIN = Parameters.NB_PEAKS_MIN();

	private static ArrayList<NativeSpectrum> first_native_list = new ArrayList<NativeSpectrum>(1000);
	private static NativeSpectrum[] native_list; // liste des spectres
	private static TransformedSpectrum[] transformed_list;
	private static int max_native_peak_count;
	private static int max_column_count; // le nombre maximum de colonnes (matrice d'alignement) generees par un spectre
	private static Pattern scan_pattern = Pattern.compile("scan=\\d+");
	private static Matcher matcher;
	
	public static void load() throws FileNotFoundException ,IOException, JMzReaderException
	{
		String path = Parameters.SPECTRA_FOLDER_PATH() + "/" + Parameters.SPECTRA_FILE();
		if (path.endsWith(".mgf"))
		{
			InputStream inputFile;
			InputStreamReader inputStream;
			BufferedReader input_buffer;
			try
			{
				inputFile = new FileInputStream(path);
				inputStream = new InputStreamReader(inputFile);
				input_buffer = new BufferedReader(inputStream);
			}
			catch (FileNotFoundException exception)
			{
				throw new FileNotFoundException("The spectra file \"" + path + "\" doesn't exists.");
			}
			try
			{
				Spectra.loadFromMGF(input_buffer);
			}
			catch (IOException e)
			{
				throw new IOException("The program stopped due to an IOException during the loading of the spectra file.");
			}			
			finally
			{
				input_buffer.close();
			}
		}
		else if (path.endsWith(".mzML"))
		{
			File inputFile;
			//MzXMLFile xmlFile;
			MzMlWrapper wrapper;
			try
			{
				inputFile = new File(path);
				//xmlFile = new MzXMLFile(inputFile);
				wrapper = new MzMlWrapper(inputFile);
				//Spectra.LoadFromMZML(xmlFile);
				Spectra.LoadFromMZML(wrapper);
			}
			catch (JMzReaderException e)
			{
				throw e;
			}
		}
		Spectra.computeTransformedSpectra();
		if (Main.commandMode)
		    System.out.println(Spectra.native_list.length + " spectra have been loaded.");
		else {
			String info = Spectra.native_list.length + " spectra have been loaded.\n";
			SpecProtGUI.LOG.append(info);
		}
	}
	
	private static void computeTransformedSpectra()
	{
		TransformedSpectrum.initialize();
		ArrayList<TransformedSpectrum> transformed_list_temp = new ArrayList<TransformedSpectrum>();
		TransformedSpectrum transformed_spectrum;
		int counter = 0;
		for (int spectrum_number = 0; spectrum_number < Spectra.first_native_list.size(); ++spectrum_number)
		{
			Spectra.first_native_list.get(spectrum_number).setID(counter);
			transformed_spectrum = Spectra.first_native_list.get(spectrum_number).firstTransformation();
			if (transformed_spectrum.isUsefull())
			{
				if ((counter%1000 == 0) && (counter >0))
				{
					System.out.println(counter+" spectra loaded");
					if (!Main.commandMode)
						SpecProtGUI.LOG.append(counter+" spectra loaded\n");
				}

				
				Spectra.first_native_list.set(counter, Spectra.first_native_list.get(spectrum_number));
				transformed_list_temp.add(transformed_spectrum);
				counter++;
				Spectra.max_column_count = Math.max(Spectra.max_column_count, transformed_spectrum.getColumnCount());
			}
		}

		
		Spectra.native_list = new NativeSpectrum[counter];
		for (int i = 0; i < counter; ++i)
		{
			Spectra.native_list[i] = Spectra.first_native_list.get(i);
		}
		Spectra.transformed_list = new TransformedSpectrum[counter];
		for (int i = 0; i < counter; ++i)
		{
			Spectra.transformed_list[i] = transformed_list_temp.get(i);
		}
	}
	
	public static NativeSpectrum[] getNatives() {return Spectra.native_list;}
	
	public static TransformedSpectrum[] getTransformed() {return Spectra.transformed_list;}
	
	public static int getMaxNativePeakCount() {return Spectra.max_native_peak_count;}
	
	public static int getMaxColumnCount() {return Spectra.max_column_count;}
	
	public static int getNumberOfSpectra() {return Spectra.native_list.length;}

	private static void loadFromMGF(BufferedReader inputBuffer) throws IOException
	{
		String line;
		String title = "";
		double precursor_mass_charge = 0.0;
		int charge = 0;
		int scan = -1;
		//double retention_time = 0.0;
		ArrayList<SimplePeak> peaks = new ArrayList<SimplePeak>(Parameters.NB_PEAKS_MAX());		
		String[] values;
		double mass;
		NativeSpectrum native_spectrum;
		while ((line = inputBuffer.readLine()) != null)
		{
			if (!line.startsWith("BEGIN IONS"))
			{
				if (line.startsWith("TITLE="))
				{
					title = line.substring(6);
				}
				else if ((line.startsWith("SCANS=")) || (line.startsWith("SCAN=")))
				{
					scan = Integer.valueOf(line.substring(6));
				}
				/*
				else if (line.startsWith("RTINSECONDS="))
				{
					retention_time = Double.valueOf(line.substring(12));
				}
				*/
				else if (line.startsWith("PEPMASS="))
				{
					precursor_mass_charge = Double.valueOf(line.substring(8).split(" ")[0]);
				}
				else if (line.startsWith("CHARGE="))
				{
					charge = Integer.valueOf(line.substring(7, 8));
				}
				else if (line.startsWith("END IONS"))
				{
					if (scan == -1)
					{
						matcher = Spectra.scan_pattern.matcher(title);
						if (matcher.find())
						{
							scan = Integer.valueOf(matcher.group().substring(5));
						}
					}
					/*
					if (counter%1000 == 0)
					{
						System.out.println(counter);
					}
					*/
					native_spectrum = new NativeSpectrum(Parameters.SPECTRA_FOLDER_PATH() + "/" + Parameters.SPECTRA_FILE() ,title, scan ,peaks ,precursor_mass_charge ,charge);
					if (peaks.size() >= NB_PEAKS_MIN)
					{
						Spectra.first_native_list.add(native_spectrum);
						Spectra.max_native_peak_count = Math.max(Spectra.max_native_peak_count, native_spectrum.getNumberOfNativePeaks());
					}

					peaks.clear();
				}
				else if (line.length() != 0 && Character.isDigit(line.charAt(0)))
				{
					values = line.split(" ");
					mass = Double.valueOf(values[0]);
					peaks.add(new SimplePeak(mass ,Double.valueOf(values[1])));
				}
			}
			else
			{
				title = "";
				scan = -1;
				precursor_mass_charge = 0.0;
				charge = 0;
			}
		}
	}
	
	//private static void LoadFromMZML(MzXMLFile file) throws MzXMLParsingException
	private static void LoadFromMZML(MzMlWrapper wrapper)
	{
		NativeSpectrum native_spectrum;
		ArrayList<SimplePeak> peaks = new ArrayList<SimplePeak>(Parameters.NB_PEAKS_MAX());		
		Iterator<Spectrum> iter = wrapper.getSpectrumIterator();
        while (iter.hasNext()) {
            Spectrum spectrum = iter.next();
            if (spectrum != null && spectrum.getMsLevel() == 2)
            {
            	for (Map.Entry<Double, Double> entry : spectrum.getPeakList().entrySet())
            	{
            		peaks.add(new SimplePeak(entry.getKey(), entry.getValue()));
            	}
            	native_spectrum = new NativeSpectrum(Parameters.SPECTRA_FOLDER_PATH() + "/" + Parameters.SPECTRA_FILE() ,spectrum.toString(), Integer.valueOf(spectrum.getId().split("scan=")[1]) ,peaks ,spectrum.getPrecursorMZ() ,spectrum.getPrecursorCharge());
				if (peaks.size() >= NB_PEAKS_MIN)
				{
					Spectra.first_native_list.add(native_spectrum);
					Spectra.max_native_peak_count = Math.max(Spectra.max_native_peak_count, native_spectrum.getNumberOfNativePeaks());
				}
				peaks.clear();
            }
            //Assert.assertNotNull(s);
            //System.out.println(s);
            
        }
	}
	
	public static void reset()
	{
		NB_PEAKS_MIN = Parameters.NB_PEAKS_MIN();	
		first_native_list.clear();
	}
}
