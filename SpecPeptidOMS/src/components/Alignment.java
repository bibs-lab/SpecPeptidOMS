package components;

import java.text.DecimalFormat;
import java.util.ArrayList;

import bioobjects.TransformedSpectrum;
import constantes.AminoAcid;
import constantes.Parameters;
import loaders.Spectra;

public class Alignment implements Comparable<Alignment> {
	
	private static int MIN_SCENARIO_SCORE = Parameters.MIN_SCENARIO_SCORE();
	private static double ACCURACY = Parameters.ACCURACY();
	
	private static final DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.00");
	
	private ArrayList<Scenario> scenarios;
	//private int[] sequence;
	private Location location;
	private TransformedSpectrum transformed_spectrum;
	private int score;
	private int nb_peaks_in_common;
	private ArrayList<Double> possible_non_aligned_masses_nter;
	private ArrayList<Double> possible_non_aligned_masses_cter;
	private ArrayList<Double> possible_non_aligned_masses;
	
	private StringBuilder hits_modified2;
	private StringBuilder real_sequence2;
	private StringBuilder not_found_plus_shift2;
	
	public Alignment()
	{
		this.scenarios = new ArrayList<Scenario>();
		this.scenarios.add(new Scenario());
		this.possible_non_aligned_masses_nter = new ArrayList<Double>();
		this.possible_non_aligned_masses_cter = new ArrayList<Double>();
		this.possible_non_aligned_masses = new ArrayList<Double>();
		this.hits_modified2 = new StringBuilder();
		this.real_sequence2 = new StringBuilder();
		this.not_found_plus_shift2 = new StringBuilder();
		this.score = MIN_SCENARIO_SCORE;
		this.nb_peaks_in_common = -1;		
		//this.score = MIN_SCENARIO_SCORE;
	}
	
    //public Alignment(InterpretationBackup interpretation_backup, ArrayList<ArrayList<Scenario>> scenarios_matrice, int[] sequence)
	/*
	public Alignment(InterpretationBackup interpretation_backup, ArrayList<ArrayList<Scenario>> scenarios_matrice, Location location)
	{
		int current_line_number = interpretation_backup.getLineNumber();
		int current_column_number = interpretation_backup.getColumnNumber();
		Scenario current_scenario = new Scenario(current_line_number ,current_column_number);
		this.scenarios = new ArrayList<Scenario>();
		this.scenarios.add(current_scenario);
		while (current_column_number > 0)
		{
			current_scenario = scenarios_matrice.get(current_line_number).get(current_column_number);
			this.scenarios.add(new Scenario(current_scenario));
			current_line_number = current_scenario.getPreviousLineNumber();
			current_column_number = current_scenario.getPreviousColumnNumber();
		}
		//this.sequence = sequence;
		this.location = location;
		this.score = interpretation_backup.getScore();
	}
	*/
	
	public void overwrite(InterpretationBackup interpretation_backup, TransformedSpectrum transformed_spectrum, ArrayList<ArrayList<Scenario>> scenarios_matrice, Location location)
	{
		int current_line_number = interpretation_backup.getLineNumber();
		int current_column_number = interpretation_backup.getColumnNumber();
		//Scenario current_scenario = new Scenario(current_line_number ,current_column_number);
		//this.scenarios = new ArrayList<Scenario>();
		//this.scenarios.add(current_scenario);
		Scenario current_scenario;
		this.scenarios.get(0).update(current_line_number, current_column_number);
		int position = 1;
		while (current_column_number > 0)
		{
			current_scenario = scenarios_matrice.get(current_line_number).get(current_column_number);
			//this.scenarios.add(new Scenario(current_scenario));
			if (position >= this.scenarios.size())
			{
				this.scenarios.add(new Scenario(current_scenario));
			}
			else
			{
				this.scenarios.get(position).update(current_scenario);
			}
			current_line_number = current_scenario.getPreviousLineNumber();
			current_column_number = current_scenario.getPreviousColumnNumber();
			++position;
		}
		this.transformed_spectrum = transformed_spectrum;
		this.location = location;
		this.score = interpretation_backup.getScore();
		this.nb_peaks_in_common = -1;
	}

	public ArrayList<Scenario> getScenarios() {return this.scenarios;}
	
	//public int[] getSequence() {return this.sequence;}
	
	public TransformedSpectrum getTransformedSpectrum() {return this.transformed_spectrum;}
	
	public Location getLocation() {return this.location;}
	
	public int getScore() {return this.score;}
	
	public int getNBPeaksInCommon() {return this.nb_peaks_in_common;}
	
	public void clear() {
		this.score = MIN_SCENARIO_SCORE;
		this.nb_peaks_in_common = -1;
		}
	/*
	public ArrayList<Double> computePossibleNonAlignedMasses(TransformedSpectrum spectrum)
	{
		ArrayList<Double> possible_non_aligned_masses = new ArrayList<Double>();
		int current_line_number = this.scenarios.get(0).getPreviousLineNumber();
		int current_column_number = this.scenarios.get(0).getPreviousColumnNumber();
		int previous_line_number = this.scenarios.get(1).getPreviousLineNumber();
		int previous_column_number = this.scenarios.get(1).getPreviousColumnNumber();
		double shift = spectrum.getLastShift(current_column_number);
		boolean found;
		if (shift > ACCURACY && current_line_number < this.sequence.length)
		{
			found = false;
			if (current_line_number + 1 < this.sequence.length)
			{
				int realignment_length = 2;
				double masses_sum = AminoAcid.get(this.sequence[current_line_number]).getMass() + AminoAcid.get(this.sequence[current_line_number + 1]).getMass();
				while (current_line_number + realignment_length + 1 <= this.sequence.length && shift - ACCURACY > masses_sum)
				{
					masses_sum += AminoAcid.get(this.sequence[current_line_number + realignment_length]).getMass();
					++realignment_length;
				}
				if (current_line_number + realignment_length > this.sequence.length || Math.abs(shift - masses_sum) > ACCURACY)
				{
					found = true;
				}
			}
			else
			{
				found = true;
			}
			if (found)
			{
				int line_expended = current_line_number;
				shift -= AminoAcid.get(this.sequence[line_expended]).getMass();
				while (shift > ACCURACY)
				{
					possible_non_aligned_masses.add(shift);
					++line_expended;
					if (line_expended == this.sequence.length)
					{
						shift = -1.0;
					}
					else
					{
						shift -= AminoAcid.get(this.sequence[line_expended]).getMass();
					}
				}
			}
		}
		int counter = 2;
		while (previous_column_number > 0)
		{
			double mass_sum = 0.0;
			for (int j = current_line_number - 1; j > previous_line_number; --j)
			{
				mass_sum += AminoAcid.get(this.sequence[j - 1]).getMass();
			}
			
			shift = spectrum.getShift(previous_column_number ,current_column_number, AminoAcid.get(this.sequence[current_line_number - 1]));
			if (shift - mass_sum > ACCURACY)
			{
				possible_non_aligned_masses.add(shift - mass_sum);
			}
			current_line_number = previous_line_number;
			current_column_number = previous_column_number;
			if (previous_column_number > 0)
			{
				previous_line_number = this.scenarios.get(counter).getPreviousLineNumber();
				previous_column_number = this.scenarios.get(counter).getPreviousColumnNumber();
				++counter;
			}
		}
		if (current_line_number == previous_line_number + 1)
		{
			shift = spectrum.getShift(previous_column_number ,current_column_number, AminoAcid.get(this.sequence[current_line_number - 1]));
			if (Math.abs(shift) > ACCURACY && current_line_number > 1)
			{
				int line_expended = current_line_number - 2;
				shift -= AminoAcid.get(this.sequence[line_expended]).getMass();
				while (shift > ACCURACY)
				{
					possible_non_aligned_masses.add(shift);
					--line_expended;
					if (line_expended < 0)
					{
						shift = -1.0;
					}
					else
					{
						shift -= AminoAcid.get(this.sequence[line_expended]).getMass();
					}
				}
			}
		}
		return possible_non_aligned_masses;
	}
	*/
	
	public ArrayList<Double> computePossibleNonAlignedMasses()
	{
		double last_shift;
		//ArrayList<Double> possible_non_aligned_masses = new ArrayList<Double>();
		this.possible_non_aligned_masses.clear();
		this.possible_non_aligned_masses_nter.clear();
		this.possible_non_aligned_masses_cter.clear();
		int current_line_number = this.scenarios.get(0).getPreviousLineNumber();
		int current_column_number = this.scenarios.get(0).getPreviousColumnNumber();
		int previous_line_number = this.scenarios.get(1).getPreviousLineNumber();
		int previous_column_number = this.scenarios.get(1).getPreviousColumnNumber();
		double shift = this.transformed_spectrum.getLastShift(current_column_number);
		boolean found;
		if (shift > ACCURACY && current_line_number < this.location.getLength())
		{
			found = false;
			if (current_line_number + 1 < this.location.getLength())
			{
				int realignment_length = 2;
				double masses_sum = AminoAcid.get(this.location.get(current_line_number)).getMass() + AminoAcid.get(this.location.get(current_line_number + 1)).getMass();
				while (current_line_number + realignment_length + 1 <= this.location.getLength() && shift - ACCURACY > masses_sum)
				{
					masses_sum += AminoAcid.get(this.location.get(current_line_number + realignment_length)).getMass();
					++realignment_length;
				}
				if (current_line_number + realignment_length > this.location.getLength() || Math.abs(shift - masses_sum) > ACCURACY)
				{
					found = true;
				}
			}
			else
			{
				found = true;
			}
			if (found) //nter
			{
				int line_expended = current_line_number;
				shift -= AminoAcid.get(this.location.get(line_expended)).getMass();
				//while (shift > ACCURACY)
				while (shift > -AminoAcid.getHMass() - ACCURACY)
				{
					//possible_non_aligned_masses.add(shift);
					possible_non_aligned_masses_nter.add(shift);
					++line_expended;
					if (line_expended == this.location.getLength())
					{
						//shift = -1.0;
						shift = -AminoAcid.getHMass() - ACCURACY - 1.0;
					}
					else
					{
						shift -= AminoAcid.get(this.location.get(line_expended)).getMass();
					}
				}
				if (possible_non_aligned_masses_nter.size() > 0)
				{
					last_shift = possible_non_aligned_masses_nter.get(possible_non_aligned_masses_nter.size() - 1);
					if (Math.abs(last_shift - AminoAcid.getHMass()) <= ACCURACY || Math.abs(-last_shift - AminoAcid.getHMass()) <= ACCURACY)
					{
						possible_non_aligned_masses_nter.clear();
						possible_non_aligned_masses_nter.add(last_shift);
					}
					else if (last_shift < 0.0)
					{
						possible_non_aligned_masses_nter.remove(possible_non_aligned_masses_nter.size() - 1);
					}
				}
			}
		}
		int counter = 2;
		while (previous_column_number > 0)
		{
			double mass_sum = 0.0;
			for (int j = current_line_number - 1; j > previous_line_number; --j)
			{
				mass_sum += AminoAcid.get(this.location.get(j - 1)).getMass();
			}
			
			shift = this.transformed_spectrum.getShift(previous_column_number ,current_column_number, AminoAcid.get(this.location.get(current_line_number - 1)));
			if (shift - mass_sum > ACCURACY)
			{
				possible_non_aligned_masses.add(shift - mass_sum);
			}
			current_line_number = previous_line_number;
			current_column_number = previous_column_number;
			if (previous_column_number > 0)
			{
				previous_line_number = this.scenarios.get(counter).getPreviousLineNumber();
				previous_column_number = this.scenarios.get(counter).getPreviousColumnNumber();
				++counter;
			}
		}
		if (current_line_number == previous_line_number + 1)
		{
			shift = this.transformed_spectrum.getShift(previous_column_number ,current_column_number, AminoAcid.get(this.location.get(current_line_number - 1)));
			if (Math.abs(shift) > ACCURACY && current_line_number > 1) //cter
			{
				int line_expended = current_line_number - 2;
				shift -= AminoAcid.get(this.location.get(line_expended)).getMass();
				//while (shift > ACCURACY)
				while (shift > -AminoAcid.getHMass() - ACCURACY)
				{
					//possible_non_aligned_masses.add(shift);
					possible_non_aligned_masses_cter.add(shift);
					--line_expended;
					if (line_expended < 0)
					{
						//shift = -1.0;
						shift = -AminoAcid.getHMass() - ACCURACY - 1.0;
					}
					else
					{
						shift -= AminoAcid.get(this.location.get(line_expended)).getMass();
					}
				}
				if (possible_non_aligned_masses_cter.size() > 0)
				{
					last_shift = possible_non_aligned_masses_cter.get(possible_non_aligned_masses_cter.size() - 1);
					if (Math.abs(last_shift - AminoAcid.getHMass()) <= ACCURACY || Math.abs(-last_shift - AminoAcid.getHMass()) <= ACCURACY)
					{
						possible_non_aligned_masses_cter.clear();
						possible_non_aligned_masses_cter.add(last_shift);
					}
					else if (last_shift < 0.0)
					{
						possible_non_aligned_masses_cter.remove(possible_non_aligned_masses_cter.size() - 1);
					}
				}			
			}
		}
		possible_non_aligned_masses.addAll(possible_non_aligned_masses_nter);
		possible_non_aligned_masses.addAll(possible_non_aligned_masses_cter);
		return possible_non_aligned_masses;
	}
	
	/*
	public String[] computeHitModifiedAndRealSequence(TransformedSpectrum spectrum)
	{
		String hits_modified = "";
		String real_sequence = "";
		DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.00");
		Scenario current_scenario = this.scenarios.get(0);
		int current_line_number = current_scenario.getPreviousLineNumber();
		int current_column_number = current_scenario.getPreviousColumnNumber();
		current_scenario = this.scenarios.get(1);
		int previous_line_number = current_scenario.getPreviousLineNumber();
		int previous_column_number = current_scenario.getPreviousColumnNumber();
		double shift = spectrum.getLastShift(current_column_number);
		if (shift > ACCURACY)
		{
			if (current_line_number + 1 < this.sequence.length)
			{
				int realignment_length = 2;
				double masses_sum = AminoAcid.get(this.sequence[current_line_number]).getMass() + AminoAcid.get(this.sequence[current_line_number + 1]).getMass();
				while (current_line_number + realignment_length + 1 <= this.sequence.length && shift - ACCURACY > masses_sum)
				{
					masses_sum += AminoAcid.get(this.sequence[current_line_number + realignment_length]).getMass();
					++realignment_length;
				}
				if (current_line_number + realignment_length <= this.sequence.length && Math.abs(shift - masses_sum) <= ACCURACY)
				{
					hits_modified = String.format("[%c]", AminoAcid.getLetter(this.sequence[current_line_number + realignment_length - 1]));
					real_sequence = Character.toString(AminoAcid.getLetter(this.sequence[current_line_number + realignment_length - 1]));
					--realignment_length;
					for (;realignment_length > 0; --realignment_length)
					{
						hits_modified = String.format("%s[%c]", hits_modified, AminoAcid.getLetter(this.sequence[current_line_number + realignment_length - 1]));
						real_sequence = String.format("%s%c", real_sequence, AminoAcid.getLetter(this.sequence[current_line_number  + realignment_length - 1]));
					}
				}
				else
				{
					hits_modified = String.format("[%s]", DECIMAL_FORMAT.format(shift));
				}
			}
			else
			{
				hits_modified = String.format("[%s]", DECIMAL_FORMAT.format(shift));
			}
		}
		int counter = 2;
		while (current_column_number > 0)
		{
			hits_modified = String.format("%s%c", hits_modified, AminoAcid.getLetter(this.sequence[current_line_number - 1]));
			real_sequence = String.format("%s%c", real_sequence, AminoAcid.getLetter(this.sequence[current_line_number - 1]));
			String not_found_plus_shift = "";
			double mass_sum = 0.0;
			for (int j = current_line_number - 1; j > previous_line_number; --j)
			{
				mass_sum += AminoAcid.get(this.sequence[j - 1]).getMass();
				not_found_plus_shift = String.format("%s[%c]", not_found_plus_shift, AminoAcid.getLetter(this.sequence[j - 1]));
				real_sequence = String.format("%s%c", real_sequence, AminoAcid.getLetter(this.sequence[j - 1]));
			}
			hits_modified = hits_modified.concat(not_found_plus_shift);
			
			shift = spectrum.getShift(previous_column_number ,current_column_number, AminoAcid.get(this.sequence[current_line_number - 1]));
			if (Math.abs(shift - mass_sum) > ACCURACY)
			{
				hits_modified = String.format("%s[%s]", hits_modified, DECIMAL_FORMAT.format(shift - mass_sum));
			}
			current_line_number = previous_line_number;
			current_column_number = previous_column_number;
			if (previous_column_number > 0)
			{
				current_scenario = this.scenarios.get(counter);
				previous_line_number = current_scenario.getPreviousLineNumber();
				previous_column_number = current_scenario.getPreviousColumnNumber();
				++counter;
			}
		}
		String[] res = {hits_modified, real_sequence};
		return res;
	}
	*/
	
	public String[] computeHitModifiedAndRealSequence()
	{
		//String hits_modified = "";
		//String real_sequence = "";
		//StringBuilder hits_modified2 = new StringBuilder();
		//StringBuilder real_sequence2 = new StringBuilder();
		//StringBuilder not_found_plus_shift2 = new StringBuilder();
		this.hits_modified2.setLength(0);
		this.real_sequence2.setLength(0);
		//DecimalFormat DECIMAL_FORMAT = new DecimalFormat("0.00");
		Scenario current_scenario = this.scenarios.get(0);
		int current_line_number = current_scenario.getPreviousLineNumber();
		int current_column_number = current_scenario.getPreviousColumnNumber();
		current_scenario = this.scenarios.get(1);
		int previous_line_number = current_scenario.getPreviousLineNumber();
		int previous_column_number = current_scenario.getPreviousColumnNumber();
		double shift = this.transformed_spectrum.getLastShift(current_column_number);
		if (shift > ACCURACY)
		{
			if (current_line_number + 1 < this.location.getLength())
			{
				int realignment_length = 2;
				double masses_sum = AminoAcid.get(this.location.get(current_line_number)).getMass() + AminoAcid.get(this.location.get(current_line_number + 1)).getMass();
				while (current_line_number + realignment_length + 1 <= this.location.getLength() && shift - ACCURACY > masses_sum)
				{
					masses_sum += AminoAcid.get(this.location.get(current_line_number + realignment_length)).getMass();
					++realignment_length;
				}
				if (current_line_number + realignment_length <= this.location.getLength() && Math.abs(shift - masses_sum) <= ACCURACY)
				{
					//hits_modified = String.format("[%c]", AminoAcid.getLetter(this.location.get(current_line_number + realignment_length - 1)));
					//real_sequence = Character.toString(AminoAcid.getLetter(this.location.get(current_line_number + realignment_length - 1)));
					hits_modified2.append('[');
					hits_modified2.append(AminoAcid.getLetter(this.location.get(current_line_number + realignment_length - 1)));
					hits_modified2.append(']');
					real_sequence2.append(AminoAcid.getLetter(this.location.get(current_line_number + realignment_length - 1)));
					--realignment_length;
					for (;realignment_length > 0; --realignment_length)
					{
						//hits_modified = String.format("%s[%c]", hits_modified, AminoAcid.getLetter(this.location.get(current_line_number + realignment_length - 1)));
						//real_sequence = String.format("%s%c", real_sequence, AminoAcid.getLetter(this.location.get(current_line_number  + realignment_length - 1)));
						hits_modified2.append('[');
						hits_modified2.append(AminoAcid.getLetter(this.location.get(current_line_number + realignment_length - 1)));
						hits_modified2.append(']');
						real_sequence2.append(AminoAcid.getLetter(this.location.get(current_line_number  + realignment_length - 1)));
					}
				}
				else
				{
					//hits_modified = String.format("[%s]", DECIMAL_FORMAT.format(shift));
					hits_modified2.append('[');
					hits_modified2.append(DECIMAL_FORMAT.format(shift));
					hits_modified2.append(']');
				}
			}
			else
			{
				//hits_modified = String.format("[%s]", DECIMAL_FORMAT.format(shift));
				hits_modified2.append('[');
				hits_modified2.append(DECIMAL_FORMAT.format(shift));
				hits_modified2.append(']');
			}
		}
		int counter = 2;
		while (current_column_number > 0)
		{
			//hits_modified = String.format("%s%c", hits_modified, AminoAcid.getLetter(this.location.get(current_line_number - 1)));
			//real_sequence = String.format("%s%c", real_sequence, AminoAcid.getLetter(this.location.get(current_line_number - 1)));
			hits_modified2.append(AminoAcid.getLetter(this.location.get(current_line_number - 1)));
			real_sequence2.append(AminoAcid.getLetter(this.location.get(current_line_number - 1)));
			//String not_found_plus_shift = "";
			not_found_plus_shift2.setLength(0);
			double mass_sum = 0.0;
			for (int j = current_line_number - 1; j > previous_line_number; --j)
			{
				mass_sum += AminoAcid.get(this.location.get(j - 1)).getMass();
				//not_found_plus_shift = String.format("%s[%c]", not_found_plus_shift, AminoAcid.getLetter(this.location.get(j - 1)));
				not_found_plus_shift2.append('[');
				not_found_plus_shift2.append(AminoAcid.getLetter(this.location.get(j - 1)));
				not_found_plus_shift2.append(']');
				//real_sequence = String.format("%s%c", real_sequence, AminoAcid.getLetter(this.location.get(j - 1)));
				real_sequence2.append(AminoAcid.getLetter(this.location.get(j - 1)));
			}
			//hits_modified = hits_modified.concat(not_found_plus_shift);
			hits_modified2.append(not_found_plus_shift2);
			
			shift = this.transformed_spectrum.getShift(previous_column_number ,current_column_number, AminoAcid.get(this.location.get(current_line_number - 1)));
			if (Math.abs(shift - mass_sum) > ACCURACY)
			{
				//hits_modified = String.format("%s[%s]", hits_modified, DECIMAL_FORMAT.format(shift - mass_sum));
				hits_modified2.append('[');
				hits_modified2.append(DECIMAL_FORMAT.format(shift - mass_sum));
				hits_modified2.append(']');
			}
			current_line_number = previous_line_number;
			current_column_number = previous_column_number;
			if (previous_column_number > 0)
			{
				current_scenario = this.scenarios.get(counter);
				previous_line_number = current_scenario.getPreviousLineNumber();
				previous_column_number = current_scenario.getPreviousColumnNumber();
				++counter;
			}
		}
		//String[] res = {hits_modified, real_sequence};
		String[] res = {hits_modified2.toString(), real_sequence2.toString()};
		return res;
	}
	
	/*
	public void computeNbPeaksInCommon(TransformedSpectrum spectrum)
	{
		//int nb_peaks_in_common = 0;
		this.nb_peaks_in_common = 0;
		int current_line_number = this.scenarios.get(0).getPreviousLineNumber();
		int current_column_number = this.scenarios.get(0).getPreviousColumnNumber();
		int previous_line_number = this.scenarios.get(1).getPreviousLineNumber();
		int previous_column_number = this.scenarios.get(1).getPreviousColumnNumber();
		double shift = spectrum.getLastShift(current_column_number);
		double masses_sum;
		int realignment_length;
		if (shift > ACCURACY)
		{
			if (current_line_number + 1 < sequence.length)
			{
				realignment_length = 2;
				masses_sum = AminoAcid.get(this.sequence[current_line_number]).getMass() + AminoAcid.get(this.sequence[current_line_number + 1]).getMass();
				while (current_line_number + realignment_length + 1 <= this.sequence.length && Math.abs(shift - masses_sum) > ACCURACY)
				{
					masses_sum += AminoAcid.get(sequence[current_line_number + realignment_length]).getMass();
					++realignment_length;
				}
				if (current_line_number + realignment_length <= this.sequence.length && Math.abs(shift - masses_sum) <= ACCURACY)
				{
					++this.nb_peaks_in_common;
					masses_sum = spectrum.getColumnMass(spectrum.getColumnCount() - 1);
					for (;realignment_length > 1; --realignment_length)
					{
						masses_sum -= AminoAcid.get(this.sequence[current_line_number - 1 + realignment_length]).getMass();
						this.nb_peaks_in_common += spectrum.massCount(masses_sum);
					}
				}
			}
		}
		int counter = 2;
		while (previous_column_number > 0)
		{
			if (spectrum.isDoubleColumn(current_column_number))
			{
				this.nb_peaks_in_common += 2;
			}
			else
			{
				++this.nb_peaks_in_common;
			}
			
			shift = spectrum.getShift(previous_column_number ,current_column_number, AminoAcid.get(this.sequence[current_line_number - 1]));			
			if (shift > ACCURACY)
			{
				masses_sum = 0.0;
				for (int j = current_line_number - 1; j > previous_line_number; --j)
				{
					masses_sum += AminoAcid.get(this.sequence[j - 1]).getMass();
				}
				this.nb_peaks_in_common += spectrum.massCount(spectrum.getColumnMass(current_column_number) - AminoAcid.get(this.sequence[current_line_number - 1]).getMass());
				if (Math.abs(shift - masses_sum) <= ACCURACY)
				{
					masses_sum = spectrum.getColumnMass(previous_column_number);
					
					for (realignment_length = 1 ;realignment_length < current_line_number - previous_line_number - 1; ++realignment_length)
					{
						masses_sum += AminoAcid.get(this.sequence[previous_line_number - 1 + realignment_length]).getMass();
						this.nb_peaks_in_common += spectrum.massCount(masses_sum);
					}
				}
			}
				
			current_line_number = previous_line_number;
			current_column_number = previous_column_number;
			if (previous_column_number > 0)
			{
				previous_line_number = this.scenarios.get(counter).getPreviousLineNumber();
				previous_column_number = this.scenarios.get(counter).getPreviousColumnNumber();
				++counter;
			}
		}
		if (spectrum.isDoubleColumn(current_column_number))
		{
			this.nb_peaks_in_common += 2;
		}
		else
		{
			++this.nb_peaks_in_common;
		}
		
		shift = spectrum.getShift(previous_column_number ,current_column_number, AminoAcid.get(this.sequence[current_line_number - 1]));			
		if (shift > ACCURACY)
		{
			masses_sum = 0.0;
			for (int j = current_line_number - 1; j > previous_line_number; --j)
			{
				masses_sum += AminoAcid.get(this.sequence[j - 1]).getMass();
			}
			this.nb_peaks_in_common += spectrum.massCount(spectrum.getColumnMass(current_column_number) - AminoAcid.get(this.sequence[current_line_number - 1]).getMass());
			if (Math.abs(shift - masses_sum) <= ACCURACY)
			{
				++this.nb_peaks_in_common;
				masses_sum = spectrum.getColumnMass(previous_column_number);
				
				for (realignment_length = 1 ;realignment_length < current_line_number - previous_line_number - 1; ++realignment_length)
				{
					masses_sum += AminoAcid.get(this.sequence[previous_line_number - 1 + realignment_length]).getMass();
					this.nb_peaks_in_common += spectrum.massCount(masses_sum);
				}
			}
		}
		else
		{
			++this.nb_peaks_in_common;
		}
		//return nb_peaks_in_common;
	}
	*/
	
	public void computeNbPeaksInCommon()
	{
		if (this.nb_peaks_in_common == -1 && this.score != MIN_SCENARIO_SCORE)
		{
			this.nb_peaks_in_common = 0;
			int current_line_number = this.scenarios.get(0).getPreviousLineNumber();
			int current_column_number = this.scenarios.get(0).getPreviousColumnNumber();
			int previous_line_number = this.scenarios.get(1).getPreviousLineNumber();
			int previous_column_number = this.scenarios.get(1).getPreviousColumnNumber();
			double shift = this.transformed_spectrum.getLastShift(current_column_number);
			double masses_sum;
			int realignment_length;
			if (shift > ACCURACY)
			{
				if (current_line_number + 1 < this.location.getLength())
				{
					realignment_length = 2;
					masses_sum = AminoAcid.get(this.location.get(current_line_number)).getMass() + AminoAcid.get(this.location.get(current_line_number + 1)).getMass();
					while (current_line_number + realignment_length + 1 <= this.location.getLength() && Math.abs(shift - masses_sum) > ACCURACY)
					{
						masses_sum += AminoAcid.get(this.location.get(current_line_number + realignment_length)).getMass();
						++realignment_length;
					}
					if (current_line_number + realignment_length <= this.location.getLength() && Math.abs(shift - masses_sum) <= ACCURACY)
					{
						++this.nb_peaks_in_common;
						masses_sum = this.transformed_spectrum.getColumnMass(this.transformed_spectrum.getColumnCount() - 1);
						for (;realignment_length > 1; --realignment_length)
						{
							masses_sum -= AminoAcid.get(this.location.get(current_line_number - 1 + realignment_length)).getMass();
							this.nb_peaks_in_common += this.transformed_spectrum.massCount(masses_sum);
						}
					}
				}
			}
			int counter = 2;
			while (previous_column_number > 0)
			{
				if (this.transformed_spectrum.isDoubleColumn(current_column_number))
				{
					this.nb_peaks_in_common += 2;
				}
				else
				{
					++this.nb_peaks_in_common;
				}
				
				shift = this.transformed_spectrum.getShift(previous_column_number ,current_column_number, AminoAcid.get(this.location.get(current_line_number - 1)));			
				if (shift > ACCURACY)
				{
					masses_sum = 0.0;
					for (int j = current_line_number - 1; j > previous_line_number; --j)
					{
						masses_sum += AminoAcid.get(this.location.get(j - 1)).getMass();
					}
					this.nb_peaks_in_common += this.transformed_spectrum.massCount(this.transformed_spectrum.getColumnMass(current_column_number) - AminoAcid.get(this.location.get(current_line_number - 1)).getMass());
					if (Math.abs(shift - masses_sum) <= ACCURACY)
					{
						masses_sum = this.transformed_spectrum.getColumnMass(previous_column_number);
						
						for (realignment_length = 1 ;realignment_length < current_line_number - previous_line_number - 1; ++realignment_length)
						{
							masses_sum += AminoAcid.get(this.location.get(previous_line_number - 1 + realignment_length)).getMass();
							this.nb_peaks_in_common += this.transformed_spectrum.massCount(masses_sum);
						}
					}
				}
					
				current_line_number = previous_line_number;
				current_column_number = previous_column_number;
				if (previous_column_number > 0)
				{
					previous_line_number = this.scenarios.get(counter).getPreviousLineNumber();
					previous_column_number = this.scenarios.get(counter).getPreviousColumnNumber();
					++counter;
				}
			}
			if (this.transformed_spectrum.isDoubleColumn(current_column_number))
			{
				this.nb_peaks_in_common += 2;
			}
			else
			{
				++this.nb_peaks_in_common;
			}
			
			shift = this.transformed_spectrum.getShift(previous_column_number ,current_column_number, AminoAcid.get(this.location.get(current_line_number - 1)));			
			if (shift > ACCURACY)
			{
				masses_sum = 0.0;
				for (int j = current_line_number - 1; j > previous_line_number; --j)
				{
					masses_sum += AminoAcid.get(this.location.get(j - 1)).getMass();
				}
				this.nb_peaks_in_common += this.transformed_spectrum.massCount(this.transformed_spectrum.getColumnMass(current_column_number) - AminoAcid.get(this.location.get(current_line_number - 1)).getMass());
				if (Math.abs(shift - masses_sum) <= ACCURACY)
				{
					++this.nb_peaks_in_common;
					masses_sum = this.transformed_spectrum.getColumnMass(previous_column_number);
					
					for (realignment_length = 1 ;realignment_length < current_line_number - previous_line_number - 1; ++realignment_length)
					{
						masses_sum += AminoAcid.get(this.location.get(previous_line_number - 1 + realignment_length)).getMass();
						this.nb_peaks_in_common += this.transformed_spectrum.massCount(masses_sum);
					}
				}
			}
			else
			{
				++this.nb_peaks_in_common;
			}
		}
	}

	@Override
	public int compareTo(Alignment alignment) {
		if (this.score > alignment.score)
		{
			return -1;
		}
		else if (this.score < alignment.score)
		{
			return 1;
		}
		else
		{
			this.computeNbPeaksInCommon();
			alignment.computeNbPeaksInCommon();
			if (this.nb_peaks_in_common > alignment.nb_peaks_in_common)
			{
				return -1;
			}
			else if (this.nb_peaks_in_common < alignment.nb_peaks_in_common)
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
	}
	
	public static void reset()
	{
		MIN_SCENARIO_SCORE = Parameters.MIN_SCENARIO_SCORE();
		ACCURACY = Parameters.ACCURACY();
	}
}
