package components;

import java.text.DecimalFormat;

public final class ComplementaryPeak implements Comparable<ComplementaryPeak> {
	
	private static final DecimalFormat Decimal_Format = new DecimalFormat("0.00");

	private double mass;
	private boolean is_double;
	private int native_peak;
	
	public ComplementaryPeak() {}
	
	/*
	public ComplementaryPeak(double mass ,boolean is_double ,int native_peak)
	{
		this.mass = mass;
		this.is_double = is_double;
		this.native_peak = native_peak;
	}
	*/
	
	public void overwrite(double mass, boolean is_double, int native_peak)
	{
		this.mass = mass;
		this.is_double = is_double;
		this.native_peak = native_peak;
	}
	
	public double get_mass() {return this.mass;}
	
	public boolean isDouble() {return this.is_double;}
	
	public int getNativePeak() {return this.native_peak;}
	
	public String toString() {return Decimal_Format.format(this.mass);}

	@Override
	public int compareTo(ComplementaryPeak peak) {return Double.compare(this.mass, peak.mass);}
}
