package comparators;

import java.util.Comparator;

import components.SimplePeak;

public class SimplePeakIntensityComparator implements Comparator<SimplePeak>{

	@Override
	public int compare(SimplePeak peak1, SimplePeak peak2) {
		return Double.compare(peak1.getIntensity(), peak2.getIntensity());
	}
}
