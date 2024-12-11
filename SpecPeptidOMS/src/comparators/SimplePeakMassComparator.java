package comparators;

import java.util.Comparator;

import components.SimplePeak;

public class SimplePeakMassComparator implements Comparator<SimplePeak>{

	@Override
	public int compare(SimplePeak peak1, SimplePeak peak2) {
		return Double.compare(peak1.getMass(), peak2.getMass());
	}
}
