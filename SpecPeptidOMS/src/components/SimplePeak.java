package components;

/*
 * Cet objet correspond simplement a un pic (intensite + masse).
 */
public final class SimplePeak
{

	private final double mass;
	private final double intensity;
	
	public SimplePeak(double mass ,double intensity) {
		this.mass = mass;
		this.intensity = intensity;
	}
	
	public double getMass()
	{
		return this.mass;
	}
	
	public double getIntensity()
	{
		return this.intensity;
	}
	
	public boolean isMoreIntense(SimplePeak peak) {return this.intensity >= peak.intensity;}
	
	public boolean isHeavier(SimplePeak peak) {return this.mass >= peak.mass;}
	
	public String toString()
	{
		return String.format("mass : %f, intensity : %f", this.mass ,this.intensity);
	}
}
