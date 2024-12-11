package components;

public final class Scenario {
	
	private int previous_line_number;
	private int previous_column_number;
	private AAPosition position;
	
	public Scenario() {}
	
	public Scenario(Scenario scenario)
	{
		this.previous_line_number = scenario.previous_line_number;
		this.previous_column_number = scenario.previous_column_number;
		this.position = scenario.position;
	}
	
	public Scenario(int previous_line_number ,int previous_column_number)
	{
		this.previous_line_number = previous_line_number;
		this.previous_column_number = previous_column_number;
	}

	public void update(Scenario scenario)
	{
		this.previous_line_number = scenario.previous_line_number;
		this.previous_column_number = scenario.previous_column_number;
		this.position = scenario.position;
	}

	public void update(int previous_line_number, int previous_column_number)
	{
		this.previous_line_number = previous_line_number;
		this.previous_column_number = previous_column_number;
	}
	
	public void update(int previous_line_number, int previous_column_number, AAPosition position)
	{
		this.previous_line_number = previous_line_number;
		this.previous_column_number = previous_column_number;
		this.position = position;
	}
	
	public int getPreviousLineNumber() {return this.previous_line_number;}
	
	public int getPreviousColumnNumber() {return this.previous_column_number;}
	
	public AAPosition getPosition() {return this.position;}
}
