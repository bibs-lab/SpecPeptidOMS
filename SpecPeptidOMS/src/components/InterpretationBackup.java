package components;

import constantes.Parameters;

public final class InterpretationBackup {
	
	private static int MIN_SCENARIO_SCORE = Parameters.MIN_SCENARIO_SCORE();

	private int line_number;
	private int column_number;
	private int score;
	
	public InterpretationBackup() {this.score = MIN_SCENARIO_SCORE;}
	
	public void reinit() {this.score = MIN_SCENARIO_SCORE;}
	
	public void tryImprove(int score ,int line_number ,int column_number)
	{
		if (score > this.score)
		{
			this.line_number = line_number;
			this.column_number = column_number;
			this.score = score;
		}
	}
	
	public int getLineNumber() {return this.line_number;}
	
	public int getColumnNumber() {return this.column_number;}
	
	public int getScore() {return this.score;}
	
	public static void reset()
	{
		MIN_SCENARIO_SCORE = Parameters.MIN_SCENARIO_SCORE();
	}
}
