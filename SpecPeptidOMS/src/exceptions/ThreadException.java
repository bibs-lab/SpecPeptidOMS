package exceptions;

public class ThreadException extends Exception{
	

	private static final long serialVersionUID = 1L;
	private Exception real_exception;
	
	public ThreadException()
	{
		super();
	}

	public ThreadException(String s)
	{
		super(s);
	}
	
	public ThreadException(Exception e)
	{
		super();
		this.real_exception = e;
	}
	
	public void printRealExcpetion()
	{
		if (this.real_exception == null)
		{
			this.printStackTrace();
		}
		else
		{
			this.real_exception.printStackTrace();
		}
	}
}
