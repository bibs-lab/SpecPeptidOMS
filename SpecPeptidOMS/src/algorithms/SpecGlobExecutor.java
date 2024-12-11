package algorithms;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import GUI.SpecProtGUI;
import bioobjects.NativeSpectrum;
import bioobjects.TransformedSpectrum;
import constantes.Parameters;
import datastructures.InterpretationsSaver;
import exceptions.ThreadException;
import loaders.Spectra;
import main.Main;

/*
 * This class manages the creation and execution of threads.
 * @author 	BENOIST Emile, TESSIER Dominique
 */
public final class SpecGlobExecutor {
	
	private static int REAL_TIME_SAVE = Parameters.REAL_TIME_SAVE();
	
	private static Thread[] threads;
	private static Thread main_thread;
	private static Exception thread_exception;
	
	private static int next_spectrum = 0;
	private static int process_step = 1;
	private static int number_of_waiting_threads = 0;
	private static long startTime;
	private static long endTime;
	
	private static NativeSpectrum[] natives_spectra;
	private static TransformedSpectrum[] transformed_spectra;
	private static InterpretationsSaver[] interpretations_saver;
	
	/**
	 * Initialize as many SpecGlobThread as the parameter 
	 * @param nbThread: the number of SpecGlobThread used (at least one)
	 */
	public static void initialise(int nbThread) throws IOException
	{
		main_thread = Thread.currentThread();
		
		/*
		 * An empty results file is created here
		 */
		FileWriter writer;
		try
		{
			writer = new FileWriter(Parameters.RESULTS_FOLDER_PATH() + "/" + Parameters.RESULTS_FILE(), false);
		}
		catch (IOException e)
		{
			throw new IOException("The program stopped due to an IOException during the creation of the results file.");
		}
		BufferedWriter buffer = new BufferedWriter(writer);
		PrintWriter output = new PrintWriter(buffer);
		output.write("spectrum title;spectrum scan;spectrum id;peptide;protein positions;alignment;score;number of shared peaks;post-processed non-aligned mass;post-processed peptide;post-processed protein positions;post-processed alignment;post-processed score;post-processed number of common peaks\n");
		output.close();
		
		SpecGlobExecutor.threads = new Thread[nbThread];
		for (int threadNumber = 0; threadNumber < nbThread; ++threadNumber)
		{
			SpecGlobExecutor.threads[threadNumber] = new Thread(new SpecGlobThread(threadNumber));
		}
		
		SpecGlobExecutor.natives_spectra = Spectra.getNatives();
		SpecGlobExecutor.transformed_spectra = Spectra.getTransformed();
		
		/*
		 * We create as many interpretations_saver as the REAL_TIME_SAVE parameter.
		 */
		SpecGlobExecutor.interpretations_saver = new InterpretationsSaver[REAL_TIME_SAVE];
		for (int i = 0; i < REAL_TIME_SAVE; ++i)
		{
			SpecGlobExecutor.interpretations_saver[i] = new InterpretationsSaver();
		}
	}
	
	/**
	 * Run each of the SpecBlogThread.
	 */
	public static void run() throws IOException, ThreadException
	{
		startTime = System.currentTimeMillis();
		int min;
		for (int threadNumber = 0; threadNumber < SpecGlobExecutor.threads.length; ++threadNumber)
		{
			threads[threadNumber].start();
		}
		
	
		/*
		 * The main thread wait until REAL_TIME_SAVE spectra are treated by the SpecGlobThreads before performing a save for these spectra and sleep again. 
		 */
		while (true)
		{
			synchronized(SpecGlobExecutor.class)
			{
				try
				{
					SpecGlobExecutor.class.wait();
		        }
				catch (InterruptedException e)
				{
					for (int threadNumber = 0; threadNumber < SpecGlobExecutor.threads.length; ++threadNumber)
					{
						threads[threadNumber].interrupt();
					}
					throw new ThreadException(thread_exception);
		        }
			}
			
			int nbSpectra = REAL_TIME_SAVE * SpecGlobExecutor.process_step;
			nbSpectra = Math.min(REAL_TIME_SAVE * SpecGlobExecutor.process_step, SpecGlobExecutor.natives_spectra.length);
			System.out.println("save..."+nbSpectra);
			if (!Main.commandMode)
				SpecProtGUI.LOG.append(nbSpectra+" spectra processed\n");
			
			endTime = System.currentTimeMillis();
			double time = ((double) (endTime - startTime))/1000;
			System.out.println("** " + time + "s");
			
			startTime=endTime;
			
			try
			{
				FileWriter writer = new FileWriter(Parameters.RESULTS_FOLDER_PATH() + "/" + Parameters.RESULTS_FILE(), true);
				BufferedWriter buffer = new BufferedWriter(writer);
				PrintWriter output = new PrintWriter(buffer);
				min = Math.min(REAL_TIME_SAVE, SpecGlobExecutor.natives_spectra.length - (SpecGlobExecutor.process_step - 1) * REAL_TIME_SAVE);
				for (int i = 0; i < min; ++i)
				{
					SpecGlobExecutor.interpretations_saver[i].save(output);
				}
				output.close();
			}
			catch (IOException e)
			{
				for (int threadNumber = 0; threadNumber < SpecGlobExecutor.threads.length; ++threadNumber)
				{
					threads[threadNumber].interrupt();
				}
				throw new IOException("The program stopped due to an IOException during a results backup. The results file may be incomplete.");
			}

			SpecGlobExecutor.number_of_waiting_threads = 0;
			++SpecGlobExecutor.process_step;
			synchronized(SpecGlobExecutor.class)
			{
				SpecGlobExecutor.class.notifyAll();
			}
			if (SpecGlobExecutor.next_spectrum == SpecGlobExecutor.natives_spectra.length)
			{
				break;
			}
		}
	}
	
	/**
	 * Return the number of the next spectrum to be treated
	 * @param thread: the SpecGlobThread asking for the next spectrum
	 */
	public static synchronized void getNext(SpecGlobThread thread)
	{
		while (SpecGlobExecutor.next_spectrum >= SpecGlobExecutor.process_step * REAL_TIME_SAVE)
		{		
			try
			{
				++SpecGlobExecutor.number_of_waiting_threads;
				if (SpecGlobExecutor.number_of_waiting_threads == SpecGlobExecutor.threads.length)
				{
					SpecGlobExecutor.class.notifyAll();
				}
				SpecGlobExecutor.class.wait();
	        }
			catch (InterruptedException e)
			{
				/*
				 * Call updateSpectrum with null, null, null parameters allows to terminated the SpecGlobThread
				 */
	            thread.updateSpectrum(null, null, null);
	            return;
	        }
		}
		if (SpecGlobExecutor.next_spectrum < SpecGlobExecutor.natives_spectra.length) 
		{
			thread.updateSpectrum(SpecGlobExecutor.natives_spectra[SpecGlobExecutor.next_spectrum], SpecGlobExecutor.transformed_spectra[SpecGlobExecutor.next_spectrum], SpecGlobExecutor.interpretations_saver[SpecGlobExecutor.next_spectrum % SpecGlobExecutor.REAL_TIME_SAVE]);
			++SpecGlobExecutor.next_spectrum;
		}
		else
		{
			++SpecGlobExecutor.number_of_waiting_threads;
			if (SpecGlobExecutor.number_of_waiting_threads == SpecGlobExecutor.threads.length)
			{
				SpecGlobExecutor.class.notify();
			}
			thread.updateSpectrum(null, null, null);
		}
	}
	
	public static void Interrupt(Exception e)
	{
		thread_exception = e;
		main_thread.interrupt();
	}
	
	public static void reset()
	{
		next_spectrum = 0;
		process_step = 1;
		number_of_waiting_threads = 0;	
		REAL_TIME_SAVE = Parameters.REAL_TIME_SAVE();
	}
}
