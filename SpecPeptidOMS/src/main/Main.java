package main;

import constantes.*;
import datastructures.InterpretationsSaver;
import exceptions.ThreadException;
import loaders.Proteins;
import loaders.Spectra;
import algorithms.PeakFilters;
import algorithms.SpecGlobExecutor;
import algorithms.SpecGlobThread;
import bioobjects.TransformedSpectrum;
import components.Alignment;
import components.InterpretationBackup;
import components.LocationBackup;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.String;
import java.nio.IntBuffer;
import java.util.Arrays;
import java.awt.EventQueue;

import GUI.SpecProtGUI;

/*
 * java version 17
 * 
 * @author : Emile Benoist / Dominique Tessier
 * 
 * There is two arguments to pass with the execution command. The first one is the relative path of the file containing spectra. The second one must be the relative path of the file containing proteins.
 */

/**
 * The program can be used through a GUI or in command mode. 
 * 	In command_mode, a --c argument should be used 
 *  Otherwise the graphical interface is displayed by default 
 * 
 * @author 	BENOIST Emile, TESSIER Dominique
 * @version 1.0
 */

/*
 * A faire :
 *  gestion des erreurs
 *  mettre dans les AAPosition l'info pic doublé et revoir tout ça
 *  voir pour plusieurs niveaux d'assurances d'acide aminé
 *  que se passe t'il si le précurseur - 17 est plus petit que certains pics ?
 *  ne pas pouvoir repasser sur les mêmes acide aminés ?
 *  
 *  remettre l'intensité dans la liste des pics de chaque spectre une fois filtré ?
 *  améliorer la sortie console
 *  mettre une ligne dans le csv meme quand le spectre n'a aucun scenario retenu
 *  revoir completement la partie AllRealignment (erreur probable + modif 2aa)
 */

/*
 * problèmes :
 *  - plusieurs fois le même acide aminé au même endroit dans le spectre -> un seul gardé, peut faire apparaitre un realignement de 1 seul aa.
 *  Ce cas peut surement être réglé dans l'algo lui même (counter_rea == 1 ou un truc comme ça)
 *  
 *  - Revoir la filtrage des pics !
 */

public class Main {
	
	public static boolean commandMode = false;

	public static void main(String[] args) throws FileNotFoundException
	{
		// checking --c argument
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("--c")) {
				commandMode = true;
			}
		    if (args[i].equals("--paramFile")){
		    	i++;
		    	Parameters.setParamFile(args[i]);
		    }
		}
		
		// if we are in command mode, we take args and use them to create an instance of
		// SpecGlob
		if (commandMode) {
			System.out.println("=====PROCESS STARTED IN COMMAND MODE=====");
			try
			{
				double x = System.currentTimeMillis();
				Parameters.load();
				Spectra.load();
				if (Spectra.getNumberOfSpectra() > 0)
				{
					Proteins.load();
					Proteins.GenerateCorrespondanceFile();
					SpecGlobExecutor.initialise(Parameters.NB_THREADS());
					SpecGlobExecutor.run();
					System.out.println("Job done in " + (System.currentTimeMillis() - x)/1000.0 + " secondes. ");
				}			
			}
			
			catch (FileNotFoundException e)
			{
				e.printStackTrace();
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}
			catch (ThreadException e)
			{
				System.err.println("\nA thread stopped due to an unknown Exception. The results file may be incomplete or even not exist. Please, send us the error message by ...\n");
				e.printRealExcpetion();
			}
			
			catch (Exception e)
			{
				System.err.println("\nThe program stopped due to an unknown Exception. The results file may be incomplete or even not exist. Please, send us the error message by ...\n");
				e.printStackTrace();
			}
		}
		 else {
			 	double x = System.currentTimeMillis();
				EventQueue.invokeLater(new Runnable() {
					public void run() {
						try {						
							Parameters.load();
							//Proteins.GenerateCorrespondanceFile();
							SpecProtGUI window = new SpecProtGUI();
							window.frmSpecglobtool.setVisible(true);
						} catch (Exception e) {
							e.printStackTrace();
						}
					}
				});
			}

	}
	
	public static void resetAll()
	{
		PeakFilters.reset();
		SpecGlobExecutor.reset();
		SpecGlobThread.reset();
		TransformedSpectrum.reset();
		Alignment.reset();
		InterpretationBackup.reset();
		LocationBackup.reset();
		AminoAcid.reset();
		InterpretationsSaver.reset();
		Proteins.reset();
		Spectra.reset();
	}
}


	
