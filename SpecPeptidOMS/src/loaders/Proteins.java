package loaders;

import constantes.Parameters;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import bioobjects.Protein;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import GUI.SpecProtGUI;

/*
 * cette classe permet le chargement de la banque de proteine au format 'fasta'. Les proteines chargees sont stockees ici
 */
public final class Proteins {

	private static Protein[] list; // la liste des proteines
	private static int max_length = 0; // la taille de la plus grande proteine stockee (sert pour la construction de(s) matrice(s) d'alignement)
	private static Pattern pattern = Pattern.compile("\\|.+\\|");
	private static Matcher matcher;
	
	public static void load() throws FileNotFoundException ,IOException
	{
		String path = Parameters.PROTEINS_FOLDER_PATH() + "/" + Parameters.PROTEINS_FILE();
		InputStream input_file;
		InputStreamReader input_stream;
		BufferedReader input_buffer;
		try
		{
			input_file = new FileInputStream(path);
			input_stream = new InputStreamReader(input_file);
			input_buffer = new BufferedReader(input_stream);
		}
		catch (FileNotFoundException e)
		{
			throw new FileNotFoundException("The proteins file \"" + path + "\" doesn't exists.");
		}
		if (path.endsWith(".fasta"))
		{
			try
			{
				Proteins.loadFromFasta(input_buffer);
			}
			catch (IOException e)
			{
				throw new IOException("The program stopped due to an IOException during the loading of the proteins file.");
			}
			finally
			{
				input_buffer.close();
			}
		}
	}
	
	public static int getMaxLength()
	{
		return Proteins.max_length;
	}
	
	public static int getNumberOfProteins()
	{
		return Proteins.list.length;
	}
	
	public static Protein[] get()
	{
		return Proteins.list;
	}
	
	public static Protein getProtein(int index) {
		return Proteins.list[index];
	}
	
	public static void GenerateCorrespondanceFile() throws IOException
	{
		FileWriter writer;
		try
		{
			writer = new FileWriter(Parameters.RESULTS_FOLDER_PATH() + "/" + "correspondanceFile.txt", false);
		}
		catch (IOException e)
		{
			throw new IOException("The program stopped due to an IOException during the creation of the proteins corresponding file.");
		}
		BufferedWriter buffer = new BufferedWriter(writer);
		PrintWriter output = new PrintWriter(buffer);
		for (Protein protein : Proteins.list)
		{
			output.write(protein.getName() + "\n");
		}
		output.close();
	}
	
	private static void loadFromFasta(BufferedReader input_buffer) throws IOException
	{
		ArrayList<Protein> list_temp = new ArrayList<Protein>();
		String line;
		String protein_sequence = "";
		String name = "";
		String protein_reference = "";
		Protein protein;
		while ((line = input_buffer.readLine()) != null)
		{
			if (!line.isEmpty())
			{
				if (line.charAt(0) == '>')
				{
					if (!protein_sequence.equals(""))
					{
						protein = new Protein(protein_sequence, name ,protein_reference);
						list_temp.add(protein);
						Proteins.max_length = Math.max(Proteins.max_length, protein_sequence.length());
						protein_sequence = "";
					}
					protein_reference = line;
					matcher = pattern.matcher(line);
					if (matcher.find())
					{
						name = matcher.group();
						name = name.substring(1, name.length() - 1);
					}
				}
				else
				{
					protein_sequence = protein_sequence.concat(line);
				}
			}
		}
		protein = new Protein(protein_sequence, name ,protein_reference);
		list_temp.add(protein);
		Proteins.max_length = Math.max(Proteins.max_length, protein_sequence.length());
		input_buffer.close();
		Proteins.list = new Protein[list_temp.size()];
		for (int i = 0; i < Proteins.list.length; ++i)
		{
			Proteins.list[i] = list_temp.get(i);
		}
		
		System.out.println(Proteins.getNumberOfProteins() + " proteins have been loaded.");
	}

	public static String foundPeptide(String peptide) {
		String answer = "";
		String protein_sequence;
		String current_positions;
		int current_position;
		int positions_counter = 0;
		int proteins_counter = 0;
		//for (int i = 0; i < list.length; ++i)
		int i = 0;
		protein_sequence = list[0].getRealSequence();
		current_position = protein_sequence.indexOf(peptide);
		if (current_position != -1)
		{
			++proteins_counter;
		}
		while (i < Proteins.list.length && proteins_counter <= 3)
		{
			current_positions = "";
			if (current_position != -1)
			{
				positions_counter = 1;
				while (current_position != -1 && positions_counter <= 3)
				{
					if (current_positions.length() == 0)
					{
						current_positions = String.format("[%d,%d]", current_position ,current_position + peptide.length() - 1);
					}
					else
					{
						current_positions = String.format("%s,[%d,%d]", current_positions ,current_position ,current_position + peptide.length() - 1);
					}
					current_position = protein_sequence.indexOf(peptide ,current_position + 1);
					if (current_position != -1)
					{
						++positions_counter;
					}
				}
				if (positions_counter > 0)
				{
					if (positions_counter > 1)
					{
						if (current_position != -1)
						{
							current_positions = String.format("(%s,...)", current_positions);
						}
						else
						{
							current_positions = String.format("(%s)", current_positions);
						}
					}
					if (answer.length() == 0)
					{
						answer = String.format("%s%s", Proteins.list[i].getName() ,current_positions);
						//answer = String.format("%s%s", i ,current_positions);
					}
					else
					{
						answer = String.format("%s,%s%s", answer, Proteins.list[i].getName() ,current_positions);
						//answer = String.format("%s,%s%s", answer, i ,current_positions);
					}
				}
			}
			++i;
			if (i < Proteins.list.length)
			{
				protein_sequence = list[i].getRealSequence();
				current_position = protein_sequence.indexOf(peptide);
				if (current_position != -1)
				{
					++proteins_counter;
				}
			}
		}
		if (i < Proteins.list.length)
		{
			answer = answer.concat(", ...");
		}
		return answer;
	}
	
	public static void reset()
	{
		max_length = 0;
	}
}
