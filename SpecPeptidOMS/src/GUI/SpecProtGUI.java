package GUI;


import java.awt.Color;
import java.awt.Cursor;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import javax.swing.DefaultComboBoxModel;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JLayeredPane;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.WindowConstants;
import javax.swing.border.EtchedBorder;
import javax.swing.filechooser.FileNameExtensionFilter;
import javax.swing.text.DefaultCaret;

import algorithms.SpecGlobExecutor;
import constantes.Parameters;
import exceptions.ThreadException;
import loaders.Proteins;
import loaders.Spectra;
import main.Main;

/**
 * This class displays a graphical interface to the user. Once the user has clicked on the "launch button",
 *  this class reads the parameters, modifies the values by default in the Parameter class, and eventually launch the application.
 *  If one parameter is incorrect, a message is sent to the user in the LOG area, and the application waits 
 *  until a new click on the launch button.
 * 
 * @author Dominique Tessier, Emile Benoist
 *
 */


public class SpecProtGUI {
	
	/**
	 * this frame contains the graphical interface
	 */
	public JFrame frmSpecglobtool;
	/**
	 * this JTextArea displays all the messages sent to the user (information or errors)
	 */
	public static JTextArea LOG;
	
	private JTextField textFieldInputMGF;
	private JTextField textFieldInputProteins;
	private JTextField textFieldOutput;
	private Color colorBIBS1 = new Color(102, 193, 191);
	private Color colorBIBS2 = new Color(0, 140, 142);
	private Color colorBIBS3 = new Color(39, 86, 98);
	private Color backGround = new Color(242,240,231);
	private File spectraFile = new File("spectra.mgf");
	private String dataType = "MGF";
	private String spectraPath;
	
	/**
	 * Creates the interface displayed to the user with callbacks
	 * 
	 */
	
	public SpecProtGUI() {
		initialize();

	}
	
	/**
	 * Initializes the contents of the frame with the following 3 Panes, the LOG area and the launch button
	 * 	<p> layeredPane manages the spectra and proteins files selection;
	 * 	<p> layeredPane1 manages the result file selection
	 * 	<p> layeredPane2 manages the other parameters
	 * 
	 */
	
	private void initialize() {

		final String fontType = "Tahoma";
		final String baseFolderPath = System.getProperty("user.home") + "\\Documents";

		ImageIcon img = new ImageIcon(getClass().getResource("/ressources/IconSpecGlobTool.png"));
		frmSpecglobtool = new JFrame();
		frmSpecglobtool.setResizable(false);
		frmSpecglobtool.setTitle("SpecPeptidOMS\r\n");
		frmSpecglobtool.setBounds(100, 100, 733, 610);
		frmSpecglobtool.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
		frmSpecglobtool.getContentPane().setLayout(null);
		frmSpecglobtool.setIconImage(img.getImage());

		JScrollPane scrollPane = new JScrollPane();
		scrollPane.setBounds(20, 370, 680, 160);
		frmSpecglobtool.getContentPane().add(scrollPane);

		LOG = new JTextArea();
		LOG.setEditable(false);
		DefaultCaret caret = (DefaultCaret) LOG.getCaret();
		caret.setUpdatePolicy(DefaultCaret.ALWAYS_UPDATE);
		scrollPane.setViewportView(LOG);

		JLabel lblNewLabel = new JLabel("Welcome to SpecPeptidOMS* !");
		lblNewLabel.setForeground(colorBIBS3);
		lblNewLabel.setFont(new Font("Tahoma", Font.BOLD, 16));
		lblNewLabel.setHorizontalAlignment(SwingConstants.CENTER);
		lblNewLabel.setBounds(180, 10, 384, 30);
		frmSpecglobtool.getContentPane().add(lblNewLabel);
		

		JLabel logoINRAE = new JLabel(new ImageIcon(getClass().getResource("/ressources/LogoINRAE.png")));
		logoINRAE.setBounds(5, 5, 100, 40);
		frmSpecglobtool.getContentPane().add(logoINRAE);

		JLabel logoUNIV = new JLabel(new ImageIcon(getClass().getResource("/ressources/LogoUNIV.png")));
		logoUNIV.setBounds(102, 5, 119, 40);
		frmSpecglobtool.getContentPane().add(logoUNIV);

		JLabel logoBIBIS = new JLabel(new ImageIcon(getClass().getResource("/ressources/LogoBIBS.png")));
		logoBIBIS.setBounds(504, 5, 180, 40);
		logoBIBIS.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		
		frmSpecglobtool.getContentPane().add(logoBIBIS);


		// INPUT Pane
		JLayeredPane layeredPane = new JLayeredPane();
		layeredPane.setForeground(Color.WHITE);
		layeredPane.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
		layeredPane.setBackground(colorBIBS2);
		layeredPane.setBounds(10, 47, 310, 170);
		frmSpecglobtool.getContentPane().add(layeredPane);
		frmSpecglobtool.getContentPane().setBackground(backGround);
		
		
		JLabel labelINPUT = new JLabel("INPUT");
		labelINPUT.setBounds(10, 10, 205, 19);
		labelINPUT.setForeground(colorBIBS3);
		layeredPane.add(labelINPUT);
		labelINPUT.setHorizontalAlignment(SwingConstants.CENTER);
		labelINPUT.setFont(new Font("Tahoma", Font.BOLD, 14));

		JLabel labelSptFile = new JLabel("Spectra file (.mgf) (.mzML)");
		labelSptFile.setFont(new Font(fontType, Font.PLAIN, 12));
		labelSptFile.setBounds(10, 32, 191, 23);
		layeredPane.add(labelSptFile);
		
		textFieldInputMGF = new JTextField();
		textFieldInputMGF.setFont(new Font(fontType, Font.PLAIN, 12));
		textFieldInputMGF.setText("spectra.mgf");
		textFieldInputMGF.setBounds(10, 54, 261, 23);
		layeredPane.add(textFieldInputMGF);
		textFieldInputMGF.setColumns(10);

		// The button to select the input spectra file
		JButton buttonInputMGF = new JButton("Folder");
		buttonInputMGF.setBounds(281, 54, 23, 23);
		
		buttonInputMGF.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFrame frame = new JFrame("Load spectra scan data");
				JFileChooser fc = new JFileChooser(new File(baseFolderPath));
				fc.setDialogTitle("Load spectra datafile (MGF or MZML)");
				fc.addChoosableFileFilter(new FileNameExtensionFilter("MGF Files", "mgf"));
				// fc.addChoosableFileFilter(new FileNameExtensionFilter("MZXML Files",
				// "mzxml"));
				fc.addChoosableFileFilter(new FileNameExtensionFilter("MZML Files", "mzML"));
				fc.setAcceptAllFileFilterUsed(false);

				int result = fc.showOpenDialog(frame);

				if (result == JFileChooser.APPROVE_OPTION) {
					File selectedFile = fc.getSelectedFile();
					String selectedFilePath = selectedFile.getAbsolutePath();
	
					spectraFile = selectedFile;
					textFieldInputMGF.setText(selectedFilePath);
					// if the input file is a MGF file :
					if (selectedFile.getAbsolutePath().endsWith(".mgf")) {
						dataType = "MGF";
						// if the input file is a MZML file :
					} else if (selectedFile.getAbsolutePath().endsWith(".mzML")) {
						dataType = "MZML";

					}
				}

			}
		});
		layeredPane.add(buttonInputMGF);

	

		JLabel labelProteinFile = new JLabel("Proteins file (.fasta)");
		labelProteinFile.setFont(new Font(fontType, Font.PLAIN, 12));
		labelProteinFile.setBounds(10, 92, 191, 23);
		layeredPane.add(labelProteinFile);
				
		textFieldInputProteins = new JTextField();
		textFieldInputProteins.setFont(new Font(fontType, Font.PLAIN, 12));
		textFieldInputProteins.setText("proteins.fasta");
		textFieldInputProteins.setBounds(10, 114, 261, 23);
		layeredPane.add(textFieldInputProteins);
		textFieldInputProteins.setColumns(10);

		// The button to select the proteins file
		JButton buttonProteins = new JButton("Folder");
		buttonProteins.setBounds(281, 114, 23, 23);
		buttonProteins.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				JFrame frame = new JFrame("Select proteins file (fasta format)");
				JFileChooser fc = new JFileChooser(new File(baseFolderPath));
				fc.setDialogTitle("Load CSV info file");
				fc.addChoosableFileFilter(new FileNameExtensionFilter("Fasta Files", "fasta"));
				fc.setAcceptAllFileFilterUsed(false);

				int result = fc.showOpenDialog(frame);
				if (result == JFileChooser.APPROVE_OPTION) {
					File selectedFile = fc.getSelectedFile();
					String selectedFilePath = selectedFile.getAbsolutePath();
					textFieldInputProteins.setText(selectedFilePath);

				}
			}
		});
		layeredPane.add(buttonProteins);
		

		// OUTPUT PANE		

		JLayeredPane layeredPane1 = new JLayeredPane();
		layeredPane1.setForeground(Color.WHITE);
		layeredPane1.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
		layeredPane1.setBackground(colorBIBS2);
		layeredPane1.setBounds(10, 240, 310, 110);
		frmSpecglobtool.getContentPane().add(layeredPane1);
		
		JLabel labelOutput = new JLabel("OUTPUT");
		labelOutput.setHorizontalAlignment(SwingConstants.CENTER);
		labelOutput.setFont(new Font("Tahoma", Font.BOLD, 14));
		labelOutput.setForeground(colorBIBS3);
		labelOutput.setBounds(10, 10, 205, 19);
		layeredPane1.add(labelOutput);


		JLabel labelAlignFile = new JLabel("Alignment result file (.csv)\r\n");
		labelAlignFile.setFont(new Font(fontType, Font.PLAIN, 12));
		labelAlignFile.setBounds(10, 32, 191, 23);
		layeredPane1.add(labelAlignFile);
		
		
		textFieldOutput = new JTextField();
		textFieldOutput.setFont(new Font(fontType, Font.PLAIN, 12));
		textFieldOutput.setText("output.csv");
		textFieldOutput.setColumns(10);
		textFieldOutput.setBounds(10, 62, 261, 23);
		layeredPane1.add(textFieldOutput);
		
		
		JButton buttonOutputCSV = new JButton("Folder");
		buttonOutputCSV.setBounds(281, 62, 23, 23);
		layeredPane1.add(buttonOutputCSV);
		buttonOutputCSV.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				
				JFrame frame = new JFrame("Select CSV output file name");
				JFileChooser fc = new JFileChooser(new File(baseFolderPath));
				fc.setDialogTitle("Select CSV output file name");
				fc.addChoosableFileFilter(new FileNameExtensionFilter("CSV Files", "csv"));
				fc.setAcceptAllFileFilterUsed(false);

				int result = fc.showSaveDialog(frame);
				if (result == JFileChooser.APPROVE_OPTION) {
					File selectedFile = fc.getSelectedFile();
					String selectedFilePath = selectedFile.getAbsolutePath();
					textFieldOutput.setText(selectedFilePath);

				}

			}
		}); 
		
		// PARAMETERS PANE

		JLayeredPane layeredPane2 = new JLayeredPane();
		layeredPane2.setForeground(Color.WHITE);
		layeredPane2.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
		layeredPane2.setBackground(colorBIBS2);
		layeredPane2.setBounds(330, 47, 350, 210);
		frmSpecglobtool.getContentPane().add(layeredPane2);

		JLabel labelParam = new JLabel("PARAMETERS");
		labelParam.setHorizontalAlignment(SwingConstants.CENTER);
		labelParam.setFont(new Font("Tahoma", Font.BOLD, 14));
		labelParam.setForeground(colorBIBS3);
		labelParam.setBounds(10, 10, 350, 19);
		layeredPane2.add(labelParam);
		
		final JComboBox<Integer> comboBoxMultiprocess = new JComboBox<>();
		comboBoxMultiprocess.setFont(new Font(fontType, Font.PLAIN, 13));
		comboBoxMultiprocess.setBackground(backGround);
		comboBoxMultiprocess.setEnabled(false);
		comboBoxMultiprocess.setModel(new DefaultComboBoxModel<>(new Integer[] { 2, 4, 8, 10, 12 }));
		comboBoxMultiprocess.setBounds(260, 39, 48, 23);
		layeredPane2.add(comboBoxMultiprocess);

		JCheckBox checkBoxMultiprocess = new JCheckBox("Multiprocessing");
		checkBoxMultiprocess.setFont(new Font(fontType, Font.PLAIN, 12));
		checkBoxMultiprocess.setBackground(backGround);
		checkBoxMultiprocess.setBounds(20, 39, 120, 23);
		checkBoxMultiprocess.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {

				comboBoxMultiprocess.setEnabled(!comboBoxMultiprocess.isEnabled());

			}
		});
		layeredPane2.add(checkBoxMultiprocess);

		JLabel labelThreads = new JLabel("Number of threads :");
		labelThreads.setFont(new Font(fontType, Font.PLAIN, 12));
		labelThreads.setBounds(140, 39, 110, 23);
		layeredPane2.add(labelThreads);
		
		JLabel labelPrecision = new JLabel("Mass accuracy on MS/MS fragments :");
		labelPrecision.setFont(new Font(fontType, Font.PLAIN, 12));
		labelPrecision.setBounds(20, 70, 272, 23);
		layeredPane2.add(labelPrecision);

		final JTextField textFieldPRECISION = new JTextField();
		textFieldPRECISION.setText("0.02");
		textFieldPRECISION.setFont(new Font(fontType, Font.PLAIN, 12));
		textFieldPRECISION.setColumns(10);
		textFieldPRECISION.setBounds(260, 70, 48, 23);
		layeredPane2.add(textFieldPRECISION);

		JLabel labelDa = new JLabel("Da");
		labelDa.setFont(new Font(fontType, Font.PLAIN, 12));
		labelDa.setBounds(320, 70, 21, 23);
		layeredPane2.add(labelDa);
		
		JLabel labelFilter = new JLabel("Number of most intense peaks filtered :");
		labelFilter.setFont(new Font(fontType, Font.PLAIN, 12));
		labelFilter.setBounds(20, 100, 210, 23);
		layeredPane2.add(labelFilter);
		
		final JTextField textFieldValueFilter = new JTextField();
		textFieldValueFilter.setText("60");
		textFieldValueFilter.setFont(new Font(fontType, Font.PLAIN, 12));
		textFieldValueFilter.setBounds(260, 100, 48, 23);
		layeredPane2.add(textFieldValueFilter);
		textFieldValueFilter.setColumns(10);
		
		JLabel labelMinScore = new JLabel("Minimal score :");
		labelMinScore.setFont(new Font(fontType, Font.PLAIN, 12));
		labelMinScore.setBounds(20, 130, 153, 23);
		layeredPane2.add(labelMinScore);
		
		final JTextField textFieldMinScore = new JTextField();
		textFieldMinScore.setText("20");
		textFieldMinScore.setFont(new Font(fontType, Font.PLAIN, 12));
		textFieldMinScore.setBounds(260, 130, 48, 23);
		layeredPane2.add(textFieldMinScore);
		textFieldMinScore.setColumns(10);
		
		JLabel labelFixedModifs = new JLabel("Fixed Modifications :");
		labelFixedModifs.setFont(new Font(fontType, Font.PLAIN, 12));
		labelFixedModifs.setBounds(20, 160, 153, 23);
		layeredPane2.add(labelFixedModifs);

		final JTextField textFieldFixedModifs = new JTextField();
		textFieldFixedModifs.setText("57.04@C");
		textFieldFixedModifs.setFont(new Font(fontType, Font.PLAIN, 12));
		textFieldFixedModifs.setBounds(200, 160, 120, 23);
		layeredPane2.add(textFieldFixedModifs);
		textFieldFixedModifs.setColumns(10);
		
		/**
		 * Launch alignment button. 
		 */
		
		JLabel infoParameters = new JLabel("*Additional parameters can be set in the execution_parameters.ini file");
		infoParameters.setFont(new Font(fontType, Font.ITALIC, 10));
		infoParameters.setBounds(10,185,340,20);
		layeredPane2.add(infoParameters);
		
		final JButton buttonLaunch = new JButton("Launch Alignments !");		
		buttonLaunch.setBackground(colorBIBS3);
		buttonLaunch.setFont(new Font("Tahoma", Font.BOLD, 12));
		buttonLaunch.setForeground(Color.WHITE);
		buttonLaunch.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));

		buttonLaunch.setBounds(420, 270, 176, 60);
		frmSpecglobtool.getContentPane().add(buttonLaunch);

		
		buttonLaunch.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				if (buttonLaunch.isEnabled()) {
					buttonLaunch.setEnabled(false);
					buttonLaunch.setBackground(backGround);
					String absolutePath = textFieldInputMGF.getText();
					absolutePath= absolutePath.replace("/", "\\");
					int folderIndex = absolutePath.lastIndexOf("\\");
					if (folderIndex != -1) {
						Parameters.setSpectra_folder_path(absolutePath.substring(0, folderIndex));
						Parameters.setSpectra_file(absolutePath.substring(folderIndex+1));				
					}
					else
					{LOG.append("The spectra file name is not correct\n");
					buttonLaunch.setEnabled(true);
					return;
					}
					
					absolutePath = textFieldInputProteins.getText();
					folderIndex = absolutePath.lastIndexOf("\\");
					if (folderIndex != -1) {
						Parameters.setProteins_folder_path(absolutePath.substring(0, folderIndex));
						Parameters.setProteins_file(absolutePath.substring(folderIndex+1));				
					}
					else
					{LOG.append("The proteins file name is not correct\n");
					buttonLaunch.setEnabled(true);
					return;
					}

					absolutePath = textFieldOutput.getText();
					folderIndex = absolutePath.lastIndexOf("\\");
					if (folderIndex != -1) {
						Parameters.setResults_folder_path(absolutePath.substring(0, folderIndex));
						Parameters.setResults_file(absolutePath.substring(folderIndex+1));				
					}
					else
					{LOG.append("The result file name is not correct\n");
					buttonLaunch.setEnabled(true);
					return;
					}					

					LOG.append("Selected spectra file: " + spectraFile.getAbsolutePath() + "\r\n");
					LOG.append("Selected csv input file: " + textFieldInputProteins.getText() + "\n");
					LOG.append("Selected file to save: " + textFieldOutput.getText() + "\n");
					
					
					try {
					String precision = textFieldPRECISION.getText();
					double accuracy = Double.parseDouble(precision);
					Parameters.setAccuracy(accuracy);
					}
					catch (Exception err) {
						LOG.append("incorrect accuracy value");
						buttonLaunch.setEnabled(true);
						return;
					}

					if (comboBoxMultiprocess.isEnabled()) {
						Parameters.setNb_threads(Integer.parseInt(comboBoxMultiprocess.getSelectedItem().toString()));
					}
					else
						Parameters.setNb_threads(1);
					
					try {
						String valueFilter = textFieldValueFilter.getText();
						int nbPeaks = Integer.parseInt(valueFilter);
						Parameters.setNb_selected_peaks(nbPeaks);
						}
						catch (Exception err) {
							LOG.append("incorrect number of peaks filter value");
							buttonLaunch.setEnabled(true);
							return;
						}
					
					try {
						String minScore = textFieldMinScore.getText();
						int score = Integer.parseInt(minScore);
						Parameters.setMin_scenario_score(score);
						}
						catch (Exception err) {
							LOG.append("incorrect score filter value");
							buttonLaunch.setEnabled(true);
							return;
						}
					
					try {				
					String ptm = textFieldFixedModifs.getText();
					String ptmList[] = ptm.split(";");
					for (int i=0; i<ptmList.length;i++) {
						String onePTM = ptmList[i];
						int index = onePTM.indexOf('@');
						if (index == -1)
						{LOG.append("No fixed modification \n");
						break;
						}
						
						String valuePTM = onePTM.substring(0, index);
						double val =  Double.parseDouble(valuePTM);
						String aa = onePTM.substring(index+1);
						switch (aa) {
						case "A": Parameters.setAModif(val);
								  break;
						case "C": Parameters.setCModif(val);
						  		  break;
						case "D": Parameters.setDModif(val);
				  		  		  break;
						case "E": Parameters.setEModif(val);
				  		  		  break;
						case "F": Parameters.setFModif(val);
				  		  		  break;
						case "G": Parameters.setGModif(val);
				  		  		  break;
						case "H": Parameters.setHModif(val);
				  		  		  break;
						case "I": Parameters.setIModif(val);
				  		  		  break;
						case "K": Parameters.setKModif(val);
		  		  		  		  break;
						case "L": Parameters.setLModif(val);
		  		  		  		  break;
						case "M": Parameters.setMModif(val);
		  		  		  		  break;
						case "N": Parameters.setNModif(val);
		  		  		  		  break;
						case "O": Parameters.setOModif(val);
		  		  		  		  break;
						case "P": Parameters.setPModif(val);
		  		  		  		  break;
						case "Q": Parameters.setQModif(val);
		  		  		  		  break;
						case "R": Parameters.setRModif(val);
		  		  		          break;
						case "S": Parameters.setSModif(val);
		  		  		  		  break;
						case "T": Parameters.setTModif(val);
		  		  		   		  break;		  		  		  
						case "U": Parameters.setUModif(val);
		  		  		  		  break;
						case "V": Parameters.setVModif(val);
		  		  		  		  break;
						case "W": Parameters.setWModif(val);
		  		  		  		  break;
						case "Y": Parameters.setYModif(val);
		  		  		  		  break;
						default:
							LOG.append("Incorrect syntax for fixed modifications \n");
							throw new IllegalArgumentException("Unexpected value: " + aa);
						}
					}
					}
					catch (Exception err) {
						LOG.append("Fixed modifications syntax incorrect. Example: 57.02@C;27.99@S");
						buttonLaunch.setEnabled(true);
						return;
					}
					


					Thread t1 = new Thread(new Runnable() { public void run(){
						try
						{
							Spectra.load();
							if (Spectra.getNumberOfSpectra() > 0)
							{
								Proteins.load();								
								SpecGlobExecutor.initialise(Parameters.NB_THREADS());
								LOG.append("proteins loaded: "+Proteins.getNumberOfProteins()+"\n");
								SpecGlobExecutor.run();
								LOG.append("Job done.\n");
								Main.resetAll();
								buttonLaunch.setEnabled(true);
								buttonLaunch.setBackground(colorBIBS3);
							}
						}
						catch (FileNotFoundException err)
						{
							//err.printStackTrace();
							LOG.append(err.getMessage());
						}
						catch (IOException err)
						{
							//err.printStackTrace();
							LOG.append(err.getMessage());
						}
						catch (ThreadException err)
						{
							//System.err.println("\nA thread stopped due to an unknown Exception. The results file may be incomplete or even not exist. Please, send us the error message by ...\n");
							LOG.append(err.getMessage());
							//err.printRealExcpetion();
						}

						catch (Exception err)
						{
							//System.err.println("\nThe program stopped due to an unknown Exception. The results file may be incomplete or even not exist. Please, send us the error message by ...\n");
							//err.printStackTrace();
							LOG.append(err.getMessage());
						}

					}


				});
					t1.start();
		}			



		JLabel lblPleaseThanksTo = new JLabel("*How to cite us : https://doi.org/10.1101/2022.05.31.494131 ");
		lblPleaseThanksTo.setHorizontalAlignment(SwingConstants.LEFT);
		lblPleaseThanksTo.setBounds(10, 540, 429, 27);
		lblPleaseThanksTo.setCursor(Cursor.getPredefinedCursor(Cursor.HAND_CURSOR));
		frmSpecglobtool.getContentPane().add(lblPleaseThanksTo);
		lblPleaseThanksTo.setForeground(new Color(0, 0, 0));
		lblPleaseThanksTo.setFont(new Font("Tahoma", Font.ITALIC, 12));

			}
		});
	}
}
			
