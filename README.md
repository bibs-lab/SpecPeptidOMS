# SpecPeptidOMS : an MS2 alignment tool that can be used in proteomics or in peptidomics

## Presentation

Welcome to the SpecPeptidOMS project !

SpecPeptidOMS is a Java program for interpreting MS2 mass spectra in peptidomics or proteomics. SpecPeptidOMS is capable of aligning several tens of thousands of experimental MS2 mass spectra to an entire proteome in just a few hours on a laptop.

## How to compile ?

Java VM is needed, java17 JDK or a more recent release is required to compile.<br>
 The main of SpecPeptidOMS is declared in the Main class.<br>
 Files in mzML format are loaded with the JmzReader parser (see pom.xml for dependencies)

We tested with success the compilation with java 17 JDK implemented by Oracle corporation (release-17.0.6)



## Use SpecPeptidOMS through its GUI Mode
A jar file can be downloaded directly from the target repository. Double-click on the jar file.<br><br>
A window appears, where you can select input and output files (spectra file, protein database file and result file). The most important parameters can be chosen via the GUI. Other parameters are stored in the *execution_parameter.ini* file that **must be in the same folder as the jar**. You can modify the parameters in execution.ini before the jar is running.<br><br>
The spectra file should be in the MGF or mzML format.<br><br>
The protein file should be in Fasta format.<br>

### Command Mode

The jar file can be launch in command mode (useful for pipeline implementations) with the option ```--c```.

Example :

	java -Xms16G -jar SpecPeptidOMS-1.0.0.jar --c
	
If a double click on the jar file generates a Java error "Java Exception has occured", please verify your java version used by default with the command (for windows): "ftype jarfile".

Parameters :

```
* --c : enable command mode
* --paramFile <paramFile name> can be used to change the parameter file by default (execution_parameters.ini)
*
## Configuration

The file **execution_parameters.ini** contains the parameters. This file must be placed in same folder as the executable jar file.

Parameters Description :


Files:
*  proteinsFolderPath: repository where the protein fasta file is
*  proteinsFile: name of the protein fasta file
*  spectraFolderPath: repository where the spectra file is 
*  spectraFile: name of the spectra file (.mzML or .mgf)
*  resultsFolderPath: repository where the result file is written 
*  resultsFile: name of the result file (.csv file)

Spectra selection:
*  nbPeaksMin (default 10): a spectrum is processed only if its number of peaks is >= nbPeaksMin
*  nbMinAAFound (default 10): a spectrum is processed only if it has at least nbMinAAFound labels

Spectra filter:
*  usedFilter (default mostIntense): currently the only filter implemented 
*  nbSelectedPeaks (default 60): number of peaks retained after filtering

General settings:
*  accuracy (default 0.02): fragment accuracy given in Daltons
*  nbThreads (default 10): number of threads once spectra and proteins are loaded

Algorithm parameters:
*  nbLssSMSaved (default 5): number of Located sub-sequence Spectrum Matches (LssSM) per spectrum
*  nbResultsReturned (default 1): number of interpretations per spectrum returned to the user
*  tolPeakMissing (default 4): maximum number of missing peaks accepted to score a realignment(in both rounds)
*  minScenarioScore (default 30): only interpretations with a score over this threshold are returned to the user
*  nbResultsAtOnce (default 500): number of spectra alignments written at a time
*  nbPeaksMax (default 500): Spectra should not exceed this number of peaks (can be adjusted with a small impact on memory usage)
*  surplus (default 5) : number of amino acids that extend the LssSM for the second alignment
*  maxRealignSizeFirstColumn (default 5): Maximum number of amino acids that can be shifted without penalty at the beginning of LssSM
*  filterLssSMOnScore (default 0.9): Filter LssSM when several are returned per spectrum after the preliminary alignment. Those whose score is below the best score multiply by this factor are excluded.

Choose algorithm:
*  shutDown3PeaksVersion=false
*  shutDownPeaksCleaning=false
*  shutDownNonAlignedMass=false
*  versionPreliminaryTreatment=1

Dynamic programming elementary scores for the preliminary alignment: Be aware, modification of the scores  is not recommended for non-experienced user.
*  certainlyFoundMain (default 10): Bonus for an aligned merged peak (original + complementary)
*  foundMain (default 7): Bonus for an aligned original/complementary peak
*  certainlyFoundWithShiftMain (default -6): Penalty for a shifted merged peak (original + complementary)
*  foundWithShiftMain (default -8): Penalty for an original/complementary peak shifted
*  notFoundMain (default -4): Penalty for missing amino acid.

Dynamic programming elementary scores for the second alignment round: Be aware, modification of the scores  is not recommended for non-experienced user.
*  certainlyFoundPost (default 10): Bonus for an aligned merged peak (original + complementary)
*  foundPost (default 7): Bonus for an aligned original/complementary peak
*  certainlyFoundWithShiftPost (default -6): Penalty for a shifted merged peak (original + complementary)
*  foundWithShiftPostn (default -8): Penalty for an original/complementary peak shifted
*  notFoundPost (default -4): Penalty for missing amino acid.

Fixed Modifications (some examples below)
* CModif= 57.021464
* AModif= 0.0
* DModif= 0.0
* NTERModif = 0.0
* CTERModif = 0.0


  ```

## Results

Results are returned under the CSV format with one or several lines per spectrum depending on the parameter <i>nbResultsReturned <i> . <br>



Detailed description of each column :

```
* The spectrum is identified by the first three columns
*   Title 
*   Scan number
*   Index in the file 
* 
* Five columns provide the best alignments after the two rounds of alignment
*   Peptide
*   List of at most 3 protein names containing the peptide, with their position in the sequence. The symbol ... indicates that the list exceeds 3 proteins, but has been truncated for readability
*   Alignment written with the same syntax we used in SpecGlobX. The syntax is provided with examples for reference below.
*   Score
*   Number of shared peaks between the peptide and the filtered spectrum
*
* Six additional columns provide a new alignment if the post-processing step improves the results
*   post-processed non-aligned mass
*   post-processed peptide
*   post-processed protein positions
*   post-processed alignment with the same syntax as above
*   post-processed score
*   post-processed number of shared peaks
```

SpecPeptidOMS uses the same syntax as SpecGlobX to express alignments as strings in the alignment columns. The aim is to summarize information about the alignment, providing a simplified fragmentation summary, highlighting stretches of detected (resp. unfound) amino acids in the alignment.

The alignment is done with the filtered experimental spectrum (not all the peaks are considered).

* When the two fragmentation peaks of an amino acid (b-ion and y-ion) are not used in the alignment, the amino acid is written between brackets
* When at least one of the fragmentation peaks of an amino acid is used in the alignment, but this alignment requires a shift, then the amino acid is written, preceded by the value of the mass shift in brackets
* When at least one the fragmentation peaks of an amino acid is used in the alignment and no shift is required for this alignment, then the amino acid is reported as such.
