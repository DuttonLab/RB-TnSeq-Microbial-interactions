# RB-TnSeq-Microbial-interactions


This pipeline has been designed to calculate and compare genome-wide gene fitness values across different growth conditions using RB-TnSeq (Wetmore *et al.,* mBio May 2015, 6 (3) e00306-15). RB-TnSeq relies on the utilization of a pooled library of barcoded transposon mutants. Each mutant has integrated a single transposon that includes a unique 20-nucleotide barcode; this barcode allows tracking of individual mutants in the population. When integrated within the 10-90% of an ORF, the transposon is expected to lead to gene disruption. RB-TnSeq libraries with high genome coverage contain multiple mutants associated with the same gene. Mutant fitness is calculated by comparing the abundance of each barcoded mutant after growth to its abundance in the T0 sample (inoculation). The mutant fitness values for mutants within a single gene are then incorporated to calculate a gene fitness value. The gene fitness value indicates the importance of the gene for fitness in the studied condition.
While the current pipeline (https://bitbucket.org/berkeleylab/feba/src/master/) developed by Wetmore et al., is designed to calculate individual gene fitness values in different conditions, nothing is implemented to statistically compare gene fitness values across conditions. Therefore, starting from the allpoolcounts file generated by BarSeqTest.pl from Wetmore *et al.,* 2015, we have developed our own pipeline to:
	-i) Generate normalized and comparable gene fitness values for each condition.
	-ii) Compare gene fitness values between two conditions. 

The pipeline is designed to implement 3 biological replicates of each condition. To obtain the normalized gene fitness value, we first calculate the normalized gene fitness values for each replicate individually (Script I: Script1_GeneFitness_Replicate.Rmd) and then average gene fitness across the three replicates (Script II: Script2_Averaging_Replicates.Rmd). Then, a third script allows you to compare gene fitness values between two conditions and identify the significant differences according to a chosen threshold of confidence and correction for multiple comparison method (Script III: Script3_2conditions_FitnessComparison.Rmd)


  ### Part I: Calculation of gene fitness values for each replicate.

Associated script: Script1_GeneFitness_Replicate.Rmd

For each replicate, a normalized gene fitness is obtained by 
i)	calculating fitness of each insertion mutant
ii)	averaging the fitness of insertion mutants associated with the same gene 
iii)	normalizing gene fitness based on gene location on the chromosome

Briefly, to calculate gene fitness values, the script uses raw counts of number of reads per barcode per condition. First, as fitness calculation is based on a log2 transformation, a pseudocount of 0.1 is added to all counts to avoid any 0 values. Then, any barcode that is associated with a non-coding region, inserted within the first 10% or the last 10% of a gene or associated with 3 or less counts in the T0 sample (inoculation sample) are filtered out. To enable comparison across all conditions, raw counts are then corrected in each condition independently using the average number of counts per neutral insertion mutants calculated thanks to chosen reference genes (ie genes whose deletion don’t have a fitness effect in any of the studied conditions). These corrected counts are then used to calculate mutant fitness and then gene fitness.
During data processing, different plots are generated to follow data modification and potentially spot some problems.



#### Inputs: 

  * **A .csv file** adapted, as described below, from the allpoolcounts.tab file generated by the Wetmore et al. perl script BarSeqTest.pl. This file contains the number of counts per barcode per condition. 
  
*Note*: as you run Script1__GeneFitness_Replicate.Rmd for each replicate individually, you need a .csv file for each replicate

IMPORTANT:
(i)the first 7 columns of this file are the first 7 columns of the allpoolcounts.tab file are the raw counts for each condition. They will be named by the condition.  
(ii)the 8th column is the T0 column named “T0”.
(iii)you will need to generate a .csv file for each replicate and maintain the same name for each condition across replicates.
See example: “Ex_Run_R1.csv”

  * The **gene.GC** file used to run the perl script TestBarSeq.pl



#### Parameters set up by user: 

In the script you have to specify a couple of parameters:
*  **org_locId** : “Num” if the locusId of the genes in your gene.GC file are numerics, “Char” if they are character. Note that when locusId’s are characters the script may take longer to run
*  **cdnb**: the number of conditions (including T0)
*  **scaffoldX** : the scaffold ID of the chromosome (plots are coded only for data on the chromosome even if the fitness values are calculated for insertions on chromosome and plasmids) 
*  **ref** : a vector containing the locusId of the reference genes chosen for raw counts normalization

#### Outputs:

An RData object containing: 
1. *All_data_Replicate*: table containing the final normalized gene fitness values for all the genes in each condition
2. *Data_norm_loc*: similar to All_data_Replicate but formatted differently
3. *Data_ref_corrected*: table containing the number of reads per barcode after filtering and normalization by the number of reads per neutral insertion mutant 
4. *Data_original*: input data
5. *genes.tab*: gene.GC file

#### Functions:

  *  **Data_prep_viz10KB()**: 
This function processes the number of reads per barcode so that we can visualize the number of reads per 10kb along the chromosome for each condition. 

 *  **Ref_counts_CorrandNorm ()**: 
For each condition, this function normalizes the number of reads per barcode according to the number of reads associated with insertion mutants of reference genes. This is performed after adding the pseudocount and filtering out the barcodes with low counts or inappropriate insertion location.

  *  **Gene_fitness ()**:
This function produces the unnormalized fitness value for each gene in each condition. It returns a list containing one table with insertion mutant fitness values, one table with unnormalized gene fitness values, and one table with gene fitness variance.
First, we calculate for each insertion mutant the strain fitness as: 
f_s=〖log〗_2 ((Corrected counts in Condition)/(Corrected counts in T0))
Then, gene fitness values are calculated as the average of the insertion mutant fitness values associated with that gene. Similarly, gene fitness variance is calculated as the variance of the insertion mutant fitness associated with that gene. 

  *  **Loc_smoothmed_norm()**:
This function performs a normalization based on the insertion location (as genes close to the replication fork may have higher count if cells are dividing, which would bias fitness calculation). This normalization is performed as described in Wetmore *et al.,* 2015. It is performed for each scaffold independently, and the transformation used depends on the length of the scaffold. For scaffolds with more than 251 genes: genes are normalized using the smoothed median in a window of 251 genes. For scaffolds with between 10 and 251 genes, we correct fitness values by the median fitness of the scaffold. If the scaffold has less than 10 genes, no normalization is performed.


  ### Part II: Averaging gene fitness values across replicates

Associated script: Script2_Averaging_Replicates.R

This script allows you to investigate correlation between your biological replicates and to average the normalized genes fitness values across the replicates.

Averaged genes fitness values are calculated as the weighted average of replicate gene fitness. 


#### Inputs: 
•	*If you have the same T0 for all your conditions*:

  *  **the .RData generated for each replicate in Script I** --> it should be 3 different .RData files
  
•	*If you have used different T0s for your conditions*:
You first need to run the .RData files generated in Script 1 for each replicates and sets of T0 in Multiple_T0s.R

  * **the .RData generated in Multiple_T0s.R**


#### Parameters: 

In this script you have to specify 1 parameter:

* **Org_locId**: “Num” if the locusId of the genes in your gene.GC file are numerics, “Char” if they are character.
* **multiT**: “TRUE” if you had different T0s for different conditions of the same replicate and had to use Multiple_T0s.R; “FALSE” if you had the same T0s for all the conditions in the replicate and directly load the .RData from Script1

#### Outputs:

An RData object containing (along with all the plots): 
1. *Final_gene_Fitness*: table containing the final averaged and normalized gene fitness values for all the genes in each condition
2. *Average_fitness*: table containing the inverse-variance weighted average gene fitness across replicates and associated standard deviation
3. *Correlation_table_fit*: table containing for each condition the Pearson correlation coefficient for all pairwise comparisons of replicates as well as the Spearman and Lin's coefficient
5. *AllReplicate*: single table containing the normalized gene fitness values for each replicate
6. *genes.tab*: gene.GC file

#### Functions:

* **Correlation_Rep()**: 
This function calculates the Pearson, Spearman and Lin’s correlation coefficients for all the replicate combinations for each condition (Replicate 1 versus Replicate 2, Replicate 1 versus Replicate 3 and Replicate 2 versus Replicate 3). This function produces a correlation table and automatically saves correlation graphs with Pearson R squared as a pdf document called “Correlation_plots_fit.pdf”

* **Weighted_average()**: 
This function averages each gene fitness value across replicates. More specifically, the gene fitness value is the inverse-variance weighted average of gene fitness across replicates. Similarly, for each gene, this function calculates the associated standard deviation. 

* **Mean_Correction()**:
Once average gene fitness has been calculated, each fitness value is corrected by the Mean of the replicates’ means of distribution. 

* **Mode_Centered()**:
This function concludes the calculation of gene fitness values. It centers the fitness distribution around the mode of the distribution, based on the assumption that most of the genes are associated with a neutral fitness.


  ### Part III: Comparing gene fitness values between two conditions

Associated script: Script3_2conditions_FitnessComparison.R

The focus of part III is the comparison of gene fitness values between two conditions. In the script, you determine a reference condition against which all the other conditions are going to be tested. 

The comparison of gene fitness values between two conditions relies on an unpaired two-sample Student test for samples with equal variance. Before running the two-sided Student test, we test for equal variance using a Fisher test. If variances are unequal between samples for a given gene, we do not perform the t-test.

The Fisher test calculates the ratio of the variance of the gene fitness in both conditions. The greatest variance has to be the numerator. If the ratio is EQUAL or GREATER than the F statistic, then variances are not equal. For the Fisher test, you can choose for a significant threshold to be alphaF=0.05 or alphaF=0.002.

In the case of equal variance, a t-statistic is calculated. In the script, you can choose to run correction for multiple comparison (parameter multi=1 or 0) and the value of alphaT, to screen for significant comparisons. This alphaT value is mostly used in the second function of the script (Category_Definition.R) which assigns labels to each gene comparison, whether it is significant (“Sig”), not significant (“Not_Sig”) or not tested (“Not_tested”).


#### Inputs: 

  * the final **.RData generated from Script II**

#### Parameters: 

In the script you have to specify these parameters:

* **Org_locId**: “Num” if the locusId of the genes in your gene.GC file are numerics, “Char” if they are character.
* **Condition1**: the condition against which you want to test all the other conditions
* **multi**: option to run correction for multiple comparison testing. 1=run correction for multiple comparison testing and calculated an adjusted-pvalue. The default method is “fdr” but can be changed in the function Comparison_test.R
* **AlphaF**: alpha you want to use for the Fisher test to reject H0 (H0=variance are equal)
* **AlphaT**: alpha you want to use for your student test to reject H0 (H0=gene fitness are equal). If multi=1, the screen is performed on the adjusted p-value

#### Outputs:

A list containing two elements: 
1. *Table_all*: the table containing, for each comparison, the gene fitness values, compared conditions, F-statistics, T-statistics and label for comparison status (“Sig”, “Not_sig”, “Not_tested”)
2. *List_plots*: list containing all the comparison plots
	

#### Functions:

* **Comparison_test()**: 
Performs for each comparison the Fisher test and the Student test. While producing the final table, it also produces a scatter plot for each comparison highlighting the genes with significantly different fitness values. 

* **Category_definition()**: 
This function is called within the Comparison_test() function. Based on the chosen alpha and the t-statistic calculated by the student test, this function assigns a category label to each test. “Sig”: the gene fitness values between the two conditions are significantly different based on your alpha choice; “Not_Sig”: the gene fitness values between the two conditions are not significantly different based on your alpha choice; “Not_tested”: the student test was not performed for that gene because the variances were unequal according to the Fisher test.


   ### Example: *E. coli* RB-TnSeq 5 conditions + T0 (Triplicate) from Piercer *et al*., 2020
   
Example resources are provided to run the pipeline using the example data from Pierce *et al*., 2020.

In this example, we use triplicates of 5 RB-TnSeq experiments (Condition 1 to 5) sharing the same T0 and performed in triplicate. They correspond to the *E. coli* K12 RB-TnSeq library (Wetmore *et al*., 2015) grown in 5 different conditions.

The pipeline is divided into 3 main scripts. We provide an associated .html file for each script:
- Script1_ GeneFitness_Replicate_Example.html
- Script2_Averaging_Replicates_Example.html
- Script3_2conditions_Fitness_Comparison.html

To open and visualize these .html files, please use the solution provided by htmlpreview.github.com:
	append https://htmlpreview.github.io/? to the URL of the .html file


**Script 1** (Script1_GeneFitness_Replicate.Rmd): Calculates gene fitness for of one replicate of a set of RB-TnSeq experiments of the same library and relying on the same T0 sample. 
Therefore, this script has to be run independently for each replicate of RB-TnSeq experiment using the same T0.
Here, the file Script1_ GeneFitness_Replicate_Example.html illustrates the run for Replicate 1.

To run that example, you need the following files:
- gene.GC
- Ex_Run_R1.csv

We also provide the files Ex_Run_R2.csv and Ex_Run_R3.csv to run that script for replicates 2 and 3 of the example.


**Script 2** (Script2_Averaging_Replicates.Rmd): Averages gene fitness values across replicates.

To run this script, you need the output of Script 1 for each replicate of the experiment. 
Then the script combines all the replicates of the 5 example conditions. 

You should be able to produce the following files from Script 1, that you can then use to run Script 2: 
- Ex_Run_R1.RData
- Ex_Run_R2.RData
- Ex_Run_R3.RData


**Script 3** (Script3_2conditions_Fitness_Comparison.Rmd): Performs the gene fitness comparisons between a given reference condition and the other conditions and identifies interaction fitness. 

That Script uses the output of Script2: All_Fitness_Values_Exemple.Rdata

Here, we used a run that compared the final gene fitness of Conditions 2 to 5 to gene fitness of reference condition, Condition 1 

**Provided files to reproduce and compare this example:**

*Initial data for Script1*:
- Replicate 1: Ex_Run_R1.csv
- Replicate 2: Ex_Run_R2.csv
- Replicate 3: Ex_Run_R3.csv
- gene.GC

*Output Script 1 and Input Script 2*:
- Replicate 1: Ex_Run_R1.RData
- Replicate 2: Ex_Run_R2.RData
- Replicate 3: Ex_Run_R3.RData

*Output Script 2*:
- All_Fitness_Values_Exemple.Rdata (input Script 3)
- All_Fitness_Values_Exemple.csv

*Output Script 3*:
	- Comparison_Ecoli_versusAlone.RData
	- Comparison_Ecoli_versusAlone.csv
	






