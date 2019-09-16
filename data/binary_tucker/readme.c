1. Nations (Nickel et al., 2011, dnations.mat): This is a 14 × 14 × 56 binary tensor consisting of 56 political relations of 14 countries between 1950 and 1965. The tensor entry indicates the presence or absence of a political action, such as “treaties”, “sends tourists to”, between the nations. Please set the diagonal elements Y(i,i,k) = 0 in the analysis.

I have pre-selected ~10 covariates for each nation. These covariates describe a few important country attributes, e.g. whether a nation is actively involved in medicine NGO, law NGO, or belongs to a catholic nation, etc. 

2. HCP (Wang et al., 2017a, Count_VSPLOT.mat): This is a 68 × 68 × 136 binary tensor consisting of structural connectivity patterns among 68 brain regions for 136 individuals from Human Connectome Project (HCP). All the individual images were preprocessed following a standard pipeline (Zhang et al., 2018), and the brain was parcellated to 68 regions-of-interest following the Desikan atlas (Desikan et al., 2006). The tensor entries encode the presence or absence of fiber connections between those 68 brain regions for each of the 136 individuals.

There are 573 covariates available for each individual. I have pre-selected ~20 important covariates. Please see the listed below. You may start with them first and then add a few new covariates of your interest. 
The full list of covariates can be found at:
https://wiki.humanconnectome.org/display/PublicData/HCP+Data+Dictionary+Public-+Updated+for+the+1200+Subject+Release?preview=/53444663/113377285/HCP_S1200_DataDictionary_April_20_2018.xlsx#HCPDataDictionaryPublic-Updatedforthe1200SubjectRelease-Excelversion:


Below is the dictinary for the important covariates:
Gender: binary, M or F
Age: categorical, 22-25, 26-30, 31-35, >35
MMSE_Score: numerical, Mini Mental Status Exam Total Score
PSQI_Score: numerical, Pittsburgh Sleep Questionnaire Total Score
PicSeq_Unadj: Measure of Episoidc Memory. NIH Toolbox Picture Sequence Memory Test: Unadjusted Scale Score, 
CardSort_Unadj: Measure of Executive Function/Cognitive Flexibility. NIH Toolbox Dimensional Change Card Sort Test: Unadjusted Scale Score,
Flanker_Unadj: Measure of Executive Function/Inhibition. NIH Toolbox Flanker Inhibitory Control and Attention Test: Unadjusted Scale Score
PMAT24_A_CR: Measure of Fluid Intelligence. Penn Progressive Matrices: Number of Correct Responses (PMAT24_A_CR)
PMAT24_A_RTCR: Measure of Fluid Intelligence. Penn Progressive Matrices: Median Reaction Time for Correct Responses (PMAT24_A_RTCR)
ReadEng_Unadj: Measure of Language/Reading Decoding. NIH Toolbox Oral Reading Recognition Test: Unadjusted Scale Score
PicVocab_Unadj: Measure of Language/Vocabulary Comprehension. NIH Toolbox Picture Vocabulary Test: Unadjusted Scale Score
ProcSpeed_Unadj: Measure of Processing Speed. NIH Toolbox Pattern Comparison Processing Speed Test: Unadjusted Scale Score 
VSPLOT_TC: Measure of Spatial Orientation. Variable Short Penn Line Orientation: Total Number Correct (VSPLOT_TC)
VSPLOT_CRTE: Measure of Spatial Orientation. Variable Short Penn Line Orientation: Median Reaction Time Divided by Expected Number of Clicks for Correct (VSPLOT_CRTE)
SCPT_FP: Measure of Sustained Attention. Short Penn Continuous Performance Test: Open Positives.
SCPT_FN: Measure of Sustained Attention. Short Penn Continuous Performance Test: Open Negatives
IWRD_TOT: Measure of Verbal Episodic Memory. Penn Word Memory Test:  Total Number of Correct Responses (IWRD_TOT)
IWRD_RTC: Measure of Verbal Episodic Memory. Penn Word Memory Test:  Median Reaction Time for Correct Responses (IWRD_RTC)
ListSort_Unadj: Measure of Working Memory. NIH Toolbox List Sorting Working Memory Test: Unadjusted Scale Score
ER40_CR: Measure of Emotion. Penn Emotion Recognition Test: Number of Correct Responses (ER40_CR)
ER40_CRT: Measure of Emotion. Penn Emotion Recognition Test: Correct Responses Median Response Time (ms) (ER40_CRT)

For more information on the covariates, please see:
https://wiki.humanconnectome.org/display/PublicData/HCP+Data+Dictionary+Public-+Updated+for+the+1200+Subject+Release?preview=/53444663/113377285/HCP_S1200_DataDictionary_April_20_2018.xlsx#HCPDataDictionaryPublic-Updatedforthe1200SubjectRelease-Excelversion:




