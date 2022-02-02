# Install the IsoformSwitchAnalyzeR package through BiocManager
# I recommend always running this as the package is updated fairly frequently to fix issues.
if (!requireNamespace("devtools", quietly = TRUE)){ install.packages("devtools") }
devtools::install_github("kvittingseerup/IsoformSwitchAnalyzeR", build_vignettes = TRUE)


# Load IsoformSwitchAnalyzeR to the environment
library(IsoformSwitchAnalyzeR)
packageVersion('IsoformSwitchAnalyzeR')

### Import quantifications:
# Quantifications are generated following the differential expression analysis steps outlined in 
# the StringTie manual:
# HISAT2 -> read alignedment -> StringTie -> assembled transcripts -> StringTie --merge 
# -> merged transcripts -> StringTie -e -B -G merged_transcripts.gtf
# Read length is set to 150bp based on FastQC results

stringtieQuant <- importIsoformExpression(parentDir = "PATH_TO_PARENT_DIRECTORY", addIsofomIdAsColumn = TRUE, readLength = 150)

# Build the deisgn matrix for 3 conditions, with 2 samples each

sample_IDs = c("42MGBA_1", "42MGBA_2", "42MGBA_TMZres_1", "42MGBA_TMZres_2", "T98G_1", "T98G_2")
conditions = c("sensitive", "sensitive", "acquired_resistant", "acquired_resistant", "intrinsically_resistant", "intrinsically_resistant")
designMatrix <- data.frame(sample_IDs, conditions)
colnames(designMatrix) <- c("sampleID", "condition")
designMatrix



# Importing the data into R
# For exon annotation, we use the .gtf file generated after the StringTie --merge step
# For transcript fasta file, we use gffread with the merged.gtf as input to get the sequences for all transcripst
# Because we use annotation (gtf) from StringTie, we set fixStringTieAnnotationProblem = TRUE

mySwitchList <- importRdata(
  isoformCountMatrix   = stringtieQuant$counts,
  isoformRepExpression = stringtieQuant$abundance,
  designMatrix         = designMatrix,
  isoformExonAnnoation = "/PATH/stringtie_merged.gtf",
  isoformNtFasta       = "/PATH/all_transcripts_merged.fa",
  fixStringTieAnnotationProblem = TRUE,
  showProgress = TRUE
)

# We predict ORF regions in this step and overwrite any ORF information since StringTie --merge
# gtf does not contain CDS annotations.

mySwitchList <- analyzeORF(
  mySwitchList,
  genomeObject = NULL,
  minORFlength=100,
  orfMethod = "longest",
  cds = NULL,
  PTCDistance = 50,
  startCodons="ATG",
  stopCodons=c("TAA", "TAG", "TGA"),
  showProgress=TRUE,
  quiet=FALSE
)


# Conducting Analysis Part 1 to extract isoform switches and consequences.


mySwitchList <- isoformSwitchAnalysisPart1(
  switchAnalyzeRlist=mySwitchList,
  prepareForWebServers = FALSE,
  outputSequences = FALSE # keeps the function from outputting the fasta files from this example
)




#Adding Alternative Splicing and Intron Retention analysis to the external analysis run. 
#Also adding the following: CPAT, ORF, PFAM, signal peptides and IDR analyses before we run Part 2.

mySwitchList <- analyzeAlternativeSplicing(
  mySwitchList,
  onlySwitchingGenes=FALSE,
  alpha=0.05,
  dIFcutoff = 0.1,
  showProgress=TRUE,
  quiet=FALSE
)


mySwitchList <- analyzeIntronRetention(
  mySwitchList,
  onlySwitchingGenes=FALSE,
  alpha=0.05,
  dIFcutoff = 0.1,
  showProgress=TRUE,
  quiet=FALSE
)


# ADDITONAL ANALYSIS STEPS THAT YOU CAN ADD BY UNCOMMENTING BELOW (NOT USED IN THIS PUBLICATION):

# 
# mySwitchList <- analyzeCPAT(
#   mySwitchList,
#   pathToCPATresultFile = NULL,
#   codingCutoff = 0.725,
#   removeNoncodinORFs = TRUE, #Recommended if the isoforms and thus ORFs are predicted de-novo.
#   quiet=FALSE
# )
# 
# mySwitchList <- analyzePFAM(
#   mySwitchList,
#   pathToPFAMresultFile = NULL,
#   showProgress=TRUE,
#   quiet=FALSE
# )
# 
# mySwitchList <- analyzeNetSurfP2(
#   mySwitchList,
#   pathToNetSurfP2resultFile = NULL,
#   smoothingWindowSize = 5,
#   probabilityCutoff = 0.5,
#   minIdrSize = 30,
#   showProgress = TRUE,
#   quiet = FALSE
# )
# 
# mySwitchList <- analyzeSignalP(
#   mySwitchList,
#   pathToSignalPresultFile = NULL,
#   minSignalPeptideProbability = 0.5,
#   quiet=FALSE
# )


# Conducting Analysis Part 2 to conduct more external analyses (CPAT, ORF, PFAM and IDR) and 
# plot all the isoform switches with consequences. 

mySwitchList <- isoformSwitchAnalysisPart2(
  mySwitchList,
  n = Inf,
  codingCutoff = NULL,
  removeNoncodinORFs = TRUE, #Recommended if the isoforms and thus ORFs are predicted de-novo.
  pathToCPATresultFile = NULL,
  pathToCPC2resultFile = NULL,
  pathToPFAMresultFile = NULL,
  pathToIUPred2AresultFile = NULL,
  pathToNetSurfP2resultFile = NULL,
  pathToSignalPresultFile = NULL,
  # ADDITONAL ANALYSIS STEPS THAT YOU CAN ADD BY UNCOMMENTING BELOW IF YOU ADDED THEM ABOVE (NOT USED IN THIS PUBLICATION):
  # consequencesToAnalyze = c(
  #   'intron_retention',
  #   'coding_potential',
  #   'ORF_seq_similarity',
  #   'NMD_status',
  #   'domains_identified',
  #   'IDR_identified',
  #   'IDR_type',
  #   'signal_peptide_identified'
  # ),
  consequencesToAnalyze = c('ORF_seq_similarity'),
  pathToOutput = "PATH/TO/OUTPUT",
  fileType = 'pdf',
  outputPlots = FALSE,
  quiet = FALSE
)


# CUSTOM PLOTS AND TABLES, I RECOMMEND WRITING YOUR OWN 

ggplot(data=mySwitchList$isoformFeatures, aes(x=dIF, y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_hline(yintercept = -log10(0.05), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  facet_wrap( ~ condition_3) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()

ggplot(data=mySwitchList$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) + 
  #facet_wrap(~ condition_2) +
  facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='Gene log2 fold change', y='dIF') +
  theme_bw()


# Add biological consequence analysis as a separate object and 
# analyze the overlap of consequences (ORF and the list below)

bioMechanismeAnalysis <- analyzeSwitchConsequences(
  mySwitchList, 
  consequencesToAnalyze = c('tss','tts','intron_structure', 
                            'last_exon', 'exon_number', 'intron_retention', 
                            '5_utr_seq_similarity', '3_utr_seq_similarity'),
  showProgress = TRUE
)$switchConsequence 

### subset to those with differences
bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoformsDifferent),]

### extract the consequences of interest already stored in the switchAnalyzeRlist
myConsequences <- mySwitchList$switchConsequence
myConsequences <- myConsequences[which(myConsequences$isoformsDifferent),]
myConsequences$isoPair <- paste(myConsequences$isoformUpregulated, myConsequences$isoformDownregulated) # id for specific iso comparison

### Obtain the mechanisms of the isoform switches with consequences
bioMechanismeAnalysis$isoPair <- paste(bioMechanismeAnalysis$isoformUpregulated, bioMechanismeAnalysis$isoformDownregulated)
bioMechanismeAnalysis <- bioMechanismeAnalysis[which(bioMechanismeAnalysis$isoPair %in% myConsequences$isoPair),]  # id for specific iso comparison


### Create list with the isoPair ids for each consequence
AS   <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'intron_structure')]
IR <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'intron_retention')]
aTSS <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'tss'             )]
aTTS <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'tts'             )]
ExNu <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'exon_number'     )]
lastE <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == 'last_exon'      )]
five_utr <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == '5_utr_seq_similarity')]
three_utr <- bioMechanismeAnalysis$isoPair[ which( bioMechanismeAnalysis$featureCompared == '3_utr_seq_similarity')]

mechList <- list(
  AS=AS,
  IR = IR,
  aTSS = aTSS,
  aTTS = aTTS,
  ExNu = ExNu,
  lastE = lastE,
  five_utr = five_utr,
  three_utr = three_utr
)


### Create Venn diagram 
#NOTE: IF YOUR mechList HAS MORE THAN 5 ITEMS, R CANNOT GENERATE A VENN DIAGRAM FOR IT.
        #INSTEAD, PLEASE USE THE UPSET PLOT CODE BELOW.
library(VennDiagram)
#> Loading required package: grid
#> Loading required package: futile.logger
myVenn <- venn.diagram(
  x = mechList,
  col='transparent',
  alpha=0.4,
  fill=RColorBrewer::brewer.pal(n=5,name='Dark2'),
  filename=NULL
)

### Plot the venn diagram
grid.newpage() ; grid.draw(myVenn)


### Create UpSet plot

install.packages("UpSetR")
library(UpSetR)

upset(fromList(mechList), nsets = 7, order.by = "freq")



#EXTRACT CONSEQUENCES, WRITE RESULTS AND GENERATE PLOTS/TABLES:

switches_w_functional_consequences <- extractSwitchSummary(
  mySwitchList,
  filterForConsequences = TRUE
) 
top_10_switches <- extractTopSwitches(mySwitchList, filterForConsequences = TRUE, n=10)

tmp_list <- subset(
  extractTopSwitches(
    mySwitchList,
    filterForConsequences = TRUE,
    n=9999999999999999999999999,
    inEachComparison = TRUE
  )[,c('gene_name','condition_1','condition_2','gene_switch_q_value','Rank')],
  condition_1 == 'acquired_resistant' & condition_2 == 'sensitive'
)

write.csv(tmp_list, "PATH/TO/CSV", row.names=FALSE)

png("image.png", width = 800, height = 600)
switchPlot(
  mySwitchList,
  gene='FAM118A',
  condition1 = 'intrinsically_resistant',
  condition2 = 'sensitive',
  localTheme = theme_light(base_size = 18) # making text sightly larger for vignette
)

extractSwitchOverlap(
  mySwitchList,
  filterForConsequences=TRUE,
  plotIsoforms = FALSE
) + theme(text = element_text(size=20))

extractConsequenceEnrichmentComparison(
  mySwitchList,
  consequencesToAnalyze=c('intron_structure','intron_retention','exon_number'),
  analysisOppositeConsequence = TRUE,
  returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)

extractConsequenceSummary(
  mySwitchList,
  consequencesToAnalyze='all',
  plotGenes = FALSE,           # enables analysis of genes (instead of isoforms)
  asFractionTotal = FALSE      # enables analysis of fraction of significant features
)


tmp_list_2 <- subset(
  extractTopSwitches(
    mySwitchList,
    filterForConsequences = TRUE,
    n=9999999999999999999999999,
    inEachComparison = TRUE
  )[,c('gene_name','condition_1','condition_2','gene_switch_q_value','Rank')],
  condition_1 == 'acquired_resistant' & condition_2 == 'intrinsically_resistant'
)
write.csv(tmp_list_2, "/PATH/TO/CSV", row.names=FALSE)

extractSwitchOverlap(
  mySwitchList,
  filterForConsequences=TRUE,
  plotIsoforms = FALSE
  
)
extractConsequenceEnrichment(
  mySwitchList,
  consequencesToAnalyze='all',
  analysisOppositeConsequence = TRUE,
  returnResult = FALSE # if TRUE returns a data.frame with the summary statistics
)


ggplot(data=mySwitchList$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) + 
  #facet_wrap(~ condition_2) +
  facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','royalblue2')) +
  labs(x='Gene log2 fold change', y='dIF') +
  theme_bw() +
  theme(text = element_text(size=20))


extractSplicingSummary(
  mySwitchList,
  asFractionTotal = TRUE,
  plotGenes=TRUE
)

extractSplicingGenomeWide(
  mySwitchList,
  featureToExtract = 'all',                 # all isoforms stored in the switchAnalyzeRlist
  splicingToAnalyze = c('A3', 'A5', 'ATSS', 'ATTS', 'ES', 'IR', 'MEE', 'MES'), # Splice types of interestcodecode
  plot=TRUE,
  returnResult=FALSE  # Preventing the summary statistics to be returned as a data.frame
)
mySwitchList$isoformFeatures$gene_log2_fold_change

write.csv(mySwitchList[,c("gene_id","condition_1", "condition_2", "gene_name", "gene_log2_fold_change")], file="outfile.csv",row.names=FALSE)
