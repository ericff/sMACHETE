use POSIX;
use Cwd;
## wrapper script to run secondary analysis
## ASSUMES THAT bowtie2 CALLS bowtie2-2.2.1 (and it likely works with later versions, although
## that has not been verified)

## FILES (NOT CODE) THAT MUST BE IN WORKING DIR (i.e. in cwd()):
## cancer_gene_census.csv, knife_badFJ, refseq_genes_for_frames_selected, ucscGenePfam,
## COSMIC_fusions.csv, onc2014406x2.csv, ChimerDB3.0_ChimerSeq.csv, 
## ChimerDB3.0_ChimerSeq.txt, ChimerDB3.0_ChimerPub.txt, NEJM_LAML_comparison.csv,
## rca_tcga_pub2015.data_mutations_extended.txt,  ov_tcga_pub.data_mutations_extended.txt
## gbm_tcga_pub2013.data_mutations_extended.txt,   paad.tcga.data_mutations_extended.txt
## laml_tcga_pub.data_mutations_extended.txt,      prad_tcga_pub.data_mutations_extended.txt
## luad_tcga_pub.data_mutations_extended.txt

## A DIRECTORY tempfiles MUST BE IN WORKING DIR (i.e. in cwd()):

## Need R libraries data.table, calibrate, binom

my $countsfile = $ARGV[0];
my $querydir = $ARGV[1];
# querydir should contain files with names containing "summary.bt.results" and contain
# rows of the form
# GARNL3-APBB3-999999999-999999999_full	0	0.0	0	0.0	0.0	0.0
# that result from processing of SBT queries

my $sequencespresentlamlfile = $ARGV[2];
# sequencespresentlamlfile is a file that gives one line for each (fusion)-(tcga case) pair
# and a 1 or 0 if the fusion is detected by the SBT query as present in that TCGA sample
# format:
# sequence.name,sample.id,is.sequence.present
# AAAA-BBBB-9999999-9999999_full,TCGA-AB-9999-99A,0
# CCCC-DDDD-9999999-9999999_full,TCGA-AB-9999-99A,0

# report csv file with info about paths of data, MERGED with metadata
my $reportpathswithmetafile = $ARGV[3];

# results of particular SBT query for bodymap 
my $sbtresultsbodymapchimera = $ARGV[4];

# directory for certain files giving one line for each sequence and sample combination
#  and 0/1 if sequence is found to be present by SBT queries
my $nametosampleidsdir = $ARGV[5];

# date string associated with queries on SB; e.g. "jan6" or "jun12"
my $querydate = $ARGV[6];



# working directory; expect files to be there; this should be the same as $codedir
my $wdir = cwd();

chdir("SB");
#
system("Rscript R_scripts_run.r ".$countsfile." ".$querydir." ".$sequencespresentlamlfile." ".$reportpathswithmetafile." ".$sbtresultsbodymapchimera." ".$nametosampleidsdir." ".$querydate);
# will automatically run system("Rscript analysis.r");
#("FigX_Mach_fusion_per_cancer.r");
#("summary_fig_plots.r");

