use POSIX;
use Cwd;
## wrapper script to run secondary analysis
## ASSUMES THAT bowtie2 CALLS bowtie2-2.2.1 (and it likely works with later versions, although
## that has not been verified)

## FILES (NOT CODE) THAT MUST BE IN WORKING DIR (i.e. in cwd()):
## cancer_gene_census.csv, knife_badFJ, refseq_genes_for_frames_selected, ucscGenePfam,
## COSMIC_fusions.csv, onc2014406x2.csv, ChimerDB3.0_ChimerSeq.csv, 
## ChimerDB3.0_ChimerSeq.txt, ChimerDB3.0_ChimerPub.txt, NEJM_LAML_comparison.csv,
## hg19_crude_geneann

## Need R libraries data.table, calibrate, binom

# code directory; expect all code that this calls to be there; this should also be the cwd()
my $codedir = $ARGV[0];
# report csv file with info about paths of data
my $reportpathsfile = $ARGV[1];

# file for metadata
my $metadatafile = $ARGV[2];

my $hg19exonsfile = $ARGV[3];
my $hg19genome = $ARGV[4];
my $countsfile = $ARGV[5];

# working directory; expect files to be there; this should be the same as $codedir
my $wdir = cwd();

# directory for output files, within working directory:
my $outputfilesdir = "${wdir}/output_files/";

unless(-e $outputfilesdir){
  mkdir $outputfilesdir;
}

open (F,"$reportpathsfile");

###########

my $myheader = <F>;  # just for throwing it away so we don't read the header

while ($l=<F>) {
  chomp($l);
  @s=split(/,/,$l);
  $nicename=substr($s[4],1,length($s[4])-2);
  $fileout="${outputfilesdir}/".$nicename."__all_merged_w_anomalies";
  if ((-e $fileout)) {		## if exists, we don't re-run
    print "$fileout exists\n";
  }
  if (!(-e $fileout)) {
    print "did not find $fileout\n";
    print "Rscript tcga_process.r $s[0] $s[1] $s[2] $wdir \n";
    system("Rscript tcga_process.r ".$s[0]." ".$s[1]." ".$s[2]." ".$wdir);
    print "Rscript append_anomalies.r $s[0] $s[1] $s[2] $s[3] $s[4] $outputfilesdir \n";
    system("Rscript append_anomalies.r ".$s[0]." ".$s[1]." ".$s[2]." ".$s[3]." ".$s[4]." ".$outputfilesdir);
  }
}
################


## vicinity check
open(F,"ucscGenePfam");# |grep kinase|grep 18086
$rf=150000;

use List::Util qw[min max];
my %protein;
while ($l=<F>){
    chomp($l);
    @s=split(/\t/,$l);
    $chr=$s[1]; #s2 and s3 are start and stop
    $start=min($s[3],$s[2]);
    $stop=max($s[2],$s[3]);
    
    $rp1=$chr."_".$rf*ceil($s[2]/$rf);    
    $rp2=$chr."_".$rf*ceil($s[3]/$rf);
    $addin="$s[1]-$s[2]-$s[3]-$s[4]";
    if (!defined($protein{$rp1})){$protein{$rp1}="";}
    if (!defined($protein{$rp2})){$protein{$rp2}="";}
        $protein{$rp1}=$protein{$rp1}.$start."_".$stop."_".$addin;
        $protein{$rp2}=$protein{$rp2}.$start."_".$stop."_".$addin;
}
close(F);

###################################################
### find frame of exon which is not used in paper #
###################################################
open (F,"refseq_genes_for_frames_selected");
my %myframe;
$l=<F>;
$debug=0;


while ($l=<F>){
    @s=split(/\t/,$l);
    chomp($l);
    @frames=split(/,/,$s[7]);
    $count=0;
    my @cumframe;
$cumframe[0]=0;

    foreach $fr (@frames){
	if ($fr>(-1)){ ## -1 is the convention for a non-existent frame
	$cumframe[$count]=$cumframe[($count-1)]+$fr;
#	print "$l\n$cumframe is now $cumframe[$count]\t$count\n";
	}
	    $count=$count+1;
}
    
    @starts=split(/,/,$s[4]);
        @stops=split(/,/,$s[5]);
    $nadd=length(@frames);
    $chr=$s[2];
    #    print "debug $nadd fields starts @starts and @frames\n";
    $jc=0;
    foreach $start (@starts){
	$ustart=$start+1;
	if ($start==1741428){print "TACC\t$l";}
	$key1=$s[2]."_".$ustart;
	$key2=$s[2]."_".$stops[$jc];
	$val=$frames[$jc];
	$val2=$cumframe[$jc];
	if ($jc>0){	$pval=$cumframe[($jc-1)];}
		if ($jc==0){	$pval=0;}
	
	$myframe{$key}=$val;
#	print "adding $key1 and $key2 and $val	and adding $val2\n";
	$myframe{$key1}=$val."_".$pval."_".$val2;
		
	$myframe{$key2}=$val."_".$pval."_".$val2;

	$jc=$jc+1;
    }
}
open (F,"hg19_exons.fa");

# OPEN and append the start and stop of exonic sequences to reports

$h=<F>;
chomp($h);
	@s=split(/ /,$h);
	$tinfo=substr($s[1],6,40);
	($chr,$range)=split(/:/,$tinfo);
	($start,$stop)=split(/-/,$range);
$stopkey=$chr."_".$stop;
$startkey=$chr."_".$start;
my %hstartseq;
my %hendseq;
$seq="";
## strand is - results in reversing the exon sequence so everything is wrt + convention
while ($h=<F>){
    chomp($h);
    if (index ($h,"known")<0){
	$seq=$seq.$h;}
if (index ($h,"known")>0){
    ## add in sequence there
    ## add sequence into stat and end

	$startseq=substr($seq,0,20);
	$endseq=substr($seq,length($seq)-20,length($seq));
 
       $hstartseq{$startkey}=$startseq."_".$endseq;
	$hendseq{$stopkey}=$startseq."_".$endseq;
	## after this is assigned, reset everything
	$seq="";
	##  show example of a single exon identified
	
	if ((index ($startkey,"139703300")>(-1) )){# | (index ($stopkey,"139703300")>(-1) )){#139703300
	    print "mykey is $startkey and $stopkey and tinfo $tinfo and $chr and $stop and seqs aren $hendseq{$stopkey} and $hstartseq{$startkey} \n ";
	}
	@s=split(/ /,$h);
	$tinfo=substr($s[1],6,40);
	($chr,$range)=split(/:/,$tinfo);
	($start,$stop)=split(/-/,$range);
		#print "added $startkey, $stopkey,$startseq $endseq\n";
	$seq="";
	## define new keys
	$stopkey=$chr."_".$stop;
	$startkey=$chr."_".$start;
           }

}

close(F);
opendir (DIR, $outputfilesdir) or die $!;

while (my $file = readdir(DIR)) {
#print "file is $file\n"; 
if ((index($file,"appended")<0) && (index($file,"merged") > 0) ){
    open (F, $outputfilesdir.$file);
    
    $h=<F>;
@header=split(/\t/,$h);

$i=0;
    ## automatically assign headers
    
    foreach $s (@header){
	if ($s eq "chr1"){$chr1i=$i;}
	if ($s eq "chr2"){$chr2i=$i;}
	if ($s eq "pos1"){$pos1i=$i;}
	if ($s eq "pos2"){$pos2i=$i;}
		$i=$i+1;
}

    ### end automatic
    
    ## apended files have protein frame and exonic sequence information, only for meta-processing as described in methods
    
    open (FO, ">".$file."_appended" );
    ## describe added fields:
    chomp($h);
print FO $h."\tseq1\tseq2\tframe1\tframe2\n";    
while ($l=<F>){
    chomp($l);
    @s=split(/\t/,$l);


        $key2= $s[$chr2i]."_".$s[$pos2i];
## begin protein syc
## pairs of fusions are 2,4 and 3,5

    $chr1=$s[$chr1i];
    $rp1=$chr1."_".$rf*ceil($s[$pos1i]/$rf);
    $chr2=$s[$chr2i];
    $rp2=$chr2."_".$rf*ceil($s[$pos2i]/$rf);
    $prot1="NONE\tNONE\tNONE\tNONE";
    $prot2="NONE\tNONE\tNONE\tNONE";
    
    if (defined($protein{$rp1})){
	($coord1,$coord2,$rest)=split(/_/,$protein{$rp1});
	if (($s[4]<$coord2) & ($s[4]>$coord1)){  
	    $prot1=$protein{$rp1};
}
}
    if (defined($protein{$rp2})){
	($coord1,$coord2,$rest)=split(/_/,$protein{$rp2});
	if (($s[3]<$coord2) & ($s[5]>$coord1)){  
	$prot2=$protein{$rp2};
	}
}
    
    $key1= $s[$chr1i]."_".$s[$pos1i];
    #print "key 1 $key1 key2 $key2\n";
    $s1="";
        $s2="";
    if (defined($hstartseq{$key1})){
	$s1=$hstartseq{$key1};
	$frame1=$myframe{$key1};
#	print "found frame $frame1 with $key1 \n";
    }
    if (defined($hendseq{$key1})){
	$s1=$hendseq{$key1};
	$frame1=$myframe{$key1};
#	print "found frame $frame1 with $key1 \n";
   }
    
    if (defined($hstartseq{$key2})){
	$s2=$hstartseq{$key2};
	$frame2=$myframe{$key2};
#	print "found frame $frame2 with $key2 \n";
    }
    if (defined($hendseq{$key2})){
	$s2=$hendseq{$key2};
		$frame2=$myframe{$key2};
#	print "found frame $frame2 with $key1 \n";
   }
    
#       print  "$l\t$s1\t$s2\t$frame1\t$frame2\t$prot1\t$prot2";
    print FO "$l\t$s1\t$s2\t$frame1\t$frame2\n"; #t$prot1\t$prot2\n";

}
    }}
## after protein domains and sequences are added, secondary processing is run
system("Rscript secondary_process.r ".$metadatafile." ".$wdir);
system("perl fusion_mapping_qual.pl ".$hg19exonsfile." ".$hg19genome);
system("Rscript meta_process.r");

