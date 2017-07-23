#### need to modify to take reverse compliment!


use POSIX;
## vicinity check

my $hg19exonsfile = $ARGV[0];
my $hg19genome = $ARGV[1];

use List::Util qw[min max];
open (F,"$hg19exonsfile");
# boundary:

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
   # if (index($h,"+")<0){ # on minus strand
	   if (1==0){ # on minus strand
#print "$h\treversing$seq\t";
  my $revcomp = reverse($seq);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  $seq=$revcomp;
#print "$seq\n";
}
# if ( index($startkey,32588618)>(-1)){ print "DEBUTGGGGGHERE $startkey and $seq \n";}    
       $hstartseq{$startkey}=$seq;
	$hendseq{$stopkey}=$seq;
## reset everything
    	if ((index ($startkey,"139703300")>(-1) ) | (index ($stopkey,"139703300")>(-1) )){
	print "mykey is $startkey and $stopkey and tinfo $tinfo and $chr and $stop and seqs aren $hendseq{$stopkey} and $hstartseq{$startkey} \n "
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
if ($start==32588618){ print "DEBUTGGGGG $startkey\n";}    
       }

}

close(F);


   
open (F,"num.poly_leq7_nGenesum_geq1.tab");
#sType	fusionInfo	junction	numReads	productPhat	junction_cdf	dirname	sampleType	chr1	gene1	pos1	strand1	chr2	gene2	pos2
    
$h=<F>;

open (FO, ">"."sb_queries_20_all" );


while ($l=<F>){
    chomp($l);
  #  print "$l\n";
    @s=split(/_/,$l);
    @s=split(/\t/,$l); # for sb format
    # $chr1key=2;
    $chr1key=2;
    $pos2key=4; 
    $chr2key=3;
    $pos2key=5; 
    
    $chr1key=8;
    $pos1key=10; 
    $chr2key=12;
    $pos2key=14; 
    $gene1key=9;
    $gene2key=13;
	
    $key1= $s[$chr1key]."_".$s[$pos1key];
    $key2= $s[$chr2key]."_".$s[$pos2key];

    print "keys $key1\t$key2\n";
    if (defined($hstartseq{$key1})){
	$s1=$hstartseq{$key1};
    }
    if (defined($hendseq{$key1})){
	$s1=$hendseq{$key1};
    }
    
    if (defined($hstartseq{$key2})){
	$s2=$hstartseq{$key2};
    }
    if (defined($hendseq{$key2})){
	$s2=$hendseq{$key2};
    }

 #   if ($s[5]==32588618){print "DEBUG $key2\n";}
    #for seven b
    
    $nflank=20;
$s120=substr($s1,length($s1)-$nflank,length($s1));
        $s220=substr($s2,0,$nflank);
    $s110=substr($s1,length($s1)-10,length($s1));
         $s210=substr($s2,0,10);
if (!($s120.$s220 eq "")){     $key1= $s[$chr1key]."_".$s[$pos1key];
    $key2= $s[$chr2key]."_".$s[$pos2key];
    
print FO ">"."$s[$gene1key]-$s[$gene2key]-$s[$pos1key]-$s[$pos2key]_full\n$s120".""."$s220\n";
}}
    
system("bowtie2 -f  --trim5 20 -p 8  --no-sq --no-head   -x ".$hg19genome." -f sb_queries_20_all > test.bowtie.out.3prime");

system("bowtie2 -f  --trim3 20 -p 8  --no-sq --no-head  -x ".$hg19genome." -f sb_queries_20_all > test.bowtie.out.5prime");
