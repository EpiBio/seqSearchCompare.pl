#!/usr/bin/perl -w
#seqSearchCompare.pl by_Jessilyn_Dunn

#This script takes in 1) a fasta file containing sequences of interest (i.e. promoters: sequences +/-1kb around gene TSS's), and genomic coordinates for features of interest (i.e. a .BED file). It looks for a specific sequence motif in 1) and grabs the Gene accession # and genomic coordinates for each occurrence of the sequence motif. It then finds coordinates from 2) that fall within regions from 1) and records each overlap occurrence as well as the metadata from 2). Finally, it incorporates metadata from a 3rd file that has corresponding gene accession #s (i.e. microarray data for gene expression).

use strict; use Getopt::Long;

############# Variable declaration ############
my $inputFastaFile = ""; #sequences for genomic regions of interest from chosen genome assembly
my $seqsearch = ""; #sequence motif of interest
my $outFile = ""; #sequence motif will be appended to start of outFile name
my $help    = "";
###############################################

# Usage subroutine
sub usage{
    print STDERR "This pipeline takes in 1) The sequences for the genomic regions of interest (-i: fasta format)\n 2).\n";
    print STDERR "$0 -i=<genomic_region_sequence_file> -m=<sequence_motif> -o=<output_filename> [-h] \n";
}

GetOptions ("i=s" => \$inputFastaFile, #sequence file for genomic region of interest (mm9txStart1kbaroundv5.fa is +/-1kb around TSS for mm9 RefSeq)
"m=s" => \$seqsearch,
"o=s" => \$outFile,
"h"   => \$help
)
or die "Invalid commandline argument: $!\n";

if ($help) {usage}

##Part 0. Create the mm9txStart fasta file

#############################################################################################################################
# Part I. Get Gene Accession Number (e.g. NM_x or NR_x) and Genomic Location (e.g. chr, sequence motif start, sequence motif end) for Genes in the chosen assembly (e.g. mm9) that Contain the Sequence Motif in their genomic region of interest.
#############################################################################################################################

my $filename = "$seqsearch".'_InRegions.fa';

#create array to store gene where the match occurs
my $gene = "";                    #gene accession number stored when the CRE exists in the promoter
my @previous_line = ();        #storage of current line
my @results = ();                        #start coord   end coord for each match
my @ranges = ();
my ($match_start, $match_end, $promoter_start, $promoter_end, $CRE_start, $CRE_end) = 0;

open(IN, "<$inputFastaFile") or die "error reading $inputFastaFile: $!\n";
open(OUT, ">$filename") or die "error writing to $filename: $!\n";

while(<IN>){
    if ($_ =~ m/$seqsearch/g)
    {   $gene = substr($previous_line[0], 23);
        @ranges = split /[:=-]+/, $previous_line[1];
        
        $promoter_start = $ranges[2];
        $match_start = pos($_) - length($&);
        $match_end = pos($_);
        $CRE_start = $promoter_start + $match_start-1;
        $CRE_end = $promoter_start + $match_end -1;
        print OUT "$ranges[1]\t$CRE_start\t$CRE_end\t$gene\n";
    }
    @previous_line = split(' ', $_);
}
close IN;
close OUT;

#### Part II. subroutine combine lines; dont need any new input files here.
my $newFileName = "packed_"."$seqsearch".'_GenePromoter_test.fa';
open(IN,"<$filename") or die "cannot open $filename: $!\n";
open(OUT, ">$newFileName") or die "cannot write to $newFileName: $!\n";

my @arr = ();
my @previous_line = ();
my $gene = "";
while (<IN>){
    @arr = split('\t',$_);
    $gene = $arr[3];
    if ($arr[1] !~ $previous_line[1] && $arr[2] !~ $previous_line[2])
    {print OUT "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]";}
    elsif ($arr[0] eq $previous_line[0] && $arr[1] eq $previous_line[1] && $arr[2] eq $previous_line[2])
    {
        seek(OUT, -1, 1);
        print OUT "\t$gene";
    }
    @previous_line = @arr;
}
close IN;
close OUT;

######### then there is the manual addition of the first line. I'll fix this! #####

##Part III:
my $outfile = "$seqsearch"."outFile";

my ($LminusR, $aLminusaR, $L, $R, $aL, $aR, $idx) = 0;
my @genes = ();

######### I will fix this script with the better sorting method so it doesnt take FOREVER to run!
open(RRBS, "<RRBS_intersection_table.txt")  or die "error reading RRBS_intersection_table.txt";
my @rrbs = <RRBS>;
close RRBS;

open(OUT, ">$outfile") or die "error writing to $outfile";
open(IN, "<$newFileName") or die "error reading $newFileName: $!\n";
my @CRE = <IN>;
close IN;

foreach my $region (@CRE){
    my @arr = split('\t', $region);
    $idx = scalar(@arr);
    @genes = @arr[3..($idx-1)];
    print "@arr\t$idx\n";
    
    foreach my $cg (@rrbs){
        my @meth = split('\t', $cg);
        
        if ($meth[0] eq $arr[0]){
            if ($meth[1] >= $arr[1]){
                if ($arr[2] >= $meth[2]){
                    $L = $meth[4];
                    $R = $meth[6];
                    $aL = $meth[8];
                    $aR = $meth[10];
                    $LminusR = 100*($meth[4]-$meth[6]);
                    $aLminusaR = 100*($meth[8]-$meth[10]);
                    print OUT "$meth[0]\t$meth[1]\t$meth[2]\t$LminusR\t$L\t$R\t$aL\t$aR\t$aLminusaR\t$meth[5]\t$meth[7]\t$meth[9]\t$meth[11]\t@genes\n";
                }
            }
        }
        
    }
}

close RRBS;
close OUT;

##### switch to R for visualization. The outfile here is the infile for the R code.

