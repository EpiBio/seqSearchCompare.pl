# seqSearchCompare.pl

This script takes in 1) a FASTA file containing sequences of interest (i.e. promoters: sequences +/-1kb around gene TSS's), and genomic coordinates for features of interest (i.e. a .BED file). It looks for a specific sequence motif in 1) and grabs the Gene accession # and genomic coordinates for each occurrence of the sequence motif. It then finds coordinates from 2) that fall within regions from 1) and records each overlap occurrence as well as the metadata from 2). Finally, it incorporates metadata from a 3rd file that has corresponding gene accession #s (i.e. microarray data for gene expression).

######################## seqSearchCompare.pl by Jessilyn Dunn 2015-01-06 ########################

I.  Download promoter sequences from UCSC (+/- 1kb around TSS)


1. Download UCSC genes for mm9 (UCSC genome browser -> table browser -> download RefSeq genome)

2. Make a new file for the "address" of the promoters (containing the chr and coordinates of the upstream and downstream portion of the TSS) - update this file to change the txStart for the - strand to be the txEnd

3. Upload the "address" file as a bed file (here, TSS_address.BED) into UCSC browser, and use the table browser to download the sequence of the defined regions.

4. The downloaded sequence is in FASTA format (mm9txStart1kbaround.fa). Search for the TFBS of interest: grep -o 'motif' mm9txStart1kbaroundv2.fa | wc -l

4.5 reformat mm9txStart1kbaround.fa so that each entry occupies one line:

cat mm9txStart1kbaround.fa | sed 's/\(^>.*\)/!\1!/' | tr '\n' '@' | tr '!' '\n' | sed 's/@//g' > mm9txStart1kbaroundv2.fa

5. Use this file as your "regions" file as an input into CRE_search.pl

Note: there is still a small bug in Part II. of CRE_draft.pl where the first line of the file doesnt insert, will fix this when I have a moment


II. In command line and in R:

6. threshold at 10 reads/site for sequencing accuracy: awk '$8 >= 10 && $9 >= 10 && $10 >= 10 && $11 >=10 {print $0}' sorted_TGACGTCA_final.fa > 10reads_sorted_TGACGTCA_final.fa

The final files contain: 
Gene Accession # | methL | methR | methaL | methaR | Accession# | expL | expR | expaL | expaR |

Import into R, and create scatterplots to determine if there is correlation between CRE methylation and gene expression (probe intensity): 
$2 vs $7
$3 vs $8
$4 vs $9
$5 vs $10
