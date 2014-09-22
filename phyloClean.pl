#!/usr/bin/perl

#Author: JP Gallagher
#Last updated: 12/Aug/2014
#Title: phyloClean.pl
#Keywords: phylogenetics, gap, ambiguous sites
#Some days you get up thinking you are almost done with something and then you
#realize how much more s*** you have to do.  This is one of those days.  This
#script takes in an alignment, eliminates all positions with a gap, then
#removes all sequences with greater than xx% ambiguous sites.  The alignment
#can be in any standard alignment format.  

#Note: "Using a hash as a reference is deprecated at phyloClean.pl line 75." is
#not an error that affects the functioning of this script

#standard pragmas for good perl coding
use warnings;
use strict;

#requires specific BioPerl modules
use Bio::AlignIO;#alignment input/output
use Bio::SimpleAlign;#alignment handling
use Bio::Seq;#sequence handling
use Bio::Tools::SeqStats;#basic sequence statistics methods
use Bio::LocatableSeq;#locatable sequence object

#Input of parameters and warning message
#Usage: phyloClean.pl <raw_msa> <percent_Ns_upper_limit> <output_name>
unless (@ARGV == 3){
	print 	"phyloClean.pl - a basic clean up tool for phylogenetic analysis\n",
			"Usage: phyloClean.pl <raw_msa> <percent_Ns_upper_limit> <output_name>\n",
			"This script uses BioPerl. For info see http://www.bioperl.org/ \n";
	exit 1;
}
#set input to specific scalars
my ($msa, $percent, $newMsa) = @ARGV;

#read in alignment using AlignIO
#use default endings to specify format.  These can be found at
#http://www.bioperl.org/wiki/HOWTO:AlignIO_and_SimpleAlign
my $in = Bio::AlignIO->new(-file => $msa);

#set up output file, default endings can be found at the website above

my $aln;
my $degapAln;
my $workingAln;
my $newAln;
#read in alignments
while( $aln = $in->next_aln ){
	#degap the alignment
	$degapAln = $aln->remove_gaps;
	#remove seqs with greater than xx% ambiguous characters
	if($percent == 100){
		#if you set the upper limit to 100% Ns, skips checking seqs
		#also, why don't you just use trimal or Gblocks in that case?
		$newAln = $degapAln;
	}
	else{
		#create working SimpleAlign object
		$workingAln = $degapAln;
		#determine missing data character for querying the hash later
		my $missing = $workingAln->missing_char();
		my @seqlist = ();
		for(my $i = 1; $i < $workingAln->num_sequences(); $i++){
			my $seq = $workingAln->get_seq_by_pos($i);
			my $seqlength = $seq->length();
			#calculate sequence base statistics
			my $seqstats = Bio::Tools::SeqStats->new(-seq => $seq);
			my $basestats = $seqstats->count_monomers();
			my $ACGT_total = 0; 
			#calculate the total number of recognized bases
			#gaps should be gone already so those should not affect the calculation
			foreach my $base (sort keys %$basestats) {
         			if($base =~ /[ACGT]/){
					#limits this to only definitive bases
					$ACGT_total += %$basestats->{$base};
				}
			}
			my $missingNo = $seqlength - $ACGT_total;
			#removes sequences from alignment that have more Ns then allowed
			if(($missingNo / $seqlength) <= $percent/100){
				push(@seqlist,$i);	
			}
		}
		#sets alignment after removing sequences
		if(@seqlist != 0){
			$newAln = $workingAln->select_noncont(@seqlist);
		}
	}
	if(defined($newAln) && $newAln->num_sequences() > $aln->num_sequences()/2){
		#writes out the new alignment
		my $out = Bio::AlignIO->new(-file => ">$newMsa");
		$out->write_aln($newAln);
	}
	else{
		#doesn't print alignment if bad
		print STDERR "Alignment $msa did not meet criteria.\n";
	}
}

exit 0;



