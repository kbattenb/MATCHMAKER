#!/usr/bin/perl

##############################################################################
#                                                                            #
#                              "MATCHMAKER.pl"                               #
#                                                                            #
#     Given a motif sequence(s) and a reference dataset, runs                #
#     an HMM search iteratively until the search converges.                  #
#                                                                            #
#     Written by Kai Battenberg                                              #
#     RIKEN CSRS Plant Symbiosis Research Team                               #
#                                                                            #
#    Copyright (C) 2022 Kai Battenberg                                       #
#                                                                            #
#    This program is free software: you can redistribute it and/or modify    #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    This program is distributed in the hope that it will be useful,         #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.    #
#                                                                            #
##############################################################################



use strict;
use warnings;
use Getopt::Long;



######Options#####
#script dir
my $script_dir = "<SET/PATH/TO/WHERE/YOU/HAVE/THIS/SAVED>";
my $overview = $script_dir."/Pipeline_overview.jpg";
my $version = "v2022.05.17";

#default values.
my $element_name = ""; #Set an ID for the element being searched
my $seq_file = ""; #Strand that will be subject to search [positive, negative, or both]
my $database = ""; #Path to the FASTA file used as the database
my $out_dir = ""; #Path to the outpud directory

my $maxround = 30; #Maxmum number of cycles run without convergence before terminating the process
my $maxwindow = 100; #Maximum window size allowed for the search in bps
my $strand = "positive"; #Strand that will be subject to search [positive, negative, or both]
my $parsemode = "ALL"; #Hits to keep for each database sequence for each round [ALL or BEST]
my $flank_size = 5; #Size of the flanking sequence (in bps) taken around the aligned segment to be kept
my $alignment = "mafft"; #Method to align parsed hits [mafft or mhhalign]
my $minbit = 1; #Minimum bit value for a position within the profile to be kept during alignment trimming
my $minweight = 0.8; #Minimum weight for a position within the profile to be kept during alignment trimming

my $threads = 2; #Number of threads passed to any of the tools called within MATCHMAKER.pl
my $help;

#help
my $help_display = "
     MATCHMAKER.pl ($version)
     Usage: MATCHMAKER.pl --element_name <ELEMENT NAME> --query <PATH/TO/FASTA> --database <PATH/TO/FASTA> --out_dir <PATH/TO/OUTPUT/DIR>
     
     --help             Display this help and exit
     
     --element_name     Set an ID for the element being searched
     --query            Path to the FASTA file used as the query
     --database         Path to the FASTA file used as the database
     --out_dir          Path to the outpud directory
     
     --maxround         Maxmum number of cycles run without convergence before terminating the process (default: 30)
     --maxwindow        Maximum window size allowed for the search in bps (default: 100)
     --strand           Strand that will be subject to search [positive, negative, or both] (default: positive)
     --parsemode        Hits to keep for each database sequence for each round [ALL or BEST] (default: ALL)
     --flank_size       Size of the flanking sequence (in bps) taken around the aligned segment to be kept (default: 5)
     --alignment        Method to align parsed hits [mafft or mhhalign] (default: mafft)
     --minbit           Minimum bit value for a position within the profile to be kept during alignment trimming (default: 1)
     --minweight        Minimum weight for a position within the profile to be kept during alignment trimming (default: 0.8)
     
     --threads          Number of threads passed to any of the tools called within MATCHMAKER.pl (default: 2)
     
";
##########



#####Checking for tools#####
#check for conda
my $conda_test = `which conda`;
if ($conda_test eq "") {
	$conda_test = "FAIL";
}
else {
	$conda_test = "PASS";
}
my $conda_version = "(not found)";
if ($conda_test eq "PASS") {
	$conda_version = `conda --version | cut -f2 -d' ' | perl -pe 'chomp'`;
}

#ckeck for perl module POSIX
my $posix_test = "FAIL";
eval "use POSIX";
if ($@ ne "") {
	print "ERROR: Perl module POSIX is required but not installed.\n";
	if ($conda_test eq "FAIL") {
		print "***conda was also not found in path. All tools required by this pipeline is available through conda. It is not required but highly recommended.***\n";
	}
	elsif ($conda_test eq "PASS") {
		print "Install POSIX with the following command:\n";
		print "\n\tconda install -c bioconda perl-posix\n\n";
	}
	die "install POSIX and please try again.\n";
}
else {
	$posix_test = "PASS";
}
require POSIX;

#check for perl module Math CDF
my $cdf_test = "FAIL";
eval "use Math::CDF";
if ($@ ne "") {
	print "ERROR: Perl module Math::CDF is required but not installed.\n";
	if ($conda_test eq "FAIL") {
		print "***conda was also not found in path. All tools required by this pipeline is available through conda. It is not required but highly recommended.***\n";
	}
	elsif ($conda_test eq "PASS") {
		print "Install Math::CDF with the following command:\n";
		print "\n\tconda install -c bioconda perl-math-cdf\n\n";
	}
	die "install Math::CDF and please try again.\n";
}
else {
	$posix_test = "PASS";
}
require Math::CDF;

#check for MAFFT
my $mafft_test = `which mafft`;
if ($mafft_test eq "") {
	print "ERROR: MAFFT is required but not found in path.\n";
	if ($conda_test eq "FAIL") {
		print "***conda was also not found in path. All tools required by this pipeline is available through conda. It is not required but highly recommended.***\n";
	}
	elsif ($conda_test eq "PASS") {
		print "Install MAFFT with the following command:\n";
		print "\n\tconda install -c bioconda mafft\n\n";
	}
	die "install MAFFT and please try again.\n";
}
else {
	$mafft_test = "PASS";
}
my $mafft_version = `mafft --help 2>&1 | grep 'MAFFT' | head -n 1 | cut -f4 -d' ' | perl -pe 'chomp'`;

#check HMMER
my $hmmer_test = `which nhmmscan`;
if ($hmmer_test eq "") {
	print "ERROR: HMMER is required but not found in path.\n";
	if ($conda_test eq "FAIL") {
		print "***conda was also not found in path. All tools required by this pipeline is available through conda. It is not required but highly recommended.***\n";
	}
	elsif ($conda_test eq "PASS") {
		print "Install HMMER with the following command:\n";
		print "\n\tconda install -c bioconda hmmer\n\n";
	}
	die "install HMMER and please try again.\n";
}
else {
	$hmmer_test = "PASS";
}
my $hmmer_version = `hmmbuild -h | grep HMMER | cut -f3 -d' ' | perl -pe 'chomp'`;

#check WebLogo
my $weblogo_test = `which weblogo`;
if ($weblogo_test eq "") {
	print "ERROR: WebLogo is required but not found in path.\n";
	if ($conda_test eq "FAIL") {
		print "***conda was also not found in path. All tools required by this pipeline is available through conda. It is not required but highly recommended.***\n";
	}
	elsif ($conda_test eq "PASS") {
		print "Install weblogo with the following command:\n";
		print "\n\tconda install -c bioconda weblogo\n\n";
	}
	die "install WebLogo and please try again.\n";
}
else {
	$weblogo_test = "PASS";
}
my $weblogo_version = `weblogo --version | cut -f2 -d' ' | perl -pe 'chomp'`;
##########



#####Setting up options#####
#making the options into external arguments.
GetOptions (
	'element_name=s' => \$element_name,
	'query=s' => \$seq_file,
	'database=s' => \$database,
	'out_dir=s' => \$out_dir,
	'maxround=s' => \$maxround,
	'maxwindow=s' => \$maxwindow,
	'strand=s' => \$strand,
	'parsemode=s' => \$parsemode,
	'flank_size=s' => \$flank_size,
	'alignment=s' => \$alignment,
	'minbit=s' => \$minbit,
	'minweight=s' => \$minweight,
	'threads=s' => \$threads,
	'help!' => \$help
	);
##########



#####Return help#####
if (defined $help) {
    die "$help_display";
}
##########



#####Checking Options#####
print "\n";
print "#####BEGIN MATCHMAKER $version#####\n";
print "checking options.\n";

#checking for required options.
if (!$element_name) {
	print "USAGE: option --element_name <ELEMENT NAME> is required.\n";
	die "##########\n\n";
}
elsif (!$seq_file) {
	print "USAGE: option --query <QUERY FASTA> is required.\n";
	die "##########\n\n";
}
elsif (!$database) {
	print "USAGE: option --database <DATABASE FASTA> is required.\n";
	die "##########\n\n";
}
elsif (!$out_dir) {
	print "USAGE: option --out_dir is required.\n";
	die "##########\n\n";
}

#checking option database
if (-e $database) {
	print " Database checked.\n";
}
else {
	print "Error: Selected database $database does not exist.\n";
	die "##########\n\n";
}

#checking option out_dir
if (-e $out_dir && -d $out_dir) {
	print " Output directory checked.\n";
}
else {
	print "Error: Selected output directory $out_dir does not exist.\n";
	die "##########\n\n";
}

#checking window
if ( $maxwindow !~ /^-?\d+\.?\d*$/ ) {
	print "Error: option --maxwindow needs to be numeric\n";
	die "##########\n\n";
}
elsif ( $maxwindow < 1 || $maxwindow !~ /^-?\d+$/ ) {
	print "Error: option --maxwindow needs to be a positive integer\n";
	die "##########\n\n";
}

#checking strand
if ( $strand !~ m/positive/ && $strand !~ m/negative/ && $strand !~ m/both/ ) {
	print "Error: option --strand needs to be positive, negative, or both\n";
	die "##########\n\n";
}

#checking parsemode
if ( $parsemode !~ m/^BEST$/ && $parsemode !~ m/^ALL$/ ) {
	print "Error: option --parsemode needs to be BEST, or ALL\n";
	die "##########\n\n";
}

#checking nat
if ( $minbit !~ /^-?\d+\.?\d*$/ ) {
	print "Error: option --minbit needs to be numeric\n";
	die "##########\n\n";
}
elsif ( $minbit <= 0 || $minbit >= 2 ) {
	print "Error: option --minbit needs to be within 0 and 2\n";
	die "##########\n\n";
}

#checking minweight
if ( $minweight !~ /^-?\d+\.?\d*$/ ) {
	print "Error: option --minweight needs to be numeric\n";
	die "##########\n\n";
}
elsif ( $minweight <= 0 || $minweight >= 1 ) {
	print "Error: option --minweight needs to be within 0 and 1\n";
	die "##########\n\n";
}

#checking flank_size
if ( $flank_size !~ /^-?\d+\.?\d*$/ ) {
	print "Error: option --flank_size needs to be numeric\n";
	die "##########\n\n";
}
elsif ( $flank_size < 0 || $flank_size !~ /^-?\d+$/ ) {
	print "Error: option --flank_size needs to be a positive integer\n";
	die "##########\n\n";
}

#checking alignment
if ( $alignment !~ m/mafft/ && $alignment !~ m/hmmalign/ ) {
	print "Error: option --alignment needs to be mafft or hmmalign\n";
	die "##########\n\n";
}

print "check complete\n\n";
##########



#####Setting output file names#####
my $out = $out_dir."/".$element_name."_MATCHMAKER";
my $input_dir = $out."/INPUT";
my $log = $out."/".$element_name."_MATCHMAKER_log.txt";
##########



#####Setup input#####
#set up all directories
system "rm -rf $out";
system "mkdir $out";
system "mkdir $input_dir";

#format query and database
my $query = $input_dir."/QUERY.fas";
&singler ($seq_file, $query);
my $subject = $input_dir."/DATABASE.fas";
&singler ($database, $subject);

#convert bit to nat
my $nat = &bit2nat($minbit);

#entering parameters to log
open (LOG, ">$log") or die "cannot open $log.\n";
print LOG "#####BEGIN MATCHMAKER $version#####\n";
print "INPUT:\n";
print LOG "INPUT:\n";
print " Element:\t$element_name\n";
print LOG " Element:\t$element_name\n";
print " Query:\t$seq_file\n";
print LOG " Query:\t$seq_file\n";
print " Database:\t$database\n";
print LOG " Database:\t$database\n";
print " Window size:\t$maxwindow\n";
print LOG " Window size:\t$maxwindow\n";
print " Keep strand: $strand\n";
print LOG " Keep strand: $strand\n";
print " Keep hits per gene: $parsemode\n";
print LOG " Keep hits per gene: $parsemode\n";
print " Minimum bits per position: $minbit\n";
print LOG " Minimum bits per position: $minbit\n";
print " Minimum weight per position: $minweight\n";
print LOG " Minimum weight per position: $minweight\n";
print " Tools:\tHMMER $hmmer_version, MAFFT $mafft_version, WebLogo $weblogo_version\n";
print LOG " Tools:\tHMMER $hmmer_version\tMAFFT $mafft_version\tWebLogo $weblogo_version\n";
print " Alignment method: $alignment\n";
print LOG " Alignment method: $alignment\n";
##########



#####MAIN Search for genes with motif#####
my $test = "GO";
my $round = 0;
my $final_seq_count = 0;
my $final_evalue = 0;

#loop until search converges
while ($test =~ m/GO/) { #Loop to keep the search running
	###set up ###
	#set up round
	$round++;
	my $round_name = "ROUND".sprintf("%.3d", $round);
	my $round_dir = $out."/".$round_name;
	system "mkdir $round_dir";
	print "Round $round:\n";
	print LOG "Round $round:\n";

	#setup files
	my $q_file = $round_dir."/01_query.fas";
	my $prof_file = $round_dir."/02_profile.hmm";
	my $prof_logo_file = $round_dir."/02_profile_logo.pdf";
	my $hit_file = $round_dir."/03_hit.txt";
	my $besthit_file = $round_dir."/04_hit_parsed.txt";
	my $besthit_fas = $round_dir."/04_hit_parsed.fas";
	my $extract_file = $round_dir."/05_hit_segments.fas";
	my $aln_file = $round_dir."/06_segments_aligned.fas";
	my $aln_logo_file = $round_dir."/06_segments_aligned_logo.txt";
	my $aln_cut_file = $round_dir."/07_segments_aligned_cut.fas";
	######



	###STEP1: Copy query file###
	if ($round < 2) {
		system "cp $query $q_file";
	}
	else {
		my $ref_round = $round - 1;
		my $ref_round_name = "ROUND".sprintf("%.3d", $ref_round);
		my $ref_file = $out."/".$ref_round_name."/07_segments_aligned_cut.fas";
		system "cp $ref_file $q_file";
	}
	#set seq count
	my $seq_count = `wc -l $q_file | sed 's/^ *//g' | cut -f1 -d' ' | perl -pe 'chomp'`;
	$seq_count = $seq_count / 2;
	$final_seq_count = $seq_count;
	my @seqlist_in = split (/ /, `cat $q_file | grep '>' | perl -pe 'chomp' | sed 's/>/ /g' | sed 's/^ //' | perl -pe 'chomp'`);
	print " Query seq count:\t$seq_count\n";
	print LOG " Query seq count:\t$seq_count\n";
	######



	###STEP2: Build HMM profile###
	my $hmm_name = $element_name."_".$round_name;
	my $hmmbuild_cmd = &hmmbuilder ($q_file, $hmm_name, $prof_file, $threads);
	my $weblogo_cmd = &weblogger ($q_file, $prof_logo_file, "pdf");
	system "hmmpress $prof_file > /dev/null";
	print " HMM profile name:\t$hmm_name\n";
	print LOG " HMM profile name:\t$hmm_name\n";
	print LOG "  CMD:\t$hmmbuild_cmd\n";
	print LOG "  CMD:\thmmpress $prof_file\n";
	print LOG "  CMD:\t$weblogo_cmd\n";
	my $seg_length = `grep 'LENG' $prof_file | rev | cut -f1 -d' ' | rev | perl -pe 'chomp'`;
	if ($seg_length >= $maxwindow) {
		die "Error: Sequence length at Rount $round exceeds the window length. Increase windowlength with --maxwindow.\n";
	}
	my $slide = $maxwindow - $seg_length + 1;
	if ($maxwindow/2 < $slide) {
		$slide = floor($maxwindow/2);
	}
	print " Slide width:\t$slide\n";
	print LOG " Slide width:\t$slide\n";
	######



	###STEP3: Search the database###
	my $subject_temp = $out."/TEMPSUBJECT_".$maxwindow."_".$slide.".fas";
	if (! -e $subject_temp) {
		&databasesegmenter ($subject, $maxwindow, $slide, $subject_temp);
	}
	my $nhmmscan_cmd = &nhmmscaner ($prof_file, $subject_temp, $hit_file, $threads);
	print "  nhmmscan completed.\n";
	print LOG "  CMD: $nhmmscan_cmd\n";
	######



	###STEP4: Parse out the hits###
	my @evals = (0.001, 0.005, 0.01, 0.05, 0.1, 0.5);
	my $evalue;
	foreach my $eval (@evals) {
		$evalue = $eval;
		&hitparser ($hit_file, $strand, $eval, $parsemode, $besthit_file);
		my $passcount = `wc -l $besthit_file | rev | cut -f2 -d' ' | rev | perl -pe 'chomp'`;
		if ($passcount > 2) {
			last;
		}
	}
	$final_evalue = $evalue;
	my $besthit_test = `wc -l $besthit_file | rev | cut -f2 -d' ' | rev | perl -pe 'chomp'`;
	if ($besthit_test < 2) {
		system "rm $out/TEMPSUBJECT_*";
		system "cp $overview $out";
		print "NO HIT FOUND (process complete)\n";
		print LOG "NO HIT FOUND (process complete)\n";
		print LOG "##########\n";
		close (LOG);
		die "##########\n\n";
	}

	my %alignments = %{ &hits2hash($besthit_file, $subject) };
	open (BESTFAS, ">", $besthit_fas) or die "cannot open $besthit_fas.\n";
	foreach my $gene (sort keys %alignments) {
		print BESTFAS ">$gene\n";
		print BESTFAS "$alignments{$gene}\n";
	}
	close (BESTFAS);

	print " E-value threshold:\t$evalue\n";
	print LOG " E-value threshold:\t$evalue\n";
	######



	###STEP5: Extract the aligned portions for the next round###
	open (EXTR, ">", $extract_file) or die "cannot open $extract_file.\n";
	foreach my $gene (sort keys %alignments) {
		my @gene_id = split (/_/, $gene);
		my $seg_strand = "+";
		if ($gene_id[6] > $gene_id[7]) {
			$seg_strand = "-";
		}
		my @seg_range = sort {$a <=> $b} ($gene_id[6], $gene_id[7]);
		$seg_range[0] = $gene_id[4] + $seg_range[0] - $flank_size - 1;
		if ($seg_range[0] < 1) {
			$seg_range[0] = 1;
		}
		$seg_range[1] = $gene_id[4] + $seg_range[1] + $flank_size - 1;
		if ($seg_range[1] > length($alignments{$gene})) {
			$seg_range[1] = length($alignments{$gene});
		}

		my $seg_header = $gene_id[0]."_".$gene_id[1]."_".$gene_id[2]."_".$gene_id[3]."_".$seg_range[0]."-".$seg_range[1]."(".$seg_strand.")";
		my $seg = substr ($alignments{$gene}, $seg_range[0] - 1, $seg_range[1] - $seg_range[0] + 1);

		print EXTR ">$seg_header\n";
		print EXTR "$seg\n";
	}
	close (EXTR);
	######



	###STEP6: Align extracted segments to the profile or on its own###
	if ($alignment =~ m/hmmalign/) {
		my $hmmalign_cmd = &hmmaligner($extract_file, $prof_file, $aln_file);
		print LOG "  CMD:\t$hmmalign_cmd\n";
	}
	elsif ($alignment =~ m/mafft/) {
		my $hit_and_ref_temp = $round_dir."/query_and_hits.fas";
		&fastamerger ($q_file, $extract_file, $hit_and_ref_temp);
		my $hit_and_ref_aln_temp = $round_dir."/query_and_hits_aln.fas";
		my $mafft_cmd = &maffter ($hit_and_ref_temp, $hit_and_ref_aln_temp, "DNA", $threads, $element_name);
		my @hits = split(/>/, `grep '>' $extract_file | perl -pe 'chomp' | sed 's/^>//' | perl -pe 'chomp'`);
		&fastaparser ($hit_and_ref_aln_temp, \@hits, $aln_file);
		system "rm $hit_and_ref_temp";
		system "rm $hit_and_ref_aln_temp";

		print LOG "  CMD:\t$mafft_cmd\n";
	}
	my $weblogo_cmd2 = &weblogger ($aln_file, $aln_logo_file, "logodata");
	print LOG "  CMD:\t$weblogo_cmd2\n";
	######



	###STEP7: Cut the alignment based on the logo data###
	my $cutoff = $nat;
	my $cutoff_bit = $minbit;
	if ($seq_count < 2 && $cutoff_bit > 0.6) {
		$cutoff_bit = 0.6;
		$cutoff = &bit2nat($cutoff_bit);
		print "***WARNING: Minimum bit threshold ($minbit) is too high with only one query sequence. Threshold lowered to $cutoff_bit for this round***\n";
	}

	my @range = @{ &logodata2range ($aln_logo_file, $cutoff, $minweight) };
	print " Minimum information (bit):\t$cutoff_bit\n";
	print LOG " Minimum information (bit):\t$cutoff_bit\n";
	if ($range[0] eq 0 && $range[1] eq 0) {
		system "rm $out/TEMPSUBJECT_*";
		system "cp $overview $out";
		print "NO POSITION WITH SUFFICIENT INFORMATION (process complete)\n";
		print LOG "NO POSITION WITH SUFFICIENT INFORMATION (process complete)\n";
		print LOG "##########\n";
		close (LOG);
		die "##########\n\n";
	}
	else {
		print " Positions kept:\t$range[0]-$range[1]\n";
		print LOG " Positions kept:\t$range[0]-$range[1]\n";
	}
	&fastatrimmer ($aln_file, $range[0], $range[1], $aln_cut_file);
	######



	###STEP8: Test for the next round###
	my @seqlist_out = split (/ /, `cat $aln_cut_file | grep '>' | perl -pe 'chomp' | sed 's/>/ /g' | sed 's/^ //' | perl -pe 'chomp'`);
	my $array_identity_test = &array_is_identical(\@seqlist_in, \@seqlist_out);
	if ($array_identity_test eq "TRUE") {
		$test = "STOP";
	}
	if ($round >= $maxround) {
		system "rm $out/TEMPSUBJECT_*";
		system "cp $overview $out";
		print "SEARCH DID NOT CONVERGE AFTER MAXIMUM NUMBER OF ITERATIONS (process complete)\n";
		print LOG "SEARCH DID NOT CONVERGE AFTER MAXIMUM NUMBER OF ITERATIONS (process complete)\n";
		print LOG "##########\n";
		close (LOG);
		die "##########\n\n";
	}
	######
}
system "rm $out/TEMPSUBJECT_*";
system "cp $overview $out";
print "Search converged!\n";
print LOG "Search converged!\n";

#####STEP9: Calculate likelihood of based on binomial distribution#####
my $trials = int(`grep -v '>' $subject | tr -d "\n" | wc -c | perl -pe 'chomp'`/$maxwindow + 0.5);
if ($strand eq "both") {
	$trials = $trials * 2;
}
my $binomial_p = 1 - Math::CDF::pbinom($final_seq_count - 1, $trials, $final_evalue);

print "\n";
print LOG "\n";
print LOG "ENRICHMENT TEST\n";
print "Nunber of hits: $final_seq_count\n";
print LOG "Nunber of hits:\t$final_seq_count\n";
print "Binomial probability: $binomial_p (x: $final_seq_count, n: $trials, p: $final_evalue)\n";
print LOG "Binomial probability:\t$binomial_p\t(x: $final_seq_count, n: $trials, p: $final_evalue)\n";
print LOG "##########\n";
close (LOG);
print "##########\n\n";
##########



#####SUB#####
#Given two arrays, tests if they are identical or not
sub array_is_identical {
	if (@_ != 2) {
		die "USAGE: module array_is_identical requires 2 arguments: <\@array1> <\@array2>.\n";
	}

	my ($array1, $array2) = @_;
	my @array1 = @{ $array1 };
	my @array2 = @{ $array2 };

	my %hash1 = map { $_ => "EXISTS" } @array1;
	my %hash2 = map { $_ => "EXISTS" } @array2;

	my $test="TRUE";
	foreach my $key1 (keys %hash1) {
		if (!exists $hash2{$key1}) {
			$test="FALSE";
			last;
		}
	}
	if ($test eq "TRUE") {
		foreach my $key2 (keys %hash2) {
			if (!exists $hash1{$key2}) {
				$test="FALSE";
				last;
			}
		}
	}
	return ($test);
}

#Given a database, window size, and slide width, generates a new fragmented database
sub databasesegmenter {
	if (@_ != 4) {
		die "USAGE: module databasesegmenter requires 4 arguments: <\$database_file> <\$windowsize> <\$slidewidth> <\$output_file>.\n";
	}

	my ($database_file, $windowsize, $slidewidth, $output_file) = @_;

	open (IN, "<", "$database_file") or die "cannot open $database_file.\n";
	open (OUT, ">", "$output_file") or die "cannot open $output_file.\n";
	while (my $line = <IN>) {
		if ($line =~ m/^>/) {
			my $header = $line;
			$header =~ s/\r//sig;
			$header =~ s/\n//sig;
			$header =~ s/^>//;

			my $seq = <IN>;
			$seq =~ s/\r//sig;
			$seq =~ s/\n//sig;

			my $seq_length = length ($seq);

			my $tail = "N" x $windowsize;
			for (my $i = 1; $i < $seq_length - $slidewidth; $i = $i + $slidewidth) {
				my $frag = substr ($seq, $i - 1, $windowsize);
				my $frag_length = length ($frag);
				if ($frag_length < $windowsize) {
					$frag = $frag.$tail;
					$frag = substr ($frag, 0, $windowsize);
				}
				my $end = $i - 1 + $windowsize;
				if ($end > $seq_length) {
					$end = $seq_length;
				}

				my $frag_header = $header."_".$i."_".$end;
				my $test = length ($frag);

				if ($frag ne $tail) {
					print OUT ">$frag_header\n";
					print OUT "$frag\n";
				}
			}
		}
		else {
			next;
		}
	}
	close (IN);
	close (OUT);
}

#Given a bit value, converts this to nat
sub bit2nat {
	if (@_ != 1) {
		die "USAGE: module bit2nat requires 1 argument: <\$bit>.\n";
	}

	my ($bit) = @_;
	my $nat = log(2)*$bit;

	return $nat;
}

#Given the population size and the number of draws, returns the number of combinations
sub combination{
	if (@_ != 2) {
		die "USAGE: module combination requires 2 arguments: <\$n> <\$r>.\n";
	}

	my ($n, $r) = @_;

	#check
	if(($n < $r) || ($n < 0) || ($r < 0)){
		die "ERROR: Parameters entered into combination is not applicable.\n";
	}

	my $combination = 1;
	my $temp = 1;

	#calculation
	my $limit = (($n - $r) > $r ? $r : ($n - $r));
	for(my $i = 1; $i <= $limit; $i++){
		$temp = $temp * (($n - ($i - 1)) / $i);
		$combination = $temp;
	}

	return $combination;
}

#Given two FASTA files, merges the two FASTA files
sub fastamerger {
	if (@_ != 3) {
		die "USAGE: module fastaparser requires 3 arguments: <\$fasta_file1> <\$fasta_file2> <\$out_file>.\n";
	}

	my ($fasta_file1, $fasta_file2, $out_file) = @_;

	my %fasta;
	open (IN1, "<", $fasta_file1) or die "cannot open $fasta_file1\n";
	while (my $line = <IN1>) {
		if ($line =~ m/^>/) {
			my $header = $line;
			chomp $header;
			$header =~ s/^>//;
			my $seq = <IN1>;
			chomp $seq;
			$fasta{$header} = $seq;
		}
		else {
			next;
		}
	}
	close (IN1);
	open (IN2, "<", $fasta_file2) or die "cannot open $fasta_file2\n";
	while (my $line = <IN2>) {
		if ($line =~ m/^>/) {
			my $header = $line;
			chomp $header;
			$header =~ s/^>//;
			my $seq = <IN2>;
			chomp $seq;
			$fasta{$header} = $seq;
		}
		else {
			next;
		}
	}
	close (IN2);

	open (OUT, ">", $out_file) or die "cannot open $out_file\n";
	foreach my $header (sort keys %fasta) {
		print OUT ">$header\n";
		print OUT "$fasta{$header}\n";
	}
	close (OUT);
}

#Given a FASTA file, list of headers to keep, and an output file name, parses the FASTA file
sub fastaparser {
	if (@_ != 3) {
		die "USAGE: module fastaparser requires 3 arguments: <\$fasta_file> <\@header_list> <\$out_file>.\n";
	}

	my ($in_file, $headers, $out_file) = @_;
	my @headers = @{ $headers };
	my %headers;
	foreach my $header (@headers) {
		$headers{$header} = "KEEP";
	}

	open (IN, "<", $in_file) or die "cannot open $in_file.\n";
	open (OUT, ">", $out_file) or die "cannot open $out_file.\n";
	while (my $line = <IN>) {
		if ($line =~ m/^>/) {
			my $header = $line;
			$header =~ s/\r//sig;
			$header =~ s/\n//sig;
			$header =~ s/^>//;

			my $seq = <IN>;
			$seq =~ s/\r//sig;
			$seq =~ s/\n//sig;

			if (exists $headers{$header}) {
				print OUT ">$header\n";
				print OUT "$seq\n";
			}
		}
		else {
			next;
		}
	}
	close (IN);
	close (OUT);
}

#Given an aligned FASTA file and the range, trims each sequence
sub fastatrimmer {
	if (@_ != 4) {
		die "USAGE: module hitparser requires 4 arguments: <\$in_file> <\$from> <\$to> <\$out_file>.\n";
	}

	my ($in_file, $from, $to, $out_file) = @_;

	open (IN, "<", $in_file) or die "cannot open $in_file.\n";
	open (OUT, ">", $out_file) or die "cannot open $out_file.\n";
	while (my $line = <IN>) {
		if ($line =~ m/^>/) {
			my $header = $line;
			chomp $header;

			my $seq = <IN>;
			chomp $seq;

			my $segment = substr ($seq, $from - 1, $to - $from + 1);

			print OUT "$header\n";
			print OUT "$segment\n";
		}
	}
	close (IN);
	close (OUT);
}

#Given the nhmmscan output, strand to keep, and e-value threshold, parses the nhmmscan file
sub hitparser {
	if (@_ != 5) {
		die "USAGE: module hitparser requires 5 arguments: <\$hit_file> <\$strand> <\$eval> <\$mode> <\$parsed_file>.\n";
	}

	my ($full_file, $strand, $eval, $mode, $parsed_file) = @_;

	#check options
	if ($strand ne "positive" && $strand ne "negative" && $strand ne "both") {
		die "USAGE: module hitparser needs strand to be positive, negative, or both.\n";
	}
	if ($mode ne "BEST" && $mode ne "ALL") {
		die "USAGE: module hitparser needs mode to be BEST or ALL.\n";
	}

	#get header;
	my $header = `head -n 1 $full_file | perl -pe 'chomp'`;

	#store all hits
	my %ref;
	my %genes;
	my $hitcount = 0;
	open (IN, "<", $full_file) or die "cannot open $full_file";
	while (my $line = <IN>) {
		$hitcount++;
		if ($hitcount < 2) {
			next;
		}

		my $hit = $line;
		$hit =~ s/\r//sig;
		$hit =~ s/\n//sig;

		my @hit = split (/\t/, $hit);
		my $segment_id = $hit[2];
		my $gene_id = (split (/_/, $segment_id))[0]."_".(split(/_/, $segment_id))[1]."_".(split(/_/, $segment_id))[2]."_".(split(/_/, $segment_id))[3];
		my @aln_range = ((split (/_/, $segment_id))[4] + $hit[6] - 1, (split (/_/, $segment_id))[4] + $hit[7] - 1);
		@aln_range = sort {$a <=> $b} @aln_range;
		my $aln_range = join ("-", @aln_range);
		my $strand = $hit[11];
		my $eval = $hit[12];

		my %hit;
		$hit{SEGID} = $segment_id;
		$hit{GENEID} = $gene_id;
		$hit{ALNRANGE} = $aln_range;
		$hit{STRAND} = $strand;
		$hit{EVAL} = $eval;
		$hit{FULL} = $hit;
		$ref{$hitcount} = \%hit;

		$genes{$gene_id} = "EXISTS";
	}
	close (IN);

	#remove any hit with e-value below the threashold
	foreach my $hit (sort {$a <=> $b} keys %ref) {
		my %hit = %{ $ref{$hit} };
		if ($hit{EVAL} > $eval) {
			delete $ref{$hit}
		}
	}

	#remove any hit on the wrong strand
	foreach my $hit (sort {$a <=> $b} keys %ref) {
		my %hit = %{ $ref{$hit} };
		if ($strand eq "positive" && $hit{STRAND} eq "-") {
			delete $ref{$hit};
		}
		elsif ($strand eq "negative" && $hit{STRAND} eq "+") {
			delete $ref{$hit};
		}
	}

	#remove hits with overlapping ranges except for the best hit
	foreach my $gene (sort keys %genes) {
		#select subset of hits to each gene
		my %hit_subset1;
		foreach my $hit (sort {$a <=> $b} keys %ref) {
			my %hit = %{ $ref{$hit} };
			if ($gene eq $hit{GENEID}) {
				$hit_subset1{$hit} = \%hit;
			}
		}

		#skip genes with only one hit
		if (scalar keys %hit_subset1 < 2) {
			next;
		}

		#get a list of loci (ranges) with one or more hit(s)
		my %hit_loci;
		foreach my $hit (sort {$a <=> $b} keys %hit_subset1) {
			my %hit = %{ $ref{$hit} };
			my $aln_head = (split (/-/, $hit{ALNRANGE}))[0];
			my $aln_tail = (split (/-/, $hit{ALNRANGE}))[1];
			for (my $i=$aln_head; $i <= $aln_tail; $i++) {
				$hit_loci{$i} = "HIT";
			}
		}
		my @hit_loci_clusters;
		my $locus_head = "";
		foreach my $locus (sort {$a <=> $b} keys %hit_loci) {
			if (!$locus_head) {
				$locus_head = $locus;
				next;
			}
			elsif (!$hit_loci{$locus + 1}) {
				my $range = $locus_head."-".$locus;
				push (@hit_loci_clusters, $range);
				$locus_head = "";
				next;
			}
		}

		#remove all but the best member of each cluster
		foreach my $cluster (@hit_loci_clusters) {
			my $c_head = (split (/-/, $cluster))[0];
			my $c_tail = (split (/-/, $cluster))[1];

			#collect all hits within a cluster
			my %hit_subset2;
			foreach my $hit (sort {$a <=> $b} keys %hit_subset1) {
				my %hit = %{ $hit_subset1{$hit} };
				my $q_head = (split (/-/, $hit{ALNRANGE}))[0];
				my $q_tail = (split (/-/, $hit{ALNRANGE}))[1];

				if ($c_head <= $q_head && $q_tail <= $c_tail) {
					$hit_subset2{$hit} = \%hit;
				}
			}

			#remove all but the best hit within a cluster
			my $best_hit = "";
			my $best_eval = "";
			foreach my $hit (sort {$a <=> $b} keys %hit_subset2) {
				my %hit = %{ $hit_subset2{$hit} };
				if (!$best_hit) {
					$best_hit = $hit;
					$best_eval = $hit{EVAL};
					next;
				}
				elsif ($hit{EVAL} < $best_eval) {
					delete $ref{$best_hit};
					$best_hit = $hit;
					$best_eval = $hit{EVAL};
					next;
				}
				else {
					delete $ref{$hit};
				}
			}
		}
	}

	#keep only the best hit for each gene
	if ($mode eq "BEST") {
		foreach my $gene (sort keys %genes) {
			#select subset of hits to each gene
			my %hit_subset1;
			foreach my $hit (sort {$a <=> $b} keys %ref) {
				my %hit = %{ $ref{$hit} };
				if ($gene eq $hit{GENEID}) {
					$hit_subset1{$hit} = \%hit;
				}
			}

			#skip genes with only one hit
			if (scalar keys %hit_subset1 < 2) {
				next;
			}

			#remove all but the best hit within a cluster
			my $best_hit = "";
			my $best_eval = "";
			foreach my $hit (sort {$a <=> $b} keys %hit_subset1) {
				my %hit = %{ $hit_subset1{$hit} };
				if (!$best_hit) {
					$best_hit = $hit;
					$best_eval = $hit{EVAL};
					next;
				}
				elsif ($hit{EVAL} < $best_eval) {
					delete $ref{$best_hit};
					$best_hit = $hit;
					$best_eval = $hit{EVAL};
					next;
				}
				else {
					delete $ref{$hit};
				}
			}
		}
	}

	#print out
	open (PARSED, ">", $parsed_file) or die "cannot open $parsed_file.\n";
	print PARSED "$header\n";
	foreach my $hit (sort keys %ref) {
		my %hit = %{ $ref{$hit} };
		print PARSED "$hit{FULL}\n";
	}
	close (PARSED);
}

#Given the nhmmscan output and the reference database, generates a hash of the database with the aligned segments indicated
sub hits2hash {
	if (@_ != 2) {
		die "USAGE: module hits2hash requires 2 arguments: <\$hit_file> <\$ref_database>.\n";
	}

	my ($hit_file, $ref_database) = @_;

	my %ref;
	open (REF, "<", $ref_database) or die "cannot open $ref_database.\n";
	while (my $line = <REF>) {
		if ($line =~ m/^>/) {
			my $header = $line;
			$header =~ s/\r//sig;
			$header =~ s/\n//sig;
			$header =~ s/^>//;
			my $seq = <REF>;
			$seq =~ s/\r//sig;
			$seq =~ s/\n//sig;
			$seq = lc ($seq);
			$ref{$header} = $seq;
		}
		else {
			next;
		}
	}
	close (REF);

	my %hits;
	open (HIT, "<", $hit_file) or die "cannot open $hit_file.\n";
	my $linecount = 0;
	while (my $line = <HIT>) {
		$linecount++;
		if ($linecount < 2) {
			next;
		}
		my $data = $line;
		$data =~ s/\r//sig;
		$data =~ s/\n//sig;
		my @data = split (/\t/, $data);
		my @gene_name = reverse (split (/_/, $data[2]));
		my $hit_id = $data[2]."_".$data[6]."_".$data[7];

		my $frag_tail = shift @gene_name;
		my $frag_head = shift @gene_name;
		my $gene_name = join ("_", reverse @gene_name);

		my @ali_range = sort {$a <=> $b} ($frag_head + $data[6] - 1, $frag_head + $data[7] - 1);
		my $ali_id = join("-", @ali_range)."_".$data[11];

		if (!$ref{$gene_name}) {
			die "ERROR: gene ($gene_name) is not found in the given reference ($ref_database)\n";
		}

		my $seq_mod = substr ($ref{$gene_name}, 0, $ali_range[0] - 1).uc(substr ($ref{$gene_name}, $ali_range[0] - 1, $ali_range[1] - $ali_range[0] + 1)).substr ($ref{$gene_name}, $ali_range[1]);

		$hits{$hit_id} = $seq_mod;
	}
	close (HIT);

	return (\%hits);
}

#Reads an aligned single-lined FASTA file and generates an HMM profile
sub hmmaligner {
	if (@_ != 3) {
		die "USAGE: module hmmaligner requires 3 arguments: <\$sequence_file> <\$hmm_file> <\$aligned_file>.\n";
	}

	#
	#Checking if HMMER is running.
	#
	my $hmmer_check = `which hmmalign`;
	if (!$hmmer_check) {
		die "USAGE: module hmmaligner requires HMMER to be running.\n";
	}

	my ($sequence_file, $hmm_file, $aligned_file) = @_;

	my $align = `hmmalign --dna --informat FASTA --outformat Stockholm --trim $hmm_file $sequence_file | grep -v '#' | grep -v '/' | sed 's/^/ /g' | perl -pe 'chomp' | sed 's/  */ /g' | sed 's/^ //' | perl -pe 'chomp'`;
	my %align = split (/ /, $align);

	open (ALN, ">", $aligned_file) or die "cannot open $aligned_file.\n";
	foreach my $header (sort keys %align) {
		print ALN ">$header\n";
		print ALN "$align{$header}\n";
	}
	close (ALN);

	my $hmmalign_cmd = "hmmalign --dna --trim --informat FASTA --outformat Stockholm $hmm_file $sequence_file";
	return ($hmmalign_cmd);
}

#Reads an aligned single-lined FASTA file and generates an HMM profile
sub hmmbuilder {
	if (@_ != 4) {
		die "USAGE: module hmmbuilder requires 4 arguments: <\$sequence_file> <\$hmm_name> <\$hmm_file> <\$threads>.\n";
	}

	#
	#Checking if HMMER is running.
	#
	my $hmmer_check = `which hmmbuild`;
	if (!$hmmer_check) {
		die "USAGE: module hmmbuilder requires HMMER to be running.\n";
	}

	my ($sequence_file, $hmm_name, $hmm_file, $threads) = @_;
	my $seq_count = `grep '>' $sequence_file | wc -l | sed 's/ //g' | perl -pe 'chomp'`;

	my $log;
	my $hmmbuild_cmd;
	if ($seq_count == 1) {
		$log = `hmmbuild -n $hmm_name --dna --fast --wpb --eent --pnone --singlemx --cpu $threads $hmm_file $sequence_file`;
		$hmmbuild_cmd = "hmmbuild -n $hmm_name --dna --fast --wpb --eent --pnone --singlemx $hmm_file $sequence_file";
	}
	else {
		$log = `hmmbuild -n $hmm_name --dna --fast --wpb --eent --pnone --cpu $threads $hmm_file $sequence_file`;
		$hmmbuild_cmd = "hmmbuild -n $hmm_name --dna --fast --wpb --eent --pnone $hmm_file $sequence_file";
	}

	return ($hmmbuild_cmd);
}

#Reads the output file of weblogo and selects a range that has information above a given threshold
sub logodata2range {
	if (@_ != 3) {
		die "USAGE: module logodata2range requires 3 arguments: <\$weblogo_file> <\$nat_threshold> <\$min_weight>.\n";
	}

	my ($weblogo_file, $nat_threshold, $min_weight) = @_;

	#get field for entropy
	my $entropy_field = (split (/Entropy/, `grep 'Entropy' $weblogo_file | tail -n 1`))[0];
	$entropy_field = (() = $entropy_field =~ /\t/gi) + 1;
	#get field for weight
	my $weight_field = (split (/Weight/, `grep 'Weight' $weblogo_file | tail -n 1`))[0];
	$weight_field = (() = $weight_field =~ /\t/gi) + 1;

	#get hash of position and entropy
	my $positions = `grep -v '#' $weblogo_file | cut -f1,$entropy_field,$weight_field | sed 's/ //g'`;
	chomp $positions;

	my @positions = split (/\n/, $positions);
	my %positions;
	foreach my $position (@positions) {
		my @data = split (/\t/, $position);
		my %data;
		$data{Entropy} = $data[1];
		$data{Weight} = $data[2];
		$positions{$data[0]} = \%data;
	}

	#get positions above threshold
	my @pass;
	foreach my $position (sort {$a <=> $b} keys %positions) {
		my %position = %{ $positions{$position} };
		if ($position{Entropy} > $nat_threshold && $position{Weight} > $min_weight) {
			push (@pass, $position);
		}
	}
	if (@pass == 0) {
		push (@pass, 0);
	}

	#get range that passes the threshold
	my @range = ($pass[0], $pass[-1]);

	return (\@range);
}

#Reads a single-lined FASTA file and runs MAFFT to generate single-lined aligned file
sub maffter {
	if (@_ != 5) {
		die "USAGE: module maffter requires 5 arguments: <\$sequence_file> <\$mafft_file> <\$sequence_type> <\$threads> <\$id>.\n";
	}

	#
	#Checking if MAFFT is running.
	#
	my $mafft_check = `which mafft`;
	if (!$mafft_check) {
		die "USAGE: module maffter requires MAFFT to be running.\n";
	}

	my ($sequence_file, $mafft_file, $sequence_type, $threads, $id) = @_;

	#setting parameters for mafft
	my $seqtype;
	if ($sequence_type =~ m/DNA/) {
		$seqtype = "--nuc";
	}
	elsif ($sequence_type =~ m/AA/) {
		$seqtype = "--amino";
	}
	else {
		die "ERROR: Sequence type needs to be DNA or AA.\n";
	}

	#testing for alignment methods
	my $seq_count = `grep '>' $sequence_file | wc -l | sed 's/ //g' | perl -pe 'chomp'`;

	#mafft conditions
	my $max_iteration = 1000;
	my $retree = 2;
	my $aln_method;
	if ($seq_count < 200) {
		$aln_method = "--localpair";
	}
	else {
		$aln_method = "--6merpair";
	}

	#temp files
	my $temp1 = "MAFFTER_01_".$id.".fas";
	my $temp2 = "MAFFTER_02_".$id.".fas";

	#running MAFFT
	system "mafft $seqtype $aln_method --retree $retree --maxiterate $max_iteration --anysymbol --quiet --thread $threads $sequence_file > $temp1";
	my $mafft_cmd = "mafft $seqtype $aln_method --retree $retree --maxiterate $max_iteration --anysymbol $sequence_file >$mafft_file";

	&singler ($temp1, $temp2);

	open (PRE, "<$temp2") or die "cannot open $temp2.\n";
	open (POST, ">$mafft_file") or die "cannot open $mafft_file.\n";
	while (my $line = <PRE>) {
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		if ($line =~ m/^>/) {
			print POST "$line\n";
			next;
		}
		else {
			$line =~ s/U/X/g;
			print POST "$line\n";
		}
	}
	close (PRE);
	close (POST);

	unlink ($temp1);
	unlink ($temp2);

	my $temp_log;
	$temp_log = "Stop codons\(\*\) are not included in the alignement.\n";
	$temp_log = $temp_log."Selenocysteine(U) is replaced by X.\n";

	return ($mafft_cmd);
}

#Given a FASTA file and a pHMM, runs nhmmscan on each sequence independently
sub nhmmscaner {
	if (@_ != 4) {
		die "USAGE: module nhmmscanner requires 4 arguments: <\$phmm> <\$subject> <\$output> <\$threads>.\n";
	}

	#
	#Checking if HMMER is running.
	#
	my $hmmer_check = `which nhmmscan`;
	if (!$hmmer_check) {
		die "USAGE: module nhmmscaner requires HMMER to be running.\n";
	}

	my ($phmm, $subject, $output, $threads) = @_;
	my $seq_count = `grep '>' $subject | wc -l | sed 's/ //g' | perl -pe 'chomp'`;
	my $id = `grep '^NAME' $phmm | rev | cut -f1 -d' ' | rev | perl -pe 'chomp'`;

	#run hmmsearch for each sequence in subject
	my $tempsub = $id."_temp_sub.fas";
	my $temphit = $id."_temp_hit.txt";
	system "echo 'target name\taccession\tquery name\taccession\thmm from\thmm to\tali from\tali to\tenv from\tenv to\tmodlen\tstrand\tE-value\tscore\tbias\tdescription of target' >> $output";
	open (SUBJECT, "<", "$subject") or die "cannot open $subject";
	my $search_count = 0;
	while (my $line = <SUBJECT>) {
		$search_count++;

		my $header = $line;
		my $seq = <SUBJECT>;

		open (SUB, ">", $tempsub);
		print SUB "$header";
		print SUB "$seq";
		close (SUB);

		#run nhmmscan
		system "nhmmscan -E 0.5 --max --tblout $temphit --noali --cpu $threads $phmm $tempsub > /dev/null";
		system "cat $temphit | grep -v '^#' | sed 's/  */\t/g' >>$output";

		if ($search_count % 1000 == 0) {
			my $percent = 100 * $search_count / $seq_count;
			$percent = sprintf ("%.2f", $percent);

			print "  nhmmscan ${percent}% completed.\n";
		}
	}
	system "rm $tempsub";
	system "rm $temphit";

	my $nhmmscan_cmd = "nhmmscan -E 0.5 --max --tblout $output --noali $phmm $subject(one segment at a time)";
	return ($nhmmscan_cmd);
}

#Reads a multi-lined FASTA file and generates a single-lined FASTA file
sub singler {
	if (@_ != 2) {
		die "USAGE: module singler requires 2 arguments: <\$multi_file> <\$single_file>.\n";
	}

	my ($multi_file, $single_file) = @_;

	open (MULTI, "<$multi_file") or die "cannot open $multi_file.\n";
	open (SINGLE, ">$single_file") or die "cannot open $single_file.\n";
	my $line_count = 0;
	while (my $line = <MULTI>) {
		$line =~ s/\r//sig;
		$line =~ s/\n//sig;
		$line =~ s/ /_/g; #specific to this pipeline and not for general use.
		$line_count++;
		if ($line_count == 1 and $line =~ m/^>/) {
			print SINGLE "$line\n";
			next;
		}

		if ($line =~ m/^>/) {
			print SINGLE "\n$line\n";
			next;
		}
		else {
			#$line = uc($line);
			print SINGLE "$line";
		}
	}
	print SINGLE "\n";
	close (MULTI);
	close (SINGLE);
}

#Geven a FASTA file, output file, and format, generates the logo data as text or as image
sub weblogger {
	if (@_ != 3) {
		die "USAGE: module weblogger requires 3 arguments: <\$fasta_file> <\$output_file> <\$format>.\n";
	}

	#
	#Checking if MAFFT is running.
	#
	my $weblogo_check = `which weblogo`;
	if (!$weblogo_check) {
		die "USAGE: module weblogger requires weblogo to be running.\n";
	}

	my ($fasta_file, $output_file, $format) = @_;

	my $weblogo_cmd;
	if ($format =~ m/pdf/) {
		system "weblogo -f $fasta_file -o $output_file -D fasta -A dna -U bits -E YES -c classic -F pdf";
		$weblogo_cmd = "weblogo -f $fasta_file -o $output_file -D fasta -A dna -U bits -E YES -c classic -F pdf";
	}
	elsif ($format =~ m/logodata/) {
		system "weblogo -f $fasta_file -o $output_file -D fasta -A dna -U nats -F logodata";
		$weblogo_cmd = "weblogo -f $fasta_file -o $output_file -D fasta -A dna -U nats -F logodata";
	}
	else {
		die "ERROR: Format for weblogger must be 'pdf' or 'logodata'.";
	}
	return ($weblogo_cmd);
}
##########
__END__
