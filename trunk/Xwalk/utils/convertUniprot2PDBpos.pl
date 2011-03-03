#!/usr/bin/perl -w
# Author: Abdullah Kahraman
# Date: 03.03.2011

################################################################################
################################################################################
### Converts UniProt residue numbers into PDB residue numbers.               ###
################################################################################
################################################################################

use strict;
use warnings;

use Getopt::Long;

use LWP::UserAgent;

my (
    # variable for parameters which are read in from commandline
    $help,
    $xlFile,
    $pdbFile,
    $fastaFile
   );
##############################################################################
### read all needed parameters from commandline ##############################

&GetOptions(
    "xl=s"    => \$xlFile,    # file with cross-links in distance file format 
    "pdb=s"   => \$pdbFile,   # PDB file
    "fasta=s" => \$fastaFile, # FASTA file
    "help!"   => \$help,      # print this help
) or die "\nTry \"$0 -h\" for a complete list of options\n\n";

##############################################################################

# help
if ($help) {printHelp(); exit}

################################################################################
# SETTINGS
###############################################################################

my $needle = "~/bin/EMBOSS-6.2.0/emboss/needle";
my $warningMsg = "";
my %aa3 = ("GLY" => "G",
	   "ALA" => "A", 
	   "VAL" => "V",
	   "LEU" => "L",
	   "ILE" => "I",
	   "PHE" => "F",
	   "PRO" => "P",
	   "TRP" => "W",
	   "ASN" => "N",
	   "GLN" => "Q",
	   "MET" => "M",
	   "CYS" => "C",
	   "THR" => "T",
	   "TYR" => "Y",
	   "SER" => "S",
	   "ARG" => "R",
	   "LYS" => "K",
	   "HIS" => "H",
	   "ASP" => "D",
	   "GLU" => "E",
	   "ASX" => "B",
	   "GLX" => "Z");

################################################################################
# SUBROUTINES
###############################################################################

###############################################################################
sub printHelp {
    # prints a help about the using and parameters of this scripts 
    # (execute if user types commandline parameter -h)
    # param:  no paramaters
    # return: no return value

    my (
	$usage,
	$sourceCode,
	@rows,
	$row,
	$option,
	$scriptInfo,
	$example,
       );

    $usage = "$0 -xl 6PGD_ECOLI_xl.txt -pdb 2zyaA.pdb -fasta 6PGD_ECOLI.fasta\n";

    print "\nUsage: " .  $usage . "\n";

    print "Valid options are:\n\n";
    open(MYSELF, "$0") or
      die "Cannot read source code file $0: $!\n";
    $sourceCode .= join "", <MYSELF>;
    close MYSELF;
    $sourceCode =~ s/^.+?\&GetOptions\(\n//s;
    $sourceCode =~ s/\n\).+$//s;
    @rows = split /\n/, $sourceCode;
    foreach $row (@rows){
        $option = $row;
	$option =~ s/\s+\"//g;
	$option =~ s/\"\s.+\#/\t\#/g;
	$option =~ s/=./\t<value> [required]/;
	$option =~ s/:./\t<value> [optional]/;
	$option =~ s/!/\t<non value> [optional]/;

	$row =~ s/^.*//;
	print "\t";
	printf("%-1s%-30s%-30s\n", "-",$option,$row);

    } # end of foreach $row (@rows)
    print "\n";
    print "Options may be abreviated, e.g. -h for --help\n\n";

    $example  = "$0";
}
################################################################################
sub createResiduePairTable{
    my $xlTable = "";

    open(X, $xlFile) or die "Failed to open $xlFile";
    while(<X>){
	next if(/^#/);
	(my $index, my $fileName, my $atom1, my $atom2, my @rest) = split(/\t/);
	(my $resName1, my $resNo1, my $chainId1, my $atomName1) = split(/-/, $atom1);
	(my $resName2, my $resNo2, my $chainId2, my $atomName2) = split(/-/, $atom2);
	
	$xlTable .= "$resNo1\t$resNo2\n";
    }
    close(X);
    return $xlTable;
}
################################################################################
sub mapUniprot2PDBseqNumber{
    (my $xlTable) = @_;

    my $newXlTable = "";
    (my $m1_1, my $m1_2, my $m2_1, my $m2_2) = &makeAlignment();
    my %aln2PDB = %$m1_1;
    my %aln2uniprot = %$m2_1;
    my %pdb2aln = %$m1_2;
    my %uniprot2aln = %$m2_2;
    my %mapPDBnum = %{&getPDBresNumbering()};

    my $newDistFile = "";
    open(X, $xlFile) or die "Failed to open $xlFile";
    while(<X>){
	chomp($_);
	next if(/^#/);
	(my $index, my $fileName, my $atom1, my $atom2, my @rest) = split(/\t/);
	(my $resName1, my $resNo1, my $chainId1, my $atomName1) = split(/-/, $atom1);
	(my $resName2, my $resNo2, my $chainId2, my $atomName2) = split(/-/, $atom2);

	my $alnPos1 = $uniprot2aln{$resNo1};
	my $alnPos2 = $uniprot2aln{$resNo2};
	if($aln2PDB{$alnPos1} != -1 && $aln2PDB{$alnPos2} != -1){
	    $resNo1 = $mapPDBnum{$aln2PDB{$alnPos1}};
	    $resNo2 = $mapPDBnum{$aln2PDB{$alnPos2}};

	    $newDistFile .= "$index\t$fileName\t$resName1-$resNo1-$chainId1-$atomName1\t".
		            "$resName2-$resNo2-$chainId2-$atomName2";
	    foreach my $t (@rest){
		$newDistFile .= "\t$t";
	    }
	    $newDistFile .= "\n";
	}
	else {
	    if($aln2PDB{$alnPos1} == -1 && $aln2PDB{$alnPos2} == -1){
		$warningMsg .= "Both UniProt residue numbers $alnPos1 and $alnPos2 could ".
	                       "not be found in the PDB structure.\n";
	    }
	    elsif($aln2PDB{$alnPos1} == -1){
		$warningMsg .= "The UniProt residue number $alnPos1 could ".
                               "not be found in the PDB structure.\n";
	    }
	    elsif($aln2PDB{$alnPos2} == -1){
		$warningMsg .= "The UniProt residue number $alnPos2 could ".
                               "not be found in the PDB structure.\n";
	    }
	}	
    }
    close(X);
    return $newDistFile;
}
################################################################################
sub makeAlignment(){
    my $fasta1 = &createFastaFile4PDBinput();
    my $fasta2 = $fastaFile;
    my $aln = $pdbFile;
    $aln =~ s/.*\///;
    $aln =~ s/(.*)\..*/$1.aln/;
    my $command = "$needle $fasta1 $fasta2 -gapopen 10 -gapextend 0.5 $aln";
    print STDERR "$command\n";
    system(`$command`) or 
	die("Error while performing Needleman-Wunsch alignment: $!");

    return &readAlignment($aln);
}
################################################################################
sub createFastaFile4PDBinput(){
    my $id = $pdbFile;
    $id =~ s/.*\///;
    $id =~ s/(.*)\..*/$1/;
    my $fasta = "$id.fasta";
    my @a = split(/\n/,`grep ^ATOM $pdbFile`);
    open(F,">$fasta") or 
        die ("Couldn't create sequence file \"$fasta\" from PDB file: $!");
    print F ">$id.pdb\n";
    my %uniqAA;
    foreach my $l (@a){
	if($l=~/^ATOM/){
	    my $resName = substr($l, 17, 3);
	    my $resNo   = substr($l, 22, 4);
	    my $chainId = substr($l, 21, 1);
	    my $aa = "$resName-$resNo-$chainId";
	    if(!exists $uniqAA{$aa}){
		$uniqAA{$aa} = $aa;
		if(exists $aa3{$resName}){
		    print F $aa3{$resName};
		}
		else{
		    print F "X";
		}
	    }
	}
    }
    close(F);

    return $fasta;
}
################################################################################
sub readAlignment(){
    (my $aln) = @_;

    my $id = $pdbFile;
    $id =~ s/.*\///;
    $id =~ s/(.*)\..*/$1/;

    open(F,$aln) or die("Failed to open alignment file \"$aln\": $!");

    my $h1;
    my $h2;
    while(<F>){
	chomp($_);
	if(/^[A-Za-z0-9]/){
	    if(/^$id.pdb/){
		$_=~s/.*\d\s([A-Za-z\-])/$1/;
		$_=~s/\s+\d+//;
		$h1 .= $_;
	    }
	    else{
		$_=~s/.*\d\s([A-Za-z\-])/$1/;
		$_=~s/\s+\d+//;
		$h2 .= $_;
	    }
	}
    }
    close(F);
    if(length($h1) != length($h2)){
	die("Error while reading alignment. Alignment sequence lengths do not match!\n");
    }

    my %map1_1;
    my %map1_2;
    my %map2_2;
    my %map2_1;
    my $j1=0;
    my $j2=0;
    for(my $i=0; $i < length($h1); $i++){
	if(substr($h1, $i, 1) ne "-"){
	    $j1++;
	    $map1_1{$i+1} = $j1;
	    $map1_2{$j1} = $i+1;
	} else {
	    $map1_1{$i+1} = -1;
	}
	
	if(substr($h2, $i, 1) ne "-"){
	    $j2++; 
	    $map2_1{$i+1} = $j2;
	    $map2_2{$j2} = $i+1;
	} else {
	    $map2_1{$i+1} = -1;
	}
    }
    return(\%map1_1, \%map1_2, \%map2_1, \%map2_2);
}
################################################################################
sub getPDBresNumbering(){
    my @a = split(/\n/,`grep ^ATOM $pdbFile`);
    my %uniqAA;
    my %map;
    foreach my $l (@a){
	if($l=~/^ATOM/){
	    my $resName = substr($l, 17, 3);
	    $resName =~ s/\s//g;
	    my $resNo   = substr($l, 22, 4);
	    $resNo =~ s/\s//g;
	    my $chainId = substr($l, 21, 1);
	    $chainId =~ s/\s//g;
	    my $aa = "$resName-$resNo-$chainId";
	    if(!exists $uniqAA{$aa}){
		$uniqAA{$aa} = $aa;
		$map{keys (%uniqAA)} = $resNo;
	    }
	}
    }
    return \%map;
}
################################################################################
# MAIN #########################################################################
################################################################################

# check whether all parameter are real files
if(defined $xlFile and defined $pdbFile and defined $fastaFile){
    die "$xlFile does not exist\n" unless(-e $xlFile);
    die "$pdbFile does not exist\n" unless(-e $pdbFile);
    die "$fastaFile does not exist\n" unless(-e $fastaFile);
    
    my $xlTable = &createResiduePairTable();

# create fasta file from pdb file
    my $newDistFile = &mapUniprot2PDBseqNumber($xlTable);
    print STDERR $warningMsg;
    print $newDistFile;
    
}
else {
    die "\nTry $0 -help to get a full list of parameters.\n\n";
}

