use warnings;
use strict;
use Data::Dumper;
use lib './';
use MODULES::FDRCalculator;
use MODULES::PSMFileIO;
use MODULES::ConcatenatedFDR;
use Term::Spinner;

my $spinner = Term::Spinner->new(clear_on_destruct => 1, output_handle => \*STDOUT,);
####################################################################################################################################
# This module integrates the multiple search engine search results , this module takes proteostats output as input for intrgration
# Improving sensitivity in proteome studies by analysis of false discovery rates for multiple search engines.
# Jones, A. R., Siepen, J. A., Hubbard, S. J., and Paton, N. W. (2009) 
# Written by : Puneet Kumar Kadimi, 26/05/2020
# Input: Concatenated FDR performed by Proteostats or Proteome Analyst on multialgorithmic Concatenated search with decoy proteins tagged as DECOY_ or decoy_
# Output: Integrated results file from multiple algorithms with Combined FDR score
#####################################################################################################################################

=start
Steps to execute:

1. Proteostats reads concatenated search result in a matrix having Q value for each file (including decoy)
2. Integrator Module 1 reads Proteostats output from step 1 calculates FDR score and reports for each file (including decoy)
3. Integrator Module 2 reads all possible search engines search results PSM pairs sums up FDR score makes Combined FDR Score reports combined result (including decoy)
4. Proteostats reads Integrator result from step 3. performs FDR on Combined FDR score reports Q Value for combined file  (excluding decoy)

Gradient and Inbtercept calculation:

1. The slope or Gradient (m) of a line passing through the two points: P=(x1,y1) and Q=(x2,y2) can be calculated as m=y2−y1/x2−x1
2. The Intercept can be calculated as : the y-intercept is b=y1−m*x1
3. The equation of the line can be written in the form y=mx+b.
=cut


# Step 1: Reading each files provided in each folder and calculating FDR score based on QValue
# Steps One expects getting a hash of {Algo name} => 'Directory address' which will be given via UI

#my $OutPutFiles=FDRSCoreForEachFile(\%AlgoDir,1); # This hash contains algo number and folder name for that algo as key => value , 1 is for directory mode
#my $OutPutFiles=FDRSCoreForEachFile(\%AlgoDirFiles,0); # This hash contains algo number and individual file for that algo as key => value , 0 is file single file mode 

# TODO: 
# 1: Read each algo dir and store their files in hash as {File name}=[Algo1,Algo2,Algon]
# Abondoned: 2: Read each algo DIR assume a table as Row: Set 1 to Set n, Column: Algo 1 to Algo n (Max 6) +1 (Total Files)
# 2: Prepare for each file from 1: in structure {Algo Number}="File with full path" and pass it to AFSandFDR

# my $MappedFiles=FileMapper(\%AlgoDir);
# print Dumper($MappedFiles);

sub MainCallerFDRSCore
{
	my $Hash=shift;
	my $Type=shift;
	
	my $OutPutHash;
	
	my @OutPutFile;
	
	if($Type == 1)
	{
		# This is hash containing Algo => Directory
		#print Dumper($Hash); exit;
		print "Step 1: FDR Score calculation..\t ";
		$OutPutHash=FDRSCoreForEachFile($Hash,1); # $OutPutHash contains data s {Algo No}=[File 1, File 2 .... File N] currently not in use the out put variable
		print "Done..\n";
		
		print "Step 2: Pairing Files..\t ";
		# Taking in the above hash and rearranging in {File name}[Algo 1, Algo 2, Algo N]
		$OutPutHash=FileMapper($Hash);
		#print Dumper($OutPutHash); exit;
		print "Done..\n";
		$spinner->clear();
		
		print "Step 3: AFS Calculation and FDR Estimation starts..\n";
		# Running for each File and rearranging it in hash as Algo => File name
		foreach my $File(keys %{$OutPutHash})
		{
			my %FDRScoreFilevsAlgo;
			
			foreach my $AlgoNo(@{${$OutPutHash}{$File}})
			{
				my $Path=${$Hash}{$AlgoNo};
				$FDRScoreFilevsAlgo{$AlgoNo}="$Path/$File";
			}
			
			# print Dumper(\%FDRScoreFilevsAlgo); exit;
			
			if(scalar(keys(%FDRScoreFilevsAlgo)) >= 2)
			{
				# Pass here for processing and check if it has atleat 2 algo else skip and move to
				# This is hash containing Algo => File name
				my $IntegratedFile=StartFDRSCoreCalculation(\%FDRScoreFilevsAlgo);
				push(@OutPutFile,$IntegratedFile);
			}
			else
			{
				print STDOUT "Can not combine less than 2 Files.. looking next set.. !!!\n";
				next;
			}
		}
		print "Step 3: AFS Calculation and FDR Estimation ends..\n";
	}
	else
	{
		if(scalar(keys(%{$Hash})) >= 2)
		{
			print "Step 1: FDR Score calculation..\t ";
			# This is hash containing Algo => File name
			$OutPutHash=FDRSCoreForEachFile($Hash,0);
			print "Done..\n";
			
			print "Step 2: AFS Calculation and FDR Estimation starts..\t ";
			my $IntegratedFile=StartFDRSCoreCalculation($OutPutHash);
			print "Step 2: AFS Calculation and FDR Estimation ends..\t";
			push(@OutPutFile,$IntegratedFile);
		}
		else
		{
			print STDOUT "Can not combine less than 2 Files.. stopping.. !!!\n";
		}
	}
	# Returning all of the output files at once
	return (\@OutPutFile);
}

=start
sub HashMapperFileToAlgo
{
my $data=shift;
my %AlgoDir=%{$data};
my %FDRScoreFilevsAlgo;

foreach my $Algo(keys %AlgoDir)
{
	foreach(@{$AlgoDir{$Algo}})
	{
		if(exists($FDRScoreFilevsAlgo{$_}))
			{
				push(@{$FDRScoreFilevsAlgo{$_}},$Algo);
			}
			else
			{
				# Assigning an array variable anynomously
				$FDRScoreFilevsAlgo{$_}=[];
				push(@{$FDRScoreFilevsAlgo{$_}},$Algo);
			}
	}
}
return (\%FDRScoreFilevsAlgo);
}
=cut

sub FileMapper
{
# Reads each algo dir arranges data in hash as and returns {File name}=[Algo1,Algo2,Algon]
my $data=shift;
my %AlgoDir=%{$data};
my %FDRScoreFilevsAlgo;

foreach my $Algo(keys %AlgoDir)
{
	# This subroutine Reads each algo directory and stores files in structure as {File name}=[Algo 1,Algo 2...Algo n]
	my $Dir=$AlgoDir{$Algo};
	opendir(DIR,$Dir);
	my @FilesOfDirectory=readdir(DIR);
	@FilesOfDirectory=grep(/\.FDRScore\.csv$/g,@FilesOfDirectory);
	close DIR;

	foreach(@FilesOfDirectory)
	{
		if($_=~m/IntegratedFDRScore/g)
		{
			next;
		}
		else
		{
			if(exists($FDRScoreFilevsAlgo{$_}))
			{
				push(@{$FDRScoreFilevsAlgo{$_}},$Algo);
			}
			else
			{
				# Assigning an array variable anynomously
				$FDRScoreFilevsAlgo{$_}=[];
				push(@{$FDRScoreFilevsAlgo{$_}},$Algo);
			}
		}
	}
}
return (\%FDRScoreFilevsAlgo);
}

# Step 2 : For each set of files reading in combining multiple PSM , taking root of scores, then FDR and Finally Combined FDR Score
sub StartFDRSCoreCalculation
{
	my $data=shift;
	my %AlgoDirFiles=%{$data}; # Expects: {Algo no} = File name with full path and pass to AFSandFDR
	my $IntegratedFDRFile=CalculateAFS(\%AlgoDirFiles);
	my $FDRoutfile=ConAscNoRankSepFDR($IntegratedFDRFile);
	print "Processing : $FDRoutfile begins Combined FDR Score calculation ..\n";
	my ($Matrix,$Header,$ModuDigit)=CalculateFDRScore($IntegratedFDRFile,9); # 9 => Integrated score multiple files FDRScore averaged	
	my $CombinedFDRSCoreFile=WriteResults($IntegratedFDRFile,$Matrix,'NOUNMATCHMATRIX',$Header,4); # Since there is no unmatch matrix its passed as NOUNMATCHMATRIX its only for Integration of files 
	print "Processing : $FDRoutfile ends Combined FDR Score calculation ..\n";
	return($CombinedFDRSCoreFile);
}


#exit;
# 1.Start finding pairs of files from each algo and arrange them in matrix as Set vs Algo names and at the end the total count of file names for each set for eaxmple 1 to 6
# 2.Run for each pair of files found o

#n different algo find pairs of PSM data can be stored as in hash {Spectrum, Peptide}[[whole line],[fdr score]]
# 3.Multiply FDR Score for same set of PSM from different algo and write a single file with comman format
# 4.Files which did not even have a single pair file from different algo search should be written in text file that these were not a have pair file to combine (User can do so manually for single option)
# 5.Read the file generated at step 3 perform concatenated FDR since it has not fixed coulmn for Spectrum, Peptide , Protein, Charge, Mass, Theoretical, Mass, Mod position,Modification, FDRSCore
# 6. Read the file generated at step 5 calculate Combined FDR Score using FDR score and new QValue

# For now write a subroutine which takes 6 algo 6 files and combines them

sub CalculateAFS
{

	# - Write the result file and perform FDR
	# - Read the file and calculate FDR Score
	# - Write the integrated result 

	my $AlgoDirFiles=shift;
	my %AlgoDirFiles=%{$AlgoDirFiles};
	my %HOHO=();
	my $SetFileName;

foreach my $AlgoNo(keys %AlgoDirFiles)
{
	# Now read the file in a hash and store it in hash
	my ($DataHash,$PsmCount)=ReadinHashFDRScoreFile($AlgoNo,$AlgoDirFiles{$AlgoNo});
	$HOHO{$PsmCount}=$DataHash;
	$SetFileName=$AlgoDirFiles{$AlgoNo}; 
	# {scan\tpeptide}=[scan,peptide,protein,experimentalmass,theoreticalmass,charge,modification,algoname,rank,fdrscore]
}
	# Here taking only the file name for ex if its "./folder1/folder2/folder3/file.csv" will take only file.csv"
	$SetFileName=~m/^\.\/.+\/(.+)/gi;
	$SetFileName=$1;
	print "Processing : $SetFileName\t";
	my @ARHO=HashreftoArray(\%HOHO);
	my ($Match,$Unmatch)=CalculateAFSforAll(@ARHO);
	my $IntegrateFileHeader="Scan,Peptide,Protein,Experimental-Mass,Theoretical-Mass,Charge,Modification,Algoname,FDRScore,AlgoCount";
	my $Integratedfile=WriteResults($SetFileName,$Match,$Unmatch,$IntegrateFileHeader,3); # $SetFileName is a name picked from any of the input set randomly
	print "AFS Done..\t";
	return ($Integratedfile);
}

sub ConAscNoRankSepFDR
{
	print "FDR begins..";
	# This is a customized set of routines from Proteostats which reads concatenated files with no Rank
	# performs FDR in Ascending order with Separate methods
	# Filters and the output with decoy tags and ranks
	# ReadTextconNoRank: integratedfile, score, scan, pep, rank,separator,prot,dectag
	# filterConAscDecTagNoRank: integratedfile, rocref, score, scan, outfile, scan, separator, prot,dectag
	# Input headers: Scan	Peptide	Protein	Experimental-Mass	Theoretical-Mass	Charge	Modification	Algoname	FDRScore	AlgoCount
	my $Integratedfile=shift;
	my ($tarref,$decref)=ReadTEXTconNoRank($Integratedfile,8,0,1,0,'comma',2,'DECOY_');
	my $ROCref=ConFDR_cutoff_asc_fdrscore($tarref,$decref);
	my $FDRoutfile=filterConAscDecTagNoRank($Integratedfile,$ROCref,8,0,$Integratedfile,0,'comma',2,'DECOY_'); # With Decoy
	print "FDR ends..\n";
	return($FDRoutfile);
}

sub HashreftoArray
{
	# Sorting the hash by values in DSC order and pusihing it in array
	my $Hash=shift;
	my @ARHO=();
	foreach (sort{ $b <=> $a }keys(%{$Hash}))
	{
		push(@ARHO,${$Hash}{$_});
	}
	return (@ARHO);
}

sub CalculateAFSforAll
{
# Data Structure: {scan\tpeptide}=[scan,peptide,protein,experimentalmass,theoreticalmass,charge,modification,algoname,fdrscore,matchflag]
my @ARHO=@_;
my $TotalAlgo=scalar(@ARHO);
my $LargestHash=shift(@ARHO);

my %Match=();
my %Unmatch=();

if($TotalAlgo == 2) # Combine 2 algo
{
	foreach my $psm (keys %{$LargestHash})
	{
		my (@modification,@algoname,@fdrscore);
		my $MatchFlag=0;
		
		if(exists(${$ARHO[0]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[0]}{$psm}[6]);
			push(@algoname,${$ARHO[0]}{$psm}[7]);
			push(@fdrscore,${$ARHO[0]}{$psm}[9]);
			${$ARHO[0]}{$psm}[9]++; # This algo spectrum got matched
			$Match{$psm}[10]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[0]}{$psm};
		}
		
		
		if($MatchFlag > 0)
		{
			# Get sum up of FDR score, Algo names, Modifications assign it to match
			push(@modification,${$LargestHash}{$psm}[6]);
			push(@algoname,${$LargestHash}{$psm}[7]);
			push(@fdrscore,${$LargestHash}{$psm}[8]);
			my ($joined_modification,$joined_algoname,$joined_fdrscore)=Sumup(\@modification,\@algoname,\@fdrscore);
			$Match{$psm}=${$LargestHash}{$psm};
			# Re assigning the information
			$Match{$psm}[6]=$joined_modification;
			$Match{$psm}[7]=$joined_algoname;
			$Match{$psm}[8]=$joined_fdrscore;
			$Match{$psm}[9]=scalar(@algoname); # Assigning the total algo counts here
		}
		else
		{
			# There is no match for the current psm hence storing it in unmatch
			$Unmatch{$psm}=${$LargestHash}{$psm};
		}
	}

	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[0]})
	{
		if(${$ARHO[0]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[0]}{$_};
		}
	}

}
elsif($TotalAlgo == 3) # Combine 3 algo
{
	foreach my $psm (keys %{$LargestHash})
	{
		my (@modification,@algoname,@fdrscore);
		my $MatchFlag=0;
		
		if(exists(${$ARHO[0]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[0]}{$psm}[6]);
			push(@algoname,${$ARHO[0]}{$psm}[7]);
			push(@fdrscore,${$ARHO[0]}{$psm}[8]);
			${$ARHO[0]}{$psm}[9]++; # This algo spectrum got matched
			$Match{$psm}[9]++; # increasing the number to reflect how many times it matched with different algo
			#$Match{$psm}=${$ARHO[0]}{$psm};
		}
		
		if(exists(${$ARHO[1]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[1]}{$psm}[6]);
			push(@algoname,${$ARHO[1]}{$psm}[7]);
			push(@fdrscore,${$ARHO[1]}{$psm}[8]);
			${$ARHO[1]}{$psm}[9]++;
			$Match{$psm}[9]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[1]}{$psm};
		}
		
		
		if($MatchFlag > 0)
		{
			# Get sum up of FDR score, Algo names, Modifications assign it to match
			push(@modification,${$LargestHash}{$psm}[6]);
			push(@algoname,${$LargestHash}{$psm}[7]);
			push(@fdrscore,${$LargestHash}{$psm}[8]);
			my ($joined_modification,$joined_algoname,$joined_fdrscore)=Sumup(\@modification,\@algoname,\@fdrscore);
			$Match{$psm}=${$LargestHash}{$psm};
			# Re assigning the information
			$Match{$psm}[6]=$joined_modification;
			$Match{$psm}[7]=$joined_algoname;
			$Match{$psm}[8]=$joined_fdrscore;
			$Match{$psm}[9]=scalar(@algoname); # Assigning the total algo counts here
		}
		else
		{
			# There is no match for the current psm hence storing it in unmatch
			$Unmatch{$psm}=${$LargestHash}{$psm};
		}
	}

	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[0]})
	{
		if(${$ARHO[0]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[0]}{$_};
		}
	}
	
	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[1]})
	{
		if(${$ARHO[1]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[1]}{$_};
		}
	}

}
elsif($TotalAlgo == 4) # Combine 4 algo
{
	#print "241: $TotalAlgo"; exit;
	foreach my $psm (keys %{$LargestHash})
	{
		my (@modification,@algoname,@fdrscore);
		my $MatchFlag=0;

		if(exists(${$ARHO[0]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[0]}{$psm}[6]);
			push(@algoname,${$ARHO[0]}{$psm}[7]);
			push(@fdrscore,${$ARHO[0]}{$psm}[8]);
			${$ARHO[0]}{$psm}[9]++; # This algo spectrum got matched
			$Match{$psm}[9]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[0]}{$psm};
		}
		
		if(exists(${$ARHO[1]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[1]}{$psm}[6]);
			push(@algoname,${$ARHO[1]}{$psm}[7]);
			push(@fdrscore,${$ARHO[1]}{$psm}[8]);
			${$ARHO[1]}{$psm}[9]++;
			$Match{$psm}[9]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[1]}{$psm};
		}
		
		if(exists(${$ARHO[2]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[2]}{$psm}[6]);
			push(@algoname,${$ARHO[2]}{$psm}[7]);
			push(@fdrscore,${$ARHO[2]}{$psm}[8]);
			${$ARHO[2]}{$psm}[9]++;
			$Match{$psm}[9]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[2]}{$psm};
		}
		
		if($MatchFlag > 0)
		{
			# Get sum up of FDR score, Algo names, Modifications assign it to match
			push(@modification,${$LargestHash}{$psm}[6]);
			push(@algoname,${$LargestHash}{$psm}[7]);
			push(@fdrscore,${$LargestHash}{$psm}[8]);
			my ($joined_modification,$joined_algoname,$joined_fdrscore)=Sumup(\@modification,\@algoname,\@fdrscore);
			
			$Match{$psm}=${$LargestHash}{$psm};
			# Re assigning the information
			$Match{$psm}[6]=$joined_modification;
			$Match{$psm}[7]=$joined_algoname;
			$Match{$psm}[8]=$joined_fdrscore;
			$Match{$psm}[9]=scalar(@algoname); # Assigning the total algo counts here
		}
		else
		{
			# There is no match for the current psm hence storing it in unmatch
			$Unmatch{$psm}=${$LargestHash}{$psm};
		}
	}

	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[0]})
	{
		if(${$ARHO[0]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[0]}{$_};
		}
	}
	
	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[1]})
	{
		if(${$ARHO[1]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[1]}{$_};
		}
	}
	
	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[2]})
	{
		if(${$ARHO[2]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[2]}{$_};
		}
	}

}
elsif($TotalAlgo == 5) # Combine 5 algo
{
	foreach my $psm (keys %{$LargestHash})
	{
		my (@modification,@algoname,@fdrscore);
		my $MatchFlag=0;
		
		if(exists(${$ARHO[0]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[0]}{$psm}[6]);
			push(@algoname,${$ARHO[0]}{$psm}[7]);
			push(@fdrscore,${$ARHO[0]}{$psm}[8]);
			${$ARHO[0]}{$psm}[9]++; # This algo spectrum got matched
			$Match{$psm}[9]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[0]}{$psm};
		}
		
		if(exists(${$ARHO[1]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[1]}{$psm}[6]);
			push(@algoname,${$ARHO[1]}{$psm}[7]);
			push(@fdrscore,${$ARHO[1]}{$psm}[8]);
			${$ARHO[1]}{$psm}[9]++;
			$Match{$psm}[9]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[1]}{$psm};
		}
		
		if(exists(${$ARHO[2]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[2]}{$psm}[6]);
			push(@algoname,${$ARHO[2]}{$psm}[7]);
			push(@fdrscore,${$ARHO[2]}{$psm}[8]);
			${$ARHO[2]}{$psm}[9]++;
			$Match{$psm}[9]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[2]}{$psm};
		}
		
		if(exists(${$ARHO[3]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[3]}{$psm}[6]);
			push(@algoname,${$ARHO[3]}{$psm}[7]);
			push(@fdrscore,${$ARHO[3]}{$psm}[8]);
			${$ARHO[3]}{$psm}[9]++;
			$Match{$psm}[9]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[3]}{$psm};
		}
		
		if($MatchFlag > 0)
		{
			# Get sum up of FDR score, Algo names, Modifications assign it to match
			push(@modification,${$LargestHash}{$psm}[6]);
			push(@algoname,${$LargestHash}{$psm}[7]);
			push(@fdrscore,${$LargestHash}{$psm}[8]);
			my ($joined_modification,$joined_algoname,$joined_fdrscore)=Sumup(\@modification,\@algoname,\@fdrscore);
			
			$Match{$psm}=${$LargestHash}{$psm};
			# Re assigning the information
			$Match{$psm}[6]=$joined_modification;
			$Match{$psm}[7]=$joined_algoname;
			$Match{$psm}[8]=$joined_fdrscore;
			$Match{$psm}[9]=scalar(@algoname); # Assigning the total algo counts here
		}
		else
		{
			# There is no match for the current psm hence storing it in unmatch
			$Unmatch{$psm}=${$LargestHash}{$psm};
		}
	}

	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[0]})
	{
		if(${$ARHO[0]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[0]}{$_};
		}
	}
	
	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[1]})
	{
		if(${$ARHO[1]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[1]}{$_};
		}
	}
	
	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[2]})
	{
		if(${$ARHO[2]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[2]}{$_};
		}
	}
	
	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[3]})
	{
		if(${$ARHO[3]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[3]}{$_};
		}
	}
}
elsif($TotalAlgo == 6) # Combine 6 algo
{
	foreach my $psm (keys %{$LargestHash})
	{
		my (@modification,@algoname,@fdrscore);
		my $MatchFlag=0;
		#print "429: $psm"; <>;
		if(exists(${$ARHO[0]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[0]}{$psm}[6]);
			push(@algoname,${$ARHO[0]}{$psm}[7]);
			push(@fdrscore,${$ARHO[0]}{$psm}[8]);
			${$ARHO[0]}{$psm}[9]++; # This algo spectrum got matched
			#$Match{$psm}[10]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[0]}{$psm};
			#print "439:Comet: $psm => $Match{$psm}\n"; <>
		}
		
		if(exists(${$ARHO[1]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[1]}{$psm}[6]);
			push(@algoname,${$ARHO[1]}{$psm}[7]);
			push(@fdrscore,${$ARHO[1]}{$psm}[8]);
			${$ARHO[1]}{$psm}[9]++;
			#$Match{$psm}[10]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[1]}{$psm};
			#print "450:Myrimatch: $psm => $Match{$psm}\n"; <>
		}
		
		if(exists(${$ARHO[2]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[2]}{$psm}[6]);
			push(@algoname,${$ARHO[2]}{$psm}[7]);
			push(@fdrscore,${$ARHO[2]}{$psm}[8]);
			${$ARHO[2]}{$psm}[9]++;
			#$Match{$psm}[10]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[2]}{$psm};
			#print "461:Masswiz: $psm => $Match{$psm}\n"; <>
		}
		
		if(exists(${$ARHO[3]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[3]}{$psm}[6]);
			push(@algoname,${$ARHO[3]}{$psm}[7]);
			push(@fdrscore,${$ARHO[3]}{$psm}[8]);
			${$ARHO[3]}{$psm}[9]++;
			#$Match{$psm}[10]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[3]}{$psm};
			#print "472:MSGF: $psm => $Match{$psm}\n"; <>
		}
		
		if(exists(${$ARHO[4]}{$psm}))
		{
			$MatchFlag++;
			push(@modification,${$ARHO[4]}{$psm}[6]);
			push(@algoname,${$ARHO[4]}{$psm}[7]);
			push(@fdrscore,${$ARHO[4]}{$psm}[8]);
			${$ARHO[4]}{$psm}[9]++;
			#$Match{$psm}[10]++; # increasing the number to reflect how many times it matched
			#$Match{$psm}=${$ARHO[4]}{$psm};
			#print "483:Tandem: $psm => $Match{$psm}\n"; <>
		}
		
		
		if($MatchFlag > 0)
		{
			# Get sum up of FDR score, Algo names, Modifications assign it to match
			push(@modification,${$LargestHash}{$psm}[6]);
			push(@algoname,${$LargestHash}{$psm}[7]);
			push(@fdrscore,${$LargestHash}{$psm}[8]);
			my ($joined_modification,$joined_algoname,$joined_fdrscore)=Sumup(\@modification,\@algoname,\@fdrscore);
			#print "498: @fdrscore";<>;
			# Assign the the psm here which was match
			$Match{$psm}=${$LargestHash}{$psm};
			
			# Re assigning the information
			$Match{$psm}[6]=$joined_modification;
			$Match{$psm}[7]=$joined_algoname;
			$Match{$psm}[8]=$joined_fdrscore;
			$Match{$psm}[9]=scalar(@algoname); # Assigning the total algo counts here
		}
		else
		{
			# Run for each data hash check for which the match flag is > 0 store them here later
			$Unmatch{$psm}=${$LargestHash}{$psm};
		}
	}

	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[0]})
	{
		if(${$ARHO[0]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[0]}{$_};
		}
	}
	
	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[1]})
	{
		if(${$ARHO[1]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[1]}{$_};
		}
	}
	
	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[2]})
	{
		if(${$ARHO[2]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[2]}{$_};
		}
	}
	
	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[3]})
	{
		if(${$ARHO[3]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[3]}{$_};
		}
	}
	
	# Storing unmatched psms to the unmatch hash
	foreach (keys %{$ARHO[4]})
	{
		if(${$ARHO[4]}{$_}[9] == 0)
		{
			$Unmatch{$_}=${$ARHO[4]}{$_};
		}
	}
}
else
{
	# die("Cant combine only single Algo file ");
	return 0; # There was only single file in the gives set
}

#print Dumper(%Match);
return (\%Match,\%Unmatch);
}

sub Sumup
{
	# This subroutine takes in modification, algoname, and fdrscore returns combined string and result
	my ($modification,$algoname,$fdrscore)=@_;
	
	my $joined_algoname=join("\:",@{$algoname});
	my $joined_modification='none';
	
	# Ommiting blanks
	@{$modification}=grep(/\w+/,@{$modification});
	
	if(scalar @{$modification} > 0)
	{
		$joined_modification=join("\:",@{$modification});
	}
	
	my $joined_fdrscore=0;
	my $Counter=0;
	
	for(my $i=0; $i < scalar @{$fdrscore}; $i++)
	{
		$joined_fdrscore*=${$fdrscore}[$i];
		$Counter++;
	}
	
	$joined_fdrscore=nthRoot($joined_fdrscore,$Counter);
	return ($joined_modification,$joined_algoname,$joined_fdrscore);
}

sub nthRoot
{
	# 218th root of x = x**(1/218)
	my ($fdrscore,$power)=@_;
	$fdrscore=$fdrscore**(1/$power);
	return $fdrscore;
}

sub ReadinHashFDRScoreFile
{
	my($AlgoNo,$InputFile)=@_;
	#print $InputFile; exit;
	my %DataHash=();
	# Get the column information first
	my ($AlgoName,$ScoreColumn,$ScanColumn,$PepColumn,$RankColumn,$ProtColumn,$ScoreType,$Separator,$ChargeColumn,$ExpMassColumn,$TheoMassColumn,$ModPosColumn,$ModificationColumn)=GetAlgoCols($AlgoNo);
	
	open(FILE,$InputFile) or die("FDR Score File not found $!");
	my @File=<FILE>;
	@File=grep(/\w+/,@File);
	shift @File; # Header is no longer at this stage
	close FILE;
	
	# Read line break the line get columns arrange in has and return
	my $Scan="";
	my $PsmCount=0;
	
	for(my $i=0; $i < scalar @File; $i++)
	{
		chomp($File[$i]);
		my @Line=split(/,/,$File[$i]);
		# Retrive the spectrum ids as following : 
		# ommsa+masswiz+msgf+tandem: QT20060328_Den18mix_03.02848.02848.2 => 03.02848.02848.2
		# myrimatch+comet: NEW_QTOF_all.13248.13248.3 = > 01.03248.03248.3 , algo no 6 & 7
		$Scan=$Line[$ScanColumn];
		
			#print "108: $Scan"; exit;
			# QT20060328_Den18mix_05.03319.03319.2 = 05.03319.03319.2
			$Scan=~m/(\d+\.\d+\.\d+\.\d)$/gi;
			$Scan=$1;
		
		# OMMSA have AA in lower case in peptides when modification occurs
		$Line[$PepColumn]=uc($Line[$PepColumn]);
		
		$DataHash{"$Scan\t$Line[$PepColumn]"}[0]=$Scan;
		$DataHash{"$Scan\t$Line[$PepColumn]"}[1]=$Line[$PepColumn];
		$DataHash{"$Scan\t$Line[$PepColumn]"}[2]=$Line[$ProtColumn];
		$DataHash{"$Scan\t$Line[$PepColumn]"}[3]=$Line[$ExpMassColumn];
		$DataHash{"$Scan\t$Line[$PepColumn]"}[4]=$Line[$TheoMassColumn];
		$DataHash{"$Scan\t$Line[$PepColumn]"}[5]=$Line[$ChargeColumn];
		
		# Comet and Myrimatch report Position and modification in separate columns
		if($AlgoNo == 6 || $AlgoNo == 7)
		{
			# Join position and mods with :
			$DataHash{"$Scan\t$Line[$PepColumn]"}[6]=$Line[$ModPosColumn].':'.$Line[$ModificationColumn];
		}
		else
		{
			$DataHash{"$Scan\t$Line[$PepColumn]"}[6]=$Line[$ModificationColumn];
		}
		
		$DataHash{"$Scan\t$Line[$PepColumn]"}[7]=$AlgoName;
		#$DataHash{"$Scan\t$Line[$PepColumn]"}[8]=$Line[$RankColumn];
		$DataHash{"$Scan\t$Line[$PepColumn]"}[8]=$Line[-1]; # storing here the FDR Score
		$DataHash{"$Scan\t$Line[$PepColumn]"}[9]=0; # This column is the indicator 0 if there is no other match 1 if match found 
		$PsmCount++;
	}
	# Data Structure: {scan\tpeptide}=[scan,peptide,protein,experimentalmass,theoreticalmass,charge,modification,algoname,fdrscore,matchflag]
	unlink($InputFile); # Deleting the inputfile since its data is been read in hash and file is no longer requires
	return (\%DataHash,$PsmCount); # Returning the data hash with the total PSM counts in it
}

sub FDRSCoreForEachFile
{
my $AlgoDirorFiles=shift;
my $Mode=shift;

# print "842: ".Dumper($AlgoDirorFiles); exit;

my %OutPutHash=();

if($Mode == 1) # Algo and Directory is being passed
{
	foreach my $Algo(keys %{$AlgoDirorFiles})
	{
		# Reading directory algo wise
		my $Dir=${$AlgoDirorFiles}{$Algo};
		opendir(DIR,$Dir);
		my @Directory=readdir(DIR);
		@Directory=grep(/\_FDR/g,@Directory);
		close DIR;
		#print "834: @Directory"; exit;
		# Processing Algo wise each file calculating FDR Score and storing the output name in hash

		$OutPutHash{$Algo}=[];
		
		for(my $i=0; $i < scalar @Directory; $i++)
		{
			my ($Matrix,$Header,$Module)=CalculateFDRScore("./$Dir/$Directory[$i]",$Algo);
			my $OutPutFile=WriteResults("./$Dir/$Directory[$i]",$Matrix,'NOUNMATCHMATRIX',$Header,1);
			$OutPutHash{$Algo}[$i]=$OutPutFile;
		}
	}
}
else # Algo and file name is being passed
{
	foreach my $Algo(keys %{$AlgoDirorFiles})
	{
		my ($Matrix,$Header,$Module)=CalculateFDRScore(${$AlgoDirorFiles}{$Algo},$Algo);
		my $OutPutFile=WriteResults(${$AlgoDirorFiles}{$Algo},$Matrix,'NOUNMATCHMATRIX',$Header,1);
		$OutPutHash{$Algo}=$OutPutFile;
	}
}
#print Dumper(\%OutPutHash); exit;
return (\%OutPutHash);
}

# Write the result in csv file
#print "Row: ".scalar @{$Matrix}."  ";
#print "Col: ".scalar @{${$Matrix}[0]}."\n";

sub CalculateFDRScore
{
	my ($InFile,$AlgoNo)=@_;
	#print "$InFile,$AlgoNo\n";
	my @Matrix;
	my @MatrixScoreQVal; # This is a paraller matrix storing only , Score,QVal,Gradient,Intercept
	my $Header;
	
	# $AlgoName,$ScoreColumn,$ScanColumn,$PepColumn,$RankColumn,$ProtColumn,$ScoreType,$Separator
	my ($AlgoName,$ScoreColumn,$ScanColumn,$PepColumn,$RankColumn,$ProtColumn,$ScoreType,$Separator,$ChargeColumn,$ExpMassColumn,$TheoMassColumn,$ModPosColumn,$ModificationColumn)=GetAlgoCols($AlgoNo);
		
	# Take these values from column hash later on 
	#my $ScoreColumn=0;
	
	
		if(-e $InFile)
		{
		open (FILE,$InFile) or die("file not found: $!");
		my @FileData=<FILE>;
		close FILE;
		
		$Header=shift @FileData;
		chomp($Header);

		if($Header!~m/q-value/gi)
		{
			die("$InFile is not FDR corrected.");
		}
		
		my $i=0;
		
		for(my $x=0; $x < scalar @FileData; $x++)
		{
			chomp($FileData[$x]);
			
			# Regex below are intended to remove any comma in Protein id and its description
			# These regex is written by keeping in mind to OMSSA results yet is kept open for all of the algo
			# These reg ex take protein ids from text and replace the whole description with only that
			# Proteins can contain decoy tags and [contaminants]
			if($FileData[$x]=~m/\".+\"/)
			{
				if($FileData[$x]=~m/\"(\w+\|\w+\|\w+).+?\"/g)
				{
					$FileData[$x]=~s/\"(\w+\|\w+\|\w+).+?\"/$1/g;
				}
				elsif($FileData[$x]=~m/\"(\[\w+\]\w+\:\w+).+?\"/g)
				{
					$FileData[$x]=~s/\"(\[\w+\]\w+\:\w+).+?\"/$1/g;
				}
				elsif($FileData[$x]=~m/\"(\w+\[\w+\]\w+\:\w+).+?\"/g)
				{
					$FileData[$x]=~s/\"(\w+\[\w+\]\w+\:\w+).+?\"/$1/g;
				}
				elsif($FileData[$x]=~m/\"(\[\w+\]\w+\|\w+\|\w+).+?\"/g)
				{
					$FileData[$x]=~s/\"(\[\w+\]\w+\|\w+\|\w+).+?\"/$1/g;
				}
				elsif($FileData[$x]=~m/\"(\w+\[\w+\]\w+).+?\"/g) # string like "DECOY_[Contaminant]sp|P35748|MYHB_RABIT MYOSIN HEAVY CHAIN, SMOOTH "
				{
					$FileData[$x]=~s/\"(\w+\[\w+\]\w+).+?\"/$1/g;
				}
				
				# if there are multiple modifications for e.g. "oxidation of M:5 ,phosphorylation of T:8 ,deamidation of N and Q:17"
				if($FileData[$x]=~m/\"(\w+\s+of.+?)\"/ig) 
				{
					my $Modification=$1;
					$Modification=~s/\,/\-/ig;
					$FileData[$x]=~s/\"(\w+\s+of.+?)\"/$Modification/ig;
				}
			}
			
			my @FileDataLine=split(",",$FileData[$x]);
			
			if(scalar @FileDataLine != 0) # Line before ROC begins is blank
			{
				for(my $j=0; $j < scalar @FileDataLine; $j++)
				{
					$Matrix[$i][$j]=$FileDataLine[$j];
				}
				
				$MatrixScoreQVal[$i][0]=$FileDataLine[$ScoreColumn]; # Adding Score 
				$MatrixScoreQVal[$i][1]=$FileDataLine[-1]; # Adding QValue
				push(@{$MatrixScoreQVal[$i]},0); # Blank Gradient
				push(@{$MatrixScoreQVal[$i]},0); # Blank Intercept
				$i++;
			}
			else
			{
				@FileData=();
				last;
			}
		}
		}
		else
		{
			die("file not found: $!");
		}
		
		
		
		# Sort by Q value in ascending order
		@Matrix=sort{$a->[-1]<=>$b->[-1]} @Matrix;
		@MatrixScoreQVal=sort{$a->[1]<=>$b->[1]} @MatrixScoreQVal;
		
		# identifyingasetofsteppoints,where the q-value ofthe identification changes (q-value (i) > q-value (i-1))
		my $Row=scalar @MatrixScoreQVal;
		my $Col=scalar @{$MatrixScoreQVal[0]};
				
		my $ScorePrev=$MatrixScoreQVal[0][0]; # Score prev
		my $QvaluePrev=$MatrixScoreQVal[0][1]; # Q Value prev 
		
		for(my $i=1; $i < $Row; $i++) # Begining from 1 to total size -1
		{
			my $ScoreCurr=$MatrixScoreQVal[$i][0];
			my $QvalueCurr=$MatrixScoreQVal[$i][1];
			
			if($QvalueCurr > $QvaluePrev)
			{
				
				my $X1=$ScorePrev; # Previous score
				my $Y1=$QvaluePrev; # Previous Q Value
				my $X2=$ScoreCurr; # Current score
				my $Y2=$QvalueCurr; # Current Q Value
				
				my ($m,$b)=GradientIntercept($X1,$Y1,$X2,$Y2);
				$MatrixScoreQVal[$i][2]=$m; # Storing Gradient
				$MatrixScoreQVal[$i][3]=$b; # Storing Intercept
				
				# Since Gradient and Intercept has been found replacing it with current
				($ScorePrev,$QvaluePrev)=($ScoreCurr,$QvalueCurr);
			}
		}
	
	# Traverse the whole matrix from bottom to up 
	# Check where gradient is present either > 0 or !0 calculate FDR score and assign it for each Score (e value)
	
	# Equation is QVal*Gradient+Intercept, it is to note that while traversing from bottom
	# It is not neccesry that the Gradient and Intercept is present in begining, in this case FDR score will become 0 
	# irrespective of the actual Q Value (suppose if its poor and 1) to avoid that situation first set of gradient and intercept is 1,0
	# It will yield as FDR SCore = Q Value at the begining since there is no gradient and intercept
	my $Gradient=1;
	my $Intercept=0;
	my $FDRScore=0;
	
	for(my $i=$Row-1; $i >= 0 ; $i--) # Traversing in reverse mode
	{
		if($MatrixScoreQVal[$i][2] != 0 ) # There is gradient present
		{
			($Gradient,$Intercept)=($MatrixScoreQVal[$i][2],$MatrixScoreQVal[$i][3]);
		}
		
		# Assigning FDR Score for each Score/E-Value
		$FDRScore=($MatrixScoreQVal[$i][0]*$Gradient+$Intercept);
		
		# Pushing here to the original data matrix
		push(@{$Matrix[$i]},$FDRScore);
	}
	
	# Releasing data from here matrix which was used for calculation only
	@MatrixScoreQVal=();
	return (\@Matrix,$Header,1);
}

sub GradientIntercept
{
	###############################################################
	#			The slope or Gradient (m) can be calculated as m=y2−y1/x2−x1		#
	#		The Intercept can be calculated as : the intercept is b=y1−m*x1			#
	#		The straight line euqation is y= mx+b
	#	gradient being negative or positive doesnt affects the final calculation	# 
	#			if gradient is negative intercept will be positive and vice versa			#
	###############################################################
	my ($X1,$Y1,$X2,$Y2)=@_;
	my $m=(($Y2-$Y1)/($X2-$X1));
	my $b=($Y1-$m*$X1);
	return ($m,$b);
}

=start
sub ScoreCalculator
{
	###############################################################
	#			The slope or Gradient (m) can be calculated as m=y2−y1/x2−x1		#
	#		The Intercept can be calculated as : the y-intercept is b=y1−m*x1		#
	#		Finally the Regression Euqation becomes: y=mx+b									#
	#	gradient being negative or positive doesnt affects the final calculation	# 
	#			if gradient is negative intercept will be positive and vice versa			#
	################################################################
	my ($X1,$Y1,$X2,$Y2)=@_;
	my $m=(($Y2-$Y1)/($X2-$X1));
	my $Intercept=($Y1-$m*$X1);
	my $FDRScore=($X1*$m+$Intercept);
	#print "$X1,$Y1,$X2,$Y2,$m,$Intercept,$FDRScore\n";
	return ($FDRScore);
}
=cut

# This subroutine only writes the output what ever is returned in Matrix/Hash
sub WriteResults
{
	# $SetFileName,$Match,$Unmatch,$IntegrateFileHeader,3
	my ($FileName,$Matrix,$MatrixUnmatch,$Header,$Module)=@_;
	my $OutFile=$FileName;

	if($OutFile=~m/\_FDR/ig)
	{
		$OutFile=~s/\_FDR\_{0,}\d{0,}//ig; # none of the files afterwards will have the word FDR
	}
	
	if($Module == 1)
	{
		$OutFile=~s/\.csv$/\.FDRScore\.csv/g;
		$Header=$Header.',FDRScore';
		open(OUT,">$OutFile") or die ("Unable to open the output file $!: $OutFile");
		print OUT "$Header\n";
		for(my $i=0; $i < scalar @{$Matrix}; $i++)
		{
			print OUT join(',',@{${$Matrix}[$i]});
			print OUT "\n";
		}
		close OUT;
	}
	elsif($Module == 3)
	{
		$OutFile=~s/\.FDRScore\.csv$/\.IntegratedFDRScore\.csv/g;
		$Header=$Header;
		open(OUT,">$OutFile") or die ("Unable to open the output file $!: $OutFile");
		print OUT "$Header\n";
		#print "$Matrix"; <>;
		foreach(keys %{$Matrix}) # writing matched psm
		{
			print OUT join(',',@{${$Matrix}{$_}});
			print OUT "\n";
		}

		foreach(keys %{$MatrixUnmatch}) # writing unmatched psm
		{
			print OUT join(',',@{${$MatrixUnmatch}{$_}});
			print OUT "\n";
		}

		close OUT;
		
	}
	elsif($Module == 4)
	{
		# Remove first averaged fdr score at [8] , remove decoy_ at [2], Produce output with fdr score < 0.5 or what ever is set at [-1]
		# unlink file with extension : IntegratedFDRScore.csv
		# $IntegratedFDRFile,$Matrix,'NOUNMATCHMATRIX',$Header,4
		$OutFile=~s/\.IntegratedFDRScore\.csv$/\.CombinedFDRScore\.csv/g;
		
		#$Header=$Header; # Not using this header using custom header since output is fixed
		
		open(OUT,">$OutFile") or die ("Unable to open the output file $!: $OutFile");

		print OUT "Scan,Peptide,Protein,Experimental-Mass,Theoretical-Mass,Charge,Modification,Algoname,AlgoCount,FDRScore\n";
		for(my $i=0; $i < scalar @{$Matrix}; $i++)
		{
			# Deleting the first averaged FDR score array size becomes 10 from 11
			# delete(${${$Matrix}[$i]}[8]);
			
			if(${${$Matrix}[$i]}[2]=~m/decoy/gi)
			{
				#delete(${${$Matrix}[$i]}[8]);
				next;
			}
=start
0 Scan
1 Peptide
2 Protein
3 Experimental-Mass
4 Theoretical-Mass
5 Charge
6 Modification
7 Algoname
8 FDRScore
9 AlgoCount
10 q-value
11 FDRScore
=cut
			if(${${$Matrix}[$i]}[-1] < 0.05) # Putting threshold of 0.05  user can later on use 0.01 in excel for use
			{
				print OUT ${${$Matrix}[$i]}[0].",".${${$Matrix}[$i]}[1].",".${${$Matrix}[$i]}[2].",".${${$Matrix}[$i]}[3].",".${${$Matrix}[$i]}[4].",";
				print OUT ${${$Matrix}[$i]}[5].",".${${$Matrix}[$i]}[6].",".${${$Matrix}[$i]}[7].",".${${$Matrix}[$i]}[9].",".${${$Matrix}[$i]}[11];
				print OUT "\n";
			}
		}
		
		close OUT;
		
		# Deleting here the file which have extension IntegratedFDRScore
		unlink($FileName);
	}
	return ($OutFile);
}


sub GetAlgoCols
{
	my $AlgoNumber=shift;

=start
	# Digit codes for Algo
	1 => 'MassWiz',
	2 => 'X!Tandem',
	3 => 'OMSSA',
	4 => 'MSAmanda',
	5 =>'MSGF+'
	6 => 'Mascot',
=cut

		my %Algo=(
			1 => 'MassWiz',
			2 => 'X!Tandem',
			3 => 'OMSSA',
			4 =>'MSAmanda',
			5 => 'MSGF+',
			6 => 'Mascot',
			9 => 'Integrated',
		);

		my %ScoreCol=(
		'MassWiz'=>7,
		'OMSSA'=>3,
		'X!Tandem'=>4,
		'MSGF+'=>14,
		'Comet'=>19,
		'Myrimatch'=>16,
		'Integrated'=>8,
		'Mascot'=>12,
		'Sequest'=>0,
		'MSAmanda'=>5,
	);
	
	my %ScanCol=(
		'MassWiz'=>0,
		'OMSSA'=>1,
		'X!Tandem'=>0,
		'MSGF+'=>0,
		'Comet'=>0,
		'Myrimatch'=>0,
		'Integrated'=>0,
		'Mascot'=>0,
		'Sequest'=>0,
		'MSAmanda'=>1,
	);
	
	my %PepCol=(
		'MassWiz'=>5,
		'OMSSA'=>2,
		'X!Tandem'=>2,
		'MSGF+'=>4,
		'Comet'=>5,
		'Myrimatch'=>5,
		'Integrated'=>1,
		'Mascot'=>4,
		'Sequest'=>0,
		'MSAmanda'=>2,
	);
	
	my %RankCol=(
		'MassWiz'=>10,
		'OMSSA'=>14, # There is no rank in OMMSA hence NIST score position is putted here it will allways remain 0
		'X!Tandem'=>11,
		'MSGF+'=>17,
		'Comet'=>6,
		'Myrimatch'=>6,
		'Integrated'=>0,
		'Mascot'=>13,
		'Sequest'=>0,
		'MSAmanda'=>7,
	);

	my %ProtCol=(
		'MassWiz'=>8,
		'OMSSA'=>9,
		'X!Tandem'=>6,
		'MSGF+'=>6,
		'Comet'=>7,
		'Myrimatch'=>7,
		'Integrated'=>2,
		'Mascot'=>6,
		'Sequest'=>0,
		'MSAmanda'=>4,
	);
	
	my %ScoreTyp=(
		'MassWiz'=>1,
		'OMSSA'=>0,
		'X!Tandem'=>0,
		'MSGF+'=>0,
		'Comet'=>0,
		'Myrimatch'=>1,
		'Integrated'=>0,
		'Mascot'=>0,
		'Sequest'=>0,
		'MSAmanda'=>1,
	);

	my %FileType=(
		'MassWiz'=>'comma',
		'OMSSA'=>'comma',
		'X!Tandem'=>'xml',
		'MSGF+'=>'tab',
		'Comet'=>'xml',
		'Myrimatch'=>'xml',
		'Integrated'=>'comma',
		'Mascot'=>'tab',
		'Sequest'=>'comma',
		'MSAmanda'=>'comma',
	);
	
	# Charge, Mass, Theoretical Mass, Mod position,Modification

	my %Charge=(
		'MassWiz'=>2,
		'OMSSA'=>11,
		'X!Tandem'=>10,
		'MSGF+'=>9,
		'Comet'=>2,
		'Myrimatch'=>2,
		'Integrated'=>5,
		'Mascot'=>9,
		'Sequest'=>0,
		'MSAmanda'=>9,
	);
	
	my %ExperimentalMass=(
		'MassWiz'=>3,
		'OMSSA'=>4,
		'X!Tandem'=>7,
		'MSGF+'=>8,
		'Comet'=>1,
		'Myrimatch'=>1,
		'Integrated'=>3,
		'Mascot'=>8,
		'Sequest'=>0,
		'MSAmanda'=>8,
	);
	
	my %TheoreticalMass=(
		'MassWiz'=>4,
		'OMSSA'=>12,
		'X!Tandem'=>8,
		'MSGF+'=>7,
		'Comet'=>8,
		'Myrimatch'=>8,
		'Integrated'=>4,
		'Mascot'=>7,
		'Sequest'=>0,
		'MSAmanda'=>8,
	);
	
	my %ModPosition=(
		'MassWiz'=>0,
		'OMSSA'=>0,
		'X!Tandem'=>0,
		'MSGF+'=>0,
		'Comet'=>12,
		'Myrimatch'=>12,
		'Integrated'=>0,
		'Mascot'=>0,
		'Sequest'=>0,
		'MSAmanda'=>0,
	);
	
	my %Modification=(
		'MassWiz'=>9,
		'OMSSA'=>10,
		'X!Tandem'=>5,
		'MSGF+'=>5,
		'Comet'=>13,
		'Myrimatch'=>13,
		'Integrated'=>6,
		'Mascot'=>5,
		'Sequest'=>0,
		'MSAmanda'=>3,
	);
	
	my $AlgoName=$Algo{$AlgoNumber};
	my $ScoreColumn=$ScoreCol{$AlgoName};
	my $ScanColumn=$ScanCol{$AlgoName};
	my $PepColumn=$PepCol{$AlgoName};
	my $RankColumn=$RankCol{$AlgoName};
	my $ProtColumn=$ProtCol{$AlgoName};
	my $ScoreType=$ScoreTyp{$AlgoName};
	my $Separator=$FileType{$AlgoName};
	
	
	my $ChargeColumn=$Charge{$AlgoName};
	my $ExpMassColumn=$ExperimentalMass{$AlgoName};
	my $TheoMassColumn=$TheoreticalMass{$AlgoName};
	my $ModPosColumn=$ModPosition{$AlgoName};
	my $ModificationColumn=$Modification{$AlgoName};
	
	return($AlgoName,$ScoreColumn,$ScanColumn,$PepColumn,$RankColumn,$ProtColumn,$ScoreType,$Separator,$ChargeColumn,$ExpMassColumn,$TheoMassColumn,$ModPosColumn,$ModificationColumn);
}

1;