use warnings;
use strict;
use Integrate;

=start
	# Digit codes for Algo
	1 => 'MassWiz',
	2 => 'X!Tandem',
	3 => 'OMSSA',
	4 => 'MSAmanda',
	5 =>'MSGF+'
	6 => 'Mascot',
=cut

=start
##################################################
# Demo folder for file pairing testing purpose only
my %AlgoDir=(
#6=>'./Input/DEMO/comet',
1=>'./Input/DEMO/masswiz',
8=>'./Input/DEMO/msgf',
#7=>'./Input/DEMO/myrimatch',
3=>'./Input/DEMO/ommsa',
2=>'./Input/DEMO/tandem',
);
################################################
=cut

=start
	# Digit codes for Algo
	1 => 'MassWiz',
	2 => 'X!Tandem',
	3 => 'OMSSA',
	4 => 'MSAmanda',
	5 =>'MSGF+'
	6 => 'Mascot',
=cut

my %AlgoDir=(
1=>'./Input/TEST/masswiz',
2=>'./Input/TEST/tandem',
3=>'./Input/TEST/ommsa',
4=>'./Input/TEST/msamanda',
5=>'./Input/TEST/msgf',
6=>'./Input/TEST/mascot',
);

# Prepare this structure
my %AlgoFiles=(
1=>'././Input/TEST/masswiz/NEW_QTOF_Concatenated_FDR.csv',
2=>'././Input/TEST/tandem/NEW_QTOF_Concatenated_FDR.csv',
3=>'././Input/TEST/ommsa/NEW_QTOF_Concatenated_FDR.csv',
5=>'././Input/TEST/msgf/NEW_QTOF_Concatenated_FDR.csv',
6=>'./Input/TEST/mascot/NEW_QTOF_Concatenated_FDR.csv',
);

#my $Files=MainCallerFDRSCore(\%AlgoFiles,0); # Single files
my $Files=MainCallerFDRSCore(\%AlgoDir,1); # For Multiple files
print "Integrated Files: \n".join("\n",@{$Files})."\n";
exit;