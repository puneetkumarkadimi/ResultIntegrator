use warnings;
use strict;

=start
			1 => 'MassWiz',
			2 => 'X!Tandem',
			3 => 'OMSSA',
			4 => 'Sequest',
			5=> 'Mascot',
			6=>'Comet',
			7=>'Myrimatch',
			8=>'MSGF+'
=cut

my $Algono=2;

#my $InFile="./Input/TEST/comet/NEW_QTOF_all.comet.pep_Concatenated_FDR.csv";
#my $InFile="./Input/TEST/masswiz/NEW_QTOF_all.Masswiz.Concatenated_FDR_1590765833.csv";
#my $InFile="./Input/TEST/msgf/NEW_QTOF_all.msgf.Concatenated_FDR_1590765420.csv";
#my $InFile="./Input/TEST/myrimatch/NEW_QTOF_all.Myrimatch.Concatenated_FDR_1590819084.csv";
#my $InFile="./Input/TEST/ommsa/NEW_QTOF_all.ommsa_Concatenated_FDR.csv";
my $InFile="./Input/TEST/tandem/NEW_QTOF_all.tandem.2016_11_21_15_15_33.tandem_Concatenated_FDR.csv";

#my $OutFile="NEW_QTOF_all.comet.pep_spectrum.csv";
#my $OutFile="NEW_QTOF_all.Masswiz.pep_spectrum.csv";
#my $OutFile="NEW_QTOF_all.msgf.pep_spectrum.csv";
#my $OutFile="NEW_QTOF_all.Myrimatch.pep_spectrum.csv";
#my $OutFile="NEW_QTOF_all.ommsa.pep_spectrum.csv";
my $OutFile="NEW_QTOF_all.tandem.pep_spectrum.csv";

my $Header=0;

my ($AlgoName,$ScoreColumn,$ScanColumn,$PepColumn,$RankColumn,$ProtColumn,$ScoreType,$Separator,$ChargeColumn,$ExpMassColumn,$TheoMassColumn,$ModPosColumn,$ModificationColumn)=GetAlgoCols($Algono);

open (File,$InFile);
open(OUT,">$OutFile");

print OUT "Algo,Scanid,Peptide\n";

while(<File>)
{
	chomp($_);
	
	if($Header == 0)
	{
		$Header=1;
		next;
	}
	
	my @Line=split(',',$_);
	
	if(scalar @Line != 0) # Line before ROC begins is blank
	{
	if($Algono == 6 || $Algono == 7)
		{
			$Line[$ScanColumn]=~m/\.(\d+\.\d+\.\d)$/gi;
			$Line[$ScanColumn]=$1;
			#print "$Line[$ScanColumn]"; exit;
			my @scan=split(/\./,$Line[$ScanColumn]);
			#print "@scan"; exit;
			@scan=grep(/\w+/,@scan);
			$scan[0]=~s/^\d//g;
			my $Predegit=$&;
			$scan[1]=~s/^\d//g;
			$Line[$ScanColumn]="0$Predegit\.0$scan[0]\.0$scan[1]\.$scan[2]";
			#print "$Line[$ScanColumn]"; exit;
		}
		else
		{
			$Line[$ScanColumn]=~m/(\d{2}\.\d+\.\d+\.\d)$/gi;
			$Line[$ScanColumn]=$1;
		}
		
	print OUT "$AlgoName,$Line[$ScanColumn],$Line[$PepColumn]\n";
}
else
{
	exit;
}	
}

close OUT;
close File;

sub GetAlgoCols
{
	my $AlgoNumber=shift;
		
		my %Algo=(
			1 => 'MassWiz',
			2 => 'X!Tandem',
			3 => 'OMSSA',
			4 => 'Sequest',
			5=> 'Mascot',
			6=>'Comet',
			7=>'Myrimatch',
			8=>'MSGF+'
		);

		my %ScoreCol=(
		'MassWiz'=>7,
		'OMSSA'=>3,
		'X!Tandem'=>4,
		'MSGF+'=>14,
		'Comet'=>19,
		'Myrimatch'=>16,
		'pepXML'=>0,
	);
	
	my %ScanCol=(
		'MassWiz'=>0,
		'OMSSA'=>1,
		'X!Tandem'=>0,
		'MSGF+'=>0,
		'Comet'=>0,
		'Myrimatch'=>0,
		'pepXML'=>0,
	);
	
	my %PepCol=(
		'MassWiz'=>5,
		'OMSSA'=>2,
		'X!Tandem'=>2,
		'MSGF+'=>4,
		'Comet'=>5,
		'Myrimatch'=>5,
		'pepXML'=>0,
	);
	
	my %RankCol=(
		'MassWiz'=>10,
		'OMSSA'=>14,
		'X!Tandem'=>11,
		'MSGF+'=>17,
		'Comet'=>6,
		'Myrimatch'=>6,
		'pepXML'=>0,
	);

	my %ProtCol=(
		'MassWiz'=>8,
		'OMSSA'=>9,
		'X!Tandem'=>6,
		'MSGF+'=>6,
		'Comet'=>7,
		'Myrimatch'=>7,
		'pepXML'=>0,
	);
	
	my %ScoreTyp=(
		'MassWiz'=>1,
		'OMSSA'=>0,
		'X!Tandem'=>0,
		'MSGF+'=>0,
		'Comet'=>0,
		'Myrimatch'=>1,
		'pepXML'=>0,
	);

	my %FileType=(
		'MassWiz'=>'comma',
		'OMSSA'=>'comma',
		'X!Tandem'=>'xml',
		'MSGF+'=>'tab',
		'Comet'=>'xml',
		'Myrimatch'=>'xml',
		'pepXML'=>'xml',
	);
	
	# Charge, Mass, Theoretical Mass, Mod position,Modification
	
	my %Charge=(
		'MassWiz'=>2,
		'OMSSA'=>11,
		'X!Tandem'=>10,
		'MSGF+'=>9,
		'Comet'=>2,
		'Myrimatch'=>2,
		'pepXML'=>0,
	);
	
	my %ExperimentalMass=(
		'MassWiz'=>3,
		'OMSSA'=>4,
		'X!Tandem'=>7,
		'MSGF+'=>8,
		'Comet'=>1,
		'Myrimatch'=>1,
		'pepXML'=>0,
	);
	
	my %TheoreticalMass=(
		'MassWiz'=>4,
		'OMSSA'=>12,
		'X!Tandem'=>8,
		'MSGF+'=>7,
		'Comet'=>8,
		'Myrimatch'=>8,
		'pepXML'=>0,
	);
	
	my %ModPosition=(
		'MassWiz'=>0,
		'OMSSA'=>0,
		'X!Tandem'=>0,
		'MSGF+'=>0,
		'Comet'=>12,
		'Myrimatch'=>12,
		'pepXML'=>0,
	);
	
	my %Modification=(
		'MassWiz'=>9,
		'OMSSA'=>10,
		'X!Tandem'=>5,
		'MSGF+'=>5,
		'Comet'=>13,
		'Myrimatch'=>13,
		'pepXML'=>0,
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