use warnings;
use strict;
use Term::Spinner;

my $spinner = Term::Spinner->new(clear_on_destruct => 1, output_handle => \*STDOUT,);

print "Processing 1:\t";
for(my $i=0; $i < 100000; $i++)
{
	$spinner->advance();
	#sleep(10);
}
print "\tDone\n";
$spinner->clear();
#$spinner->finish();
print "Processing 2:\t";
for(my $i=0; $i < 100000; $i++)
{
	$spinner->advance();
	#sleep(10);
}
print "\tDone\n";
$spinner->clear();
#$spinner->finish();
print "Processing 3:\t";
for(my $i=0; $i < 100000; $i++)
{
	$spinner->advance();
	#sleep(10);
}
print "\tDone\n";
$spinner->clear();
#$spinner->finish();
#undef $spinner; # clears final spinner output by default.