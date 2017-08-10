#! /usr/bin/perl

my @array_loop;

for(my $i=0; $i<80;$i++)
{
    push(@array_loop,$i*10);
}

(my $header = $ARGV[0]) =~ s/\.[^.]+$//;
print $header."\n";

my $template=$ARGV[1];
my $minEnergy=$ARGV[2];
    
foreach my $rangeEnergy (@array_loop)
{
    my $Energy = $minEnergy+$rangeEnergy;
    print "Energy:".$Energy."\n";
    
    my $newheader = $header."range".$Energy.".in";
    print $newheader."\n";
    system "cp $header.in $newheader\n";
    system "./LoopMocadi.pl $newheader $template ./4-1e/mocadiR-41e $Energy";
}
