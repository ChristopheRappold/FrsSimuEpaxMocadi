#! /usr/bin/perl

#use lib "/u/crappold/frs/epax/usr/local/share/perl";
use lib "./usr/lib/perl5";
use IO::CaptureOutput qw/capture/;

use threads;
use Config;
$Config{useithreads} or die      "Recompilez Perl avec les threads activés pour faire tourner ce programme.";

my ($stdout, $stderr);


open my $fileIn, '<', $ARGV[0] or die "error opening $filename: $!";
#my $allFile = do { local $/ = undef; <$fileIn> };

$nSplit = int($ARGV[1]);

@file_array;
@namefile_array;
for(0 .. $nSplit-1)
{
    (my $namefile_out = $ARGV[0]) =~ s/\.[^.]+$//;
    $namefile_out.="_part".$_.".dat";
    push(@namefile_array,$namefile_out);
    unlink $namefile_out;
    open my $fileOut, ">", $namefile_out or die "error opening $namefile_out : $!";
    push(@file_array,$fileOut);
}

while(my $line=<$fileIn>)
{
    my $nL = $. % $nSplit;
    #print "$nL\n";
    print {$file_array[int($nL)]} $line;
}

foreach $fileO (@file_array)
{
    close $fileO;
}


sub sub1 
{
    my @Parameters = @_;
    print "In thread\n";
    print "Parameters >", join("<>", @Parameters), "<\n";
    
    open my $fileIN, "<", $Parameters[0] or die "error opening $Parametres[0]: $!";
    print scalar <$fileIN>;
    close my $fileIN;

    my $namefileIN = $Parameters[0];
    capture sub 
    {
	system "./SystematicMocadi.pl $namefileIN para";
    } => \$stdout, \$stderr;
    return $stdout;

}

@thr_array;
foreach $namefileO (@namefile_array)
{
    my $thr = threads->new(\&sub1, $namefileO);
    push(@thr_array,$thr);
}

foreach $thr1 (@thr_array)
{
    my @datareturned = $thr1->join;
    print @datareturned; 
}
