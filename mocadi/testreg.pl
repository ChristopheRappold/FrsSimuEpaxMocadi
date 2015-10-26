#! /usr/bin/perl

$th = "5.3400e+02";

$th =~ /(-?\d.\d+)e(.\d+)/;
my ($const2, $expon2) = ($1, $2);
$target_th2 = $const2 * 10 ** $expon2;

#print $target_th2."\n";


open my $fileEpax, '<', $ARGV[0] or die "error opening $fileEpax: $!";

while(my $line=<$fileEpax>)
{
    if($line =~ /^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]*\.[0-9]*) (.+?) (.+?) [0-9]+/)
    {
	my $Af = $1;
	my $Zf = $2;
	my $Ap = $3;
	my $Zp = $4;
	my $At = $5;
	my $Zt = $6;
	my $Rho = $7;
	my $xs = $8;
	my $Ratio = $9;

#	print " line : $Af $Zf $Rho $xs and $Ratio ok !\n";
    }
}

close $fileEpax;

open my $etcRoot, "<", "../RadioNuclides.txt" or die "error !";

%table_ListStableLifeTime=();


while(my $line=<$etcRoot>)
{   if($line=~ /^\s+([0-9]+)-([A-Z][a-z]?)-([0-9]*)\s+\d+\s+\d+\s+\d+\s+[0-9]+\.?[0.9]*?\s+.*?\s+(.*?)\s+/)
    {
	($Zr,$Name,$Ar,$T) = ($1,$2,$3,$4);
	my $Ar2 = int($Ar);
	my $Name2 = "A".$Ar2."Z".$Zr;
	my $T2 = $T*2.99792458000000000e+08;
	if($Ar2<10)
	{
	    print "$Zr $Name $Ar $T\n";
	    print "$T2 $Name2\n";
	}
	if($T == -1)
	{
	    $table_ListStableLifeTime{$Name2} = $T2;
	}
	if($T2 > 5.)
	{
	    $table_ListStableLifeTime{$Name2} = $T2;
	}
    }
    # else
    # {
    # 	if($line=~ /^\s+\d+-/)
    # 	    {
    # 		print $line;
    # 	    }
    # }
}

close $etcRoot;


foreach $keyName (sort keys %table_ListStableLifeTime)
{
    if($keyName =~ /A5/)
    {
	print "$keyName ".$table_ListStableLifeTime{$keyName}." \n";
    }

}


$NameTemp = "A9Z6";
if(exists $table_ListStableLifeTime{$NameTemp})
{
    print "exist $NameTemp !\n";

}
else
{
    print "dont exist $NameTemp !\n";
}
$NameTemp = "A4Z2";
if(exists $table_ListStableLifeTime{$NameTemp})
{
    print "exist $NameTemp !\n";

}
else
{
    print "dont exist $NameTemp !\n";
}
$NameTemp = "A5Z2";
if(exists $table_ListStableLifeTime{$NameTemp})
{
    print "exist $NameTemp !\n";

}
else
{
    print "dont exist $NameTemp!\n";
}


my $Aparasist1 = 2;
$Aparasist1 -= 1;
my $Zparasist1 = 5;
my $NameTemp1 = "A".$Aparasist1."Z".$Zparasist1;
while(!exists $table_ListStableLifeTime{$NameTemp1})
{
    $Aparasist1 += 1;
    $NameTemp1 = "A".$Aparasist1."Z".$Zparasist1;
}

print "$NameTemp1\n";

my $Aparasist2 = 3;
$Aparasist2 -= 1;
my $Zparasist2 = 6;
my $NameTemp2 = "A".$Aparasist2."Z".$Zparasist2;
while(!exists $table_ListStableLifeTime{$NameTemp2})
{
    $Aparasist2 += 1;
    $NameTemp2 = "A".$Aparasist2."Z".$Zparasist2;
}

print "$NameTemp2\n";

my $Aparasist3 = 3;
$Aparasist3 -= 1;
my $Zparasist3 = 7;
my $NameTemp3 = "A".$Aparasist3."Z".$Zparasist3;
while(!exists $table_ListStableLifeTime{$NameTemp3})
{
    $Aparasist3 += 1;
    $NameTemp3 = "A".$Aparasist3."Z".$Zparasist3;
}

print "$NameTemp3 $Aparasist3 $Zparasist3\n";
