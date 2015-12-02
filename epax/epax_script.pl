#! /usr/bin/perl

use lib "/u/crappold/frs/epax/usr/local/share/perl";
use IO::CaptureOutput qw/capture/;
my ($stdout, $stderr);


$At = $ARGV[0];
$Zt = $ARGV[1];

$Ap = $ARGV[2];
$Zp = $ARGV[3];

print "Target : A".$At."Z".$Zt." / Proj : A".$Ap."Z".$Zp."\n";

open IN, "/usr/local/pub/debian5.0/i686/gcc432-11/rootmgr/534-00-rc1/etc/RadioNuclides.txt" or die "impossible d'ouvrir !";

%table_AZ=();
%table_Z_minA=();
%table_ListStable=();

while($line=<IN>)
{
    if(($Zr,$Name,$Ar) = ($line=~ /^.*?([0-9]*)-([A-Z][a-z]*)-([0-9]*).*?[0-9]*/))
	{
	    #print OUTTEMP "$Name A".$Ar."Z".$Zr."\n";
	    $Zr=int($Zr);
	    $Ar=int($Ar);
	    if(!exists $table_Z_minA{$Zr})
	    {
		#print "$Zr $Z $Ar $Name\n";
		$table_Z_minA{$Zr}=$Ar;
	    }
	    $table_ZName{$Zr}=$Name;
	    
	    $Name=$Name.$Ar;
	    $inside="A".$Ar."Z".$Zr;
	    $table_AZ{$inside}=$Name;

	}
}
$table_AZ{"A23Z14"}="23Si";



open IN2, "ListBeam.dat" or die "impossible d'ouvrir !";

while($line=<IN2>)
{
    if(($NameS,$As,$Zs) = ($line=~ /^([0-9]*[A-Z][a-z]?) ([0-9]*) ([0-9]*)$/))
	{
	    #print "$NameS A".$As."Z".$Zs."\n";
	    $Zr=int($Zs);
	    $Ar=int($As);
	    if(!exists $table_ListStable{$Zs})
	    {
		#print "Stable $Zs $As $NameS\n";
		$inside="A".$As."Z".$Zs;
		$table_ListStable{$inside}=$NameS;
	    }
	}
}


# foreach $key (sort keys %table_ListStable)
# {
#     #if(($aa,$zz)=($value=~/A([0-9]*)Z([0-9]*)$/))
#     #{
# #	if($value<10)
# #	    {
# 		print $key.", ".$table_ListStable{$key}."\n";
# #	    }
# #    }
# }
# print "table_AZ\n";
# foreach $key (sort keys %table_AZ)
# {
#     #if(($aa,$zz)=($value=~/A([0-9]*)Z([0-9]*)$/))
#     #{
# #	if($value<10)
# #	    {
# 		print $key.", ".$table_AZ{$key}."\n";
# #	    }
# #    }
# }
# print "table_Z_minA\n"; 
# foreach $key (sort keys %table_Z_minA)
# {
#     my $val = int($table_Z_minA{$key}) - int($key);
#     print $key.", ".$table_Z_minA{$key}." ".$val."\n";
# }



my $MaxProj = "A".$Ap."Z".$Zp;
my $MaxTarget = "A".$At."Z".$Zt;
my $MaxProjName = $table_AZ{$MaxProj};
my $MaxTargetName = $table_AZ{$MaxTarget};

print "name : $MaxProj $MaxTarget $MaxProjName $MaxTargetName \n";
print "$table_AZ{$MaxProj}\n";
#$count = 0;
my $barn = 0.1;
my $InMb = 1000.;

my $nextstop_T = 0;
my $temp_nameT = "A1Z1";

%table_result=();

while($nextstop_T == 0)
{

    my ($tA_t,$tZ_t)=($temp_nameT=~/A([0-9]*)Z([0-9]*)$/);
    #print "Target : $temp_nameT \n";
    my $temp_nameP = "A1Z1";
    my $nextstop_P = 0;

    if($temp_nameT eq $MaxTarget)
    {
	$nextstop_T = 1;
	#print "next target exits !\n";
    }
    
    my $isStableT = 0;
    if(exists $table_ListStable{$temp_nameT})
    {
	$isStableT = 1;
    }
    #print "is stable $isStableT target $temp_nameT beam $temp_nameP\n";
    if($isStableT == 1) 
	{
	    while($nextstop_P == 0)
	    {
		my $nextstep = 0;
		my $counter = 0;
		my $temp_name = "A1Z0";
		
		if($temp_nameP eq $MaxProj)
		{
		    $nextstop_P = 1;
		    #print "next beam stop !\n";
		}
		
		my ($tA_p,$tZ_p)=($temp_nameP=~/A([0-9]*)Z([0-9]*)$/);
		#print "Projectile : $temp_nameP | Target : $temp_nameT\n";
		
		my $isStableP = 0;
		if(exists $table_ListStable{$temp_nameP})
		{
		    $isStableP = 1;
		}
		#print "is stable $isStableP beam $temp_nameP\n";
		if($isStableP == 1) 
		{
		
		    while ($nextstep==0)
		    {
			$counter = $counter + 1;
			if($temp_name eq $temp_nameP)
			{
			    $nextstep = 1;
			}
			if($counter > 1000)
			{
			    $nextstep = 1;
			    #print "stop counter $temp_nameT $temp_nameP \n";
			}
			
			my ($tempA,$tempZ)=($temp_name=~/A([0-9]*)Z([0-9]*)$/);
			#print "Fragment $tempA $tempZ\n";
			
		    
			capture sub {
			    system "./epax_v31_www $tA_p $tZ_p $tA_t $tZ_t $tempA $tempZ"; 
			    #system($command);
			} => \$stdout, \$stderr;
			
			#open IN2, "temp_epax.html" or die "impossible d'ouvrir !";
			($barn2) = ($stdout =~ /sigma = (.*?) b \<br/);
			$barn = 1000*$barn2;
			#print "$barn \n";
			#print "$barn2 \n";
			
			#print "A".$tempA."Z".$tempZ." = $barn mb | $barn2 b \n";
			
			$table_result{$temp_name}{$temp_nameP}{$temp_nameT}=$barn;
			
			if($barn=~/nan/)
			{
			    $nextstep=1;
			}
			
			$tempA=$tempA+1;
			$temp_name="A".$tempA."Z".$tempZ;
		    
			if(! exists $table_AZ{$temp_name})
			{
			    #print "fragment A$tempA Z$tempZ does not exist, next is";
			    $tempZ=$tempZ+1;
			    $tempA=$table_Z_minA{$tempZ};
			    #print "A$tempA Z$tempZ\n";
			    $temp_name="A".$tempA."Z".$tempZ;
			}
		    
			
		    }
		}
		#print "#end [$temp_nameP] [$temp_nameT]\n\n";
		
		$tA_p=$tA_p+1;
		$temp_nameP="A".$tA_p."Z".$tZ_p;
		
		if(! exists $table_AZ{$temp_nameP})
		{
		    #print "projectile A$tA_p Z$tZ_p does not exist, next is";
		    $tZ_p=$tZ_p+1;
		    $tA_p=$table_Z_minA{$tZ_p};
		    #print "A$tA_p Z$tZ_p\n";
		    $temp_nameP="A".$tA_p."Z".$tZ_p;
		}
	
	    }
	}

    $tA_t=$tA_t+1;
    $temp_nameT="A".$tA_t."Z".$tZ_t;
    
    if(! exists $table_AZ{$temp_nameT})
    {
	#print "target A$tA_t Z$tZ_t does not exist, next is";
	$tZ_t=$tZ_t+1;
	$tA_t=$table_Z_minA{$tZ_t};
	#print "A$tA_t Z$tZ_t\n";
	$temp_nameT="A".$tA_t."Z".$tZ_t;
    }

}


$nameout = "File_".$MaxTargetName."_".$MaxProjName.".dat";

open (OUT,"> $nameout") || die "impossible d'ecrire: $!";

for my $k1 ( sort keys %table_result ) 
{
    #print "fragment: $k1 \n";
    for my $k2 ( keys %{ $table_result{$k1} } ) 
    {
	for my $k3 (keys %{ $table_result{ $k1 }{ $k2 } })
	{
	    my $temp_res = $table_result{ $k1 }{ $k2 }{ $k3 };
	    
	    my ($tA_f,$tZ_f)=($k1=~/A([0-9]*)Z([0-9]*)$/);
	    my ($tA_p,$tZ_p)=($k2=~/A([0-9]*)Z([0-9]*)$/);
	    my ($tA_t,$tZ_t)=($k3=~/A([0-9]*)Z([0-9]*)$/);
	    
	    print OUT "$tA_f $tZ_f $tA_p $tZ_p $tA_t $tZ_t $temp_res\n";

	    #print "proj: $k2 target: $k3  = $temp_res\n";
	}
    }
    #print "\n";
}

