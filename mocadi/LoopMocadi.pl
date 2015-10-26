#! /usr/bin/perl

use lib "/u/crappold/frs/epax/usr/local/share/perl";
use IO::CaptureOutput qw/capture/;
my ($stdout, $stderr);


(my $file_in = $ARGV[0]) =~ s/\.[^.]+$//;

open my $fileheader, '<', $ARGV[0] or die "error opening $filename: $!";
my $header = do { local $/ = undef; <$fileheader> };

close $fileheader;

open my $fileheader, '<', $ARGV[0] or die "error opening $filename: $!";
my $Af =0;
my $Zf =0;
while(my $lineheader=<$fileheader>)
{
    if($lineheader =~ /^.*?([0-9]*)\,.*?([0-9]*)\,.*?Sollfragment$/)
    {
	$Af = $1;
	$Zf = $2;
    }
}
close $fileheader;
print "fragment :A".$Af."Z".$Zf."\n";


# capture sub {
#     system "./source3.5/mocadi-35 $file_in"; 
# } => \$stdout, \$stderr;


#open IN, $ARGV[0] or die "impossible d'ouvrir !";

open IN2, "<", "template_heb.in" or die " impossible d'ouvrir !";

#open OUT, ">test.in" or die " impossible d'ouvrir !";

@array_save;
while(my $line2=<IN2>)
{
    my $savestate ;
    if(($savestate) = ($line2=~ /^ SAVE save(.*?)$/))
    {
	push(@array_save,$savestate); 
    }
}

close IN2;

print @array_save;
print "\n";

my $previous_save = "#0";
my $new_header;
my $previous_out = "0";
foreach $save (@array_save)
{
    print $save."\n";
    $namesave=substr $save, 1;
    $file_in_mocadi = $file_in."_temp_".$namesave.".in";

    open OUT, ">", $file_in_mocadi or die " impossible d'ouvrir out";
    print "outfile : $file_in_mocadi \n";

    open IN2, "<", "template_heb.in" or die " impossible d'ouvrir !";

    %table_fragment=();

    if($previous_save =~ "#0")
    {
	print "first header\n";
	print OUT $header;
    }
    else
    {
	print "update header\n";
	print OUT $new_header;

	(my $file_out = $previous_out) =~ s/\.in/\.out/;
	(my $file_out_hbk = $previous_out) =~ s/\.in/\.hbk/;
	my $file_out_hbk2 = lc($file_out_hbk);
	
	print "process outfile :".$file_out."\n";

	#my $file_in = "C12B10_2GeV_heb.out";
	print "in file ".$file_out."\n";
	open IN,"<", $file_out or die "impossible d'ouvrir !";
	
	my $frag = 0;
	my $Energy = 0;
	my $step = 0;
	
	while ($line =<IN>)
	{
	    if($line=~ /expected value.*?([0-9]*) \*\*/)
	    {
		$step = $1;
	    }
	    if($line=~ /i_fragment.*?([0-9]*)$/)
	    {
		$frag = $1;
	    }
	    
	    if(($string_in) = ($line=~ / \<  energy   \>= (.*?) MeV/))
	    {
		$string_in =~ /(-?\d.\d+)e(.\d+)/;
		my ($const, $expon) = ($1, $2);
		$Energy2 = $const * 10 ** $expon;
		#print "step #".$step." fragment ".$frag." = ".$Energy."\n";
		#print $Energy2."\n";
		$table_fragment{$frag}{$step}=$Energy2;
	    }
	}
	close IN;

	unlink $file_out;
	unlink $previous_out;
#	unlink $file_out_hbk2;

	for my $k1 ( sort keys %table_fragment ) 
	{
		for my $k2 ( sort keys %{ $table_fragment{$k1} } ) 
		{
		    my $temp_res = $table_fragment{ $k1 }{ $k2 };
		    print "step #".$k2." fragment ".$k1." = ".$temp_res."\n";
		}
	}
    }

    # for my $k1 ( sort keys %table_fragment ) 
    # {
    # 	for my $k2 ( sort keys %{ $table_fragment{$k1} } ) 
    # 	{
    # 	    my $temp_res = $table_fragment{ $k1 }{ $k2 };
    # 	    print "step #".$k2." fragment ".$k1." = ".$temp_res."\n";
    # 	}
    # }

    print "previous save :".$previous_save."\n";
    my $to_cat = 0;
    my $Energy;
    my $headersave = substr $previous_save, 1;
    $Energy = $table_fragment{ "1" }{  $headersave };
    print $Energy."\n";
    while(my $line2=<IN2>)
    {
	$line2 =~ s/ TempA / $Af /;
	$line2 =~ s/ TempZ / $Zf /;
	if($line2=~/^ SAVE save$previous_save$/)
	{
	    $to_cat =1;
	    next;
	}
	if($to_cat==1)
	{
	    # if($line2 =~ /^ SAVE save#([0-9]*)$/)
	    # {
	    # 	my $headersave = $1;
	    # 	$Energy = $table_fragment{ "1" }{  $headersave };
	    # 	print $Energy."\n";
	    # }
	    $line2 =~ s/ $Af  , $Zf , 2000./ $Af  , $Zf , $Energy/;
	    print OUT $line2;
	}
	if($line2=~ /^ SAVE save$save$/)
	{
	    last; 
	}
    }
    print OUT " END";
    close OUT;
    close IN2;

    $previous_save = $save;

    capture sub 
    {
     	system "./source3.5/mocadi-35 $file_in_mocadi";
    } => \$stdout, \$stderr;

    $previous_out = $file_in_mocadi;

    
    open IN3, "<", $file_in_mocadi or die " impossible d'ouvrir !";
    my $temp_new_header;
    my $Energy;
    while($linenewheader=<IN3>)
    {
	if($linenewheader !~ /^ END$/)
	{
	    $temp_new_header .= $linenewheader;
	}
    }
    $new_header = $temp_new_header;
#print $new_header;
    
}


