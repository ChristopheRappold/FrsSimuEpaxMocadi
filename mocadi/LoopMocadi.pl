#! /usr/bin/perl



# Author : Christophe Rappold, Dr < c.rappold [()] gsi.de >
# Copyright 2013 Christophe Rappold
# License : GPLv3


# The script works without any internal modifications. Everything are passed by command line.
# Example : ./LoopMocadi.pl header_file.in template_file.in mocadi_exe
# In this way the script can used for any version of mocadi that is available in the system. The input files can adjusted also externally.

# Optional : Please install the corresponding package of perl : 
# Optional : IO::CaptureOutput from CPAN or your linux distribution package manager (ex: debian package is libio-captureoutput-perl)

use strict;
use warnings;

#use lib "/u/crappold/frs/epax/usr/local/share/perl";
#use IO::CaptureOutput qw/capture/;
#use lib "./usr/local/share/perl";
use lib "./usr/lib/perl5";
my $noCapture = 0;
eval
{
# If IO::CaptureOutput is installed 
    require IO::CaptureOutput;
    IO::CaptureOutput->import(qw/capture_exec/);
};
if($@)
{
    $noCapture = 1;
    die " No IO::CaptureOutput \n";
    
}
else
{
    $noCapture = 0;
    print " IO::CaptureOutput Loaded \n";
}

my ($stdout, $stderr);

die "Usage: $0 HEADER_FILE [TEMPLATE_FILE] [MOCADI_EXE] [MEANENERGY]\n" if($#ARGV != 3 and $#ARGV != 2 and $#ARGV != 1 and $#ARGV != 0) ;

my $template_mocadi = "template_heb_3.in";
if($#ARGV == 1 or $#ARGV == 2 or $#ARGV ==3)
{
    $template_mocadi = $ARGV[1];
}

die "TEMPLATEFILE not readable !" if(!-e $template_mocadi);
die "HEADERFILE not readable !" if(!-e $ARGV[0]);
die "MOCADI_EXE not executable !" if ($#ARGV == 2 and !-X $ARGV[2]);

my $mocadi_exe = "./source3.5/mocadi-35";
if($#ARGV == 2 or $#ARGV == 3)
{
    $mocadi_exe = $ARGV[2];
}

my $MeanEnergy =0;
if($#ARGV == 3)
{
    $MeanEnergy = $ARGV[3];
}


(my $file_in = $ARGV[0]) =~ s/\.[^.]+$//;

## -> Load the header file in the variable $header for later concatenation

open my $fileheader, '<', $ARGV[0] or die "error opening $ARGV[0]: $! \n";
my $header = do { local $/ = undef; <$fileheader> };

close $fileheader;

## -> Parse the header file to set the fragment of interest $Af $Zf

open my $fileheader2, '<', $ARGV[0] or die "error opening $ARGV[0]: $! \n";
my $Af =0;
my $Zf =0;
while(my $lineheader=<$fileheader2>)
{
#    print $lineheader;

    if($lineheader =~ /^\*?\s*([0-9]*\.[0-9]*|[0-9]*)\,\s*([0-9]*)\,\s*Sollfragment/)
    {	
#	print $lineheader;
	$Af = $1;
	$Zf = $2;
    }
}
if($Af==0 or $Zf==0)
{
    die "No valid fragment A=$Af Z=$Zf !\n";
}
close $fileheader2;
print "fragment :A".$Af."Z".$Zf."\n";

## -> 1) Open the second file of the mocadi setup : Matrices of the magnets and SAVE statements at the wanted place
## -> 2) Look for the save statements and save the position for the loop of the SAVE points

open IN2, "<", $template_mocadi or die " error opening $template_mocadi: $! \n";

my @array_save;
while(my $line2=<IN2>)
{
    my $savestate ;
    if(($savestate) = ($line2=~ /^\s*SAVE save(.*?)\s?$/))
    {
	push(@array_save,$savestate); 
    }
}

close IN2;

print @array_save;
print "\n";


## Loop over the save statements position 
## 1) Open the file of the structure of the separator 
## 2) Prepare the temporary mocadi input file to be run  
## 2 a) Add the prepared header file 
## 2 b) Open the file of results of the previous mocadi run
## 2 c) Extract the mean energy of each fragments
## 2 d) Remove from the filesystem the temporary files of the previous run
## 2 e) Update the structure of the separator with the new information from the previous mocadi run
## 3) Run mocadi with the prepared temporary input file of the current iteration
## 4) Parse the current temporary input file and remove the END statement to save it for the next loop iteration
## This until the last SAVE statement. The last files will remain as the final updated input file for the mocadi run and the output files of this run.


my $previous_save = "#0";
my $new_header;
my $previous_out = "0";
foreach my $save (@array_save)
{
    print $save."\n";
    my $namesave=substr $save, 1;
    my $file_in_mocadi = $file_in."_temp_".$namesave.".in";

    open OUT, ">", $file_in_mocadi or die " error to open OUT $file_in_mocadi ! \n";
    print "outfile : $file_in_mocadi \n";

    open IN2, "<", $template_mocadi or die " error opening $template_mocadi : $! \n";

    my %table_fragment=();

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
	(my $file_out_root = $previous_out) =~ s/\.in/\.root/;
	
	
	print "process outfile :".$file_out."\n";
	if($MeanEnergy>0 && $previous_save =~ "#1")
	{
	    $table_fragment{0}{1}=$MeanEnergy;
	}
	else
	{
	    print "in file ".$file_out."\n";
	    open IN,"<", $file_out or die "error opening $file_out: $! \n";
	    
	    my $frag = 0;
	    my $Energy = 0;
	    my $step = 0;
	    
	    while (my $line =<IN>)
	    {
		if($line=~ /expected value.*?([0-9]*) \*\*/)
		{
		    $step = $1;
		}
		if($line=~ /i_fragment.*?([0-9]*)$/)
		{
		    $frag = $1;
		}
		
		if(my ($string_in) = ($line=~ / \<  energy   \>= (.*?) MeV/))
		{
		    $string_in =~ /(-?\d.\d+)e(.\d+)/;
		    my ($const, $expon) = ($1, $2);
		    my $Energy2 = $const * 10 ** $expon;
		    #print "step #".$step." fragment ".$frag." = ".$Energy."\n";
		    #print $Energy2."\n";
		    $table_fragment{$frag}{$step}=$Energy2;
		}
	    }
	    close IN;
	    
	    ## this is to remove the temporary file of previous iterration. You can comment them if you want to keep them all
	    #
	    unlink $file_out;
	    unlink $previous_out;
	    unlink $file_out_hbk2;
	    unlink $file_out_root;
	    #
	    
	    for my $k1 ( sort keys %table_fragment ) 
	    {
		for my $k2 ( sort keys %{ $table_fragment{$k1} } ) 
		{
		    my $temp_res = $table_fragment{ $k1 }{ $k2 };
		    print "step #".$k2." fragment ".$k1." = ".$temp_res."\n";
		}
	    }
	}
    }

   # print "---- Debug:\n";
    # for my $k1 ( sort keys %table_fragment ) 
    # {
    # 	for my $k2 ( sort keys %{ $table_fragment{$k1} } ) 
    # 	{
    # 	    my $temp_res = $table_fragment{ $k1 }{ $k2 };
    # 	    print "step #".$k2." fragment ".$k1." = ".$temp_res."\n";
    # 	}
    # }
    # print " ----\n ";


    print "previous save :".$previous_save."\n";
    my $to_cat = 0;
    my $Energy;
    my $headersave = substr $previous_save, 1;
    if(exists( $table_fragment{"1"} ) )
    {
	$Energy = $table_fragment{ "1" }{  $headersave }; #"1" normally
    }
    else
    {
	$Energy = $table_fragment{ "0" }{  $headersave };
    }
    print "set E:".$Energy."\n" if(defined $Energy);
    while(my $line2=<IN2>)
    {
	$line2 =~ s/TempA/$Af/;
	$line2 =~ s/TempZ/$Zf/;

	if($line2=~/^\s*SAVE save$previous_save\s?/)
	{
	    $to_cat =1;
	    next;
	}
	if($to_cat==1)
	{
	    $line2 =~ s/\s*$Af.*?,.*?$Zf.*?,.*?TempEnergy.*?/ $Af , $Zf , $Energy/;
	    print OUT $line2;
	}
	if($line2=~ /^\s*SAVE save$save\s?/)
	{
	    last; 
	}
    }
    print OUT " END";
    close OUT;
    close IN2;

    $previous_save = $save;

    # If IO::CaptureOuput is installed or not
    if($noCapture == 0)
    {
	print "Exe: $mocadi_exe,$file_in_mocadi\n";
	$stdout = capture_exec($mocadi_exe,$file_in_mocadi);
    }
    else
    {
	#print "Backup Solution !\n";
	system "$mocadi_exe $file_in_mocadi >| log_LoopMocadi.log";
    }


    $previous_out = $file_in_mocadi;

    
    open IN3, "<", $file_in_mocadi or die " error opening $file_in_mocadi : $! \n";
    my $temp_new_header;
    while(my $linenewheader=<IN3>)
    {
	if($linenewheader !~ /^ END$/)
	{
	    $temp_new_header .= $linenewheader;
	}
    }
    $new_header = $temp_new_header;
    close IN3;
#print $new_header;
    
}
