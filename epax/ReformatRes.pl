#! /usr/bin/perl

my $databaseFile = $ARGV[0];

my $FileOutDatabase = $databaseFile."UpETP.dat";

open (INDATABASE, $databaseFile) || die "impossible to open the database file : $! !\n";

open( OUTDATABASE, ">", $FileOutDatabase) || die "impossible to open the outfile for database : $! !\n";

while(my $lineData=<INDATABASE>)
{
    #                                                                                                                 Af       Zf      Ap       Zp       At       Zt          Rho          xs   Ratio  para1 para2   Z*Z    namfileRes  namefileParaRes
    (my $Af, $Zf, $Ap, $Zp, $At, $Zt, $Rho, $xs, $Ratio,$par1,$par2,$par3, $par4,$nameRes,$nameParaRes) = ($lineData =~ /^([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]*\.[0-9]*) (.+?) (.+?) ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+) (.*?\.dat) (.*?\.dat)$/ );
#    print "A$Af.Z$Zf A$Ap.Z$Zp A$At.Z$Zt Rho:$Rho / CX: $xs, $Ratio / $par1, $par2, $par3, $par4 / $nameRes + $nameParaRes\n"; 

    if($Af == $Ap and $Zf == $Zp)
    {
	next;
    }
    
    my $nameRes2 = "../mocadi/".$nameRes;
    my $nameParaRes2 = "../mocadi/".$nameParaRes;

    (my $nameResBase) = ( $nameRes2 =~ /(.*?)_result.dat/);
#    print $nameResBase."\n";
    
    
    my @fileOuts = `ls $nameResBase*out`;
    my %table_Tck_ES;
    foreach my $fileOut (@fileOuts)
    {
	chomp($fileOut);
	(my $Tck) = ($fileOut =~ /_Tck(.+?)_/);
	#print "$fileOut : Tck= $Tck ";

	open (INMOCADIRUN, $fileOut) || die "impossible to open Outfile of MOCADI : $! !\n";
	my $SurvivalRatio = 0;
	my $MeanEnergy = 0;
	my $SigmaEnergy = 0;
	my $First = 0;
	while (my $lineMocadi = <INMOCADIRUN>)
	{
	    if($lineMocadi =~ /survival ratio/)
	    {
		($SurvivalRatio) = ($lineMocadi =~ /survival ratio = ([0-9]+\.[0-9]+?) /);
		#print "Survival $SurvivalRatio\n";
	    }
	    if($First == 1)
	    {	       
		#print "in first : $lineMocadi \n";
		if((my $TempMeanEnergy, my $TempSigmaEnergy) = ($lineMocadi =~ /\<\s+energy\s+\>= (.*?) MeV\/u\s+sigma  energy\s+= (.*?) MeV/))
		{
		    $TempMeanEnergy =~ /(-?\d.\d+)e(.\d+)/;
		    my ($const, $expon) = ($1, $2);
		    $MeanEnergy = $const * 10 ** $expon;

		    $TempSigmaEnergy =~ /(-?\d.\d+)e(.\d+)/;
		    my ($const2, $expon2) = ($1, $2);
		    $SigmaEnergy = $const2 * 10 ** $expon2;		    
		    #$Frist = 0;
		}
	    }
	    
	    if($lineMocadi =~ /i_fragment.*?([0-9]*)$/)
	    {
		#print "First SELECTED\n";
		$First = $1;
	    }
	}
	close INMOCADIRUN;

	#print " Energy : $MeanEnergy +- $SigmaEnergy \n";
	$table_Tck_ES{$Tck} = [$SurvivalRatio, $MeanEnergy, $SigmaEnergy];
    }
    
    
    open (INRES, $nameRes2) || die "impossible to open Res file : $! \n";
    my %table_Tck_T ;
    while(my $lineRes = <INRES>)
    {
	(my $Tck, $Trans, $TRatio) = ($lineRes =~ /\s+([0-9]+\.)[0-9]*\s+.+?\s+.+?\s+.+?\s+.+?\s+.+?\s+.+?\s+.+?\s+.+?\s+.+?\s+.+?\s+.+?\s+.+?\s(.+?)\s(.+?)\s/);
	#print "Tck : $Tck -> $Trans | $TRatio \n";
	# my $TRatio2;
	# if($TRatio =~ /(-?\d.\d+)e(.\d+)/);
	# {
	#     my ($const, $expon) = ($1, $2);
	#     $TRatio2 = $const * 10 ** $expon;
	# }

	$table_Tck_T{$Tck} = [$Trans, $TRatio];	
    }
    close INRES;


    open (INPARARES, $nameParaRes2) || die "impossible to open Res file : $! \n";
    my %table_Tck_P ;
    #print "$nameParaRes2\n";
    while(my $lineParaRes = <INPARARES>)
    {
	if((my $Tck, $Para1, $Para2, $Para3) = ($lineParaRes =~ /([0-9]+\.*)[0-9]* .+? 1 .+? [0-9]\s(.+?)\s[0-9]\s(.+?)\s[0-9]\s(.+?)\s/))
	    {
		#print "Tck : $Tck --> $Para1 $Para2 $Para3 \n";
		my $Tck2 = "$Tck\.";
		$table_Tck_P{$Tck2} = [$Para1, $Para2, $Para3];	
	    }
    }
    close INPARARES;

    
    foreach my $key_Tck (keys(%table_Tck_ES))
    {
	print OUTDATABASE "$Af $Zf $Ap $Zp $At $Zt $Rho $xs $Ratio $par1 $par2 $par3 $par4 $key_Tck ";
	
	if($table_Tck_ES{$key_Tck}>0)
	{
	    #print "$key_Tck : Energy : $table_Tck_ES{$key_Tck}[0] $table_Tck_ES{$key_Tck}[1] $table_Tck_ES{$key_Tck}[2] |";
	    print OUTDATABASE "$table_Tck_ES{$key_Tck}[0] $table_Tck_ES{$key_Tck}[1] $table_Tck_ES{$key_Tck}[2] ";
	}
	if($table_Tck_T{$key_Tck}>0)
	{
	    #print " Trans : $table_Tck_T{$key_Tck}[0] $table_Tck_T{$key_Tck}[1] |"; 
	    print OUTDATABASE "$table_Tck_T{$key_Tck}[0] $table_Tck_T{$key_Tck}[1] ";
	}
	if($table_Tck_P{$key_Tck}>0)
	{
	    #print " Para : $table_Tck_P{$key_Tck}[0] $table_Tck_P{$key_Tck}[1] $table_Tck_P{$key_Tck}[2] |"; 
	    print OUTDATABASE "$table_Tck_P{$key_Tck}[0] $table_Tck_P{$key_Tck}[1] $table_Tck_P{$key_Tck}[2]";
	}
	#print "\n";
	print OUTDATABASE "\n";
    }
    
}

close OUTDATABASE;
close INDATABASE;
