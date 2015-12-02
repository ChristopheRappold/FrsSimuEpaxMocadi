#! /usr/bin/perl

#use lib "/u/crappold/frs/epax/usr/local/share/perl";
use lib "./usr/lib/perl5";
use IO::CaptureOutput qw/capture/;
my ($stdout, $stderr);


open my $etcRoot, "<", "../RadioNuclides.txt" or die "error !";

%table_ListStableLifeTime=();

while(my $line=<$etcRoot>)
{
    if($line=~ /^\s+([0-9]+)-([A-Z][a-z]?)-([0-9]*)\s+\d+\s+\d+\s+\d+\s+[0-9]+\.?[0.9]*?\s+.*?\s+(.*?)\s+/)
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
	if($T2 > 10.)
	{
	    $table_ListStableLifeTime{$Name2} = $T2;
	}
    }
}

close $etcRoot;


%table_Z=();
$table_Z{"0"}="n";
$table_Z{"1"}="H";$table_Z{"2"}="He";$table_Z{"3"}="Li";$table_Z{"4"}="Be";$table_Z{"5"}="B";$table_Z{"6"}="C";$table_Z{"7"}="N";$table_Z{"8"}="O";$table_Z{"9"}="F";
$table_Z{"10"}="Ne";$table_Z{"11"}="Na";$table_Z{"12"}="Mg";$table_Z{"13"}="Al";$table_Z{"14"}="Si";$table_Z{"15"}="P";$table_Z{"16"}="S";$table_Z{"17"}="Cl";
$table_Z{"18"}="Ar";$table_Z{"19"}="K";$table_Z{"20"}="Ca";$table_Z{"21"}="Sc";$table_Z{"22"}="Ti";$table_Z{"23"}="V";$table_Z{"24"}="Cr";$table_Z{"25"}="Mn";
$table_Z{"26"}="Fe";$table_Z{"27"}="Co";$table_Z{"28"}="Ni";$table_Z{"29"}="Cu";$table_Z{"30"}="Zn";$table_Z{"31"}="Ga";$table_Z{"32"}="Ge";$table_Z{"33"}="As";
$table_Z{"34"}="Se";$table_Z{"35"}="Br";$table_Z{"36"}="Kr";$table_Z{"37"}="Rb";$table_Z{"38"}="Sr";$table_Z{"39"}="Y";$table_Z{"40"}="Zr";$table_Z{"41"}="Nb";
$table_Z{"42"}="Mo";$table_Z{"43"}="Tc";$table_Z{"44"}="Ru";$table_Z{"45"}="Rh";$table_Z{"46"}="Pd";$table_Z{"47"}="Ag";$table_Z{"48"}="Cd";$table_Z{"49"}="In";
$table_Z{"50"}="Sn";$table_Z{"51"}="Sb";$table_Z{"52"}="Te";$table_Z{"53"}="I";$table_Z{"54"}="Xe";$table_Z{"55"}="Cs";$table_Z{"56"}="Ba";$table_Z{"57"}="La";
$table_Z{"58"}="Ce";$table_Z{"59"}="Pr";$table_Z{"60"}="Nd";$table_Z{"61"}="Pm";$table_Z{"62"}="Sm";$table_Z{"63"}="Eu";$table_Z{"64"}="Gd";
$table_Z{"65"}="Tb";$table_Z{"66"}="Dy";$table_Z{"67"}="Ho";$table_Z{"68"}="Er";$table_Z{"69"}="Tm";$table_Z{"70"}="Yb";$table_Z{"71"}="Lu";$table_Z{"72"}="Hf";
$table_Z{"73"}="Ta";$table_Z{"74"}="W";$table_Z{"75"}="Re";$table_Z{"76"}="Os";$table_Z{"77"}="Ir";$table_Z{"78"}="Pt";$table_Z{"79"}="Au";$table_Z{"80"}="Hg";
$table_Z{"81"}="Tl";$table_Z{"82"}="Pb";$table_Z{"83"}="Bi";$table_Z{"84"}="Po";$table_Z{"85"}="At";$table_Z{"86"}="Rn";$table_Z{"87"}="Fr";$table_Z{"88"}="Ra";
$table_Z{"89"}="Ac";$table_Z{"90"}="Th";$table_Z{"91"}="Pa";$table_Z{"92"}="U";$table_Z{"93"}="Np";$table_Z{"94"}="Pu";$table_Z{"95"}="Am";$table_Z{"96"}="Cm";
$table_Z{"97"}="Bk";$table_Z{"98"}="Cf";$table_Z{"99"}="Es";$table_Z{"100"}="Fm";$table_Z{"101"}="Md";$table_Z{"102"}="No";$table_Z{"103"}="Lr";$table_Z{"104"}="Rf";
$table_Z{"105"}="Db";$table_Z{"106"}="Sg";$table_Z{"107"}="Bh";$table_Z{"108"}="Hs";$table_Z{"109"}="Mt";$table_Z{"110"}="Ds"; 



open my $fileheader, '<', "template_header.in" or die "error opening $filename: $!";
my $header = do { local $/ = undef; <$fileheader> };

close $fileheader;

@array_thickness;
push(@array_thickness,"1.");
push(@array_thickness,"2.");
push(@array_thickness,"5.");
push(@array_thickness,"8.");
push(@array_thickness,"10.");
push(@array_thickness,"13.");
push(@array_thickness,"15.");
push(@array_thickness,"20.");

print @array_thickness;
print "\n";

open my $fileEpax, '<', $ARGV[0] or die "error opening $fileEpax: $!";

(my $new_database_file = $ARGV[0]) =~ s/\.[^.]+$//;
$new_database_file.="_mocadi.dat";

open my $fileEpaxMocadiOut, '>', $new_database_file or die "error opening $new_database_file : $!";

$numberARG = $#ARGV;
print "number ARGV : $numberARG\n";
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

	my $AoverZ=$Af/$Zf;
	
	
	my $Aparasist1 = int($AoverZ*($Zf-1));
	$Aparasist1 -= 1;
	my $Zparasist1 = int($Zf-1);
	my $NameTemp1 = "A".$Aparasist1."Z".$Zparasist1;
	while(!exists $table_ListStableLifeTime{$NameTemp1})
	{
	    $Aparasist1 += 1;
	    $NameTemp1 = "A".$Aparasist1."Z".$Zparasist1;
	}

	my $Aparasist2 = int($AoverZ*($Zf-2));
	$Aparasist2 -= 1;
	my $Zparasist2 = int($Zf-2);
	my $NameTemp2 = "A".$Aparasist2."Z".$Zparasist2;
	while(!exists $table_ListStableLifeTime{$NameTemp2})
	{
	    $Aparasist2 += 1;
	    $NameTemp2 = "A".$Aparasist2."Z".$Zparasist2;
	}


	my $Aparasist3 = int($AoverZ*($Zf-3));
	$Aparasist3 -= 1;
	my $Zparasist3 = int($Zf-3);
	my $NameTemp3 = "A".$Aparasist3."Z".$Zparasist3;
	while(!exists $table_ListStableLifeTime{$NameTemp3})
	{
	    $Aparasist3 += 1;
	    $NameTemp3 = "A".$Aparasist3."Z".$Zparasist3;
	}
	   

	$Rho = 1000*$Rho; 

	my %table_result=();
	my %table_result_th=();
	
	my $name_file_result = "header_files/".$table_Z{$Zp}.$Ap."_".$table_Z{$Zt}.$At."_Frag_".$table_Z{$Zf}.$Af."_result.dat";
	my $name_file_result_para = "header_files/".$table_Z{$Zp}.$Ap."_".$table_Z{$Zt}.$At."_Frag_".$table_Z{$Zf}.$Af."_result_para.dat";

	foreach $Thick (@array_thickness)
	{
	    my $new_header = $header;

	    $new_header =~ s/TempBeamA/$Ap/g;
	    $new_header =~ s/TempBeamZ/$Zp/g;
	    $new_header =~ s/TempTargetA/$At/g;
	    $new_header =~ s/TempTargetZ/$Zt/g;
	    $new_header =~ s/TempFragA/$Af/g;
	    $new_header =~ s/TempFragZ/$Zf/g;
	    $new_header =~ s/TempTargetRho/$Rho/g;
	    $new_header =~ s/TempThCm/$Thick/g;

	    if($numberARG == "1")
	    {
	    	$new_header =~ s/\*FRAGMENT/FRAGMENT/g;
	    	$new_header =~ s/\*2, 4/2, 3/g;
	    	$new_header =~ s/\*  TempParaZ1/  $Zparasist1/g;
	    	$new_header =~ s/TempParaA1/$Aparasist1/g;
	    	$new_header =~ s/\*  TempParaZ2/  $Zparasist2/g;
	    	$new_header =~ s/TempParaA2/$Aparasist2/g;
	    	$new_header =~ s/\*  TempParaZ3/  $Zparasist3/g;
	    	$new_header =~ s/TempParaA3/$Aparasist3/g;
	    }

	    my $name_file = "header_files/".$table_Z{$Zp}.$Ap."_".$table_Z{$Zt}.$At."_Frag_".$table_Z{$Zf}.$Af."_Tck".$Thick."_header.in";
	    
	    print "process for fragment : A".$Af."Z".$Zf." with beam A".$Ap."Z".$Zp." and A".$At."Z".$Zt." (Rho : ".$Rho.") and Thickness ".$Thick.": ".$name_file."\n";

	    #$file_out = 
	    open OUT, ">", $name_file or die "error opening outfile  $!";
	    print OUT $new_header;
	    close OUT;
	    
	    capture sub 
	    {
		system "./LoopMocadi.pl $name_file";
	    } => \$stdout, \$stderr;
	    
	    #print $stdout;
	    
	    (my $file_out = $name_file) =~ s/\.[^.]+$//;
	    $file_out_mocadi = $file_out."_temp_6.out";
	    
	    open IN2, "<", $file_out_mocadi or die "error opening outfile  $!";
	    my $target_th;
	    my $target_th2;
	    my $step;
	    my $frag;
	    my $trans;
	    while (my $line_out =<IN2>)
	    {
		if($line_out=~ /thickness   =.*? (.*?)  cm/)
		{
		    $target_th = $1;
		    #print "Thickness ".$target_th."\n";
		}
		if(($string_in2) = ($line_out=~ /thickness   =\s*(.*?)mg/))
		{
		    $string_in2 =~ /(-?\d.\d+)e(.\d+)/;
		    my ($const2, $expon2) = ($1, $2);
		    $target_th2 = $const2 * 10 ** $expon2;
		    $table_result_th{$target_th}=$target_th2;
		    #print "Thickness2 ".$target_th2."\n";
		}
		if($line_out=~ /expected value.*?([0-9]*) \*\*/)
		{
		    $step = $1;
		}
		if($line_out=~ /i_fragment.*?([0-9]*)$/)
		{
		    $frag = $1;
		}
		
		if($line_out=~ /tr: opt    =\s*(.*?)\s*tr: total/)
		{
		    $trans = $1;
		    $table_result{$target_th}{$frag}{$step}=$trans;
		    #print $Thick." ".$target_th." ".$frag." ".$step." ".$trans."\n";
		}
	    }
	    close IN2;
	    print "Analysis of the Mocadi file $file_out_mocadi done !\n";
	}
	
	open OUT2, ">", $name_file_result or die "error opening outfile  $!";
	for my $k1 ( sort keys %table_result ) 
	{
	    print OUT2 $k1." ";
	    print OUT2 $table_result_th{$k1}." ";
	    for my $k2 ( sort keys %{ $table_result{$k1} } ) 
	    {
		if($k2=="1")
		{
		    #print OUT2 $k2." ";
		    my $size = scalar(keys %{ $table_result{$k1}{$k2} });
		    print OUT2 $size." "; 
		    for my $k3 (sort keys %{ $table_result{$k1}{$k2} } )
		    {
			#my $temp_res = (sort keys %{ $table_result{$k1}{$k2} })[-1];
			my $temp_res = $table_result{ $k1 }{ $k2 }{$k3};
			my $temp_res2 = $temp_res * $Ratio * $k1;
			print OUT2 $temp_res." ".$temp_res2." ";
		    }
		}
	    }
	    print OUT2 "\n";

	}
	close OUT2;
	print "Result Parasist file $name_file_result done !\n";

	if($numberARG == "1")
	{
	    open OUT3, ">", $name_file_result_para or die "error opening outfile  $!";
	    print OUT3 "#Length(cm) Tickness(g/cm2) id A".$Af."Z".$Zf." id A".$Aparasist1."Z".$Zparasist1." id A".$Aparasist2."Z".$Zparasist2." id A".$Aparasist3."Z".$Zparasist3." \n";
	    for my $k1 ( sort keys %table_result ) 
	    {
		print OUT3 $k1." ";
		print OUT3 $table_result_th{$k1}." ";
		for my $k2 ( sort keys %{ $table_result{$k1} } ) 
		{
		    if($k2!="0")
		    {
			print OUT3 $k2." ";
			my @temp_array_res;
			for my $k3 (sort keys %{ $table_result{$k1}{$k2} })
			{
			    my $temp_res =  $table_result{$k1}{$k2}{$k3};
			    my $temp_res2 = $temp_res * $Ratio * $k1;
			    push(@temp_array_res,$temp_res2);
			    #print "$k1 $k2 $k3 $temp_res2 \n";
			}
			my $temp_res3 = @temp_array_res[-1];
			print OUT3 $temp_res3." ";
		    }
		}
		print OUT3 "\n";
		
	    }
	    close OUT3;
	    print "Result Parasist file $name_file_result_para done !\n";
	}

	my $line_new = $line;
	chomp $line_new;
	print $fileEpaxMocadiOut $line_new." ".$name_file_result;
	if($numberARG == "1")
	{
	    print $fileEpaxMocadiOut " ".$name_file_result_para;
	}
	print $fileEpaxMocadiOut "\n";
    }
}

close $fileEpax;
close $fileEpaxMocadiOut;
