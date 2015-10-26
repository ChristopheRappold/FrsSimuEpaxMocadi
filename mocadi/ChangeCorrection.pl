#! /usr/bin/perl

use File::Copy;

my @file_old = `ls ./header_files`;

foreach my $file (@file_old)
{
    if($file =~ /_result_para.dat$/)
    {
	
	chomp($file);
	my $old_file = "./header_files/";
	$old_file.=$file;
	print "$old_file ";
	my $tempfile = $old_file;
	$tempfile =~ s/\.[^.]+$//;
	$tempfile.="_old.dat";
	print "$tempfile \n";
	move($old_file,$tempfile);
	#system "mv $file $tempfile";
    }
}

foreach my $file (@file_old)
{
    if($file =~ /_result_para_new.dat$/)
    {
	
	chomp($file);
	my $old_file = "./header_files/";
	$old_file.=$file;
	print "$old_file ";
	my $tempfile = $old_file;
	$tempfile =~ s/\.[^.]+$//;
	$tempfile =~ s/_new//;
	$tempfile.=".dat";
	print "$tempfile \n";
	move($old_file,$tempfile);
	#system "mv $file $tempfile";
    }
}

