#!/usr/bin/perl -I/home/mzheng/mzheng-data/perl_prog/hus6 -I/home/mzheng/mzheng-data/perl_prog/common

use strict;
use FileHandle;
use commonFunctions;
use DBI;

#############################################
# ENVIRONMENT VARIABLES#
#############################################
$ENV {'ORACLE_BASE'} = '/export/app/oracle';
$ENV{'ORACLE_SID'} = 'lseq';
$ENV{'ORACLE_TERM'} = 'vt100';
$ENV {'ORACLE_HOME'} = '/export/app/oracle/product/11.1.0';


my $flanking = 10000;

my $tableSNPvsMMgene = "B36_PERLEGEN_SNP_VS_MMGENE";
my $tableSNPcodons = "B36_PERLEGEN_SNP_CODONS";


#### input SNP_vs_mmgene
my @row_array3;
my $dbh3 = DBI->connect('DBI:Oracle:PELTZ', 'peltz_prod','test')
                or die "Couldn't connect to database: " . DBI->errstr;

my $sql3 = "select * from $tableSNPvsMMgene where SNP_id = \'NES08846476\'";
my $sth3 = $dbh3->prepare($sql3);
$sth3->execute;

my ($i, $j,  $key);

my $count = 0;

while( @row_array3 = $sth3->fetchrow_array)
{
    for($i = 0; $i <= $#row_array3; $i++)
    {
	print ">$row_array3[$i]<\t";
    }
    print "\n";

}
$sth3->finish;
$dbh3->disconnect;


exit(0);

my $outputFile = "/home/mzheng/mzheng-data/projects/hapmapper/SNP_info/perlegen/snp_gene.txt";

#### input SNP_vs_mmgene
my @row_array2;
my $dbh2 = DBI->connect('DBI:Oracle:PELTZ', 'peltz_prod','test')
                or die "Couldn't connect to database: " . DBI->errstr;

my $sql2 = "select * from $tableSNPcodons";
my $sth2 = $dbh2->prepare($sql2);
$sth2->execute;

my ($i, $j,  $key);

my $count = 0;
my %hasCodonChange = ();

while( @row_array2 = $sth2->fetchrow_array)
{
    print "$row_array2[0]\n" if $row_array2[0] eq "NES08846476";
    $hasCodonChange{$row_array2[0]} = 1;
    $count++;
    if($count % 100000 == 0)
    {
	print "$count\n";
    }
}
$sth2->finish;
$dbh2->disconnect;

my $output = new FileHandle;
$output->open("> $outputFile") or die "$outputFile can't be opened\n";


#### input SNP_vs_mmgene
my @row_array1;
my $dbh1 = DBI->connect('DBI:Oracle:PELTZ', 'peltz_prod','test')
                or die "Couldn't connect to database: " . DBI->errstr;

my $sql1 = "select * from $tableSNPvsMMgene";
my $sth1 = $dbh1->prepare($sql1);
$sth1->execute;

my ($i, $j,  $key);

my $count = 0;

while( @row_array1 = $sth1->fetchrow_array)
{
    next if(!$row_array1[3]);
    
    if($row_array1[7] < $flanking)
    {
	if($hasCodonChange{$row_array1[1]})
	{
	    print $output "$row_array1[1],$row_array1[3],1\n";
	}
	else
	{
	    print $output "$row_array1[1],$row_array1[3],0\n";
	}
    }

    $count++;
    if($count % 100000 == 0)
    {
	print "$count\n";
    }
}
$sth1->finish;
$dbh1->disconnect;

close($output);
print "succeeded.\n";


