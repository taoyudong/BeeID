## This script is designed to take as input .vcf and output diplotized .vcf # ## file
#usr/bin/perl -w
$infile1 = $ARGV[0]; # .vcf file 
$outfile = $ARGV[1]; # output file name

# Read .vcf.gz file 
#open(READ, "gunzip -c $infile1 |") || die "cannot open pipe to $infile1";
open(READ, "<$infile1") or die "cannot open the file $infile1\n";
open(WR,">$outfile");

my @temp=();
my @temp2=();
my $id=();
my $line = ();
my %coords_hash=();
my %vals_hash=();
my @sampleids = ();
my @snpids = ();
              
while($line=<READ>) # read line by line .vcf.gz file 
{
	chomp($line);
	if($line =~ /^\#CHROM/) # line that contains all sample ids
	{
		print WR "$line\n";
		#@temp = split(/\t/,$line);
		#@sampleids = @temp[9..@temp];     #print "@sampleids\n";		
	}
	elsif($line !~ /^\#/) # SNP lines
	{
		#print  "INside $line\n";	
		@temp = split(/\t/,$line);
		$chr = $temp[0];
		$crd = $temp[1];
		for($i=9; $i <@temp; $i++)
		{
			# capture the (genotype)value just before the semicolon
			($gtype) = ($temp[$i] =~ /(.*?)\:.*/);
			# replace the existing value with the phased diploid genotype
			$temp[$i] = $gtype."|".$gtype;
		}
		$line2 = join("\t", @temp);
		print WR "$line2\n";
	}	
	else # print the header
	{
		print  WR "$line\n";
	}
}

