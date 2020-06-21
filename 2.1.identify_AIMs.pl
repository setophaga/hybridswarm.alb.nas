use strict ; 
use warnings ; 

### input file p1.p2. is the first command line argument
my $input = $ARGV[0] ; 

### give the mpileup as the second command line argument
my $mpileup = $ARGV[1] ; 

### minimum distance between SNPs to include them in the analysis. 
my $min_distance = 100 ; 

### recombination rate per site per generations in Morgans per site
my $rec = 3e-8 ;  ### this value corresponds to 3cm/mb. 

### data object to store ancestry informative sites
my %aim ; 

print STDERR "Reading Reference Panel Input\n" ; 
my $total_aims = 0 ; 

open IN, "<$input" ;
<IN> ; 
while (<IN>) { 
	chomp ; 
	$_ =~ s/\"//g ;
	my @split = split ( /\,/, $_ ) ; 

	### require a minimum of 6 chromosomes and a fixation 
	if ( $split[2] > 3 && $split[7] > 3 && abs( $split[4] - $split[9] ) == 1 ) { 
		if ( $split[4] == 1 ) { 
			$aim{$split[0]}{$split[1]} = $split[3]."\t".$split[5]."\t".$split[2]*$split[4]."\t".$split[2]*(1-$split[4])."\t".$split[7]*$split[9]."\t".$split[7]*(1-$split[9]) ;
		}
		else {
			$aim{$split[0]}{$split[1]} = $split[5]."\t".$split[3]."\t".$split[2]*(1-$split[4])."\t".$split[2]*$split[4]."\t".$split[7]*(1-$split[9])."\t".$split[7]*$split[9] ;
		}
		$total_aims ++ ; 
	}
}
close IN ; 
print STDERR "Found ", $total_aims, " AIMs\n" ; 


### read from mpileup file, do not give samtools a reference when you make this
print STDERR "Populating SNP matrix and outputting\n" ;
my $last = -1 * $min_distance ; 
my $chrom = "NA" ; 
open IN, "<$mpileup" ; 
while (<IN>) { 

	chomp $_ ; 
	my @split = split ( /\t/, $_ ) ; 

	### update to the next chromosome
	if ( $split[0] ne $chrom ) { 
		$chrom = $split[0] ; 
		$last = -1 * $min_distance ; 
	}

	### check if the site exists and is sufficiently far from the last site
	if ( !exists( $aim{$split[0]}{$split[1]} ) || $split[1] - $last < $min_distance ) { 
		next ;
	}

	### get the site information 
	my @a = split ( /\t/, $aim{$split[0]}{$split[1]} ) ;

	### this uses a uniform recombination rate, 1e-8 per site, or 1cm/mb 
	print $split[0], "\t", $split[1], "\t", $a[2], "\t", $a[3], "\t", $a[4], "\t", $a[5], "\t", ( $split[1] - $last ) * $rec ; 
	$last = $split[1] ;

	### iterate through the mpileup and add sites to our data input file 
	for ( my $i = 4 ; $i < $#split + 1 ; $i += 3 ) { 
		my $c1 = 0 ; my $c2 = 0 ; 
		for ( my $p = 0 ; $p < length( $split[$i] ) ; $p ++ ) { 
			if ( substr( $split[$i], $p, 1 ) =~ m/$a[0]/i ) { 
				$c1 ++ ; 
			}
			elsif ( substr( $split[$i], $p, 1 ) =~ m/$a[1]/i ) { 
        	                $c2 ++ ; 
			}
		}
		print "\t", $c1, "\t", $c2 ;		
	}
	print "\n" ;
}

