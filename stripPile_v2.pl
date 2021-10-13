#!/usr/bin/perl

$noIndel = $ARGV[0];

use Getopt::Long;

GetOptions (
#	'noIndel' => \$noIndel,	 #filter out sites with any indel component
	'cutIndel' => \$cutIndel,	#just remove the indel-containing call, but allow the site
	'callHap' => \$callHap,  #call pseudohaplotypes
	'makeCns' => \$makeCns,	 #dump out a fasta file, wrapped and according to the pileup ref names
	'mincvg=s' => \$mincvg,
	'maxcvg=s' => \$maxcvg,
	'minsup=s' => \$minsup,	#minimum support to consider an allele
	'gzip' => \$gzip,	#gzip the output
	'out=s' => \$out,	#out stem
	'seed=s' => \$seed,	#srand seed
);

if ($seed) {
	srand($seed);
}

$minsup = 1 unless ($minsup);
$mincvg = 1 unless ($mincvg);
$maxcvg = 1e9 unless ($maxcvg);

if ($makeCns) {
	$callHap = 1;
	#$ct = 0;
	if ($gzip) {
		open OUT, "| gzip -c > $out.fa.gz";
	}
	else {
		open OUT, ">$out.fa";
	}
}
elsif ($callHap) {
	if ($gzip) {
		open OUT, "| gzip -c > $out.hap.txt.gz";
	}
	else {
		open OUT, ">$out.hap.txt";
	}
	print OUT "$out\n";
}

while ($line = <STDIN>) {
	chomp $line;
	@d = split /\s+/, $line;
	$ping = $d[4];
	$ref = $d[2];
	$ping =~ s/[^acgtACGTnN\-\+]//g;
	if ($noIndel) {
		$checkIndel = 0;
		if ($ping =~ /[\-\+]/) {
			$checkIndel = 1;
		}
	}
	if ($cutIndel and $ping =~ /[\-\+]/) {
		$checkIndel = 0;
		@id = split(/(?=[\-\+])/, $ping);
		if (@id>0) {
			@new = ();
			foreach $seg (@id) {
				@s = split(/(?<=[0-9])/, $seg);
				$correct = $s[0];
				$safe = $s[1];
				$correct =~ s/[\-\+]//;
				$s[1] =~ s/^.{$correct}//;
				push @new, $s[1];
				$d[3]--;
			}
		}
		$ping = join "", @new;
#		print $ping."\n";
	}
	$ping =~ s/[\.\,]/$ref/g;
	$ping =~ s/[^acgtACGT]//g;
	$ping =~ tr/[acgt]/[ACGT]/;
	$len = length $ping;
	$ping =~ s/A//g;
	$a = $len - (length $ping);
	$len = length $ping;
	$ping =~ s/C//g;
	$c = $len - (length $ping);
	$len = length $ping;
	$ping =~ s/G//g;
	$g = $len - (length $ping);
	$t = length $ping;
	$tot = $a + $c + $g + $t;
	pop @d;
	pop @d;
	pop @d;
	pop @d;
	$set = join "\t", @d;
	$chr = $d[0];
	if ($checkIndel==1) {
		$tot = 0;
		$a = 0;
		$c = 0;
		$g = 0;
		$t = 0;
	}
	if ($callHap) {
		if ($tot >= $mincvg and $tot <= $maxcvg) {
			undef @o;
			if ($a >= $minsup) {
				push @o, "A";
			}
			if ($c >= $minsup) {
				push @o, "C";
			}
			if ($g >= $minsup) {
				push @o, "G";
			}
			if ($t >= $minsup) {
				push @o, "T";
			}
			if (@o == 0) {
				$pick = 0;
			}
			else {
				$pick = int rand @o;
				$pick = $o[$pick];
			}
		}
		else {
			$pick = 0;
		}
	}
	if ($makeCns) {
		if ($chr ne $lastchr) {
			if ($. == 1) {
				print OUT ">$chr\n";
				$ct = -1;
			}
			else {
				print OUT "\n>$chr\n";
				$ct = -1;
			}
		}
		$lastchr = $chr;
		$ct++;
		if ($ct==100) {
			print OUT "\n";
			$ct = 0;
		}
		else {
		}
		$pick =~ s/0/N/;
		print OUT $pick;
	}
	elsif ($callHap) {
		print OUT $pick."\n";
	}
	else {
		print "$set\t$tot\t$a\t$c\t$g\t$t\n";
	}
}

=pod
#!/usr/bin/perl
#
#Takes a bed file of SNP sites (chr, beg, fin) and a pileup file run through stripPile.pl, and adds any missing sites from the bed to make a complete pileup
#

($pile, $bed, $out) = @ARGV;

die "\nUsage: perl $0 infile.pile markers.bed outfile\n\n" unless (@ARGV==3);

open IN, $pile;
while ($line = <IN>) {
	chomp $line;
	($chr, $pos, $n, $a, $c, $g, $t) = split /\s+/, $line;
	$ping = "$chr.$pos";
	$ct = $a + $c + $g + $t;
	$line = join "\t", ($chr, $pos, $ct, $a, $c, $g, $t);
	$hold{$ping} = $line;
}
close IN;

open BED, $bed;
open OUT, ">$out";
while ($line = <BED>) {
	chomp $line;
	($chr, $n, $pos) = split /\s+/, $line;
	$ping = "$chr.$pos";
	$hold{$ping} = "$chr\t$pos\t0\t.\t.\t.\t." unless (exists ($hold{$ping}));
	print OUT $hold{$ping}."\n";
}
