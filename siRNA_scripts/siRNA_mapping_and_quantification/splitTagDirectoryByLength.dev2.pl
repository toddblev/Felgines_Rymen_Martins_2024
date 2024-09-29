#!/usr/bin/perl -w
#
#
if (@ARGV < 1) {
	print STDERR "\n\tsplitTagDirectoryByLength.pl <tag directory> [options]\n";
	print STDERR "\n\tOptions:\n";
	print STDERR "\t\t-minLength <#> (default: 20)\n";
	print STDERR "\t\t-maxLength <#> (default: 25)\n";
	print STDERR "\n";
	exit;
}


my $tagDir = $ARGV[0];
$tagDir =~ s/\/+$//;

my $minLength = 20;
my $maxLength = 25;
for (my $i=1;$i<@ARGV;$i++) {
	if ($ARGV[$i] eq '-minLength') {
		$minLength = $ARGV[++$i];
	} elsif ($ARGV[$i] eq '-maxLength') {
		$maxLength = $ARGV[++$i];
	} else {
		print STDERR "!Wrong option\n";
		exit;
	}
}

my $ogTotal = 0;
open IN, "$tagDir/tagInfo.txt";
my $c = 0;
while (<IN>) {
	$c++;
	next if ($c < 2);
	chomp;
	my @line = split /\t/;
	$ogTotal = $line[2];

	last;
}
close IN;

print STDERR "Total reads= $ogTotal\n";

my $logFileName = $tagDir . ".log";
open LOG, ">$logFileName";

for (my $i=$minLength;$i<=$maxLength;$i++) {
	my $newDir = "$tagDir-$i" . "nt";
	`cp -r "$tagDir" "$newDir"`;

	my @tagFiles = ();
	`ls "$newDir"/*.tags.tsv > .ls`;
	open IN, ".ls";
	while (<IN>) {
		chomp;
		push(@tagFiles, $_);
	}
	close IN;

	my $total = 0;
	my $n = 0;
	foreach(@tagFiles) {
		my $file = $_;
		open IN, $file;
		open OUT, ">.tmp";
		while (<IN>) {
			my $og = $_;
			chomp;
			my @line = split /\t/;
			my $v = $line[4];
			$total+=$v;
			if ($line[5] eq $i) {
				print OUT $og;
				$n+=$v;
			}
		}
		close IN;
		close OUT;
		`mv ".tmp" $file`;
	}
	my $fraction = sprintf("%.2lf",100*$n/$total);
	my $fraction2 = sprintf("%.2lf",100*$n/$ogTotal);
	print LOG "\tProcessing $tagDir  $i  ($minLength to $maxLength): reads $n of $total ($fraction%) [total mappable: $n or $ogTotal ($fraction2%)]\n";
	print STDERR "\tProcessing $tagDir  $i  ($minLength to $maxLength): reads $n of $total ($fraction%) [total mappable: $n or $ogTotal ($fraction2%)]\n";
}
