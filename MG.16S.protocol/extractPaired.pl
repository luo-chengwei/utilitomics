## this script takes blastn results from read-vs-reference 16S (I usually use E.coli 16s)
## and find the reads with both ends mapped to 16S.

## usage: perl extractPaired.pl [infile in blast m8 format] 
## the output file would be the .PairedEnd.bls file.

open(IN,$ARGV[0]);
$outfile=$ARGV[0];
$outfile=~s/\.bls$/PairedEnd.bls/;
open(OUT,">$outfile");

@tags=();
%lib=();
while(<IN>){
  chomp;
  $line=$_;
  @col=split("\t",$line);
  $col[0]=~s/\/\d{1}$//;
  push(@tags,$col[0]);
  if(!exists $lib{$col[0]}){
    $lib{$col[0]}=();
  }
  push(@{$lib{$col[0]}},$line);
}

foreach $tag (@tags){
  next if $#{$lib{$tag}}!=1;
  print OUT $lib{$tag}[0],"\n",$lib{$tag}[1],"\n";
  delete $lib{$tag};
}