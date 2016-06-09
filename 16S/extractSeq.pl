open(LIB,$ARGV[0]);
open(SEQ,$ARGV[1]);
open(OUT,">$ARGV[2]");

%lib=();
while(<LIB>){
  chomp;
  $line=$_;
  @col=split("\t",$_);
  $lib{$col[0]}=1;
}
close(LIB);

while(<SEQ>){
  chomp;
  $tag=$_;
  $tag=~s/^\>//;
  $seq=<SEQ>;
  chomp($seq);
  if(exists $lib{$tag}){
    print OUT '>',$tag,"\n",$seq,"\n";
  }
}
close(SEQ);
close(OUT);
