open(LIB,$ARGV[0]);
open(BLS,$ARGV[1]);
open(OUT,">$ARGV[2]");

%lib=();
while(<LIB>){
  chomp;
  @col=split(/\s+/,$_);
  $id=$col[0];
  $desc=$col[1];
  $lib{$id}=$desc;
}
close(LIB);

while(<BLS>){
  chomp;
  @col=split(/\s+/,$_);
  $id=$col[1];
  $desc=$lib{$id};
  push(@col,$desc);
  $line=join("\t",@col);
  print OUT $line,"\n";
}