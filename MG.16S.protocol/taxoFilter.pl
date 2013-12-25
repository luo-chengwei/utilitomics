open(IN,$ARGV[0]);
open(OUT,">$ARGV[1]");
while(<IN>){
    chomp;
    $line1=$_;
    @col1=split("\t",$_);
    @ele1=split(';',$col1[-1]);
    $line2=<IN>;
    chomp($line2);
    @col2=split("\t",$line2);
    @ele2=split(';',$col2[-1]);
    $order1='NA';
    $order2='NA';
    foreach $e (@ele1){
      if($e=~/^o__/){
        $order1=$e;
        $order1=~s/^o\_\_//;
      }
    }
    foreach $e (@ele2){
      if($e=~/^o__/){
        $order2=$e;
        $order2=~s/^o\_\_//;
      }
    }
    #print $order1,"\t",$order2,"\n";
    next if $col1[3]<40 or $col2[3]<40;
    next if $col1[2]<90 or $col2[2]<90;
    next if $order1 ne $order2;
    print OUT $line1,"\n",$line2,"\n";
}
