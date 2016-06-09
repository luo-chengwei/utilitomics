$infile=$ARGV[0];
$outprefix=$ARGV[1];

$tmpFaa=$outprefix.'.faa.tmp';
$tmpFfn=$outprefix.'.ffn.tmp';
$faa=$outprefix.'.faa';
$ffn=$outprefix.'.ffn';
$gbk=$outprefix.'.gbk';

system("prodigal -a $tmpFaa -d $tmpFfn -f gbk -i $infile -o $gbk -p meta");

open(FATMP,$tmpFaa);
open(FNTMP,$tmpFfn);
open(FAA,">$faa");
open(FFN,">$ffn");

while(<FATMP>){
  chomp;
  if(/^\>/){
    $_=~s/\s\#\s(\d+)\s\#\s(\d+)\s\#\s(\-?\d+)\s\#\s/\|$1\-$2\|$3\|/;
  }
  print FAA $_,"\n";
}
close(FATMP);
close(FAA);

while(<FNTMP>){
  chomp;
  if(/^\>/){
    $_=~s/\s\#\s(\d+)\s\#\s(\d+)\s\#\s(\-?\d+)\s\#\s/\|$1\-$2\|$3\|/;
  }
  print FFN $_,"\n";
}
close(FNTMP);
close(FNN);

system("rm $tmpFaa");
system("rm $tmpFfn");