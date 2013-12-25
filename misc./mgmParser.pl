#usage perl mgmParser.pl [infile] [output prefix]

open(IN,$ARGV[0]);
$faa=$ARGV[1].'.faa';
$ffn=$ARGV[1].'.ffn';
open(FAA,">$faa");
open(FFN,">$ffn");

@genes=();
@dna=();
@prot=();

while(<IN>){
  chomp;
  $line=$_;
  if(not $line=~/^\#/ and $line){
    @col=split("\t",$line);
    if($col[6] eq '+'){
      $gene=$col[0].'|'.$col[3].'-'.$col[4];
    }else{
      $gene=$col[0].'|c'.$col[4].'-'.$col[3];
    }
    push(@genes,$gene);
    
    while(1){
    	$line=<IN>;
    	chomp($line);
    	if(not $line){
    	  last;
    	}
    	if(not $line=~/^\#/){
    	  @col=split("\t",$line);
          if($col[6] eq '+'){
      	  $gene=$col[0].'|'.$col[3].'-'.$col[4];
    	  }else{
      		$gene=$col[0].'|c'.$col[4].'-'.$col[3];
    	  }
    	  push(@genes,$gene);
    	}else{
    	  if($line=~/^\#\#Protein/){
    	  	$prot_seq='';
    	  	while(1){
    	  		$seq=<IN>;
    	  		chomp($seq);
    	  		if($seq=~/^\#\#end-Protein/){
    	  			last;
    	  		}
    	  		$seq=~s/#//g;
    	  		$prot_seq.=$seq;
    	  	}
    	  	push(@prot,$prot_seq);
    	  }elsif($line=~/^\#\#DNA/){
    	    $dna_seq='';
    	  	while(1){
    	  		$seq=<IN>;
    	  		chomp($seq);
    	  		if($seq=~/^\#\#end-DNA/){
    	  			last;
    	  		}
    	  		$seq=~s/#//g;
    	  		$dna_seq.=$seq;
    	  	}
    	  	push(@dna,$dna_seq);
    	  }
    	}
    	
    }
  }
}

foreach $gene (@genes){
  $prot_seq=shift(@prot);
  $dna_seq=shift(@dna);
  print FAA ">$gene\n$prot_seq\n";
  print FFN ">$gene\n$dna_seq\n";
}

close(FAA);
close(FFN);
