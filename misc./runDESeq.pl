use Statistics::R;

sub createLib{
	$libfile=shift;
	$outfile=shift;
	$line=shift;
	@index=split("\t",$line);
	%needlib=();
	foreach $i (@index){
		$needlib{$i}=1;
	}
	open(LIB,$libfile);
	open(OUT,">$outfile");
	
	$header=<LIB>;
	chomp($header);
	
	@ele=split("\t",$header);
	%hitIndex=();
	$newHeader=$ele[0];
	foreach $j (0..$#ele){
		if(exists $needlib{$ele[$j]}){
			$hitIndex{$j}=1;
			$newHeader.="\t".$ele[$j];
		}
	}
	print OUT $newHeader,"\n";
	 
	while(<LIB>){
		chomp;
		@col=split("\t");
		$line=$col[0];
		foreach $j (1..$#col){
			if(exists $hitIndex{$j}){
				$line.="\t".$col[$j];
			}
		}
		print OUT $line,"\n";
	}
	close(LIB);
	close(OUT);
}
=head
sub runDESeq{
	$infile=shift;
	$outfile=shift;
	$R=Statistics::R->new();
	$R->run(
		q`library(DESeq)`,
		qq`countsTable<-read.delim("./$infile", header=TRUE, stringsAsFactors=TRUE)`,
		q`rownames(countsTable)<-countsTable$subsystem`,
		q`countsTable<-countsTable[,-1]`,
		q`conds<-c("C","C","C","T","T","T")`,
		q`cds<-newCountDataSet(countsTable, conds)`,
		q`cds<-estimateSizeFactors(cds)`,
		q`cds<-estimateVarianceFunctions(cds)`,
		q`res<-nbinomTest(cds,"C","T")`,
		qq`write.table(res, file="./$outfile", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=TRUE, col.names=TRUE)`,
		#q`dev.off()`
		);
	$R->stop();
}
=cut

open(COM,$ARGV[0]);
$lib=$ARGV[1];
$outdir=$ARGV[2];

$i=0;
$R=Statistics::R->new();
while(<COM>){
	chomp;
	$i++;
	print $i,"\n";
	$line=$_;
	$tempfile="temp.tab";
	createLib($lib,$tempfile,$line);
	$outfile=$outdir.'/'.$i.'.res';
	$R->run(
		q`library(DESeq)`,
		qq`countsTable<-read.delim("./$tempfile", header=TRUE, stringsAsFactors=TRUE)`,
		q`rownames(countsTable)<-countsTable$subsystem`,
		q`countsTable<-countsTable[,-1]`,
		q`conds<-c("C","C","C","T","T","T")`,
		q`cds<-newCountDataSet(countsTable, conds)`,
		q`cds<-estimateSizeFactors(cds)`,
		q`cds<-estimateVarianceFunctions(cds)`,
		q`res<-nbinomTest(cds,"C","T")`,
		qq`write.table(res, file="./$outfile", append=FALSE, quote=FALSE, sep="\t", eol="\n", na="NA", dec=".", row.names=TRUE, col.names=TRUE)`,
		#q`dev.off()`
	);
	system("rm temp.tab");
	#last;
}
$R->stop();
close(COM);