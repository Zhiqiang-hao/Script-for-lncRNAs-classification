#!/usr/bin/perl
open$in1,"/home/li/Other_studys/lncRNA_VVI/Genome/annotation/Vvinifera_457_Genoscope.12X.fa"or die;
while($in2=<$in1>){
chomp$in2;
if($in2=~/>(\S+)/){
$chr=$1;next;
}
$seq{$chr}+=length$in2;
}
open$in1,"/home/li/Other_studys/lncRNA_VVI/Genome/annotation/Vvinifera_457_v2.1.gene_exons_only.gff3"or die;
while($in2=<$in1>){
if($in2=~/^(\S+)\s+\S+\s+mRNA\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+.*Name=(\S+?);.*longest=1/){
$nam=$5;$chr=$1;$sta=$2;$end=$3;$str=$4;
$end{$chr}{$sta}{$nam}=$end;
$str{$nam}=$str;}}

$num=1;
sub sn{$sta_int{$chr}{$a}<=>$sta_int{$chr}{$b}}
for$chr(keys%end){
for$sta(keys%{$end{$chr}}){
for(keys%{$end{$chr}{$sta}}){
$end=$end{$chr}{$sta}{$_};
for($n=$sta;$n<=$end;$n++){
$gen{$chr}{$n}=1;
$nam{$chr}{$n}.=$_.$str{$_};
}
}
}
for($n=1;$n<=$seq{$chr};$n++){
if(!$gen{$chr}{$n}){
$sta_int{$chr}{$num}=$n and print$chr."\t".$n."\t".$num."\n" unless $sta_int{$chr}{$num};
#$nam_sta{$chr}{$num}=$nam{$chr}{$n-1} unless $sta_int{$chr}{$num} and $nam{$chr}{$n-1};
}
if ($gen{$chr}{$n} and !$gen{$chr}{$n-1}){
$end_int{$chr}{$num}=$n-1;
$nam_end{$chr}{$num}=$nam{$chr}{$n};
$num++;
}
$end_int{$chr}{$num}=$n if $n==$seq{$chr} and !$gen{$chr}{$n};
}

for(sort sn keys%{$sta_int{$chr}}){$n=$sta_int{$chr}{$_};
$c.=$_."\t".$chr."\t".$sta_int{$chr}{$_}."\t".$end_int{$chr}{$_}."\t".$nam{$chr}{$n-1}."_".$nam_end{$chr}{$_}."\n";
}
%gen=();%nam=();%sta_int=();%end_int=();%nam_sta=();%nam_end=();
}
open$out,'>',"/home/li/Other_studys/lncRNA_VVI/V5/classify/5intergenic_region"or die;
print$out($c);
