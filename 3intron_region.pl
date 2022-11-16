#!/usr/bin/perl
open$in1,"/home/li/Other_studys/lncRNA_VVI/Genome/annotation/Vvinifera_457_v2.1.gene_exons_only.gff3"or die;
while($in2=<$in1>){
#if($in2=~/^(\S+)\s+\S+\s+mRNA\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+\S+\s+.*Name=(\S+?);.*longest=1/){
#$t2{$5}=1;print$5;
#}
if($in2=~/^(\S+)\s+\S+\s+exon\s+(\S+)\s+(\S+)\s+\S+\s+(\S+).*Parent=(\S+?)\.v2/){
#next unless $t2{$5};
print$in2;
$nam=$5;$chr=$1;$sta=$2;$end=$3;$str=$4;
if($nam1 eq $nam){
if($str eq '+'){
$c.=$nam."\t".$chr."\t".$str."\t".($end1+1)."\t".($sta-1)."\n";
print $in2 if $send1>$sta;
}
elsif($str eq '-'){
$c.=$nam."\t".$chr."\t".$str."\t".($end+1)."\t".($sta1-1)."\n";
print $in2 if $sta1<$end;
}
}
$nam1=$nam;$sta1=$sta;$end1=$end;
}
}
open$out,'>',"/home/li/Other_studys/lncRNA_VVI/Genome/annotation/vvi_v2.1_intron_location"or die;
print$out($c);
