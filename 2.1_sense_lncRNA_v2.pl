#!/usr/bin/perl
open$in1,"/home/li/Other_studys/lncRNA_VVI/Genome/annotation/Vvinifera_457_v2.1.gene_exons.gff3"or die;
while($in2=<$in1>){
if($in2=~/^(\S+)\s+\S+\s+mRNA\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+.*Name=(\S+?);.*?longest=1/){
$long{$5}=1;}
if($in2=~/^(\S+)\s+\S+\s+exon\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+.*ID=(\S+)\.v2\.1\.exon\.(\d+)/){
next unless $long{$5};
$sta{$1}{$5.'.exon'.$6}=$2;$end{$1}{$5.'.exon'.$6}=$3;$str{$1}{$5.'.exon'.$6}=$4;
}
}
open$in1,"/home/li/Other_studys/lncRNA_VVI/V5/pipline_out/6static_map_non_rfam2"or die;
while($in2=<$in1>){
if($in2=~/^(\S+)\s+\S+(\S+)\s+(MSTRG\S+)/){print$in2; ##### version 3
$non{$3}=$1;$non_t2{$3}=1;
}
}
open$in1,"/home/li/Other_studys/lncRNA_VVI/V5/strand/1_str"or die;
while($in2=<$in1>){
if($in2=~/^(MSTRG\S+)/){print$in2; ##### version 3
$non_str{$1}=1;
}
}
open$in1,"/home/li/Other_studys/lncRNA_VVI/V5/48merged.gtf"or die;
while($in2=<$in1>){

if($in2=~/^(\S+)\s+\S+\s+exon\s+(\S+)\s+(\S+)\s+\S+\s+(\S+).*transcript_id\s+\"(\S+?)\"/){
next unless$non_t2{$5};
next unless$non_str{$5}; ### version3
#print$5."\n";
$chr=$1;$sta=$2;$end=$3;$str=$4;$nam=$5;
$mid=($sta+$end)/2; $r=abs($end-$sta+1)/2;
$str1{$nam}=$str;$chr1{$nam}=$chr;

for$t2(keys%{$sta{$chr}}){

$sta2=$sta{$chr}{$t2};$end2=$end{$chr}{$t2};$str2=$str{$chr}{$t2};
$mid2=($sta2+$end2)/2; $r2=abs($end2-$sta2+1)/2;

if(abs($mid-$mid2)<($r2+$r)){
$exon{$nam}=1;
$info{$nam}.="\tq:".$sta."-".$end.';;'."$t2:".$sta2.'-'.$end2.":".$str{$chr}{$t2};

if(($str2 eq $str or $str2 eq '.')){
$sen{$nam}+=1;
}
elsif(($str2 ne $str and $str2 ne '.')){
$ant{$nam}+=1;
}
}}}}

for(keys%sen){
$c.=$non{$_}."\t".$_."\t".$sen{$_}."\t".$chr1{$_}."\t".$str1{$_}.$info{$_}."\n";
}

for(keys%ant){
$d.=$non{$_}."\t".$_."\t".$ant{$_}."\t".$chr1{$_}."\t".$str1{$_}.$info{$_}."\n" unless $sen{$_};
}

for(keys%exon){
$e.=$_."\n";
}
open$out,'>',"/home/li/Other_studys/lncRNA_VVI/V5/classify/2_non_sense_v3"or die;
print$out($c);

open$out,'>',"/home/li/Other_studys/lncRNA_VVI/V5/classify/2_non_anti_v3"or die;
print$out($d);
open$out,'>',"/home/li/Other_studys/lncRNA_VVI/V5/classify/2_non_exon_overlap_v3"or die;
print$out($e);
