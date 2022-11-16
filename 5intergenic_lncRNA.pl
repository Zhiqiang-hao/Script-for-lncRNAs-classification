#!/usr/bin/perl
open$in1,"/home/li/Other_studys/lncRNA_VVI/V5/classify/5intergenic_region"or die;
while($in2=<$in1>){chomp$in2;
@b=split"\t",$in2;
$sta=$b[2];$end=$b[3];$chr=$b[1];$nam=$b[4];
#print (join"\t",@b);
#print$chr."\n";
#print$sta."\t".$end."\n";
$end{$chr}{$nam}{$sta}=$end;
}
open$in1,"/home/li/Other_studys/lncRNA_VVI/V5/pipline_out/6static_map_non_rfam2"or die;
while($in2=<$in1>){
if($in2=~/^(\S+)\s+\S+(\S+)\s+(MSTRG\S+)/){
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
if($in2=~/^(\S+)\s+\S+\s+exon\s+(\S+)\s+(\S+)\s+\S+\s+(\S+).*transcript_id\s+\"(\S+?)\";\s+exon_number\s+\"(\d+)\"/){
next unless$non_t2{$5};
#next unless$non_str{$5}; ### version3
$chr=$1;$sta=$2;$end=$3;$str=$4;$nam=$5;$num=$6;$exon{$nam}=$num;
#print$sta."\t".$end."\n";
$str1{$nam}=$str;$chr1{$nam}=$chr;

for$t2(keys%{$end{$chr}}){
for$sta2(keys%{$end{$chr}{$t2}}){
$end2=$end{$chr}{$t2}{$sta2};#print$t2."\t".$sta2."\t".$end2."\n";
if($sta2<=$sta and $end<=$end2 and $sta<$end){
$inter{$nam}{$num}=1;
$info{$nam}.="\tq$num:".$sta.'-'.$end.";;$t2:".$sta2.'-'.$end2;
}
}}

}}

for$nam(keys%inter){
for(keys%{$inter{$nam}}){
if($inter{$nam}{$_}){
$tot{$nam}++;
}
}
}
for(keys%tot){
if($tot{$_} == $exon{$_}){
$c.=$non{$_}."\t".$_."\t".$tot{$_}."\t".$chr1{$_}."\t".$str1{$_}.$info{$_}."\n";
}
}
open$out,'>',"/home/li/Other_studys/lncRNA_VVI/V5/classify/5inter_lncRNA4"or die;
print$out($c);
