#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);

my $data_dir = "$Bin/../Data";

mmy $ref_file = '';

my $simple_repeat = '';

my $TE_fasta = '';

my $gap_bed = '';

my $bam_file = '';

my $samtool_path = '';

my $cores = 1;

my $out_prefix = '';

my $min_mapQ = 1;
my $min_mapQ_idup = 20;

my $max_SAR = 0.6;

my $min_ins_reads = 2;
my $min_del_reads = 2;
my $min_str_reads = 3;

my $min_VRR = 0.05;
my $min_str_vrr = 0.15;

my $max_depth_fold = 15;

my $min_maplen = 500;

my $chromium = 0;

my $include_secalign = 0;

my $targeted_seq = 0;

my $max_mismatch = 15;

my $min_indel_size = 50;
my $min_str_indel_size = 50;
my $min_ins_str_mei = 200;

my $indel_rate = 10;

my $min_coverage = 80;
my $min_me_coverage = 60;

my $min_str_identity = 70;

my $dup_find = 1;

my $intersperse_dup_find = 1;

my $mei_find = 1;

my $bp_diff0 = 20;
my $bp_diff = 100;
my $bp_diff2 = 200;
my $bp_diff3 = 150;
my $ins_bp_diff = 50;
my $max_bp_diff = 300;
my $inv_bp_diff = 6000;
my $min_overlap_rate = 0.7;
my $min_overlap_rate_eval = 0.5;
my $max_dist = 10;
my $max_dist2 = 15;
my $max_match_size = 200;
my $min_clip_len = 50;
my $max_insbp_diff = 200;
my $min_inv_size = 100;
my $max_del_size = 50000000;
my $max_inv_size = 50000000;
my $max_dup_len = 10000000;
my $min_dup_len = 50;
my $max_bp_sd_small_ins = 5;

my $Mbin_size = 1000000;
my $Mbin_size2 = 100000;

my $config = shift @ARGV;

my $target_chr = shift @ARGV;

my $step2_out = shift @ARGV;

open (FILE, $config) or die "$config is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($arg, $value) = split (/\t/, $line);
    if ($arg eq 'bam_file'){
        $bam_file = $value;
    }
    elsif ($arg eq 'ref_file'){
        $ref_file = $value;
    }
    elsif ($arg eq 'repeat_bed'){
        $simple_repeat = $value;
    }
    elsif ($arg eq 'te_fasta'){
        $TE_fasta = $value;
    }
    elsif ($arg eq 'gap_bed'){
        $gap_bed = $value;
    }
    elsif ($arg eq 'out_prefix'){
        $out_prefix = $value;
    }
    elsif ($arg eq 'min_len'){
        $min_indel_size = $value;
    }
    elsif ($arg eq 'min_str_len'){
        $min_str_indel_size = $value;
    }
    elsif ($arg eq 'min_ins_read'){
        $min_ins_reads = $value;
    }
    elsif ($arg eq 'min_del_read'){
        $min_del_reads = $value;
    }
    elsif ($arg eq 'min_tr_read'){
        $min_str_reads = $value;
    }
    elsif ($arg eq 'min_vrr'){
        $min_VRR = $value;
    }
    elsif ($arg eq 'min_tr_vrr'){
        $min_str_vrr = $value;
    }
    elsif ($arg eq 'max_dpf'){
        $max_depth_fold = $value;
    }
    elsif ($arg eq 'min_mapq'){
        $min_mapQ = $value;
    }
    elsif ($arg eq 'max_sar'){
        $max_SAR = $value;
    }
    elsif ($arg eq 'mismatch_rate'){
        $max_mismatch = $value;
    }
    elsif ($arg eq 'indel_rate'){
        $indel_rate = $value;
    }
    elsif ($arg eq 'incl_sec'){
        $include_secalign = $value;
    }
    elsif ($arg eq 'targeted'){
        $targeted_seq = $value;
    }
    elsif ($arg eq 'samtool_path'){
        $samtool_path = $value;
    }
}
close (FILE);

$max_dist = $min_indel_size * 0.9 if ($max_dist > $min_indel_size * 0.9);
$max_dist = $min_str_indel_size * 0.9 if ($max_dist > $min_str_indel_size * 0.9);

if ($samtool_path ne ''){
    $ENV{PATH} = "$samtool_path:" . $ENV{PATH};
}

my $temp_dir = "$out_prefix.temp";
system ("mkdir $temp_dir") if (!-d $temp_dir);

my %STR;
my %step2;
my %step2_type;
my %SAR;

open (FILE, $simple_repeat) or die "$simple_repeat is not found: $!\n" if ($simple_repeat !~ /\.gz$/);
open (FILE, "gzip -dc $simple_repeat |") or die "$simple_repeat is not found: $!\n" if ($simple_repeat =~ /\.gz$/);
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#|^$/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    next if ($chr ne $target_chr);
    my $pos = $line[1];
    my $end = $line[2];
    $STR{$pos} = $end;
}
close (FILE);

open (FILE, $step2_out) or die "$step2_out is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($chr, $pos, $type, $len, $info) = split (/\t/, $line);
    next if (!defined $info);
    my $read = $1 if ($info =~ /RN-(\d+)/);
    my $bp = 0;
    $bp = $1 if ($info =~ /BP-(\d+)/);
    if ($info =~ /RN1-(\d+)/){
        $read = $1;
        if ($info =~ /RN2-(\d+)/){
            $read += $1;
            $read = int ($read * 0.5 + 0.5);
        }
    }
    my $Mbin1 = int ($pos / $Mbin_size);
    my $end = $pos;
    if ($type =~ /TR/){
        $end = $STR{$pos};
    }
    elsif ($type =~ /INS/){
        if (($len == 0) and ($info =~ /BPLEN-(\d+)/)){
            $len = $1;
            $end = $pos + $len - 1;
        }
    }
    else{
        $end = $pos + $len - 1;
    }
    ${$step2{$Mbin1}}{$pos} = $end;
    ${${$step2_type{$Mbin1}}{$pos}}{$type} = $len;
    my $Mbin2 = int ($end / $Mbin_size);
    if ($Mbin2 > $Mbin1){
        ${${$step2_type{$Mbin2}}{$pos}}{$type} = $len;
        if (!exists ${$step2{$Mbin2}}{$pos}){
            ${$step2{$Mbin2}}{$pos} = $end;
        }
        else{
            my $end2 = ${$step2{$Mbin2}}{$pos};
            if ($end > $end2){
                ${$step2{$Mbin2}}{$pos} = $end;
            }
        }
    }
}
close (FILE);

open (FILE, "samtools view $bam_file $target_chr |") or die "$bam_file is not found:$!\n" if ($bam_file =~ /\.bam$/);
open (FILE, "samtools view --reference $ref_file $bam_file $target_chr |") or die "$bam_file is not found:$!\n" if ($bam_file =~ /\.cram$/);
while (my $sline = <FILE>){
    chomp $sline;
    my @sline = split (/\t/, $sline);
    next if ($sline[4] != 0);
    my $spos = $sline[3];
    my $send = $spos;
    my $cigar = $sline[5];
    my $flag = 0;
    my %hit_pos;
    while ($cigar =~ /(\d+)([MD=])/g){
        $send += $1;
    }
    my $Mbin1 = int ($spos / $Mbin_size);
    my $Mbin2 = int ($send / $Mbin_size);
    next if (!exists $step2{$Mbin1}) and (!exists $step2{$Mbin2});
    foreach my $pos (sort {$a <=> $b} keys %{$step2{$Mbin1}}){
        last if ($pos > $send);
        my $end = ${$step2{$Mbin1}}{$pos};
        next if ($end < $spos);
        $hit_pos{$pos} = 1;
        $flag = 1;
    }
    if (($flag == 0) and ($Mbin2 > $Mbin1)){
        foreach my $pos (sort {$a <=> $b} keys %{$step2{$Mbin2}}){
            last if ($pos > $send);
            my $end = ${$step2{$Mbin2}}{$pos};
            next if ($end < $spos);
            $hit_pos{$pos} = 1;
            $flag = 2;
        }
    }
    next if ($flag == 0);
    my $Mbin = $Mbin1;
    $Mbin = $Mbin2 if ($flag == 2);
    foreach my $pos (sort {$a <=> $b} keys %hit_pos){
        foreach my $type (keys %{${$step2_type{$Mbin}}{$pos}}){
            my $len = ${${$step2_type{$Mbin}}{$pos}}{$type};
            my $sar = 0;
            my $send2 = $spos;
            my $cigar2 = $cigar;
            if ($type =~ /TR/){
                my $start = $pos - $bp_diff0;
                my $end = $STR{$pos} + $bp_diff0;
                my $type2 = 'INS' if ($type =~ /ins/);
                $type2 = 'DEL' if ($type =~ /del/);
                my $str_indel = 0;
                my $str_bp = 0;
                my %indel;
                my $count = 0;
                while ($cigar2 =~ /(\d+)([SHMIDX=])/g){
                    my $sublen = $1;
                    my $tag = $2;
                    $count ++;
                    if (($tag =~ /S|H/) and ($count == 1)){
                        if (($sublen >= $min_clip_len) and ($spos >= $start) and ($spos <= $pos + $bp_diff0)){
                            $str_bp ++ if ($type2 eq 'INS');
                        }
                        next;
                    }
                    if (($tag eq 'M') or ($tag eq '=')){
                        $send2 += $sublen;
                    }
                    elsif ($tag eq 'D'){
                        $indel{$send2} = -$sublen if ($sublen >= $max_dist) and ($send2 >= $start) and ($send2 <= $end);
                        $send2 += $sublen;
                    }
                    elsif ($tag eq 'I'){
                        $indel{$send2} = $sublen if ($sublen >= $max_dist) and ($send2 >= $start) and ($send2 <= $end);
                    }
                    elsif ($tag eq 'X'){
                        $send2 += $sublen;
                    }
                    elsif (($tag =~ /S|H/) and ($sublen >= $min_clip_len)){
                        if (($send2 >= $end - $bp_diff * 2) and ($send2 <= $end)){
                            $str_bp ++ if ($type2 eq 'INS');
                        }
                    }
                    last if ($send2 > $end);
                }
                my $sum_indel = 0;
                foreach my $ipos (sort {$a <=> $b} keys %indel){
                    $sum_indel += $indel{$ipos};
                }
                if (($type2 eq 'DEL') and ($sum_indel < 0)){
                    $sum_indel = 0 - $sum_indel;
                    if (($sum_indel / $len >= 0.67) and ($sum_indel / $len <= 1.5)){
                        $str_indel ++;
                    }
                }
                elsif (($type2 eq 'INS') and ($sum_indel > 0)){
                    if (($sum_indel / $len >= 0.67) and ($sum_indel / $len <= 1.5)){
                        $str_indel ++;
                    }
                }
                $str_indel += int ($str_bp * 0.5) if ($type2 eq 'INS');
                ${$SAR{$pos}}{$type} += $str_indel;
            }
            elsif (($type =~ /INS/) and ($len > 0)){
                my $start = $pos - $bp_diff2;
                $start = $pos - $max_bp_diff if ($len >= 500);
                my $end = $pos + $bp_diff2;
                $end = $pos + $max_bp_diff if ($len >= 500);
                my %ins;
                my $ins = 0;
                my $count = 0;
                while ($cigar2 =~ /(\d+)([SHMIDX=])/g){
                    my $sublen = $1;
                    my $tag = $2;
                    $count ++;
                    if (($tag =~ /S|H/) and ($count == 1)){
                        next;
                    }
                    if (($tag eq 'M') or ($tag eq '=')){
                        $send2 += $sublen;
                    }
                    elsif ($tag eq 'D'){
                        $send2 += $sublen;
                    }
                    elsif ($tag eq 'I'){
                        $ins{$send2} = $sublen if ($sublen >= $max_dist) and ($send2 >= $start) and ($send2 <= $end);
                    }
                    elsif ($tag eq 'X'){
                        $send2 += $sublen;
                    }
                    elsif (($tag =~ /S|H/) and ($count > 1)){
                        next;
                    }
                    last if ($send2 > $end);
                }
                my $pre_ipos = 0;
                my $pre_ilen = 0;
                foreach my $ipos (sort {$a <=> $b} keys %ins){
                    my $ilen = $ins{$ipos};
                    my $distance = $ipos - $pre_ipos + 1;
                    my $flag = 0;
                    if (($pre_ipos > 0) and ($distance <= 50) and ($distance < $pre_ilen) and ($distance < $ilen)){
                        $flag = 1;
                    }
                    elsif (($pre_ipos > 0) and ($distance <= 100) and ($distance < $pre_ilen * 0.2) and ($distance < $ilen * 0.2)){
                        $flag = 1;
                    }
                    if ($flag == 1){
                        if ($ilen >= $pre_ilen){
                            $ins{$ipos} += $pre_ilen;
                            delete $ins{$pre_ipos};
                            $pre_ilen = $ilen + $pre_ilen;
                            $pre_ipos = $ipos;
                            next;
                        }
                        else{
                            $ins{$pre_ipos} += $ilen;
                            delete $ins{$ipos};
                            $pre_ilen = $ilen + $pre_ilen;
                            next;
                        }
                    }
                    $pre_ipos = $ipos;
                    $pre_ilen = $ilen;
                }
                foreach my $ipos (sort {$a <=> $b} keys %ins){
                    my $ilen = $ins{$ipos};
                    if (($ilen / $len >= 0.67) and ($ilen / $len <= 1.5)){
                        $ins ++;
                    }
                }
                ${$SAR{$pos}}{$type} += $ins;
            }
            elsif ($type eq 'DEL'){
                my $start = $pos - $bp_diff2;
                $start = int ($pos - $len * 0.5) if ($len >= 500);
                my $end = $pos + $len + $bp_diff2;
                $end = int ($pos + $len + $len * 0.5) if ($len >= 500);
                my %del;
                my $del = 0;
                my $count = 0;
                while ($cigar2 =~ /(\d+)([SHMIDX=])/g){
                    my $sublen = $1;
                    my $tag = $2;
                    $count ++;
                    if (($tag =~ /S|H/) and ($count == 1)){
                        next;
                    }
                    if (($tag eq 'M') or ($tag eq '=')){
                        $send2 += $sublen;
                    }
                    elsif ($tag eq 'D'){
                        $del{$send2} = $sublen if ($sublen >= $max_dist) and ($send2 >= $start) and ($send2 <= $end);
                        $send2 += $sublen;
                    }
                    elsif ($tag eq 'I'){
                        next;
                    }
                    elsif ($tag eq 'X'){
                        $send2 += $sublen;
                    }
                    elsif (($tag =~ /S|H/) and ($count > 1)){
                        next;
                    }
                    last if ($send2 > $end);
                }
                my $pre_dpos = 0;
                my $pre_dlen = 0;
                foreach my $dpos (sort {$a <=> $b} keys %del){
                    my $dlen = $del{$dpos};
                    my $pre_dend = $pre_dpos + $dlen - 1;
                    my $distance = $dpos - $pre_dend + 1;
                    my $flag = 0;
                    if (($pre_dpos > 0) and ($distance <= 500) and ($distance < $pre_dlen * 0.5) and ($distance < $dlen * 0.5)){
                        $flag = 1;
                    }
                    if ($flag == 1){
                        if ($dlen >= $pre_dlen){
                            $del{$dpos} = $dpos + $dlen - $pre_dpos;
                            delete $del{$pre_dpos};
                            $pre_dlen = $dpos + $dlen - $pre_dpos;
                            $pre_dpos = $dpos;
                            next;
                        }
                        else{
                            $del{$pre_dpos} = $dpos + $dlen - $pre_dpos;
                            delete $del{$dpos};
                            $pre_dlen = $dpos + $dlen - $pre_dpos;
                            next;
                        }
                    }
                    $pre_dpos = $dpos;
                    $pre_dlen = $dlen;
                }
                foreach my $dpos (sort {$a <=> $b} keys %del){
                    my $dlen = $del{$dpos};
                    if (($dlen / $len >= 0.67) and ($dlen / $len <= 1.5)){
                        $del ++;
                    }
                }
                ${$SAR{$pos}}{$type} += $del;
            }
            else{
                my $pos2 = $pos + $len - 1;
                my $start1 = $pos - $bp_diff;
                my $end1 = $pos + $bp_diff;
                my $start2 = $pos2 - $bp_diff;
                my $end2 = $pos2 + $bp_diff;
                my $count = 0;
                my $bp = 0;
                while ($cigar2 =~ /(\d+)([SHMIDX=])/g){
                    my $sublen = $1;
                    my $tag = $2;
                    $count ++;
                    if (($tag =~ /S|H/) and ($count == 1) and ($type =~ /INS-BP|^DUP|INV/)){
                        if (($sublen >= $min_clip_len) and ($spos >= $start1) and ($spos <= $end1)){
                            $bp ++;
                        }
                        next;
                    }
                    elsif (($tag =~ /S|H/) and ($count == 1) and ($type =~ /DEL|INV/)){
                        if (($sublen >= $min_clip_len) and ($spos >= $start2) and ($spos <= $end2)){
                            $bp ++;
                        }
                        next;
                    }
                    elsif (($tag eq 'M') or ($tag eq '=')){
                        $send2 += $sublen;
                    }
                    elsif ($tag eq 'D'){
                        $send2 += $sublen;
                    }
                    elsif ($tag eq 'I'){
                        next;
                    }
                    elsif ($tag eq 'X'){
                        $send2 += $sublen;
                    }
                    elsif (($tag =~ /S|H/) and ($sublen >= $min_clip_len) and ($type =~ /DEL|INV/)){
                        if (($send2 >= $start1) and ($send2 <= $end1)){
                            $bp ++;
                        }
                    }
                    elsif (($tag =~ /S|H/) and ($sublen >= $min_clip_len) and ($type =~ /INS-BP|^DUP|INV/)){
                        if (($send2 >= $start2) and ($send2 <= $end2)){
                            $bp ++;
                        }
                    }
                    last if ($send2 > $end2);
                }
                ${$SAR{$pos}}{$type} += $bp;
            }
        }
    }
}
close (FILE);

my $step3_out = $step2_out . 2;
open (OUT, "> $step3_out");
open (FILE, $step2_out) or die "$step2_out is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($chr, $pos, $type, $len, $info) = split (/\t/, $line);
    next if (!defined $info);
    my $read = $1 if ($info =~ /RN-(\d+)/);
    my $bp = 0;
    $bp = $1 if ($info =~ /BP-(\d+)/);
    if ($info =~ /RN1-(\d+)/){
        $read = $1;
        if ($info =~ /RN2-(\d+)/){
            $read += $1;
            $read = int ($read * 0.5 + 0.5);
        }
    }
    my $sar = 0;
    if (exists ${$SAR{$pos}}{$type}){
        my $var = ${$SAR{$pos}}{$type};
        if ($type =~ /TR/){
            $sar = int ($var / ($read + $var) * 100 + 0.5) / 100;
        }
        elsif (($type =~ /INS/) and ($len > 0)){
            $sar = int ($var / ($read - $bp + $var) * 100 + 0.5) / 100 if ($read - $bp + $var > 0);
        }
        elsif ($type eq 'DEL'){
            $sar = int ($var / ($read - $bp + $var) * 100 + 0.5) / 100 if ($read - $bp + $var> 0);
        }
        else{
            $sar = int ($var / ($read * 2 + $var) * 100 + 0.5) / 100;
        }
    }
    if ($line !~ /SAR-/){
        $line .= ";SAR-$sar";
    }
    else{
        $line =~ s/SAR-.+/SAR-$sar/;
    }
    print OUT "$line\n";
}
close (FILE);
close (OUT);


my $step2_wc = `wc -l $step2_out`;
chomp $step2_wc;
$step2_wc = $1 if ($step2_wc =~ /^(\d+)/);
my $step3_wc = `wc -l $step3_out`;
chomp $step3_wc;
$step3_wc = $1 if ($step3_wc =~ /^(\d+)/);

if ($step2_wc != $step3_wc){
    die "$target_chr: step2 out wc ($step2_wc) not equal to step3 out wc ($step3_wc)\n";
}
else{
    system ("rm $step2_out");
    system ("mv $step3_out $step2_out");
}
