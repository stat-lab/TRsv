#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);

my $data_dir = "$Bin/../Data";

my $ref_file = '';

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
my $min_ins_str_mei = 200;

my $indel_rate = 10;

my $min_coverage = 80;
my $min_me_coverage = 60;

my $min_str_identity = 70;

my $dup_find = 1;

my $intersperse_dup_find = 1;

my $mei_find = 1;

my $bp_diff0 = 50;
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

if ($samtool_path ne ''){
    $ENV{PATH} = "$samtool_path:" . $ENV{PATH};
}

my $temp_dir = "$out_prefix.temp";
system ("mkdir $temp_dir") if (!-d $temp_dir);

my %STR;

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
    if ($type =~ /TR/){
        my $start = $pos - $bp_diff;
        my $end = $STR{$pos} + $bp_diff;
        my @sam = `samtools view $bam_file $target_chr:$start-$end | awk \'\$5==0\'` if ($bam_file =~ /\.bam$/);
        @sam = `samtools view --reference $ref_file $bam_file $target_chr:$start-$end | awk \'\$5==0\'` if ($bam_file =~ /\.cram$/);
        $type = 'INS' if ($type =~ /ins/);
        $type = 'DEL' if ($type =~ /del/);
        my $str_indel = 0;
        my $str_bp = 0;
        foreach (@sam){
            my @sline = split (/\t/, $_);
            my $spos = $sline[3];
            my $cigar = $sline[5];
            my $send = $spos;
            my %indel;
            my $count = 0;
            while ($cigar =~ /(\d+)([SHMIDX=])/g){
                my $sublen = $1;
                my $tag = $2;
                $count ++;
                if (($tag =~ /S|H/) and ($count == 1)){
                    if (($sublen >= $min_clip_len) and ($spos >= $start) and ($spos <= $pos + $bp_diff)){
                        $str_bp ++ if ($type eq 'INS');
                    }
                    next;
                }
                if (($tag eq 'M') or ($tag eq '=')){
                    $send += $sublen;
                }
                elsif ($tag eq 'D'){
                    $indel{$end} = -$sublen if ($sublen >= $max_dist) and ($send >= $start) and ($send <= $end);
                    $send += $sublen;
                }
                elsif ($tag eq 'I'){
                    $indel{$end} = $sublen if ($sublen >= $max_dist) and ($send >= $start) and ($send <= $end);
                }
                elsif ($tag eq 'X'){
                    $send += $sublen;
                }
                elsif (($tag =~ /S|H/) and ($sublen >= $min_clip_len)){
                    if (($send >= $end - $bp_diff * 2) and ($send <= $end)){
                        $str_bp ++ if ($type eq 'INS');
                    }
                }
                last if ($send > $end);
            }
            my $sum_indel = 0;
            foreach my $ipos (sort {$a <=> $b} keys %indel){
                $sum_indel += $indel{$ipos};
            }
            if (($type eq 'DEL') and ($sum_indel < 0)){
                $sum_indel = 0 - $sum_indel;
                if (($sum_indel / $len >= 0.67) and ($sum_indel / $len <= 1.5)){
                    $str_indel ++;
                }
            }
            elsif (($type eq 'INS') and ($sum_indel > 0)){
                if (($sum_indel / $len >= 0.67) and ($sum_indel / $len <= 1.5)){
                    $str_indel ++;
                }
            }
        }
        $str_indel += int ($str_bp * 0.5) if ($type eq 'INS');
        $sar = int ($str_indel / ($read + $str_indel) * 100 + 0.5) / 100;
    }
    elsif (($type =~ /INS/) and ($len > 0)){
        my $start = $pos - $bp_diff2;
        $start = $pos - $max_bp_diff if ($len >= 500);
        my $end = $pos + $bp_diff2;
        $end = $pos + $max_bp_diff if ($len >= 500);
        my @sam = `samtools view $bam_file $target_chr:$start-$end | awk \'\$5==0\'` if ($bam_file =~ /\.bam$/);
        @sam = `samtools view --reference $ref_file $bam_file $target_chr:$start-$end | awk \'\$5==0\'` if ($bam_file =~ /\.cram$/);
        my $ins = 0;
        foreach (@sam){
            my @sline = split (/\t/, $_);
            my $spos = $sline[3];
            my $cigar = $sline[5];
            my $send = $spos;
            my %ins;
            my $count = 0;
            while ($cigar =~ /(\d+)([SHMIDX=])/g){
                my $sublen = $1;
                my $tag = $2;
                $count ++;
                if (($tag =~ /S|H/) and ($count == 1)){
                    next;
                }
                if (($tag eq 'M') or ($tag eq '=')){
                    $send += $sublen;
                }
                elsif ($tag eq 'D'){
                    $send += $sublen;
                }
                elsif ($tag eq 'I'){
                    $ins{$end} = $sublen if ($sublen >= $max_dist) and ($send >= $start) and ($send <= $end);
                }
                elsif ($tag eq 'X'){
                    $send += $sublen;
                }
                elsif (($tag =~ /S|H/) and ($count > 1)){
                    next;
                }
                last if ($send > $end);
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
        }
        $sar = int ($ins / ($read - $bp + $ins) * 100 + 0.5) / 100 if ($read - $bp + $ins > 0);
    }
    elsif ($type eq 'DEL'){
        my $start = $pos - $bp_diff2;
        $start = int ($pos - $len * 0.5) if ($len >= 500);
        my $end = $pos + $len + $bp_diff2;
        $end = int ($pos + $len + $len * 0.5) if ($len >= 500);
        my @sam = `samtools view $bam_file $target_chr:$start-$end | awk \'\$5==0\'` if ($bam_file =~ /\.bam$/);
        @sam = `samtools view --reference $ref_file $bam_file $target_chr:$start-$end | awk \'\$5==0\'` if ($bam_file =~ /\.cram$/);
        my $del = 0;
        foreach (@sam){
            my @sline = split (/\t/, $_);
            my $spos = $sline[3];
            my $cigar = $sline[5];
            my $send = $spos;
            my %del;
            my $count = 0;
            while ($cigar =~ /(\d+)([SHMIDX=])/g){
                my $sublen = $1;
                my $tag = $2;
                $count ++;
                if (($tag =~ /S|H/) and ($count == 1)){
                    next;
                }
                if (($tag eq 'M') or ($tag eq '=')){
                    $send += $sublen;
                }
                elsif ($tag eq 'D'){
                    $del{$end} = $sublen if ($sublen >= $max_dist) and ($send >= $start) and ($send <= $end);
                    $send += $sublen;
                }
                elsif ($tag eq 'I'){
                    next;
                }
                elsif ($tag eq 'X'){
                    $send += $sublen;
                }
                elsif (($tag =~ /S|H/) and ($count > 1)){
                    next;
                }
                last if ($send > $end);
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
        }
        $sar = int ($del / ($read - $bp + $del) * 100 + 0.5) / 100 if ($read - $bp + $del > 0);
    }
    else{
        if (($len == 0) and ($info =~ /BPLEN-(\d+)/)){
            $len = $1;
        }
        my $pos2 = $pos + $len - 1;
        my $start1 = $pos - $bp_diff;
        my $end1 = $pos + $bp_diff;
        my $start2 = $pos2 - $bp_diff;
        my $end2 = $pos2 + $bp_diff;
        my @sam1 = `samtools view $bam_file $target_chr:$start1-$end1 | awk \'\$5==0\'` if ($bam_file =~ /\.bam$/);
        @sam1 = `samtools view --reference $ref_file $bam_file $target_chr:$start1-$end1 | awk \'\$5==0\'` if ($bam_file =~ /\.cram$/);
        my @sam2 = `samtools view $bam_file $target_chr:$start2-$end2 | awk \'\$5==0\'` if ($bam_file =~ /\.bam$/);
        @sam2 = `samtools view --reference $ref_file $bam_file $target_chr:$start2-$end2 | awk \'\$5==0\'` if ($bam_file =~ /\.cram$/);
        my $bp1 = 0;
        my $bp2 = 0;
        foreach (@sam1){
            my @sline = split (/\t/, $_);
            my $spos = $sline[3];
            my $cigar = $sline[5];
            my $send = $spos;
            my $count = 0;
            while ($cigar =~ /(\d+)([SHMIDX=])/g){
                my $sublen = $1;
                my $tag = $2;
                $count ++;
                if (($tag =~ /S|H/) and ($count == 1)){
                    if (($sublen >= $min_clip_len) and ($spos >= $start1) and ($spos <= $end1)){
                        $bp1 ++ if ($type =~ /INS-BP|^DUP|INV/);
                    }
                    next;
                }
                if (($tag eq 'M') or ($tag eq '=')){
                    $send += $sublen;
                }
                elsif ($tag eq 'D'){
                    $send += $sublen;
                }
                elsif ($tag eq 'I'){
                    next;
                }
                elsif ($tag eq 'X'){
                    $send += $sublen;
                }
                elsif (($tag =~ /S|H/) and ($sublen >= $min_clip_len)){
                    if (($send >= $start1) and ($send <= $end1)){
                        $bp1 ++ if ($type =~ /DEL|INV/);
                    }
                }
                last if ($send > $end1);
            }
        }
        foreach (@sam2){
            my @sline = split (/\t/, $_);
            my $spos = $sline[3];
            my $cigar = $sline[5];
            my $send = $spos;
            my $count = 0;
            while ($cigar =~ /(\d+)([SHMIDX=])/g){
                my $sublen = $1;
                my $tag = $2;
                $count ++;
                if (($tag =~ /S|H/) and ($count == 1)){
                    if (($sublen >= $min_clip_len) and ($spos >= $start2) and ($spos <= $end2)){
                        $bp2 ++ if ($type =~ /DEL|INV/);
                    }
                    next;
                }
                if (($tag eq 'M') or ($tag eq '=')){
                    $send += $sublen;
                }
                elsif ($tag eq 'D'){
                    $send += $sublen;
                }
                elsif ($tag eq 'I'){
                    next;
                }
                elsif ($tag eq 'X'){
                    $send += $sublen;
                }
                elsif (($tag =~ /S|H/) and ($sublen >= $min_clip_len)){
                    if (($send >= $start2) and ($send <= $end2)){
                        $bp2 ++ if ($type =~ /INS-BP|^DUP|INV/);
                    }
                }
                last if ($send > $end2);
            }
        }
        $sar = int (($bp1 + $bp2) / ($read * 2 + $bp1 + $bp2) * 100 + 0.5) / 100;
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


