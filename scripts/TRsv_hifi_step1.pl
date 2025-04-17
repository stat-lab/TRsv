#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);

my $data_dir = "$Bin/../Data";

my $ref_file = '';
my $ref_index = '';

my $simple_repeat = '';

my $TE_fasta = '';

my $gap_bed = '';

my $bam_file = '';

my $samtool_path = '';

my $cores = 1;

my $out_prefix = '';

my $min_mapQ = 1;
my $min_mapQ_idup = 20;

my $min_ins_reads = 2;
my $min_del_reads = 2;
my $min_str_reads = 3;

my $min_VRR = 0.05;
my $min_str_vrr = 0.15;

my $str_max_len_rate = 1.5;

my $max_depth_fold = 15;

my $min_maplen = 200;   # ONT/CLR: 500

my $chromium = 0;

my $include_secalign = 0;

my $targeted_seq = 0;

my $exclude_bed = '';

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

my $bp_diff0 = 50;
my $bp_diff = 100;
my $bp_diff2 = 200;
my $bp_diff3 = 150;
my $ins_bp_diff2 = 20;
my $ins_bp_diff = 50;
my $max_bp_diff = 300;
my $inv_bp_diff = 6000;
my $min_overlap_rate = 0.5;
my $min_overlap_rate2 = 0.7;
my $max_dist = 10;
my $max_dist2 = 15;
my $max_match_size = 200;
my $min_clip_len = 20;      # ONT/CLR: 50
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
    elsif ($arg eq 'ref_index'){
        $ref_index = $value;
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
    elsif ($arg eq 'min_tr_len'){
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
    elsif ($arg eq 'exclude'){
        $exclude_bed = $value;
    }
    elsif ($arg eq 'max_tr_rate'){
        $str_max_len_rate = $value;
    }
    elsif ($arg eq 'samtool_path'){
        $samtool_path = $value;
    }
}
close (FILE);

if ($samtool_path ne ''){
    $ENV{PATH} = "$samtool_path:" . $ENV{PATH};
}

if ($max_dist > $min_indel_size){
    $max_dist = int ($min_indel_size * 0.8);
    $max_dist = $min_indel_size - 1 if ($min_indel_size < 5);
}
if ($max_dist > $min_str_indel_size){
    $max_dist = int ($min_str_indel_size * 0.8);
    $max_dist = $min_str_indel_size - 1 if ($min_str_indel_size < 5);
}
if ($max_dist2 > $min_str_indel_size){
    $max_dist2 = int ($min_str_indel_size * 0.8);
    $max_dist2 = $min_str_indel_size - 1 if ($min_str_indel_size < 5);
}

my $str_min_len_rate = int (1 / $str_max_len_rate * 100 + 0.5) / 100;

my $temp_dir = "$out_prefix.temp";
system ("mkdir $temp_dir") if (!-d $temp_dir);

my %chr_len;
my %STR;
my %STR2;
my %STRovl;
my %repeat;
my %STR_ME;
my %STR_motif;
my %gap_bp;
my %exclude_region;
my $hg19_flag = 0;
my $total_STR_len = 0;
my $exclude_flag = 0;

open (FILE, $ref_index) or die "$ref_index is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($chr, $len) = split (/\t/, $line);
    $chr_len{$chr} = $len;
    if ($chr =~ /^chr/){
        $hg19_flag = 1;
    }
}
close (FILE);

if ($gap_bed ne ''){
    open (FILE, "gzip -dc $gap_bed |") or die "$gap_bed is not found:$!\n" if ($gap_bed =~ /\.gz$/);
    open (FILE, $gap_bed) or die "$gap_bed is not found:$!\n" if ($gap_bed !~ /\.gz$/);
    while (my $line = <FILE>){
        chomp $line;
        my @line = split (/\s+/, $line);
        my $chr = $line[0];
        my $pos = $line[1];
        my $end = $line[2];
        for (my $i = $pos - 10; $i <= $pos + 2; $i++){
            ${$gap_bp{$chr}}{$i} = 1;
        }
        for (my $i = $end - 2; $i <= $end + 10; $i++){
            ${$gap_bp{$chr}}{$i} = 1;
        }
    }
    close (FILE);
}

if (-f $exclude_bed){
    open (FILE, $exclude_bed) or die "$exclude_bed is not found: $!\n" if ($exclude_bed !~ /\.gz$/);
    open (FILE, "gzip -dc $exclude_bed |") or die "$exclude_bed is not found: $!\n" if ($exclude_bed =~ /\.gz$/);
    while (my $line = <FILE>){
        chomp $line;
        my ($chr, $start, $end) = split (/\s+/, $line);
        next if ($chr ne $target_chr);
        $exclude_region{$start} = $end;
    }
    close (FILE);
    $exclude_flag = 1 if (scalar keys %exclude_region > 0);
}

my $pre_sid = 0;
my $pre_send = 0;

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
    my $len = $end - $pos + 1;
    if ($exclude_flag == 1){
        my $flag = 0;
        foreach my $xpos (sort {$a <=> $b} keys %exclude_region){
            last if ($xpos > $end);
            my $xend = $exclude_region{$xpos};
            next if ($xend < $pos);
            my $overlap = 0;
            if (($xpos <= $pos) and ($xend >= $end)){
                $flag = 1;
                last;
            }
            elsif (($xpos >= $pos) and ($xpos <= $end)){
                $overlap = $end - $xpos;
            }
            elsif (($xend >= $pos) and ($xend <= $end)){
                $overlap = $xend - $pos;
            }
            elsif (($xpos <= $pos) and ($xend >= $end)){
                $overlap = $xend - $xpos;
            }
            if ($overlap >= $len * 0.5){
                $flag = 1;
                last;
            }
        }
        if ($flag == 1){
            next;
        }
    }
    my $id = $line[3];
    my $me = $line[6];
    my $motif = $line[7];
    my $motif_size = $line[4];
    $STR{$pos} = $id;
    $STR2{$id} = "$pos=$end=$motif_size";
    $STR_ME{$id} = $me if ($me ne '-');
    $STR_motif{$id} = $motif;
    $total_STR_len += $end - $pos + 1;
    my $Mbin_start = int ($pos / $Mbin_size);
    my $Mbin_end = int ($end / $Mbin_size);
    if ($Mbin_start == $Mbin_end){
        ${$repeat{$Mbin_start}}{$pos} = $end;
    }
    else{
        my $Mbin = $Mbin_start;
        while ($Mbin <= $Mbin_end){
            ${$repeat{$Mbin}}{$pos} = $end;
            $Mbin ++;
        }
    }
    if ($pos <= $pre_send){
        if (!exists $STRovl{$pre_sid}){
            $STRovl{$pre_sid} = $id;
        }
        else{
            $STRovl{$pre_sid} .= ",$id";
        }
        if (!exists $STRovl{$id}){
            $STRovl{$id} = $pre_sid;
        }
        else{
            $STRovl{$id} .= ",$pre_sid";
        }
    }
    $pre_sid = $id;
    $pre_send = $end;
}
close (FILE);

my $chr = $target_chr;
my %bam_indel;
my %bam_bp1;
my %bam_bp2;
my %bam_2bp;
my %read_bp1;
my %read_bp2;
my %read_2bp;
my %read_2bp_rpos;
my %read_align5;
my %read_align3;
my %bam_ins;
my %bam_del;
my %bam_dup;
my %bam_inv;
my %bam_rep;
my %bam_str;
my %bam_str2;
my %bam_str_pos;
my %bam_str_cov;
my %bam_str_cov2;
my %bam_str_insinfo;
my %bam_str_insinfo2;
my %bam_indel_str;
my %bam_ins_str_read;
my %bam_str_bp1;
my %bam_str_bp2;
my %bam_ins_Q0;
my %overlap_str_opt;

my %gap;
my %DP;
my %DP_Q0;
my %bam_indel_Q0cov;
my %shortmap;
my $max_sublen = 100;
my $chr_len = $chr_len{$chr};
my $total_align = 0;
my $total_realign = 0;
my %align_reads;
my $sec_align_num = 0;
my $total_del_it = 0;
my $delete_mapq = 0;
my @type_str = ('DEL', 'INS', 'BP');
my $sum_align = 0;

open (FILE, "samtools view $bam_file $chr |") or die "$bam_file is not found: $!\n" if ($bam_file =~ /\.bam$/);
open (FILE, "samtools view -S $bam_file $chr |") or die "$bam_file is not found: $!\n" if ($bam_file =~ /\.sam$/);
open (FILE, "samtools view --reference $ref_file $bam_file $chr |") or die "$bam_file is not found: $!\n" if ($bam_file =~ /\.cram$/);
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^\@/);
    my @line = split (/\t/, $line);
    my $read_id = $line[0];
    my $pos = $line[3];
    my $cigar = $line[5];
    next if ($cigar eq '*');
    $total_align ++;
    $align_reads{$read_id} ++;
    my $tag = $line[1];
    my $mapQ = $line[4];
    my $read_pos = 0;
    my $map_len = 0;
    my %ins;
    my %del;
    my %rep;
    my %del_rpos;
    my %Dpos;
    last if ($chr !~ /^c*h*r*[\dXY]+$/);
    my $end = $pos;
    my $count = 0;
    my $S5_flag = 0;
    my $S3_flag = 0;
    my $flag = 0;
    my $strand = 'F';
    my $sec_align = 0;
    my $exclude_overlap_flag = 0;
    if ($tag >= 2048){
        $tag -= 2048;
    }
    if ($tag >= 1024){
        $tag -= 1024;
    }
    if ($tag >= 512){
        $tag -= 512;
    }
    if ($tag >= 256){
        $tag -= 256;
        $sec_align = 1;
    }
    if ($tag >= 128){
        $tag -= 128;
    }
    if ($tag >= 64){
        $tag -= 64;
    }
    if ($tag >= 32){
        $tag -= 32;
    }
    if ($tag >= 16){
        $tag -= 16;
        $strand = 'R';
    }
    if ($sec_align == 1){
        $sec_align_num ++;
        $delete_mapq ++ if ($mapQ < $min_mapQ);
        next if ($include_secalign == 0);
    }
    if ($mapQ < $min_mapQ){
        $delete_mapq ++;
        while ($cigar =~ /(\d+)([SHMIDX=])/g){
            my $sublen = $1;
            my $tag = $2;
            if ($tag =~ /S|H/){
                $read_pos += $sublen;
            }
            elsif (($tag eq 'M') or ($tag eq '=')){
                $end += $sublen;
                $read_pos += $sublen;
            }
            elsif ($tag eq 'D'){
                $end += $sublen;
            }
            elsif ($tag eq 'I'){
                my $end2 = int ($end / 10 + 0.5) * 10;
                push @{$bam_ins_Q0{$end2}}, "$read_id=$strand=$sublen=$read_pos=$end" if ($sublen >= 1000);
                $read_pos += $sublen;
            }
        }
        $read_align5{$read_id} = $pos;
        $read_align3{$read_id} = $end;
        next;
    }
    my $sublen_5S = 0;
    my $last_cig = '';
    while ($cigar =~ /(\d+)[MIX=]/g){
        $map_len += $1;
    }
    next if ($map_len < $min_maplen);
    $sum_align += $map_len;

    while ($cigar =~ /(\d+)([SHMIDX=])/g){
        my $sublen = $1;
        my $tag = $2;
        $count ++;
        if (($tag =~ /S|H/) and ($count == 1)){
            if ($sublen >= $min_clip_len){
                $S5_flag = $sublen;
                $sublen_5S = $sublen;
                push @{$bam_bp2{$pos}}, "$read_id=$strand=$sublen=$tag=$map_len" if (!exists ${$gap_bp{$chr}}{$pos});
                ${$read_bp2{$read_id}}{$pos} = "$map_len=$sublen=$cigar=$strand" if (!exists ${$gap_bp{$chr}}{$pos});
                ${$shortmap{$read_id}}{$pos} = 1 if ($map_len < 1000);
            }
            $read_pos += $sublen;
            next;
        }
        if (($tag eq 'M') or ($tag eq '=')){
            $end += $sublen;
            $read_pos += $sublen;
        }
        elsif ($tag eq 'D'){
            $del{$end} = $sublen if ($sublen >= $max_dist);
            $del_rpos{$end} = $read_pos if ($sublen >= $max_dist);
            if ($sublen >= 40){
                my $end2 = int ($end / 50 + 0.5) * 50;
                my $end3 = $end + $sublen;
                $end3 = $end2 if ($end2 > $end3);
                while ($end2 <= $end3){
                    $Dpos{$end2} = 1;
                    $end2 += 50;
                }
            }
            $end += $sublen;
        }
        elsif ($tag eq 'I'){
            $ins{$end} = "$read_id=$strand=$sublen=$read_pos" if ($sublen >= $max_dist);
            $read_pos += $sublen;
        }
        elsif ($tag eq 'X'){
            if ($sublen >= $max_dist){
                my $prefix = $`;
                my $suffix = $';
                if (($prefix =~ /[M=]$/) and ($suffix =~ /^\d+[M=]/)){
                    $rep{$end} = $sublen;
                }
            }
            $end += $sublen;
            $read_pos += $sublen;
        }
        elsif (($tag =~ /S|H/) and ($sublen >= $min_clip_len)){
            $S3_flag = $sublen;
            if (($map_len >= $min_maplen) and (!exists ${$gap_bp{$chr}}{$end})){
                push @{$bam_bp1{$end}}, "$read_id=$strand=$sublen=$tag";
                ${$read_bp1{$read_id}}{$end} = "$read_pos=$sublen=$cigar=$strand";
            }
            ${$shortmap{$read_id}}{$end} = 1 if ($map_len < 1000);
            if (($S5_flag >= $max_sublen) and ($S3_flag >= $max_sublen) and (!exists ${$gap_bp{$chr}}{$end})){
                my $len = $end - $pos + 1;
                ${${$read_2bp{$read_id}}{$pos}}{$end} = 1;
                push @{${$read_2bp_rpos{$read_id}}{$end}}, "$S5_flag=$S3_flag=$read_pos=$strand=$pos";
                push @{$bam_2bp{$pos}}, "$end=$read_id=$strand=$S5_flag=$S3_flag";
            }
        }
        $last_cig = $sublen . $tag;
    }
    if ($exclude_flag == 1){
        foreach my $xpos (sort {$a <=> $b} keys %exclude_region){
            last if ($xpos > $end);
            my $xend = $exclude_region{$xpos};
            next if ($xend < $pos);
            if (($pos >= $xpos) and ($end <= $xend)){
                $exclude_overlap_flag = 3;
            }
            elsif (($pos >= $xpos) and ($pos <= $xend)){
                $exclude_overlap_flag = 1;
            }
            elsif (($end >= $xpos) and ($end <= $xend)){
                $exclude_overlap_flag = 2;
            }
        }
        if (($exclude_overlap_flag == 1) or ($exclude_overlap_flag == 3)){
            delete $bam_bp2{$pos};
            delete ${$read_bp2{$read_id}}{$pos};
        }
        if (($exclude_overlap_flag == 2) or ($exclude_overlap_flag == 3)){
            delete $bam_bp1{$end};
            delete ${$read_bp1{$read_id}}{$end};
        }
        if ($exclude_overlap_flag > 0){
            if (exists $bam_2bp{$pos}){
                delete $bam_2bp{$pos};
                delete ${$read_2bp_rpos{$read_id}}{$end};
                delete ${$read_2bp{$read_id}}{$pos};
                if (scalar keys %{$read_2bp_rpos{$read_id}} == 0){
                    delete $read_2bp_rpos{$read_id};
                }
                if (scalar keys %{$read_2bp{$read_id}} == 0){
                    delete $read_2bp{$read_id};
                }
            }
        }
        next if ($exclude_overlap_flag == 3);
    }
    $read_align5{$read_id} = $pos;
    $read_align3{$read_id} = $end;
    my $pos2 = int ($pos / 50 + 0.5) * 50;
    while ($pos2 <= $end){
        my $Mbin = int ($pos2 / $Mbin_size);
        my $Mbin_res = $pos2 % $Mbin_size;
        if (!exists $Dpos{$pos2}){
            ${$DP{$Mbin}}{$Mbin_res} ++;
        }
        $pos2 += 50;
    }
    my $Mbin_s = int ($pos / $Mbin_size);
    my $Mbin_e = int ($end / $Mbin_size);
    my %BP_str;
    my %ovl_str;
    my %str_bp_flag;
    my %str_bp_flag2;
    foreach my $sbin (sort {$a <=> $b} keys %repeat){      # determine read coverage and the number of break ends within TRs
        last if ($sbin > $Mbin_e);
        next if ($sbin < $Mbin_s);
        foreach my $spos (sort {$a <=> $b} keys %{$repeat{$sbin}}){
            last if ($spos > $end + 20);
            my $send = ${$repeat{$sbin}}{$spos};
            next if ($send < $pos - 20);
            my $sid = $STR{$spos};
            $ovl_str{$sid} = 1;
            my $bp5_flag = 0;
            my $bp3_flag = 0;
            if (($pos >= $spos - 20) and ($pos <= $send + 20)){
                $BP_str{$sid} = 5;
                if (($S5_flag > 0) and (!exists ${$gap_bp{$chr}}{$pos})){
                    push @{${$bam_str{$sid}}{'BP5'}}, $pos;
                    $bam_str_bp2{$pos} = $sid;
                    $bp5_flag = 1;
                }
            }
            if (($end >= $spos - 20) and ($end <= $send + 20)){
                $BP_str{$sid} = 3;
                if (($S3_flag > 0) and (!exists ${$gap_bp{$chr}}{$end})){
                    push @{${$bam_str{$sid}}{'BP3'}}, $end;
                    $bam_str_bp1{$end} = $sid;
                    $bp3_flag = 1;
                }
            }
            if (($bp5_flag == 1) and ($bp3_flag == 1)){
                push @{${$bam_str{$sid}}{'BP53'}}, $pos if (!exists ${$gap_bp{$chr}}{$pos}) and (!exists ${$gap_bp{$chr}}{$end});
            }
            if (($bp5_flag == 1) or ($bp3_flag == 1)){
                $str_bp_flag{$sid} = 1;
            }
            if (($pos < $spos - 20) and ($end > $send + 20)){
                $str_bp_flag2{$sid} = 1;
            }
            elsif (($pos >= $spos - 20) and ($pos <= $send + 20) and (exists ${$gap_bp{$chr}}{$pos})){
                $str_bp_flag2{$sid} = 1;
            }
            elsif (($end >= $spos - 20) and ($end <= $send + 20) and (exists ${$gap_bp{$chr}}{$end})){
                $str_bp_flag2{$sid} = 1;
            }
        }
    }
    if (scalar keys %ovl_str > 0){
        foreach my $sid (keys %ovl_str){
            my ($spos, $send) = split (/=/, $STR2{$sid});
            my $slen = $send - $spos + 1;
            my $overlap = 0;
            if (($pos <= $spos) and ($end >= $send)){
                $overlap = $slen;
            }
            elsif (($pos >= $spos) and ($end <= $send)){
                $overlap = $end - $pos + 1;
            }
            elsif (($pos >= $spos) and ($pos <= $send)){
                $overlap = $send - $pos + 1;
            }
            elsif (($end >= $spos) and ($end <= $send)){
                $overlap = $end - $spos + 1;
            }
            my $ovl_rate = int ($overlap / $slen * 10 + 0.5) / 10;
            $ovl_rate = 1 if ($ovl_rate >= 0.8);
            if ($ovl_rate > 0){
                $bam_str_cov{$sid} += $ovl_rate;
                $bam_str_cov2{$sid} += $ovl_rate if (exists $str_bp_flag2{$sid});
            }
        }
    }
    my %STRins;
    my %STRins_id;
    my %STRdel;
    my %STRdel_id;
    my %STR_id;
    my %ipos_str;
    if (scalar keys %ins > 0){
        if ($chromium == 1){        # delete INSs whose BPs are flanked with a stretch of N (N6)
            my $rseq = $line[9];
            foreach my $ipos (keys %ins){
                my ($id1, $strand1, $ilen, $rpos) = split (/=/, $ins{$ipos});
                if ($cigar =~ /^(\d+)H/){
                    $rpos -= $1;
                }
                my $rend = $rpos + $ilen;
                my $rseq_1 = '';
                my $rseq_2 = '';
                $rseq_1 = substr ($rseq, $rpos - 19, 18) if ($rpos >= 19);
                $rseq_2 = substr ($rseq, $rend, 18) if (length ($rseq) - $rend > 18);
                if (($rseq_1 =~ /N{6}/) or ($rseq_2 =~ /N{6}/)){
                    delete $ins{$ipos};
                }
            }
        }
        foreach my $ipos (sort {$a <=> $b} keys %ins){      # assign INSs (>= 20 bp) within TRs
            if ($exclude_overlap_flag > 0){
                my $ovl_flag = 0;
                foreach my $xpos (sort {$a <=> $b} keys %exclude_region){
                    last if ($xpos > $ipos);
                    my $xend = $exclude_region{$xpos};
                    next if ($xend < $ipos);
                    $ovl_flag = 1;
                    last;
                }
                if ($ovl_flag == 1){
                    delete $ins{$ipos};
                    next;
                }
            }
            my $Mbin_start = int ($ipos / $Mbin_size);
            if (exists $repeat{$Mbin_start}){
                my ($id1, $strand1, $ilen, $rpos) = split (/=/, $ins{$ipos});
                my %match_str;
                foreach my $spos (sort {$a <=> $b} keys %{$repeat{$Mbin_start}}){
                    last if ($spos > $ipos + 20);
                    my $send = ${$repeat{$Mbin_start}}{$spos};
                    next if ($send < $ipos - 20);
                    my $str = $STR{$spos};
                    $match_str{$str} = 1;
                }
                if (scalar keys %match_str > 0){
                    if (scalar keys %match_str > 1){
                        my %str_within;
                        foreach my $strid (keys %match_str){
                            my ($spos, $send) = split (/=/, $STR2{$strid});
                            if (($ipos >= $spos) and ($ipos < $send)){
                                $str_within{$strid} = 1;
                            }
                        }
                        if (scalar keys %str_within > 0){
                            foreach my $strid (keys %match_str){
                                delete $match_str{$strid} if (!exists $str_within{$strid});
                            }
                        }
                        if (scalar keys %match_str > 1){
                            my %str_diff;
                            foreach my $strid (keys %match_str){
                                my ($spos, $send) = split (/=/, $STR2{$strid});
                                my $distance = abs ($ipos - $spos);
                                $distance = abs ($send - $ipos + 1) if (abs ($send - $ipos + 1) < $distance);
                                $str_diff{$strid} = $distance;
                            }
                            my $opt_str = '';
                            my $opt_str2 = '';
                            if (scalar keys %str_within == 0){
                                foreach my $strid (sort {$str_diff{$a} <=> $str_diff{$b}} keys %str_diff){
                                    if ($opt_str eq ''){
                                        $opt_str = $strid;
                                        next;
                                    }
                                    if ($opt_str2 eq ''){
                                        $opt_str2 = $strid;
                                        last;
                                    }
                                }
                            }
                            else{
                                foreach my $strid (sort {$str_diff{$b} <=> $str_diff{$a}} keys %str_diff){
                                    if ($opt_str eq ''){
                                        $opt_str = $strid;
                                        next;
                                    }
                                    if ($opt_str2 eq ''){
                                        $opt_str2 = $strid;
                                        last;
                                    }
                                }
                            }
                            $overlap_str_opt{$opt_str} = $opt_str2;
                            foreach my $strid (keys %match_str){
                                if ($strid ne $opt_str){
                                    delete $match_str{$strid};
                                }
                            }
                        }
                    }
                    foreach my $match_str (keys %match_str){
                        $STR_id{$match_str} = 1;
                        $STRins{$ipos} = $match_str;
                        ${$STRins_id{$match_str}}{$ipos} = $ilen;
                        if (!exists $ipos_str{$match_str}){
                            $ipos_str{$match_str} = "$rpos-$ilen-$strand";
                        }
                        else{
                            $ipos_str{$match_str} .= ",$rpos-$ilen-$strand";
                        }
                    }
                }
            }
        }
    }
    my %overlap_str_del_opt;
    if (scalar keys %del > 0){
        if ($chromium == 1){        # delete DELs whose BPs are flanked with a stretch of N (N6)
            my $rseq = $line[9];
            foreach my $dpos (keys %del){
                my $rpos = $del_rpos{$dpos};
                if ($cigar =~ /^(\d+)H/){
                    $rpos -= $1;
                }
                my $rseq_1 = '';
                my $rseq_2 = '';
                $rseq_1 = substr ($rseq, $rpos - 19, 18) if ($rpos >= 19);
                $rseq_2 = substr ($rseq, $rpos, 18) if (length ($rseq) - $rpos > 18);
                if (($rseq_1 =~ /N{6}/) or ($rseq_2 =~ /N{6}/)){
                    delete $del{$dpos};
                }
            }
        }
        foreach my $dpos (sort {$a <=> $b} keys %del){      # assign DELs (>= 20 bp) within TRs
            my $dlen = $del{$dpos};
            my $dend = $dpos + $dlen - 1;
            if ($exclude_overlap_flag > 0){
                my $ovl_flag = 0;
                foreach my $xpos (sort {$a <=> $b} keys %exclude_region){
                    last if ($xpos > $dend);
                    my $xend = $exclude_region{$xpos};
                    next if ($xend < $dpos);
                    my $overlap = 0;
                    if (($dpos >= $xpos) and ($dend <= $xend)){
                        $ovl_flag = 1;
                        last;
                    }
                    elsif (($xpos >= $dpos) and ($xpos <= $dend)){
                        $overlap = $dend - $xpos + 1;
                    }
                    elsif (($xend >= $dpos) and ($xend <= $dend)){
                        $overlap = $xend - $dpos + 1;
                    }
                    if ($overlap >= $dlen * 0.5){
                        $ovl_flag = 1;
                        last;
                    }
                }
                if ($ovl_flag == 1){
                    delete $del{$dpos};
                    next;
                }
            }
            my $Mbin_start = int ($dpos / $Mbin_size);
            if (exists $repeat{$Mbin_start}){
                my %match_str;
                foreach my $spos (sort {$a <=> $b} keys %{$repeat{$Mbin_start}}){
                    last if ($spos > $dend);
                    my $send = ${$repeat{$Mbin_start}}{$spos};
                    next if ($send < $dpos);
                    my $slen = $send - $spos + 1;
                    my $str = $STR{$spos};
                    my $overlap = 0;
                    if (($dpos <= $spos) and ($dend >= $send)){
                        $overlap = $slen;
                    }
                    elsif (($dpos >= $spos) and ($dpos <= $send)){
                        $overlap = $send - $dpos + 1;
                        $overlap = $dlen if ($dend < $send);
                    }
                    elsif (($dend >= $spos) and ($dend <= $send)){
                        $overlap = $dend - $spos + 1;
                        $overlap = $dlen if ($dpos > $spos);
                    }
                    if ($overlap >= $dlen * 0.1){
                        $match_str{$str} = $overlap;
                    }
                }
                if (scalar keys %match_str > 0){
                    my $match_num = scalar keys %match_str;
                    my $opt_str = '';
                    my $opt_str2 = '';
                    foreach my $str (sort {$match_str{$b} <=> $match_str{$a}} keys %match_str){
                        my $overlap = $match_str{$str};
                        if (($match_num == 1) and ($overlap >= $dlen * $min_overlap_rate)){
                            $STRdel{$dpos} = $str;
                            ${$STRdel_id{$str}}{$dpos} = $dlen;
                            $STR_id{$str} = 1;
                        }
                        elsif ($match_num > 1){
                            if (($opt_str eq '') and ($overlap >= $dlen * $min_overlap_rate)){
                                $STRdel{$dpos} = $str;
                                ${$STRdel_id{$str}}{$dpos} = $dlen;
                                $STR_id{$str} = 1;
                                $opt_str = $str;
                                next;
                            }
                            if (($opt_str ne '') and ($opt_str2 eq '')){
                                $opt_str2 = $str;
                            }
                        }
                        last;
                    }
                    if ($opt_str2 ne ''){
                        $overlap_str_del_opt{$opt_str} = $opt_str2;
                    }
                }
            }
        }
    }
    foreach my $sid1 (sort keys %overlap_str_opt){           # merge TR-INS overlapping with neighboring TRs
        next if (!exists $STRins_id{$sid1});
        my $sid2 = $overlap_str_opt{$sid1};
        my $ovl_inslen1 = 0;
        my $ovl_inslen2 = 0;
        next if (!exists $STRins_id{$sid2});
        my @ovl_ipos1;
        my @ovl_ipos2;
        my @ipos1;
        my @ipos2;
        my ($spos1, $send1) = split (/=/, $STR2{$sid1});
        my ($spos2, $send2) = split (/=/, $STR2{$sid2});
        foreach my $ipos (sort {$a <=> $b} keys %{$STRins_id{$sid1}}){
            my $ilen = ${$STRins_id{$sid1}}{$ipos};
            next if ($ilen < $max_dist2);
            if (($ipos >= $spos2) and ($ipos <= $send2)){
                $ovl_inslen1 += $ilen;
                push @ovl_ipos1, $ipos;
            }
            else{
                push @ipos1, $ipos;
            }
        }
        foreach my $ipos (sort {$a <=> $b} keys %{$STRins_id{$sid2}}){
            my $ilen = ${$STRins_id{$sid2}}{$ipos};
            next if ($ilen < $max_dist2);
            if (($ipos >= $spos1) and ($ipos <= $send1)){
                $ovl_inslen2 += $ilen;
                push @ovl_ipos2, $ipos;
            }
            else{
                push @ipos2, $ipos;
            }
        }
        next if ($ovl_inslen1 == 0) or ($ovl_inslen2 == 0);
        if ((@ipos1 == 0) and (@ipos2 > 0)){
            delete $STRins_id{$sid1};
        }
        elsif ((@ipos1 > 0) and (@ipos2 == 0)){
            delete $STRins_id{$sid2};
        }
        elsif ((@ipos1 == 0) and (@ipos2 == 0)){
            my $str_overlap = 0;
            if ($spos1 > $spos2){
                $str_overlap = $send1 - $spos2 + 1;
            }
            else{
                $str_overlap = $send2 - $spos1 + 1;
            }
            my $ovl_rate1 = int ($str_overlap / ($send1 - $spos1 + 1) * 100) / 100;
            my $ovl_rate2 = int ($str_overlap / ($send2 - $spos2 + 1) * 100) / 100;
            if ($ovl_rate1 <= $ovl_rate2){
                delete $STRins_id{$sid2};
            }
            else{
                delete $STRins_id{$sid1};
            }
        }
        else{
            my $sum_ovl_pos1 = 0;
            my $sum_ovl_pos2 = 0;
            my $sum_pos1 = 0;
            my $sum_pos2 = 0;
            map{$sum_ovl_pos1 += $_} @ovl_ipos1;
            map{$sum_ovl_pos2 += $_} @ovl_ipos2;
            map{$sum_pos1 += $_} @ipos1;
            map{$sum_pos2 += $_} @ipos2;
            my $ave_ovl_pos1 = int ($sum_ovl_pos1 / @ovl_ipos1);
            my $ave_ovl_pos2 = int ($sum_ovl_pos2 / @ovl_ipos2);
            my $ave_pos1 = int ($sum_pos1 / @ipos1);
            my $ave_pos2 = int ($sum_pos2 / @ipos2);
            my $diff1 = abs ($ave_ovl_pos1 - $ave_pos1);
            my $diff2 = abs ($ave_ovl_pos2 - $ave_pos2);
            if ($diff1 <= $diff2){
                foreach my $ipos (@ipos2){
                    delete ${$STRins_id{$sid2}}{$ipos}
                }
            }
            else{
                foreach my $ipos (@ipos1){
                    delete ${$STRins_id{$sid1}}{$ipos}
                }
            }
        }
    }
    foreach my $sid1 (sort keys %overlap_str_del_opt){       # merge TR-DEL overlapping with neighboring TRs
        next if (exists $STRdel_id{$sid1});
        my $sid2 = $overlap_str_del_opt{$sid1};
        my $dellen1 = 0;
        my $dellen2 = 0;
        next if (!exists $STRdel_id{$sid2});
        my ($spos1, $send1) = split (/=/, $STR2{$sid1});
        my ($spos2, $send2) = split (/=/, $STR2{$sid2});
        my $slen1 = $send1 - $spos1 + 1;
        my $slen2 = $send2 - $spos2 + 1;
        foreach my $dpos (sort {$a <=> $b} keys %{$STRdel_id{$sid1}}){
            my $dlen = ${$STRdel_id{$sid1}}{$dpos};
            next if ($dlen < $max_dist2);
            my $dend = $dpos + $dlen - 1;
            my $overlap1 = 0;
            my $overlap2 = 0;
            if (($dpos >= $spos2) and ($dpos <= $send2)){
                $overlap2 = $send2 - $dpos + 1;
                $overlap2 = $dlen if ($dend < $send2);
            }
            elsif (($dend >= $spos2) and ($dend <= $send2)){
                $overlap2 = $dend - $spos2 + 1;
                $overlap2 = $dlen if ($dpos > $spos2);
            }
            elsif (($dpos < $spos2) and ($dend > $send2)){
                $overlap2 = $slen2;
            }
            next if ($overlap2 == 0);
            if (($dpos >= $spos1) and ($dpos <= $send1)){
                $overlap1 = $send1 - $dpos + 1;
                $overlap1 = $dlen if ($dend < $send1);
            }
            elsif (($dend >= $spos1) and ($dend <= $send1)){
                $overlap1 = $dend - $spos1 + 1;
                $overlap1 = $dlen if ($dpos > $spos1);
            }
            elsif (($dpos < $spos1) and ($dend > $send1)){
                $overlap1 = $slen1;
            }
            if ($overlap1 >= $overlap2){
                delete ${$STRdel_id{$sid2}}{$dpos} if (exists ${$STRdel_id{$sid2}}{$dpos});
            }
            else{
                delete ${$STRdel_id{$sid1}}{$dpos};
                if (!exists ${$STRdel_id{$sid2}}{$dpos}){
                    ${$STRdel_id{$sid2}}{$dpos} = $dlen;
                }
            }
        }
        if (scalar keys %{$STRdel_id{$sid1}} == 0){
            delete $STRdel_id{$sid1};
        }
        if (scalar keys %{$STRdel_id{$sid2}} == 0){
            delete $STRdel_id{$sid2};
        }
    }
    foreach my $sid (sort keys %STR_id){                         # Take either INSs or DELs, that were mapped in a read aligned within a TR, and estimate a merged size of fragmented INSs mapped within a read within a TR
        next if (!exists $STRins_id{$sid}) and (!exists $STRdel_id{$sid});
        my $inslen = 0;
        my $dellen = 0;
        my $max_ipos = 0;
        my $max_ilen = 0;
        my ($spos, $send, $motif_size) = split (/=/, $STR2{$sid});
        my $slen = $send - $spos + 1;
        my @ipos;
        foreach my $ipos (sort {$a <=> $b} keys %{$STRins_id{$sid}}){
            my $ilen = ${$STRins_id{$sid}}{$ipos};
            next if ($ilen < $max_dist2);
            if ($ilen > $max_ilen){
                $max_ilen = $ilen;
                $max_ipos = $ipos;
            }
            $inslen += $ilen;
            push @ipos, $ipos if ($ilen >= $min_str_indel_size);
        }
        foreach my $dpos (sort {$a <=> $b} keys %{$STRdel_id{$sid}}){
            my $dlen = ${$STRdel_id{$sid}}{$dpos};
            next if ($dlen < $max_dist2);
            $dellen += $dlen;
        }
        if ($inslen >= $dellen){
            my $ilen = $inslen - $dellen;
            if (($ilen >= $min_str_indel_size) and (!exists $str_bp_flag{$sid})){
                push @{${$bam_str{$sid}}{'INS'}}, $ilen;
                if (($slen > 1000) and (exists $str_bp_flag2{$sid})){
                    push @{${$bam_str2{$sid}}{'INS'}}, $ilen;
                }
                if ((!exists $BP_str{$sid}) and ($mapQ > 0) and (exists $ipos_str{$sid})){
                    ${$bam_str_insinfo{$sid}}{$read_id} = "$ilen=$ipos_str{$sid}=$max_ipos";
                }
                elsif (($mapQ > 0) and (exists $ipos_str{$sid})){
                    ${$bam_str_insinfo2{$sid}}{$read_id} = "$ilen=$ipos_str{$sid}=$max_ipos";
                }
                if ((exists $STRovl{$sid}) and (@ipos > 0)){
                    my $sum_ipos = 0;
                    map{$sum_ipos += $_} @ipos;
                    my $ave = int ($sum_ipos / @ipos + 0.5);
                    push @{${$bam_str_pos{$sid}}{$ilen}}, $ave;
                }
                delete $STRdel_id{$sid};
            }
            else{
                delete $STRins_id{$sid};
                delete $STRdel_id{$sid};
            }
        }
        else{
            my $dlen = $dellen - $inslen;
            if (($dlen >= $min_str_indel_size) and (!exists $str_bp_flag{$sid})){
                push @{${$bam_str{$sid}}{'DEL'}}, $dlen;
                if (($slen > 1000) and (exists $str_bp_flag2{$sid})){
                    push @{${$bam_str2{$sid}}{'DEL'}}, $dlen;
                }
                delete $STRins_id{$sid};
            }
            else{
                delete $STRins_id{$sid};
                delete $STRdel_id{$sid};
            }
        }
    }

    if (scalar keys %ins > 0){
        my $pre_pos = 0;
        my $pre_len = 0;
        my %merge;
        foreach my $ipos (sort {$a <=> $b} keys %ins){
            next if (exists $STRins{$ipos});
            my ($id1, $strand1, $ilen, $rpos) = split (/=/, $ins{$ipos});
            if ($pre_pos == 0){
                $pre_len = $ilen;
                $pre_pos = $ipos;
                $merge{$ipos} = $ilen;
                next;
            }
            my $distance = $ipos - $pre_pos + 1;
            my $flag = 0;
            if (($distance <= 50) and ($distance < $ilen) and ($distance < $pre_len)){
                $flag = 1;
            }
            elsif (($distance <= 200) and ($distance < $ilen * 0.2) and ($distance < $pre_len * 0.2)){
                $flag = 1;
            }
            elsif (($pre_pos > 0) and ($distance <= 200) and ($ilen >= 200) and ($pre_len >= 200) and (($distance < $pre_len * 0.5) or ($distance < $ilen * 0.5))){
                $flag = 1;
            }
            elsif (($pre_pos > 0) and ($distance <= 200) and (($ilen / $pre_len >= 10) or ($pre_len / $ilen >= 10)) and (($distance < $pre_len * 0.2) or ($distance < $ilen * 0.2))){
                $flag = 1;
            }
            if ($flag == 1){
                $merge{$ipos} = $ilen;
                $pre_len += $ilen;
                $pre_pos = $ipos;
                next;
            }
            else{
                if (scalar keys %merge > 1){
                    my $max_pos = 0;
                    my $max_len = 0;
                    my $sum_len = 0;
                    my $first_pos = 0;
                    my $last_pos = 0;
                    foreach my $ipos2 (sort {$a <=> $b} keys %merge){
                        my ($id2, $strand2, $ilen2, $rpos2) = split (/=/, $ins{$ipos2});
                        if ($ilen2 > $max_len){
                            $max_len = $ilen2;
                            $max_pos = $ipos2;
                        }
                        if ($first_pos == 0){
                            $first_pos = $ipos2;
                        }
                        $last_pos = $ipos2;
                    }
                    my ($id2, $strand2, $ilen2, $rpos2) = split (/=/, $ins{$max_pos});
                    foreach my $ipos2 (sort {$a <=> $b} keys %ins){
                        last if ($ipos2 > $last_pos);
                        next if ($ipos2 < $first_pos);
                        my ($id2, $strand2, $ilen2, $rpos2) = split (/=/, $ins{$ipos2});
                        $sum_len += $ilen2;
                        delete $ins{$ipos2};
                    }
                    $ins{$max_pos} = "$id2=$strand2=$sum_len=$rpos2";
                }
                %merge = ();
                $merge{$ipos} = $ilen;
            }
            $pre_len = $ilen;
            $pre_pos = $ipos;
        }
        if (scalar keys %merge > 1){
            my $max_pos = 0;
            my $max_len = 0;
            my $sum_len = 0;
            my $first_pos = 0;
            my $last_pos = 0;
            foreach my $ipos2 (sort {$a <=> $b} keys %merge){
                my ($id2, $strand2, $ilen2, $rpos2) = split (/=/, $ins{$ipos2});
                if ($ilen2 > $max_len){
                    $max_len = $ilen2;
                    $max_pos = $ipos2;
                }
                if ($first_pos == 0){
                    $first_pos = $ipos2;
                }
                $last_pos = $ipos2;
            }
            my ($id2, $strand2, $ilen2, $rpos2) = split (/=/, $ins{$max_pos});
            foreach my $ipos2 (sort {$a <=> $b} keys %ins){
                last if ($ipos2 > $last_pos);
                next if ($ipos2 < $first_pos);
                my ($id2, $strand2, $ilen2, $rpos2) = split (/=/, $ins{$ipos2});
                $sum_len += $ilen2;
                delete $ins{$ipos2};
            }
            $ins{$max_pos} = "$id2=$strand2=$sum_len=$rpos2";
        }
        foreach my $ipos (keys %ins){
            next if (exists $STRins{$ipos});
            my ($id, $istrand, $ilen, $rpos) = split (/=/, $ins{$ipos});
            if ($ilen >= $max_dist){
                if (exists $bam_ins{$ipos}){
                    my $sumlen = 0;
                    foreach (@{$bam_ins{$ipos}}){
                        my ($id1, $strand1, $ilen1) = split (/=/, $_);
                        $sumlen += $ilen1;
                    }
                    my $ave = int ($sumlen / @{$bam_ins{$ipos}});
                    if (($ave / $ilen <= 1.7) and ($ave / $ilen >= 0.59)){
                        push @{$bam_ins{$ipos}}, "$ins{$ipos}=$ipos";
                    }
                    else{
                        my $ipos2 = $ipos;
                        while (1){
                            $ipos2 ++;
                            if (exists $bam_ins{$ipos2}){
                                my ($id2, $strand2, $ilen2) = split (/=/, ${$bam_ins{$ipos2}}[0]);
                                if (($ilen2 / $ilen <= 1.7) and ($ilen2 / $ilen >= 0.59)){
                                    push @{$bam_ins{$ipos2}}, "$ins{$ipos}=$ipos2";
                                    last;
                                }
                            }
                            else{
                                push @{$bam_ins{$ipos2}}, "$ins{$ipos}=$ipos2";
                                last;
                            }
                        }
                    }
                }
                else{
                    push @{$bam_ins{$ipos}}, "$ins{$ipos}=$ipos";
                }
            }
        }
    }
    if (scalar keys %del > 0){
        my %removed;
        foreach my $dpos1 (sort {$a <=> $b} keys %del){
            next if (exists $removed{$dpos1});
            next if (exists $STRdel{$dpos1});
            my $dlen1 = $del{$dpos1};
            my $dend1 = $dpos1 + $dlen1 - 1;
            while (1){
                my $flag = 0;
                my @dpos2 = ();
                foreach my $dpos2 (sort {$a <=> $b} keys %del){
                    next if ($dpos2 <= $dpos1);
                    next if (exists $removed{$dpos2});
                    next if (exists $STRdel{$dpos2});
                    last if ($dpos2 > $dend1 + 500);
                    my $dlen2 = $del{$dpos2};
                    my $dend2 = $dpos2 + $dlen2 - 1;
                    my $distance = $dpos2 - $dend1 + 1;
                    if (($distance < $dlen1 * 0.5) and ($distance < $dlen2 * 0.5) and ($distance < 500)){
                        $flag = $dpos2;
                    }
                    push @dpos2, $dpos2;
                }
                if ($flag > 0){
                    my $dlen2 = $del{$flag};
                    my $dend2 = $flag + $dlen2 - 1;
                    foreach my $dpos3 (@dpos2){
                        last if ($dpos3 > $flag);
                        next if (!exists $del{$dpos3});
                        delete $del{$dpos3};
                        $removed{$dpos3} = 1;
                    }
                    my $new_len = $dend2 - $dpos1 + 1;
                    $del{$dpos1} = $new_len;
                    $dlen1 = $new_len;
                    $dend1 = $dpos1 + $new_len - 1;
                }
                else{
                    last;
                }
            }
        }
        foreach my $dpos (sort {$a <=> $b} keys %del){
            next if (exists $STRdel{$dpos});
            my $dlen = $del{$dpos};
            if ($dlen >= $max_dist){
                if (exists $bam_del{$dpos}){
                    my $sumlen = 0;
                    foreach (@{$bam_del{$dpos}}){
                        my ($dlen2) = split (/=/, $_);
                        $sumlen += $dlen2;
                    }
                    my $ave = int ($sumlen / @{$bam_del{$dpos}});
                    if (($ave / $dlen <= 1.7) and ($ave / $dlen >= 0.59)){
                        push @{$bam_del{$dpos}}, "$dlen=$dpos=$read_id";
                    }
                    else{
                        my $dpos2 = $dpos;
                        while (1){
                            $dpos2 ++;
                            if (exists $bam_del{$dpos2}){
                                my ($dlen3) = split (/=/, ${$bam_del{$dpos2}}[0]);
                                if (($dlen3 / $dlen <= 1.7) and ($dlen3 / $dlen >= 0.59)){
                                    push @{$bam_del{$dpos2}}, "$dlen=$dpos2=$read_id";
                                    last;
                                }
                            }
                            else{
                                push @{$bam_del{$dpos2}}, "$dlen=$dpos2=$read_id";
                                last;
                            }
                        }
                    }
                }
                else{
                    push @{$bam_del{$dpos}}, "$dlen=$dpos=$read_id";
                }
            }
        }
    }
    if (scalar keys %rep > 0){
        my $pre_pos = 0;
        my $pre_end = 0;
        my $pre_len = 0;
        foreach my $dpos (sort {$a <=> $b} keys %rep){
            my $dlen = $rep{$dpos};
            my $dend = $dpos + $dlen - 1;
            my $distance = $dpos - $dend + 1;
            if (($pre_pos > 0) and ($distance < $pre_len) and ($distance < $dlen) and ($distance < 500)){
                my $new_len = $dend - $pre_pos + 1;
                $rep{$pre_pos} = $new_len;
                delete $rep{$dpos};
                $pre_len = $new_len;
                $pre_end = $dend;
                next;
            }
            $pre_pos = $dpos;
            $pre_len = $dlen;
            $pre_end = $dend;
        }
        foreach my $rpos (keys %rep){
            my $xlen = $rep{$rpos};
            if ($xlen >= $min_indel_size){
                push @{$bam_rep{$rpos}}, $xlen;
            }
        }
    }
}
close (FILE);

my $ave_depth = int ($sum_align / $chr_len{$chr} * 10 + 0.5) / 10;
$ave_depth = 20 if ($ave_depth < 20);
my $max_depth = $ave_depth * $max_depth_fold;

my %highDP;
my %high_depth_range;
my $window = 1000;
foreach my $mbin (sort {$a <=> $b} keys %DP){
    foreach my $res (sort {$a <=> $b} keys %{$DP{$mbin}}){
        next if ($res % 100 > 0);
        next if (${$DP{$mbin}}{$res} <= $max_depth);
        my $dp = ${$DP{$mbin}}{$res};
        my $pos = $Mbin_size * $mbin + $res;
        $highDP{$pos} = $dp;
    }
}

if ($targeted_seq == 0){
    print STDERR "$chr: Mean depth: $ave_depth\n";
    print STDERR "$chr: Max threshold depth: $max_depth ($ave_depth ", 'x ', "$max_depth_fold)\n";
    print STDERR "$chr: High-coverage regions\n";

    if (scalar keys %highDP > 0){
        my $pre_dppos = 0;
        my @dp;
        my $start = 0;
        my $sum_drange = 0;
        print STDERR "#CHR\tSTART\tEND\tDEPTH\n";
        foreach my $pos (sort {$a <=> $b} keys %highDP){
            my $dp = $highDP{$pos};
            if ($pre_dppos == 0){
                $pre_dppos = $pos;
                push @dp, $dp;
                $start = $pos;
                next;
            }
            if (($pre_dppos > 0) and ($pos - $pre_dppos <= $window)){
                push @dp, $dp;
            }
            else{
                my $end = $pre_dppos + 100;
                $high_depth_range{$start} = $end;
                my $sumdp = 0;
                map{$sumdp += $_} @dp;
                my $avedp = int ($sumdp / @dp + 0.5);
                print STDERR "$chr\t$start\t$end\t$avedp\n";
                $sum_drange += $end - $start;
                $start = $pos;
                @dp = ();
                push @dp, $dp;
            }
            $pre_dppos = $pos;
        }
        if ($start != $pre_dppos){
            $high_depth_range{$start} = $pre_dppos;
            my $sumdp = 0;
            map{$sumdp += $_} @dp;
            my $avedp = int ($sumdp / @dp + 0.5);
            print STDERR "$chr\t$start\t$pre_dppos\t$avedp\n";
            $sum_drange += $pre_dppos - $start;
        }
        %highDP = ();
        print STDERR "$chr: Total length of high-coverage regions (bp): $sum_drange\n";
    }

    if (scalar keys %high_depth_range > 0){
        foreach my $dpos (sort {$a <=> $b} keys %high_depth_range){
            my $dend = $high_depth_range{$dpos};
            foreach my $pos (sort {$a <=> $b} keys %bam_ins){
                last if ($pos > $dend);
                next if ($pos < $dpos);
                delete $bam_ins{$pos};
            }
            foreach my $pos (sort {$a <=> $b} keys %bam_del){
                last if ($pos > $dend - 50);
                next if ($pos < $dpos);
                delete $bam_del{$pos};
            }
            foreach my $pos (sort {$a <=> $b} keys %bam_bp1){
                last if ($pos > $dend);
                next if ($pos < $dpos);
                delete $bam_bp1{$pos};
            }
            foreach my $pos (sort {$a <=> $b} keys %bam_bp2){
                last if ($pos > $dend);
                next if ($pos < $dpos);
                delete $bam_bp2{$pos};
            }
            foreach my $readid (keys %read_bp1){
                foreach my $pos (sort {$a <=> $b} keys %{$read_bp1{$readid}}){
                    last if ($pos > $dend);
                    next if ($pos < $dpos);
                    delete ${$read_bp1{$readid}}{$pos};
                    delete $read_bp1{$readid} if (scalar keys %{$read_bp1{$readid}} == 0);
                }
            }
            foreach my $readid (keys %read_bp2){
                foreach my $pos (sort {$a <=> $b} keys %{$read_bp2{$readid}}){
                    last if ($pos > $dend);
                    next if ($pos < $dpos);
                    delete ${$read_bp2{$readid}}{$pos};
                    delete $read_bp2{$readid} if (scalar keys %{$read_bp2{$readid}} == 0);
                }
            }
        }
    }
}

my %str_used;
foreach my $sid1 (sort keys %bam_str){        # check and reaasign TR-INS between neighboring TRs with <= 20 bp distance
    if ((exists ${$bam_str{$sid1}}{'INS'}) and (exists $STRovl{$sid1})){
        my @sid2 = split (/,/, $STRovl{$sid1});
        foreach my $sid2 (@sid2){
            next if (!exists ${$bam_str{$sid2}}{'INS'});
            next if (exists ${$str_used{$sid1}}{$sid2});
            ${$str_used{$sid2}}{$sid1} = 1;
#            next if (!exists $bam_str_ins_ovl{$sid1}) and (!exists $bam_str_ins_ovl{$sid2});
            my $max1 = 0;
            my $max2 = 0;
            my $min1 = 1000000000;
            my $min2 = 1000000000;
            my @max1;
            my @max2;
            my @min1;
            my @min2;
            my $sid1_len1 = 0;
            my $sid1_len2 = 0;
            my $sid2_len1 = 0;
            my $sid2_len2 = 0;
            my $num1 = scalar @{${$bam_str{$sid1}}{'INS'}};
            my $num2 = scalar @{${$bam_str{$sid2}}{'INS'}};
            foreach (@{${$bam_str{$sid1}}{'INS'}}){
                if ($_ > $max1){
                    $max1 = $_;
                }
                if ($_ < $min1){
                    $min1 = $_;
                }
            }
            foreach (@{${$bam_str{$sid2}}{'INS'}}){
                if ($_ > $max2){
                    $max2 = $_;
                }
                if ($_ < $min2){
                    $min2 = $_;
                }
            }
            if ($max1 / $min1 <= 1.5){
                @max1 = sort {$a <=> $b} @{${$bam_str{$sid1}}{'INS'}};
                my $hn = int (@max1 * 0.5);
                $sid1_len1 = $max1[$hn];
            }
            else{
                foreach (@{${$bam_str{$sid1}}{'INS'}}){
                    my $diff1 = $max1 - $_;
                    my $diff2 = $_ - $min1;
                    if ($diff1 <= $diff2){
                        push @max1, $_;
                    }
                    else{
                        push @min1, $_;
                    }
                }
                my $hn1 = int (@max1 * 0.5);
                @max1 = sort {$a <=> $b} @max1;
                $sid1_len1 = $max1[$hn1];
                my $hn2 = int (@min1 * 0.5);
                @min1 = sort {$a <=> $b} @min1;
                $sid1_len2 = $min1[$hn2];
                if (($sid1_len1 / $sid1_len2 <= 1.5) and ($sid1_len1 / $sid1_len2 >= 0.67)){
                    push @max1, @min1;
                    @min1 = ();
                }
            }
            if ($max2 / $min2 <= 1.5){
                @max2 = sort {$a <=> $b} @{${$bam_str{$sid2}}{'INS'}};
                my $hn = int (@max2 * 0.5);
                $sid2_len1 = $max2[$hn];
            }
            else{
                foreach (@{${$bam_str{$sid2}}{'INS'}}){
                    my $diff1 = $max2 - $_;
                    my $diff2 = $_ - $min2;
                    if ($diff1 <= $diff2){
                        push @max2, $_;
                    }
                    else{
                        push @min2, $_;
                    }
                }
                my $hn1 = int (@max2 * 0.5);
                @max2 = sort {$a <=> $b} @max2;
                $sid2_len1 = $max2[$hn1];
                my $hn2 = int (@min2 * 0.5);
                @min2 = sort {$a <=> $b} @min2;
                $sid2_len2 = $min2[$hn2];
                if (($sid2_len1 / $sid2_len2 <= 1.5) and ($sid2_len1 / $sid2_len2 >= 0.67)){
                    push @max2, @min2;
                    @min2 = ();
                }
            }
            my $max1_pos = 0;
            my $max2_pos = 0;
            my $min1_pos = 0;
            my $min2_pos = 0;
            my @max1_pos;
            my @max2_pos;
            my @min1_pos;
            my @min2_pos;
            if (@max1 > 0){
                my $sum1 = 0;
                foreach my $ilen (@max1){
                    next if (!exists ${$bam_str_pos{$sid1}}{$ilen});
                    push @max1_pos, @{${$bam_str_pos{$sid1}}{$ilen}};
                }
                if (@max1_pos > 0){
                    map{$sum1 += $_} @max1_pos;
                    $max1_pos = int ($sum1 / @max1_pos + 0.5);
                }
            }
            if (@max2 > 0){
                my $sum2 = 0;
                foreach my $ilen (@max2){
                    next if (!exists ${$bam_str_pos{$sid2}}{$ilen});
                    push @max2_pos, @{${$bam_str_pos{$sid2}}{$ilen}};
                }
                if (@max2_pos > 0){
                    map{$sum2 += $_} @max2_pos;
                    $max2_pos = int ($sum2 / @max2_pos + 0.5);
                }
            }
            if (@min1 > 0){
                my $sum1 = 0;
                foreach my $ilen (@min1){
                    next if (!exists ${$bam_str_pos{$sid1}}{$ilen});
                    push @min1_pos, @{${$bam_str_pos{$sid1}}{$ilen}};
                }
                if (@min1_pos > 0){
                    map{$sum1 += $_} @min1_pos;
                    $min1_pos = int ($sum1 / @min1_pos + 0.5);
                }
            }
            if (@min2 > 0){
                my $sum2 = 0;
                foreach my $ilen (@min2){
                    next if (!exists ${$bam_str_pos{$sid2}}{$ilen});
                    push @min2_pos, @{${$bam_str_pos{$sid2}}{$ilen}};
                }
                if (@min2_pos > 0){
                    map{$sum2 += $_} @min2_pos;
                    $min2_pos = int ($sum2 / @min2_pos + 0.5);
                }
            } 
            if (($sid1_len1 / $sid2_len1 <= 1.4) and ($sid1_len1 / $sid2_len1 >= 0.71)){
                my $flag = 0;
                if (abs ($max1_pos - $max2_pos) <= 150){
                    $flag = 1;
                }
                elsif (((abs ($max1_pos - $max2_pos) <= 300) or ($max1_pos == 0) or ($max2_pos == 0)) and ($sid1_len1 / $sid2_len1 <= 1.25) and ($sid1_len1 / $sid2_len1 >= 0.8)){
                    $flag = 1;
                }
                if ($flag == 1){
                    if ($num1 >= $num2){
                        push @max1, @max2;
                        @max2 = ();
                    }
                    else{
                        push @max2, @max1;
                        @max1 = ();
                    }
                }
            }
            elsif ((@min2 > 0) and ($sid1_len1 / $sid2_len2 <= 1.4) and ($sid1_len1 / $sid2_len2 >= 0.71)){
                my $flag = 0;
                if (abs ($max1_pos - $min2_pos) <= 150){
                    $flag = 1;
                }
                elsif (((abs ($max1_pos - $min2_pos) <= 300) or ($max1_pos == 0) or ($min2_pos == 0)) and ($sid1_len1 / $sid2_len2 <= 1.25) and ($sid1_len1 / $sid2_len2 >= 0.8)){
                    $flag = 1;
                }
                if ($flag == 1){
                    if ($num1 >= $num2){
                        push @max1, @min2;
                        @min2 = ();
                    }
                    else{
                        push @min2, @max1;
                        @max1 = ();
                    }
                }
            }
            elsif ((@min1 > 0) and ($sid2_len1 / $sid1_len2 <= 1.4) and ($sid2_len1 / $sid1_len2 >= 0.71)){
                my $flag = 0;
                if (abs ($max2_pos - $min1_pos) <= 150){
                    $flag = 1;
                }
                elsif (((abs ($max2_pos - $min1_pos) <= 300) or ($max2_pos == 0) or ($min1_pos == 0)) and ($sid2_len1 / $sid1_len2 <= 1.25) and ($sid2_len1 / $sid1_len2 >= 0.8)){
                    $flag = 1;
                }
                if ($flag == 1){
                    if ($num1 >= $num2){
                        push @min1, @max2;
                        @max2 = ();
                    }
                    else{
                        push @max2, @min1;
                        @min1 = ();
                    }
                }
            }
            if ((@min1 > 0) and (@min2 > 0) and ($sid1_len2 / $sid2_len2 <= 1.4) and ($sid1_len2 / $sid2_len2 >= 0.71)){
                my $flag = 0;
                if (abs ($min1_pos - $min2_pos) <= 150){
                    $flag = 1;
                }
                elsif (((abs ($min1_pos - $min2_pos) <= 300) or ($min1_pos == 0) or ($min2_pos == 0)) and ($sid1_len2 / $sid2_len2 <= 1.25) and ($sid1_len2 / $sid2_len2 >= 0.8)){
                    $flag = 1;
                }
                if ($flag == 1){
                    if ($num1 >= $num2){
                        push @min1, @min2;
                        @min2 = ();
                    }
                    else{
                        push @min2, @min1;
                        @min1 = ();
                    }
                }
            }
            my $new_num1 = @max1 + @min1;
            my $new_num2 = @max2 + @min2;
            if ($num1 != $new_num1){
                if ($new_num2 > 0){
                    @{${$bam_str{$sid2}}{'INS'}} = (@max2, @min2);
                }
                else{
                    delete ${$bam_str{$sid2}}{'INS'};
                    delete $bam_str{$sid2} if (scalar keys %{$bam_str{$sid2}} == 0);
                    if (exists ${$bam_str2{$sid2}}{'INS'}){
                        delete ${$bam_str2{$sid2}}{'INS'};
                        delete $bam_str2{$sid2} if (scalar keys %{$bam_str{$sid2}} == 0);
                    }
                }
                if ($new_num1 > 0){
                    @{${$bam_str{$sid1}}{'INS'}} = (@max1, @min1);
                }
                else{
                    delete ${$bam_str{$sid1}}{'INS'};
                    delete $bam_str{$sid1} if (scalar keys %{$bam_str{$sid1}} == 0);
                    if (exists ${$bam_str2{$sid1}}{'INS'}){
                        delete ${$bam_str2{$sid1}}{'INS'};
                        delete $bam_str2{$sid1} if (scalar keys %{$bam_str2{$sid1}} == 0);
                    }
                    last;
                }
            }
        }
    }
}
%str_used = ();
foreach my $sid1 (sort keys %bam_str2){        # check and reaasign TR-INS between neighboring TRs with <= 20 bp distance
    if ((exists ${$bam_str2{$sid1}}{'INS'}) and (exists $STRovl{$sid1})){
        my @sid2 = split (/,/, $STRovl{$sid1});
        foreach my $sid2 (@sid2){
            next if (!exists ${$bam_str2{$sid2}}{'INS'});
            next if (exists ${$str_used{$sid1}}{$sid2});
            ${$str_used{$sid2}}{$sid1} = 1;
#            next if (!exists $bam_str_ins_ovl2{$sid1}) and (!exists $bam_str_ins_ovl2{$sid2});
            my $max1 = 0;
            my $max2 = 0;
            my $min1 = 1000000000;
            my $min2 = 1000000000;
            my @max1;
            my @max2;
            my @min1;
            my @min2;
            my $sid1_len1 = 0;
            my $sid1_len2 = 0;
            my $sid2_len1 = 0;
            my $sid2_len2 = 0;
            my $num1 = scalar @{${$bam_str2{$sid1}}{'INS'}};
            my $num2 = scalar @{${$bam_str2{$sid2}}{'INS'}};
            foreach (@{${$bam_str2{$sid1}}{'INS'}}){
                if ($_ > $max1){
                    $max1 = $_;
                }
                if ($_ < $min1){
                    $min1 = $_;
                }
            }
            foreach (@{${$bam_str2{$sid2}}{'INS'}}){
                if ($_ > $max2){
                    $max2 = $_;
                }
                if ($_ < $min2){
                    $min2 = $_;
                }
            }
            if ($max1 / $min1 <= 1.5){
                @max1 = sort {$a <=> $b} @{${$bam_str2{$sid1}}{'INS'}};
                my $hn = int (@max1 * 0.5);
                $sid1_len1 = $max1[$hn];
            }
            else{
                foreach (@{${$bam_str2{$sid1}}{'INS'}}){
                    my $diff1 = $max1 - $_;
                    my $diff2 = $_ - $min1;
                    if ($diff1 <= $diff2){
                        push @max1, $_;
                    }
                    else{
                        push @min1, $_;
                    }
                }
                my $hn1 = int (@max1 * 0.5);
                @max1 = sort {$a <=> $b} @max1;
                $sid1_len1 = $max1[$hn1];
                my $hn2 = int (@min1 * 0.5);
                @min1 = sort {$a <=> $b} @min1;
                $sid1_len2 = $min1[$hn2];
                if (($sid1_len1 / $sid1_len2 <= 1.5) and ($sid1_len1 / $sid1_len2 >= 0.67)){
                    push @max1, @min1;
                    @min1 = ();
                }
            }
            if ($max2 / $min2 <= 1.5){
                @max2 = sort {$a <=> $b} @{${$bam_str2{$sid2}}{'INS'}};
                my $hn = int (@max2 * 0.5);
                $sid2_len1 = $max2[$hn];
            }
            else{
                foreach (@{${$bam_str2{$sid2}}{'INS'}}){
                    my $diff1 = $max2 - $_;
                    my $diff2 = $_ - $min2;
                    if ($diff1 <= $diff2){
                        push @max2, $_;
                    }
                    else{
                        push @min2, $_;
                    }
                }
                my $hn1 = int (@max2 * 0.5);
                @max2 = sort {$a <=> $b} @max2;
                $sid2_len1 = $max2[$hn1];
                my $hn2 = int (@min2 * 0.5);
                @min2 = sort {$a <=> $b} @min2;
                $sid2_len2 = $min2[$hn2];
                if (($sid2_len1 / $sid2_len2 <= 1.5) and ($sid2_len1 / $sid2_len2 >= 0.67)){
                    push @max2, @min2;
                    @min2 = ();
                }
            }
            my $max1_pos = 0;
            my $max2_pos = 0;
            my $min1_pos = 0;
            my $min2_pos = 0;
            if (@max1 > 0){
                my @max1_pos;
                my $sum1 = 0;
                foreach my $ilen (@max1){
                    next if (!exists ${$bam_str_pos{$sid1}}{$ilen});
                    push @max1_pos, @{${$bam_str_pos{$sid1}}{$ilen}};
                }
                if (@max1_pos > 0){
                    map{$sum1 += $_} @max1_pos;
                    $max1_pos = int ($sum1 / @max1_pos + 0.5);
                }
            }
            if (@max2 > 0){
                my @max2_pos;
                my $sum2 = 0;
                foreach my $ilen (@max2){
                    next if (!exists ${$bam_str_pos{$sid2}}{$ilen});
                    push @max2_pos, @{${$bam_str_pos{$sid2}}{$ilen}};
                }
                if (@max2_pos > 0){
                    map{$sum2 += $_} @max2_pos;
                    $max2_pos = int ($sum2 / @max2_pos + 0.5);
                }
            }
            if (@min1 > 0){
                my @min1_pos;
                my $sum1 = 0;
                foreach my $ilen (@min1){
                    next if (!exists ${$bam_str_pos{$sid1}}{$ilen});
                    push @min1_pos, @{${$bam_str_pos{$sid1}}{$ilen}};
                }
                if (@min1_pos > 0){
                    map{$sum1 += $_} @min1_pos;
                    $min1_pos = int ($sum1 / @min1_pos + 0.5);
                }
            }
            if (@min2 > 0){
                my @min2_pos;
                my $sum2 = 0;
                foreach my $ilen (@min2){
                    next if (!exists ${$bam_str_pos{$sid2}}{$ilen});
                    push @min2_pos, @{${$bam_str_pos{$sid2}}{$ilen}};
                }
                if (@min2_pos > 0){
                    map{$sum2 += $_} @min2_pos;
                    $min2_pos = int ($sum2 / @min2_pos + 0.5);
                }
            }
            if (($sid1_len1 / $sid2_len1 <= 1.4) and ($sid1_len1 / $sid2_len1 >= 0.71)){
                my $flag = 0;
                if (abs ($max1_pos - $max2_pos) <= 150){
                    $flag = 1;
                }
                elsif (((abs ($max1_pos - $max2_pos) <= 300) or ($max1_pos == 0) or ($max2_pos == 0)) and ($sid1_len1 / $sid2_len1 <= 1.25) and ($sid1_len1 / $sid2_len1 >= 0.8)){
                    $flag = 1;
                }
                if ($flag == 1){
                    if ($num1 >= $num2){
                        push @max1, @max2;
                        @max2 = ();
                    }
                    else{
                        push @max2, @max1;
                        @max1 = ();
                    }
                }
            }
            elsif ((@min2 > 0) and ($sid1_len1 / $sid2_len2 <= 1.4) and ($sid1_len1 / $sid2_len2 >= 0.71)){
                my $flag = 0;
                if (abs ($max1_pos - $min2_pos) <= 150){
                    $flag = 1;
                }
                elsif (((abs ($max1_pos - $min2_pos) <= 300) or ($max1_pos == 0) or ($min2_pos == 0)) and ($sid1_len1 / $sid2_len2 <= 1.25) and ($sid1_len1 / $sid2_len2 >= 0.8)){
                    $flag = 1;
                }
                if ($flag == 1){
                    if ($num1 >= $num2){
                        push @max1, @min2;
                        @min2 = ();
                    }
                    else{
                        push @min2, @max1;
                        @max1 = ();
                    }
                }
            }
            elsif ((@min1 > 0) and ($sid2_len1 / $sid1_len2 <= 1.4) and ($sid2_len1 / $sid1_len2 >= 0.71)){
                my $flag = 0;
                if (abs ($max2_pos - $min1_pos) <= 150){
                    $flag = 1;
                }
                elsif (((abs ($max2_pos - $min1_pos) <= 300) or ($max2_pos == 0) or ($min1_pos == 0)) and ($sid2_len1 / $sid1_len2 <= 1.25) and ($sid2_len1 / $sid1_len2 >= 0.8)){
                    $flag = 1;
                }
                if ($flag == 1){
                    if ($num1 >= $num2){
                        push @min1, @max2;
                        @max2 = ();
                    }
                    else{
                        push @max2, @min1;
                        @min1 = ();
                    }
                }
            }
            if ((@min1 > 0) and (@min2 > 0) and ($sid1_len2 / $sid2_len2 <= 1.4) and ($sid1_len2 / $sid2_len2 >= 0.71)){
                my $flag = 0;
                if (abs ($min1_pos - $min2_pos) <= 150){
                    $flag = 1;
                }
                elsif (((abs ($min1_pos - $min2_pos) <= 300) or ($min1_pos == 0) or ($min2_pos == 0)) and ($sid1_len2 / $sid2_len2 <= 1.25) and ($sid1_len2 / $sid2_len2 >= 0.8)){
                    $flag = 1;
                }
                if ($flag == 1){
                    if ($num1 >= $num2){
                        push @min1, @min2;
                        @min2 = ();
                    }
                    else{
                        push @min2, @min1;
                        @min1 = ();
                    }
                }
            }
            my $new_num1 = @max1 + @min1;
            my $new_num2 = @max2 + @min2;
            if ($num1 != $new_num1){
                if ($new_num2 > 0){
                    @{${$bam_str2{$sid2}}{'INS'}} = (@max2, @min2);
                }
                else{
                    delete ${$bam_str2{$sid2}}{'INS'};
                    delete $bam_str2{$sid2} if (scalar keys %{$bam_str2{$sid2}} == 0);
                }
                if ($new_num1 > 0){
                    @{${$bam_str2{$sid1}}{'INS'}} = (@max1, @min1);
                }
                else{
                    delete ${$bam_str2{$sid1}}{'INS'};
                    delete $bam_str2{$sid1} if (scalar keys %{$bam_str2{$sid1}} == 0);
                    last;
                }
            }
        }
    }
}

foreach my $sid (keys %bam_str){        # assign the sizes and counts (genotype) of INSs and DELs within TRs
    my $sum_del = 0;
    my $sum_ins = 0;
    next if (!exists $bam_str_cov{$sid}) or ($bam_str_cov{$sid} == 0);
    my ($spos, $send, $motif_size) = split (/=/, $STR2{$sid});
    my $slen = $send - $spos + 1;
    my $sum_read = $bam_str_cov{$sid};
    my $sum_read_all = $sum_read;
    my $flag2 = 0;
    $sum_del = scalar @{${$bam_str{$sid}}{'DEL'}} if (exists ${$bam_str{$sid}}{'DEL'});
    $sum_ins = scalar @{${$bam_str{$sid}}{'INS'}} if (exists ${$bam_str{$sid}}{'INS'});
    if ($slen > 1000){
        my $sum_del2 = 0;
        my $sum_ins2 = 0;
        $sum_del2 = scalar @{${$bam_str2{$sid}}{'DEL'}} if (exists ${$bam_str2{$sid}}{'DEL'});
        $sum_ins2 = scalar @{${$bam_str2{$sid}}{'INS'}} if (exists ${$bam_str2{$sid}}{'INS'});
        if (($sum_ins >= $min_str_reads) and ($sum_del >= $min_str_reads)){
            if (($sum_del2 >= $min_str_reads) and ($sum_ins2 >= $min_str_reads)){
                $flag2 = 1;
            }
            elsif (($sum_ins2 >= $min_str_reads * 2) or ($sum_del2 >= $min_str_reads * 2)){
                $flag2 = 1;
            }
        }
        elsif ($sum_ins >= $min_str_reads){
            if ($sum_ins2 >= $min_str_reads){
                $flag2 = 1;
            }
        }
        elsif ($sum_del >= $min_str_reads){
            if ($sum_del2 >= $min_str_reads){
                $flag2 = 1;
            }
        }
        if ($flag2 == 1){
            $sum_del = $sum_del2;
            $sum_ins = $sum_ins2;
            $sum_read = $bam_str_cov2{$sid} if (exists $bam_str_cov2{$sid});
        }
    }
    my $sum_bp5 = 0;
    my $sum_bp3 = 0;
    my $sum_bp35 = 0;
    my $term_bp5 = 0;
    my $term_bp3 = 0;
    my $sum_bp = 0;
    my $max_ins = 0;
    my @bp5;
    my @bp3;
    my $med_bp5 = 0;
    my $med_bp3 = 0;
    if (exists ${$bam_str{$sid}}{'BP5'}){
        $sum_bp5 = scalar @{${$bam_str{$sid}}{'BP5'}};
        foreach (@{${$bam_str{$sid}}{'BP5'}}){
            if ((abs ($spos - $_) <= 100) and ($send - $_ > 200)){
                $term_bp5 ++;
            }
            push @bp5, $_;
        }
        @bp5 = sort {$a <=> $b} @bp5;
        my $hn = int (@bp5 * 0.5);
        $med_bp5 = $bp5[$hn];
    }
    if (exists ${$bam_str{$sid}}{'BP3'}){
        $sum_bp3 = scalar @{${$bam_str{$sid}}{'BP3'}};
        foreach (@{${$bam_str{$sid}}{'BP3'}}){
            if ((abs ($send - $_) <= 100) and ($_ - $spos > 200)){
                $term_bp3 ++;
            }
            push @bp3, $_;
        }
        @bp3 = sort {$a <=> $b} @bp3;
        my $hn = int (@bp3 * 0.5);
        $med_bp3 = $bp3[$hn];
    }
    if (($term_bp5 == 0) and ($sum_bp5 > 0) and ($sum_bp3 == 0)){
        $sum_bp5 = 0;
    }
    if (($term_bp3 == 0) and ($sum_bp3 > 0) and ($sum_bp5 == 0)){
        $sum_bp3 = 0;
    }
    if (($sum_bp5 > 0) and ($sum_bp3 > 0) and ($term_bp5 == 0) and ($term_bp3 == 0) and (abs ($med_bp5 - $med_bp3) > 50)){
        $sum_bp5 = 0;
        $sum_bp3 = 0;
    }
    $sum_bp35 = scalar @{${$bam_str{$sid}}{'BP53'}} if (exists ${$bam_str{$sid}}{'BP53'});
    $sum_bp5 -= int ($sum_bp35 * 0.5) if ($sum_bp5 > 0);
    $sum_bp3 -= int ($sum_bp35 * 0.5) if ($sum_bp3 > 0);
    $sum_bp = $sum_bp5 + $sum_bp3;
    if ($sum_ins > 0){
        foreach my $ilen (@{${$bam_str{$sid}}{'INS'}}){
            $max_ins = $ilen if ($ilen > $max_ins);
        }
        if (($max_ins < 200) and ($sum_bp < $sum_read * 0.2)){
            $sum_bp = 0;
        }
    }
    my $sum_bp_2 = $sum_bp;
    my $ins_rate = int ($sum_ins / $sum_read * 100 + 0.5) / 100;
    my $del_rate = int ($sum_del / $sum_read * 100 + 0.5) / 100;
    my $bp_rate = int ($sum_bp / $sum_read * 100 + 0.5) / 100;
    my ($med_len1, $med_len2) = (0, 0);
    my ($var_num1, $var_num2) = (0, 0);
    my $type1 = 'NA';
    my $type2 = 'NA';
    my $dprate = 0;
    if (($del_rate > 0) and ($ins_rate + $bp_rate > 0) and ($del_rate + $ins_rate + $bp_rate >= $min_str_vrr * 0.3)){
        my @del_len = @{${$bam_str{$sid}}{'DEL'}};
        @del_len = @{${$bam_str2{$sid}}{'DEL'}} if ($flag2 == 1);
        my ($med_lenD1, $med_lenD2, $var_numD1, $var_numD2) = &clust_indel (\@del_len, "$sid-DEL");
        my ($med_lenI1, $med_lenI2, $var_numI1, $var_numI2) = (0, 0, 0, 0);
        if ($ins_rate > 0){
            my @ins_len = @{${$bam_str{$sid}}{'INS'}};
            @ins_len = @{${$bam_str2{$sid}}{'INS'}} if ($flag2 == 1);
            ($med_lenI1, $med_lenI2, $var_numI1, $var_numI2) = &clust_indel (\@ins_len, "$sid-INS");
        }
        if (($slen >= 200) and ($term_bp5 >= $sum_read_all * 0.15) and ($term_bp3 >= $sum_read_all * 0.15)){
            my ($flank_dp) = &calc_dp ($spos, $send, \%DP, 'STR');
            $dprate = int ($sum_read_all / $flank_dp * 100 + 0.5) / 100 if ($flank_dp > 0);
            my $dup_len = 0;
            if ($dprate >= 1.2){
                if ($bp_rate < 0.7){
                    $dup_len = int ($slen * ($dprate * 2 - 2));
                }
                else{
                    $dprate = 2 if ($dprate < 2);
                    $dup_len = int ($slen * ($dprate - 1));
                }
            }
            if ($dup_len > 0){
                my $temp_len = $med_lenI1;
                my $temp_num = $var_numI1;
                if (($med_lenI1 > 0) and ($med_lenI2 > 0)){
                    if ($med_lenI1 >= $med_lenI2){
                        $med_lenI1 = $dup_len if ($dup_len >= $med_lenI1);
                        if ($temp_len < $dup_len * $str_min_len_rate){
                            $med_lenI2 = $temp_len;
                            $var_numI2 = $temp_num;
                            $var_numI1 = 0;
                        }
                    }
                    else{
                        $med_lenI1 = $dup_len if ($dup_len > $med_lenI2);
                        $med_lenI1 = $med_lenI2 if ($dup_len < $med_lenI2);
                        $med_lenI2 = $temp_len;
                        $var_numI1 = $var_numI2;
                        $var_numI2 = $temp_num;
                        if (($med_lenI1 / $med_lenI2 >= $str_min_len_rate) and ($med_lenI1 / $med_lenI2 <= $str_max_len_rate)){
                            $var_numI1 += $var_numI2;
                            $var_numI2 = 0;
                            $med_lenI2 = 0;
                        }
                    }
                }
                elsif ($med_lenI1 > 0){
                    $med_lenI1 = $dup_len if ($dup_len > $med_lenI1);
                    if ($temp_len < $dup_len * $str_min_len_rate){
                        $med_lenI2 = $temp_len;
                        $var_numI2 = $temp_num;
                        $var_numI1 = 0;
                    }
                }
            }
        }
        if ($var_numI1 > $var_numI2 * 4){
            $var_numI2 = 0;
            $med_lenI2 = 0;
        }
        if ($med_lenI1 == 0){
            $var_numI1 = 0;
            if (exists ${$bam_str{$sid}}{'BP5'}){
                foreach (@{${$bam_str{$sid}}{'BP5'}}){
                    delete $bam_str_bp2{$_};
                }
            }
            if (exists ${$bam_str{$sid}}{'BP3'}){
                foreach (@{${$bam_str{$sid}}{'BP3'}}){
                    delete $bam_str_bp1{$_};
                }
            }
        }
        if ($med_lenI1 > $med_lenI2){
            $var_numI1 += $sum_bp if ($med_lenI1 >= 500);
        }
        else{
            if ($med_lenI2 >= 500){
                $var_numI2 += $sum_bp;
                if ($var_numI2 > $var_numI1){
                    my $temp_len = $med_lenI1;
                    my $temp_num = $var_numI1;
                    $med_lenI1 = $med_lenI2;
                    $var_numI1 = $var_numI2;
                    $med_lenI2 = $temp_len;
                    $var_numI2 = $temp_num;
                }
            }
        }
        if ($var_numI1 == 0){
            $med_len1 = $med_lenD1;
            $var_num1 = $var_numD1;
            $type1 = 'DEL';
            if ($var_numD2 > 0){
                $med_len2 = $med_lenD2;
                $var_num2 = $var_numD2;
                $type2 = 'DEL';
            }
            else{
                $med_len2 = 0;
                $var_num2 = 0;
                $type2 = 'NA';
            }
        }
        elsif ($var_numD2 > $var_numI1){
            if ($var_numD2 > $var_numI1 * $str_max_len_rate){
                $med_len1 = $med_lenD1;
                $med_len2 = $med_lenD2;
                $var_num1 = $var_numD1;
                $var_num2 = $var_numD2;
                $type1 = 'DEL';
                $type2 = 'DEL';
            }
            else{
                $med_len1 = $med_lenD1;
                $med_len2 = $med_lenI1;
                $var_num1 = $var_numD1;
                $var_num2 = $var_numI1;
                $type1 = 'DEL';
                $type2 = 'INS';
            }
        }
        elsif ($var_numI2 > $var_numD1){
            if ($var_numI2 > $var_numD1 * $str_max_len_rate){
                $med_len1 = $med_lenI1;
                $med_len2 = $med_lenI2;
                $var_num1 = $var_numI1;
                $var_num2 = $var_numI2;
                $type1 = 'INS';
                $type2 = 'INS';
                if (($med_len1 / $med_len2 >= $str_min_len_rate) and ($med_len1 / $med_len2 <= $str_max_len_rate)){
                    $var_num1 += $var_num2;
                    $var_num2 = $var_numD1;
                    $med_len2 = $med_lenD1;
                    $type2 = 'DEL';
                }
            }
            else{
                $med_len1 = $med_lenI1;
                $med_len2 = $med_lenD1;
                $var_num1 = $var_numI1;
                $var_num2 = $var_numD1;
                $type1 = 'INS';
                $type2 = 'DEL';
            }
        }
        elsif ($var_numD1 >= $var_numI1){
            $med_len1 = $med_lenD1;
            $med_len2 = $med_lenI1;
            $var_num1 = $var_numD1;
            $var_num2 = $var_numI1;
            $type1 = 'DEL';
            $type2 = 'INS';
        }
        else{
            $med_len1 = $med_lenI1;
            $med_len2 = $med_lenD1;
            $var_num1 = $var_numI1;
            $var_num2 = $var_numD1;
            $type1 = 'INS';
            $type2 = 'DEL';
        }
        if (($var_num1 > $var_num2 * 4) or (($var_num1 > $var_num2 * 3) and ($var_num2 < 5))){
            $var_num2 = 0;
            $med_len2 = 0;
            $type2 = 'NA';
        }
    }
    elsif ($del_rate >= $min_str_vrr * 0.3){
        my @del_len = @{${$bam_str{$sid}}{'DEL'}};
        @del_len = @{${$bam_str2{$sid}}{'DEL'}} if ($flag2 == 1);
        ($med_len1, $med_len2, $var_num1, $var_num2) = &clust_indel (\@del_len, "$sid-DEL");
        $type1 = 'DEL';
        $type2 = 'DEL' if ($med_len2 > 0);
    }
    elsif ($ins_rate + $bp_rate >= $min_str_vrr * 0.3){
        if ($sum_ins == 1){
            $med_len1 = ${${$bam_str{$sid}}{'INS'}}[0];
            $sum_bp = 0 if ($med_len1 <= 500);
            $var_num1 = $sum_ins + $sum_bp;
        }
        elsif ($sum_ins > 1){
            my @ins_len = @{${$bam_str{$sid}}{'INS'}};
            @ins_len = @{${$bam_str2{$sid}}{'INS'}} if ($flag2 == 1);
            ($med_len1, $med_len2, $var_num1, $var_num2) = &clust_indel (\@ins_len, "$sid-INS");
            if (($slen >= 200) and ($term_bp5 >= $sum_read_all * 0.15) and ($term_bp3 >= $sum_read_all * 0.15)){     # assign long INS overlapping TR DUP
                my ($flank_dp) = &calc_dp ($spos, $send, \%DP, 'STR');
                $dprate = int ($sum_read_all / $flank_dp * 100 + 0.5) / 100 if ($flank_dp > 0);
                if ($dprate >= 1.2){
                    my $duplen = 0;
                    if ($bp_rate < 0.7){
                        $duplen = int ($slen * ($dprate * 2 - 2));
                    }
                    else{
                        $dprate = 2 if ($dprate < 2);
                        $duplen = int ($slen * ($dprate - 1));
                    }
                    my $temp_len = $med_len1;
                    my $temp_num = $var_num1;
                    if (($med_len1 > 0) and ($med_len2 > 0)){
                        if ($med_len1 > $med_len2){
                            $med_len1 = $duplen if ($duplen > $med_len1);
                            if ($temp_len < $duplen * $str_min_len_rate){
                                $med_len2 = $temp_len;
                                $var_num2 = $temp_num;
                                $var_num1 = $sum_bp;
                            }
                            else{
                                $var_num1 += $sum_bp;
                            }
                        }
                        else{
                            $med_len1 = $duplen if ($duplen >= $med_len2);
                            $med_len1 = $med_len2 if ($duplen < $med_len2);
                            $var_num1 = $sum_bp;
                            $med_len2 = $temp_len;
                            $var_num2 = $temp_num;
                            if (($med_len1 / $med_len2 >= $str_min_len_rate) and ($med_len1 / $med_len2 <= $str_max_len_rate)){
                                $var_num1 += $var_num2;
                                $var_num2 = 0;
                                $med_len2 = 0;
                            }
                        }
                    }
                    elsif ($med_len1 > 0){
                        $med_len1 = $duplen if ($duplen > $med_len1);
                        if ($temp_len < $duplen * $str_min_len_rate){
                            $med_len2 = $temp_len;
                            $var_num2 = $temp_num;
                            $var_num1 = $sum_bp;
                        }
                        else{
                            $var_num1 += $sum_bp;
                        }
                    }
                    $sum_bp = 0;
                }
            }
            if ($var_num1 > $var_num2 * 4){
                $var_num2 = 0;
                $med_len2 = 0;
            }
            if ($sum_bp > 0){
                if (($med_len1 > 500) and ($med_len1 > $med_len2 * 2)){
                    $var_num1 += $sum_bp;
                }
                elsif (($med_len2 > 500) and ($med_len2 > $med_len1 * 2)){
                    $var_num2 += $sum_bp;
                }
                elsif (($med_len1 > 0) and ($med_len2 > 0)){ 
                    if (($med_len1 > 500) and ($med_len2 > 500)){
                        my $sum_hbp = int ($sum_bp * 0.5 + 0.5);
                        $var_num1 += $sum_hbp;
                        $var_num2 += $sum_bp - $sum_hbp;
                    }
                    elsif ($med_len1 > 500){
                        $var_num1 += $sum_bp
                    }
                    elsif ($med_len2 > 500){
                        $var_num2 += $sum_bp
                    }
                }
                else{
                    $var_num1 += $sum_bp if ($med_len1 > 500);
                }
            }
            $type2 = 'INS' if ($med_len2 > 0);
            if (($med_len1 > 0) and ($med_len2 > 0)){
                if ($var_num2 > $var_num1){
                    my $len1 = $med_len1;
                    my $num1 = $var_num1;
                    $med_len1 = $med_len2;
                    $var_num1 = $var_num2;
                    $med_len2 = $len1;
                    $var_num2 = $num1;
                }
            }
        }
        elsif ($sum_ins == 0){
            if (($slen >= 200) and ($term_bp5 >= $sum_read_all * 0.2) and ($term_bp3 >= $sum_read * 0.2)){
                my ($flank_dp) = &calc_dp ($spos, $send, \%DP, 'STR');
                $dprate = int ($sum_read_all / $flank_dp * 100 + 0.5) / 100 if ($flank_dp > 0);
                if ($dprate < 1.2){
                    if (exists ${$bam_str{$sid}}{'BP5'}){
                        foreach (@{${$bam_str{$sid}}{'BP5'}}){
                            delete $bam_str_bp2{$_};
                        }
                    }
                    if (exists ${$bam_str{$sid}}{'BP3'}){
                        foreach (@{${$bam_str{$sid}}{'BP3'}}){
                            delete $bam_str_bp1{$_};
                        }
                    }
                    next;
                }
                $var_num1 = $sum_bp;
                if ($var_num1 / $sum_read < 0.7){
                    $med_len1 = int ($slen * ($dprate * 2 - 2));
                }
                else{
                    $dprate = 2 if ($dprate < 2);
                    $med_len1 = int ($slen * ($dprate - 1));
                }
            }
            else{
                if (exists ${$bam_str{$sid}}{'BP5'}){
                    foreach (@{${$bam_str{$sid}}{'BP5'}}){
                        delete $bam_str_bp2{$_};
                    }
                }
                if (exists ${$bam_str{$sid}}{'BP3'}){
                    foreach (@{${$bam_str{$sid}}{'BP3'}}){
                        delete $bam_str_bp1{$_};
                    }
                }
                next;
            }
        }
        $type1 = 'INS';
    }
    next if ($var_num1 <= 1);
    $bam_indel_str{$spos} = "$sid=$send=$motif_size=$sum_read=$type1=$type2=$med_len1=$med_len2=$var_num1=$var_num2=$sum_bp_2=$dprate" if ($type1 ne 'NA');
    my $hit_read1 = '';
    my $hit_read2 = '';
    my $hit_ipos1 = '';
    my $hit_ipos2 = '';
    my $hit_pos1 = 0;
    my $hit_pos2 = 0;
    if ($type1 eq 'INS'){
        my %ipos;
        my %ipos_info;
        my $opt_ipos = 0;
        my $opt_num = 0;
        foreach my $readid (keys %{$bam_str_insinfo{$sid}}){
            my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo{$sid}}{$readid});
            if (($ilen / $med_len1 <= 2) and ($ilen / $med_len1 >= 0.5)){
                $ipos{$ipos} ++;
                $ipos_info{$ipos} = "$readid=$ipos_str";
            }
        }
        foreach my $readid (keys %{$bam_str_insinfo2{$sid}}){
            my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo2{$sid}}{$readid});
            if (($ilen / $med_len1 <= 1.5) and ($ilen / $med_len1 >= 0.65)){
                $ipos{$ipos} ++;
                $ipos_info{$ipos} = "$readid=$ipos_str";
            }
        }
        foreach my $ipos (sort {$ipos{$b} <=> $ipos{$a}} keys %ipos){
            $opt_ipos = $ipos;
            $opt_num = $ipos{$ipos};
            last;
        }
        if ($opt_num >= 3){
            $hit_pos1 = $opt_ipos;
            ($hit_read1, $hit_ipos1) = split (/=/, $ipos_info{$opt_ipos});
        }
        else{
            foreach my $readid (keys %{$bam_str_insinfo{$sid}}){
                my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo{$sid}}{$readid});
                if (($ilen / $med_len1 <= 1.2) and ($ilen / $med_len1 >= 0.8)){
                    $hit_read1 = $readid;
                    $hit_ipos1 = $ipos_str;
                    $hit_pos1 = $ipos;
                    last;
                }
            }
            if ($hit_read1 eq ''){
                foreach my $readid (keys %{$bam_str_insinfo{$sid}}){
                    my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo{$sid}}{$readid});
                    if (($ilen / $med_len1 <= 1.5) and ($ilen / $med_len1 >= 0.65)){
                        $hit_read1 = $readid;
                        $hit_ipos1 = $ipos_str;
                        $hit_pos1 = $ipos;
                        last;
                    }
                }
            }
            if ($hit_read1 eq ''){
                foreach my $readid (keys %{$bam_str_insinfo{$sid}}){
                    my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo{$sid}}{$readid});
                    if (($ilen / $med_len1 <= 2) and ($ilen / $med_len1 >= 0.5)){
                        $hit_read1 = $readid;
                        $hit_ipos1 = $ipos_str;
                        $hit_pos1 = $ipos;
                        last;
                    }
                }
            }
            if ($hit_read1 eq ''){
                foreach my $readid (keys %{$bam_str_insinfo2{$sid}}){
                    my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo2{$sid}}{$readid});
                    if (($ilen / $med_len1 <= 1.2) and ($ilen / $med_len1 >= 0.8)){
                        $hit_read1 = $readid;
                        $hit_ipos1 = $ipos_str;
                        $hit_pos1 = $ipos;
                        last;
                    }
                }
            }
            if ($hit_read1 eq ''){
                foreach my $readid (keys %{$bam_str_insinfo2{$sid}}){
                    my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo2{$sid}}{$readid});
                    if (($ilen / $med_len1 <= 1.5) and ($ilen / $med_len1 >= 0.65)){
                        $hit_read1 = $readid;
                        $hit_ipos1 = $ipos_str;
                        $hit_pos1 = $ipos;
                        last;
                    }
                }
            }
        }
    }
    if (($type2 eq 'INS') and ($med_len2 > 0)){
        my %ipos;
        my %ipos_info;
        my $opt_ipos = 0;
        my $opt_num = 0;
        foreach my $readid (keys %{$bam_str_insinfo{$sid}}){
            my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo{$sid}}{$readid});
            if (($ilen / $med_len2 <= 2) and ($ilen / $med_len2 >= 0.5)){
                $ipos{$ipos} ++;
                $ipos_info{$ipos} = "$readid=$ipos_str";
            }
        }
        foreach my $readid (keys %{$bam_str_insinfo2{$sid}}){
            my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo2{$sid}}{$readid});
            if (($ilen / $med_len2 <= 1.5) and ($ilen / $med_len2 >= 0.65)){
                $ipos{$ipos} ++;
                $ipos_info{$ipos} = "$readid=$ipos_str";
            }
        }
        foreach my $ipos (sort {$ipos{$b} <=> $ipos{$a}} keys %ipos){
            $opt_ipos = $ipos;
            $opt_num = $ipos{$ipos};
            last;
        }
        if ($opt_num >= 3){
            $hit_pos2 = $opt_ipos;
            ($hit_read2, $hit_ipos2) = split (/=/, $ipos_info{$opt_ipos});
        }
        else{
            foreach my $readid (keys %{$bam_str_insinfo{$sid}}){
                my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo{$sid}}{$readid});
                if (($ilen / $med_len2 <= 1.2) and ($ilen / $med_len2 >= 0.8)){
                    $hit_read2 = $readid;
                    $hit_ipos2 = $ipos_str;
                    $hit_pos2 = $ipos;
                    last;
                }
            }
            if ($hit_read2 eq ''){
                foreach my $readid (keys %{$bam_str_insinfo{$sid}}){
                    my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo{$sid}}{$readid});
                    if (($ilen / $med_len2 <= 1.5) and ($ilen / $med_len2 >= 0.65)){
                        $hit_read2 = $readid;
                        $hit_ipos2 = $ipos_str;
                        $hit_pos2 = $ipos;
                        last;
                    }
                }
            }
            if ($hit_read2 eq ''){
                foreach my $readid (keys %{$bam_str_insinfo{$sid}}){
                    my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo{$sid}}{$readid});
                    if (($ilen / $med_len2 <= 2) and ($ilen / $med_len2 >= 0.5)){
                        $hit_read2 = $readid;
                        $hit_ipos2 = $ipos_str;
                        $hit_pos2 = $ipos;
                        last;
                    }
                }
            }
            if ($hit_read2 eq ''){
                foreach my $readid (keys %{$bam_str_insinfo2{$sid}}){
                    my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo2{$sid}}{$readid});
                    if (($ilen / $med_len2 <= 1.2) and ($ilen / $med_len2 >= 0.8)){
                        $hit_read2 = $readid;
                        $hit_ipos2 = $ipos_str;
                        $hit_pos2 = $ipos;
                        last;
                    }
                }
            }
            if ($hit_read2 eq ''){
                foreach my $readid (keys %{$bam_str_insinfo2{$sid}}){
                    my ($ilen, $ipos_str, $ipos) = split (/=/, ${$bam_str_insinfo2{$sid}}{$readid});
                    if (($ilen / $med_len2 <= 1.5) and ($ilen / $med_len2 >= 0.65)){
                        $hit_read2 = $readid;
                        $hit_ipos2 = $ipos_str;
                        $hit_pos2 = $ipos;
                        last;
                    }
                }
            }
        }
    }
    if ($hit_read1 ne ''){
        ${$bam_ins_str_read{$spos}}{1} = "$hit_read1=$hit_ipos1=$hit_pos1";
    }
    if ($hit_read2 ne ''){
        ${$bam_ins_str_read{$spos}}{2} = "$hit_read2=$hit_ipos2=$hit_pos2";
    }
}

foreach my $sid (keys %bam_str){
    my ($spos, $send) = split (/=/, $STR2{$sid});
    my $bp5 = 0;
    my $bp3 = 0;
    foreach my $tag (keys %{$bam_str{$sid}}){
        if ($tag eq 'BP5'){
            $bp5 += @{${$bam_str{$sid}}{$tag}};
        }
        elsif ($tag eq 'BP3'){
            $bp3 += @{${$bam_str{$sid}}{$tag}};
        }
    }
    if (($bp5 > 3) and ($bp3 == 0)){
        delete $bam_str_bp2{$spos};
    }
    if (($bp3 > 3) and ($bp5 == 0)){
        delete $bam_str_bp1{$spos};
    }
}

# DEL and INS finding

my %delete;
my %bam_ins_len;

if ($min_indel_size <= 3){
    my $pre_pos = 0;
    foreach my $pos (sort {$a <=> $b} keys %bam_ins){
        if (($pre_pos > 0) and ($pos - $pre_pos <= 1)){
            if (@{$bam_ins{$pre_pos}} >= @{$bam_ins{$pos}}){
                push @{$bam_ins{$pre_pos}}, @{$bam_ins{$pos}};
                delete $bam_ins{$pos};
                next;
            }
            else{
                push @{$bam_ins{$pos}}, @{$bam_ins{$pre_pos}};
                delete $bam_ins{$pre_pos};
            }
        }
        $pre_pos = $pos;
    }
    foreach my $pos (sort {$a <=> $b} keys %bam_ins){
        my @len;
        foreach (@{$bam_ins{$pos}}){
            my ($id, $strand, $ilen) = split (/=/, $_);
            push @len, $ilen;
        }
        my $hn = int (@len * 0.5);
        @len = sort {$a <=> $b} @len;
        my $med_len = $len[$hn];
        if ($med_len < $min_indel_size){
            delete $bam_ins{$pos};
        }
        elsif (($med_len <= 3) and ((@len < 3) or (@len < $min_ins_reads))){
            delete $bam_ins{$pos};
        }
    }
}

my $pre_pos = 0;
foreach my $pos (sort {$a <=> $b} keys %bam_ins){
    if (($pre_pos > 0) and ($pos - $pre_pos <= $ins_bp_diff2)){
        if (@{$bam_ins{$pre_pos}} >= @{$bam_ins{$pos}}){
            push @{$bam_ins{$pre_pos}}, @{$bam_ins{$pos}};
            delete $bam_ins{$pos};
            next;
        }
        else{
            push @{$bam_ins{$pos}}, @{$bam_ins{$pre_pos}};
            delete $bam_ins{$pre_pos};
        }
    }
    $pre_pos = $pos;
}
$pre_pos = 0;
foreach my $pos (sort {$a <=> $b} keys %bam_ins){
    if (($pre_pos > 0) and ($pos - $pre_pos <= $ins_bp_diff)){
        if (@{$bam_ins{$pre_pos}} >= @{$bam_ins{$pos}}){
            push @{$bam_ins{$pre_pos}}, @{$bam_ins{$pos}};
            delete $bam_ins{$pos};
            next;
        }
        else{
            push @{$bam_ins{$pos}}, @{$bam_ins{$pre_pos}};
            delete $bam_ins{$pre_pos};
        }
    }
    $pre_pos = $pos;
}

foreach my $pos (sort {$a <=> $b} keys %bam_ins){       # devide or delete highly deviated size of INSs
    next if (!exists $bam_ins{$pos});
    my $max_len = 0;
    my $min_len = 10000000000;
    my @len;
    foreach (@{$bam_ins{$pos}}){
        my ($id, $strand, $ilen) = split (/=/, $_);
        if ($ilen > $max_len){
            $max_len = $ilen;
        }
        if ($ilen < $min_len){
            $min_len = $ilen;
        }
        push @len, $ilen;
    }
    if (($max_len == 0) or (@len == 0)){
        delete $bam_ins{$pos};
        next;
    }
    my $diff = int ($max_len / $min_len * 100 + 0.5) / 100;
    if ($diff >= 1.3){
        my $mergin = $min_len;
        $mergin = 500 if ($mergin < 500);
        my %pos2_ids;
        foreach my $pos2 (sort {$a <=> $b} keys %bam_ins){
            last if ($pos2 > $pos + $mergin);
            next if ($pos2 < $pos - $mergin);
            next if ($pos2 == $pos);
            next if (!exists $bam_ins{$pos2});
            foreach (@{$bam_ins{$pos2}}){
                my ($id, $strand, $ilen) = split (/=/, $_);
                ${$pos2_ids{$id}}{$pos2} = $ilen;
            }
        }
        if (scalar keys %pos2_ids > 0){
            my %remove_ins;
            my @new_info;
            foreach (@{$bam_ins{$pos}}){
                my ($id, $strand, $ilen, $rpos, $ipos) = split (/=/, $_);
                if (exists $pos2_ids{$id}){
                    my $ilen3 = $ilen;
                    foreach my $pos2 (sort {$a <=> $b} keys %{$pos2_ids{$id}}){
                        my $ilen2 = ${$pos2_ids{$id}}{$pos2};
                        $ilen3 += $ilen2;
                    }
                    if (($ilen3 / $max_len >= 0.9) and ($ilen3 / $max_len <= 1.11)){
                        foreach my $pos2 (keys %{$pos2_ids{$id}}){
                            push @{$remove_ins{$pos2}}, $id;
                        }
                        push @new_info, "$id=$strand=$ilen3=$rpos=$ipos";
                    }
                    else{
                        push @new_info, $_;
                    }
                }
                else{
                    push @new_info, $_;
                }
            }
            if (scalar keys %remove_ins > 0){
                @{$bam_ins{$pos}} = (@new_info);
                foreach my $pos2 (keys %remove_ins){
                    if (@{$remove_ins{$pos2}} == @{$bam_ins{$pos2}}){
                        delete $bam_ins{$pos2};
                    }
                    else{
                        my %ids;
                        my @info2;
                        foreach my $id (@{$remove_ins{$pos2}}){
                            $ids{$id} = 1;
                        }
                        foreach (@{$bam_ins{$pos2}}){
                            my ($id) = split (/=/, $_);
                            if (!exists $ids{$id}){
                                push @info2, $_;
                            }
                        }
                        @{$bam_ins{$pos2}} = (@info2);
                    }
                }
                $max_len = 0;
                $min_len = 10000000000;
                @len = ();
                foreach (@{$bam_ins{$pos}}){
                    my ($id, $strand, $ilen) = split (/=/, $_);
                    if ($ilen > $max_len){
                        $max_len = $ilen;
                    }
                    if ($ilen < $min_len){
                        $min_len = $ilen;
                    }
                    push @len, $ilen;
                }
                $diff = int ($max_len / $min_len * 100 + 0.5) / 100;
            }
        }
    }
    if ($diff <= 1.5){
        my $hn = int (@len * 0.5);
        @len = sort {$a <=> $b} @len;
        my $med_len = $len[$hn];
        if ($med_len < $min_indel_size){
            delete $bam_ins{$pos};
        }
        else{
            $bam_ins_len{$pos} = $med_len;
        }
        next;
    }
    my @max;
    my @min;
    my @minlen;
    my @maxlen;
    my @minpos;
    my @maxpos;
    foreach (@{$bam_ins{$pos}}){
        my ($id, $strand, $ilen, $rpos, $ipos) = split (/=/, $_);
        my $max_diff = $max_len - $ilen;
        my $min_diff = $ilen - $min_len;
        if ($max_diff <= $min_diff){
            push @max, $_;
            push @maxlen, $ilen;
            push @maxpos, $ipos;
        }
        else{
            push @min, $_;
            push @minlen, $ilen;
            push @minpos, $ipos;
        }
    }
    my $min_hn = int (@minlen * 0.5);
    my $max_hn = int (@maxlen * 0.5);
    @minlen = sort {$a <=> $b} @minlen;
    @maxlen = sort {$a <=> $b} @maxlen;
    my $med_minlen = $minlen[$min_hn];
    my $med_maxlen = $maxlen[$max_hn];
    my $max_num = scalar @max;
    my $min_num = scalar @min;
    if ($med_minlen < $min_indel_size){
        delete $bam_ins{$pos};
        push @{$bam_ins{$pos}}, @max;
        if ($med_maxlen < $min_indel_size){
            delete $bam_ins{$pos};
            next;
        }
        $bam_ins_len{$pos} = $med_maxlen;
        $min_num = 0;
    }
    
    if (($max_num > 0) and ($min_num > 0)){
        @maxpos = sort {$a <=> $b} @maxpos;
        @minpos = sort {$a <=> $b} @minpos;
        my $med_minpos = $minpos[$min_hn];
        my $med_maxpos = $maxpos[$max_hn];
        if ($med_minpos == $med_maxpos){
            $med_maxpos ++;
        }
        delete $bam_ins{$pos};
        next if ($med_maxlen < $min_indel_size);
        if ($med_minlen < $min_indel_size){
            push @{$bam_ins{$med_maxpos}}, @max;
            $bam_ins_len{$med_maxpos} = $med_maxlen;
            next;
        }
        push @{$bam_ins{$med_minpos}}, @min;
        push @{$bam_ins{$med_maxpos}}, @max;
        $bam_ins_len{$med_minpos} = $med_minlen;
        $bam_ins_len{$med_maxpos} = $med_maxlen;
    }
}

foreach my $pos1 (sort {$a <=> $b} keys %bam_ins){       # merge similar size of INS with distance 100-500 bp
    next if (exists $delete{$pos1});
    my $num1 = scalar @{$bam_ins{$pos1}};
    my %rid;
    foreach (@{$bam_ins{$pos1}}){
        my ($id, $strand, $ilen) = split (/=/, $_);
        $rid{$id} = 1;
    }
    my $len1 = 0;
    $len1 = $bam_ins_len{$pos1} if (exists $bam_ins_len{$pos1});
    if (!exists $bam_ins_len{$pos1}){
        my @len1;
        foreach (@{$bam_ins{$pos1}}){
            my ($id, $strand, $ilen) = split (/=/, $_);
            push @len1, $ilen;
        }
        next if (@len1 == 0);
        @len1 = sort {$a <=> $b} @len1;
        my $hn = int (@len1 * 0.5);
        $len1 = $len1[$hn];
    }
    if ((!defined $len1) or ($len1 == 0)){
        delete $bam_ins{$pos1};
        next;
    }
    my $max_distance = $len1;
    if ($len1 > 1000){
        $max_distance = 700;
    }
    elsif ($len1 > 500){
        $max_distance = 500;
    }
    elsif ($len1 > 200){
        $max_distance = 300;
    }
    elsif ($len1 > 100){
        $max_distance = 200;
    }
    my %hit;
    foreach my $pos2 (sort {$a <=> $b} keys %bam_ins){
        next if (exists $delete{$pos2});
        next if ($pos2 <= $pos1);
        last if ($pos2 > $pos1 + $max_distance);
        my $num2 = scalar @{$bam_ins{$pos2}};
        my $dupid_flag = 0;
        my @id2;
        foreach (@{$bam_ins{$pos2}}){
            my ($id, $strand, $ilen) = split (/=/, $_);
            if (exists $rid{$id}){
                $dupid_flag ++;
            }
            push @id2, $id;
        }
        my $len2 = 0;
        $len2 = $bam_ins_len{$pos2} if (exists $bam_ins_len{$pos2});
        if (!exists $bam_ins_len{$pos2}){
            my @len2;
            foreach (@{$bam_ins{$pos2}}){
                my ($id, $strand, $ilen) = split (/=/, $_);
                push @len2, $ilen;
            }
            next if (@len2 == 0);
            @len2 = sort {$a <=> $b} @len2;
            my $hn = int (@len2 * 0.5);
            $len2 = $len2[$hn];
        }
        if ((!defined $len2) or ($len2 == 0)){
            delete $bam_ins{$pos2};
            next;
        }
        my $rate = int ($len1 / $len2 * 100 + 0.5) / 100;
        if (($rate >= 0.8) and ($rate <= 1.2) and ($dupid_flag < $num1 * 0.5) and ($dupid_flag < $num2 * 0.5)){
            $hit{$pos2} = $num2;
            map{$rid{$_} = 1} @id2; 
        }
    }
    if (scalar keys %hit > 0){
        my $sum_pos = $pos1 * $num1;
        my $total_num = $num1;
        foreach my $pos2 (keys %hit){
            my $num2 = $hit{$pos2};
            $sum_pos += $pos2 * $num2;
            $total_num += $num2;
        }
        my $ave_pos = int ($sum_pos / $total_num + 0.5);
        my $med_pos = 0;
        my $min_diff = 1000000000;
        foreach my $pos (keys %hit){
            my $diff = abs ($pos - $ave_pos);
            if ($diff < $min_diff){
                $min_diff = $diff;
                $med_pos = $pos;
            }
        }
        if (exists $delete{$med_pos}){
            while (1){
                $med_pos ++;
                last if (!exists $delete{$med_pos}) and (!exists $bam_ins{$med_pos});
            }
        }
        $hit{$pos1} = $num1;
        my @len;
        my @pos2;
        foreach my $pos (keys %hit){
            if ($pos != $med_pos){
                push @{$bam_ins{$med_pos}}, @{$bam_ins{$pos}};
                delete $bam_ins{$pos};
                push @pos2, $pos;
                $delete{$pos} = 1;
            }
            foreach (@{$bam_ins{$pos}}){
                my ($id, $strand, $ilen) = split (/=/, $_);
                push @len, $ilen;
            }
        }
        my $hn = int (@len * 0.5);
        @len = sort {$a <=> $b} @len;
        my $med_len = $len[$hn];
        next if ($med_len == 0);
        $bam_ins_len{$med_pos} = $med_len;
    } 
}

for (my $i = 1; $i <= 2; $i++){
    my $pre_pos = 0;
    my $pre_len = 0;
    foreach my $pos (sort {$a <=> $b} keys %bam_ins){       # merge partially overlapped inss and closely neighbored inss
        my $len1 = $bam_ins_len{$pos};
        my $num = scalar @{$bam_ins{$pos}};
        if ($num == 0){
            delete $bam_ins{$pos};
            next;
        }
        my $pre_num = 0;
        $pre_num = scalar @{$bam_ins{$pre_pos}} if (exists $bam_ins{$pre_pos});
        my $len_rate = int ($len1 / $pre_len * 100 + 0.5) / 100 if ($pre_len > 0);
        my $distance = $pos - $pre_pos + 1;
        my $match_flag = 0;
        if (($pre_pos > 0) and ($distance <= $bp_diff) and ($distance < $len1) and ($distance < $pre_len) and ($len_rate <= 1.7) and ($len_rate >= 0.59)){
            $match_flag = 1;
        }
        elsif (($pre_pos > 0) and ($distance <= 200) and ($len1 > $distance) and ($pre_len > $distance) and ($len_rate <= 1.2) and ($len_rate >= 0.8) ){
            $match_flag = 1;
        }
        elsif (($pre_pos > 0) and ($distance < $len1) and ($distance < $pre_len) and ($len_rate <= 1.1) and ($len_rate >= 0.9)){
            $match_flag = 1;
        }
        if ($match_flag == 1){
            my %ids;
            my $share_id = 0;
            foreach (@{$bam_ins{$pos}}){
                my ($id) = split (/=/, $_);
                $ids{$id} = 1;
            }
            foreach (@{$bam_ins{$pre_pos}}){
                my ($id) = split (/=/, $_);
                $share_id ++ if (exists $ids{$id});
            }
            if (($share_id >= $num * 0.5) or ($share_id >= $pre_num * 0.5)){
                $pre_pos = $pos;
                $pre_len = $len1;
                next;
            }
            my $len2 = int (($len1 * $num + $pre_len * $pre_num) / ($num + $pre_num) + 0.5);
            if (($pre_num == 1) and ($num == 1)){
                my ($pre_id, $pre_strand, $pre_ilen) = split (/=/, ${$bam_ins{$pre_pos}}[0]);
                if ($len1 >= $pre_ilen){
                    push @{$bam_ins{$pos}}, @{$bam_ins{$pre_pos}};
                    delete $bam_ins{$pre_pos};
                    $len1 = $len2;
                }
                else{
                    push @{$bam_ins{$pre_pos}}, @{$bam_ins{$pos}};
                    delete $bam_ins{$pos};
                    $pre_len = $len2;
                    next;
                }
            }
            elsif ($pre_num > $num){
                push @{$bam_ins{$pre_pos}}, @{$bam_ins{$pos}};
                delete $bam_ins{$pos};
                $pre_len = $len2;
                next;
            }
            elsif ($pre_num < $num){
                push @{$bam_ins{$pos}}, @{$bam_ins{$pre_pos}};
                delete $bam_ins{$pre_pos};
                $len1 = $len2;
            }
            else{
                my $sum_prelen = 0;
                foreach (@{$bam_ins{$pre_pos}}){
                    my ($pre_id, $pre_strand, $pre_ilen) = split (/=/, $_);
                    $sum_prelen += $pre_ilen;
                }
                my $ave_prelen = int ($sum_prelen / $pre_num);
                if ($len1 >= $ave_prelen){
                    push @{$bam_ins{$pos}}, @{$bam_ins{$pre_pos}};
                    delete $bam_ins{$pre_pos};
                    $len1 = $len2;
                }
                else{
                    push @{$bam_ins{$pre_pos}}, @{$bam_ins{$pos}};
                    delete $bam_ins{$pos};
                    $pre_len = $len2;
                    next;
                }
            }
        }
        $pre_pos = $pos;
        $pre_len = $len1;
    }
}

foreach my $pos (sort {$a <=> $b} keys %bam_ins){
    my $num = scalar @{$bam_ins{$pos}};
    if ($num == 0){
        delete $bam_ins{$pos};
        next;
    }
    my $sum_len = 0;
    my @pos;
    my @len;
    my %ids;
    foreach (@{$bam_ins{$pos}}){
        my ($id, $strand, $ilen, $rpos, $ipos) = split (/=/, $_);
        push @len, $ilen;
        push @pos, $ipos;
        ${$ids{$id}}{$rpos} = $ilen;
    }
    my $hnum = int ($num * 0.5);
    @len = sort {$a <=> $b} @len;
    my @len2 = (@len);
    @pos = sort {$a <=> $b} @pos;
    my $med_ilen = $len[$hnum];
    my $med_ipos = $pos[$hnum];
    if ($med_ilen == 0){
        delete $bam_ins{$pos};
        next;
    }
    $bam_ins_len{$pos} = $med_ilen;
    my %keep_id;
    foreach my $id (keys %ids){
        if (scalar keys %{$ids{$id}} > 1){
            my %diff;
            foreach my $rpos (keys %{$ids{$id}}){
                my $ilen = ${$ids{$id}}{$rpos};
                my $diff = abs ($ilen - $med_ilen);
                $diff{$rpos} = $diff;
            }
            foreach my $rpos (sort {$diff{$a} <=> $diff{$b}} keys %diff){
                ${$keep_id{$id}}{$rpos} = 1;
                last;
            }
        }
    }
    if (scalar keys %keep_id > 0){
        my @info = ();
        @pos = ();
        @len = ();
        foreach (@{$bam_ins{$pos}}){
            my ($id, $strand, $ilen, $rpos, $ipos) = split (/=/, $_);
            if (exists $keep_id{$id}){
                if (exists ${$keep_id{$id}}{$rpos}){
                    push @info, $_;
                    push @pos, $ipos;
                    push @len, $ilen;
                }
            }
            else{
                push @info, $_;
                push @pos, $ipos;
                push @len, $ilen;
            }
        }
        if (@info == 0){
            delete $bam_ins{$pos};
            next;
        }
        $num = scalar @info;
        $hnum = int ($num * 0.5);
        @len = sort {$a <=> $b} @len;
        @pos = sort {$a <=> $b} @pos;
        $med_ilen = $len[$hnum];
        $med_ipos = $pos[$hnum];
        @{$bam_ins{$pos}} = (@info);
    }
    map{$sum_len += $_} @len;
    my $ave_len = int ($sum_len / @len + 0.5);
    $ave_len = int ($ave_len / 10 + 0.5) * 10 if ($ave_len < $min_indel_size);
    if ($ave_len < $min_indel_size){
        delete $bam_ins{$pos};
        next;
    }
    if ($pos != $med_ipos){
        push @{$bam_ins{$med_ipos}}, @{$bam_ins{$pos}};
        delete $bam_ins{$pos};
        delete $bam_ins_len{$pos};
        $bam_ins_len{$med_ipos} = $med_ilen;
    }
}

foreach my $pos1 (sort {$a <=> $b} keys %bam_del){
    my @len1;
    next if (!exists $bam_del{$pos1});
    foreach (@{$bam_del{$pos1}}){
        my ($dlen) = split (/=/, $_);
        push @len1, $dlen;
    }
    if (@len1 == 0){
        delete $bam_del{$pos1};
        next;
    }
    my $hn1 = int (@len1 * 0.5);
    @len1 = sort {$a <=> $b} @len1;
    my $med_len1 = $len1[$hn1];
    my %merge;
    foreach my $pos2 (sort {$a <=> $b} keys %bam_del){
        last if ($pos2 > $pos1 + $ins_bp_diff);
        next if ($pos2 <= $pos1);
        next if (!exists $bam_del{$pos2});
        my $distance = $pos2 - $pos1 + 1;
        my @len2;
        foreach (@{$bam_del{$pos2}}){
            my ($dlen) = split (/=/, $_);
            push @len2, $dlen;
        }
        if (@len2 == 0){
            delete $bam_del{$pos2};
            next;
        }
        my $hn2 = int (@len2 * 0.5);
        @len2 = sort {$a <=> $b} @len2;
        my $med_len2 = $len2[$hn2];
        if (($distance < $med_len1) and ($distance < $med_len2) and ($med_len1 / $med_len2 <= 1.5) and ($med_len1 / $med_len2 >= 0.67)){
            $merge{$pos2} = scalar @len2;
        }
    }
    if (scalar keys %merge > 0){
        $merge{$pos1} = scalar @len1;
        my $top_pos = 0;
        foreach my $dpos (sort {$merge{$b} <=> $merge{$a}} keys %merge){
            $top_pos = $dpos;
            last;
        }
        foreach my $dpos (keys %merge){
            if ($dpos != $top_pos){
                push @{$bam_del{$top_pos}}, @{$bam_del{$dpos}};
                delete $bam_del{$dpos};
            }
        }
    }
}

foreach my $pos (sort {$a <=> $b} keys %bam_del){       # devide or delete highly deviated size of INSs
    my $max_len = 0;
    my $min_len = 10000000000;
    my @len;
    foreach (@{$bam_del{$pos}}){
        my ($dlen) = split (/=/, $_);
        if ($dlen > $max_len){
            $max_len = $dlen;
        }
        if ($dlen < $min_len){
            $min_len = $dlen;
        }
        push @len, $dlen;
    }
    my $diff = int ($max_len / $min_len * 100 + 0.5) / 100;
    if ($diff <= 1.5){
        my $hn = int (@len * 0.5);
        @len = sort {$a <=> $b} @len;
        my $med_len = $len[$hn];
        if ($med_len < $min_indel_size){
            delete $bam_del{$pos};
        }
        next;
    }
    my @max;
    my @min;
    my @minlen;
    my @maxlen;
    my @minpos;
    my @maxpos;
    foreach (@{$bam_del{$pos}}){
        my ($dlen, $dpos) = split (/=/, $_);
        my $max_diff = $max_len - $dlen;
        my $min_diff = $dlen - $min_len;
        if ($max_diff <= $min_diff){
            push @max, $_;
            push @maxlen, $dlen;
            push @maxpos, $dpos;
        }
        else{
            push @min, $_;
            push @minlen, $dlen;
            push @minpos, $dpos;
        }
    }
    my $min_hn = int (@minlen * 0.5);
    my $max_hn = int (@maxlen * 0.5);
    @minlen = sort {$a <=> $b} @minlen;
    @maxlen = sort {$a <=> $b} @maxlen;
    my $med_minlen = $minlen[$min_hn];
    my $med_maxlen = $maxlen[$max_hn];
    my $max_num = scalar @max;
    my $min_num = scalar @min;
    if ($med_minlen < $min_indel_size){
        delete $bam_del{$pos};
        push @{$bam_del{$pos}}, @max;
        if ($med_maxlen < $min_indel_size){
            delete $bam_del{$pos};
            next;
        }
    }
    elsif (($max_num > 0) and ($min_num > 0)){
        @maxpos = sort {$a <=> $b} @maxpos;
        @minpos = sort {$a <=> $b} @minpos;
        my $med_minpos = $maxpos[$max_hn];
        my $med_maxpos = $minpos[$min_hn];
        if ($med_minpos == $med_maxpos){
            $med_maxpos ++;
        }
        delete $bam_del{$pos};
        push @{$bam_del{$med_minpos}}, @max;
        push @{$bam_del{$med_maxpos}}, @min;
    }
}

my %clust_del;
my %used_del;

foreach my $pos (sort {$a <=> $b} keys %bam_del){   # merge clustered dels for each read
    next if (exists $used_del{$pos});
    my $sumlen = 0;
    foreach (@{$bam_del{$pos}}){
        my ($dlen) = split (/=/, $_);
        $sumlen += $dlen;
    }
    my $num = scalar @{$bam_del{$pos}};
    my $len = int ($sumlen / $num + 0.5);
    my $end = $pos + $len - 1;
    foreach my $pos2 (sort {$a <=> $b} keys %bam_del){
        next if (exists $used_del{$pos2});
        last if ($pos2 >= $end);
        next if ($pos2 <= $pos);
        my $sumlen2 = 0;
        foreach (@{$bam_del{$pos2}}){
            my ($dlen2) = split (/=/, $_);
            $sumlen2 += $dlen2;
        }
        my $num2 = scalar @{$bam_del{$pos2}};
        my $len2 = int ($sumlen2 / $num2 + 0.5);
        my $end2 = $pos2 + $len2 - 1;
        if ($pos2 < $end){
            my $overlap = $end - $pos2 + 1;
            $overlap = $len2 if ($end2 < $end);
            if (($overlap >= $len * $min_overlap_rate2) and ($overlap >= $len2 * $min_overlap_rate2)){
                if (!exists $clust_del{$pos}){
                    push @{$clust_del{$pos}}, @{$bam_del{$pos}};
                    push @{$clust_del{$pos}}, @{$bam_del{$pos2}};
                    $used_del{$pos} = 1;
                    $used_del{$pos2} = 1;
                }
                else{
                    push @{$clust_del{$pos}}, @{$bam_del{$pos2}};
                    $used_del{$pos2} = 1;
                }
            }
        }
    }
}

foreach my $pos (keys %bam_del){
    delete $bam_del{$pos} if (exists $used_del{$pos});
}

foreach my $pos (sort {$a <=> $b} keys %clust_del){
    my $sum_pos = 0;
    foreach (@{$clust_del{$pos}}){
        my ($dlen, $dpos) = split (/=/, $_);
        $sum_pos += $dpos;
    }
    my $ave_pos = int ($sum_pos / @{$clust_del{$pos}});
    if (!exists $bam_del{$ave_pos}){
        push @{$bam_del{$ave_pos}}, @{$clust_del{$pos}};
    }
    else{
        while (1){
            $ave_pos ++;
            next if (exists $bam_del{$ave_pos});
            push @{$bam_del{$ave_pos}}, @{$clust_del{$pos}};
            last;
        }
    }
}
%clust_del = ();
%used_del = ();

my %clust_del2;
my %used_del2;

foreach my $pos (sort {$a <=> $b} keys %bam_del){   # merge clustered dels for each read
    next if (exists $used_del2{$pos});
    my $sumlen = 0;
    foreach (@{$bam_del{$pos}}){
        my ($dlen) = split (/=/, $_);
        $sumlen += $dlen;
    }
    my $num = scalar @{$bam_del{$pos}};
    my $len = int ($sumlen / $num + 0.5);
    my $end = $pos + $len - 1;
    foreach my $pos2 (sort {$a <=> $b} keys %bam_del){
        next if (exists $used_del2{$pos2});
        last if ($pos2 >= $end);
        next if ($pos2 <= $pos);
        my $sumlen2 = 0;
        foreach (@{$bam_del{$pos2}}){
            my ($dlen2) = split (/=/, $_);
            $sumlen2 += $dlen2;
        }
        my $num2 = scalar @{$bam_del{$pos2}};
        my $len2 = int ($sumlen2 / $num2 + 0.5);
        my $end2 = $pos2 + $len2 - 1;
        if ($pos2 < $end){
            my $overlap = $end - $pos2 + 1;
            $overlap = $len2 if ($end2 < $end);
            if (($overlap >= $len * $min_overlap_rate) and ($overlap >= $len2 * $min_overlap_rate)){
                if (!exists $clust_del2{$pos}){
                    push @{$clust_del2{$pos}}, @{$bam_del{$pos}};
                    push @{$clust_del2{$pos}}, @{$bam_del{$pos2}};
                    $used_del2{$pos} = 1;
                    $used_del2{$pos2} = 1;
                }
                else{
                    push @{$clust_del2{$pos}}, @{$bam_del{$pos2}};
                    $used_del2{$pos2} = 1;
                }
            }
        }
    }
}

foreach my $pos (keys %bam_del){
    delete $bam_del{$pos} if (exists $used_del2{$pos});
}

foreach my $pos (sort {$a <=> $b} keys %clust_del2){
    my $sum_pos = 0;
    foreach (@{$clust_del2{$pos}}){
        my ($dlen, $dpos) = split (/=/, $_);
        $sum_pos += $dpos;
    }
    my $ave_pos = int ($sum_pos / @{$clust_del2{$pos}});
    if (!exists $bam_del{$ave_pos}){
        push @{$bam_del{$ave_pos}}, @{$clust_del2{$pos}};
    }
    else{
        while (1){
            $ave_pos ++;
            next if (exists $bam_del{$ave_pos});
            push @{$bam_del{$ave_pos}}, @{$clust_del2{$pos}};
            last;
        }
    }
}
%clust_del2 = ();
%used_del2 = ();

foreach my $pos (sort {$a <=> $b} keys %bam_del){   # merge clustered dels for each read
    next if (exists $used_del2{$pos});
    my $sumlen = 0;
    foreach (@{$bam_del{$pos}}){
        my ($dlen) = split (/=/, $_);
        $sumlen += $dlen;
    }
    my $num = scalar @{$bam_del{$pos}};
    my $len = int ($sumlen / $num + 0.5);
    my $end = $pos + $len - 1;
    foreach my $pos2 (sort {$a <=> $b} keys %bam_del){
        next if (exists $used_del2{$pos2});
        last if ($pos2 >= $end);
        next if ($pos2 <= $pos);
        my $sumlen2 = 0;
        foreach (@{$bam_del{$pos2}}){
            my ($dlen2) = split (/=/, $_);
            $sumlen2 += $dlen2;
        }
        my $num2 = scalar @{$bam_del{$pos2}};
        my $len2 = int ($sumlen2 / $num2 + 0.5);
        my $end2 = $pos2 + $len2 - 1;
        if ($pos2 < $end){
            my $len_rate = int ($len / $len2 * 100 + 0.5) / 100;
            if (($len_rate <= 1.1) and ($len_rate >= 0.9)){
                if (!exists $clust_del2{$pos}){
                    push @{$clust_del2{$pos}}, @{$bam_del{$pos}};
                    push @{$clust_del2{$pos}}, @{$bam_del{$pos2}};
                    $used_del2{$pos} = 1;
                    $used_del2{$pos2} = 1;
                }
                else{
                    push @{$clust_del2{$pos}}, @{$bam_del{$pos2}};
                    $used_del2{$pos2} = 1;
                }
            }
        }
    }
}

foreach my $pos (keys %bam_del){
    delete $bam_del{$pos} if (exists $used_del2{$pos});
}

foreach my $pos (sort {$a <=> $b} keys %clust_del2){
    my $sum_pos = 0;
    foreach (@{$clust_del2{$pos}}){
        my ($dlen, $dpos) = split (/=/, $_);
        $sum_pos += $dpos;
    }
    my $ave_pos = int ($sum_pos / @{$clust_del2{$pos}});
    if (!exists $bam_del{$ave_pos}){
        push @{$bam_del{$ave_pos}}, @{$clust_del2{$pos}};
    }
    else{
        while (1){
            $ave_pos ++;
            next if (exists $bam_del{$ave_pos});
            push @{$bam_del{$ave_pos}}, @{$clust_del2{$pos}};
            last;
        }
    }
}
%clust_del2 = ();
%used_del2 = ();

my %bam_del_len;

foreach my $pos (sort {$a <=> $b} keys %bam_del){       # devide or delete highly deviated size of DELs
    my $max_len = 0;
    my $min_len = 10000000000;
    foreach (@{$bam_del{$pos}}){
        my ($dlen) = split (/=/, $_);
        if ($dlen > $max_len){
            $max_len = $dlen;
        }
        if ($dlen < $min_len){
            $min_len = $dlen;
        }
    }
    my $diff = int ($max_len / $min_len * 100 + 0.5) / 100;
    if ((@{$bam_del{$pos}} == 2) and ($diff > 3)){
        delete $bam_del{$pos};
        next;
    }
    next if ($diff <= 1.7);
    my @max;
    my @min;
    my @minlen;
    my @maxlen;
    foreach (@{$bam_del{$pos}}){
        my ($dlen) = split (/=/, $_);
        my $max_diff = $max_len - $dlen;
        my $min_diff = $dlen - $min_len;
        if ($max_diff <= $min_diff){
            push @max, $_;
            push @maxlen, $dlen;
        }
        else{
            push @min, $_;
            push @minlen, $dlen;
        }
    }
    my $sumlen = 0;
    map{$sumlen += $_} @minlen;
    my $avelen = int ($sumlen / @minlen + 0.5);
    $avelen = int ($avelen / 10 + 0.5) * 10 if ($avelen < $min_indel_size);
    my $max_num = scalar @max;
    my $min_num = scalar @min;
    my @max_dpos;
    my @min_dpos;
    if ($avelen < $min_indel_size){
        delete $bam_del{$pos};
        push @{$bam_del{$pos}}, @max;
        my $sumlen2 = 0;
        map{$sumlen2 += $_} @maxlen;
        my $avelen2 = int ($sumlen2 / @maxlen + 0.5);
        $avelen2 = int ($avelen2 / 10 + 0.5) * 10 if ($avelen2 < $min_indel_size);
        if ($avelen2 < $min_indel_size){
            delete $bam_del{$pos};
        }
    }
    elsif ($max_num > $min_num * 2.5){
        delete $bam_del{$pos};
        push @{$bam_del{$pos}}, @max;
    }
    elsif ($min_num > $max_num * 2.5){
        delete $bam_del{$pos};
        push @{$bam_del{$pos}}, @min;
    }
    elsif (($max_num > 0) and ($min_num > 0)){
        foreach (@max){
            my ($dlen, $dpos) = split (/=/, $_);
            push @max_dpos, $dpos;
        }
        foreach (@min){
            my ($dlen, $dpos) = split (/=/, $_);
            push @min_dpos, $dpos;
        }
        my $max_hnum = int ($max_num * 0.5);
        my $min_hnum = int ($min_num * 0.5);
        @max_dpos = sort {$a <=> $b} @max_dpos;
        @min_dpos = sort {$a <=> $b} @min_dpos;
        my $max_med_pos = $max_dpos[$max_hnum];
        my $min_med_pos = $min_dpos[$min_hnum];
        if ($max_med_pos == $min_med_pos){
            $min_med_pos ++;
        }
        delete $bam_del{$pos};
        push @{$bam_del{$max_med_pos}}, @max;
        push @{$bam_del{$min_med_pos}}, @min;
    }
}

foreach my $pos (sort {$a <=> $b} keys %bam_del){
    my $num = scalar @{$bam_del{$pos}};
    if ($num == 1){
        my ($dlen) = split (/=/, ${$bam_del{$pos}}[0]);
        $bam_del_len{$pos} = $dlen;
        next;
    }
    my $hnum = int ($num * 0.5);
    my @pos = ();
    my @len = ();
    my %ids;
    foreach (@{$bam_del{$pos}}){
        my ($dlen, $dpos, $id) = split (/=/, $_);
        push @pos, $dpos;
        push @len, $dlen;
        ${$ids{$id}}{$dpos} = $dlen;
    }
    my $med_dlen = $len[$hnum];
    my %keep_id;
    foreach my $id (keys %ids){
        if (scalar keys %{$ids{$id}} > 1){
            my %diff;
            foreach my $dpos (keys %{$ids{$id}}){
                my $dlen = ${$ids{$id}}{$dpos};
                my $diff = abs ($dlen - $med_dlen);
                $diff{$dpos} = $diff;
            }
            foreach my $dpos (sort {$diff{$a} <=> $diff{$b}} keys %diff){
                ${$keep_id{$id}}{$dpos} = 1;
                last;
            }
        }
    }
    if (scalar keys %keep_id > 0){
        my @info = ();
        @pos = ();
        @len = ();
        foreach (@{$bam_del{$pos}}){
            my ($dlen, $dpos, $id) = split (/=/, $_);
            if (exists $keep_id{$id}){
                if (exists ${$keep_id{$id}}{$dpos}){
                    push @info, $_;
                    push @pos, $dpos;
                    push @len, $dlen;
                }
            }
            else{
                push @info, $_;
                push @pos, $dpos;
                push @len, $dlen;
            }
        }
        if (@info == 0){
            delete $bam_del{$pos};
            next;
        }
        @{$bam_del{$pos}} = (@info);
        $num = scalar @info;
        $hnum = int ($num * 0.5);
    }
    @pos = sort {$a <=> $b} @pos;
    my $med_pos = $pos[$hnum];
    my $sum_len = 0;
    map{$sum_len += $_} @len;
    my $ave_len = int ($sum_len / @len + 0.5);
    $ave_len = int ($ave_len / 10 + 0.5) * 10 if ($ave_len < $min_indel_size);
    if ($ave_len < $min_indel_size){
        delete $bam_del{$pos};
        next;
    }
    if ($pos != $med_pos){
        push @{$bam_del{$med_pos}}, @{$bam_del{$pos}};
        delete $bam_del{$pos};
        $bam_del_len{$med_pos} = $ave_len;
    }
    else{
        $bam_del_len{$pos} = $ave_len;
    }
}

if (scalar keys %bam_rep > 0){
    my $pre_pos = 0;
    foreach my $pos (sort {$a <=> $b} keys %bam_rep){
        if (($pre_pos > 0) and ($pos - $pre_pos > $max_dist)){
            if (@{$bam_rep{$pre_pos}} >= @{$bam_rep{$pos}}){
                push @{$bam_rep{$pre_pos}}, @{$bam_rep{$pos}};
                delete $bam_rep{$pos};
                next;
            }
            else{
                push @{$bam_rep{$pos}}, @{$bam_rep{$pre_pos}};
                delete $bam_rep{$pre_pos};
            }
        }
        $pre_pos = $pos;
    }
}

# select clustered double-clipped reads (>= 3) for DUP and INV
my %unsplit_subread1;
my %unsplit_subread2; 
my %clust_2bp_hcov;
my %clust_2bp_hcov2;
my %clust_2bp_hcov_bp1;
my %clust_2bp_hcov_bp2;
my %clust_2bp_bps;
my %clust_2bp;
my %dup_dp;

foreach my $pos1 (sort {$a <=> $b} keys %bam_2bp){
    my $hit_flag = 0;
    foreach (@{$bam_2bp{$pos1}}){
        my ($end1, $id) = split (/=/, $_);
        if ((exists ${$unsplit_subread1{$id}}{$pos1}) or (exists ${$unsplit_subread1{$id}}{$end1})){
            $hit_flag = 1;
            last;
        }
    }
    next if ($hit_flag == 1);
    push @{$clust_2bp{$pos1}}, @{$bam_2bp{$pos1}};
}
$pre_pos = 0;
foreach my $pos1 (sort {$a <=> $b} keys %clust_2bp){
    if (($pre_pos > 0) and ($pos1 - $pre_pos <= $ins_bp_diff)){
        if (@{$clust_2bp{$pos1}} > @{$clust_2bp{$pre_pos}}){
            push @{$clust_2bp{$pos1}}, @{$clust_2bp{$pre_pos}};
            push @{$clust_2bp_bps{$pos1}}, $pre_pos;
            delete $clust_2bp{$pre_pos};
        }
        elsif (@{$bam_2bp{$pos1}} <= @{$bam_2bp{$pre_pos}}){
            push @{$clust_2bp{$pre_pos}}, @{$clust_2bp{$pos1}};
            push @{$clust_2bp_bps{$pre_pos}}, $pos1;
            delete $clust_2bp{$pos1};
            next;
        }
    }
    $pre_pos = $pos1;
}
foreach my $pos1 (keys %clust_2bp){
    if (@{$clust_2bp{$pos1}} < 2){
        delete $clust_2bp{$pos1};
        next;
    }
    my $sum = 0;
    foreach (@{$clust_2bp{$pos1}}){
        my ($end) = split (/=/, $_);
        $sum += $end;
    }
    my $ave_end = int ($sum / @{$clust_2bp{$pos1}} + 0.5);
    my @ends;
    foreach (@{$clust_2bp{$pos1}}){
        my ($end) = split (/=/, $_);
        if (abs ($ave_end - $end) <= $ins_bp_diff * 0.5){
            push @ends, $_;
        }
    }
    my $bp_distance = $ave_end - $pos1;
    if (@ends >= 2){
        @{$clust_2bp{$pos1}} = (@ends);
        my $dp_rate = 0;
        ($dp_rate) = &calc_dprate ($pos1, $ave_end, \%DP, 'DUP');
        if ($dp_rate >= 2){                 # for DUP
            push @{$clust_2bp_hcov{$pos1}}, @{$clust_2bp{$pos1}};
        }
        if (($dp_rate >= 1.2) and ($bp_distance >= 400)){                 # for DUP
            push @{$clust_2bp_hcov2{$pos1}}, @{$clust_2bp{$pos1}};
        }
    }
}
foreach my $pos1 (keys %clust_2bp_hcov){
    foreach (@{$clust_2bp_hcov{$pos1}}){
        my ($end) = split (/=/, $_);
        $clust_2bp_hcov_bp2{$end} = 1;
    }
}
foreach my $pos1 (keys %clust_2bp_hcov){
    if (exists $clust_2bp_bps{$pos1}){
        $clust_2bp_hcov_bp1{$pos1} = 1;
        foreach my $pos1_2 (@{$clust_2bp_bps{$pos1}}){
            $clust_2bp_hcov_bp1{$pos1_2} = 1;
        }
    }
    else{
        $clust_2bp_hcov_bp1{$pos1} = 1;
    }
}
%clust_2bp_hcov = ();
%clust_2bp = ();
%clust_2bp_bps = ();

# DUP finding

my %dup_clust;
my %clust_dup_bps;
my %used_dupinv_bp;

foreach my $read_id (keys %read_bp1){       # find reads whose clipped sequence at bp1 is mapped at upstream regions, whose clipped sequence at bp2 is mapped at the downstram site
    if (exists $read_bp2{$read_id}){
        my $total_bp_align = scalar keys %{$read_bp1{$read_id}};
        $total_bp_align += scalar keys %{$read_bp2{$read_id}};
        foreach my $bp1 (sort {$a <=> $b} keys %{$read_bp1{$read_id}}){
            my $bp1_p = int ($bp1 / 10 + 0.5) * 10;
            my %bp2p;
            my $flag = 0;
            foreach my $bp2 (sort {$a <=> $b} keys %{$read_bp2{$read_id}}){
                if (($bp2 <= $bp1 + 20000) and ($bp2 >= $bp1 + $bp_diff)){      # when the read has 3'-clip supplementary alignment within 20 Kb dounstream of bp1
                    $flag = 1;
                }
                last if ($bp2 >= $bp1 + 20000);
                my $distance = $bp1 - $bp2 + 1;
                next if ($distance > $max_dup_len);
                next if ($distance < 100);
                next if (exists $bam_str_bp1{$bp1}) and (exists $bam_str_bp2{$bp2} and ($bam_str_bp1{$bp1} eq $bam_str_bp2{$bp2}));
                my ($readpos1, $cliplen1, $cigar1, $strand1) = split (/=/, ${$read_bp1{$read_id}}{$bp1});
                my ($maplen2, $cliplen2, $cigar2, $strand2) = split (/=/, ${$read_bp2{$read_id}}{$bp2});
                next if ($cigar1 eq $cigar2);
#                next if ($strand1 ne $strand2);
                my $dup_bp_diff = 150;
                $dup_bp_diff = $distance * 0.01 if ($distance > 15000);
                my $double_clip_reads = 0;
                foreach my $cpos1 (sort {$a <=> $b} keys %bam_2bp){     # count the dup-supporting reads < double-clipped alignments restricted to the dup region
                    last if ($cpos1 > $bp1 + $dup_bp_diff);
                    next if ($cpos1 < $bp1 - $dup_bp_diff);
                    foreach (@{$bam_2bp{$cpos1}}){
                        my ($cpos2) = split (/=/, $_);
                        if (($cpos2 - $cpos1 + 1 >= $distance * 0.9) and ($cpos2 - $cpos1 + 1 <= $distance * 1.1)){
                            $double_clip_reads ++;
                        }
                    }
                }
                my $diff = $readpos1 - $cliplen2;
                next if ($total_bp_align - $double_clip_reads > 3);
                if ((abs ($diff) > $distance * 0.6) and ($double_clip_reads == 0)){           # check whether the aligned length at the bp1 and the clipped length at the bp2 of the same read alignment are similar
                    next;
                }
                my $ins_match = 0;
                foreach my $ipos (sort {$a <=> $b} keys %bam_ins){
                    last if ($ipos > $bp1 + $bp_diff);
                    next if ($ipos < $bp2 - $bp_diff);
                    my $ilen = $bam_ins_len{$ipos};
                    next if ($ilen < $distance);
                    my $rnum = scalar @{$bam_ins{$ipos}};
                    next if ($rnum < $min_ins_reads);
                    if (($ipos >= $bp1 - $bp_diff) and ($ipos <= $bp2 + $bp_diff)){
                        $ins_match = 1;
                        last;
                    }
                }
                next if ($ins_match == 1);
                my $bp2_p = int ($bp2 / 10 + 0.5) * 10;
                push @{${$dup_clust{$bp2_p}}{$bp1_p}}, "$read_id=$strand2=$cliplen2=$strand1=$cliplen1=$double_clip_reads";
                ${$bp2p{$bp2_p}}{$bp1_p} = 1;
                if (!exists $clust_dup_bps{$bp2_p}){
                    push @{$clust_dup_bps{$bp2_p}}, $bp2;
                }
                else{
                    push @{$clust_dup_bps{$bp2_p}}, $bp2 if (@{$clust_dup_bps{$bp2_p}} < 10);
                }
            }
            if ($flag == 1){
                foreach my $b2 (keys %bp2p){
                    foreach my $b1 (keys %{$bp2p{$b2}}){
                        delete ${$dup_clust{$b2}}{$b1};
                    }
                    delete $dup_clust{$b2} if (scalar keys %{$dup_clust{$b2}} == 0);
                }
                next;
            }
            if (!exists $clust_dup_bps{$bp1_p}){
                push @{$clust_dup_bps{$bp1_p}}, $bp1;
            }
            else{
                push @{$clust_dup_bps{$bp1_p}}, $bp1 if (@{$clust_dup_bps{$bp1_p}} < 10);
            }
        }
    }
}

foreach my $read_id (keys %read_2bp_rpos){      # detect clustered double-clipped reads as DUP
    my $pre_bp1 = 0;
    my %dup_2bp;
    foreach my $bp1 (sort {$a <=> $b} keys %{$read_2bp_rpos{$read_id}}){
        my $bp1_p = int ($bp1 / 10 + 0.5) * 10;
        next if (exists $dup_clust{$bp1_p});
        if (($pre_bp1 > 0) and ($bp1 - $pre_bp1 <= $bp_diff)){
            if (!exists $dup_2bp{$pre_bp1}){
                push @{$dup_2bp{$pre_bp1}}, @{${$read_2bp_rpos{$read_id}}{$pre_bp1}};
            }
            push @{$dup_2bp{$pre_bp1}}, @{${$read_2bp_rpos{$read_id}}{$bp1}};
            next;
        }
        else{
            if (@{${$read_2bp_rpos{$read_id}}{$bp1}} > 1){
                push @{$dup_2bp{$bp1}}, @{${$read_2bp_rpos{$read_id}}{$bp1}};
            }
        }
        $pre_bp1 = $bp1;
    }
    foreach my $bp1 (sort {scalar @{$dup_2bp{$b}} <=> scalar @{$dup_2bp{$a}}} keys %dup_2bp){
        my %rpos;
        my %bp2;
        my @bp2;
        my $strand = '';
        my $clip3_len = 0;
        my $clip5_len = 0;
        foreach (@{$dup_2bp{$bp1}}){
            my ($rstart, $clip3len, $rend, $strand2, $bp2) = split (/=/, $_);
            $rpos{$rstart} = $rend;
            $bp2{$rstart} = $bp2;
            $clip3_len = $clip3len if ($clip3_len == 0);
            $clip5_len = $rstart if ($clip5_len == 0);
            $strand = $strand2 if ($strand eq '');
        }
        my $pre_rpos = 0;
        my $pre_rend = 0;
        my $cn = 0;
        foreach my $rpos (sort {$a <=> $b} keys %rpos){
            my $rend = $rpos{$rpos};
            my $rlen = $rend - $rpos + 1;
            if (($pre_rpos > 0) and (abs ($rpos - $pre_rend) <= $rlen)){
                push @bp2, $bp2{$rpos};
                $cn ++;
            }
            $pre_rpos = $rpos;
            $pre_rend = $rend;
        }
        if ($cn >= 1){
            my $sum_bp2 = 0;
            map{$sum_bp2 += $_} @bp2;
            my $bp2 = int ($sum_bp2 / @bp2 + 0.5);
            my $bp2_p = int ($bp2 / 10 + 0.5) * 10;
            my $bp1_p = int ($bp1 / 10 + 0.5) / 10;
            push @{${$dup_clust{$bp2_p}}{$bp1_p}}, "$read_id=$strand=$clip5_len=$strand=$clip3_len=$cn=2BP";
            push @{$clust_dup_bps{$bp1_p}}, $bp1;
            push @{$clust_dup_bps{$bp2_p}}, @bp2;
            my $len = $bp1 - $bp2 + 1;
        }
    }
}

foreach my $bp1 (sort {$a <=> $b} keys %clust_2bp_hcov2){   # detect >= 400 bp distance double-clipped reads as DUP
    my $bp1_p = int ($bp1 / 10 + 0.5) * 10;
    next if (exists $dup_clust{$bp1_p}); 
    my @bp2;
    foreach (@{$clust_2bp_hcov2{$bp1}}){
        my ($bp2) = split (/=/, $_);
        push @bp2, $bp2;
    }
    my $sum_bp2 = 0;
    map{$sum_bp2 += $_} @bp2;
    my $ave_bp2 = int ($sum_bp2 / @bp2);
    next if ($ave_bp2 - $bp1 < 400);
    foreach (@{$clust_2bp_hcov2{$bp1}}){
        my ($bp2, $read_id, $strand, $clip5_len, $clip3_len) = split (/=/, $_);
        my $bp2_p = int ($bp2 / 10 + 0.5) * 10;
        push @{${$dup_clust{$bp1_p}}{$bp2_p}}, "$read_id=$strand=$clip5_len=$strand=$clip3_len=0=2BP";
        push @{$clust_dup_bps{$bp1_p}}, $bp1;
        push @{$clust_dup_bps{$bp2_p}}, @bp2;
        my $bplen = $bp2 - $bp1 + 1;
    }
}

my $pre_bp1 = 0;
my $pre_num = 0;
foreach my $bp1 (sort {$a <=> $b} keys %dup_clust){            # merge DUP BP1 clusters
    if (($pre_bp1 > 0) and ($bp1 - $pre_bp1 <= $bp_diff)){
        my $num = 0;
        my $bp2_str = '';
        foreach my $bp2 (keys %{$dup_clust{$bp1}}){
            $num += scalar @{${$dup_clust{$bp1}}{$bp2}};
            $bp2_str .= "$bp2,";
        }
        if ($pre_num == 0){
            foreach my $bp2 (keys %{$dup_clust{$pre_bp1}}){
                $pre_num += scalar @{${$dup_clust{$pre_bp1}}{$bp2}};
            }
        }
        if ($pre_num > $num){
            foreach my $bp2 (keys %{$dup_clust{$bp1}}){
                push @{${$dup_clust{$pre_bp1}}{$bp2}}, @{${$dup_clust{$bp1}}{$bp2}};
            }
            $pre_num += $num;
            push @{$clust_dup_bps{$pre_bp1}}, @{$clust_dup_bps{$bp1}} if (@{$clust_dup_bps{$bp1}} < 10);
            delete $dup_clust{$bp1};
            next;
        }
        else{
            foreach my $bp2 (keys %{$dup_clust{$pre_bp1}}){
                push @{${$dup_clust{$bp1}}{$bp2}}, @{${$dup_clust{$pre_bp1}}{$bp2}};
            }
            $pre_num += $num;
            push @{$clust_dup_bps{$bp1}}, @{$clust_dup_bps{$pre_bp1}} if (@{$clust_dup_bps{$pre_bp1}} < 10);
            delete $dup_clust{$pre_bp1};
        }
    }
    else{
        $pre_num = 0;
    }
    $pre_bp1 = $bp1;
}

foreach my $bp1 (sort {$a <=> $b} keys %dup_clust){            # merge DUP BP2 clusters
    my $pre_bp2 = 0;
    my $pre_len = 0;
    foreach my $bp2 (sort {$a <=> $b} keys %{$dup_clust{$bp1}}){
        my $len = $bp2 - $bp1 + 1;
        my $dup_bp_diff = $bp_diff;
        $dup_bp_diff = 50 if ($len < 1000);
        $dup_bp_diff = 150 if ($len > 10000);
        if ($pre_bp2 > 0){
            my $diff = $bp2 - $pre_bp2 + 1;
            if ($diff <= $dup_bp_diff){
                push @{${$dup_clust{$bp1}}{$bp2}}, @{${$dup_clust{$bp1}}{$pre_bp2}};
                push @{$clust_dup_bps{$bp2}}, @{$clust_dup_bps{$pre_bp2}} if (@{$clust_dup_bps{$pre_bp2}} < 10);
                my $pre_num = scalar @{${$dup_clust{$bp1}}{$pre_bp2}};
                my $pre_bp_num = scalar @{$clust_dup_bps{$pre_bp2}};
                my $num = scalar @{${$dup_clust{$bp1}}{$bp2}};
                my $bp_num = scalar @{$clust_dup_bps{$bp2}};
                delete ${$dup_clust{$bp1}}{$pre_bp2};
            }
        }
        $pre_bp2 = $bp2;
        $pre_len = $len;
    }
}

foreach my $bp1 (sort {$a <=> $b} keys %dup_clust){            # delete DUP with shared BP1
    if (scalar keys %{$dup_clust{$bp1}} == 0){
        delete $dup_clust{$bp1};
        next;
    }
    if (scalar keys %{$dup_clust{$bp1}} > 1){
        my $top_num = 0;
        my $sec_num = 0;
        my $top_bp2 = 0;
        my $sec_bp2 = 0;
        foreach my $bp2 (sort {scalar @{${$dup_clust{$bp1}}{$b}} <=> scalar @{${$dup_clust{$bp1}}{$a}}} keys %{$dup_clust{$bp1}}){
            if ($top_num == 0){
                $top_num = scalar @{${$dup_clust{$bp1}}{$bp2}};
                $top_bp2 = $bp2;
                next;
            }
            elsif ($sec_num == 0){
                $sec_num = scalar @{${$dup_clust{$bp1}}{$bp2}};
                $sec_bp2 = $bp2;
                last;
            }
        }
        if ($top_num > $sec_num){
            foreach my $bp2 (keys %{$dup_clust{$bp1}}){
                next if ($bp2 == $top_bp2);
                delete ${$dup_clust{$bp1}}{$bp2};
            }
        }
        else{
            foreach my $bp2 (keys %{$dup_clust{$bp1}}){
                my $num = scalar @{${$dup_clust{$bp1}}{$bp2}};
                next if ($num >= $top_num);
                delete ${$dup_clust{$bp1}}{$bp2};
            }
        }
    }
}   
 
$pre_pos = 0;
my $pre_len = 0;
my @clust_bp1;
foreach my $pos1 (sort {$a <=> $b} keys %dup_clust){        # filter clustered DUPs < 1.5 Kb, where the distance between neighboring DUPs is < DUP length and contain at least 3 DUPs
    my @pos2;
    my $sum = 0;
    if (scalar keys %{$dup_clust{$pos1}} == 0){
        delete $dup_clust{$pos1};
        next;
    }
    foreach my $pos2 (keys %{$dup_clust{$pos1}}){
        push @pos2, $pos2;
    }
    map{$sum += $_} @pos2;
    my $ave_pos2 = int ($sum / @pos2 + 0.5);
    my $len = $ave_pos2 - $pos1;
    next if ($len >= 1500);
    my $distance = $pos1 - $pre_pos;
    if (($pre_pos > 0) and ($distance <= $len) or ($distance <= $pre_len)){
        push @clust_bp1, $pre_pos if (@clust_bp1 == 0);
        push @clust_bp1, $pos1;
    }
    else{
        if (@clust_bp1 >= 3){
            foreach (@clust_bp1){
                delete $dup_clust{$_};
            }
        }
        @clust_bp1 = ();
    }
    $pre_pos = $pos1;
    $pre_len = $len;
}

foreach my $pos1 (sort {$a <=> $b} keys %dup_clust){
    foreach my $pos2 (sort {$a <=> $b} keys %{$dup_clust{$pos1}}){
        next if (exists ${$dup_dp{$pos1}}{$pos2});
        my $dp_rate = 0;                # calc DP rate between flanking and DUP regions
        ($dp_rate) = &calc_dprate ($pos1, $pos2, \%DP, 'DUP');
        if ($dp_rate <= 1.1){
            delete ${$dup_clust{$pos1}}{$pos2};
            delete $dup_clust{$pos1} if (scalar keys %{$dup_clust{$pos1}} == 0);
            next;
        }
        ${$dup_dp{$pos1}}{$pos2} = $dp_rate;
    }
}
my %dup_clust2;
foreach my $pos1 (sort {$a <=> $b} keys %dup_clust){        # select the most reliable DUP from DUPs with the same BP1
    next if (scalar keys %{$dup_clust{$pos1}} <= 1);
    my %bp2_info;
    foreach my $pos2 (keys %{$dup_clust{$pos1}}){
        my $num = scalar @{${$dup_clust{$pos1}}{$pos2}};
        $bp2_info{$pos2} = $num;
    }
    my $top_pos2 = 0;
    my $second_pos2 = 0;
    my $top_num = 0;
    my $second_num = 0;
    foreach my $bp2 (sort {$bp2_info{$b} <=> $bp2_info{$b}} keys %bp2_info){
        if ($top_pos2 == 0){
            $top_pos2 = $bp2;
            $top_num = $bp2_info{$bp2};
        }
        elsif ($second_pos2 == 0){
            $second_pos2 = $bp2;
            $second_num = $bp2_info{$bp2};
            last;
        }
    }
    if (($top_num >= 3) and ($top_num > $second_num)){
    }
    else{
        my $dprate1 = ${$dup_dp{$pos1}}{$top_pos2};
        my $dprate2 = ${$dup_dp{$pos1}}{$second_pos2};
        if ($dprate1 < $dprate2){
            $top_pos2 = $second_pos2;
        }
    }
    foreach my $pos2 (keys %{$dup_clust{$pos1}}){
        if ($pos2 == $top_pos2){
            ${$dup_clust2{$pos2}}{$pos1} = scalar @{${$dup_clust{$pos1}}{$pos2}};
            next;
        }
        delete ${$dup_clust{$pos1}}{$pos2};
        delete $dup_clust{$pos1} if (scalar keys %{$dup_clust{$pos1}} == 0);
    }
}

foreach my $pos2 (sort {$a <=> $b} keys %dup_clust2){        # select the most reliable DUP from DUPs with the same BP2
    next if (scalar keys %{$dup_clust2{$pos2}} == 1);
    my %bp1_info;
    foreach my $pos1 (keys %{$dup_clust2{$pos2}}){
        my $num = ${$dup_clust2{$pos2}}{$pos1};
        $bp1_info{$pos1} = $num;
    }
    my $top_pos1 = 0;
    my $second_pos1 = 0;
    my $top_num = 0;
    my $second_num = 0;
    foreach my $bp1 (sort {$bp1_info{$b} <=> $bp1_info{$b}} keys %bp1_info){
        if ($top_pos1 == 0){
            $top_pos1 = $bp1;
            $top_num = $bp1_info{$bp1};
        }
        elsif ($second_pos1 == 0){
            $second_pos1 = $bp1;
            $second_num = $bp1_info{$bp1};
            last;
        }
    }
    if (($top_num >= 3) and ($top_num > $second_num)){
    }
    else{
        my $dprate1 = ${$dup_dp{$top_pos1}}{$pos2};
        my $dprate2 = ${$dup_dp{$second_pos1}}{$pos2};
        if ($dprate1 < $dprate2){
            $top_pos1 = $second_pos1;
        }
    }
    foreach my $pos1 (keys %{$dup_clust2{$pos2}}){
        if ($pos1 == $top_pos1){
            next;
        }
        delete ${$dup_clust{$pos1}}{$pos2};
        delete $dup_clust{$pos1} if (scalar keys %{$dup_clust{$pos1}} == 0);
    }
}
 
$pre_pos = 0;
my $pre_end = 0;
foreach my $pos1 (sort {$a <=> $b} keys %dup_clust){        # filter overlapped DUPs
    my ($pos2) = keys %{$dup_clust{$pos1}};
    my $len = $pos2 - $pos1 + 1;
    if (($pre_pos > 0) and ($pos1 < $pre_end)){
        my $overlap = $pre_end - $pos1;
        $overlap = $len if ($pos2 < $pre_end);
        my $pre_len = $pre_end - $pre_pos + 1;
        if (($pos2 < $pre_end) or (($overlap >= $len * $min_overlap_rate) and ($overlap >= $pre_len * $min_overlap_rate))){
            my $num = scalar @{${$dup_clust{$pos1}}{$pos2}};
            my $pre_num = scalar @{${$dup_clust{$pre_pos}}{$pre_end}};
            if (($num >= 2) and ($num > $pre_num)){
                delete ${$dup_clust{$pre_pos}}{$pre_end};
                delete $dup_clust{$pre_pos} if (scalar keys %{$dup_clust{$pre_pos}} == 0);
                ${$used_dupinv_bp{$pre_pos}}{$pre_end} = 1 if ($pre_num >= $min_ins_reads);
            }
            else{
                delete ${$dup_clust{$pos1}}{$pos2};
                delete $dup_clust{$pos1} if (scalar keys %{$dup_clust{$pos1}} == 0);
                ${$used_dupinv_bp{$pos1}}{$pos2} = 1 if ($num >= $min_ins_reads);
                next;
            }
        }
    }
    $pre_pos = $pos1;
    $pre_end = $pos2;
}
     
foreach my $bp1 (sort {$a <=> $b} keys %dup_clust){
    my ($bp2) = keys %{$dup_clust{$bp1}};
    my %clust_bp1;
    my %clust_bp2;
    my $pre_sbp1 = 0;
    my $pre_sbp2 = 0;
    my $top_sbp1 = 0;
    my $top_sbp2 = 0;
    my $second_sbp1 = 0;
    my $second_sbp2 = 0;
    my $top_sbp1_num = 0;
    my $top_sbp2_num = 0;
    my $second_sbp1_num = 0;
    my $second_sbp2_num = 0;
    foreach my $sbp1 (sort {$a <=> $b} @{$clust_dup_bps{$bp1}}){
        if ($pre_sbp1 == 0){
            push @{$clust_bp1{$sbp1}}, $sbp1;
            $pre_sbp1 = $sbp1;
            next;
        }
        if ($sbp1 - $pre_sbp1 <= 30){
            push @{$clust_bp1{$pre_sbp1}}, $sbp1;
            next;
        }
        else{
            push @{$clust_bp1{$sbp1}}, $sbp1;
        }
        $pre_sbp1 = $sbp1;
    }
    foreach my $sbp2 (sort {$a <=> $b} @{$clust_dup_bps{$bp2}}){
        if ($pre_sbp2 == 0){
            push @{$clust_bp2{$sbp2}}, $sbp2;
            $pre_sbp2 = $sbp2;
            next;
        }
        if ($sbp2 - $pre_sbp2 <= 30){
            push @{$clust_bp2{$pre_sbp2}}, $sbp2;
            next;
        }
        else{
            push @{$clust_bp2{$sbp2}}, $sbp2;
        }
        $pre_sbp2 = $sbp2;
    }
    foreach my $sbp1 (sort {scalar @{$clust_bp1{$b}} <=> scalar @{$clust_bp1{$a}}} keys %clust_bp1){
        if ($top_sbp1 == 0){
            $top_sbp1 = $sbp1;
            $top_sbp1_num = scalar @{$clust_bp1{$sbp1}};
        }
        elsif ($second_sbp1 == 0){
            $second_sbp1 = $sbp1;
            $second_sbp1_num = scalar @{$clust_bp1{$sbp1}};
            last;
        }
    }
    foreach my $sbp2 (sort {scalar @{$clust_bp2{$b}} <=> scalar @{$clust_bp2{$a}}} keys %clust_bp2){
        if ($top_sbp2 == 0){
            $top_sbp2 = $sbp2;
            $top_sbp2_num = scalar @{$clust_bp2{$sbp2}};
        }
        elsif ($second_sbp2 == 0){
            $second_sbp2 = $sbp2;
            $second_sbp2_num = scalar @{$clust_bp2{$sbp2}};
            last;
        }
    }
    my $select_sbp1 = 0;
    my $select_sbp2 = 0;
    my @sbp1;
    my @sbp2;
    if (($top_sbp1_num == $second_sbp1_num) and ($top_sbp1_num > 1)){
        @sbp1 = (@{$clust_bp1{$top_sbp1}}, @{$clust_bp1{$second_sbp1}});
    }
    elsif ($top_sbp1_num > 1){
        @sbp1 = (@{$clust_bp1{$top_sbp1}});
    }
    else{
        @sbp1 = (@{$clust_dup_bps{$bp1}});
        
    }
    if (($top_sbp2_num == $second_sbp2_num) and ($top_sbp2_num > 1)){
        @sbp2 = (@{$clust_bp2{$top_sbp2}}, @{$clust_bp2{$second_sbp2}});
    }
    elsif ($top_sbp2_num > 1){
        @sbp2 = (@{$clust_bp2{$top_sbp2}});
    }
    else{
        @sbp2 = (@{$clust_dup_bps{$bp2}});
        
    }
    my $bp1_num = scalar @sbp1;
    my $bp1_h = int ($bp1_num / 2);
    @sbp1 = sort {$a <=> $b} @sbp1;
    $select_sbp1 = $sbp1[$bp1_h];
    my $bp2_num = scalar @sbp2;
    my $bp2_h = int ($bp2_num / 2);
    @sbp2 = sort {$a <=> $b} @sbp2;
    $select_sbp2 = $sbp2[$bp2_h];

    my $info2 = '';
    foreach my $info (@{${$dup_clust{$bp1}}{$bp2}}){
        $info2 .= "$info|";
    }
    $info2 =~ s/\|$//;
    ${$bam_dup{$select_sbp1}}{$select_sbp2} = $info2;
}

# INV finding start

my %inv_clust_bp1;
my %inv_clust_bp2;
my %clust_inv_bps1;
my %clust_inv_bps2;

foreach my $id1 (keys %read_bp1){               # find reads whose clipped sequences mapped at the bp1 or bp2 with reverse orientation
    foreach my $bp1 (sort {$a <=> $b} keys %{$read_bp1{$id1}}){
        my ($rend1, $cliplen1, $cigar1, $strand1) = split (/=/, ${$read_bp1{$id1}}{$bp1});
        my $bp1_p = int ($bp1 / 10 + 0.5) * 10;
        foreach my $bp2 (sort {$b <=> $a} keys %{$read_bp1{$id1}}){
            last if ($bp2 < $bp1 + $min_inv_size);
            next if ($bp2 > $bp1 + $max_inv_size);
            next if (exists $bam_str_bp1{$bp1}) and (exists $bam_str_bp2{$bp2} and ($bam_str_bp1{$bp1} eq $bam_str_bp2{$bp2}));
            my ($rend2, $cliplen2, $cigar2, $strand2) = split (/=/, ${$read_bp1{$id1}}{$bp2});
            my $distance = $bp2 - $bp1 + 1;
            next if ($strand1 eq $strand2);
            my $bp2_p = int ($bp2 / 10 + 0.5) * 10;
            push @{${$inv_clust_bp1{$bp1_p}}{$bp2_p}}, "bp1=$id1=$cliplen1=$strand1=$cliplen2=$strand2";
            if (!exists $clust_inv_bps1{$bp2_p}){
                push @{$clust_inv_bps1{$bp2_p}}, $bp2;
            }
            else{
                push @{$clust_inv_bps1{$bp2_p}}, $bp2 if (@{$clust_inv_bps1{$bp2_p}} < 10);
            }
        }
        if (exists $inv_clust_bp1{$bp1_p}){
            if (!exists $clust_inv_bps1{$bp1_p}){
                push @{$clust_inv_bps1{$bp1_p}}, $bp1;
            }
            else{
                push @{$clust_inv_bps1{$bp1_p}}, $bp1 if (@{$clust_inv_bps1{$bp1_p}} < 10);
            }
        }
    }
}
foreach my $id1 (keys %read_bp2){
    foreach my $bp1 (sort {$a <=> $b} keys %{$read_bp2{$id1}}){
        my ($rend1, $cliplen1, $cigar1, $strand1) = split (/=/, ${$read_bp2{$id1}}{$bp1});
        my $bp1_p = int ($bp1 / 10 + 0.5) * 10;
        foreach my $bp2 (sort {$b <=> $a} keys %{$read_bp2{$id1}}){
            last if ($bp2 < $bp1 + $min_inv_size);
            next if ($bp2 > $bp1 + $max_inv_size);
            next if (exists $bam_str_bp1{$bp1}) and (exists $bam_str_bp2{$bp2} and ($bam_str_bp1{$bp1} eq $bam_str_bp2{$bp2}));
            my ($rend2, $cliplen2, $cigar2, $strand2) = split (/=/, ${$read_bp2{$id1}}{$bp2});
            my $distance = $bp2 - $bp1 + 1;
            next if ($strand1 eq $strand2);
            my $rend2_2 = $rend2;
            $rend2_2 += $1 if ($cigar2 =~ /(\d+)[SH]$/);
            my $rend1_2 = $rend1;
            $rend1_2 += $1 if ($cigar1 =~ /(\d+)[SH]$/);
            my $rate2 = int ($rend2_2 / $cliplen1 * 10) / 10;
            my $diff2 = $rend2_2 - $cliplen1;
            my $bp2_p = int ($bp2 / 10 + 0.5) * 10;
            push @{${$inv_clust_bp2{$bp1_p}}{$bp2_p}}, "bp2=$id1=$cliplen1=$strand1=$cliplen2=$strand2";
            if (!exists $clust_inv_bps2{$bp2_p}){
                push @{$clust_inv_bps2{$bp2_p}}, $bp2;
            }
            else{
                push @{$clust_inv_bps2{$bp2_p}}, $bp2 if (@{$clust_inv_bps2{$bp2_p}} < 10);
            }
        }
        if (exists $inv_clust_bp2{$bp1_p}){
            if (!exists $clust_inv_bps2{$bp1_p}){
                push @{$clust_inv_bps2{$bp1_p}}, $bp1;
            }
            else{
                push @{$clust_inv_bps2{$bp1_p}}, $bp1 if (@{$clust_inv_bps2{$bp1_p}} < 10);
            }
        }
    }
}

my %invclust1;
my %invclust2;
my %inv_clust1;
my %inv_clust2;
my %inv_clust;
my $clust_count = 1;
$pre_bp1 = 0;

foreach my $bp1 (sort {$a <=> $b} keys %inv_clust_bp1){     # merge inv bp1 cluster
    if (($pre_bp1 > 0) and ($bp1 - $pre_bp1 <= $bp_diff0)){
        push @{$invclust1{$clust_count}}, $bp1;
        next;
    }
    elsif ($pre_bp1 > 0){
        $clust_count ++;
        push @{$invclust1{$clust_count}}, $bp1;
    }
    push @{$invclust1{$clust_count}}, $bp1 if ($pre_bp1 == 0);
    $pre_bp1 = $bp1;
}
foreach my $clust (sort {$a <=> $b} keys %invclust1){
    my %clust1;
    my %clust2;
    foreach my $bp1 (@{$invclust1{$clust}}){
        foreach my $bp2 (sort {$a <=> $b} keys %{$inv_clust_bp1{$bp1}}){
            $clust1{$bp1} += @{${$inv_clust_bp1{$bp1}}{$bp2}};
        }
    }
    my $top_bp1 = 0;
    my $top_bp2 = 0;
    my %select_bp2;
    foreach my $bp1 (sort {$clust1{$b} <=> $clust1{$a}} keys %clust1){
        $top_bp1 = $bp1 if ($top_bp1 == 0);
        foreach my $bp2 (sort {$a <=> $b} keys %{$inv_clust_bp1{$bp1}}){
            $clust2{$bp2} += @{${$inv_clust_bp1{$bp1}}{$bp2}};
        }
    }
    foreach my $bp2 (sort {$clust2{$b} <=> $clust2{$a}} keys %clust2){
        if ($top_bp2 == 0){
            $top_bp2 = $bp2;
            $select_bp2{$bp2} = 1;
            next;
        }
        if (abs ($bp2 - $top_bp2) <= $bp_diff0){
            $select_bp2{$bp2} = 1;
        }
    }
    
    foreach my $bp1 (@{$invclust1{$clust}}){
        foreach my $bp2 (sort {$a <=> $b} keys %{$inv_clust_bp1{$bp1}}){
            next if (!exists $select_bp2{$bp2});
            push @{${$inv_clust1{$top_bp1}}{$top_bp2}}, @{${$inv_clust_bp1{$bp1}}{$bp2}};
        }
    }
}

$clust_count = 1;
$pre_bp1 = 0;

foreach my $bp1 (sort {$a <=> $b} keys %inv_clust_bp2){     # merge inv bp2 cluster
    if (($pre_bp1 > 0) and ($bp1 - $pre_bp1 <= $bp_diff0)){
        push @{$invclust2{$clust_count}}, $bp1;
        next;
    }
    elsif ($pre_bp1 > 0){
        $clust_count ++;
        push @{$invclust2{$clust_count}}, $bp1;
    }
    push @{$invclust2{$clust_count}}, $bp1 if ($pre_bp1 == 0);
    $pre_bp1 = $bp1;
}
foreach my $clust (sort {$a <=> $b} keys %invclust2){
    my %clust1;
    my %clust2;
    foreach my $bp1 (@{$invclust2{$clust}}){
        foreach my $bp2 (sort {$a <=> $b} keys %{$inv_clust_bp2{$bp1}}){
            $clust1{$bp1} += @{${$inv_clust_bp2{$bp1}}{$bp2}};
        }
    }
    my $top_bp1 = 0;
    my $top_bp2 = 0;
    my %select_bp2;
    foreach my $bp1 (sort {$clust1{$b} <=> $clust1{$a}} keys %clust1){
        $top_bp1 = $bp1 if ($top_bp1 == 0);
        foreach my $bp2 (sort {$a <=> $b} keys %{$inv_clust_bp2{$bp1}}){
            $clust2{$bp2} += @{${$inv_clust_bp2{$bp1}}{$bp2}};
        }
    }
    foreach my $bp2 (sort {$clust2{$b} <=> $clust2{$a}} keys %clust2){
        if ($top_bp2 == 0){
            $top_bp2 = $bp2;
            $select_bp2{$bp2} = 1;
            next;
        }
        if (abs ($bp2 - $top_bp2) <= $bp_diff0){
            $select_bp2{$bp2} = 1;
        }
    }
    
    foreach my $bp1 (@{$invclust2{$clust}}){
        foreach my $bp2 (sort {$a <=> $b} keys %{$inv_clust_bp2{$bp1}}){
            next if (!exists $select_bp2{$bp2});
            push @{${$inv_clust2{$top_bp1}}{$top_bp2}}, @{${$inv_clust_bp2{$bp1}}{$bp2}};
        }
    }
}

%invclust1 = ();
%invclust2 = ();
%inv_clust_bp1 = ();
%inv_clust_bp2 = ();

foreach my $bp1 (sort {$a <=> $b} keys %inv_clust1){        # assign inv bp1 and bp2 with <= 5 Kb distance
    foreach my $bp2 (sort {$a <=> $b} keys %inv_clust2){
        last if ($bp2 > $bp1 + $bp_diff);
        next if ($bp2 < $bp1 - $inv_bp_diff);
        my @info;
        my ($bp1_2) = keys %{$inv_clust1{$bp1}};
        push @info, @{${$inv_clust1{$bp1}}{$bp1_2}};
        my ($bp2_2) = keys %{$inv_clust2{$bp2}};
        next if ($bp2_2 - $bp1 < 100);
        push @info, @{${$inv_clust2{$bp2}}{$bp2_2}};
        next if (@info < 3);
        push @{${$inv_clust{$bp1}}{$bp2_2}}, @info;
    }
}

$pre_pos = 0;
$pre_end = 0;
foreach my $pos1 (sort {$a <=> $b} keys %inv_clust){        # filter overlapped INVs
    my ($pos2) = keys %{$inv_clust{$pos1}};
    my $len = $pos2 - $pos1 + 1;
    my $distance = $pre_pos - $pos1;
    if (($pre_pos > 0) and ($pos1 < $pre_end)){
        my $overlap = $pre_end - $pos1;
        $overlap = $len if ($pos2 < $pre_end);
        my $pre_len = $pre_end - $pre_pos + 1;
        if (($pos2 < $pre_end) or (($overlap >= $len * $min_overlap_rate) or ($overlap >= $pre_len * $min_overlap_rate))){
            my $num = scalar @{${$inv_clust{$pos1}}{$pos2}};
            my $pre_num = scalar @{${$inv_clust{$pre_pos}}{$pre_end}};
            if (($num >= 3) and ($num >= $pre_num)){
                delete ${$inv_clust{$pre_pos}}{$pre_end};
                delete $inv_clust{$pre_pos} if (scalar keys %{$inv_clust{$pre_pos}} == 0);
                ${$used_dupinv_bp{$pre_pos}}{$pre_end} = 1 if ($pre_num >= $min_ins_reads);
            }
            else{
                delete ${$inv_clust{$pos1}}{$pos2};
                delete $inv_clust{$pos1} if (scalar keys %{$inv_clust{$pos1}} == 0);
                ${$used_dupinv_bp{$pos1}}{$pos2} = 1 if ($num >= $min_ins_reads);
                next;
            }
        }
    }
    $pre_pos = $pos1;
    $pre_end = $pos2;
}
$pre_pos = 0;
$pre_len = 0;
my %clust_bp1;
my %clust_bp1_len;
foreach my $pos1 (sort {$a <=> $b} keys %inv_clust){        # filter overlapped INVs
    if (scalar keys %{$inv_clust{$pos1}} == 0){
        delete $inv_clust{$pos1};
        next;
    }
    my ($pos2) = keys %{$inv_clust{$pos1}};
    my $len = $pos2 - $pos1 + 1;
    my $distance = $pos1 - $pre_pos;
    next if ($pre_pos > 0) and (scalar keys %{$inv_clust{$pre_pos}} == 0);
    if (($pre_pos > 0) and ($distance <= $len) and ($distance <= $pre_len)){
        my ($pre_pos2) = keys %{$inv_clust{$pre_pos}};
        my $num = scalar @{${$inv_clust{$pos1}}{$pos2}};
        my $pre_num = scalar @{${$inv_clust{$pre_pos}}{$pre_pos2}};
        $clust_bp1{$pre_pos} = $pre_num if (scalar keys %clust_bp1 == 0);
        $clust_bp1_len{$pre_pos} = $pre_len if (scalar keys %clust_bp1_len == 0);
        $clust_bp1{$pos1} = $num;
        $clust_bp1_len{$pos1} = $len;
    }
    else{
        if (scalar keys %clust_bp1 >= 3){
            my $top_bp1 = 0;
            my $second_bp1 = 0;
            my $top_num = 0;
            my $second_num = 0;
            foreach my $bp1 (sort {$clust_bp1{$b} <=> $clust_bp1{$a}} keys %clust_bp1){
                if ($top_bp1 == 0){
                    $top_bp1 = $bp1;
                    $top_num = $clust_bp1{$bp1};
                }
                elsif ($second_bp1 == 0){
                    $second_bp1 = $bp1;
                    $second_num = $clust_bp1{$bp1};
                    last;
                }
            }
            my $top_len = $clust_bp1_len{$top_bp1};
            my $second_len = $clust_bp1_len{$second_bp1};
            if ($top_num == $second_num){
                if ($top_len > $second_len){
                    $top_bp1 = $second_bp1;
                    $top_len = $second_len;
                }
            }
            foreach my $bp1 (keys %clust_bp1){
                if ($bp1 != $top_bp1){
                    delete $inv_clust{$bp1};
                }
            }
            if ($top_bp1 + $top_len - 1 > $pos2){
                $pre_pos = $top_bp1 + $top_len - 1;
                $pre_len = $top_len;
                next;
            }
        }
        @clust_bp1 = ();
    }
    $pre_pos = $pos1;
    $pre_len = $len;
}

my %inv_bps;
foreach my $bp1 (sort {$a <=> $b} keys %inv_clust){
    my ($bp2) = keys %{$inv_clust{$bp1}};
    next if (!exists $clust_inv_bps1{$bp1}) or (!exists $clust_inv_bps2{$bp2});
    my %clust_bp1;
    my %clust_bp2;
    my $pre_sbp1 = 0;
    my $pre_sbp2 = 0;
    my $top_sbp1 = 0;
    my $top_sbp2 = 0;
    my $second_sbp1 = 0;
    my $second_sbp2 = 0;
    my $top_sbp1_num = 0;
    my $top_sbp2_num = 0;
    my $second_sbp1_num = 0;
    my $second_sbp2_num = 0;
    foreach my $sbp1 (sort {$a <=> $b} @{$clust_inv_bps1{$bp1}}){
        if ($pre_sbp1 == 0){
            push @{$clust_bp1{$sbp1}}, $sbp1;
            $pre_sbp1 = $sbp1;
            next;
        }
        if ($sbp1 - $pre_sbp1 <= 30){
            push @{$clust_bp1{$pre_sbp1}}, $sbp1;
            next;
        }
        else{
            push @{$clust_bp1{$sbp1}}, $sbp1;
        }
        $pre_sbp1 = $sbp1;
    }
    foreach my $sbp2 (sort {$a <=> $b} @{$clust_inv_bps2{$bp2}}){
        if ($pre_sbp2 == 0){
            push @{$clust_bp2{$sbp2}}, $sbp2;
            $pre_sbp2 = $sbp2;
            next;
        }
        if ($sbp2 - $pre_sbp2 <= 30){
            push @{$clust_bp2{$pre_sbp2}}, $sbp2;
            next;
        }
        else{
            push @{$clust_bp2{$sbp2}}, $sbp2;
        }
        $pre_sbp2 = $sbp2;
    }
    foreach my $sbp1 (sort {scalar @{$clust_bp1{$b}} <=> scalar @{$clust_bp1{$a}}} keys %clust_bp1){
        if ($top_sbp1 == 0){
            $top_sbp1 = $sbp1;
            $top_sbp1_num = scalar @{$clust_bp1{$sbp1}};
        }
        elsif ($second_sbp1 == 0){
            $second_sbp1 = $sbp1;
            $second_sbp1_num = scalar @{$clust_bp1{$sbp1}};
            last;
        }
    }
    foreach my $sbp2 (sort {scalar @{$clust_bp2{$b}} <=> scalar @{$clust_bp2{$a}}} keys %clust_bp2){
        if ($top_sbp2 == 0){
            $top_sbp2 = $sbp2;
            $top_sbp2_num = scalar @{$clust_bp2{$sbp2}};
        }
        elsif ($second_sbp2 == 0){
            $second_sbp2 = $sbp2;
            $second_sbp2_num = scalar @{$clust_bp2{$sbp2}};
            last;
        }
    }
    my $select_sbp1 = 0;
    my $select_sbp2 = 0;
    my @sbp1;
    my @sbp2;
    if (($top_sbp1_num == $second_sbp1_num) and ($top_sbp1_num > 1)){
        @sbp1 = (@{$clust_bp1{$top_sbp1}}, @{$clust_bp1{$second_sbp1}});
    }
    elsif ($top_sbp1_num > 1){
        @sbp1 = (@{$clust_bp1{$top_sbp1}});
    }
    else{
        @sbp1 = (@{$clust_inv_bps1{$bp1}});
        
    }
    if (($top_sbp2_num == $second_sbp2_num) and ($top_sbp2_num > 1)){
        @sbp2 = (@{$clust_bp2{$top_sbp2}}, @{$clust_bp2{$second_sbp2}});
    }
    elsif ($top_sbp2_num > 1){
        @sbp2 = (@{$clust_bp2{$top_sbp2}});
    }
    else{
        @sbp2 = (@{$clust_inv_bps2{$bp2}});
        
    }
    my $bp1_num = scalar @sbp1;
    my $bp1_h = int ($bp1_num / 2);
    @sbp1 = sort {$a <=> $b} @sbp1;
    $select_sbp1 = $sbp1[$bp1_h];
    my $bp2_num = scalar @sbp2;
    my $bp2_h = int ($bp2_num / 2);
    @sbp2 = sort {$a <=> $b} @sbp2;
    $select_sbp2 = $sbp2[$bp2_h];

    my $info2 = '';
    foreach my $info (@{${$inv_clust{$bp1}}{$bp2}}){
        $info2 .= "$info|";
    }
    $info2 =~ s/\|$//;
    if (@{${$inv_clust{$bp1}}{$bp2}} >= $min_ins_reads){
        foreach my $bp_1 (@{$clust_inv_bps1{$bp1}}){
            $inv_bps{$bp_1} = 1;
        }
        foreach my $bp_2 (@{$clust_inv_bps2{$bp2}}){
            $inv_bps{$bp_2} = 1;
        }
    }
    ${$bam_inv{$select_sbp1}}{$select_sbp2} = $info2;
}
%clust_inv_bps1 = ();
%clust_inv_bps2 = ();

# For DEL and INS finding

my %DEL;
my %used;
my %bam_del_bp;
my %bam_ins_bp;
my %bam_rep_bp;
my %del_bp;
my %rep_bp;
my %ins_bp_inread;
my %clust_del_bps;
my %clust_rep_bps;
my %clust_ins_bps;
my %rep_bp_info;
my %pseudo_2bp;

foreach my $id (keys %read_bp1){                
    foreach my $bp1 (sort {$a <=> $b} keys %{$read_bp1{$id}}){ 
        if (exists $read_bp2{$id}){     # assign DELs/REPs with two clustered breakpoints (Minialign2 maps also clipped read regions with the same read id)
            my $bp1_p = int ($bp1 / 10 + 0.5) * 10;
            my ($readpos1, $cliplen1, $cigar1, $strand1) = split (/=/, ${$read_bp1{$id}}{$bp1});
            if (exists $read_2bp{$id}){         # detect bp1 corresponding to clustered double-clipped ends
                foreach my $bp_1 (keys %{$read_2bp{$id}}){
                    my ($bp_2) = keys %{${$read_2bp{$id}}{$bp_1}};
                    if ($bp_2 == $bp1){
                        my $bp_1_p = int ($bp_1 / 50 + 0.5) * 50;
                        ${$pseudo_2bp{$bp1_p}}{$bp_1_p} ++;
                    }
                }
            }
            foreach my $bp2 (sort {$a<=> $b} keys %{$read_bp2{$id}}){
                last if ($bp2 > $bp1 + $max_del_size);
                last if ($bp2 > $bp1 - 20000) and ($bp2 < $bp1 - $bp_diff2); # when the read has 5'-clip supplementary alignment end within 20 Kb downstream of bp1
                next if ($bp2 < $bp1 - $bp_diff2);
                next if (exists $bam_str_bp1{$bp1}) and (exists $bam_str_bp2{$bp2} and ($bam_str_bp1{$bp1} eq $bam_str_bp2{$bp2}));
                my ($maplen2, $cliplen2, $cigar2, $strand2) = split (/=/, ${$read_bp2{$id}}{$bp2});
                next if ($strand1 ne $strand2);
                next if (exists $inv_bps{$bp1}) or (exists $inv_bps{$bp2});
                next if (exists ${$shortmap{$id}}{$bp1}) and (exists ${$shortmap{$id}}{$bp2});
                next if (exists $clust_2bp_hcov_bp1{$bp2}) or (exists $clust_2bp_hcov_bp2{$bp1});
                my $bp2_p = int ($bp2 / 10 + 0.5) * 10;
                my $distance = $bp2 - $bp1 + 1;
                my $hit_flag = 0;
                my $diff = $readpos1 - $cliplen2;
                if ($distance >= $bp_diff2){   # check whether the aligned length at the bp1 and the clipped length at the bp2 of the same read alignment are similar
                    if ($diff < 0){
                        push @{${$bam_indel{$bp1_p}}{$bp2_p}}, -$diff ;
                    }
                    ${$del_bp{$bp1_p}}{$bp2_p} ++;
                    if (!exists $clust_del_bps{$bp2_p}){
                        push @{$clust_del_bps{$bp2_p}}, $bp2;
                    }
                    else{
                        push @{$clust_del_bps{$bp2_p}}, $bp2 if (@{$clust_del_bps{$bp2_p}} < 10);
                    }
                    $hit_flag = 1;
                }
                elsif (abs ($distance) < 100){           # INSs with BP1 and BP2 in the same read
                    my $maplen2_2 = $maplen2;
                    $maplen2_2 += $1 if ($cigar2 =~ /(\d+)[SH]$/);
                    my $ins_len = $cliplen1 - $maplen2_2 - $distance;
                    if (($ins_len >= 200) and ($ins_len < 5000)){
                        my $bp_m = int (($bp1 + $bp2) / 2 / 10 + 0.5) * 10;
                        if ((!exists $bam_str_bp1{$bp1}) and (!exists $bam_str_bp2{$bp2})){
                            push @{$ins_bp_inread{$bp_m}}, "$id=$strand1=$ins_len=$readpos1";
                            push @{$clust_ins_bps{$bp_m}}, $bp1, $bp2;
                            $hit_flag = 1;
                        }
                    }
                }
                if (($hit_flag == 0) and ($distance >= $ins_bp_diff) and ($distance < 1000)){
                    my $indel_len = $cliplen2 - $readpos1;
                    if (($indel_len >= $ins_bp_diff) and (($indel_len > $distance * 0.5) and ($indel_len < $distance * 1.1))){
                        ${$rep_bp{$bp1_p}}{$bp2_p} ++;      # conversion/replacement (indel) event
                        push @{${$rep_bp_info{$bp1_p}}{$bp2_p}}, "$id~$strand1~$cliplen1~$cliplen2";
                        if (!exists $clust_rep_bps{$bp2_p}){
                            push @{$clust_rep_bps{$bp2_p}}, $bp2;
                        }
                        else{
                            push @{$clust_rep_bps{$bp2_p}}, $bp2 if (@{$clust_rep_bps{$bp2_p}} < 10);
                        }
                    }
                }
                last;
            }
            if (!exists $clust_del_bps{$bp1_p}){
                push @{$clust_del_bps{$bp1_p}}, $bp1;
            }
            else{
                push @{$clust_del_bps{$bp1_p}}, $bp1 if (@{$clust_del_bps{$bp1_p}} < 10) and (exists $del_bp{$bp1_p});
            }
            if (!exists $clust_rep_bps{$bp1_p}){
                push @{$clust_rep_bps{$bp1_p}}, $bp1;
            }
            else{
                push @{$clust_rep_bps{$bp1_p}}, $bp1 if (@{$clust_rep_bps{$bp1_p}} < 10) and (exists $rep_bp{$bp1_p});
            }
        }
    }
}

if (scalar keys %pseudo_2bp > 0){
    foreach my $bp1 (keys %pseudo_2bp){
        my $pseudo_flag = 0;
        foreach my $bp2 (keys %{$pseudo_2bp{$bp1}}){
            my $readnum = 0;
            foreach my $bp2b (keys %{$del_bp{$bp1}}){
                $readnum += ${$del_bp{$bp1}}{$bp2b};
            }
            if (${$pseudo_2bp{$bp1}}{$bp2} >= $readnum * 0.9){
                $pseudo_flag = 1;
                last;
            }
        }
        if ($pseudo_flag == 1){
            delete $del_bp{$bp1};
        }
    }
}

$pre_bp1 = 0;
$pre_num = 0;
foreach my $bp1 (sort {$a <=> $b} keys %del_bp){            # merge DEL BP1 clusters
    if (($pre_bp1 > 0) and ($bp1 - $pre_bp1 <= $bp_diff)){
        my $num = 0;
        foreach my $bp2 (keys %{$del_bp{$bp1}}){
            $num += ${$del_bp{$bp1}}{$bp2};
        }
        if ($pre_num == 0){
            foreach my $bp2 (keys %{$del_bp{$pre_bp1}}){
                $pre_num += ${$del_bp{$pre_bp1}}{$bp2};
            }
        }
        if ($pre_num > $num){
            foreach my $bp2 (keys %{$del_bp{$bp1}}){
                ${$del_bp{$pre_bp1}}{$bp2} += ${$del_bp{$bp1}}{$bp2};
            }
            $pre_num += $num;
            push @{$clust_del_bps{$pre_bp1}}, @{$clust_del_bps{$bp1}} if (@{$clust_del_bps{$bp1}} < 10);
            delete $del_bp{$bp1};
            next;
        }
        else{
            foreach my $bp2 (keys %{$del_bp{$pre_bp1}}){
                ${$del_bp{$bp1}}{$bp2} += ${$del_bp{$pre_bp1}}{$bp2};
            }
            $pre_num += $num;
            push @{$clust_del_bps{$bp1}}, @{$clust_del_bps{$pre_bp1}} if (@{$clust_del_bps{$pre_bp1}} < 10);
            delete $del_bp{$pre_bp1};
        }
    }
    else{
        $pre_num = 0;
    }
    $pre_bp1 = $bp1;
}

foreach my $bp1 (sort {$a <=> $b} keys %del_bp){            # merge DEL BP2 clusters
    my $pre_bp2 = 0;
    my $pre_len = 0;
    foreach my $bp2 (sort {$a <=> $b} keys %{$del_bp{$bp1}}){
        my $len = $bp2 - $bp1 + 1;
        my $del_bp_diff = $bp_diff;
        $del_bp_diff = 50 if ($len < 1000);
        $del_bp_diff = 150 if ($len > 10000);
        if ($pre_bp2 > 0){
            my $diff = $bp2 - $pre_bp2 + 1;
            if ($diff <= $del_bp_diff){
                ${$del_bp{$bp1}}{$bp2} += ${$del_bp{$bp1}}{$pre_bp2};
                push @{$clust_del_bps{$bp2}}, @{$clust_del_bps{$pre_bp2}} if (@{$clust_del_bps{$pre_bp2}} < 10);
                delete ${$del_bp{$bp1}}{$pre_bp2};
            }
        }
        $pre_bp2 = $bp2;
        $pre_len = $len;
    }
}

foreach my $bp1 (sort {$a <=> $b} keys %del_bp){            # delete DEL with shared BP1
    if (scalar keys %{$del_bp{$bp1}} == 0){
        delete $del_bp{$bp1};
        next;
    }
    if (scalar keys %{$del_bp{$bp1}} > 1){
        my $top_num = 0;
        my $sec_num = 0;
        my $top_bp2 = 0;
        my $sec_bp2 = 0;
        foreach my $bp2 (sort {${$del_bp{$bp1}}{$b} <=> ${$del_bp{$bp1}}{$a}} keys %{$del_bp{$bp1}}){
            if ($top_num == 0){
                $top_num = ${$del_bp{$bp1}}{$bp2};
                $top_bp2 = $bp2;
                next;
            }
            elsif ($sec_num == 0){
                $sec_num = ${$del_bp{$bp1}}{$bp2};
                $sec_bp2 = $bp2;
                last;
            }
        }
        if ($top_num > $sec_num){
            foreach my $bp2 (keys %{$del_bp{$bp1}}){
                next if ($bp2 == $top_bp2);
                delete ${$del_bp{$bp1}}{$bp2};
            }
        }
        else{
            foreach my $bp2 (keys %{$del_bp{$bp1}}){
                my $num = ${$del_bp{$bp1}}{$bp2};
                next if ($num >= $top_num);
                delete ${$del_bp{$bp1}}{$bp2};
            }
        }
    }
}

my %del_clust2;
foreach my $pos1 (sort {$a <=> $b} keys %del_bp){        # select the most reliable DEL from DELs with the same BP1
    if (scalar keys %{$del_bp{$pos1}} == 0){
        delete $del_bp{$pos1};
        next;
    }
    next if (scalar keys %{$del_bp{$pos1}} <= 1);
    my %bp2_info;
    foreach my $pos2 (keys %{$del_bp{$pos1}}){
        my $num = ${$del_bp{$pos1}}{$pos2};
        $bp2_info{$pos2} = $num;
    }
    my $top_pos2 = 0;
    my $second_pos2 = 0;
    my $top_num = 0;
    my $second_num = 0;
    foreach my $bp2 (sort {$bp2_info{$b} <=> $bp2_info{$b}} keys %bp2_info){
        if ($top_pos2 == 0){
            $top_pos2 = $bp2;
            $top_num = $bp2_info{$bp2};
        }
        elsif ($second_pos2 == 0){
            $second_pos2 = $bp2;
            $second_num = $bp2_info{$bp2};
            last;
        }
    }
    if ($top_num == $second_num){
        if ($second_pos2 < $top_pos2){
            $top_pos2 = $second_pos2;
        }
    }
    foreach my $pos2 (keys %{$del_bp{$pos1}}){
        if ($pos2 == $top_pos2){
            ${$del_clust2{$pos2}}{$pos1} = ${$del_bp{$pos1}}{$pos2};
            next;
        }
        delete ${$del_bp{$pos1}}{$pos2};
        delete $del_bp{$pos1} if (scalar keys %{$del_bp{$pos1}} == 0);
    }
}
foreach my $pos2 (sort {$a <=> $b} keys %del_clust2){        # select the most reliable DEL from DELs with the same BP2
    next if (scalar keys %{$del_clust2{$pos2}} == 1);
    my %bp1_info;
    foreach my $pos1 (keys %{$del_clust2{$pos2}}){
        my $num = ${$del_clust2{$pos2}}{$pos1};
        $bp1_info{$pos1} = $num;
    }
    my $top_pos1 = 0;
    my $second_pos1 = 0;
    my $top_num = 0;
    my $second_num = 0;
    foreach my $bp1 (sort {$bp1_info{$b} <=> $bp1_info{$b}} keys %bp1_info){
        if ($top_pos1 == 0){
            $top_pos1 = $bp1;
            $top_num = $bp1_info{$bp1};
        }
        elsif ($second_pos1 == 0){
            $second_pos1 = $bp1;
            $second_num = $bp1_info{$bp1};
            last;
        }
    }
    if ($top_num == $second_num){
        if ($second_pos1 > $top_pos1){
            $top_pos1 = $second_pos1;
        }
    }
    foreach my $pos1 (keys %{$del_clust2{$pos2}}){
        if ($pos1 == $top_pos1){
            next;
        }
        delete ${$del_bp{$pos1}}{$pos2};
        delete $del_bp{$pos1} if (scalar keys %{$del_bp{$pos1}} == 0);
    }
}
$pre_pos = 0;
$pre_end = 0;
foreach my $pos1 (sort {$a <=> $b} keys %del_bp){        # filter overlapped DELs
    my ($pos2) = keys %{$del_bp{$pos1}};
    my $len = $pos2 - $pos1 + 1;
    if (($pre_pos > 0) and ($pos1 < $pre_end)){
        my $overlap = $pre_end - $pos1;
        $overlap = $len if ($pos2 < $pre_end);
        my $pre_len = $pre_end - $pre_pos + 1;
        if (($pos2 < $pre_end) or (($overlap >= $len * $min_overlap_rate) and ($overlap >= $pre_len * $min_overlap_rate))){
            my $num = ${$del_bp{$pos1}}{$pos2};
            my $pre_num = ${$del_bp{$pre_pos}}{$pre_end};
            if (($num >= 2) and ($num >= $pre_num)){
                delete ${$del_bp{$pre_pos}}{$pre_end};
                delete $del_bp{$pre_pos} if (scalar keys %{$del_bp{$pre_pos}} == 0);
            }
            else{
                delete ${$del_bp{$pos1}}{$pos2};
                delete $del_bp{$pos1} if (scalar keys %{$del_bp{$pos1}} == 0);
                next;
            }
        }
    }
    $pre_pos = $pos1;
    $pre_end = $pos2;
}

foreach my $bp1 (sort {$a <=> $b} keys %del_bp){
    my ($bp2) = keys %{$del_bp{$bp1}};
    my $read_num = ${$del_bp{$bp1}}{$bp2};
    my %clust_bp1;
    my %clust_bp2;
    my %clust_bp1_num;
    my %clust_bp2_num;
    my $pre_sbp1 = 0;
    my $pre_sbp2 = 0;
    my $top_sbp1 = 0;
    my $top_sbp2 = 0;
    my $second_sbp1 = 0;
    my $second_sbp2 = 0;
    my $top_sbp1_num = 0;
    my $top_sbp2_num = 0;
    my $second_sbp1_num = 0;
    my $second_sbp2_num = 0;
    foreach my $sbp1 (sort {$a <=> $b} @{$clust_del_bps{$bp1}}){
        if ($pre_sbp1 == 0){
            push @{$clust_bp1{$sbp1}}, $sbp1;
            $clust_bp1_num{$sbp1} += ${$del_bp{$bp1}}{$bp2};
            $pre_sbp1 = $sbp1;
            next;
        }
        if ($sbp1 - $pre_sbp1 <= 30){
            push @{$clust_bp1{$pre_sbp1}}, $sbp1;
            $clust_bp1_num{$pre_sbp1} += ${$del_bp{$bp1}}{$bp2};
            next;
        }
        else{
            push @{$clust_bp1{$sbp1}}, $sbp1;
            $clust_bp1_num{$sbp1} += ${$del_bp{$bp1}}{$bp2};
        }
        $pre_sbp1 = $sbp1;
    }
    foreach my $sbp2 (sort {$a <=> $b} @{$clust_del_bps{$bp2}}){
        if ($pre_sbp2 == 0){
            push @{$clust_bp2{$sbp2}}, $sbp2;
            $clust_bp2_num{$sbp2} += ${$del_bp{$bp1}}{$bp2};
            $pre_sbp2 = $sbp2;
            next;
        }
        if ($sbp2 - $pre_sbp2 <= 30){
            push @{$clust_bp2{$pre_sbp2}}, $sbp2;
            $clust_bp2_num{$pre_sbp2} += ${$del_bp{$bp1}}{$bp2};
            next;
        }
        else{
            push @{$clust_bp2{$sbp2}}, $sbp2;
            $clust_bp2_num{$sbp2} += ${$del_bp{$bp1}}{$bp2};
        }
        $pre_sbp2 = $sbp2;
    }
    foreach my $sbp1 (sort {$clust_bp1_num{$b} <=> $clust_bp1_num{$a}} keys %clust_bp1_num){
        if ($top_sbp1 == 0){
            $top_sbp1 = $sbp1;
            $top_sbp1_num = $clust_bp1_num{$sbp1};
        }
        elsif ($second_sbp1 == 0){
            $second_sbp1 = $sbp1;
            $second_sbp1_num = $clust_bp1_num{$sbp1};
            last;
        }
    }
    foreach my $sbp2 (sort {$clust_bp2_num{$b} <=> $clust_bp2_num{$a}} keys %clust_bp2_num){
        if ($top_sbp2 == 0){
            $top_sbp2 = $sbp2;
            $top_sbp2_num = $clust_bp2_num{$sbp2};
        }
        elsif ($second_sbp2 == 0){
            $second_sbp2 = $sbp2;
            $second_sbp2_num = $clust_bp2_num{$sbp2};
            last;
        }
    }
    my $select_sbp1 = 0;
    my $select_sbp2 = 0;
    my @sbp1;
    my @sbp2;
    if (($top_sbp1_num == $second_sbp1_num) and ($top_sbp1_num > 1)){
        @sbp1 = (@{$clust_bp1{$top_sbp1}}, @{$clust_bp1{$second_sbp1}});
    }
    elsif ($top_sbp1_num > 1){
        @sbp1 = (@{$clust_bp1{$top_sbp1}});
    }
    else{
        @sbp1 = (@{$clust_del_bps{$bp1}});
        
    }
    if (($top_sbp2_num == $second_sbp2_num) and ($top_sbp2_num > 1)){
        @sbp2 = (@{$clust_bp2{$top_sbp2}}, @{$clust_bp2{$second_sbp2}});
    }
    elsif ($top_sbp2_num > 1){
        @sbp2 = (@{$clust_bp2{$top_sbp2}});
    }
    else{
        @sbp2 = (@{$clust_del_bps{$bp2}});
        
    }
    my $bp1_num = scalar @sbp1;
    my $bp1_h = int ($bp1_num / 2);
    @sbp1 = sort {$a <=> $b} @sbp1;
    $select_sbp1 = $sbp1[$bp1_h];
    my $bp2_num = scalar @sbp2;
    my $bp2_h = int ($bp2_num / 2);
    @sbp2 = sort {$a <=> $b} @sbp2;
    $select_sbp2 = $sbp2[$bp2_h];

    my $insbp_num = 0;
    for (my $i = $select_sbp1 - 1000; $i <= $select_sbp1 + $bp_diff; $i++){    # check the number of 5'-clipped reads within 1 Kb upstream of DEL BP1
        if (exists $bam_bp2{$i}){
            $insbp_num += @{$bam_bp2{$i}};
        }
    }
    next if ($insbp_num >= $read_num);
    
    my $delsize = $select_sbp2 - $select_sbp1 + 1;
    my $inslen = 0;
    if (exists ${$bam_indel{$bp1}}{$bp2}){
        my $sum = 0;
        map{$sum += $_} @{${$bam_indel{$bp1}}{$bp2}};
        $inslen = int ($sum / @{${$bam_indel{$bp1}}{$bp2}});
    }
    $bam_del_bp{$select_sbp1} = "$delsize=$read_num=$read_num=$inslen";
    $DEL{$select_sbp1} = 1;
    $used{$select_sbp1} = 1;
}

foreach my $bp1 (sort {$a <=> $b} keys %bam_del_bp){        # remove BP-DELs overlapping assigned DELs
    my ($blen) = split (/=/, $bam_del_bp{$bp1});
    foreach my $dpos (sort {$a <=> $b} keys %bam_del){
        last if ($dpos > $bp1 + 20);
        next if ($dpos < $bp1 - 500);
        my $dlen = $bam_del_len{$dpos};
        if (($blen / $dlen >= 0.8) and ($blen / $dlen <= 1.2)){
            delete $bam_del_bp{$bp1};
            last;
        }
    }
}

$pre_bp1 = 0;
$pre_num = 0;
foreach my $bp1 (sort {$a <=> $b} keys %rep_bp){            # merge REP BP1 clusters
    if (($pre_bp1 > 0) and ($bp1 - $pre_bp1 <= $bp_diff)){
        my $num = 0;
        foreach my $bp2 (keys %{$rep_bp{$bp1}}){
            $num += ${$rep_bp{$bp1}}{$bp2};
        }
        if ($pre_num == 0){
            foreach my $bp2 (keys %{$rep_bp{$pre_bp1}}){
                $pre_num += ${$rep_bp{$pre_bp1}}{$bp2};
            }
        }
        if ($pre_num > $num){
            foreach my $bp2 (keys %{$rep_bp{$bp1}}){
                ${$rep_bp{$pre_bp1}}{$bp2} += ${$rep_bp{$bp1}}{$bp2};
                push @{${$rep_bp_info{$pre_bp1}}{$bp2}}, @{${$rep_bp_info{$bp1}}{$bp2}};
            }
            $pre_num += $num;
            push @{$clust_rep_bps{$pre_bp1}}, @{$clust_rep_bps{$bp1}} if (@{$clust_rep_bps{$bp1}} < 10);
            delete $rep_bp{$bp1};
            next;
        }
        else{
            foreach my $bp2 (keys %{$rep_bp{$pre_bp1}}){
                ${$rep_bp{$bp1}}{$bp2} += ${$rep_bp{$pre_bp1}}{$bp2};
                push @{${$rep_bp_info{$bp1}}{$bp2}}, @{${$rep_bp_info{$pre_bp1}}{$bp2}};
            }
            $pre_num += $num;
            push @{$clust_rep_bps{$bp1}}, @{$clust_rep_bps{$pre_bp1}} if (@{$clust_rep_bps{$pre_bp1}} < 10);
            delete $rep_bp{$pre_bp1};
        }
    }
    else{
        $pre_num = 0;
    }
    $pre_bp1 = $bp1;
}
 
foreach my $bp1 (sort {$a <=> $b} keys %rep_bp){            # merge REP BP2 clusters
    my $pre_bp2 = 0;
    my $pre_len = 0;
    foreach my $bp2 (sort {$a <=> $b} keys %{$rep_bp{$bp1}}){
        my $len = $bp2 - $bp1 + 1;
        my $rep_bp_diff = $bp_diff;
        $rep_bp_diff = 50 if ($len < 1000);
        $rep_bp_diff = 150 if ($len > 10000);
        if ($pre_bp2 > 0){
            my $diff = $bp2 - $pre_bp2 + 1;
            if ($diff <= $rep_bp_diff){
                ${$rep_bp{$bp1}}{$bp2} += ${$rep_bp{$bp1}}{$pre_bp2};
                push @{$clust_rep_bps{$bp2}}, @{$clust_rep_bps{$pre_bp2}} if (@{$clust_rep_bps{$pre_bp2}} < 10);
                delete ${$rep_bp{$bp1}}{$pre_bp2};
            }
        }
        $pre_bp2 = $bp2;
        $pre_len = $len;
    }
}

foreach my $bp1 (sort {$a <=> $b} keys %rep_bp){            # delete REP with shared BP1
    if (scalar keys %{$rep_bp{$bp1}} == 0){
        delete $rep_bp{$bp1};
        next;
    }
    if (scalar keys %{$rep_bp{$bp1}} > 1){
        my $top_num = 0;
        my $sec_num = 0;
        my $top_bp2 = 0;
        my $sec_bp2 = 0;
        foreach my $bp2 (sort {${$rep_bp{$bp1}}{$b} <=> ${$rep_bp{$bp1}}{$a}} keys %{$rep_bp{$bp1}}){
            if ($top_num == 0){
                $top_num = ${$rep_bp{$bp1}}{$bp2};
                $top_bp2 = $bp2;
                next;
            }
            elsif ($sec_num == 0){
                $sec_num = ${$rep_bp{$bp1}}{$bp2};
                $sec_bp2 = $bp2;
                last;
            }
        }
        if ($top_num > $sec_num){
            foreach my $bp2 (keys %{$rep_bp{$bp1}}){
                next if ($bp2 == $top_bp2);
                delete ${$rep_bp{$bp1}}{$bp2};
            }
        }
        else{
            foreach my $bp2 (keys %{$rep_bp{$bp1}}){
                my $num = ${$rep_bp{$bp1}}{$bp2};
                next if ($num >= $top_num);
                delete ${$rep_bp{$bp1}}{$bp2};
            }
        }
    }
}

foreach my $bp1 (sort {$a <=> $b} keys %rep_bp){
    my $select_bp2 = 0;
    my $pre_num = 0;
    my $pre_len = 0;
    if (scalar keys %{$rep_bp{$bp1}} == 0){
        delete $rep_bp{$bp1};
        next;
    }
    foreach my $bp2 (sort {${$rep_bp{$bp1}}{$b} <=> ${$rep_bp{$bp1}}{$a}} keys %{$rep_bp{$bp1}}){
        my $len = $bp2 - $bp1 + 1;
        my $num = ${$rep_bp{$bp1}}{$bp2};
        last if ($pre_num > $num);
        if ($pre_num == $num){
            next if ($pre_len < $len);
        }
        $select_bp2 = $bp2;
        $pre_num = $num;
        $pre_len = $len;
    }
    if ($select_bp2 != 0){
        my $read_num = ${$rep_bp{$bp1}}{$select_bp2};
        my $bp1_med = 0;
        my $bp2_med = 0;
        my $bp1_num = scalar @{$clust_rep_bps{$bp1}};
        my $bp2_num = scalar @{$clust_rep_bps{$select_bp2}};
        my $bp1_hn = int ($bp1_num / 2);
        my $bp2_hn = int ($bp2_num / 2);
        my @rep_bp1 = sort {$a <=> $b} @{$clust_rep_bps{$bp1}};
        my @rep_bp2 = sort {$a <=> $b} @{$clust_rep_bps{$select_bp2}};
        $bp1_med = $rep_bp1[$bp1_hn];
        $bp2_med = $rep_bp2[$bp2_hn];
        my $repsize = $bp2_med - $bp1_med + 1;
        if (!exists ${$rep_bp_info{$bp1}}{$select_bp2}){
        }
        my $info = '';
        $info = join ('|', @{${$rep_bp_info{$bp1}}{$select_bp2}}) if (exists ${$rep_bp_info{$bp1}}{$select_bp2});
        $bam_rep_bp{$bp1_med} = "$repsize=$read_num=$read_num=$info";
        $DEL{$bp1_med} = 1;
        $used{$bp1_med} = 1;
    }
}

# INS BP finding
my %ins_bp1;

foreach my $bp1 (sort {$a <=> $b} keys %ins_bp_inread){
    my $delete_flag1 = 0;
    foreach (@{$clust_ins_bps{$bp1}}){
        if ((exists $DEL{$_}) or (exists $inv_bps{$_})){
            $delete_flag1 = 1;
            last;
        }
    }
    next if ($delete_flag1 == 1);
    my %idist;
    foreach my $bp2 (sort {$a <=> $b} keys %bam_ins){
        last if ($bp2 > $bp1 + 50);
        next if ($bp2 < $bp1 - 50);
        my $dist = abs ($bp1 - $bp2);
        $idist{$dist} = $bp2;
    }
    if (scalar keys %idist > 0){
    }
    else{
        my $sum_bpins = 0;
        my $bp_info = '';
        my %bpins_id;
        my $bpnum = scalar @{$ins_bp_inread{$bp1}};
        next if ($bpnum < $min_ins_reads);
        foreach (@{$ins_bp_inread{$bp1}}){
            my ($id1, $strand1, $ilen1, $rpos1) = split (/=/, $_);
            $sum_bpins += $ilen1;
            $bp_info .= "$id1=$strand1=$ilen1=$rpos1=$bp1|";
            $bpins_id{$id1} = 1;
        }
        my $bpins_num = 0;
        my $ins_num = 0;
        $bpins_num = scalar keys %bpins_id;
        my $ave_bplen = int ($sum_bpins / $bpins_num + 0.5);
        $bp_info =~ s/\|$//;
        $bam_ins_bp{$bp1} = "$bp_info~$bpins_num~$ave_bplen";
    }
}

my $pre_bp = 0;
foreach my $bp1 (sort {$a <=> $b} keys %bam_ins_bp){
    if (($pre_bp > 0) and ($bp1 - $pre_bp <= $bp_diff3)){
        my ($pre_info, $pre_num, $pre_len) = split (/~/, $bam_ins_bp{$pre_bp});
        my ($info, $num, $len) = split (/~/, $bam_ins_bp{$bp1});
        if (($len / $pre_len > 0.5) or ($len / $pre_len < 2)){
            $info .= "|$pre_info";
            $num += $pre_num;
            my $sum_len = $len * $num + $pre_len * $pre_num;
            my $new_len = int ($sum_len / $num + 0.5);
            delete $bam_ins_bp{$pre_bp};
            $bam_ins_bp{$bp1} = "$info~$num~$new_len";
        }
        else{
            if ($num >= $pre_num){
                delete $bam_ins_bp{$pre_bp};
            }
            else{
                delete $bam_ins_bp{$bp1};
            }
        }
    }
    $pre_bp = $bp1;
}


# delete INS BP overlapping defined DEL-BP, DUP-BP. INV-BP or >2 double-clipped reads

my %assign_bp;
foreach my $bp1 (keys %bam_del_bp){
    my ($dlen, $bp1_num, $bp2_num) = split (/=/, $bam_del_bp{$bp1});
    my $read_num = $bp1_num;
    $read_num = $bp2_num if ($bp2_num > $bp1_num);
    next if ($read_num < $min_del_reads);
    $assign_bp{$bp1} = 1;
    my $bp2 = $bp1 + $dlen - 1;
    $assign_bp{$bp2} = 1;
}
foreach my $bp1 (keys %bam_dup){
    my $del_hit = 0;
    foreach my $dpos (sort {$a <=> $b} keys %bam_del){
        next if (!exists $bam_del_len{$dpos});
        last if ($dpos > $bp1 + 30);
        my $dlen = $bam_del_len{$dpos};
        next if ($dlen < 500);
        my $dend = $dpos + $dlen - 1;
        next if ($dend < $bp1 - 30);
        if (abs ($dend - $bp1) <= 30){
            $del_hit = 1;
            last;
        }
    }
    if ($del_hit == 1){
        delete $bam_dup{$bp1};
        next;
    }
    foreach my $dpos (sort {$a <=> $b} keys %bam_del_bp){
        last if ($dpos > $bp1 + 30);
        my ($dlen) = split (/=/, $bam_del_bp{$dpos});
        my $dend = $dpos + $dlen - 1;
        next if ($dend < $bp1 - 30);
        if (abs ($dend - $bp1) <= 30){
            $del_hit = 1;
            last;
        }
    }
    if ($del_hit == 1){
        delete $bam_dup{$bp1};
        next;
    }
    foreach my $bp2 (keys %{$bam_dup{$bp1}}){
        foreach my $dpos (sort {$a <=> $b} keys %bam_del){
            next if (!exists $bam_del_len{$dpos});
            last if ($dpos > $bp2 + 30);
            my $dlen = $bam_del_len{$dpos};
            next if ($dlen < 500);
            my $dend = $dpos + $dlen - 1;
            next if ($dpos < $bp2 - 30);
            if (abs ($dpos - $bp2) <= 30){
                $del_hit = 1;
                last;
            }
        }
        if ($del_hit == 1){
            delete ${$bam_dup{$bp1}}{$bp2};
            delete $bam_dup{$bp1} if (scalar keys %{$bam_dup{$bp1}} == 0);
            last;
        }
        foreach my $dpos (sort {$a <=> $b} keys %bam_del_bp){
            last if ($dpos > $bp2 + 30);
            my ($dlen) = split (/=/, $bam_del_bp{$dpos});
            my $dend = $dpos + $dlen - 1;
            next if ($dpos < $bp2 - 30);
            if (abs ($dpos - $bp2) <= 30){
                $del_hit = 1;
                last;
            }
        }
        if ($del_hit == 1){
            delete ${$bam_dup{$bp1}}{$bp2};
            delete $bam_dup{$bp1} if (scalar keys %{$bam_dup{$bp1}} == 0);
            last;
        }
        my @info = split (/\|/, ${$bam_dup{$bp1}}{$bp2});
        my $read_num = scalar @info;
        my $dp_rate = 0;
        my $flank_dp = 0;
        ($dp_rate, $flank_dp) = &calc_dprate ($bp1, $bp2, \%DP, 'DUP');
        next if ($dp_rate < 1.2);
        my $dup_rate = 0;
        $dup_rate = int ($read_num / $flank_dp * 100) / 100 if ($flank_dp > 0);
        if (($read_num >= $min_ins_reads) and ($dup_rate >= $min_VRR) and ($dp_rate >= 1.2)){
            $assign_bp{$bp1} = 1;
            $assign_bp{$bp2} = 1;
        }
    }
}
foreach my $bp1 (keys %bam_inv){
    foreach my $bp2 (keys %{$bam_inv{$bp1}}){
        my @info = split (/\|/, ${$bam_inv{$bp1}}{$bp2});
        my $read_num = scalar @info;
        my $dp = 0;
        ($dp) = &calc_dp ($bp1, $bp2, \%DP, 'INV');
        my $inv_rate = 0;
        $inv_rate = int ($read_num / $dp * 100) / 100 if ($dp > 0);
        if (($read_num >= $min_ins_reads) and ($inv_rate >= $min_VRR)){
            $assign_bp{$bp1} = 1;
            $assign_bp{$bp2} = 1;
        }
    }
}
foreach my $bp1 (keys %used_dupinv_bp){
    $assign_bp{$bp1} = 1;
    foreach my $bp2 (keys %{$used_dupinv_bp{$bp1}}){
        $assign_bp{$bp2} = 1;
    }
}
foreach my $bp (sort {$a <=> $b} keys %bam_ins_bp){
    my $sbp_match = 0;
    foreach my $sbp (sort {$a <=> $b} keys %assign_bp){
        last if ($sbp > $bp + $bp_diff);
        if (abs ($bp - $sbp) <= $bp_diff){
            $sbp_match = 1;
            last;
        }
    }
    if ($sbp_match == 1){
        delete $bam_ins_bp{$bp};
    }
    else{
        my $dclip_match = 0;
        foreach my $dbp1 (sort {$a <=> $b} keys %bam_2bp){
            next if ($bp < $dbp1 - $bp_diff);
            if (abs ($bp -$dbp1) <= $bp_diff){
                $dclip_match ++;
            }
            foreach (@{$bam_2bp{$dbp1}}){
                my ($dbp2) = split (/=/, $_);
                if (abs ($bp -$dbp2) <= $bp_diff){
                    $dclip_match ++;
                }
            }
        }
        if ($dclip_match >= 3){
            delete $bam_ins_bp{$bp};
        }
    }
}
foreach my $pos (sort {$a <=> $b} keys %bam_del){
    my $sum_len = 0;
    my $max_len = 0;
    my $read_num = scalar @{$bam_del{$pos}};
    foreach (@{$bam_del{$pos}}){
        my ($dlen) = split (/=/, $_);
        $sum_len += $dlen;
        $max_len = $dlen if ($max_len < $dlen);
    }
    my $ave_len = int ($sum_len / $read_num + 0.5);
    next if ($ave_len < 300);
    for (my $i = $pos - 20; $i <= $pos + 20; $i ++){
        $assign_bp{$i} = 1;
    }
    my $end = $pos + $max_len - 1;
    for (my $i = $end - 50; $i <= $end + 50; $i ++){
        $assign_bp{$i} = 1;
    }
}
foreach my $pos (sort {$a <=> $b} keys %bam_ins){
    my $ilen = $bam_ins_len{$pos};
    next if ($ilen < 200);
    for (my $i = $pos - 20; $i <= $pos + 20; $i ++){
        $assign_bp{$i} = 1;
    }
}

# collecting INS BPs from different reads

$pre_bp1 = 0;
foreach my $bp1 (sort {$a <=> $b} keys %bam_bp1){
    next if (exists $bam_str_bp1{$bp1});
    next if (exists $assign_bp{$bp1});
    if (($pre_bp1 > 0) and ($bp1 - $pre_bp1 <= 20)){
        push @{$bam_bp1{$bp1}}, @{$bam_bp1{$pre_bp1}};
        delete $bam_bp1{$pre_bp1};
    }
    $pre_bp1 = $bp1;
}
my $pre_bp2 = 0;
foreach my $bp2 (sort {$a <=> $b} keys %bam_bp2){
    next if (exists $bam_str_bp2{$bp2});
    next if (exists $assign_bp{$bp2});
    if (($pre_bp2 > 0) and ($bp2 - $pre_bp2 <= 20)){
        push @{$bam_bp2{$bp2}}, @{$bam_bp2{$pre_bp2}};
        delete $bam_bp2{$pre_bp2};
    }
    $pre_bp2 = $bp2;
}

my %ins_bp12;
my %ins_bp3;
my %ins_bp3_tag;
foreach my $bp1 (sort {$a <=> $b} keys %bam_bp1){
    my $assign_flag = 0;
    $assign_flag = 1 if (exists $bam_str_bp1{$bp1}) or (exists $assign_bp{$bp1});
    my $read_num1 = scalar @{$bam_bp1{$bp1}};
    next if ($read_num1 == 1);
    my $read_num2 = 0;
    my @bp2 = ();
    my $top_bp2 = 0;
    my $top_num2 = 0;
    foreach my $bp2 (sort {$a <=> $b} keys %bam_bp2){
        last if ($bp2 > $bp1 + 1000);
        next if ($bp2 < $bp1 - 10000);
        $assign_flag ++ if (exists $bam_str_bp2{$bp2}) or (exists $assign_bp{$bp2});
        push @bp2, $bp2;
        my $num = scalar @{$bam_bp2{$bp2}};
        if ($num > $top_num2){
            $top_num2 = $num;
            $top_bp2 = $bp2;
        }
    }
    if ($top_bp2 > 0){
        next if ($assign_flag > 0);
        $read_num2 = scalar @{$bam_bp2{$top_bp2}};
        next if ($read_num2 <= 1) or ($read_num1 + $read_num2 <= 3);
        if ($top_bp2 <= $bp1){
            if (exists $ins_bp12{$top_bp2}){
                my $pre_rnum1 = scalar @{${$ins_bp12{$top_bp2}}{1}};
                if ($read_num1 > $pre_rnum1){
                    @{${$ins_bp12{$top_bp2}}{1}} = (@{$bam_bp1{$bp1}});
                    @{${$ins_bp12{$top_bp2}}{2}} = (@{$bam_bp2{$top_bp2}});
                    $ins_bp3{$top_bp2} = $bp1;
                    $ins_bp3_tag{$top_bp2} = '5T';
                }
            }
            else{
                push @{${$ins_bp12{$top_bp2}}{1}}, @{$bam_bp1{$bp1}};
                push @{${$ins_bp12{$top_bp2}}{2}}, @{$bam_bp2{$top_bp2}};
                $ins_bp3{$top_bp2} = $bp1;
                $ins_bp3_tag{$top_bp2} = '5T';
            }
        }
        else{
            if (exists $ins_bp12{$bp1}){
                my $pre_rnum2 = scalar @{${$ins_bp12{$bp1}}{2}};
                if ($read_num2 > $pre_rnum2){
                    @{${$ins_bp12{$bp1}}{1}} = (@{$bam_bp1{$bp1}});
                    @{${$ins_bp12{$bp1}}{2}} = (@{$bam_bp2{$top_bp2}});
                    $ins_bp3{$bp1} = $top_bp2;
                    $ins_bp3_tag{$bp1} = '3T';
                }
            }
            else{
                push @{${$ins_bp12{$bp1}}{1}}, @{$bam_bp1{$bp1}};
                push @{${$ins_bp12{$bp1}}{2}}, @{$bam_bp2{$top_bp2}};
                $ins_bp3{$bp1} = $top_bp2;
                $ins_bp3_tag{$bp1} = '3T';
            }
        }
    }
}

$pre_bp1 = 0;
$pre_bp2 = 0;
$pre_num = 0;
foreach my $bp1 (sort {$a <=> $b} keys %ins_bp3){
    my $bp2 = $ins_bp3{$bp1};
    my $num = scalar @{${$ins_bp12{$bp1}}{1}};
    $num += scalar @{${$ins_bp12{$bp1}}{2}};
    if (($pre_bp1 > 0) and ($pre_bp2 >= $bp1)){
        if ($num >= $pre_num){
            delete $ins_bp12{$pre_bp1};
        }
        else{
            delete $ins_bp12{$bp1};
            next;
        }
    }
    $pre_bp1 = $bp1;
    $pre_bp2 = $bp2;
    $pre_num = $num;
}

my %bam_ins_bp2;

foreach my $pos (sort {$a <=> $b} keys %ins_bp12){      # assign INS to INS-BP with BPLEN or INS-BP to DUP
    my $read_num = 0;
    my $S5_clip = 0;
    my $S3_clip = 0;
    my @pos2;
    foreach my $pos2 (sort {$a <=> $b} keys %{$ins_bp12{$pos}}){
        my $clip_tag = '3T' if ($pos2 == 1);
        $clip_tag = '5T' if ($pos2 == 2);
        foreach my $info (@{${$ins_bp12{$pos}}{$pos2}}){
            my ($readid, $strand, $cliplen, $tag) = split (/=/, $info);
            $S5_clip ++ if ($clip_tag eq '5T') and ($tag eq 'S');
            $S3_clip ++ if ($clip_tag eq '3T') and ($tag eq 'S');
        }
    }
    if (($S5_clip <= 1) or ($S3_clip <= 1)){
        delete $ins_bp12{$pos};
        next;
    }
    if (($S5_clip > $S3_clip * 3)){
        delete $ins_bp12{$pos};
        next;
    }
    if ($S3_clip > $S5_clip * 3){
        delete $ins_bp12{$pos};
        next;
    }
    my $BP = int (($S5_clip + $S3_clip) * 0.5 + 0.5);
    my $pos2 = $ins_bp3{$pos};
    my $bplen = $pos2 - $pos + 1;
    if ($bplen > 100){
        my $dup_hit = 0;
        foreach my $dpos1 (sort {$a <=> $b} keys %bam_dup){     # match betweem the breakpoints of INS-BP and DUP
            last if ($dpos1 > $pos2 + 50);
            my $match_dlen = 0;
            foreach my $dpos2 (sort {$a <=> $b} keys %{$bam_dup{$dpos1}}){
                next if ($dpos2 < $pos - 50);
                my $duplen = $dpos2 - $dpos1 + 1;
                if ((abs ($pos - $dpos1) <= 50) and (abs ($pos2 - $dpos2) <= 50)){
                    $match_dlen = $duplen;
                    last;
                }
            }
            if ($match_dlen > 0){
                delete $ins_bp12{$pos};
                $dup_hit = 1;
                last;
            }
        }
        next if ($dup_hit == 1);

        my %ins_info;
        my $max_len = 0;
        foreach my $ipos (sort {$a <=> $b} keys %bam_ins){
            last if ($ipos > $pos2 + 20);
            next if ($ipos < $pos - 20);
            my $ilen = $bam_ins_len{$ipos} if (exists $bam_ins_len{$ipos});
            if (!exists $bam_ins_len{$ipos}){
                my @len;
                foreach (@{$bam_ins{$ipos}}){
                    my ($id, $strand, $ilen2) = split (/=/, $_);
                    push @len, $ilen2;
                }
                my $hn = int (@len * 0.5);
                @len = sort {$a <=> $b} @len;
                $ilen = $len[$hn];
            }
            if (($ilen >= 200) and ($ilen > $bplen * 0.5)){
                push @{$ins_info{$ipos}}, @{$bam_ins{$ipos}};
                if ($ilen > $max_len){
                    $max_len = $ilen;
                }
            }
        }
        if ($max_len > 0){
            my @select_info;
            my %ins_ids;
            my %ins_ids_dup;
            foreach my $ipos (sort {$a <=> $b} keys %ins_info){
                foreach (@{$ins_info{$ipos}}){
                    my ($id, $strand, $ilen, $rpos) = split (/=/, $_);
                    if (!exists $ins_ids{$id}){
                        $ins_ids{$id} = $_;
                    }
                    else{
                        my ($id2, $strand2, $ilen2, $rpos2, $ipos2) = split (/=/, $ins_ids{$id});
                        $ilen2 += $ilen;
                        my $new_inf = "$id2=$strand2=$ilen2=$rpos2=$ipos2";
                        $ins_ids{$id} = $new_inf;
                        $ins_ids_dup{$id} = 1;
                    }
                    my $bp_flag = 0;
                    my $read_start = $read_align5{$id};
                    my $read_end = $read_align3{$id};
                    if (($read_start > $pos + 50) and ($read_start < $pos2)){
                        $bp_flag = 1;
                    }
                    elsif (($read_end > $pos) and ($read_end < $pos2 - 50)){
                         $bp_flag = 1;
                    }
                    if ($bp_flag == 0){
                        foreach my $bpos (sort {$a <=> $b} keys %{$read_bp2{$id}}){
                            last if ($bpos > $pos2);
                            next if ($bpos < $pos - 50);
                            $bp_flag = 1;
                        }
                    }
                    if ($bp_flag == 0){
                        foreach my $bpos (sort {$a <=> $b} keys %{$read_bp1{$id}}){
                            last if ($bpos > $pos2 + 50);
                            next if ($bpos < $pos);
                            $bp_flag = 1;
                        }
                    }
                    if ($bp_flag == 0){
                        push @select_info, $_;
                    }
                }
            }
            if (@select_info == 0){
                my $dp_rate = 0;
                my $flank_dp = 0;
                ($dp_rate, $flank_dp) = &calc_dprate ($pos, $pos2, \%DP, 'DUP');
                if ($dp_rate >= 1.2){
                    my $dup_info = '';
                    foreach my $ipos (keys %ins_info){
                        foreach (@{$ins_info{$ipos}}){
                            my ($id, $strand, $ilen) = split (/=/, $_);
                            $dup_info .= "$id=$strand=1-$strand=1=0|";
                        }
                    }
                    $dup_info =~ s/\|$//;
                    ${$bam_dup{$pos}}{$pos2} = $dup_info;
                    delete $ins_bp12{$pos};
                }
                foreach my $ipos (keys %ins_info){
                    delete $bam_ins{$ipos};
                }
            }
            else{
                delete $bam_ins{$pos} if (exists $bam_ins{$pos});
                $bam_ins_bp2{$pos} = "$BP=$S5_clip=$S3_clip=$pos=$pos2";
                foreach my $ipos (keys %ins_info){
                    delete $bam_ins{$ipos};
                }
                my @select_len;
                my $max_len2 = 0;
                my $min_len2 = 1000000000;
                foreach (@select_info){
                    my ($id, $strand, $ilen) = split (/=/, $_);
                    if ($ilen > $max_len2){
                        $max_len2 = $ilen;
                    }
                    if ($ilen < $min_len2){
                        $min_len2 = $ilen;
                    }
                }
                if ($max_len2 / $min_len2 <= 1.5){
                    foreach (@select_info){
                        my ($id, $strand, $ilen) = split (/=/, $_);
                        if (exists $ins_ids_dup{$id}){
                            $_ = $ins_ids{$id};
                        }
                        push @select_len, $ilen;
                        push @{$bam_ins{$pos}}, $_;
                    }
                    @select_len = sort {$a <=> $b} @select_len;
                    my $hn = int (@select_len * 0.5);
                    my $med_ilen = $select_len[$hn];
                    $bam_ins_len{$pos} = $med_ilen;
                }
                else{
                    my @max;
                    my @min;
                    my @minlen2;
                    foreach (@select_info){
                        my ($id, $strand, $ilen) = split (/=/, $_);
                        my $diff1 = $max_len2 - $ilen;
                        my $diff2 = $ilen - $min_len2;
                        if ($diff1 <= $diff2){
                            push @max, $_;
                        }
                        else{
                            push @min, $_;
                            push @minlen2, $ilen;
                        }
                    }
                    foreach (@max){
                        my ($id, $strand, $ilen) = split (/=/, $_);
                        if (exists $ins_ids_dup{$id}){
                            $_ = $ins_ids{$id};
                        }
                        push @select_len, $ilen;
                        push @{$bam_ins{$pos}}, $_;
                    }
                    $bam_ins_len{$pos} = $max_len2;
                    my $min_hn = int (@minlen2 * 0.5);
                    @minlen2 = sort {$a <=> $b} @minlen2;
                    my $med_minlen = $minlen2[$min_hn];
                    if ($med_minlen >= $bplen * 0.5){
                        my $pos_a = $pos + 1;
                        $bam_ins_len{$pos_a} = $med_minlen;
                        foreach (@min){
                            my ($id, $strand, $ilen) = split (/=/, $_);
                            if (exists $ins_ids_dup{$id}){
                                $_ = $ins_ids{$id};
                            }
                            push @{$bam_ins{$pos_a}}, $_;
                        }
                        $bam_ins_bp2{$pos_a} = "$BP=$S5_clip=$S3_clip=$pos_a=$pos2";
                        $bam_ins_len{$pos_a} = $med_minlen;
                    }
                    else{
                        foreach my $ipos (keys %ins_info){
                            next if (!exists $bam_ins_len{$ipos});
                            my $ilen = $bam_ins_len{$ipos};
                            if (($ilen / $med_minlen >= 0.8) and ($ilen / $med_minlen <= 1.25)){
                                push @{$bam_ins{$ipos}}, @{$ins_info{$ipos}} if ($ipos != $pos);
                            }
                        }
                    }
                }
                delete $ins_bp12{$pos};
            }
        }
    }
}

my %bam_ins_BP;
foreach my $pos (sort {$a <=> $b} keys %bam_ins){           # assign >200 bp INS to the separated BPs of 5'- and 3'-ends
    my $read_num = scalar @{$bam_ins{$pos}};
    next if ($read_num == 0);
    next if (exists $bam_ins_bp2{$pos});
    my $ilen = $bam_ins_len{$pos} if (exists $bam_ins_len{$pos});
    if (!exists $bam_ins_len{$pos}){
        my @len;
        foreach (@{$bam_ins{$pos}}){
            my ($id, $strand, $ilen2) = split (/=/, $_);
            push @len, $ilen2;
        }
        my $hn = int (@len * 0.5);
        @len = sort {$a <=> $b} @len;
        $ilen = $len[$hn];
    }
    next if ($ilen < 200);
    my %bp1;
    my %bp2;
    foreach my $bp1 (sort {$a <=> $b} keys %bam_bp2){
        last if ($bp1 > $pos + $ilen);
        next if ($bp1 < $pos - $ilen);
        my $num = scalar @{$bam_bp2{$bp1}};
        $bp1{$bp1} = $num;
    }
    foreach my $bp2 (sort {$a <=> $b} keys %bam_bp1){
        last if ($bp2 > $pos + $ilen);
        next if ($bp2 < $pos - $ilen);
        my $num = scalar @{$bam_bp1{$bp2}};
        $bp2{$bp2} = $num;
    }
    my $pre_bp1 = 0;
    my $pre_bp2 = 0;
    foreach my $bp1 (sort {$a <=> $b} keys %bp1){
        if (($pre_bp1 > 0) and ($bp1 - $pre_bp1 <= 20)){
            $bp1{$pre_bp1} += $bp1{$bp1};
            delete $bp1{$bp1};
            next;
        }
        $pre_bp1 = $bp1;
    }
    foreach my $bp2 (sort {$a <=> $b} keys %bp2){
        if (($pre_bp2 > 0) and ($bp2 - $pre_bp2 <= 20)){
            $bp2{$pre_bp2} += $bp2{$bp2};
            delete $bp1{$bp2};
            next;
        }
        $pre_bp2 = $bp2;
    }
    next if (scalar keys %bp1 == 0) or (scalar keys %bp2 == 0);
    my $top_bp1 = 0;
    my $top_bp2 = 0;
    foreach my $bp1 (sort {$bp1{$b} <=> $bp1{$a}} keys %bp1){
        $top_bp1 = $bp1;
        last;
    }
    foreach my $bp2 (sort {$bp2{$b} <=> $bp2{$a}} keys %bp2){
        next if ($bp2 < $top_bp1 - 100);
        $top_bp2 = $bp2;
        last;
    }
    if (($top_bp1 > 0) and ($top_bp2 > 0)){
        if ((($top_bp1 < $pos - 50) and ($top_bp2 < $pos - 50)) or (($top_bp1 > $pos + 50) and ($top_bp2 > $pos + 50))){
            my $bp1_keys = scalar keys %bp1;
            my $bp2_keys = scalar keys %bp2;
            foreach my $bp1 (sort {$bp1{$b} <=> $bp1{$a}} keys %bp1){
                next if ($bp1 == $top_bp1) and ($bp1_keys > 1);
                $top_bp1 = $bp1;
                last;
            }
            foreach my $bp2 (sort {$bp2{$b} <=> $bp2{$a}} keys %bp2){
                next if ($bp2 < $top_bp1 - 100);
                next if ($bp2 == $top_bp2) and ($bp2_keys > 1);
                $top_bp2 = $bp2;
                last;
            }
        }
    }
    if (($top_bp2 - $top_bp1 > $ilen * 2) or ($top_bp1 == 0) or ($top_bp2 == 0)){
        $top_bp1 = 0;
        $top_bp2 = 0;
        foreach my $bp2 (sort {$bp2{$b} <=> $bp2{$a}} keys %bp2){
            $top_bp2 = $bp2;
            last;
        }
        foreach my $bp1 (sort {$bp1{$b} <=> $bp1{$a}} keys %bp1){
            next if ($bp1 > $top_bp2 + 100);
            $top_bp1 = $bp1;
            last;
        }
    }
    next if ($top_bp1 == 0) and ($top_bp2 == 0);
    my $bp1_read = 0;
    my $bp2_read = 0;
    $bp1_read = $bp1{$top_bp1} if (exists $bp1{$top_bp1});
    $bp2_read = $bp2{$top_bp2} if (exists $bp2{$top_bp2});
    if (($top_bp1 == 0) and (abs ($top_bp2 - $pos) <= 50) and ($bp2_read >= 2)){
        $top_bp1 = $pos;
    }
    elsif (($top_bp2 == 0) and (abs ($top_bp1 - $pos) <= 50) and ($bp1_read >= 2)){
        $top_bp1 = $pos;
    }
    my $bplen = 0;
    $bplen = $top_bp2 - $top_bp1 + 1 if ($top_bp1 > 0) and ($top_bp2 > 0);
    next if ($bplen > $ilen * 2);
    my $bp1_3T = 0;
    my $bp2_5T = 0;
    if ($bplen >= 100){
        foreach my $bp1 (sort {$a <=> $b} keys %bam_bp1){
            last if ($bp1 > $top_bp1 + 40);
            next if ($bp1 < $top_bp1 - 40);
            $bp1_3T ++;
        }
        foreach my $bp2 (sort {$a <=> $b} keys %bam_bp2){
            last if ($bp2 > $top_bp2 + 40);
            next if ($bp2 < $top_bp2 - 40);
            $bp2_5T ++;
        }
        next if ($bp1_3T >= $bp1_read * 0.1);
        next if ($bp2_5T >= $bp2_read * 0.1);
    }
    if (($top_bp1 > 0) and ($top_bp2 > 0)){
        next if ($top_bp1 < $pos - 50) and ($top_bp2 < $pos - 50);
        next if ($top_bp1 > $pos + 50) and ($top_bp2 > $pos + 50);
    }
    my $BP = int (($bp1_read + $bp2_read) * 0.5 + 0.5);
    my $sum_clip = 0;
    foreach my $bp1 (keys %bp1){
        $sum_clip += $bp1{$bp1};
    }
    foreach my $bp2 (keys %bp2){
        $sum_clip += $bp2{$bp2};
    }
    next if ($sum_clip > ($bp1_read + $bp2_read) * 2);
    $bam_ins_bp2{$pos} = "$BP=$bp1_read=$bp2_read=$top_bp1=$top_bp2";
    if (($BP > 1) and ($bp1_read >= 2) and ($bp2_read >= 2) and ($bplen >= 100)){
        my $bp_range = "$top_bp1-$top_bp2";
        push @{$bam_ins_BP{$bp_range}}, "$pos=$ilen";
        $bam_ins_bp2{$top_bp1} = "$BP=$bp1_read=$bp2_read=$top_bp1=$top_bp2";
    }
}

foreach my $bp_range (keys %bam_ins_BP){
    my %max;
    my $max_len = 0;
    my $min_len = 100000000;
    my %read_ilen;
    my %read_ipos;
    my %select_read;
    my %bp_info;
    my ($bp1, $bp2) = split (/-/, $bp_range);
    foreach (@{$bam_ins_BP{$bp_range}}){
        my ($ipos, $ilen) = split (/=/, $_);
        next if (!exists $bam_ins{$ipos});
        foreach (@{$bam_ins{$ipos}}){
            my ($id, $strand, $ilen, $rpos) = split (/=/, $_);
            $read_ilen{$id} += $ilen;
            if (!exists $read_ipos{$id}){
                $read_ipos{$id} = $ipos;
            }
            else{
                $read_ipos{$id} = $ipos if ($ipos < $read_ipos{$id});
            }
        }
    }
    foreach my $id (keys %read_ilen){
        my $ilen = $read_ilen{$id};
        my $bp_flag = 0;
        my $read_start = $read_align5{$id};
        my $read_end = $read_align3{$id};
        if (($read_start > $bp1 + 50) and ($read_start < $bp2)){
            $bp_flag = 1;
        }
        elsif (($read_end > $bp1) and ($read_end < $bp2 - 50)){
            $bp_flag = 1;
        }
        if ($bp_flag == 0){
            foreach my $bpos (sort {$a <=> $b} keys %{$read_bp2{$id}}){
                last if ($bpos > $bp2);
                next if ($bpos < $bp1 - 50);
                $bp_flag = 1;
                if (($bpos >= $bp1 - 50) and ($bpos <= $bp1 + 50)){
                    my ($mlen, $cliplen, $cigar, $strand) = split (/=/, ${$read_bp2{$id}}{$bpos});
                    $bp_info{$id} = "5T=$cliplen=$strand";
                }
            }
        }
        if ($bp_flag == 0){
            foreach my $bpos (sort {$a <=> $b} keys %{$read_bp1{$id}}){
                last if ($bpos > $bp2 + 50);
                next if ($bpos < $bp1);
                $bp_flag = 1;
                if (($bpos >= $bp2 - 50) and ($bpos <= $bp2 + 50)){
                    my ($rpos, $cliplen, $cigar, $strand) = split (/=/, ${$read_bp1{$id}}{$bpos});
                    $bp_info{$id} = "3T=$cliplen=$strand";
                }
            }
        }
        if ($bp_flag == 0){
            $select_read{$id} = $ilen;
            if ($ilen > $max_len){
                $max_len = $ilen;
            }
            if ($ilen < $min_len){
                $min_len = $ilen;
            }
        }
    }
    if ($max_len == 0){
        my $dp_rate = 0;
        my $flank_dp = 0;
        ($dp_rate, $flank_dp) = &calc_dprate ($bp1, $bp2, \%DP, 'DUP');
        if ($dp_rate >= 1.2){
            my $dup_info = '';
            foreach my $id (keys %read_ipos){
                my $ipos = $read_ipos{$id};
                foreach (@{$bam_ins{$ipos}}){
                    my ($id, $strand, $ilen) = split (/=/, $_);
                    $dup_info .= "$id=$strand=1-$strand=1=0|";
                }
            }
            $dup_info =~ s/\|$//;
            ${$bam_dup{$bp1}}{$bp2} = $dup_info;
        }
        else{
            my @info_5T;
            my @info_3T;
            foreach my $id (keys %read_ipos){
                if (exists $bp_info{$id}){
                    my ($tag, $cliplen, $strand) = split (/=/, $bp_info{$id});
                    push @info_5T, "$id=$strand=$cliplen=S" if ($tag eq '5T');
                    push @info_3T, "$id=$strand=$cliplen=S" if ($tag eq '3T');
                }
            }
            next if ((@info_5T == 0) or (@info_3T == 0));
            push @{${$ins_bp12{$bp1}}{1}}, @info_3T;
            push @{${$ins_bp12{$bp1}}{2}}, @info_5T;
            $ins_bp3{$bp1} = $bp2;
            $ins_bp3_tag{$bp1} = '5T';
        }
        foreach (@{$bam_ins_BP{$bp_range}}){
            my ($ipos, $ilen) = split (/=/, $_);
            delete $bam_ins{$ipos};
        }
    }
    else{
        my $min_pos = 0;
        my %ipos;
        my %ipos2;
        my $bplen = $bp2 - $bp1 + 1;
        if ($max_len / $min_len >= 1.5){
            my @min_len;
            foreach my $rid (keys %select_read){
                my $ilen = $select_read{$rid};
                my $diff1 = $max_len - $ilen;
                my $diff2 = $ilen - $min_len;
                my $ipos = $read_ipos{$rid};
                if ($diff1 <= $diff2){
                    $ipos{$ipos} = 1;
                }
                else{
                    $ipos2{$ipos} = 1;
                    push @min_len, $ilen;
                }
            }
            my $hn = int (@min_len * 0.5);
            @min_len = sort {$a <=> $b} @min_len;
            my $med_minlen = $min_len[$hn];
            if ($med_minlen >= $bplen * 0.5){
                $min_pos = $bp1 + 1;
                $min_pos ++ if (exists $bam_ins{$min_pos});
                $min_pos ++ if (exists $bam_ins{$min_pos});
                delete $bam_ins{$min_pos};
                my ($BP, $bp1_read, $bp2_read, $top_bp1, $top_bp2) = split (/=/, $bam_ins_bp2{$bp1});
                $bam_ins_bp2{$min_pos} = "$BP=$bp1_read=$bp2_read=$min_pos=$top_bp2";
                foreach my $ipos (keys %ipos2){
                    next if (!exists $bam_ins{$ipos});
                    push @{$bam_ins{$min_pos}}, @{$bam_ins{$ipos}};
                    delete $bam_ins{$ipos};
                }
                $bam_ins_len{$min_pos} = $med_minlen;
                %ipos2 = ();
            }
        }
        else{
            foreach my $rid (keys %select_read){
                my $ilen = $select_read{$rid};
                my $ipos = $read_ipos{$rid};
                push @{$ipos{$ipos}}, @{$bam_ins{$ipos}};
            }
        }
        my @info;
        foreach my $ipos (keys %ipos){
            next if (!exists $bam_ins{$ipos});
            push @info, @{$bam_ins{$ipos}};
        }
        delete $bam_ins{$bp1} if (exists $bam_ins{$bp1});
        push @{$bam_ins{$bp1}}, @info;
        foreach (@{$bam_ins_BP{$bp_range}}){
            my ($ipos, $ilen) = split (/=/, $_);
            next if ($ipos == $bp1);
            next if ($ipos == $min_pos);
            next if (exists $ipos2{$ipos});
            delete $bam_ins{$ipos};
        }
        $bam_ins_len{$bp1} = $max_len;
    }
}
%bam_ins_BP = ();

my %bam_dup_ins;
foreach my $pos1 (sort {$a <=> $b} keys %bam_dup){              # assign INS within DUP to DUP
    foreach my $pos2 (sort {$a <=> $b} keys %{$bam_dup{$pos1}}){
        my $duplen = $pos2 - $pos1 + 1;
        my $end = $pos1 + $duplen - 1;
        my @info = split (/\|/, ${$bam_dup{$pos1}}{$pos2});
        my $orig_read = scalar @info;
        my $read_num = 0;
        my $bp5 = 0;
        my $bp3 = 0;
        my $filt_flag = 0;
        for (my $i = $pos1 - 20; $i <= $pos1 + 20; $i++){
            if (exists $bam_bp2{$i}){
                $bp5 += scalar @{$bam_bp2{$i}};
            }
        }
        for (my $i = $end - 50; $i <= $end + 50; $i++){
            if (exists $bam_bp1{$i}){
                $bp3 += scalar @{$bam_bp1{$i}};
            }
        }
        if (@info < 5){
            $filt_flag = 1 if ($bp3 > $bp5 * 3) or ($bp5 > $bp3 * 3);
        }
        $read_num = int (($bp5 + $bp3) * 0.5 + 0.5);
        $read_num = scalar @info if ($read_num < @info);
        $filt_flag = 1 if ($read_num < $min_ins_reads);
        $filt_flag = 1 if ($read_num < $min_ins_reads + 1) and ($duplen >= 1000000);
        $filt_flag = 1 if (($orig_read / $read_num <= 0.1) or ($orig_read == 2)) and ($duplen >= 1000000);
        if ($filt_flag == 1){
            delete ${$bam_dup{$pos1}}{$pos1};
            delete $bam_dup{$pos1} if (scalar keys %{$bam_dup{$pos1}} == 0);
            next;
        }
        my $dp_rate = 0;
        my $flank_dp = 0;
        my $incons_rate = 0;
        ($dp_rate, $flank_dp, $incons_rate) = &calc_dprate ($pos1, $end, \%DP, 'DUP');
        $filt_flag = 1 if ($flank_dp <= $ave_depth * 0.2);
        $filt_flag = 1 if ($dp_rate < 1.1);
        $filt_flag = 1 if ($incons_rate >= 0.5);
        if ($filt_flag == 1){
            delete ${$bam_dup{$pos1}}{$pos1};
            delete $bam_dup{$pos1} if (scalar keys %{$bam_dup{$pos1}} == 0);
            next;
        }
        next if ($duplen > 20000);
        my $gt = 'HT';
        $gt = 'HM' if ($read_num / $flank_dp > 0.8);
        my %read_len;
        my %read_pos;
        my %ins_pos;
        for (my $i = $pos1 - 50; $i <= $pos1 + 50; $i++){
            if ((exists $bam_ins_bp2{$i}) and (exists $bam_ins{$i})){
                my ($BP, $bp1_read, $bp2_read, $ins_bp1, $ins_bp2) = split (/=/, $bam_ins_bp2{$i});
                if (abs ($ins_bp2 - $pos2) <= 50){
                    foreach (@{$bam_ins{$i}}){
                        my ($id, $strand, $ilen) = split (/=/, $_);
                        $read_len{$id} += $ilen;
                        $read_pos{$id} = $i if (!exists $read_pos{$id});
                        $ins_pos{$i} = $id;
                    }
                }
            }
        }
        my $bp_flag = 0;
        my %dup_ins;
        if (scalar keys %read_pos > 0){
            foreach my $id (keys %read_pos){
                my $ilen = $read_len{$id};
                foreach my $bpos (sort {$a <=> $b} keys %{$read_bp2{$id}}){
                    last if ($bpos > $pos1 + 50);
                    next if ($bpos < $pos1 - 50);
                    $bp_flag = 1;
                }
                if ($bp_flag == 0){
                    foreach my $bpos (sort {$a <=> $b} keys %{$read_bp1{$id}}){
                        last if ($bpos > $pos2 + 50);
                        next if ($bpos < $pos2 - 50);
                        $bp_flag = 1;
                    }
                }
                if ($bp_flag == 0){
                    $dup_ins{$id} = $ilen;
                }
            }
            my $ins_read = scalar keys %dup_ins;
            my $max_inslen = 0;
            foreach my $id (sort {$dup_ins{$b} <=> $dup_ins{$a}} keys %dup_ins){
                $max_inslen = $dup_ins{$id};
                last;
            }
            $bam_dup_ins{$pos1} = "$ins_read=$max_inslen" if ($ins_read > 0);
            foreach my $ipos (keys %ins_pos){
                delete $bam_ins{$ipos};
            }
        }
        %read_len = ();
        %read_pos = ();
        %ins_pos = ();
        $bp_flag = 0;
        %dup_ins = ();
        if ($gt eq 'HM'){
            foreach my $ipos (sort {$a <=> $b} keys %bam_ins){
                last if ($ipos > $pos2 + 50);
                next if ($ipos < $pos1 - 50);
                if (exists $bam_ins{$ipos}){
                    foreach (@{$bam_ins{$ipos}}){
                        my ($id, $strand, $ilen) = split (/=/, $_);
                        $read_len{$id} += $ilen;
                        $read_pos{$id} = $ipos if (!exists $read_pos{$id});
                        $ins_pos{$ipos} = $id;
                    }
                }
            }
            if (scalar keys %read_pos > 0){
                foreach my $id (keys %read_pos){
                    my $ilen = $read_len{$id};
                    foreach my $bpos (sort {$a <=> $b} keys %{$read_bp2{$id}}){
                        last if ($bpos > $pos1 + 50);
                        next if ($bpos < $pos1 - 50);
                        $bp_flag = 1;
                    }
                    if ($bp_flag == 0){
                        foreach my $bpos (sort {$a <=> $b} keys %{$read_bp1{$id}}){
                            last if ($bpos > $pos2 + 50);
                            next if ($bpos < $pos2 - 50);
                            $bp_flag = 1;
                        }
                    }
                    if ($bp_flag == 0){
                        $dup_ins{$id} = $ilen;
                    }
                }
                my $ins_read = scalar keys %dup_ins;
                my $max_inslen = 0;
                foreach my $id (sort {$dup_ins{$b} <=> $dup_ins{$a}} keys %dup_ins){
                    $max_inslen = $dup_ins{$id};
                    last;
                }
                $bam_dup_ins{$pos1} = "$ins_read=$max_inslen" if ($ins_read > 0);
                foreach my $ipos (keys %ins_pos){
                    delete $bam_ins{$ipos};
                }
            }
        }
    }
}

my $out_chr = "$temp_dir/$out_prefix.chr$chr.discov.txt";
$out_chr = "$temp_dir/$out_prefix.$chr.discov.txt" if ($chr !~ /^[\dXY]+$/);
if (scalar keys %align_reads == 0){
    system ("rm $out_chr") if (-f $out_chr);
    die "Total alignment: 0 in $chr:\n";
}
open (OUT, "> $out_chr");
foreach my $pos (sort {$a <=> $b} keys %bam_ins){
    my $sum_len = 0;
    my $read_info = '';
    my %clust_id;
    my @len;
    my $read_num = scalar @{$bam_ins{$pos}};
    next if ($read_num == 0);
    foreach (@{$bam_ins{$pos}}){
        my ($id, $strand, $len, $rpos) = split (/=/, $_);
        $read_info .= "$_|";
        $clust_id{$id} = 1;
        if ($strand ne 'BP'){
            push @len, $len;
        }
    }
    map{$sum_len += $_} @len;
    my $BP = 0;
    my $pos1 = 0;
    my $pos2 = 0;
    my $bp1 = 0;
    my $bp2 = 0;
    my $bp1_3T = 0;
    my $bp2_5T = 0;
    my $sum_clip = 0;
    ($BP, $bp1, $bp2, $pos1, $pos2) = split (/=/, $bam_ins_bp2{$pos}) if (exists $bam_ins_bp2{$pos});
    $BP = 0 if (($bp1 == 0) or ($bp2 == 0)) and ($bp1 <= 1) and ($bp2 <= 1);
    my $bplen = 0;
    $bplen = $pos2 - $pos + 1 if ($pos2 > 0);
    my $read_num2 = scalar keys %clust_id;
    $read_num2 += $BP;
    my $ave_len = int ($sum_len / @len + 0.5);
    next if ($read_num2 < $min_ins_reads) and ($ave_len < 1000);
    $read_info =~ s/\|$//;
    my $hnum = int (@len * 0.5);
    @len = sort {$a <=> $b} @len;
    my $med_len = $len[$hnum];
    my $dp = 0;
    my $dp_rate = 0;   
    if ($bplen >= 200){
        my $end = $pos + $bplen - 1;
        ($dp_rate, $dp) = &calc_dprate ($pos, $end, \%DP, 'DUP');
    }
    else{
        ($dp) = &calc_dp ($pos, $pos, \%DP, 'INS');
    }
    my $ins_rate = 1;
    $ins_rate = int ($read_num2 / $dp * 100) / 100 if ($dp > 0);
    $ins_rate = 1 if ($ins_rate > 1);
    next if ($read_num2 < $min_ins_reads);
    next if ($read_num2 < $min_ins_reads + 1) and ($med_len <= 10);
    next if ($read_num2 < $min_ins_reads + 2) and ($med_len <= 5);
    next if ($ins_rate < $min_VRR);
    next if ($ins_rate < $min_VRR * 1.5) and ($med_len <= 10);
    next if ($ins_rate < $min_VRR * 2) and ($med_len <= 5);
    print OUT "$chr\t$pos\tINS\t$med_len\t$read_info\t$read_num2,$dp,$ins_rate,$BP,$bplen,$dp_rate\n";
}
foreach my $pos (sort {$a <=> $b} keys %bam_ins_bp){    # INSs supported by reads with both 5'-clipped and 3'-clipped alignments at the BPs
    my ($read_info1, $bpins_num, $inslen) = split (/~/, $bam_ins_bp{$pos});
    my $dp = 0;
    ($dp) = &calc_dp ($pos, $pos, \%DP, 'INS');
    my $read_num = $bpins_num;
    my $ins_rate = 1;
    $ins_rate = int ($read_num / $dp * 100) / 100 if ($dp > 0);
    next if ($read_num < $min_ins_reads + 1);
    next if ($ins_rate < $min_VRR);
    print OUT "$chr\t$pos\tINS-BP\t$inslen\t$read_info1\t$read_num,$dp,$ins_rate,$read_num,0,0\n";
}
foreach my $pos (sort {$a <=> $b} keys %ins_bp12){      # INSs supported by reads with either 5'-clipped or 3'-clipped alignment at the BPs
    my $read_info = '';
    my $read_num = 0;
    my $S5_clip = 0;
    my $S3_clip = 0;
    foreach my $pos2 (sort {$a <=> $b} keys %{$ins_bp12{$pos}}){
        my $clip_tag = '3T' if ($pos2 == 1);
        $clip_tag = '5T' if ($pos2 == 2);
        foreach my $info (@{${$ins_bp12{$pos}}{$pos2}}){
            my ($readid, $strand, $cliplen, $tag) = split (/=/, $info);
            $S5_clip ++ if ($clip_tag eq '5T') and ($tag eq 'S');
            $S3_clip ++ if ($clip_tag eq '3T') and ($tag eq 'S');
            $read_info .= "$readid=$strand=$cliplen=$clip_tag=$tag|";
        }
    }
    next if ($S5_clip <= 1) or ($S3_clip <= 1);
    next if ($S5_clip > $S3_clip * 3);
    next if ($S3_clip > $S5_clip * 3);
    $read_num = int (($S5_clip + $S3_clip) * 0.5 + 0.5);
    $read_info =~ s/\|$//;
    my $flag = 0;
    for (my $i = $pos - 50; $i <= $pos + 50; $i++){
        if (exists $bam_ins{$i}){
            my @len;
            my $sumlen = 0;
            foreach (@{$bam_ins{$i}}){
                my ($id, $strand, $len) = split (/=/, $_);
                if ($strand ne 'BP'){
                    push @len, $len;
                }
            }
            if (@len > 0){
                map{$sumlen += $_} @len;
                my $avelen = int ($sumlen / @len);
                if ($avelen >= 100){
                    $flag = 1;
                    last;
                }
            }
        }
        if (exists $bam_del{$i}){
            my @len;
            my $sumlen = 0;
            foreach (@{$bam_del{$i}}){
                my ($dlen) = split (/=/, $_);
                push @len, $dlen;
            }
            if (@len > 0){
                map{$sumlen += $_} @len;
                my $avelen = int ($sumlen / @len);
                if ($avelen >= 200){
                    $flag = 1;
                    last;
                }
            }
        }
    }
    next if ($flag == 1);
    my @insQ0_info;
    foreach my $ipos (sort {$a <=> $b} keys %bam_ins_Q0){
        last if ($ipos > $pos + 20);
        next if ($ipos < $pos - 20);
        push @insQ0_info, @{$bam_ins_Q0{$ipos}};
    }
    my $read_num2 = $read_num;
    $read_num2 += scalar @insQ0_info;
    my $pos2 = $ins_bp3{$pos};
    my $pos3 = $pos2;
    $pos3 = $pos + 200 if ($pos2 < $pos + 100);
    my $bp_tag = $ins_bp3_tag{$pos};
    my $dprate = 0;
    my $flank_dp = 0;
    ($dprate, $flank_dp) = &calc_dprate ($pos, $pos3, \%DP, 'INS');
    my $ins_rate = 1;
    $ins_rate = int ($read_num2 / $flank_dp * 100) / 100 if ($flank_dp > 0);
    $ins_rate = 1 if ($ins_rate > 1);
    next if ($read_num2 < $min_ins_reads + 1);
    next if ($ins_rate < $min_VRR);
    my $ins_info = 'NA';
    my $bp_distance = $pos2 - $pos + 1;
    if (@insQ0_info > 0){
        $read_info = '';
        my @len;
        foreach (@insQ0_info){
            my ($id, $strand, $len, $rpos) = split (/=/, $_);
            my $bp_flag = 0;
            my $read_start = $read_align5{$id};
            my $read_end = $read_align3{$id};
            if (($read_start > $pos + 50) and ($read_start < $pos2)){
                $bp_flag = 1;
            }
            elsif (($read_end > $pos) and ($read_end < $pos2 - 50)){
                $bp_flag = 1;
            }
            if ($bp_flag == 0){
                foreach my $bpos (sort {$a <=> $b} keys %{$read_bp2{$id}}){
                    last if ($bpos > $pos2);
                    next if ($bpos < $pos - 50);
                    $bp_flag = 1;
                }
            }
            if ($bp_flag == 0){
                foreach my $bpos (sort {$a <=> $b} keys %{$read_bp1{$id}}){
                    last if ($bpos > $pos2 + 50);
                    next if ($bpos < $pos);
                    $bp_flag = 1;
                }
            }
            if ($bp_flag == 0){
                $read_info .= "$_|";
                push @len, $len;
            }
        }
        if (@len > 0){
            $read_info =~ s/\|$//;
            my $hnum = int (@len * 0.5);
            @len = sort {$a <=> $b} @len;
            my $med_len = $len[$hnum];
            print OUT "$chr\t$pos\tINS\t$med_len\t$read_info\t$read_num2,$flank_dp,$ins_rate,$read_num,$bp_distance,$dprate\n";
        }
        else{
            @insQ0_info = ();
        }
    }
    if (@insQ0_info == 0){
        if ($bp_distance >= 150){
            if ($bp_tag eq '5T'){
                next if ($dprate <= 1);
                $ins_info = "DUP-$bp_distance";
            }
            elsif ($bp_tag eq '3T'){
                next if ($dprate >= 0.9);
                $ins_info = "DEL-$bp_distance";
            }
        }
        print OUT "$chr\t$pos\tINS-BP2\t0\t$read_info\t$read_num,$dprate,$ins_rate,$S5_clip,$S3_clip\t$ins_info\n";
    }
}
my $max_reads = 5;
foreach my $pos (sort {$a <=> $b} keys %bam_del){
    my $sum_len = 0;
    my $max_len = 0;
    my $read_num = scalar @{$bam_del{$pos}};
    my $del_reads = '';
    my $read_count = 0;
    my @len;
    foreach (@{$bam_del{$pos}}){
        my ($dlen, $dpos, $id) = split (/=/, $_);
        push @len, $dlen;
        $max_len = $dlen if ($max_len < $dlen);
        $read_count ++;
        $del_reads .= "$id," if ($read_count <= $max_reads);
    }
    $del_reads =~ s/,$//;
    map{$sum_len += $_} @len;
    my $hnum = int (@len * 0.5);
    @len = sort {$a <=> $b} @len;
    my $med_len = $len[$hnum];
    my $ave_len = int ($sum_len / $read_num + 0.5);
    my $end = $pos + $max_len - 1;
    my $bp_len = 0;
    my $bp_num_str = '';
    my $bp_num = 0;
    my $ave_end = $pos + $ave_len;
    my $BP = 0;
    if ($max_len >= 500){
        my $bp5 = 0;
        my $bp3 = 0;
        my $end = $pos + $ave_len - 1;
        for (my $i = $pos - 20; $i <= $pos + 20; $i++){
            if (exists $bam_bp1{$i}){
                $bp3 += scalar @{$bam_bp1{$i}};
            }
        }
        for (my $i = $end - 50; $i <= $end + 50; $i++){
            if (exists $bam_bp2{$i}){
                $bp5 += scalar @{$bam_bp2{$i}};
            }
        }
        $BP = int (($bp5 + $bp3) * 0.5 + 0.5) if ($bp5 > 0) and ($bp3 > 0);
    }
    $read_num += $BP;
    next if ($read_num < $min_del_reads);
    my $dp_rate = 0;
    my $flank_dp = 0;
    ($dp_rate, $flank_dp) = &calc_dprate ($pos, $ave_end, \%DP, 'DEL') if ($ave_len >= 50);
    ($flank_dp) = &calc_dp ($pos, $ave_end, \%DP, 'DEL') if ($ave_len < 50);
    my $del_rate = 1;
    $del_rate = int ($read_num / $flank_dp * 100) / 100 if ($flank_dp > 0);
    $del_rate = 1 if ($del_rate > 1);
    next if ($read_num < $min_del_reads + 1) and ($ave_len <= 10);
    next if ($read_num < $min_del_reads + 2) and ($ave_len <= 5);
    next if ($del_rate < $min_VRR);
    next if ($del_rate < $min_VRR * 1.5) and ($ave_len <= 10);
    next if ($del_rate < $min_VRR * 2) and ($ave_len <= 5);
    print OUT "$chr\t$pos\tDEL\t$ave_len\t$read_num\t$del_reads\t$dp_rate,$flank_dp,$del_rate,$BP\n";
    $total_del_it ++;
}
foreach my $pos (sort {$a <=> $b} keys %bam_del_bp){
    my ($dlen, $bp1_num, $bp2_num, $inslen) = split (/=/, $bam_del_bp{$pos});
    next if ($bp1_num < $min_del_reads) or ($bp2_num < $min_del_reads);
    my $end = $pos + $dlen - 1;
    my $bp5 = 0;
    my $bp3 = 0;
    my $read_num = 0;
    for (my $i = $pos - 20; $i <= $pos + 20; $i++){
        if (exists $bam_bp1{$i}){
            $bp3 += scalar @{$bam_bp1{$i}};
        }
    }
    for (my $i = $end - 50; $i <= $end + 50; $i++){
        if (exists $bam_bp2{$i}){
            $bp5 += scalar @{$bam_bp2{$i}};
        }
    }
    next if ($bp3 > $bp5 * 3) and ($bp2_num < 4);
    next if ($bp5 > $bp3 * 3) and ($bp1_num < 4);
    $read_num = int (($bp5 + $bp3) * 0.5 + 0.5);
    if ($bp5 + $bp3 < $bp1_num + $bp2_num){
        $read_num = int (($bp1_num + $bp2_num) * 0.5 + 0.5);
    }
    if ($bp1_num < $bp3){
        $bp1_num = $bp3;
    }
    if ($bp2_num < $bp5){
        $bp2_num = $bp5;
    }
    next if ($read_num < $min_del_reads);
    next if ($read_num < $min_del_reads + 1) and ($dlen >= 1000000);
    next if (($bp1_num / $read_num <= 0.1) or ($bp1_num == 2)) and ($dlen >= 1000000);
    my $dp_rate = 0;
    my $flank_dp = 0;
    my $incons_rate = 0;
    ($dp_rate, $flank_dp, $incons_rate) = &calc_dprate ($pos, $end, \%DP, 'DEL');
    my $del_rate = 1;
    $del_rate = int ($read_num / $flank_dp * 100) / 100 if ($flank_dp > 0);
    $del_rate = 1 if ($del_rate > 1);
    next if ($flank_dp >= $ave_depth * 5);
    next if ($del_rate < $min_VRR);
    next if ($incons_rate >= 0.5);
    print OUT "$chr\t$pos\tDEL-BP\t$dlen\t$bp1_num,$bp2_num\t$inslen\t$dp_rate,$flank_dp,$del_rate\n"; 
}

foreach my $pos (sort {$a <=> $b} keys %bam_indel_str){
    my ($strid, $send, $unit_len, $read_cov, $type1, $type2, $len1, $len2, $num1, $num2, $insbp, $dpr) = split (/=/, $bam_indel_str{$pos});
    my $len = $len1;
    my $num = $num1;
    my $sec_strid = 'NA';
    $sec_strid = $overlap_str_opt{$strid} if (exists $overlap_str_opt{$strid});
    if (($type1 ne 'NA') and ($type2 ne 'NA')){
        if (($type1 eq 'DEL') and ($type2 eq 'INS')){
            my $read_ipos = 'NA=NA=NA';
            $read_ipos = ${$bam_ins_str_read{$pos}}{2} if (exists ${$bam_ins_str_read{$pos}}{2});
            print OUT "$chr\t$pos\tTR-del\t$len1\t$num1\t$read_cov=$unit_len\n";
            print OUT "$chr\t$pos\tTR-ins\t$len2\t$num2\t$read_cov=$unit_len=$read_ipos=$insbp=$dpr=$sec_strid\n";
        }
        elsif (($type1 eq 'INS') and ($type2 eq 'DEL')){
            my $read_ipos = 'NA=NA=NA';
            $read_ipos = ${$bam_ins_str_read{$pos}}{1} if (exists ${$bam_ins_str_read{$pos}}{1});
            print OUT "$chr\t$pos\tTR-ins\t$len1\t$num1\t$read_cov=$unit_len=$read_ipos=$insbp=$dpr=$sec_strid\n";
            print OUT "$chr\t$pos\tTR-del\t$len2\t$num2\t$read_cov=$unit_len\n";
        }
        elsif ($type1 eq 'DEL'){
            print OUT "$chr\t$pos\tTR-del1\t$len1\t$num1\t$read_cov=$unit_len\n" if ($num1 >= $min_str_reads);
            if ($num2 >= $min_del_reads){
                print OUT "$chr\t$pos\tTR-del2\t$len2\t$num2\t$read_cov=$unit_len\n";
            }
        }
        elsif ($type1 eq 'INS'){
            my $read_ipos1 = 'NA=NA=NA';
            my $read_ipos2 = 'NA=NA=NA';
            $read_ipos1 = ${$bam_ins_str_read{$pos}}{1} if (exists ${$bam_ins_str_read{$pos}}{1});
            $read_ipos2 = ${$bam_ins_str_read{$pos}}{2} if (exists ${$bam_ins_str_read{$pos}}{2});
            print OUT "$chr\t$pos\tTR-ins1\t$len1\t$num1\t$read_cov=$unit_len=$read_ipos1=$insbp=$dpr=$sec_strid\n" if ($num1 >= $min_str_reads);
            if ($num2 >= $min_del_reads){
                print OUT "$chr\t$pos\tTR-ins2\t$len2\t$num2\t$read_cov=$unit_len=$read_ipos2=$insbp=1=$sec_strid\n";
            }
        }
    }
    else{
        my $read_ipos = 'NA=NA=NA';
        $read_ipos = ${$bam_ins_str_read{$pos}}{1} if (exists ${$bam_ins_str_read{$pos}}{1});
        print OUT "$chr\t$pos\tTR-del\t$len\t$num\t$read_cov=$unit_len\n" if ($type1 eq 'DEL');
        print OUT "$chr\t$pos\tTR-ins\t$len\t$num\t$read_cov=$unit_len=$read_ipos=$insbp=$dpr=$sec_strid\n" if ($type1 eq 'INS');
    }
}

foreach my $pos (sort {$a <=> $b} keys %bam_rep){
    my $sum_len = 0;
    my $read_num = scalar @{$bam_rep{$pos}};
    foreach my $dlen (@{$bam_rep{$pos}}){
        $sum_len += $dlen;
    }
    my $ave_len = int ($sum_len / $read_num + 0.5);
    my $end = $pos + $ave_len - 1;
    my $bp_len = 0;
    my $bp_num_str = '';
    my $bp_num = 0;
    my $bp_num2 = 0;
    foreach my $bpos (sort {$a <=> $b} keys %bam_rep_bp){
        last if ($bpos > $end);
        my ($rlen, $bp1_num, $bp2_num) = split (/=/, $bam_rep_bp{$bpos});
        $bp_num = $bp1_num;
        $bp_num = $bp2_num if ($bp1_num > $bp2_num);
        $bp_num2 += $bp1_num + $bp2_num;
        my $bend = $bpos + $rlen - 1;
        next if ($bend < $pos);
        next if ($bpos < $pos - $bp_diff);
        if ($bpos < $end){
            my $overlap = $end - $bpos;
            $overlap = $rlen if ($bend < $end);
            if (($overlap >= $ave_len * $min_overlap_rate) and ($overlap >= $rlen * $min_overlap_rate)){
                $bp_len = $rlen;
                $bp_num_str = "$bp1_num,$bp2_num";
                delete $bam_rep_bp{$bpos};
                last;
            }
        }
    }
    my $dp = 0;
    ($dp) = &calc_dp ($pos, $end, \%DP, 'REP');
    $read_num += $bp_num if ($bp_num_str ne '');
    my $rep_rate = 1;
    $rep_rate = int ($read_num / $dp * 100) / 100 if ($dp > 0);
    next if ($read_num < $min_ins_reads);
    next if ($read_num < $min_ins_reads + 1) and ($ave_len <= 10);
    next if ($read_num < $min_ins_reads + 2) and ($ave_len <= 5);
    next if ($rep_rate < $min_VRR);
    print OUT "$chr\t$pos\tREP\t$ave_len\t$read_num\tNA\t$dp,$rep_rate\n";
}
foreach my $pos (sort {$a <=> $b} keys %bam_rep_bp){
    my ($rlen, $bp1_num, $bp2_num, $info) = split (/=/, $bam_rep_bp{$pos});
    next if ($bp1_num > $bp2_num * 3);
    next if ($bp2_num > $bp1_num * 3);
    $info = 'NA' if (!defined $info);
    my $end = $pos + $rlen;
    my $dp = 0;
    ($dp) = &calc_dp ($pos, $end, \%DP, 'REP');
    my $read_num = $bp1_num;
    $read_num = $bp2_num if ($bp2_num > $bp1_num);
    my $rep_rate = 1;
    $rep_rate = int ($read_num / $dp * 100) / 100 if ($dp > 0);
    next if ($read_num < $min_ins_reads);
    next if ($rep_rate < $min_VRR);
    print OUT "$chr\t$pos\tREP-BP\t$rlen\t$bp1_num,$bp2_num\t$info\t$dp,$rep_rate\n";
}

foreach my $pos1 (sort {$a <=> $b} keys %bam_dup){
    foreach my $pos2 (sort {$a <=> $b} keys %{$bam_dup{$pos1}}){
        my $duplen = $pos2 - $pos1 + 1;
        my $end = $pos1 + $duplen - 1;
        my @info = split (/\|/, ${$bam_dup{$pos1}}{$pos2});
        my $orig_read = scalar @info;
        my $read_num = 0;
        my $bp5 = 0;
        my $bp3 = 0;
        for (my $i = $pos1 - 20; $i <= $pos1 + 20; $i++){
            if (exists $bam_bp2{$i}){
                $bp5 += scalar @{$bam_bp2{$i}};
            }
        }
        for (my $i = $end - 50; $i <= $end + 50; $i++){
            if (exists $bam_bp1{$i}){
                $bp3 += scalar @{$bam_bp1{$i}};
            }
        }
        if (@info < 5){
            next if ($bp3 > $bp5 * 3);
            next if ($bp5 > $bp3 * 3);
        }
        $read_num = int (($bp5 + $bp3) * 0.5 + 0.5);
        $read_num = scalar @info if ($read_num < @info);
        my $dup_ins = 0;
        my $dup_ins_len = 0;
        ($dup_ins, $dup_ins_len) = split (/=/, $bam_dup_ins{$pos1}) if (exists $bam_dup_ins{$pos1});
        $read_num += $dup_ins if ($dup_ins_len >= 200);
        $dup_ins_len = 0 if ($dup_ins_len < 200);
        next if ($read_num < $min_ins_reads);
        next if ($read_num < $min_ins_reads + 1) and ($duplen >= 1000000);
        next if (($orig_read / $read_num <= 0.1) or ($orig_read == 2)) and ($duplen >= 1000000);
        my $dp_rate = 0;
        my $flank_dp = 0;
        my $incons_rate = 0;
        ($dp_rate, $flank_dp, $incons_rate) = &calc_dprate ($pos1, $end, \%DP, 'DUP');
        next if ($flank_dp <= $ave_depth * 0.2);
        next if ($dp_rate < 1.1);
        next if ($incons_rate >= 0.5);
        my $dup_rate = 1;
        $dup_rate = int ($read_num / $flank_dp * 100) / 100 if ($flank_dp > 0);
        $dup_rate = 1 if ($dup_rate  > 1);
        next if ($dup_rate < $min_VRR);
        print OUT "$chr\t$pos1\tDUP\t$duplen\t${$bam_dup{$pos1}}{$pos2}\t$dp_rate,$flank_dp,$dup_rate,$read_num,$bp5,$bp3,$dup_ins,$dup_ins_len\n";
    }
}

foreach my $pos1 (sort {$a <=> $b} keys %bam_inv){
    foreach my $pos2 (sort {$a <=> $b} keys %{$bam_inv{$pos1}}){
        my $invlen = $pos2 - $pos1 + 1;
        my @info = split (/\|/, ${$bam_inv{$pos1}}{$pos2});
        my $read_num = scalar @info;
        next if ($read_num < $min_ins_reads);
        next if ($read_num < $min_ins_reads + 1) and ($invlen >= 1000000);
        my $dp = 0;
        ($dp) = &calc_dp ($pos1, $pos2, \%DP, 'INV');
        my $inv_rate = 1;
        $inv_rate = int ($read_num / $dp * 100) / 100 if ($dp > 0);
        next if ($inv_rate < $min_VRR);
        print OUT "$chr\t$pos1\tINV\t$invlen\t${$bam_inv{$pos1}}{$pos2}\t$dp,$inv_rate\n";
    }
}
close (OUT);

my $align_reads = scalar keys %align_reads;
my $rate_secalign = int ($sec_align_num / $align_reads * 1000) / 10;

my $chr3 = $chr;
$chr3 = 'chr' . $chr if ($chr !~ /^chr/i);
print STDERR "Total alignments in bam $chr3:    $total_align\n";
print STDERR "Total aligned reads in bam $chr3: $align_reads\n";
print STDERR "Total secondary alignments in $chr3: $sec_align_num ($rate_secalign%)\n";
print STDERR "Removed low MAPQ (< $min_mapQ) reads  $chr3: $delete_mapq\n";


sub calc_dp {
    my ($pos1, $pos2, $ref_DP, $type) = @_;
    my $flank_len = 500;
    my $ave_dp = 0;
    if ($type eq 'INS'){
        my $bp1 = int ($pos1 / 50) * 50;
        my @flank_dp;
        for (my $i = $bp1 - $flank_len; $i <= $bp1 + $flank_len; $i += 50){
            my $Mbin = int ($i / $Mbin_size);
            my $Mbin_res = $i % $Mbin_size;
            my $dp = 0;
            $dp = ${${$ref_DP}{$Mbin}}{$Mbin_res} if (exists ${${$ref_DP}{$Mbin}}{$Mbin_res});
            push @flank_dp, $dp if ($dp > 0);
        }
        my $sum = 0;
        map{$sum += $_} @flank_dp;
        $ave_dp = int ($sum / @flank_dp * 100 + 0.5) / 100 if (@flank_dp > 0);
    }
    else{
        my $bp1 = int ($pos1 / 50) * 50;
        my $bp2 = int ($pos2 / 50) * 50 + 100;
        my @flank_dp;
        my @flank_dp2;
        for (my $i = $bp1 - $flank_len; $i <= $bp1; $i += 50){
            my $Mbin = int ($i / $Mbin_size);
            my $Mbin_res = $i % $Mbin_size;
            my $dp = 0;
            $dp = ${${$ref_DP}{$Mbin}}{$Mbin_res} if (exists ${${$ref_DP}{$Mbin}}{$Mbin_res});
            push @flank_dp, $dp if ($dp > 0);
        }
        for (my $i = $bp2; $i <= $bp2 + $flank_len; $i += 50){
            my $Mbin = int ($i / $Mbin_size);
            my $Mbin_res = $i % $Mbin_size;
            my $dp = 0;
            $dp = ${${$ref_DP}{$Mbin}}{$Mbin_res} if (exists ${${$ref_DP}{$Mbin}}{$Mbin_res});
            if ($type eq 'STR'){
                push @flank_dp2, $dp if ($dp > 0);
            }
            else{
                push @flank_dp, $dp if ($dp > 0);
            }
        }
        if ($type eq 'STR'){
            my $sum1 = 0;
            my $sum2 = 0;
            my $ave_dp1 = 0;
            my $ave_dp2 = 0;
            map{$sum1 += $_} @flank_dp;
            map{$sum2 += $_} @flank_dp2;
            $ave_dp1 = int ($sum1 / @flank_dp * 100 + 0.5) / 100 if (@flank_dp > 0);
            $ave_dp2 = int ($sum2 / @flank_dp2 * 100 + 0.5) / 100 if (@flank_dp2 > 0);
            $ave_dp = $ave_dp1;
            $ave_dp = $ave_dp2 if ($ave_dp2 > $ave_dp1);
        }
        else{
            my $sum = 0;
            map{$sum += $_} @flank_dp;
            $ave_dp = int ($sum / @flank_dp * 100 + 0.5) / 100 if (@flank_dp > 0);
        }
    }
    return ($ave_dp);
}

sub calc_dprate {
    my ($pos1, $pos2, $ref_DP, $type) = @_;
    my $len = $pos2 - $pos1 + 1;
    my $flank_len = 500;
    my $bp_distance = $pos2 - $pos1 + 1;
    if ($bp_distance > 100000){
        $flank_len = 5000;
    }
    elsif ($bp_distance > 10000){
        $flank_len = int ($bp_distance * 0.05 / 50) * 50;
    }
    my $bp1_f2 = int (($pos1 - 100) / 50) * 50;
    my $bp1_f1 = $bp1_f2 - $flank_len;
    $bp1_f1 = 100 if ($bp1_f1 < 100);
    my $bp1 = int (($pos1 + 100) / 50) * 50 if ($len >= 100);
    if ($len < 100){
        $bp1 = int ($pos1 / 50) * 50;
        $bp1 += 50 if ($bp1 < $pos1);
    }
    my $bp2 = int (($pos2 - 100) / 50) * 50 if ($len >= 100);
    if ($len < 100){
        $bp2 = int ($pos2 / 50) * 50;
        $bp2 -= 50 if ($bp2 > $pos2);
    }
    $bp2 = $bp1 if ($bp1 > $bp2);
    my $bp2_f1 = int (($pos2 + 100) / 50) * 50;
    my $bp2_f2 = $bp1_f1 + $flank_len;
    my @flank_dp;
    my @dup_dp;
    for (my $i = $bp1_f1; $i <= $bp1_f2; $i += 50){
        my $Mbin = int ($i / $Mbin_size);
        my $Mbin_res = $i % $Mbin_size;
        my $dp = 0;
        $dp = ${${$ref_DP}{$Mbin}}{$Mbin_res} if (exists ${${$ref_DP}{$Mbin}}{$Mbin_res});
        push @flank_dp, $dp if ($dp > 0);
    }
    for (my $i = $bp2_f1; $i <= $bp2_f2; $i += 50){
        my $Mbin = int ($i / $Mbin_size);
        my $Mbin_res = $i % $Mbin_size;
        my $dp = 0;
        $dp = ${${$ref_DP}{$Mbin}}{$Mbin_res} if (exists ${${$ref_DP}{$Mbin}}{$Mbin_res});
        push @flank_dp, $dp if ($dp > 0);
    }
    my $step_size = 50;
    $step_size = 100 if ($len >= 10000);
    $step_size = 1000 if ($len >= 100000);
    for (my $i = $bp1; $i <= $bp2; $i += $step_size){
        my $Mbin = int ($i / $Mbin_size);
        my $Mbin_res = $i % $Mbin_size;
        my $dp = 0;
        $dp = ${${$ref_DP}{$Mbin}}{$Mbin_res} if (exists ${${$ref_DP}{$Mbin}}{$Mbin_res});
        push @dup_dp, $dp;
    }
    if (@flank_dp > 0){
        my $incons_num = 0;
        my @dup_dp2;
        my $sum_flank_dp = 0;
        my $sum_dp = 0;
        my $ave_flank_dp = 0;
        my $ave_dp = 0;
        map{$sum_flank_dp += $_} @flank_dp;
        $ave_flank_dp = int ($sum_flank_dp / @flank_dp * 100 + 0.5) / 100;
        foreach (@dup_dp){
            if ($type eq 'DEL'){
                if ($_ > $ave_flank_dp * 2){
                    $incons_num ++;
                }
                else{
                    push @dup_dp2, $_;
                }
            }
            elsif ($type eq 'DUP'){
                if ($_ < $ave_flank_dp * 0.5){
                    $incons_num ++;
                }
                else{
                    push @dup_dp2, $_;
                }
            }
            else{
                push @dup_dp2, $_;
            }
        }
        my $incons_rate = int ($incons_num / @dup_dp * 100 + 0.5) / 100;
        if ($incons_rate < 0.3){
            map{$sum_dp += $_} @dup_dp2;
            $ave_dp = int ($sum_dp / @dup_dp2 * 100 + 0.5) / 100 if (@dup_dp2 > 0);
        }
        else{
            map{$sum_dp += $_} @dup_dp;
            $ave_dp = int ($sum_dp / @dup_dp * 100 + 0.5) / 100 if (@dup_dp > 0);
        }
        my $dp_rate = int ($ave_dp / $ave_flank_dp * 100 + 0.5) / 100;
        return ($dp_rate, $ave_flank_dp, $incons_rate);
    }
    else{
        return (0, 0, 0);
    }
}

sub clust_indel{
    my ($ref_indel, $str) = @_;
    my @len = sort {$a <=> $b} @{$ref_indel};
    my $all_num = scalar @len;
    my %clust;
    my $clust_num = 1;
    my $pre_len = 0;
    foreach my $len (@len){
        if ($pre_len > 0){
            my $num = scalar @{$clust{$clust_num}};
            my $hnum = int ($num * 0.5);
            my $med = ${$clust{$clust_num}}[$hnum];
            if ($len / $med <= $str_max_len_rate){
                push @{$clust{$clust_num}}, $len;
            }
            else{
                $clust_num ++;
                push @{$clust{$clust_num}}, $len;
            }
        }
        else{
            push @{$clust{$clust_num}}, $len;
        }
        $pre_len = $len;
    }
    my $top_clust = 0;
    my $sec_clust = 0;
    my $thr_clust = 0;
    my $top_num = 0;
    my $sec_num = 0;
    my $thr_num = 0;
    my $sec_len = 0;
    my $thr_len = 0;
    foreach my $clust (sort {scalar @{$clust{$b}} <=> scalar @{$clust{$a}}} keys %clust){
        if ($top_clust == 0){
            $top_clust = $clust;
            $top_num = scalar @{$clust{$clust}};
        }
        elsif ($sec_clust == 0){
            $sec_clust = $clust;
            $sec_num = scalar @{$clust{$clust}};
            my $sum_len = 0;
            foreach (@{$clust{$clust}}){
                $sum_len += $_;
            }
            $sec_len = int ($sum_len / @{$clust{$clust}} + 0.5);
            next;
        }
        elsif ($thr_clust == 0){
            my $sum_len = 0;
            my $thr_num = scalar @{$clust{$clust}};
            foreach (@{$clust{$clust}}){
                $sum_len += $_;
            }
            $thr_len = int ($sum_len / $thr_num + 0.5);
            if (($thr_len >= 1000) and ($thr_num > 1) and ($thr_len > $sec_len)){
                $sec_clust = $clust;
                $sec_num = $thr_num;
                $sec_len = $thr_len;
                last;
            }
        }
    }
    if ($sec_num < $top_num * 0.5){
        my $top_hnum = int ($top_num * 0.5);
        @{$clust{$top_clust}} = sort {$a <=> $b} @{$clust{$top_clust}};
        my $top_med = ${$clust{$top_clust}}[$top_hnum];
        if ($sec_clust > 0){
            my $sec_hnum = int ($sec_num * 0.5);
            @{$clust{$sec_clust}} = sort {$a <=> $b} @{$clust{$sec_clust}};
            my $sec_med = ${$clust{$sec_clust}}[$sec_hnum];
            if (($top_med / $sec_med <= $str_max_len_rate) and ($top_med / $sec_med >= $str_min_len_rate)){
                push @{$clust{$top_clust}}, @{$clust{$sec_clust}};
                $top_num += $sec_num;
                delete $clust{$sec_clust};
                $sec_clust = 0;
                $sec_num = 0;
            }
        }
    }
    if (($clust_num > 2) and ($sec_clust > 0)){
        my $pre_clust = 0;
        my $pre_num = 0;
        my $pre_med = 0;
        foreach my $clust (sort {$a <=> $b} keys %clust){
            my $num = scalar @{$clust{$clust}};
            my $hnum = int ($num * 0.5);
            @{$clust{$clust}} = sort {$a <=> $b} @{$clust{$clust}};
            my $med = ${$clust{$clust}}[$hnum];
            my $next_clust = $clust + 1;
            while ($next_clust <= $clust_num){
                last if (exists $clust{$next_clust});
                $next_clust ++;
            }
            my $next_num = 0;
            my $next_hnum = 0;
            my $next_med = 0;
            if (exists $clust{$next_clust}){
                $next_num = scalar @{$clust{$next_clust}};
                $next_hnum = int ($next_num * 0.5);
                @{$clust{$next_clust}} = sort {$a <=> $b} @{$clust{$next_clust}};
                $next_med = ${$clust{$next_clust}}[$next_hnum];
            }
            if ($pre_clust > 0){
                if ((($pre_clust eq $top_clust) and ($next_clust eq $sec_clust) and (exists $clust{$next_clust})) or (($pre_clust eq $sec_clust) and ($next_clust eq $top_clust) and (exists $clust{$next_clust}))){
                    if (($med / $pre_med <= $str_max_len_rate) and ($med / $pre_med >= $str_min_len_rate) and ($next_med / $med <= $str_max_len_rate) and ($next_med / $med >= $str_min_len_rate)){
                        if ($pre_num > $next_num){
                            push @{$clust{$next_clust}}, @{$clust{$clust}};
                            delete $clust{$clust};
                            $pre_clust = 0;
                            next;
                        }
                        else{
                            push @{$clust{$pre_clust}}, @{$clust{$clust}};
                            delete $clust{$clust};
                            $pre_clust = 0;
                            next;
                        }
                    }
                    elsif (($med / $pre_med <= $str_max_len_rate) and ($med / $pre_med > $str_min_len_rate)){
                        push @{$clust{$pre_clust}}, @{$clust{$clust}};
                        delete $clust{$clust};
                        $pre_clust = 0;
                        next;
                    }
                    elsif (($next_med / $med <= $str_max_len_rate) and ($next_med / $med > $str_min_len_rate)){
                        push @{$clust{$next_clust}}, @{$clust{$clust}};
                        delete $clust{$clust};
                        $pre_clust = 0;
                        next;
                    }
                }
                elsif ((($pre_clust eq $top_clust) and ($clust ne $sec_clust)) or (($pre_clust eq $sec_clust) and ($clust ne $top_clust))){
                    if (($med / $pre_med <= $str_max_len_rate) and ($med / $pre_med > $str_min_len_rate)){
                        push @{$clust{$pre_clust}}, @{$clust{$clust}};
                        delete $clust{$clust};
                        $pre_clust = 0;
                        next;
                    }
                }
            }
            $pre_clust = $clust;
            $pre_num = $num;
            $pre_med = $med;
        }
        $top_num = scalar @{$clust{$top_clust}};
        $sec_num = scalar @{$clust{$sec_clust}};
    }
    my $top_med = 0;
    my $sec_med = 0;
    my $top_hnum = int ($top_num * 0.5);
    my $sec_hnum = int ($sec_num * 0.5);
    @{$clust{$top_clust}} = sort {$a <=> $b} @{$clust{$top_clust}};
    $top_med = ${$clust{$top_clust}}[$top_hnum];
    @{$clust{$sec_clust}} = sort {$a <=> $b} @{$clust{$sec_clust}} if (exists $clust{$sec_clust});
    $sec_med = ${$clust{$sec_clust}}[$sec_hnum] if (exists $clust{$sec_clust});
    if ($sec_num >= $min_str_reads){
        return ($top_med, $sec_med, $top_num, $sec_num);
    }
    else{
        return ($top_med, 0, $top_num, 0);
    }
}

