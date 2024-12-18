#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use FindBin qw($Bin);

my $data_dir = "$Bin/../Data";

my $ref_file = '';

my $simple_repeat = '';

my $TE_fasta = '';

my $gap_bed = '';

my $bam_file = '';

my $samtool_path = '';
my $trf_path = '';
my $yass_path = '';
my $multalin_path = '';

my $cores = 1;

my $out_prefix = '';

my $min_mapQ = 1;
my $min_mapQ_idup = 20;

my $min_ins_reads = 2;
my $min_del_reads = 2;
my $min_str_reads = 3;

my $min_VRR = 0.05;
my $min_str_vrr = 0.15;

my $min_maplen = 500;

my $chromium = 0;

my $include_secalign = 0;

my $max_mismatch = 15;

my $min_indel_size = 50;
my $min_str_indel_size = 50;
my $min_ins_str_mei = 200;
my $min_str_len_rate = 0.5;
my $min_str_cn = 0.02;
my $str_max_len_rate = 1.1;

my $indel_rate = 10;

my $min_coverage = 80;
my $min_me_coverage = 60;

my $min_str_identity = 70;

my $dup_find = 1;

my $intersperse_dup_find = 1;

my $mei_find = 1;

my $bp_diff = 100;
my $bp_diff2 = 200;
my $bp_diff3 = 150;
my $ins_bp_diff = 50;
my $max_bp_diff = 300;
my $min_overlap_rate = 0.7;
my $min_overlap_rate_eval = 0.5;
my $max_dist = 10;
my $max_dist2 = 15;
my $max_match_size = 200;
my $min_clip_len = 50;
my $max_insbp_diff = 200;
my $min_inv_size = 100;
my $max_del_size = 5000000;
my $max_inv_size = 500000;
my $max_dup_len = 10000000;
my $min_dup_len = 50;
my $max_bp_sd_small_ins = 5;

my $Mbin_size = 1000000;
my $Mbin_size2 = 100000;

my $config = shift @ARGV;

my $target_chr = shift @ARGV;

my $step1_out = shift @ARGV;

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
    elsif ($arg eq 'min_tr_len'){
        $min_str_indel_size = $value;
    }
    elsif ($arg eq 'min_tr_lrate'){
        $min_str_len_rate = $value;
    }
    elsif ($arg eq 'min_tr_cn'){
        $min_str_cn = $value;
    }
    elsif ($arg eq 'min_tr_ident'){
        $min_str_identity = $value;
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
    elsif ($arg eq 'max_tr_rate'){
        $str_max_len_rate = $value;
    }
    elsif ($arg eq 'samtool_path'){
        $samtool_path = $value;
    }
    elsif ($arg eq 'trf_path'){
        $trf_path = $value;
    }
    elsif ($arg eq 'yass_path'){
        $yass_path = $value;
    }
    elsif ($arg eq 'multalin_path'){
        $multalin_path = $value;
    }
}
close (FILE);

my $max_mismatch2 = $max_mismatch - 5;
my $max_mismatch3 = $max_mismatch - 10;

if ($samtool_path ne ''){
    $ENV{PATH} = "$samtool_path:" . $ENV{PATH};
}
if ($trf_path ne ''){
    $ENV{PATH} = "$trf_path:" . $ENV{PATH};
}
if ($yass_path ne ''){
    $ENV{PATH} = "$yass_path:" . $ENV{PATH};
}
if ($multalin_path ne ''){
    $ENV{PATH} = "$multalin_path:" . $ENV{PATH};
}

my $str_min_len_rate = int (1 / $str_max_len_rate * 100 + 0.5) / 100;

my $temp_dir = "$out_prefix.temp";
system ("mkdir $temp_dir") if (!-d $temp_dir);

my %chr_len;
my %STR;
my %STR2;
my %STR_ME;
my %STR_motif;
my %STR_pos;
my $hg19_flag = 0;
my $total_STR_len = 0;

my $ref_file2 = $ref_file;
$ref_file2 = $1 if ($ref_file2 =~ /(.+)\.gz$/);
my $ref_index = "$ref_file2.fai";
if (!-f $ref_index){
    my $ref_base = $ref_file2;
    $ref_base = $1 if ($ref_file2 =~ /\/(.+?)$/);
    if ($ref_file =~ /\.gz$/){
        system ("gzip -dc $ref_file > $ref_base");
    }
    else{
        system ("ln -s $ref_file");
    }
    system ("samtools faidx $ref_base");
    $ref_index = "$ref_base.fai";
}

open (FILE, $ref_index) or die "$ref_index is not found: $!\n";while (my $line = <FILE>){
    chomp $line;
    my ($chr, $len) = split (/\t/, $line);
    $chr_len{$chr} = $len;
    if ($chr =~ /^chr/){
        $hg19_flag = 1;
    }
}
close (FILE);

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
    my $id = $line[3];
    my $me = $line[6];
    my $motif = $line[7];
    my $motif_size = $line[4];
    $STR{$pos} = $id;
    $STR2{$id} = "$pos=$end=$motif_size=$motif";
    $STR_ME{$id} = $me if ($me ne '-');
    $STR_motif{$id} = $motif;
    my $Mbin1 = int (($pos - 50) / $Mbin_size);
    my $Mbin2 = int (($end + 50) / $Mbin_size);
    ${$STR_pos{$Mbin1}}{$pos} = $end;
    if ($Mbin2 > $Mbin1){
        ${$STR_pos{$Mbin2}}{$pos} = $end;
    }
}
close (FILE);

my %read;
my %TE;
my %TE_len;
my @TE;

if ($TE_fasta ne ''){
    my $header = '';
    my $seq = '';
    open (FILE, $TE_fasta) or die "$TE_fasta is not found: $!\n" if ($TE_fasta !~ /\.gz$/);
    open (FILE, "gzip -dc $TE_fasta |") or die "$TE_fasta is not found: $!\n" if ($TE_fasta =~ /\.gz$/);
    while (my $line = <FILE>){
        chomp $line;
        if ($line =~ /^>(\S+)/){
            if ($seq ne ''){
                $TE{$header} = $seq;
                $TE_len{$header} = length $seq;
                push @TE, $header;
            }
            $seq = '';
            $header = $1;
        }
        else{
            $seq .= uc $line;
        }
    }
    if ($seq ne ''){
        $TE{$header} = $seq;
        $TE_len{$header} = length $seq;
        push @TE, $header;
        $seq = '';
    }
    close (FILE);
}

my %seq;
my $header = '';
my $seq = '';
my $Mb_count = 0;

open (FILE, $ref_file) or die "$ref_file is not found: $!\n" if ($ref_file !~ /\.gz$/);
open (FILE, "gzip -dc $ref_file |") or die "$ref_file is not found: $!\n" if ($ref_file =~ /\.gz$/);
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^>(\S+)/){
        if ($seq ne ''){
            while (length $seq > 0){
                if (length $seq >= $Mbin_size){
                    my $seq_Mb = substr ($seq, 0, $Mbin_size, '');
                    $seq{$Mb_count} = $seq_Mb if ($header eq $target_chr);
                }
                else{
                    $seq{$Mb_count} = $seq if ($header eq $target_chr);
                    $seq = '';
                }
                $Mb_count ++;
            }
        }
        $seq = '';
        $Mb_count = 0;
        $header = $1;
        $header = $1 if ($header =~ /(.+)-[pm]aternal$/);
        last if ($header eq $target_chr) and (scalar keys %seq > 0);
    }
    else{
        $seq .= uc $line;
    }
}
if ($seq ne ''){
    while (length $seq > 0){
        if (length $seq >= $Mbin_size){
            my $seq_Mb = substr ($seq, 0, $Mbin_size, '');
            $seq{$Mb_count} = $seq_Mb if ($header eq $target_chr);
        }
        else{
            $seq{$Mb_count} = $seq if ($header eq $target_chr);
            $seq = '';
        }
        $Mb_count ++;
    }
}
close (FILE);

my %ins_pos_id;
my %ins_bp_id;
my %ins_str_id;
my %ins_strID;
my %ins_strID_add;
my %hit_strID;
my %inv_id;
my %inv_seq;
my %ins_pos_seq;
my %ins_bp_seq;
my %ins_str_seq;
my %read_id;

open (FILE, $step1_out) or die "$step1_out is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($chr, $pos, $type, $len, $info, $info2) = split (/\t/, $line);
    if ($type =~ /INS/){
        my %idist;
        my %id;
        my @info = split (/\|/, $info);
        if ($type ne 'INS-BP2'){
            foreach my $info3 (@info){
                my ($id, $strand, $ilen, $rpos, $ipos) = split (/=/, $info3);
                next if ($strand eq 'BP');
                $id = $1 if ($id =~ /(.+)~[35]$/);
                $id{$id} = "$ilen=$rpos=$strand";
                my $dist = abs ($pos - $ipos);
                $ilen = int ($ilen / 2) * 2;
                push @{${$idist{$dist}}{$ilen}}, $id;
            }
            my $select_id1 = '';
            my $select_id2 = '';
            foreach my $dist (sort {$a <=> $b} keys %idist){
                foreach my $len (sort {scalar @{${$idist{$dist}}{$b}} <=> scalar @{${$idist{$dist}}{$a}}} keys %{$idist{$dist}}){
                    foreach my $id (@{${$idist{$dist}}{$len}}){
                        if ($select_id1 eq ''){
                            $select_id1 = $id;
                            last if ($len >= 100);
                        }
                        elsif ($select_id2 eq ''){
                            $select_id2 = $id;
                            last;
                        }
                    }
                }
            }
            ${$ins_pos_id{$select_id1}}{$pos} = $id{$select_id1};
            $read_id{$select_id1} = 1;
            if ($select_id2 ne ''){
                ${$ins_pos_id{$select_id1}}{$pos + 1} = $id{$select_id2};
                $read_id{$select_id2} = 1;
            }
        }
        else{
            my $clip5_id = '';
            my $clip3_id = '';
            my $clip5_len = 0;
            my $clip3_len = 0;
            my $strand5 = '';
            my $strand3 = '';
            foreach my $info3 (@info){
                my ($id, $strand, $cliplen, $cliptag, $tag) = split (/=/, $info3);
                if ($tag eq 'S'){
                    if ($cliptag eq '5T'){
                        if ($cliplen > $clip5_len){
                            $clip5_id = $id;
                            $clip5_len = $cliplen;
                            $strand5 = $strand;
                        }
                    }
                    else{
                        if ($cliplen > $clip3_len){
                            $clip3_id = $id;
                            $clip3_len = $cliplen;
                            $strand3 = $strand;
                        }
                    }
                }
            }
            ${${$ins_bp_id{$clip5_id}}{$pos}}{'5T'} = "$clip5_len=$strand5";
            ${${$ins_bp_id{$clip3_id}}{$pos}}{'3T'} = "$clip3_len=$strand3";
            $read_id{$clip5_id} = 1;
            $read_id{$clip3_id} = 1;
        }
    }
    elsif ($type =~ /TR-ins/){
        my $strid = $STR{$pos};
        $ins_strID{$strid} = 1;
        my @info2 = split (/=/, $info2);
        if (@info2 >= 4){
            my $id = $info2[2];
            next if ($id eq 'NA');
            my $inf2 = $info2[3];
            if (($type eq 'TR-ins') or ($type eq 'TR-ins1')){
                ${${$ins_str_id{$id}}{$pos}}{1} = $inf2;
            }
            elsif ($type eq 'TR-ins2'){
                ${${$ins_str_id{$id}}{$pos}}{2} = $inf2;
            }
            $read_id{$id} = 1 if ($id ne '');
        }
    }
    elsif ($type eq 'INV'){
        my @info = split (/\|/, $info);
        foreach (@info){
            my ($tag, $id, $cliplen1, $strand1, $cliplen2, $strand2) = split (/=/, $_);
            if ($tag !~ /bp[12]/){
                ${$inv_id{$id}}{$pos} = "$tag=$cliplen1=$strand1=$cliplen2=$strand2";
                $read_id{$id} = 1;
            }
        }
    }
    elsif ($type eq 'REP-BP'){
        next if ($info eq '1,1');
        my @info = split (/\|/, $info2);
        foreach (@info){
            my ($id) = split (/~/, $_);
            $read_id{$id} = 1;
        }
    }
}
close (FILE);


my $chr = $target_chr;

open (FILE, "samtools view $bam_file $chr |") or die "$bam_file is not found:$!\n" if ($bam_file =~ /\.bam$/);
open (FILE, "samtools view -S $bam_file $chr |") or die "$bam_file is not found:$!\n" if ($bam_file =~ /\.sam$/);
open (FILE, "samtools view --reference $ref_file $bam_file $chr |") or die "$bam_file is not found:$!\n" if ($bam_file =~ /\.cram$/);
while (my $line = <FILE>){
    chomp $line;
    my @line = split (/\t/, $line);
    my $read_id = $line[0];
    next if (!exists $read_id{$read_id});
    next if ($line[5] =~ /H/);
    my $tag = $line[1];
    my $sec_align = 0;
    my $strand = 'F';
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
    next if ($sec_align == 1);
    my $seq = $line[9];
    if (exists $ins_pos_id{$read_id}){
        foreach my $ipos (keys %{$ins_pos_id{$read_id}}){
            my ($ilen, $rpos, $istrand) = split (/=/, ${$ins_pos_id{$read_id}}{$ipos});
            next if ($strand ne $istrand);
            next if ($rpos + $ilen > length $seq);
            my $ins_seq = substr ($seq, $rpos, $ilen);
            $ins_pos_seq{$ipos} = $ins_seq;
        }
    }
    if (exists $ins_str_id{$read_id}){
        foreach my $ipos (keys %{$ins_str_id{$read_id}}){
            foreach my $num (keys %{${$ins_str_id{$read_id}}{$ipos}}){
                my @info = split (/,/, ${${$ins_str_id{$read_id}}{$ipos}}{$num});
                my $ins_seq = '';
                foreach my $inf (@info){
                    my ($rpos, $ilen, $istrand) = split (/-/, $inf);
                    next if ($strand ne $istrand);
                    next if ($rpos + $ilen > length $seq);
                    my $insseq = substr ($seq, $rpos, $ilen);
                    $ins_seq .= "$insseq-";
                }
                $ins_seq =~ s/-$// if ($ins_seq =~ /-$/);
                ${$ins_str_seq{$ipos}}{$num} = $ins_seq;
            }
        }
    }
    if (exists $ins_bp_id{$read_id}){
        my $cigar = $line[5];
        my $clip5len = 0;
        my $clip3len = 0;
        if ($cigar =~ /^(\d+)S/){
            $clip5len = $1;
        }
        if ($cigar =~ /(\d+)S$/){
            $clip3len = $1;
        }
        foreach my $ipos (keys %{$ins_bp_id{$read_id}}){
            foreach my $tag (keys %{${$ins_bp_id{$read_id}}{$ipos}}){
                my ($cliplen, $strand) = split (/=/, ${${$ins_bp_id{$read_id}}{$ipos}}{$tag});
                my $clipseq = '';
                if ($tag eq '5T'){
                    next if ($cliplen != $clip5len);
                    $clipseq = substr ($seq, 0, $cliplen);
                }
                elsif ($tag eq '3T'){
                    next if ($cliplen != $clip3len);
                    $clipseq = substr ($seq, -$cliplen, $cliplen);
                }
                ${$ins_bp_seq{$ipos}}{$tag} = $clipseq;
            }
        }
    }
    if (exists $inv_id{$read_id}){
        my $cigar = $line[5];
        my $clip5len = 0;
        my $clip3len = 0;
        if ($cigar =~ /^(\d+)S/){
            $clip5len = $1;
        }
        if ($cigar =~ /(\d+)S$/){
            $clip3len = $1;
        }
        foreach my $ipos (keys %{$inv_id{$read_id}}){
            my ($tag, $cliplen1, $strand1, $cliplen2, $strand2) = split (/=/, ${$inv_id{$read_id}}{$ipos});
            my $clipseq1 = '';
            my $clipseq2 = '';
            if ($tag eq 'bp1'){
                if (($strand eq $strand1) and ($cliplen1 == $clip3len)){
                     $clipseq1 = substr ($seq, -$cliplen1, $cliplen1);
                }
                if (($strand eq $strand2) and ($cliplen2 == $clip3len)){
                    $clipseq2 = substr ($seq, -$cliplen2, $cliplen2);
                }           
            }
            elsif ($tag eq 'bp2'){
                if (($strand eq $strand1) and ($cliplen1 == $clip5len)){
                    $clipseq1 = substr ($seq, 0, $cliplen1);
                }
                if (($strand eq $strand2) and ($cliplen2 == $clip5len)){
                    $clipseq2 = substr ($seq, 0, $cliplen2);
                }
            }
            if ($clipseq1 ne ''){
                ${${$inv_seq{$ipos}}{$read_id}}{1} = "$tag=$clipseq1";
            }
            if ($clipseq1 ne ''){
                ${${$inv_seq{$ipos}}{$read_id}}{2} = "$tag=$clipseq2";
            }
        }
    }
}
close (FILE);

undef %read_id;
undef %ins_pos_id;
undef %ins_bp_id;
undef %ins_str_id;
undef %inv_id;

my $out_file = $step1_out;
$out_file =~ s/\.txt$/.out/;
my $ins_out = $step1_out;
$ins_out =~ s/\.txt$//;
$ins_out .= '.INS.fa';
my $ins_undef_fasta = "$temp_dir/Z.$chr.INS-undefined.fa";
my %INS_result;
my %INS_seq;
my %INS_undef;
my $STRins_multi = 0;
my $STRins_multi_match = 0;

open (OUTD, "> $out_file");
open (FILE, $step1_out) or die "$step1_out is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    my $gt = 'NA';
    my ($chr, $pos, $type, $len, $info, $dp_str, $dp2) = split (/\t/, $line);
    if ($type =~ /INS/){
        my ($read_num, $dp, $ins_rate, $bp_read, $bplen, $dprate) = split (/,/, $dp_str);
        $bplen = 0 if (!defined $bplen);
        if ($dp > 0){
            if ($ins_rate >= 0.7){
                $gt = 'HM';
            }
            else{
                $gt = 'HT';
            }
        }
        my %ins_annot;
        my $str_flag = 0;
        my $tdup_flag = 0;
        my $ins_seq = '';
        my $ins_seq2 = '';
        my $insbp_len = '';
        my $insbp_overlap = 0;
        
        if ($type eq 'INS-BP2'){
#                next if ($ins_rate < $min_VRR);
            $dprate = $dp;
            my $clipseq_3T = '';
            my $clipseq_5T = '';
            $clipseq_3T = ${$ins_bp_seq{$pos}}{'3T'} if (exists ${$ins_bp_seq{$pos}}{'3T'});
            $clipseq_5T = ${$ins_bp_seq{$pos}}{'5T'} if (exists ${$ins_bp_seq{$pos}}{'5T'});
            if (($clipseq_3T eq '') or ($clipseq_5T eq '')){
                print STDERR "INS-BP2: $chr:$pos no-clipped sequences available\n";
                next;
            }
            my ($match, $insseq) = &find_overlap ($clipseq_3T, $clipseq_5T, $chr, $pos); # check end-to-end alignment between 3'-clipped sequence and 5'-clipped sequence
            if ($match == 1){
                $ins_seq = $insseq;
                $insbp_overlap = 1;
            }
            else{
                $ins_seq = $clipseq_3T;
                $ins_seq = $clipseq_5T if (length $clipseq_3T <= length $clipseq_5T);
            }
            $len = length $ins_seq if ($ins_seq ne '');
            $len = 0 if ($match == 0);
            if (exists $INS_result{$pos}){
                my $pre_bp = 0;
                my $pre_result = $INS_result{$pos};
                $pre_bp = $1 if ($pre_result =~ /BP-(\d+)/);
                if ($pre_bp < $read_num){
                    $pre_result =~ s/BP-\d+/BP-$read_num/;
                    my $add_read = $read_num - $pre_bp;
                    my $pre_read = $1 if ($pre_result =~ /RN-(\d+)/);
                    $pre_read += $add_read;
                    $pre_result =~ s/RN-\d+/RN-$pre_read/;
                    $INS_result{$pos} = $pre_result;
                }
                next;
            }
            if ($dp2 =~ /DUP-(\d+)/){
                $ins_annot{'BP'} = "0=$1=DUP";
            }
            elsif ($dp2 =~ /DEL-(\d+)/){
                $ins_annot{'BP'} = "0=$1=DEL";
            }
        }
        else{
            $ins_seq = $ins_pos_seq{$pos} if (exists $ins_pos_seq{$pos});
        }
        my $BP = 0;
        if ($type =~ /BP/){
            $BP = $read_num;
        }
        else{
            $BP = $bp_read if (defined $bp_read);
        }
        if ($ins_seq eq ''){
            $INS_result{$pos} = "$chr\t$pos\t$type\t$len\tRN-$read_num;VRR-$ins_rate;GT-$gt;BP-$BP";
            next;
        }
        if (($dup_find == 1) and ($len >= 20)){
            my ($dup_pos1, $dup_pos2, $direction) = &find_dup_1 ($ins_seq, 'tdup', $chr, $pos);
            if (($len < 100) and ($dup_pos1 == 0) and (exists ${$ins_pos_seq{$chr}}{$pos + 1})){
                my $ins_seq2 = ${$ins_pos_seq{$chr}}{$pos + 1};
                ($dup_pos1, $dup_pos2, $direction) = &find_dup_1 ($ins_seq2, 'tdup', $chr, $pos);
            }
            if ($dup_pos1 > 0){
                my $duplen = $dup_pos2 - $dup_pos1 + 1;
                $ins_annot{'DUP'} = "$dup_pos1=$duplen=$direction";
                $tdup_flag = 1;
            }
            
        }
        if (($mei_find == 1) and (@TE > 0) and ($len >= 90)){
            my ($me_type, $overlaplen, $cn, $sec_me) = &find_me ($ins_seq, $chr, $pos, 'NA');
            if ($overlaplen > 0){
                $ins_annot{'ME'} = "$me_type=$overlaplen=$cn=$sec_me";
            }
        }
        my $annot = '';
        my $annot2 = '';
        if (exists $ins_annot{'BP'}){
            my ($dup_pos1, $dlen, $bptype) = split (/=/, $ins_annot{'BP'});
            if ($bptype eq 'DUP'){
                $annot = "DUPBPLEN-$dlen;";
            }
            elsif ($bptype eq 'DEL'){
                $annot = "DELBPLEN-$dlen;";
            }
        }
        if (exists $ins_annot{'DUP'}){
            my ($dup_pos1, $dlen, $dir) = split (/=/, $ins_annot{'DUP'});
            if ($dlen >= $min_dup_len){
                $annot .= "DUPPOS-$dup_pos1,DUPLEN-$dlen;";
                $annot2 = 'DUP,';
                $type = 'INS(DUP)';
                $type = 'INS(DUP:R)' if ($dir eq 'R');
            }
        }
        if (exists $ins_annot{'ME'}){
            my ($me_type, $overlaplen, $me_cn, $sec_me) = split (/=/, $ins_annot{'ME'});
            $annot .= "MEI-$me_type,MEILEN-$overlaplen,MEICN-$me_cn;" if ($sec_me eq '');
            $annot .= "MEI-$me_type,MEILEN-$overlaplen,MEICN-$me_cn,MEI2-$sec_me;" if ($sec_me ne '');
            $annot2 .= "MEI($me_type),";
        }
        else{
            my ($match_cn, $match_motif) = &TRF_test1 ($ins_seq, $chr);
            if (($match_cn > 0) and ($len > 0)){
                my $ibin = int ($pos / $Mbin_size);
                my %hit_str;
                my $strid = '';
                my $hit_strid = '';
                my $ins_cn = 0;
                my $hit_str_pos = 0;
                my $str_flag = 0;
                foreach my $spos (sort {$a <=> $b} keys %{$STR_pos{$ibin}}){
                    last if ($spos > $pos + 50);
                    my $send = ${$STR_pos{$ibin}}{$spos};
                    next if ($send < $pos - 50);
                    my $strid2 = $STR{$spos};
                    my $diff = 0;
                    if (($pos >= $spos) and ($pos <= $send)){
                        $diff = $pos - $spos;
                        $diff = $send - $pos if ($diff < $send - $pos);
                    }
                    elsif ($pos > $send){
                        $diff = $send - $pos;
                    }
                    elsif ($pos < $spos){
                        $diff = $pos - $spos;
                    }
                    $hit_str{$strid2} = $diff;
                }
                if (scalar keys %hit_str > 0){
                    foreach my $sid (sort {$hit_str{$b} <=> $hit_str{$a}} keys %hit_str){
                        $strid = $sid;
                        last;
                    }
                }
                if ($strid ne ''){
                    my $match_len = 0;
                    my ($spos, $send, $motif_size, $motif) = split (/=/, $STR2{$strid});
                    if (($len >= $motif_size) and ($motif_size <= 4)){
                        ($match_len) = &motif_test ($ins_seq, $motif);
                    }
                    elsif (($len >= $motif_size * 2) and ($motif_size > 4)){
                        ($match_len, my $tmotif) = &TRF_test2 ($ins_seq, $motif, $chr);
                        my $match_cov1 = int ($match_len / $len * 10 + 0.5) / 10;
                        if ($match_cov1 >= 0.5){
                            my $motif_match_flag = 0;
                            my $motifx2 = $motif . $motif;
                            my $tmotifx2 = $tmotif . $tmotif;
                            if ($motifx2 =~ /$tmotif/){
                                $motif_match_flag = 1;
                            }
                            elsif ($tmotifx2 =~ /$motif/){
                                $motif_match_flag = 1;
                            }
                            my $tmotif_size = length $tmotif;
                            if ($motif_match_flag == 0){
                                my ($ident, $cov, $match) = &multalin_ins1 ($ins_seq, $tmotif, $chr) if ($tmotif_size >= $len);
                                ($ident, $cov, $match) = &multalin_ins1 ($tmotif, $ins_seq, $chr) if ($tmotif_size < $len);
                                if (($ident < $min_str_identity) or ($cov < 80)){
                                    $match_len = 0;
                                    $match_cov1 = 0;
                                }
                            }
                        }
                        if ($match_cov1 < 0.5){
                            my ($ident, $cov, $match) = &multalin_ins2 ($ins_seq, $motif, $chr);
                            if (($ident >= $min_str_identity) and ($match >= $len * $min_str_len_rate) and ($match > $match_len)){
                                $match_len = $match;
                            }
                        }
                    }
                    else{
                        if ($motif_size >= $len){
                            my ($ident, $cov, $match) = &multalin_ins1 ($ins_seq, $motif, $chr);
                            if (($ident >= $min_str_identity) and ($match >= $len * $min_str_len_rate)){
                                $match_len = $match;
                            }
                        }
                        else{
                            my ($ident, $cov, $match) = &multalin_ins2 ($ins_seq, $motif, $chr);
                            if (($ident >= $min_str_identity) and ($match >= $len * $min_str_len_rate)){
                                $match_len = $match;
                            }
                        }
                    }
                    my $match_cov = $match_len / $len;
                    if ($match_cov >= $min_str_len_rate){
                        if (($motif_size == 1) and ($match_cov < 0.8)){
                            $str_flag = 0;
                        }
                        else{
                            $str_flag = 1;
                            $ins_cn = int ($match_len / $motif_size * 10 + 0.5) / 10;
                            $ins_cn = int ($match_len / $motif_size * 100 + 0.5) / 100 if ($ins_cn == 0);
                            if ($ins_cn < $min_str_cn){
                                $str_flag = 0;
                            }
                            else{
                                $hit_strid = $strid;
                                $hit_str_pos = $spos;
                            }
                        }
                    }
                }
                if ($str_flag == 1){
                    my $str_line2 = "$chr\t$hit_str_pos\tTR-ins\t$len\tRN-$read_num;RD-1;GT-$gt;TRCN-$ins_cn;TR-$hit_strid;BP-$BP";
                    $annot2 = 'NA' if ($annot2 eq '');
                    my $ins_seq_info = "$chr:$hit_str_pos-$len $annot2 $ins_seq";
                    push @{$ins_strID_add{$hit_strid}}, "$str_line2==$hit_str_pos==$ins_seq_info";
                    next;
                }
                else{
                    $annot .= "MOTIF-$match_motif:$match_cn;" if ($match_cn >= 1.9);
                }
            }
        }
        if ($bplen > 0){
            $annot .= "BPLEN2-$bplen;";
        }
        if ($dprate > 0){
            $annot .= "DPR-$dprate;";
        }
        $annot =~ s/;$// if ($annot =~ /;$/);
        $annot2 =~ s/,$// if ($annot2 =~ /,$/);
        if ($annot eq ''){
            $INS_result{$pos} = "$chr\t$pos\t$type\t$len\tRN-$read_num;VRR-$ins_rate;GT-$gt;BP-$BP";
        }
        else{
            $INS_result{$pos} = "$chr\t$pos\t$type\t$len\t$annot;RN-$read_num;VRR-$ins_rate;GT-$gt;BP-$BP";
        }

        if ($ins_seq ne ''){
            $annot2 = 'unknown' if ($annot2 eq '');
            $annot2 = 'unknown,INS-BP2(partial)' if ($annot2 eq 'unknown') and ($type eq 'INS-BP2') and ($insbp_overlap == 0);
            push @{$INS_seq{$pos}}, "$chr:$pos-$len $annot2 $ins_seq";
            if ($annot2 =~ /unknown/){
                $INS_undef{$pos} = "$chr:$pos $annot2 $ins_seq";
            }
        }
    }
    
    if ($type eq 'INV'){
        my @info = split (/\|/, $info);
        my %inv_pair;
        my %inv_rev_align;
        my $inv_pair = 0;
        my $inv_rev_align = 0;
        my @info2;
        foreach my $info2 (@info){
            if ($info2 =~ /^bp[12]+=/){
                my ($tag, $id, $if1, $if2, $if3, $if4, $align_tag) = split (/=/, $info2);
                $inv_pair{$id} = 1;
                $inv_rev_align{$id} = 1 if ($tag eq 'bp12');
            }
            else{
                push @info2, $info2;
            }
        }
        $inv_pair = scalar keys %inv_pair;
        $inv_rev_align = scalar keys %inv_rev_align;
        my ($dp, $inv_rate) = split (/,/, $dp_str);
        if ($dp > 0){
            if ($inv_pair / $dp >= 0.7){
                $gt = 'HM';
            }
            else{
                $gt = 'HT';
            }
        }
        if ($info =~ /bp1=|bp2=/){
            print OUTD "$chr\t$pos\t$type\t$len\tRN-$inv_pair;VRR-$inv_rate;GT-$gt\n" if ($inv_rev_align == 0);
            print OUTD "$chr\t$pos\t$type\t$len\tRN-$inv_pair;ALIGN-$inv_rev_align;VRR-$inv_rate;GT-$gt\n" if ($inv_rev_align > 0);
        }
        else{
            if ($inv_rev_align > 0){
                print OUTD "$chr\t$pos\t$type\t$len\tRN-$inv_pair;ALIGN-$inv_rev_align;VRR-$inv_rate;GT-$gt\n";
            }
            elsif (exists $inv_seq{$pos}){
                my ($match_flag) = &check_inv ($chr, $pos, $len);
                if ($match_flag >= 1){
                    print OUTD "$chr\t$pos\t$type\t$len\tRN-$inv_pair;ALIGN-1;VRR-$inv_rate;GT-$gt\n";
                }
            }
        }
    }
    elsif ($type eq 'DUP'){
        my @info = split (/\|/, $info);
        my $ave_mTD = 0;
        my @mTD;
        foreach (@info){
            my @item = split (/=/, $_);
            if (@item > 5){
                push @mTD, $item[5];
            }
        }
        my ($dprate, $flank_dp, $dup_rate, $read_num, $bp5, $bp3, $ins_read, $insdup_len) = split (/,/, $dp_str);
        $ins_read = 0 if (!defined $ins_read);
        next if ($dprate < 1) and ($len >= 500);
        if (@mTD > 0){
            my $sum = 0;
            map{$sum += $_} @mTD;
            $ave_mTD = int ($sum / @mTD + 0.5);
            $ave_mTD += 3;
        }
        else{
            $ave_mTD = int ($dprate * 2 + 0.5);
        }
        if ($ave_mTD <= 3){
            $gt = 'HT';
            $gt = 'HM' if ($dup_rate >= 0.8);
        }
        else{
            if ($dup_rate < 0.7){
                $gt = 'HT';
            }
            else{
                $gt = 'HM';
            }
        }

        my $dup_info = '';
        my $ref_seq = '';
        if (($mei_find == 1) and (@TE > 0) and ($len < 50000)){
            if ($ref_seq eq ''){
                my $end = $pos + $len - 1;
                my $divnum1 = int ($pos / $Mbin_size);
                my $res1 = $pos % $Mbin_size;
                my $divnum2 = int ($end / $Mbin_size);
                my $res2 = $end % $Mbin_size;
                if ($divnum1 == $divnum2){
                    if ($res1 > 0){
                        $ref_seq = substr ($seq{$divnum1}, $res1 - 1, $len);
                    }
                    else{
                        $ref_seq = substr ($seq{$divnum1}, 0, $len - 1);
                    }
                }
                else{
                    my $subseq1 = '';
                    $subseq1 = substr ($seq{$divnum1}, $res1 - 1) if ($res1 > 0);
                    my $subseq2 = substr ($seq{$divnum2}, 0, $res2);
                    $ref_seq = $subseq1 . $subseq2;
                }
            }
            my ($me_type, $overlaplen, $cn, $sec_me) = &find_me ($ref_seq, $chr, $pos, 'NA');
            if ($overlaplen > 0){
                $dup_info .= "MEI-$me_type,MEILEN-$overlaplen,MEICN-$cn;" if ($sec_me eq '');
                $dup_info .= "MEI-$me_type,MEILEN-$overlaplen,MEICN-$cn,MEI2-$sec_me;" if ($sec_me ne '');
            }
        }
        if ($ins_read > 0){
            $dup_info .= "INSREAD-$ins_read;"
        }
        if ($insdup_len > 0){
            $dup_info .= "INSLEN-$insdup_len;"
        }
        $dup_info =~ s/;$// if ($dup_info =~ /;$/);
        print OUTD "$chr\t$pos\tDUP\t$len\tRN-$read_num;mTD-$ave_mTD;DPR-$dprate;VRR-$dup_rate;GT-$gt\n" if ($dup_info eq '');
        print OUTD "$chr\t$pos\tDUP\t$len\t$dup_info;RN-$read_num;mTD-$ave_mTD;DPR-$dprate;VRR-$dup_rate;GT-$gt\n" if ($dup_info ne '');
    }
    elsif ($type eq 'DEL'){
        my $read_num = $info;
        my ($dp_rate, $flank_dp, $del_rate, $bp_read)  = split (/,/, $dp2);
        if (($flank_dp > 0) and ($len < 500)){
            if ($read_num / $flank_dp >= 0.7){
                $gt = 'HM';
            }
            else{
                $gt = 'HT';
            }
            if (($len >= 100) and ($dp_rate <= 0.1)){
                $gt = 'HM';
            }
        }
        elsif ($flank_dp > 0){
            if ($dp_rate <= 0.2){
                $gt = 'HM';
            }
            elsif ($dp_rate >= 0.4){
                $gt = 'HT';
            }
            else{
                if ($read_num / $flank_dp >= 0.7){
                    $gt = 'HM';
                }
                else{
                    $gt = 'HT';
                }
            }
        }
        else{
            $gt = 'HT';
        }
        my $BP = 0;
        $BP = $bp_read if (defined $read_num);
        print OUTD "$chr\t$pos\tDEL\t$len\tRN-$read_num;DPR-$dp_rate;VRR-$del_rate;GT-$gt;BP-$BP\n" if ($len >= 100);
        print OUTD "$chr\t$pos\tDEL\t$len\tRN-$read_num;VRR-$del_rate;GT-$gt;BP-$BP\n" if ($len < 100);
    }
    elsif ($type eq 'DEL-BP'){
        my ($rnum1, $rnum2) = split (/,/, $info);
        my $dprate_dp  = $dp2;
        my ($dprate, $flank_dp, $del_rate) = split (/,/, $dprate_dp);
        my $read_num = int (($rnum1 + $rnum2) * 0.5 + 0.5);
        next if ($dprate > 0.9) and ($len >= 500);
        next if ($dprate > 0.85) and ($len >= 1000) and ($del_rate < 0.2);
        if ($flank_dp > 0){
            if ($del_rate >= 0.9){
                $gt = 'HM';
            }
            elsif ($del_rate <= 0.6){
                $gt = 'HT';
            }
        }
        if ($dprate <= 0.1){
            $gt = 'HM';
        }
        elsif ($dprate >= 0.3){
            $gt = 'HT';
        }
        print OUTD "$chr\t$pos\tDEL-BP\t$len\tRN1-$rnum1;RN2-$rnum2;INSLEN-$dp_str;DPR-$dprate;VRR-$del_rate;GT-$gt;BP-$read_num\n";
    }
    elsif ($type eq 'REP'){
        my $read_num = $info;
        my ($flank_dp, $rep_rate) = split (/,/, $dp2);
        if ($flank_dp > 0){
            if ($read_num / $flank_dp >= 0.7){
                $gt = 'HM';
            }
            else{
                $gt = 'HT';
            }
        }
        print OUTD "$chr\t$pos\tREP\t$len\tRN-$read_num;VRR-$rep_rate;GT-$gt\n";
    }
    elsif ($type eq 'REP-BP'){
        my ($rnum1, $rnum2) = split (/,/, $info);
        next if ($info eq '1,1');
        my ($flank_dp, $rep_rate) = split (/,/, $dp2);
        my $pos2 = $pos + $len - 1;
        my $read_num = $rnum1 + $rnum2;
        next if ($read_num < $flank_dp * 0.1) and ($flank_dp > 0);
        if ($flank_dp > 0){
            if ($read_num / $flank_dp >= 0.7){
                $gt = 'HM';
            }
            else{
                $gt = 'HT';
            }
        }
        print OUTD "$chr\t$pos\tREP-BP\t$len\tRN1-$rnum1;RN2-$rnum2;VRR-$rep_rate;GT-$gt\n";
    }
    elsif ($type =~ /TR-/){
        my ($read_cov, $unit_size, $readid, $ipos_str, $ipos, $insbp, $dpr, $sec_str) = split (/=/, $dp_str);
        $dpr = 0 if (!defined $dpr);
        my $strid = $STR{$pos};
        my ($spos, $send, $motif_size, $motif) = split (/=/, $STR2{$strid});
        my $slen = $send - $spos + 1;
        my %mei;
        my $ins_mei = '';
        my $mei_len = 0;
        my $mei_cn = 0;
        my $ilen = $len;
        $gt = 'HM';
        if ($info !~ /\//){
            my $rate = int ($info / $read_cov * 100 + 0.5) / 100;
            if ($rate < 0.7){
                $gt = 'HT';
            }
        }
        my $ins_num = 1;
        $ins_num = 2 if ($type eq 'TR-ins2');
        my $str_line = '';
        if ($type =~ /TR-del/){
            my $del_cn = int ($len / $motif_size * 10 + 0.5) / 10;
            $del_cn = int ($len / $motif_size * 100 + 0.5) / 100 if ($del_cn == 0);
            if (($del_cn < $min_str_cn) and ($len < $min_str_indel_size)){
                next;
            }
        }
        if (($type =~ /TR-ins/) and (exists ${$ins_str_seq{$pos}}{$ins_num}) and ($dpr < 1.2)){
            my $match_flag = 0;
            my $hit_strid = '';
            my $ins_seq = '';
            $ins_seq = ${$ins_str_seq{$pos}}{$ins_num};
            my @ins_seq = split (/-/, $ins_seq);
            my $annot = $strid;
            my $ins_cn = 0;
            foreach my $insseq (@ins_seq){
                if (length $insseq >= $min_ins_str_mei){
                    my ($me_type, $overlaplen, $cn) = &find_me ($insseq, $chr, $pos, $strid);
                    if ($overlaplen > 0){
                        my $str_motif = $STR_motif{$strid};
                        my $motif_size = length $str_motif;
                        if (($me_type eq 'ALU') and ($overlaplen >= 180) and ($motif_size <= 150)){
                            $mei{$me_type} += $overlaplen;
                        }
                        elsif (($me_type eq 'SVA') and ($overlaplen >= 1000)){
                            $mei{$me_type} += $overlaplen;
                        }
                        elsif (($me_type eq 'LINE1') and ($overlaplen >= 3000)){
                            $mei{$me_type} += $overlaplen;
                        }
                        elsif (($me_type eq 'HERVK') and ($overlaplen >= 3000)){
                            $mei{$me_type} += $overlaplen;
                        }
                    }
                }
            }
            
            if ((scalar keys %mei == 0) and ($min_str_len_rate > 0)){
                my $alt_str_info = '';
                my $sum_ins_len = 0;
                my $sum_ins_match_len = 0;
                my $sum_ins_match_len2 = 0;
                my $sec_str_match_motif = '';
                my %match_insseq;
                my %match_insseq2;
                foreach my $insseq (@ins_seq){
                    my $ins_len = length $insseq;
                    next if ($ins_len == 0);
                    next if ($motif_size <= 5) and ($ins_len < $motif_size);
                    $sum_ins_len += $ins_len;
                    my $match_len = 0;
                    if (($ins_len >= $motif_size) and ($motif_size <= 4)){
                        ($match_len) = &motif_test ($insseq, $motif);
                    }
                    elsif (($ins_len >= $motif_size * 2) and ($motif_size > 4)){
                        ($match_len, my $tmotif) = &TRF_test2 ($insseq, $motif, $chr);
                        my $match_cov1 = int ($match_len / $ins_len * 10 + 0.5) / 10;
                        if ($match_cov1 >= 0.5){
                            my $motif_match_flag = 0;
                            my $motifx2 = $motif . $motif;
                            my $tmotifx2 = $tmotif . $tmotif;
                            if ($motifx2 =~ /$tmotif/){
                                $motif_match_flag = 1;
                            }
                            elsif ($tmotifx2 =~ /$motif/){
                                $motif_match_flag = 1;
                            }
                            if ($motif_match_flag == 0){
                                my $tmotif_size = length $tmotif;
                                my ($ident, $cov, $match) = &multalin_ins1 ($insseq, $tmotif, $chr) if ($tmotif_size >= $ins_len);
                                ($ident, $cov, $match) = &multalin_ins1 ($tmotif, $insseq, $chr) if ($tmotif_size < $ins_len);
                                if (($ident < $min_str_identity) or ($cov < 80)){
                                    $match_len = 0;
                                    $match_cov1 = 0;
                                }
                            }
                        }
                        if ($match_cov1 < 0.5){
                            my ($ident, $cov, $match) = &multalin_ins2 ($insseq, $motif, $chr);
                            if (($ident >= $min_str_identity) and ($match >= $ins_len * $min_str_len_rate) and ($match > $match_len)){
                                $match_len = $match;
                            }
                        }
                    }
                    else{
                        if ($motif_size >= $ins_len){
                            my ($ident, $cov, $match) = (0, 0, 0);
                            if ($ins_len <= 5){
                                my $motifR2 = $motif . $motif;
                                if ($motifR2 =~ /$insseq/){
                                    $ident = 100;
                                    $cov = 100;
                                    $match = 5;
                                }
                                elsif (($min_str_identity <= 80) and ($ins_len == 5)){
                                    my $mflag = 0;
                                    for (my $i = 0; $i <= 4; $i++){
                                        my $minsseq = $insseq;
                                        substr ($minsseq, $i, 1, '[ACGT]');
                                        if ($motifR2 =~ /$minsseq/){
                                            $mflag = 1;
                                            last;
                                        }
                                    }
                                    if ($mflag == 1){
                                        $ident = 80;
                                        $cov = 80;
                                        $match = 4;
                                    }
                                }
                                elsif (($min_str_identity <= 75) and ($ins_len == 4)){
                                    my $mflag = 0;
                                    for (my $i = 0; $i <= 3; $i++){
                                        my $minsseq = $insseq;
                                        substr ($minsseq, $i, 1, '[ACGT]');
                                        if ($motifR2 =~ /$minsseq/){
                                            $mflag = 1;
                                            last;
                                        }
                                    }
                                    if ($mflag == 1){
                                        $ident = 75;
                                        $cov = 75;
                                        $match = 3;
                                    }
                                }
                            }
                            else{
                                ($ident, $cov, $match) = &multalin_ins1 ($insseq, $motif, $chr);
                            }
                            if (($ident >= $min_str_identity) and ($match >= $ins_len * $min_str_len_rate)){
                                $match_len = $match;
                            }
                        }
                        else{
                            my ($ident, $cov, $match) = &multalin_ins2 ($insseq, $motif, $chr);
                            if (($ident >= $min_str_identity) and ($match >= $ins_len * $min_str_len_rate)){
                                $match_len = $match;
                            }
                        }
                    }
                    if ($match_len >= $ins_len * 0.5){
                        $sum_ins_match_len += $match_len;
                        $match_insseq{$insseq} = 1 if ($match_len > 0);
                    }
                    else{
#                        $seg_align{$ins_count} = 0 if ($ins_len >= 50);
                    }
                }
                my $match_cov1 = 0;
                my $match_cov2 = 0;
                $match_cov1 = int ($sum_ins_match_len / $sum_ins_len * 10 + 0.5) / 10 if ($sum_ins_len > 0);
                if ($match_cov1 < 0.8){
                    if ($sec_str eq 'NA'){
                        my $ibin = int ($ipos / $Mbin_size);
                        my %hit_str;
                        foreach my $spos (sort {$a <=> $b} keys %{$STR_pos{$ibin}}){
                            last if ($spos > $ipos + 20);
                            my $send = ${$STR_pos{$ibin}}{$spos};
                            next if ($send < $ipos - 20);
                            my $strid2 = $STR{$spos};
                            next if ($strid2 eq $strid);
                            my $diff = 0;
                            if (($ipos >= $spos) and ($ipos <= $send)){
                                $diff = $ipos - $spos;
                                $diff = $send - $ipos if ($diff < $send - $ipos);
                            }
                            elsif ($ipos > $send){
                                $diff = $send - $ipos;
                            }
                            elsif ($ipos < $spos){
                                $diff = $ipos - $spos;
                            }
                            $hit_str{$strid2} = $diff;
                        }
                        if (scalar keys %hit_str > 0){
                            foreach my $sid (sort {$hit_str{$b} <=> $hit_str{$a}} keys %hit_str){
                                $sec_str = $sid;
                                last;
                            }
                        }
                    }
                    if ($sec_str ne 'NA'){
                        $sum_ins_len = 0;
                        foreach my $insseq (@ins_seq){
                            my $ins_len = length $insseq;
                            next if ($ins_len == 0);
                            my ($spos2, $send2, $motif_size2, $motif2) = split (/=/, $STR2{$sec_str});
                            next if ($motif_size2 <= 5) and ($ins_len < $motif_size2);
                            $sum_ins_len += $ins_len;
                            my $match_len2 = 0;
                            if (($ins_len >= $motif_size2) and ($motif_size2 <= 4)){
                                ($match_len2) = &motif_test ($insseq, $motif2);
                            }
                            elsif (($ins_len >= $motif_size2 * 2) and ($motif_size2 > 4)){
                                ($match_len2, my $tmotif2) = &TRF_test2 ($insseq, $motif2, $chr);
                                my $match_cov2 = int ($match_len2 / $ins_len * 10 + 0.5) / 10;
                                if ($match_cov2 >= 0.5){
                                    my $motif_match_flag = 0;
                                    my $motifx2 = $motif2 . $motif2;
                                    my $tmotifx2 = $tmotif2 . $tmotif2;
                                    if ($motifx2 =~ /$tmotif2/){
                                        $motif_match_flag = 1;
                                    }
                                    elsif ($tmotifx2 =~ /$motif2/){
                                        $motif_match_flag = 1;
                                    }
                                    if ($motif_match_flag == 0){
                                        my $tmotif_size = length $tmotif2;
                                        my ($ident, $cov, $match) = &multalin_ins1 ($insseq, $tmotif2, $chr) if ($tmotif_size >= $ins_len);
                                        ($ident, $cov, $match) = &multalin_ins1 ($tmotif2, $insseq, $chr) if ($tmotif_size < $ins_len);
                                        if (($ident < $min_str_identity) or ($cov < 80)){
                                            $match_len2 = 0;
                                            $match_cov2 = 0;
                                        }
                                    }
                                }
                                if ($match_cov2 < 0.5){
                                    my ($ident, $cov, $match) = &multalin_ins2 ($insseq, $motif2, $chr);
                                    if (($ident >= $min_str_identity) and ($match >= $ins_len * $min_str_len_rate) and ($match > $match_len2)){
                                        $match_len2 = $match;
                                    }
                                }
                            }
                            else{
                                if ($motif_size2 >= $ins_len){
                                    my ($ident, $cov, $match) = (0, 0, 0);
                                    if ($ins_len <= 5){
                                        my $motifR2 = $motif2 . $motif2;
                                        if ($motifR2 =~ /$insseq/){
                                            $ident = 100;
                                            $cov = 100;
                                            $match = 5;
                                        }
                                        elsif (($min_str_identity <= 80) and ($ins_len == 5)){
                                            my $mflag = 0;
                                            for (my $i = 0; $i <= 4; $i++){
                                                my $minsseq = $insseq;
                                                substr ($minsseq, $i, 1, '[ACGT]');
                                                if ($motifR2 =~ /$minsseq/){
                                                    $mflag = 1;
                                                    last;
                                                }
                                            }
                                            if ($mflag == 1){
                                                $ident = 80;
                                                $cov = 80;
                                                $match = 4;
                                            }
                                        }
                                        elsif (($min_str_identity <= 75) and ($ins_len == 4)){
                                            my $mflag = 0;
                                            for (my $i = 0; $i <= 3; $i++){
                                                my $minsseq = $insseq;
                                                substr ($minsseq, $i, 1, '[ACGT]');
                                                if ($motif2 =~ /$minsseq/){
                                                    $mflag = 1;
                                                    last;
                                                }
                                            }
                                            if ($mflag == 1){
                                                $ident = 75;
                                                $cov = 75;
                                                $match = 3;
                                            }
                                        }
                                    }
                                    else{
                                        ($ident, $cov, $match) = &multalin_ins1 ($insseq, $motif2, $chr);
                                    }
                                    if (($ident >= $min_str_identity) and ($match >= $ins_len * $min_str_len_rate)){
                                        $match_len2 = $match;
                                    }
                                }
                                else{
                                    my ($ident, $cov, $match) = &multalin_ins2 ($insseq, $motif2, $chr);
                                    if (($ident >= $min_str_identity) and ($match >= $ins_len * $min_str_len_rate)){
                                        $match_len2 = $match;
                                    }
                                }
                            }
                            if ($match_len2 >= $ins_len * 0.5){
                                $sum_ins_match_len2 += $match_len2;
                                $match_insseq2{$insseq} = 1;
                                $alt_str_info = "$sec_str:$motif2:$spos2" if ($alt_str_info eq '');
                            }
                            else{
#                                $seg_align2{$ins_count2} = 0 if ($ins_len >= 50);
                            }
                        }
                        $match_cov2 = int ($sum_ins_match_len2 / $sum_ins_len * 10 + 0.5) / 10 if ($sum_ins_len > 0);
                    }
                }
                if ($match_cov1 >= $match_cov2){
                    if ($match_cov1 >= $min_str_len_rate){
                        $match_flag = 1;
                        $sum_ins_match_len = $len if ($sum_ins_match_len > $len);
                        if (($motif_size == 1) and ($match_cov1 < 0.8)){
                            $match_flag = 0;
                        }
                        else{
                            $hit_strid = $strid;
                            $ins_cn = int ($sum_ins_match_len / $motif_size * 10 + 0.5) / 10;
                            $ins_cn = int ($sum_ins_match_len / $motif_size * 100 + 0.5) / 100 if ($ins_cn == 0);
                            if ($ins_cn < $min_str_cn){
                                $match_flag = 0;
                            }
                        }
                    }
                }
                else{
                    if ($match_cov2 >= $min_str_len_rate){
                        $match_flag = 1;
                        $sum_ins_match_len2 = $len if ($sum_ins_match_len2 > $len);
                        my ($strid2, $motif2, $pos2) = split (/:/, $alt_str_info);
                        my $motif_size2 = length $motif2;
                        if (($motif_size2 == 1) and ($match_cov2 < 0.8)){
                            $match_flag = 0;
                        }
                        else{
                            $hit_strid = $strid2;
                            $ins_cn = int ($sum_ins_match_len2 / length ($motif2) * 10 + 0.5) / 10;
                            $ins_cn = int ($sum_ins_match_len2 / length ($motif2) * 100 + 0.5) / 100if ($ins_cn == 0);
                            if ($ins_cn < $min_str_cn){
                                $match_flag = 0;
                            }
                        }
                    }
                }
                if ($match_flag == 0){
                    $ins_seq =~ s/-//g if ($ins_seq =~ /-/);
                }
                else{
                    $ins_seq = '';
                    foreach my $insseq (@ins_seq){
                        if ($match_cov1 >= $match_cov2){
                            if (exists $match_insseq{$insseq}){
                                $ins_seq .= $insseq;
                            }
                        }
                        else{
                            if (exists $match_insseq2{$insseq}){
                                $ins_seq .= $insseq;
                            }
                        }
                    }
                }
            }
            else{
                $ins_seq =~ s/-//g if ($ins_seq =~ /-/);
                $match_flag = 2;
            }
            if (scalar keys %mei > 0){
                $ins_seq =~ s/-//g if ($ins_seq =~ /-/);
                foreach my $mei (sort {$mei{$b} <=> $mei{$a}} keys %mei){
                    $ins_mei = $mei;
                    $mei_len = $mei{$mei};
                    $mei_cn = int ($mei_len / $TE_len{$mei} * 10 + 0.5) / 10;
                    last;
                }
                $annot = "MEI($ins_mei),$strid" if ($ins_mei ne '');
                push @{$INS_seq{$ipos}}, "$chr:$ipos-$ilen $annot $ins_seq";
                my $vrr = int ($info / $read_cov * 100 + 0.5) / 100;
                my ($pos2, $end2) = split (/=/, $STR2{$strid});
                if (($ipos >= $pos2) and ($ipos <= $end2)){
                    $str_line = "$chr\t$ipos\tINS\t$ilen\tMEI-$ins_mei,MEILEN-$mei_len,MEICN-$mei_cn;RN-$info;VRR-$vrr;GT-$gt;TR-$strid";
                }
                else{
                    $str_line = "$chr\t$ipos\tINS\t$ilen\tMEI-$ins_mei,MEILEN-$mei_len,MEICN-$mei_cn;RN-$info;VRR-$vrr;GT-$gt";
                }
                $str_line .= ";BP-$insbp" if (defined $insbp);
            }
            elsif ($match_flag == 0){
                if ($len >= $min_indel_size){
                    my ($pos2, $end2) = split (/=/, $STR2{$strid});
                    my $match_flag2 = 0;
                    $match_flag2 = 1 if ($ipos >= $pos2) and ($ipos <= $end2);
                    my $vrr = int ($info / $read_cov * 100 + 0.5) / 100;
                    my $ins_len = length $ins_seq;
                    my ($match_cn, $match_motif) = &TRF_test1 ($ins_seq, $chr);
                    $match_motif = '' if ($match_cn < 1.9);
                    if ($match_flag2 == 1){
                        $annot = "Weak-TR-homology,$strid";
                        $annot .= ",$match_motif:$match_cn" if ($match_motif ne '');
                        $str_line = "$chr\t$ipos\tINS\t$ilen\tRN-$info;VRR-$vrr;GT-$gt;TR-$strid";
                        $str_line .= ";MOTIF-$match_motif:$match_cn" if ($match_motif ne '');
                    }
                    else{
                        $annot = 'unknown';
                        $annot = "$match_motif:$match_cn" if ($match_motif ne '');
                        $str_line = "$chr\t$ipos\tINS\t$ilen\tRN-$info;VRR-$vrr;GT-$gt";
                        $str_line .= ";MOTIF-$match_motif:$match_cn" if ($match_motif ne '');
                    }
                    $str_line .= ";BP-$insbp" if (defined $insbp);
                    push @{$INS_seq{$ipos}}, "$chr:$ipos-$ilen $annot $ins_seq";
                }
            }
            elsif ($match_flag == 2){
                push @{$INS_seq{$pos}}, "$chr:$pos-$len $annot $ins_seq";
                $str_line = "$chr\t$pos\t$type\t$len\tRN-$info;RD-$read_cov;GT-$gt;TR-$strid";
                $str_line .= ";BP-$insbp" if (defined $insbp);
                $hit_strID{$strid} = $str_line;
            }
            elsif (($hit_strid ne '') and ($hit_strid ne $strid)){
                my ($pos2, $end2) = split (/=/, $STR2{$hit_strid});
                $annot = $hit_strid;
                if (exists $ins_strID{$hit_strid}){
                    my $str_line2 = "$chr\t$pos2\tTR-ins\t$len\tRN-$info;RD-$read_cov;GT-$gt;TRCN-$ins_cn;TR-$hit_strid";
                    $str_line2 .= ";BP-$insbp" if (defined $insbp);
                    my $ins_seq_info = "$chr:$pos2-$len $annot $ins_seq";
                    push @{$ins_strID_add{$hit_strid}}, "$str_line2==$pos2==$ins_seq_info";
                    next;
                }
                push @{$INS_seq{$pos2}}, "$chr:$pos2-$len $annot $ins_seq";
                $str_line = "$chr\t$pos2\t$type\t$len\tRN-$info;RD-$read_cov;GT-$gt;TRCN-$ins_cn;TR-$hit_strid";
                $str_line .= ";BP-$insbp" if (defined $insbp);
                $hit_strID{$hit_strid} = $str_line;
            }
            elsif ($hit_strid ne ''){
                push @{$INS_seq{$pos}}, "$chr:$pos-$len $annot $ins_seq";
                $str_line = "$chr\t$pos\t$type\t$len\tRN-$info;RD-$read_cov;GT-$gt;TRCN-$ins_cn;TR-$hit_strid";
                $str_line .= ";BP-$insbp" if (defined $insbp);
                $hit_strID{$hit_strid} = $str_line;
            }
        }
        else{
            $str_line = "$chr\t$pos\t$type\t$len\tRN-$info;RD-$read_cov;GT-$gt;TR-$strid";
            $str_line .= ";BP-$insbp" if (defined $insbp);
            $str_line .= ";TRDUP-$dpr" if ($dpr >= 1.2);
        }
        print OUTD "$str_line\n" if ($str_line ne '');
    }
}
close (FILE);

foreach my $sid (sort keys %ins_strID_add){
    my $count = 0;
    foreach (@{$ins_strID_add{$sid}}){
        my ($str_line, $pos2, $insseq_info) = split (/==/, $_);
        if (!exists $hit_strID{$sid}){
            $count ++;
            last if ($count == 3);
            if (@{$ins_strID_add{$sid}} == 1){
                print OUTD "$str_line\n";
                push @{$INS_seq{$pos2}}, $insseq_info;
            }
            else{
                my @strline = split (/\t/, $str_line);
                $strline[2] = 'TR-ins' . $count;
                $str_line = join ("\t", @strline);
                print OUTD "$str_line\n";
                push @{$INS_seq{$pos2}}, $insseq_info;
            }
        }
        else{
            my @strline1 = split (/\t/, $hit_strID{$sid});
            my @strline2 = split (/\t/, $str_line);
            my $len1 = $strline1[3];
            my $len2 = $strline2[3];
            next if ($len1 < $min_str_indel_size);
            next if ($len2 < $min_str_indel_size);
            my $flag = 0;
            if ($len1 >= $len2){
                $flag = 1 if ($len1 / $len2 > $str_max_len_rate);
            }
            else{
                $flag = 1 if ($len1 / $len2 < $str_min_len_rate);
            }
            if ($flag == 1){
                $strline2[2] = 'TR-ins2';
                $str_line = join ("\t", @strline2);
                print OUTD "$str_line\n";
                push @{$INS_seq{$pos2}}, $insseq_info;
                last;
            }
        }
    }
}

if (scalar keys %INS_result > 0){
    if ((scalar keys %INS_undef > 0) and ($intersperse_dup_find == 1)){
        my %qlen;
        my %qpos;
        open (INS2, "> $ins_undef_fasta");
        foreach my $pos2 (sort {$a <=> $b} keys %INS_undef){
            my ($header, $annot, $seq) = split (/\s+/, $INS_undef{$pos2});
            print INS2 ">$header\n$seq\n";
            $qlen{$header} = length $seq;
            $qpos{$header} = $pos2;
        }
        close (INS2);
        
        my $minimap_command = "minimap2 -ax sr -t1 -Y --MD $ref_file $ins_undef_fasta > $temp_dir/Z.$chr.sam 2>$temp_dir/Z.$chr.minimap.log";
        system ("$minimap_command");
        
        my %qalign;
        open (FILE2, "$temp_dir/Z.$chr.sam") or die "$temp_dir/Z.$chr.sam is not found: \n";
        while (my $line2 = <FILE2>){
            chomp $line2;
            next if ($line2 =~ /^\@/);
            my @line2 = split (/\t/, $line2);
            my $query = $line2[0];
            my $chr2 = $line2[2];
            my $pos2 = $line2[3];
            my $cigar = $line2[5];
            next if ($cigar eq '*');
            my $tag = $line2[1];
            my $mapQ = $line2[4];
            next if ($mapQ < $min_mapQ_idup);
            my $strand = 'F';
            my $sec_align = 0;
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
            next if ($sec_align == 1);
            my $maplen = 0;
            my $mismatch = 0;
            while ($cigar =~ /(\d+)([MIX=])/g){
                $maplen += $1;
                $mismatch += $1 if ($2 eq 'I');
            }
            next if ($maplen < $min_dup_len);
            my $qlen = $qlen{$query};
            my $map_rate = int ($maplen / $qlen * 1000 + 0.5) / 10 if ($qlen > 0);
            next if ($map_rate < $min_coverage);
            my $end2 = $pos2;
            while ($cigar =~ /(\d+)[MDX=]/g){
                $end2 += $1;
            }
            my $tlen = $end2 - $pos2 + 1;
            if ($line2 =~ /MD:Z:(\S+)/){
                my $MD = $1;
                $MD =~ s/[\d\^]+//g;
                $mismatch += length $MD;
            }
            my $mismatch_rate = int ($mismatch / $qlen * 1000 + 0.5) / 10;
            next if ($mismatch_rate > $max_mismatch2);
            next if ($mismatch_rate > $max_mismatch3) and ($qlen <= 150);
            $maplen -= $mismatch;
            if (!exists $qalign{$query}){
                $qalign{$query} = "$chr2\t$pos2\t$tlen\t$maplen\t$strand\t$qlen";
            }
            else{
                my ($chr3, $pos3, $tlen3, $mlen3) = split (/\t/, $qalign{$query});
                if ($maplen > $mlen3){
                    $qalign{$query} = "$chr2\t$pos2\t$tlen\t$maplen\t$strand\t$qlen";
                }
            }
        }
        close (FILE2);
        
        foreach my $qname (keys %qalign){
            my $qpos = $qpos{$qname};
            my ($chr2, $pos2, $tlen2, $maplen2, $strand2, $qlen2) = split (/\t/, $qalign{$qname});
            my $dist = $qlen2 * 0.5;
            $dist = 100 if ($qlen2 <= 150);
            if (($chr2 eq $chr) and ($pos2 >= $qpos - $qlen2 - $dist) and ($pos2 <= $qpos + $qlen2 + $dist)){
                if (exists $INS_seq{$qpos}){
                    my $qcount = 0;
                    foreach (@{$INS_seq{$qpos}}){
                        my ($header, $annot, $seq) = split (/\s+/, $_);
                        $annot =~ s/unknown/DUP/;
                        ${$INS_seq{$qpos}}[$qcount] = "$header $annot $seq";
                        $qcount ++;
                    }
                }
                if (exists $INS_result{$qpos}){
                    my @ins_line = split (/\t/, $INS_result{$qpos});
                    my $type2 = 'INS(DUP)';
                    $type2 = 'INS(DUP:R)' if ($strand2 eq 'R');
                    $ins_line[2] = $type2;
                    $ins_line[4] .= ";DUPPOS-$pos2,DUPLEN-$tlen2";
                    my $new_line = join ("\t", @ins_line);
                    $INS_result{$qpos} = $new_line;
                }
            }
            else{
                if (exists $INS_seq{$qpos}){
                    my $qcount = 0;
                    foreach (@{$INS_seq{$qpos}}){
                        my ($header, $annot, $seq) = split (/\s+/, $_);
                        $annot =~ s/unknown/intDUP/ if ($annot =~ /unknown/);
                        ${$INS_seq{$qpos}}[$qcount] = "$header $annot $seq"; 
                        $qcount ++;
                    }
                }
                if (exists $INS_result{$qpos}){
                    my @ins_line = split (/\t/, $INS_result{$qpos});
                    my $type2 = 'INS(intDUP)';
                    $type2 = 'INS(intDUP:R)' if ($strand2 eq 'R');
                    $ins_line[2] = $type2;
                    $ins_line[4] .= ";DUPPOS-$chr2:$pos2,DUPLEN-$tlen2";
                    my $new_line = join ("\t", @ins_line);
                    $INS_result{$qpos} = $new_line;
                }
            }
        }
    }
    foreach my $pos2 (sort {$a <=> $b} keys %INS_result){
        print OUTD "$INS_result{$pos2}\n";
    }
    %INS_undef = ();
    %INS_result = ();
}
close (OUTD);

if (scalar keys %INS_seq > 0){
    open (INS, "> $ins_out");
    foreach my $pos2 (sort {$a <=> $b} keys %INS_seq){
        foreach (@{$INS_seq{$pos2}}){
            my ($header, $annot, $seq) = split (/\s+/, $_);
            print INS ">$header $annot\n$seq\n" if ($seq ne '');
        }
    }
    close (INS);
}

undef %INS_seq;
undef %inv_seq;
undef %ins_pos_seq;
undef %ins_bp_seq;
undef %ins_str_seq;
undef %seq;


sub find_overlap{
    my ($clip3t_seq, $clip5t_seq, $chr2, $pos) = @_;
    my ($match, $insseq) = &yass_align_overlap ($clip3t_seq, $clip5t_seq, $chr2, $pos);
    if ($match == 1){
        return (1, $insseq);
    }
    else{
        return (0, '');
    }
}

sub find_repeat{
    my ($read_seq, $chr2, $pos) = @_;
    my $seq_len = length $read_seq;
    my $seq_fasta = "$temp_dir/Z.$chr2-seq.fa";
    open (NEWFILE, "> $seq_fasta");
    print NEWFILE '>seq', "\n", $read_seq;
    close (NEWFILE);
    my @trf_result = `trf $seq_fasta 2 7 7 80 10 50 500 -ngs -h`;
    my $motif_copy = 0;
    my $motif_size = 0;
    my $motif = '';
    my $identity = 0;
    my %rep;
    if (@trf_result > 0){
        foreach my $result (@trf_result){
            chomp $result;
            my @result = split (/\s+/, $result);
            next if (@result < 15);
            $motif_size = $result[2];
            $motif_copy = $result[3];
            $identity = $result[5];
            $motif = $result[13];
            my $rep_len = $motif_size * $motif_copy;
            my $rep_coverage = int ($rep_len / $seq_len * 1000 + 0.5) / 10;
            if (($identity >= $min_str_identity) and ($rep_coverage >= $min_coverage)){
                $rep{$motif} = "$motif_copy=$identity=$rep_coverage";
            }
        }
    }
    if (scalar keys %rep == 0){
        return ('', 0, 0, 0);
    }
    else{
        foreach my $motif (sort {length $a <=> length $b} keys %rep){
            my ($motif_copy, $identity, $rep_coverage) = split (/=/, $rep{$motif});
            return ($motif, $motif_copy, $identity, $rep_coverage);
            last;
        }
    }
}

sub find_dup_1{             # find tandem DUP with an ins-containing read
    my ($read_seq, $tag, $chr2, $pos) = @_;
    my $ilen = length $read_seq;
    my $upos1 = $pos - $ilen - int ($ilen * 0.2) if ($tag eq 'tdup');
    $upos1 = $pos - $ilen - int ($ilen * 0.5) if ($tag eq 'sdup');
    $upos1 = 1 if ($upos1 <= 0);
    my $dpos2 = $pos + $ilen + int ($ilen * 0.2) if ($tag eq 'tdup');
    $dpos2 = $pos + $ilen + int ($ilen * 0.5) if ($tag eq 'sdup');
    $dpos2 = $chr_len{$chr} if ($dpos2 > $chr_len{$chr});
    my $ref_len = $dpos2 - $upos1 + 1;
    my $ref_seq = '';
    my $divnum1 = int ($upos1 / $Mbin_size);
    my $res1 = $upos1 % $Mbin_size;
    my $divnum2 = int ($dpos2 / $Mbin_size);
    my $res2 = $dpos2 % $Mbin_size;
    if ($divnum1 == $divnum2){
        if ($res1 > 0){
            $ref_seq = substr ($seq{$divnum1}, $res1 - 1, $ref_len);
        }
        else{
            $ref_seq = substr ($seq{$divnum1}, 0, $ref_len - 1);
        }
    }
    else{
        if ($divnum1 + 1 == $divnum2){
            my $subseq1 = '';
            $subseq1 = substr ($seq{$divnum1}, $res1 - 1) if ($res1 > 0);
            my $subseq2 = substr ($seq{$divnum2}, 0, $res2);
            $ref_seq = $subseq1 . $subseq2;
        }
        else{
            my $subseq1 = '';
            $subseq1 = substr ($seq{$divnum1}, $res1 - 1) if ($res1 > 0);
            my $subseq2 = substr ($seq{$divnum2}, 0, $res2);
            my $submid = '';
            my $count = $divnum1 + 1;
            while ($count < $divnum2){
                $submid .= $seq{$count};
                $count ++;
            }
            $ref_seq = $subseq1 . $submid . $subseq2;
        }
    }
    my ($dup_start, $dup_end, $direction) = &yass_align_dup ($ref_seq, $read_seq, $chr2);
    if ($dup_start > 0){
        $dup_start += $upos1;
        $dup_end += $upos1;
        return ($dup_start, $dup_end, $direction);
    }
    else{
        return (0, 0, 'F');
    }
}

sub find_dup_2{             # find tandem DUP with an ins-containing read
    my ($read_seq, $tag, $chr2, $pos) = @_;
    my $ilen = length $read_seq;
    my $upos1 = 0;
    my $dpos2 = 0;
    $upos1 = $pos if ($tag eq 'T3');
    $upos1 = $pos - $ilen - 100 if ($tag eq 'T5');
    $upos1 = 1 if ($upos1 <= 0);
    $dpos2 = $pos if ($tag eq 'T5');
    $dpos2 = $pos + $ilen + 100 if ($tag eq 'T3');
    $dpos2 = $chr_len{$chr} if ($dpos2 > $chr_len{$chr2});
    my $ref_len = $dpos2 - $upos1 + 1;
    my $ref_seq = '';
    my $divnum1 = int ($upos1 / $Mbin_size);
    my $res1 = $upos1 % $Mbin_size;
    my $divnum2 = int ($dpos2 / $Mbin_size);
    my $res2 = $dpos2 % $Mbin_size;
    if ($divnum1 == $divnum2){
        if ($res1 > 0){
            $ref_seq = substr ($seq{$divnum1}, $res1 - 1, $ref_len);
        }
        else{
            $ref_seq = substr ($seq{$divnum1}, 0, $ref_len - 1);
        }
    }
    else{
        if ($divnum1 + 1 == $divnum2){
            my $subseq1 = '';
            $subseq1 = substr ($seq{$divnum1}, $res1 - 1) if ($res1 > 0);
            my $subseq2 = substr ($seq{$divnum2}, 0, $res2);
            $ref_seq = $subseq1 . $subseq2;
        }
        else{
            my $subseq1 = '';
            $subseq1 = substr ($seq{$divnum1}, $res1 - 1) if ($res1 > 0);
            my $subseq2 = substr ($seq{$divnum2}, 0, $res2);
            my $submid = '';
            my $count = $divnum1 + 1;
            while ($count < $divnum2){
                $submid .= $seq{$count};
                $count ++;
            }
            $ref_seq = $subseq1 . $submid . $subseq2;
        }
    }
    my ($dup_start, $dup_end, $direction) = &yass_align_dup ($ref_seq, $read_seq, $chr2);
    if ($dup_start > 0){
        return (1);
    }
    else{
        return (0);
    }
}

sub find_me{
    my ($read_seq, $chr2, $pos, $sid) = @_;
    my $ilen = length $read_seq;
    my %match;
    my %match2;
    my $skip_te = '';
    $skip_te = $STR_ME{$sid} if (exists $STR_ME{$sid});
    if ($skip_te ne ''){
        my $str_motif_size = length $STR_motif{$sid};
        my $te_len = $TE_len{$skip_te};
        if ($str_motif_size < $te_len * 0.5){
            $skip_te = '';
        }
    }
    foreach my $me (@TE){
        next if ($me eq $skip_te);
        my $me_len = $TE_len{$me};
        next if ($ilen < 150) and ($me_len >= 1000);
#        next if ($ilen > $me_len * 3);
        my ($overlap, $cn, $te_cov, $query_cov) = &yass_align_me ($me, $read_seq, $chr, $pos);
        if ($overlap > 0){
            $match{$me} = $query_cov;
            $match2{$me} = "$overlap=$cn=$query_cov";
            last if ($overlap / $ilen >= 0.9);
        }
    }
    if (scalar keys %match > 0){
        my $te_type = '';
        my $te_type2 = '';
        foreach my $me (sort {$match{$b} <=> $match{$a}} keys %match){
            if ($te_type eq ''){
                $te_type = $me;
                next;
            }
            if ($te_type2 eq ''){
                $te_type2 = $me;
                last;
            }
        }
        my ($overlaplen, $te_cn, $te_cov) = split (/=/, $match2{$te_type});
        my ($overlaplen2, $te_cn2, $te_dcov2) = split (/=/, $match2{$te_type2}) if ($te_type2 ne '');
        my $te_type2_len = "$te_type2-$overlaplen2" if ($te_type2 ne '');
        if (($ilen < 400) and ($te_type2 eq 'ALU') and ($overlaplen2 > 100)){
            $te_type2_len = "$te_type-$overlaplen";
            $te_type = $te_type2;
            $te_cn = $te_cn2;
            $overlaplen = $overlaplen2;
        }
        return ($te_type, $overlaplen, $te_cn, '') if ($te_type2 eq '');
        return ($te_type, $overlaplen, $te_cn, $te_type2_len) if ($te_type2 ne '');
    }
    else{
        return ('', 0, 0, '');
    }
}

sub check_inv{             # check whether the reverse-complemented sequences of 5'- and 3'-clipped sequences of a double-clipped read are the downstream and upstream reference flanking sequences of the INV BPs, respectively
    my ($ref_info, $chr2, $pos, $len) = @_;
    my $select_id1 = '';
    my $select_id2 = '';

    my $count = 0;
    my $match_flag = 0;
    foreach my $readid (keys %{$inv_seq{$pos}}){
        my $clipseq1 = '';
        my $clipseq2 = '';
        my $ref_seq1 = '';
        my $ref_seq2 = '';
        my $ref_start1 = 0;
        my $ref_end1 = 0;
        my $ref_start2 = 0;
        my $ref_end2 = 0;
        $count ++;
        (my $tag1, $clipseq1) = split (/=/, ${${$inv_seq{$pos}}{$readid}}{1}) if (exists ${${$inv_seq{$pos}}{$readid}}{1});
        (my $tag2, $clipseq2) = split (/=/, ${${$inv_seq{$pos}}{$readid}}{2}) if (exists ${${$inv_seq{$pos}}{$readid}}{2});
        next if ($tag1 ne $tag2);
        next if (length $clipseq1 < 200) and (length $clipseq2 < 200);

        if ($tag1 eq 'bp1'){
            $ref_start1 = $pos;
            $ref_end1 = $pos + $len - 1;
            $ref_start2 = $pos - length $clipseq2;
            $ref_end2 = $pos - 1;
        }
        elsif ($tag1 eq 'bp2'){
            $ref_start1 = $pos + $len;
            $ref_end1 = $ref_start1 + length ($clipseq1) - 1;
            $ref_start2 = $pos;
            $ref_end2 = $pos + $len - 1;
        }
        my $divnum1 = int ($ref_start1 / $Mbin_size);
        my $res1 = $ref_start1 % $Mbin_size;
        my $divnum2 = int ($ref_end1 / $Mbin_size);
        my $res2 = $ref_end1 % $Mbin_size;
        if ($divnum1 == $divnum2){
            if ($res1 > 0){
                $ref_seq1 = substr ($seq{$divnum1}, $res1 - 1, $len);
            }
            else{
                $ref_seq1 = substr ($seq{$divnum1}, 0, $len);
            }
        }
        else{
            my $subseq1 = substr ($seq{$divnum1}, $res1 - 1) if ($res1 > 0);
            $subseq1 = substr ($seq{$divnum1}, 0) if ($res1 == 0);
            my $subseq2 = substr ($seq{$divnum2}, 0, $res2);
            $ref_seq1 = $subseq1 . $subseq2;
        }
        $divnum1 = int ($ref_start2 / $Mbin_size);
        $res1 = $ref_start2 % $Mbin_size;
        $divnum2 = int ($ref_end2 / $Mbin_size);
        $res2 = $ref_end2 % $Mbin_size;
        if ($divnum1 == $divnum2){
            if ($res1 > 0){
                $ref_seq2 = substr ($seq{$divnum1}, $res1 - 1, $len);
            }
            else{
                $ref_seq2 = substr ($seq{$divnum1}, 0, $len);
            }
        }
        else{
            my $subseq1 = substr ($seq{$divnum1}, $res1 - 1) if ($res1 > 0);
            $subseq1 = substr ($seq{$divnum1}, 0) if ($res1 == 0);
            my $subseq2 = substr ($seq{$divnum2}, 0, $res2);
            $ref_seq2 = $subseq1 . $subseq2;
        }
        $clipseq1 = reverse $clipseq1;
        $clipseq1 =~ tr/ACGT/TGCA/;
        $clipseq2 = reverse $clipseq2;
        $clipseq2 =~ tr/ACGT/TGCA/;
        my $match_tag1 = 0;
        my $match_tag2 = 0;
        ($match_tag1) = &yass_align_inv ($ref_seq1, $clipseq1, $chr) if ($clipseq1 ne '');
        ($match_tag2) = &yass_align_inv ($ref_seq2, $clipseq2, $chr) if ($clipseq2 ne '');
        if (($match_tag1 == 1) and ($match_tag2 == 1)){
            $match_flag = 1;
            last;
        }
        last if ($count == 3);
    }
    
    return ($match_flag);
}

sub yass_align_dup{
    my ($seq1, $seq2, $chr2) = @_;
    open (NEWFILE1, "> $temp_dir/Z.$chr2-seq1.fa");
    print NEWFILE1 '>seq1', "\n", $seq1;
    close (NEWFILE1);
                
    open (NEWFILE1, "> $temp_dir/Z.$chr2-seq2.fa");
    print NEWFILE1 '>seq2', "\n", $seq2;
    close (NEWFILE1);
    my @result = `yass -O 10 -m $max_mismatch $temp_dir/Z.$chr2-seq1.fa $temp_dir/Z.$chr2-seq2.fa 2>/dev/null`;
    my $seq2_len = length $seq2;
    my $overlap_len = 0;
    my $match_pos1 = 0;
    my $match_pos2 = 0;
    my %match_pos;
    my $match_flag = 0;
    my $direction = 'F';
    my $mismatch_count = 0;
    my $pre_end = 0;
    foreach my $line (@result){
        chomp $line;
        if (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match_flag == 0)){
            $match_flag = 1;
            if ($1 < $2){
                $match_pos1 = $1;
                $match_pos2 = $2;
            }
            else{
                $match_pos1 = $2;
                $match_pos2 = $1;
                $direction = 'R';
            }
            $match_pos{$match_pos1} = $match_pos2;
            $pre_end = $match_pos2;
        }
        elsif (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match_flag >= 1)){
            my $dir2 = 'F';
            if ($1 < $2){
                $match_pos1 = $1;
                $match_pos2 = $2;
            }
            else{
                $match_pos1 = $2;
                $match_pos2 = $1;
                $dir2 = 'R';
            }
            if ((!exists $match_pos{$match_pos1}) and ($direction eq $dir2)){
                $match_pos{$match_pos1} = $match_pos2;
            }
            if ($match_pos1 < $pre_end){
                $match_flag = 2;
            }
            else{
                $match_flag = 1;
            }
            $pre_end = $match_pos2 if ($match_pos2 > $pre_end);
        }
        elsif ((($match_flag == 1)) and ($line =~ /^[\|\:\.\s]+$/)){
            my $match_count = 0;
            $match_count ++ while ($line =~ /\|/g);
            my $bases = length $line;
            $mismatch_count += $bases - $match_count;
            if ($line =~ /(\s{10,})/){
                my $indel_len = length ($1) - 1;
                $mismatch_count -= $indel_len;
            }
        }
    }
    if ($match_flag >= 1){
        my $pre_end2 = 0;
        foreach my $pos1 (sort {$a <=> $b} keys %match_pos){
            my $pos2 = $match_pos{$pos1};
            my $matchlen = $pos2 - $pos1 + 1;
            if ($pos2 <= $pre_end2){
                next;
            }
            elsif ($pos1 < $pre_end2){
                $matchlen = $pos2 - $pre_end2 + 1;
            }
            $overlap_len += $matchlen;
            $pre_end2 = $pos2;
        }
        my $coverage = int ($overlap_len / $seq2_len * 1000 + 0.5) / 10;
        my $mmrate = int ($mismatch_count / $overlap_len * 1000 + 0.5) / 10;
        if (($mmrate <= $max_mismatch2) and ($coverage >= $min_coverage)){
            return ($match_pos1, $match_pos2, $direction);
        }
        else{
            return (0, 0, 'F');
        }
    }
    else{
        return (0, 0, 'F');
    }
}

sub yass_align_me{
    my ($me, $insseq, $chr2, $pos) = @_;
    if (!-f "$temp_dir/chr$chr.$me.fa"){
        my $me_seq = $TE{$me};
        open (NEWFILE1, "> $temp_dir/Z.$chr2.$me.fa");
        print NEWFILE1 ">$me\n$me_seq\n";
        close (NEWFILE1);
    }        
    open (NEWFILE1, "> $temp_dir/Z.$chr2-seq.fa");
    print NEWFILE1 ">seq\n$insseq\n";
    close (NEWFILE1);
    my @result = `yass -O 10 -m $max_mismatch $temp_dir/Z.$chr2.$me.fa $temp_dir/Z.$chr2-seq.fa 2>/dev/null`;
    my $me_len = $TE_len{$me};
    my $ins_len = length $insseq;
    my $overlap_tlen = 0;
    my $overlap_qlen = 0;
    my $match_tpos1 = 0;
    my $match_tpos2 = 0;
    my $match_qpos1 = 0;
    my $match_qpos2 = 0;
    my %match_tpos;
    my %match_qpos;
    my $match_flag = 0;
    my $t_direction = 'F';
    my $q_direction = 'F';
    my $mismatch_count = 0;
    my $min_me_coverage2 = $min_me_coverage;
    $min_me_coverage2 = 50 if ($ins_len < 150);
    my $cn = 0;
    my $pre_tend = 0;
    foreach my $line (@result){
        chomp $line;
        if (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match_flag == 0)){
            $match_flag = 1;
            if ($1 < $2){
                $match_tpos1 = $1;
                $match_tpos2 = $2;
            }
            else{
                $match_tpos1 = $2;
                $match_tpos2 = $1;
                $t_direction = 'R';
            }
            if ($3 < $4){
                $match_qpos1 = $3;
                $match_qpos2 = $4;
            }
            else{
                $match_qpos1 = $4;
                $match_qpos2 = $3;
                $q_direction = 'R';
            }
            $match_tpos{$match_tpos1} = $match_tpos2;
            $match_qpos{$match_qpos1} = $match_qpos2;
            $pre_tend = $match_tpos2;
        }
        elsif (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match_flag >= 1)){
            my $tdir2 = 'F';
            my $qdir2 = 'F';
            if ($1 < $2){
                $match_tpos1 = $1;
                $match_tpos2 = $2;
            }
            else{
                $match_tpos1 = $2;
                $match_tpos2 = $1;
                $tdir2 = 'R';
            }
            if ($3 < $4){
                $match_qpos1 = $3;
                $match_qpos2 = $4;
            }
            else{
                $match_qpos1 = $4;
                $match_qpos2 = $3;
                $qdir2 = 'R';
            }
            if ((!exists $match_tpos{$match_tpos1}) and ($t_direction eq $tdir2)){
                $match_tpos{$match_tpos1} = $match_tpos2;
            }
            if ((!exists $match_qpos{$match_qpos1}) and ($q_direction eq $qdir2)){
                $match_qpos{$match_qpos1} = $match_qpos2;
            }
            if ($match_tpos1 < $pre_tend){
                $match_flag = 2;
            }
            else{
                $match_flag = 1;
            }
            $pre_tend = $match_tpos2 if ($match_tpos2 > $pre_tend);
        }
        elsif ((($match_flag == 1)) and ($line =~ /^[\|\:\.\s]+$/)){
            my $match_count = 0;
            $match_count ++ while ($line =~ /\|/g);
            my $bases = length $line;
            $mismatch_count += $bases - $match_count;
            if ($line =~ /(\s{10,})/){
                my $indel_len = length ($1) - 1;
                $mismatch_count -= $indel_len;
            }
        }
    }
    if ($match_flag >= 1){
        my $pre_end = 0;
        foreach my $pos1 (sort {$a <=> $b} keys %match_tpos){
            my $pos2 = $match_tpos{$pos1};
            my $matchlen = $pos2 - $pos1 + 1;
            if ($pos2 <= $pre_end){
                next;
            }
            elsif ($pos1 < $pre_end){
                $matchlen = $pos2 - $pre_end + 1;
            }
            $overlap_tlen += $matchlen;
            $pre_end = $pos2;
        }
        $pre_end = 0;
        foreach my $pos1 (sort {$a <=> $b} keys %match_qpos){
            my $pos2 = $match_qpos{$pos1};
            my $matchlen = $pos2 - $pos1 + 1;
            if ($pos2 <= $pre_end){
                next;
            }
            elsif ($pos1 < $pre_end){
                $matchlen = $pos2 - $pre_end + 1;
            }
            $overlap_qlen += $matchlen;
            $pre_end = $pos2;
        }
        my $mmrate = int ($mismatch_count / $overlap_tlen * 1000 + 0.5) / 10;
        my $coverage_query = int ($overlap_qlen / $ins_len * 1000 + 0.5) / 10;
        my $coverage_me = int ($overlap_tlen / $me_len * 1000 + 0.5) / 10;
        if (($mmrate <= $max_mismatch) and (($coverage_query >= $min_me_coverage2) or ($coverage_me >= $min_me_coverage2))){
            $cn = int ($overlap_qlen / $me_len * 10 + 0.5) / 10;
            return ($overlap_qlen, $cn, $coverage_me, $coverage_query);
        }
        else{
            return (0, 0, 0, 0);
        }
    }
    else{
        return (0, 0, 0, 0);
    }
}

sub yass_align_overlap{
    my ($seq1, $seq2, $chr2, $pos) = @_;
    open (NEWFILE1, "> $temp_dir/Z.$chr2-seq1.fa");
    print NEWFILE1 '>seq1', "\n", $seq1;
    close (NEWFILE1);
                
    open (NEWFILE1, "> $temp_dir/Z.$chr2-seq2.fa");
    print NEWFILE1 '>seq2', "\n", $seq2;
    close (NEWFILE1);
    my @result = `yass -O 10 -m $max_mismatch -r 0 $temp_dir/Z.$chr2-seq1.fa $temp_dir/Z.$chr2-seq2.fa 2>/dev/null`;
    my $tlen = length $seq1;
    my $qlen = length $seq2;
    my $overlap_len = 0;
    my $tpos1 = 0;
    my $tpos2 = 0;
    my $qpos1 = 0;
    my $qpos2 = 0;
    my $match = 0;
    my $indel_count = 0;
    my $insseq = '';
    foreach my $line (@result){
        chomp $line;
        if (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match == 0)){
            next if ($1 > $2);
            my $match_len1 = $2 - $1 + 1;
            my $match_len2 = $4 - $3 + 1;
            $overlap_len = $match_len1;
            $overlap_len = $match_len2 if ($match_len1 < $match_len2);
            next if ($overlap_len < 30);
            $tpos1 = $1;
            $tpos2 = $2;
            $qpos1 = $3;
            $qpos2 = $4;
            if (($qpos1 <= 50) and ($tpos2 >= $tlen - 50)){
                my $subseq = substr ($seq1, 0, $tpos1 - 1);
                $insseq = $subseq . $seq2;
            }
            elsif (($qpos1 <= 50) and ($tpos2 >= $tlen - 50)){
                $insseq = substr ($seq2, 0, $qpos2 - 1);
            }
            elsif (($qpos1 <= 50) and ($qpos2 >= $qlen - 50)){
                $insseq = substr ($seq1, 0, $tpos2 - 1);
            }
            elsif (($tpos1 <= 50) and ($tpos2 >= $tlen - 50)){
                my $ilen = $qlen - $qpos1 + 1;
                $insseq = substr ($seq2, $qpos1 - 1, $ilen);
            }
            elsif (($qpos2 >= $qlen - 50) and ($tpos1 <= 50)){
                $insseq = substr ($seq1, 0, $tpos2 - 1);
            }
            $match = 1 if ($insseq ne '');
        }
        elsif (($match == 1) and ($line =~ /^[ACGTN\-]+$/)){
            $indel_count ++ while ($line =~ /-/g);
        }
        if (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match == 1)){
            last;
        }
    }
    if (($match == 1) and ($overlap_len >= 50)){
        my $indelrate = int ($indel_count / $overlap_len * 1000 + 0.5) / 10;
        return (0, $insseq);
    }
    else{
        return (0, '');
    }
}

sub yass_align_inv{
    my ($seq1, $seq2, $chr2) = @_;
    open (NEWFILE1, "> $temp_dir/Z.$chr2-seq1.fa");
    print NEWFILE1 '>seq1', "\n", $seq1;
    close (NEWFILE1);
                
    open (NEWFILE1, "> $temp_dir/Z.$chr2-seq2.fa");
    print NEWFILE1 '>seq2', "\n", $seq2;
    close (NEWFILE1);
    my @result = `yass -O 10 -m $max_mismatch -r 0 $temp_dir/Z.$chr2-seq1.fa $temp_dir/Z.$chr2-seq2.fa 2>/dev/null`;
    my $seq1_len = length $seq1;
    my $seq2_len = length $seq2;
    my $overlap_len = 0;
    my $match_pos1 = 0;
    my $match_pos2 = 0;
    my $match_flag = 0;
    my %match_pos;
    my $insseq = '';
    my $pre_end = 0;
    my $mismatch_count = 0;
    foreach my $line (@result){
        chomp $line;
#print STDERR "$line\n";# if ($pos == 64637275) or ($pos == 68357399); 
        if (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match_flag == 0)){
            $match_flag = 1;
            if ($1 < $2){
                $match_pos1 = $1;
                $match_pos2 = $2;
            }
            else{
                $match_pos1 = $2;
                $match_pos2 = $1;
            }
            $match_pos{$match_pos1} = $match_pos2;
            $pre_end = $match_pos2;
        }
        elsif (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match_flag >= 1)){
            if ($1 < $2){
                $match_pos1 = $1;
                $match_pos2 = $2;
            }
            else{
                $match_pos1 = $2;
                $match_pos2 = $1;
            }
            if (!exists $match_pos{$match_pos1}){
                $match_pos{$match_pos1} = $match_pos2;
            }
            if ($match_pos1 < $pre_end){
                $match_flag = 2;
            }
            else{
                $match_flag = 1;
            }
            $pre_end = $match_pos2 if ($match_pos2 > $pre_end);
        }
        elsif ((($match_flag == 1)) and ($line =~ /^[\|\:\.\s]+$/)){
            my $match_count = 0;
            $match_count ++ while ($line =~ /\|/g);
            my $bases = length $line;
            $mismatch_count += $bases - $match_count;
            if ($line =~ /(\s{10,})/){
                my $indel_len = length ($1) - 1;
                $mismatch_count -= $indel_len;
            }
        }
    }
    if ($match_flag >= 1){
        my $pre_end = 0;
        foreach my $pos1 (sort {$a <=> $b} keys %match_pos){
            my $pos2 = $match_pos{$pos1};
            my $matchlen = $pos2 - $pos1 + 1;
            if ($pos2 <= $pre_end){
                next;
            }
            elsif ($pos1 < $pre_end){
                $matchlen = $pos2 - $pre_end + 1;
            }
            $overlap_len += $matchlen;
            $pre_end = $pos2;
        }
        my $mmrate = int ($mismatch_count / $overlap_len * 1000 + 0.5) / 10;
        my $coverage = int ($overlap_len / $seq2_len * 1000 + 0.5) / 10;
        if (($mmrate <= $max_mismatch2) and ($coverage >= 0.7)){
            return (1, $overlap_len);
        }
        else{
            return (0, 0);
        }
    }
    else{
        return (0, 0);
    }
}

sub yass_realign{
    my ($seq1, $seq2, $chr2) = @_;
    open (NEWFILE1, "> $temp_dir/Z.$chr2-seq1.fa");
    print NEWFILE1 '>seq1', "\n", $seq1;
    close (NEWFILE1);
                
    open (NEWFILE1, "> $temp_dir/Z.$chr2-seq2.fa");
    print NEWFILE1 '>seq2', "\n", $seq2;
    close (NEWFILE1);
    my $indel_rate2 = 15;
    my @result = `yass -O 10 -m $max_mismatch -i $indel_rate2 -r 0 $temp_dir/Z.$chr2-seq1.fa $temp_dir/Z.$chr2-seq2.fa 2>/dev/null`;
    my $seq1_len = length $seq1;
    my $seq2_len = length $seq2;
    my $overlap_len = 0;
    my $merge_seq = '';
    my $match = 0;
    my $indel_count = 0;
    foreach my $line (@result){
        chomp $line;
        if (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match == 0)){
            my $match_len1 = $2 - $1 + 1;
            my $match_len2 = $4 - $3 + 1;
            $overlap_len = $match_len1;
            $overlap_len = $match_len2 if ($match_len1 < $match_len2);
            if (($overlap_len >= $seq1_len * 0.8) or (($2 >= $seq1_len - 50) and ($1 <= 50))){
                $match = 1;
            }
            else{
                next;
            }
        }
        elsif (($match == 1) and ($line =~ /^[ACGTN\-]+$/)){
            $indel_count ++ while ($line =~ /-/g);
        }
        elsif (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match == 1)){
            last;
        }
    }
    if ($match == 1){
        my $indelrate = int ($indel_count / $overlap_len * 1000 + 0.5) / 10;
        if ($indelrate > $indel_rate2){
            $match = 0;
        }
    }
    return ($match);
}

sub yass_realign2{
    my ($seq1, $seq2, $chr2) = @_;
    open (NEWFILE1, "> $temp_dir/Z.$chr2-seq1.fa");
    print NEWFILE1 '>seq1', "\n", $seq1;
    close (NEWFILE1);
                
    open (NEWFILE1, "> $temp_dir/Z.$chr2-seq2.fa");
    print NEWFILE1 '>seq2', "\n", $seq2;
    close (NEWFILE1);
    my @result = `yass -O 20 -m 10 -i 10 $temp_dir/Z.$chr2-seq1.fa $temp_dir/Z.$chr2-seq2.fa 2>/dev/null`;
    my $seq1_len = length $seq1;
    my $seq2_len = length $seq2;
    my $overlap_len = 0;
    my $merge_seq = '';
    my $match = 0;
    my $mismatch = 0;
    my $indel_count = 0;
    foreach my $line (@result){
        chomp $line;
        if (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match == 0)){
            my $match_len1 = $2 - $1 + 1;
            my $match_len2 = $4 - $3 + 1;
            $overlap_len = $match_len1;
            $overlap_len = $match_len2 if ($match_len1 < $match_len2);
            next if ($overlap_len < 20);
            $match = 1;
        }
        elsif (($match == 1) and ($line =~ /ts\s+:\s+(\d+)\s+tv\s+:\s+(\d+)/)){
            $mismatch = $1 + $2;
        }
        elsif (($match == 1) and ($line =~ /^[ACGTN\-]+$/)){
            $indel_count ++ while ($line =~ /-/g);
        }
        elsif (($line =~ /\*\((\d+)-(\d+)\)\((\d+)-(\d+)\)/) and ($match == 1)){
            last;
        }
    }
    if ($match == 1){
        my $match_len = $overlap_len - $mismatch - $indel_count;
        my $len1 = length $seq1;
        my $len2 = length $seq2;
        my $identity = int ($match_len / $overlap_len * 1000 + 0.5) / 10;
        if ($identity < 85){
            $match = 0;
        }
        elsif (($overlap_len / $len1 < 0.85) and ($overlap_len / $len2 < 0.85)){
            $match = 0;
        }
    }
    return ($match);
}

sub motif_test{
    my ($insseq, $motif) = @_;
    my @matchlen;
    my $mlen = length $motif;
    my ($match) = &motif_test_sub ($insseq, $motif);
    push @matchlen, $match;
    if ($mlen > 1){
        my $motif2 = $motif;
        while ($mlen > 1){
            my $base1 = substr ($motif2, 0, 1, '');
            $motif2 .= $base1;
            my ($match) = &motif_test_sub ($insseq, $motif2);
            push @matchlen, $match;
            $mlen --;
        }
    }
    @matchlen = sort {$b <=> $a} @matchlen;
    return ($matchlen[0]);
}

sub motif_test_sub{
    my ($insseq, $motif) = @_;
    my $mlen = length $motif;
    my $ins_len = length $insseq;
    my $cn = 0;
    my @base = ('A', 'C', 'G', 'T');
    $cn = $insseq =~ s/$motif/X/g;
    $cn = 0 if (!defined $cn);
    my $cn2 = 0;
    my $ident = 0;
    if ($mlen > 1){
        my @motifs;
        my @mbases = split (//, $motif);
        for (my $i = 0; $i < $mlen; $i++){
            foreach my $base (@base){
                next if ($base eq $mbases[$i]);
                my $motif2 = $motif;
                substr ($motif2, $i, 1, $base);
                push @motifs, $motif2;
            }
        }
        foreach my $motif1 (@motifs){
            my $cn_2 = $insseq =~ s/$motif1/XX/g;
            if (defined $cn_2){
                $cn2 += $cn_2;
            }
        }
    }
    my $mlen2 = $mlen - 1;
    my $alt_mrate = 0;
    $alt_mrate = $cn2 * ($mlen2 / $mlen) if ($mlen > 0);
    my $sum_cn = $cn + $cn2;
    $ident = int (($cn + $alt_mrate) / $sum_cn * 1000 + 0.5) / 10 if ($sum_cn > 0);
    my $match_len = $sum_cn * $mlen;
    $match_len = 0 if ($ident < $min_str_identity);
    return ($match_len);
}

sub TRF_test1{
    my ($insseq, $chr2) = @_;
    my $ins_len = length $insseq;
    my $min_score = 50;
    if ($ins_len < 30){
        $min_score = 10;
    }
    elsif ($ins_len < 50){
        $min_score = 20;
    }
    elsif ($ins_len < 100){
        $min_score = 30;
    }
    my $test_fasta = "$temp_dir/TRF.$chr2.test.fa";
    open (OUT, "> $test_fasta");
    print OUT ">test\n";
    print OUT "$insseq\n";
    close (OUT);
    my $top_motif = '';
    my $top_len = 0;
    my $top_cn = 0;
    my $top_pos = 0;
    my $top_end = 0;
    my @result = `trf $test_fasta 2 7 7 80 10 $min_score 2000 -d -l 1 -h -ngs`;
    if (@result > 0){
        foreach (@result){
            chomp $_;
            next if ($_ =~ /^\@/);
            my @line = split (/\s+/, $_);
            my $apos = $line[0];
            my $aend = $line[1];
            my $amlen = $line[2];
            my $amotif = $line[13];
            my $alen = $aend - $apos + 1;
            if ($alen > $top_len){
                $top_len = $alen;
                $top_motif = $amotif;
                $top_pos = $apos;
                $top_end = $aend;
            }
        }
        my $motif_size = length $top_motif;
        if ($top_pos - 1 < $motif_size){
            $top_len += $top_pos - 1;
        }
        if ($ins_len - $top_end < $motif_size){
            $top_len += $ins_len - $top_end;
        }
        $top_cn = int ($top_len / $motif_size * 10 + 0.5) / 10 if ($motif_size > 0);
    }
    if ($top_len >= $ins_len * $min_str_len_rate){
        return ($top_cn, $top_motif);
    }
    else{
        return (0, $top_motif);
    }
}

sub TRF_test2{
    my ($insseq, $motif, $chr2) = @_;
    my $mlen = length $motif;
    my $ins_len = length $insseq;
    my $min_score = 50;
    if ($ins_len < 30){
        $min_score = 10;
    }
    elsif ($ins_len < 50){
        $min_score = 20;
    }
    elsif ($ins_len < 100){
        $min_score = 30;
    }
    my $test_fasta = "$temp_dir/TRF.$chr2.test.fa";
    open (OUT, "> $test_fasta");
    print OUT ">test\n";
    print OUT "$insseq\n";
    close (OUT);
    my $match_len = 0;
    my $match_motif = '';
    my $top_motif = '';
    my $match_pos = 0;
    my $match_end = 0;
    my $top_len = 0;
    my $top_pos = 0;
    my $top_end = 0;
    my $sec_len = 0;
    my $sec_motif = '';
    my $motif2 = $motif . $motif;
    my @result = `trf $test_fasta 2 7 7 80 10 $min_score 2000 -d -l 1 -h -ngs`;
    if (@result > 0){
        foreach (@result){
            chomp $_;
            next if ($_ =~ /^\@/);
            my @line = split (/\s+/, $_);
            my $apos = $line[0];
            my $aend = $line[1];
            my $amlen = $line[2];
            my $amotif = $line[13];
            my $alen = $aend - $apos + 1;
            if (($amlen >= $mlen * 0.8) and ($amlen <= $mlen * 1.25)){
                $match_len = $alen;
                $match_motif = $amotif;
                $match_pos = $apos;
                $match_end = $aend;
                last;
            }
        }
        my $motif_size = length $match_motif;
        if ($match_len > 0){
            if ($match_pos - 1 < $motif_size){
                $match_len += $match_pos - 1;
            }
            if ($ins_len - $match_end < $motif_size){
                $match_len += $ins_len - $match_end;
            }
        }
        $match_len = $ins_len if ($match_len > $ins_len);
    }
    return ($match_len, $match_motif);
}

sub multalin_ins1{
    my ($ins_seq, $motif, $chr2) = @_;
    $ins_seq =~ s/-//g if ($ins_seq =~ /-/);
    my $mlen = length $motif;
    my $ins_len = length $ins_seq;
    my $motif2 = $motif . $motif;
    my $input_fasta = "$temp_dir/Z.$chr2.multalin.fasta";
    my $input_base = "$temp_dir/Z.$chr2.multalin";
    my $multalin_out = "$input_base.msf";
    my $multalin_out1 = "$input_base.clu";
    my $multalin_out2 = "$input_base.cl2";
    my $multalin_log = "$input_base.log";
    my @ins_seq;
    push @ins_seq, $ins_seq;
    my $match_len = 0;
    my $mismatch_len = 0;
    my $match_len2 = 0;
    my $mismatch_len2 = 0;
    my $cov_len = 0;
    open (OUT1, "> $input_fasta");
    print OUT1 ">STR\n";
    print OUT1 "$motif2\n";
    print OUT1 ">INS\n";
    print OUT1 "$ins_seq\n";
    close (OUT1);
    system ("multalin -q $input_fasta > $multalin_log");
    
    if (!-f $multalin_out){
        return (0, 0, 0);
    }
    my $qalign = '';
    my $talign = '';
    my $ins_cons = '';
    open (FILE1, $multalin_out);
    while (my $line = <FILE1>){
        chomp $line;
        if ($line =~ /^\s+Consensus\s+(.+)/){
            my $subseq = $1;
            $subseq =~ s/\s+//g;
            $ins_cons .= $subseq;
        }
        elsif ($line =~ /^\s+INS\s+(.+)/){
            my $subseq = $1;
            $subseq =~ s/\s+//g;
            $qalign .= $subseq;
        }
        elsif ($line =~ /^\s+STR\s+(.+)/){
            my $subseq = $1;
            $subseq =~ s/\s+//g;
            $talign .= $subseq;
        }
    }
    close (FILE1);
    
    my $match = 0;
    my $mismatch = 0;
    my $cov = 0;
    $qalign =~ s/^\.*//;
    $qalign =~ s/\.*$//;
    $talign =~ s/^\.*//;
    $talign =~ s/\.*$//;
    next if ($qalign eq '');
    my @bases = split (//, $talign);
    foreach (@bases){
        if ($_ eq '.'){
            $mismatch ++;
        }
    }
    my @qbases = split (//, $qalign);
    foreach (@qbases){
        if ($_ eq '.'){
            $mismatch ++;
        }
    }
    $ins_cons =~ s/^\.*//;
    $ins_cons =~ s/\.*$//;
    my @bases2 = split (//, $ins_cons);
    foreach (@bases2){
        if ($_ eq '.'){
        }
        elsif ($_ =~ /[acgt]/){
            $mismatch ++;
            $cov ++;
        }
        else{
            $match ++;
            $cov ++;
        }
    }
    $match_len += $match;
    $mismatch_len += $mismatch;
    my $match_rate = 0;
    $match_rate = int ($match / ($match + $mismatch) * 1000 + 0.5) / 10 if ($match + $mismatch > 0);
    if ($match_rate >= $min_str_identity){
        $cov_len += $cov;
        $match_len2 += $match;
        $mismatch_len2 += $mismatch;
    }
    system ("rm -f $input_fasta $multalin_out $multalin_out1 $multalin_out2 $multalin_log");
    my $coverage = int ($cov_len / $ins_len * 1000 + 0.5) / 10;

    return ($match_rate, $coverage, $cov_len);
}

sub multalin_ins2{
    my ($ins_seq, $str_seq, $chr2) = @_;
    $ins_seq =~ s/-//g if ($ins_seq =~ /-/);
    my $ins_len = length $ins_seq;
    my $str_len = length $str_seq;
    my $input_fasta = "$temp_dir/Z.$chr2.multalin.fasta";
    my $input_base = "$temp_dir/Z.$chr2.multalin";
    my $multalin_out = "$input_base.msf";
    my $multalin_out1 = "$input_base.clu";
    my $multalin_out2 = "$input_base.cl2";
    my $multalin_log = "$input_base.log";
    
    my $match_len = 0;
    my $mismatch_len = 0;
    my $cov_len = 0;
    my $match_pos = 0;
    my $match_end = 0;
    my $count = 0;
    my $step_size = int ($str_len * 1.6);
    my $first_cov = 0;
    my $last_3unmatch = 0;
    while ($match_end < $ins_len - $str_len * 0.5){
        $count ++;
        my $iseq = '';
        $iseq = substr ($ins_seq, $match_end, $step_size) if ($match_end < $ins_len - $step_size);
        $iseq = substr ($ins_seq, $match_end) if ($match_end >= $ins_len - $step_size);
        open (OUT1, "> $input_fasta");
        print OUT1 ">INS\n";
        print OUT1 "$iseq\n";
        print OUT1 ">STR\n";
        print OUT1 "$str_seq\n";
        close (OUT1);
        system ("multalin -q $input_fasta > $multalin_log");
        if (!-f $multalin_out){
            $match_end += $step_size;
            next;
        }
        my $qalign = '';
        my $talign = '';
        my $ins_cons = '';
        open (FILE1, $multalin_out);
        while (my $line = <FILE1>){
            chomp $line;
            if ($line =~ /^\s+Consensus\s+(.+)/){
                my $subseq = $1;
                $subseq =~ s/\s+//g;
                $ins_cons .= $subseq;
            }
            elsif ($line =~ /^\s+INS\s+(.+)/){
                my $subseq = $1;
                $subseq =~ s/\s+//g;
                $talign .= $subseq;
            }
            elsif ($line =~ /^\s+STR\s+(.+)/){
                my $subseq = $1;
                $subseq =~ s/\s+//g;
                $qalign .= $subseq;
            }
        }
        close (FILE1);
        system ("rm -f $input_fasta $multalin_out $multalin_out1 $multalin_out2 $multalin_log");
        my $match = 0;
        my $mismatch = 0;
        my $cov = 0;
        my $unmatch5 = '';
        my $unmatch3 = '';
        $unmatch5 = $1 if ($qalign =~ /^(\.+)/);
        $unmatch3 = $1 if ($qalign =~ /(\.+)$/);
        if ($count == 1){
            $match_pos = 1 + length ($unmatch5);
        }
        my $unmatch3_len = length $unmatch3;
        $last_3unmatch = $unmatch3_len;
        $match_end = $match_end + $step_size - $unmatch3_len;
        $qalign =~ s/^\.*//;
        $qalign =~ s/\.*$//;
        $talign =~ s/^\.*//;
        $talign =~ s/\.*$//;
        next if ($talign eq '');
        my @bases = split (//, $qalign);
        foreach (@bases){
            if ($_ eq '.'){
                $mismatch ++;
            }
        }
        my @tbases = split (//, $talign);
        foreach (@tbases){
            if ($_ eq '.'){
                $mismatch ++;
            }
        }
        $ins_cons =~ s/^\.*//;
        $ins_cons =~ s/\.*$//;
        my @bases2 = split (//, $ins_cons);
        foreach (@bases2){
            if ($_ eq '.'){
            }
            elsif ($_ =~ /[acgt]/){
                $mismatch ++;
                $cov ++;
            }
            else{
                $match ++;
                $cov ++;
            }
        }
        my $match_rate = 0;
        $match_rate = int ($match / ($match + $mismatch) * 1000 + 0.5) / 10 if ($match + $mismatch > 0);
        if ($match_rate >= $min_str_identity - 10){
            if ($cov / $str_len >= $min_str_len_rate){
                $cov = $str_len;
            }
            $cov_len += $cov;
            $match_len += $match;
            $mismatch_len += $mismatch;
        }
        $first_cov = $cov if ($count == 1);
    }
    my $identity = 0;
    $identity = int ($match_len / ($match_len + $mismatch_len) * 1000 + 0.5) / 10 if ($match_len + $mismatch_len > 0);
    $cov_len += $match_pos if ($first_cov < $str_len);
    $cov_len += $last_3unmatch if ($last_3unmatch < $str_len);
    $cov_len = $ins_len if ($cov_len > $ins_len);
    my $coverage = int ($cov_len / $ins_len * 1000 + 0.5) / 10;

    return ($identity, $coverage, $cov_len);
}

