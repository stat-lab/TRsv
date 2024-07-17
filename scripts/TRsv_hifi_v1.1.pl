#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use threads;
use FindBin qw($Bin);

# detect SVs from error-corrected long read alignment bam file
# define repeats in INS sequence with rtf
# when an INS flank sequence is an TR and the INS is repeat sequence, confirm whether the TR and the INS repeat motif are identical
# for non-repeat INS, the INS sequence is aligned to the flanking regions (up tp 100 Kb) with yass to define the INS to be DUP
# call INSs with single BPs (removing candidates whose 5'- and 3'-clipped sequences are aligned to the flanking ref sequences upstream and downstream of the BP)
# call DELs in a repeat-aware manner

my $data_dir = "$Bin/../Data";

my $ref_file = '';

my $simple_repeat = '';

my $simple_repeat_unanalyzed = '';

my $lowconf_TR = '';

my $exclude_bed = '';

my $TE_fasta = '';

my $gap_bed = '';

my $samtool_path = '';
my $yass_path = '';
my $trf_path = '';
my $multalin_path = '';

print "# $0 @ARGV\n";

my $bam_file = '';

my $cores = 1;

my $out_prefix = '';

my $chr_analyzed = 'ALL';
my $exclude_chr = '';

my $min_mapQ = 1;
my $min_mapQ_idup = 20;

my $max_SAR = 0.7;

my $min_ins_reads = 2;
my $min_del_reads = 2;
my $min_str_reads = 3;

my $min_VRR = 0.05;
my $min_str_vrr = 0.15;
my $min_str_vrr2 = 0.3;

my $str_max_len_rate = 1.5;

my $max_depth_fold = 15;

my $min_maplen = 500;

my $chromium = 0;

my $include_secalign = 0;

my $max_mismatch = 15;

my $min_indel_size = 50;
my $min_str_indel_size = 0;
my $min_ins_str_mei = 200;
my $min_str_len_rate = 0.5;
my $min_str_cn = 0.1;
my $min_str_identity = 70;

my $indel_rate = 10;

my $min_coverage = 80;
my $min_me_coverage = 60;

my $dup_find = 1;

my $intersperse_dup_find = 1;

my $mei_find = 1;

my $bp_diff = 100;
my $bp_diff2 = 200;
my $bp_diff3 = 150;
my $ins_bp_diff = 50;
my $max_bp_diff = 300;
my $min_overlap_rate = 0.5;
my $max_dist = 10;
my $max_dist2 = 15;
my $max_match_size = 200;
my $min_clip_len = 50;
my $max_insbp_diff = 200;
my $min_inv_size = 100;
my $max_del_size = 50000000;
my $max_inv_size = 5000000;
my $max_dup_len = 10000000;
my $min_dup_len = 50;
my $max_bp_sd_small_ins = 5;

my $Mbin_size = 1000000;
my $Mbin_size2 = 100000;

my $platform = 'hifi';

my $skip_stage1 = 0;

my $targeted_seq = 0;

my $non_human = 0;
my $build = 37;

my $skip_step = 0;

my $conf_file = '';

my $help;

GetOptions(
    'bam_file|b=s' => \$bam_file,
    'ref_file|r=s' => \$ref_file,
    'thread|n=i' => \$cores,
    'prefix|p=s' => \$out_prefix,
    'conf|c=s' => \$conf_file,
    'chr|tc=s' => \$chr_analyzed,
    'excl_chr|xc=s' => \$exclude_chr,
    'gap_bed|gb=s' => \$gap_bed,
    'repeat_bed|rep=s' => \$simple_repeat,
    'repeat_u|reu=s' => \$simple_repeat_unanalyzed,
    'lowconf_tr|lcs=s' => \$lowconf_TR,
    'exclude_bed|exc=s' => \$exclude_bed,
    'te_fasta|tf=s' => \$TE_fasta,
    'min_len|ml=i' => \$min_indel_size,
    'min_tr_len|msl=i' => \$min_str_indel_size,
    'min_tr_lrate|mslr=f' => \$min_str_len_rate,
    'min_tr_cn|msc=f' => \$min_str_cn,
    'min_tr_ident|msi=i' => \$min_str_identity,
    'min_ins_read|mir=i' => \$min_ins_reads,
    'min_del_read|mdr=i' => \$min_del_reads,
    'min_tr_read|msr=i' => \$min_str_reads,
    'min_vrr|mv=f' => \$min_VRR,
    'min_tr_vrr|msv=f' => \$min_str_vrr,
    'max_tr_rate|xsr=i' => \$str_max_len_rate,
    'max_dpf|xd=f' => \$max_depth_fold,
    'min_mapq|mq=i' => \$min_mapQ,
    'max_sar|sar=f' => \$max_SAR,
    'max_mismatch|xm=i' => \$max_mismatch,
    'max_indel|xi=i' => \$indel_rate,
    'incl_sec|insec' => \$include_secalign,
    'platform|x=s' => \$platform,
    'skip|sk=i' => \$skip_step,
    'targeted|t' => \$targeted_seq,
    'non_human|nh' => \$non_human,
    'build=s' => \$build,
    'samtool_path|sp=s' => \$samtool_path,
    'trf_path|tp=s' => \$trf_path,
    'yass_path|yp=s' => $yass_path,
    'multalin_path|mp=s' => \$multalin_path,
    'help|h' => \$help,
    
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  TRsv_hifi_v1.1.pl <option>

  Options:
   --bam_file or -b <STR>       bam file for long reads [mandatory]
   --ref_file or -r <STR>       reference fasta file [mandatory]
   --repeat_bed or -rep <STR>   simple/short repeat (STR/VNTR) bed file (TRsv-compatible format), which is automatically selected for human. [mandatory for non-human]
   --prefix or -p <STR>         output prefix [mandatory]

   --conf or -c <STR>           configure file specifying options. The specified options in conf file override separately specified options. [optional]
   --repeat_u or -reu <STR>     simple/short repeat (STR/VNTR) bed file of TRs (i.e., > 10 Kb TRs) not ppresent in the file specified with --repeat, which is automatically selected for human. [optional]
   --lowconf_tr or -lcs <STR>   list file of TR IDs in the TR file specified with --repeat_u whose TR regions could cause low-confident TR-CNV calling, which is automatically selected for human.
                                TR-CNVs within these regions are marked 'LowConf' in the FILTER field of output vcf file. [optional]
   --exclude_bed or -exc <STR>  bed file indicating chromosomal regions to be excluded [optional]
                                [default for human: centromere regions, which is automatically selected for human. Specify any non-file character if you want to include centromeres.]
   --gap_bed or -gb <STR>       gap bed file for reference gap regions, which is automatically selected for human. [optional]
   --te_fasta or -tf <STR>      retroelement sequence fasta file, which is automatically selected for human. [optional]
   
   --chr or -tc <STR>           target chromosome(s) (comma-separated chromosome name(s)) [default: ALL]
   --excl_chr or -xc <STR>      chromosome(s) to be excluded for analysis (comma-separated)
   --platform or -x <STR>       sequencing platform (hifi|ont|clr)
   --min_len or -ml <INT>       minimum size of INS/DEL [default: 50]
   --min_tr_len or -msl  <INT>  minimum size of TR-CNV [default: 50 (N specified with -ml)]
   --min_vrr or -mv <FLOAT>     minimum VRR (variant rate to read depth) [default: 0.05]
   --min_tr_vrr or -msv <FLOAT> minimum VRR (variant rate to read depth) for TR-CNVs [default: 0.15]
   --min_mapq or -mq <INT>      minimum mapping quality [default: 1]
   --incl_sec or -insec <BOOLEAN> include secondary alignments [default: false]
   --skip or -sk <INT>          skip step-1~3 (1: skip step-1, 2: skip step-2, 3: skip step-3, 12: skip steps-1 and -2, 123: skip steps-1, -2, and -3) [default: 0]
   --non_human or -nh <BOOLEAN> sample is a non-human species [default: false]
   --build <STR>                reference build (37|38|T2T) when using human sample [default: 37]
   --thread or -n <INT>         number of threads [default: 1]

   --samtool_path or -sp <STR>  path of samtools (${samtool_path}/samtools) if the corresponding path is not set in $PATH
   --trf_path or -tp <STR>      path of trf (${trf_path}/trf) if the corresponding path is not set in $PATH
   --yass_path or -yp <STR>     path of yass (${yass_path}/yass) if the corresponding path is not set in $PATH
   --multalin_path or -mp <STR> path of multalin (${multalin_path}/multalin) if the corresponding path is not set in $PATH

   --min_ins_read or -mir <INT> minimum number of reads supporting INSs/DUPs/INVs [default: 2]
   --min_del_read or -mdr <INT> minimum number of reads supporting DELs [default: 2]
   --min_tr_read or -msr <INT>  minimum number of reads supporting TR-CNVs [default: 3]
   --max_tr_rate or -xsr <FLOAT> maximum rate of length for descriminating different TR-CNV alleles at a site [deafult: 1.5, when the ratio of two different TR-CNV length is larger than 1.5 or smaller than 1/1.5, these are considered two alleles]
   --min_tr_lrate or -mslr <FLOAT> minimum rate of TR unit content in TR-INS sequence [default: 0.5] (TR-INS that a content of TR unit/motif is smaller than the specified rate is regarded as a normal INS)
                                If 0 is specified, all INSs within a TR region are regarded as TR-INS without checking the homology between the INS sequence and TR motif
   --min_tr_cn or -msc <FLOAT>  minimum copy number of TR-CNV [default: 0.1] (TR-CNV with a smaller copy number than the specified value and smaller than the value specified with --min_len is discarded)
   --min_tr_ident or -msi <INT> minimum identity (%) of TR-INS sequence with TR motif [default: 70]
   --max_dpf or -xd <FLOAT>     maximum fold of mean read depth to consider SV calling (not call in high-depth regions with specified value x mean depth) [default: 15]
   --max_sar or -sar <FLOAT>    maximum rate of SVs supported by alignments with maping quality 0, including secondary alignmnents (SVs exceeding this value are marked as 'LowQual' in the FILTER field) [default: 0.7]
   --max_mismatch or -xm <INT>  maximum percentage of mismatch for yass aligner to search INS homology to TE or flanking regions [default: 15]
   --max_indel or -xi <INT>     maximum percentage of indel for yass aligner to search INS homology to TE or flanking regions [default: 10]
   --targeted or -t <BOOLEAN>   the data is targeted sequencing data [default: false]
   --help or -h                 output help message
   
=cut

if ($conf_file ne ''){
    open (FILE, $conf_file) or die "$conf_file is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#|^$/);
        $line = $1 if ($line =~ /(.+?)#/);
        my ($arg, $value) = split (/\s+/, $line);
        next if ($value eq '') or ($value eq '#') or ($value eq '-');
        if ($arg eq 'bam_file'){
            $bam_file = $value;
        }
        elsif ($arg eq 'ref_file'){
            $ref_file = $value;
        }
        elsif ($arg eq 'repeat_bed'){
            $simple_repeat = $value;
        }
        elsif ($arg eq 'repeat_u'){
            $simple_repeat_unanalyzed = $value;
        }
        elsif ($arg eq 'lowconf_tr'){
            $lowconf_TR = $value;
        }
        elsif ($arg eq 'exclude_bed'){
            $exclude_bed = $value;
        }
        elsif ($arg eq 'te_fasta'){
            $TE_fasta = $value;
        }
        elsif ($arg eq 'gap_bed'){
            $gap_bed = $value;
        }
        elsif ($arg eq 'prefix'){
            $out_prefix = $value;
        }
        elsif ($arg eq 'min_len'){
            $min_indel_size = $value;
        }
        elsif ($arg eq 'min_tr_len'){
            $min_str_indel_size = $value;
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
        elsif ($arg eq 'incl_sec'){
            $include_secalign = $value;
        }
        elsif ($arg eq 'non_human'){
            $non_human = $value;
        }
        elsif ($arg eq 'build'){
            $build = $value;
        }
        elsif ($arg eq 'skip'){
            $skip_step = $value;
        }
        elsif ($arg eq 'thread'){
            $cores = $value;
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
        elsif ($arg eq 'min_tr_lrate'){
            $min_str_len_rate = $value;
        }
        elsif ($arg eq 'min_tr_cn'){
            $min_str_cn = $value;
        }
        elsif ($arg eq 'min_tr_ident'){
            $min_str_identity = $value;
        }
        elsif ($arg eq 'max_tr_rate'){
            $str_max_len_rate = $value;
        }
        elsif ($arg eq 'max_tr_rate'){
            $str_max_len_rate = $value;
        }
        elsif ($arg eq 'max_dpf'){
            $max_depth_fold = $value;
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
        elsif ($arg eq 'targeted'){
            $targeted_seq = $value;
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
}

die "-bam_file option not specified:\n" if ($bam_file eq '');
die "-ref option not specified:\n" if ($ref_file eq '');
die "-prefix option not specified:\n" if ($out_prefix eq '');

my $temp_dir = "$out_prefix.temp";
system ("mkdir $temp_dir") if (!-d $temp_dir);

if ($samtool_path eq ''){
    my $Spath = `which samtools`;
    chomp $Spath;
    if (($Spath =~ /\s/) or (!-f $Spath)){
        die "samtools path is not specified with --samtool_path or not in PATH:\n";
    }
}
else{
    if (!-f "$samtool_path/samtools"){
        die "samtools does not exist in the specified path: $samtool_path:\n";
    }
    $ENV{PATH} = "$samtool_path:" . $ENV{PATH};
}
if ($trf_path eq ''){
    my $Tpath = `which trf`;
    chomp $Tpath;
    if (($Tpath =~ /\s/) or (!-f $Tpath)){
        die "trf path is not specified with --trf_path or not in PATH:\n";
    }
}
else{
    if (!-f "$samtool_path/samtools"){
        die "trf does not exist in the specified path: $trf_path:\n";
    }
    $ENV{PATH} = "$trf_path:" . $ENV{PATH};
}
if ($yass_path eq ''){
    my $Ypath = `which yass`;
    chomp $Ypath;
    if (($Ypath =~ /\s/) or (!-f $Ypath)){
        die "yass path is not specified with --yass_path or not in PATH:\n";
    }
}
else{
    if (!-f "$yass_path/yass"){
        die "yass does not exist in the specified path: $yass_path:\n";
    }
    $ENV{PATH} = "$yass_path:" . $ENV{PATH};
}
if ($multalin_path eq ''){
    my $Mpath = `which multalin`;
    chomp $Mpath;
    if (($Mpath =~ /\s/) or (!-f $Mpath)){
        die "multalin path is not specified with --multalin_path or not in PATH:\n";
    }
}
else{
    if (!-f "$multalin_path/yass"){
        die "multalin does not exist in the specified path: $multalin_path:\n";
    }
    $ENV{PATH} = "$multalin_path:" . $ENV{PATH};
}

if ($non_human == 0){
    $TE_fasta = "$data_dir/TE.fa";
    if ($build eq '37'){
        $simple_repeat = "$data_dir/simpleRepeat.b37.20-10000.addHip.ME.bed.gz" if ($simple_repeat eq '');
        $simple_repeat_unanalyzed = "$data_dir/simpleRepeat.b37.ov10K.bed.gz" if ($simple_repeat_unanalyzed eq '');
        $lowconf_TR = "$data_dir/Low-confidence_TR-ov10Kb_region_b37.txt" if ($lowconf_TR eq '');
        $gap_bed = "$data_dir/gap.bed" if ($gap_bed eq '');
    }
    if ($build eq '38'){
        $simple_repeat = "$data_dir/simpleRepeat.b38.20-10000.addHip.ME.bed.gz" if ($simple_repeat eq '');
        $simple_repeat_unanalyzed = "$data_dir/simpleRepeat.b38.ov10K.bed.gz" if ($simple_repeat_unanalyzed eq '');
        $lowconf_TR = "$data_dir/Low-confidence_TR-ov10Kb_region_b38.txt" if ($lowconf_TR eq '');
        $exclude_bed = "$data_dir/hg38.centromere.bed" if ($exclude_bed eq '');
        $gap_bed = "$data_dir/gap.b38.bed" if ($gap_bed eq '');
    }
    elsif ($build eq 'T2T'){
        $simple_repeat = "$data_dir/simpleRepeat.t2t-chm13.20-10000.ME.bed.gz" if ($simple_repeat eq '');
        $simple_repeat_unanalyzed = "$data_dir/simpleRepeat.t2t-chm13.v2.0.ov10K.bed.gz" if ($simple_repeat_unanalyzed eq '');
        $lowconf_TR = "$data_dir/Low-confidence_TR-ov10Kb_region_T2T.txt" if ($lowconf_TR eq '');
        $exclude_bed = "$data_dir/chm13.v2.0.centromere.bed" if ($exclude_bed eq '');
    }
}
else{
    die "repeat_bed file is not specified or not present:\n" if ($simple_repeat eq '') or (!-f $simple_repeat);
}

if ($min_str_indel_size == 0){
    $min_str_indel_size = $min_indel_size;
}

my $min_VRR2 = $min_VRR * 2;
$min_VRR2 = 0.01 if ($min_VRR2 > 0.01);

my $ref_base = $ref_file;
$ref_base = $1 if ($ref_file =~ /\/(.+?)$/);
$ref_base = $1 if ($ref_base =~ /(.+)\.gz$/);
if ($ref_file =~ /(.+)\.gz$/){
    system ("gzip -dc $ref_file > $ref_base");
    $ref_file = $ref_base;
}
my $ref_index = "$ref_file.fai";
if (!-f $ref_index){
    if (!-f $ref_base){
        system ("ln -s $ref_file");
    }
    system ("samtools faidx $ref_base");
    $ref_index = "$ref_base.fai";
    $ref_file = $ref_base;
}

my $config_file = "$out_prefix.config.txt";
open (OUT, "> $config_file");
print OUT "bam_file\t$bam_file\n";
print OUT "ref_file\t$ref_file\n";
print OUT "repeat_bed\t$simple_repeat\n";
print OUT "te_fasta\t$TE_fasta\n";
print OUT "gap_bed\t$gap_bed\n";
print OUT "out_prefix\t$out_prefix\n";
print OUT "min_len\t$min_indel_size\n";
print OUT "min_tr_len\t$min_str_indel_size\n";
print OUT "min_tr_lrate\t$min_str_len_rate\n";
print OUT "min_tr_cn\t$min_str_cn\n";
print OUT "min_tr_ident\t$min_str_identity\n";
print OUT "min_ins_read\t$min_ins_reads\n";
print OUT "min_del_read\t$min_del_reads\n";
print OUT "min_tr_read\t$min_str_reads\n";
print OUT "max_tr_rate\t$str_max_len_rate\n";
print OUT "min_vrr\t$min_VRR\n";
print OUT "min_tr_vrr\t$min_str_vrr\n";
print OUT "max_dpf\t$max_depth_fold\n";
print OUT "min_mapq\t$min_mapQ\n";
print OUT "max_sar\t$max_SAR\n";
print OUT "mismatch_rate\t$max_mismatch\n";
print OUT "indel_rate\t$indel_rate\n";
print OUT "incl_sec\t$include_secalign\n";
print OUT "targeted\t$targeted_seq\n";
print OUT "exclude\t$exclude_bed\n";
print OUT "non_human\t$non_human\n";
print OUT "samtool_path\t$samtool_path\n" if ($samtool_path ne '');
print OUT "trf_path\t$trf_path\n" if ($trf_path ne '');
print OUT "yass_path\t$yass_path\n" if ($yass_path ne '');
print OUT "multalin_path\t$multalin_path\n" if ($multalin_path ne '');

close (OUT);

my %target_chr;
my @tchr = split (/,/, $chr_analyzed) if ($chr_analyzed ne 'ALL');
map{$target_chr{$_} = 1} @tchr;

my %excl_chr;
my @xchr = split (/,/, $exclude_chr) if ($exclude_chr ne '');
map{$excl_chr{$_} = 1} @xchr;


my @chr;
open (FILE, $ref_index) or die "$ref_index is not found:$!\n";
while (my $line = <FILE>){
    chomp $line;
    my ($chr) = split (/\t/, $line);
    if ($chr_analyzed eq 'ALL'){
        next if ($non_human == 0) and ($chr !~ /^c*h*r*[\dXY]+$/);
    }
    elsif (!exists $target_chr{$chr}){
        next;
    }
    if (exists $excl_chr{$chr}){
        next;
    }
    push @chr, $chr;
}
close (FILE);

my %out_files;

my $count = 0;
my @jobs;
foreach my $chr (@chr){
    my $out_file = "$temp_dir/$out_prefix.chr$chr.discov.txt";
    $out_file = "$temp_dir/$out_prefix.$chr.discov.txt" if ($chr !~ /^[\dXY]+$/);
    if ((-f $out_file) and ($skip_step =~ /1/)){
        $out_files{$chr} = $out_file;
        next;
    }
    $count ++;
    my ($thread_t) = threads->new(\&run_step1, $config_file, $chr);
    push @jobs, $thread_t;
    if ($count == $cores){
        foreach (@jobs){
            my ($file1, $chr2, $align_info) = $_->join;
print STDERR "Finished 1st step: $chr2\n";
            $out_files{$chr2} = $file1;
        }
        $count = 0;
        undef @jobs;
    }
}
if (@jobs > 0){
    foreach (@jobs){
        my ($file1, $chr2, $align_info) = $_->join;
print STDERR "Finished 1st step: $chr2\n";
        $out_files{$chr2} = $file1;
    }
    $count = 0;
    undef @jobs;
}

my %out_files_2;
my @jobs2;
my $count2 = 0;

foreach my $chr (@chr){
    next if (!exists $out_files{$chr});
    my $out_chr = $out_files{$chr};

    my $out_file = "$temp_dir/$out_prefix.chr$chr.discov.out";
    $out_file = "$temp_dir/$out_prefix.$chr.discov.out" if ($chr !~ /^[\dXY]+$/);
    if ((-f $out_file) and ($skip_step =~ /2/)){
        $out_files_2{$chr} = $out_file;
        next;
    }
    if (!-f $out_chr){
        print STDERR "Warning: $out_chr from 1st step is not found:\n";
        next;
    }

    $count2 ++;
    my ($thread_t) = threads->new(\&run_step2, $config_file, $chr, $out_chr);
    push @jobs2, $thread_t;
    if ($count2 == $cores){
        foreach (@jobs2){
            my ($file1, $chr2) = $_->join;
            $out_files_2{$chr2} = $file1;
print STDERR "Finished 2nd step: $chr2\n";
        }
        $count2 = 0;
        undef @jobs2;
    }
}
if (@jobs2 > 0){
    foreach (@jobs2){
        my ($file1, $chr2) = $_->join;
        $out_files_2{$chr2} = $file1;
print STDERR "Finished 2nd step: $chr2\n";
    }
    $count2 = 0;
    undef @jobs2;
}

my %out_files_3;
my @jobs3;
my $count3 = 0;

foreach my $chr (@chr){
    next if (!exists $out_files_2{$chr});
    my $out_chr = $out_files_2{$chr};

    my $out_file = "$temp_dir/$out_prefix.chr$chr.discov.out";
    $out_file = "$temp_dir/$out_prefix.$chr.discov.out" if ($chr !~ /^[\dXY]+$/);
    if ((-f $out_file) and ($skip_step =~ /3/)){
        $out_files_3{$chr} = $out_file;
        next;
    }
    if (!-f $out_chr){
        print STDERR "Warning: $out_chr from 2nd step is not found:\n";
        next;
    }

    $count3 ++;
    my ($thread_t) = threads->new(\&run_step3, $config_file, $chr, $out_chr);
    push @jobs3, $thread_t;
    if ($count3 == $cores){
        foreach (@jobs3){
            my ($file1, $chr2) = $_->join;
            $out_files_3{$chr2} = $file1;
print STDERR "Finished 3rd step: $chr2\n";
        }
        $count3 = 0;
        undef @jobs3;
    }
}
if (@jobs3 > 0){
    foreach (@jobs3){
        my ($file1, $chr2) = $_->join;
        $out_files_3{$chr2} = $file1;
print STDERR "Finished 3rd step: $chr2\n";
    }
    $count3 = 0;
    undef @jobs3;
}

my %eval_call;
my %sv;
my %INS;
my $file_str = '';

foreach my $chr (@chr){
    next if (!exists $out_files_3{$chr});
    my $out_chr = $out_files_3{$chr};
    my $low_hom_ins = "$temp_dir/low-hom-INS.$chr.txt";
    $file_str .= "$low_hom_ins " if (-f $low_hom_ins);
    open (FILE, $out_chr) or die "$out_chr is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        my ($chr2, $pos, $type, $len, $info) = split (/\t/, $line);
print STDERR "$chr\n" if (!defined $info);
        my $type2 = $type;
        $type2 =~ s/-BP$// if ($type2 =~ /-BP$/);
        $type2 =~ s/-BP2$// if ($type2 =~ /-BP2$/);
        $type2 = 'INS' if ($type =~ /INS/);
        my $tag = 'IT';
        $tag = 'BP' if ($type eq 'DUP') or ($type eq 'INV') or ($type =~ /-BP/);
        $tag = 'DUP' if ($type eq 'INS(DUP)');
        $tag = 'DUP:R' if ($type eq 'INS(DUP:R)');
        $tag = 'intDUP' if ($type eq 'INS(intDUP)');
        $tag = 'intDUP:R' if ($type eq 'INS(intDUP:R)');
        $tag = 'TR' if ($type =~ /TR/);
        next if ($type2 eq 'DEL') and ($len > $max_del_size);
        next if ($type2 eq 'DUP') and ($len > $max_dup_len);
        next if (($type2 eq 'INV') or ($type2 eq 'REP')) and ($len > $max_inv_size);
        if (!exists ${${$eval_call{$type2}}{$chr}}{$pos}){
            ${${$eval_call{$type2}}{$chr}}{$pos} = "$tag=$len=$info";
        }
        else{
            my ($pre_tag) = split (/=/, ${${$eval_call{$type2}}{$chr}}{$pos});
            if (($pre_tag eq 'IT') and ($tag eq 'BP')){
            }
            elsif (($type =~ /INS\(DUP/) or ($tag eq 'BP')){
                ${${$eval_call{$type2}}{$chr}}{$pos} = "$tag=$len=$info";
            }
        }
    }
    close (FILE);
    my $out_ins = $out_chr;
    $out_ins =~ s/\.out$//;
    $out_ins .= '.INS.fa';
    if (-f $out_ins){
        my $header = '';
        open (FILE, $out_ins);
        while (my $line = <FILE>){
            chomp $line;
            if ($line =~ /^>(.+)/){
                $header = $1;
            }
            else{
                my ($chrposlen) = split (/\s+/, $header);
                my ($chrpos, $len) = split (/-/, $chrposlen);
                my ($chr, $pos) = split (/:/, $chrpos);
                my $chr02d = $chr;
                $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
                if (!exists ${$INS{$chr02d}}{$pos}){
                    ${$INS{$chr02d}}{$pos} = ">$header\n$line\n";
                }
                else{
                    ${$INS{$chr02d}}{$pos} .= ">$header\n$line\n";
                }
            }
        }
        close (FILE);
    }
}
if ($file_str ne ''){
    system ("cat $file_str > $temp_dir/low-hom-INS.txt");
    system ("rm $file_str") if (!-z "$temp_dir/low-hom-INS.txt");
}

my %STR;
my %STR2;
my %repeat;
my %uSTR;
my %uSTR2;
my %uSTR_id;
my %lowconf_STR;
my %lowconf_region;

open (FILE, $simple_repeat) or die "$simple_repeat is not found: $!\n" if ($simple_repeat !~ /\.gz$/);
open (FILE, "gzip -dc $simple_repeat |") or die "$simple_repeat is not found: $!\n" if ($simple_repeat =~ /\.gz$/);
while (my $line = <FILE>){
    chomp $line;
    next if ($line =~ /^#|^$/);
    my @line = split (/\t/, $line);
    my $chr = $line[0];
    my $pos = $line[1];
    my $end = $line[2];
    my $id = $line[3];
    my $me = $line[6];
    my $motif = $line[7];
    my $motif_size = $line[4];
    ${$STR{$chr}}{$pos} = $id;
    $STR2{$id} = "$pos=$end=$motif_size";
    my $Mbin_start = int ($pos / $Mbin_size);
    my $Mbin_end = int ($end / $Mbin_size);
    if ($Mbin_start == $Mbin_end){
        ${${$repeat{$chr}}{$Mbin_start}}{$pos} = $end;
    }
    else{
        my $Mbin = $Mbin_start;
        while ($Mbin <= $Mbin_end){
            ${${$repeat{$chr}}{$Mbin}}{$pos} = $end;
            $Mbin ++;
        }
    }
}
close (FILE);

if (-f $simple_repeat_unanalyzed){
    open (FILE, $simple_repeat_unanalyzed) or die "$simple_repeat_unanalyzed is not found: $!\n" if ($simple_repeat_unanalyzed !~ /\.gz$/);
    open (FILE, "gzip -dc $simple_repeat_unanalyzed |") or die "$simple_repeat_unanalyzed is not found: $!\n" if ($simple_repeat_unanalyzed =~ /\.gz$/);
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#|^$/);
        my @line = split (/\t/, $line);
        my $chr = $line[0];
        my $pos = $line[1];
        my $end = $line[2];
        my $id = $line[3];
        ${$uSTR{$chr}}{$pos} = $end;
        ${$uSTR2{$chr}}{$pos} = $id;
        $uSTR_id{$id} = "$chr=$pos=$end";
    }
    close (FILE);
}

if (-f $lowconf_TR){
    my %lowconfSTR_region;
    open (FILE, $lowconf_TR) or die "$lowconf_TR is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        if ($line =~ /(.+)-(.+)/){
            my $str1 = $1;
            my $str2 = $2;
            my $str_prefix = $1 if ($str1 =~ /^([^\d]+)\d+/);
            my $str1_num = $1 if ($str1 =~ /(\d+)/);
            my $str2_num = $1 if ($str2 =~ /(\d+)/);
            for (my $i = $str1_num; $i <= $str2_num; $i++){
                my $strid = $str_prefix . $i;
                $lowconf_STR{$strid} = 1;
                my ($chr, $pos, $end) = split (/=/, $uSTR_id{$strid});
                if (!exists ${$lowconfSTR_region{$chr}}{$pos}){
                    ${$lowconfSTR_region{$chr}}{$pos} = $end;
                }
                else{
                    my $pre_end = ${$lowconfSTR_region{$chr}}{$pos};
                    if ($end > $pre_end){
                        ${$lowconfSTR_region{$chr}}{$pos} = $end;
                    }
                }
            }
        }
        else{
            my @str = split (/\//, $line);
            foreach (@str){
                $lowconf_STR{$_} = 1;
                my ($chr, $pos, $end) = split (/=/, $uSTR_id{$_});
                if (!exists ${$lowconfSTR_region{$chr}}{$pos}){
                    ${$lowconfSTR_region{$chr}}{$pos} = $end;
                }
                else{
                    my $pre_end = ${$lowconfSTR_region{$chr}}{$pos};
                    if ($end > $pre_end){
                        ${$lowconfSTR_region{$chr}}{$pos} = $end;
                    }
                }
            }
        }
    }
    close (FILE);
    foreach my $chr (keys %lowconfSTR_region){
        my $pre_pos = 0;
        my $pre_end = 0;
        foreach my $pos (sort {$a <=> $b} keys %{$lowconfSTR_region{$chr}}){
            my $end = ${$lowconfSTR_region{$chr}}{$pos};
            if (($pre_pos > 0) and ($pos <= $pre_end)){
                if ($end <= $pre_end){
                    delete ${$lowconfSTR_region{$chr}}{$pos};
                    next;
                }
                else{
                    delete ${$lowconfSTR_region{$chr}}{$pos};
                    ${$lowconfSTR_region{$chr}}{$pre_pos} = $end;
                    $pre_end = $end;
                    next;
                }
            }
            $pre_pos = $pos;
            $pre_end = $end;
        }
    }
    foreach my $chr (keys %lowconfSTR_region){
        my $pre_end = 0;
        my $pre_len = 0;
        foreach my $pos (sort {$a <=> $b} keys %{$lowconfSTR_region{$chr}}){
            my $end = ${$lowconfSTR_region{$chr}}{$pos};
            my $len = $end - $pos + 1;
            my $distance = $pos - $pre_end + 1;
            if (($distance <= $len * 0.5) and ($distance <= $pre_len * 0.5)){
                ${$lowconf_region{$chr}}{$pre_end} = $pos;
            }
            $pre_end = $end;
            $pre_len = $len;
        }
    }
}

foreach my $type (keys %eval_call){
    next if ($type =~ /TR/);
    foreach my $chr (keys %{$eval_call{$type}}){
        my $chr02d = $chr;
        $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
        my %removed;
        foreach my $pos1 (sort {$a <=> $b} keys %{${$eval_call{$type}}{$chr}}){
            next if (exists $removed{$pos1});
            my $end1 = $pos1;
            my ($tag1, $len1, $info1) = split (/=/, ${${$eval_call{$type}}{$chr}}{$pos1});
            next if ($len1 < 20) and ($len1 > 0);
            $end1 = $pos1 + $len1 - 1 if ($type ne 'INS');
            my $read_num1 = 0;
            if ($info1 =~ /RN-(\d+)/){
                $read_num1 = $1;
            }
            elsif ($info1 =~ /RN1-(\d+)/){
                $read_num1 = $1;
                if ($info1 =~ /RN2-(\d+)/){
                    $read_num1 = $1 if ($read_num1 < $1);
                }
            }
            my @pos2;
            foreach my $pos2 (sort {$a <=> $b} keys %{${$eval_call{$type}}{$chr}}){
                next if ($pos2 <= $pos1);
                last if ($pos2 > $end1 + 1000);
                my $end2 = $pos2;
                my ($tag2, $len2, $info2) = split (/=/, ${${$eval_call{$type}}{$chr}}{$pos2});
                next if ($len2 < 20) and ($len2 > 0);
                $end2 = $pos2 + $len2 - 1 if ($type ne 'INS');
                my $read_num2 = 0;
                if ($info2 =~ /RN-(\d+)/){
                    $read_num2 = $1;
                }
                elsif ($info2 =~ /RN1-(\d+)/){
                    $read_num2 = $1;
                    if ($info2 =~ /RN2-(\d+)/){
                        $read_num2 = $1 if ($read_num2 < $1);
                    }
                }
                my $distance = $pos2 - $end1 + 1;
                my $flag = 0;
                if ($type eq 'INS'){
                    if (($distance <= $bp_diff) and (($len1 == 0) or ($len2 == 0))){
                        if (($len1 == 0) and ($len2 > 500)){
                            ${${$eval_call{$type}}{$chr}}{$pos2} .= "|${${$eval_call{$type}}{$chr}}{$pos1}";
                            delete ${${$eval_call{$type}}{$chr}}{$pos1};
                            delete ${$INS{$chr02d}}{$pos1};
                            $removed{$pos1} = 1;
                            $read_num2 += $read_num1;
                            $flag = 2;
                        }
                        elsif (($len2 == 0) and ($len1 > 500)){
                            ${${$eval_call{$type}}{$chr}}{$pos1} .= "|${${$eval_call{$type}}{$chr}}{$pos2}";
                            delete ${${$eval_call{$type}}{$chr}}{$pos2};
                            delete ${$INS{$chr02d}}{$pos2};
                            $removed{$pos2} = 1;
                            $read_num1 += $read_num2;
                            $flag = 1;
                        }
                        elsif (($len1 == 0) and ($len2 == 0)){
                            ${${$eval_call{$type}}{$chr}}{$pos1} .= "|${${$eval_call{$type}}{$chr}}{$pos2}";
                            delete ${${$eval_call{$type}}{$chr}}{$pos2};
                            delete ${$INS{$chr02d}}{$pos2};
                            $removed{$pos2} = 1;
                            $read_num1 += $read_num2;
                            $flag = 1;
                        }
                    }
                }
                else{
                    if ($distance < 0){
                        my $overlap = $end1 - $pos2 + 1;
                        $overlap = $len2 if ($end2 < $end1);
                        if (($overlap >= $len1 * $min_overlap_rate) and ($overlap >= $len2 * $min_overlap_rate)){
                            if (($tag2 !~ /BP/) and ($tag1 =~ /BP/)){
                                delete ${${$eval_call{$type}}{$chr}}{$pos1};
                                $removed{$pos1} = 1;
                                $flag = 2;
                            }
                            elsif (($tag2 =~ /BP/) and ($tag1 !~ /BP/)){
                                delete ${${$eval_call{$type}}{$chr}}{$pos2};
                                $removed{$pos2} = 1;
                                $flag = 1;
                            }
                            elsif (($type ne 'DEL') or (($tag2 =~ /BP/) and ($tag1 =~ /BP/))){
                                if ($read_num1 <= $read_num2){
                                    ${${$eval_call{$type}}{$chr}}{$pos2} .= "|${${$eval_call{$type}}{$chr}}{$pos1}";
                                    delete ${${$eval_call{$type}}{$chr}}{$pos1};
                                    $removed{$pos1} = 1;
                                    $read_num2 += $read_num1;
                                    $flag = 2;
                                }
                                elsif ($read_num1 > $read_num2){
                                    ${${$eval_call{$type}}{$chr}}{$pos1} .= "|${${$eval_call{$type}}{$chr}}{$pos2}";
                                    delete ${${$eval_call{$type}}{$chr}}{$pos2};
                                    $removed{$pos2} = 1;
                                    $read_num1 += $read_num2;
                                    $flag = 1;
                                }
                            }
                            elsif (($type eq 'DEL') and ((($len2 / $len1 <= 1.5) and ($len2 / $len1 >= 0.67)) or (($pos2 - $pos1 <= 200) and (abs ($end2 - $end1) <= 200)))){
                                if ($read_num1 <= $read_num2){
                                    ${${$eval_call{$type}}{$chr}}{$pos2} .= "|${${$eval_call{$type}}{$chr}}{$pos1}";
                                    delete ${${$eval_call{$type}}{$chr}}{$pos1};
                                    $removed{$pos1} = 1;
                                    $read_num2 += $read_num1;
                                    $flag = 2;
                                }
                                else{
                                    ${${$eval_call{$type}}{$chr}}{$pos1} .= "|${${$eval_call{$type}}{$chr}}{$pos2}";
                                    delete ${${$eval_call{$type}}{$chr}}{$pos2};
                                    $removed{$pos2} = 1;
                                    $read_num1 += $read_num2;
                                    $flag = 1;
                                }
                            }
                        }
                    }
                }
                if ($flag >= 1){
                    if (@pos2 > 0){
                        foreach my $pos3 (@pos2){
                            next if ($pos3 == $pos2);
                            next if (!exists ${${$eval_call{$type}}{$chr}}{$pos3});
                            my ($tag3, $len3) = split (/=/, ${${$eval_call{$type}}{$chr}}{$pos3});
                            if (($len3 < $len1) and ($len3 < $len2)){
                                delete ${${$eval_call{$type}}{$chr}}{$pos3};
                                $removed{$pos3} = 1;
                            }
                        }
                        @pos2 = ();
                    }
                    last if ($flag == 2);
                }
                push @pos2, $pos2;
            }
        }
    }
}

foreach my $type (keys %eval_call){
    foreach my $chr (keys %{$eval_call{$type}}){
        my $chr02d = $chr;
        $chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
        foreach my $pos (sort {$a <=> $b} keys %{${$eval_call{$type}}{$chr}}){
            my @info = split (/\|/, ${${$eval_call{$type}}{$chr}}{$pos});
            ${${$sv{$chr02d}}{$pos}}{$type} = ${${$eval_call{$type}}{$chr}}{$pos};
        }
    }
}

my %svtype;
my %sv2;
my %del_pos;

foreach my $chr (sort keys %sv){
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$sv{$chr}}){
        my $pos2 = $pos;
        $pos2 -- if (exists ${$STR{$chr}}{$pos});
        foreach my $type (keys %{${$sv{$chr}}{$pos}}){
            next if ($type =~ /TR/);
            my $sum_len = 0;
            my %tags;
            my $inv_align = 0;
            my @indel_len = ();
            my @dprate = ();
            my @mTD = ();
            my $insbplen_tag = '';
            my %gt;
            my @info_str;
            my @info = split (/\|/, ${${$sv{$chr}}{$pos}}{$type});
            my %rnum_len;
            my @IT_len;
            my @vrr;
            my @insbplen = ();
            my @duppos = ();
            my @duplen = ();
            my @sar = ();
            my $iduppos = '';
            my $iduplen = 0;
            my $dupbplen = 0;
            my $delbplen = 0;
            my $bplen2 = 0;
            my $insread = 0;
            my @insduplen = ();
            my $insduplen = 0;
            my $motif_info = '';
            my %MEI;
            my $type2 = $type;
            my $TAG = '';
            my $BP = 0;
            my $strid = '';
            foreach (@info){
                my ($tag, $len, $info) = split (/=/, $_);
                push @info_str, $info;
                my $read_num = 0;
                if ($info =~ /RN-(\d+)/){
                    $read_num = $1;
                }
                elsif ($info =~ /RN1-(\d+)/){
                    $read_num = $1;
                    if ($info =~ /RN2-(\d+)/){
                        $read_num = $1 if ($read_num < $1);
                    }
                }
                $tags{$tag} += $read_num;
                next if ($len eq 'NA');
                push @{$rnum_len{$read_num}}, $len;
                push @IT_len, $len if ($tag =~ /IT/);
                $BP += $1 if ($info =~ /BP-(\d+)/);
            }
            my $ave_len = 0;
            if (@IT_len > 0){
                @IT_len = sort {$b <=> $a} @IT_len;
                $ave_len = $IT_len[0];
            }
            elsif (scalar keys %rnum_len > 0){
                my @len;
                foreach my $rnum (sort {$b <=> $a} keys %rnum_len){
                    push @len, @{$rnum_len{$rnum}};
                    last;
                }
                @len = sort {$b <=> $a} @len;
                $ave_len = $len[0];
            }
            foreach my $tg (sort {$tags{$b} <=> $tags{$a}} keys %tags){
                $TAG = $tg if ($TAG eq '');
            }
            my $end = $pos;
            $end = $pos + $ave_len - 1 if ($type ne 'INS');
            my $readnum = 0;
            my $readnum2 = 0;
            foreach my $inf (@info_str){
                if ($inf =~ /RN-(\d+)/){
                    $readnum += $1;
                }
                elsif ($inf =~ /RN1-|RN2-/){
                    my $rn1 = 0;
                    my $rn2 = 0;
                    if ($inf =~ /RN1-(\d+)/){
                        $rn1 = $1;
                    }
                    if ($inf =~ /RN2-(\d+)/){
                        $rn2 = $1;
                    }
                    $readnum += int (($rn1 + $rn2) * 0.5 + 0.5);
                }
                if ($inf =~ /INSBPLEN-(\d+)/){
                    push @insbplen, $1;
                    $insbplen_tag = $1 if ($inf =~ /INSBPLEN-\d+\(([PF])/);
                }
                if ($inf =~ /ALIGN-(\d+)/){
                    $inv_align += $1;
                }
                if ($inf =~ /DPR-([\d\.]+)/){
                    push @dprate, $1;
                }
                if ($inf =~ /mTD-(\d+)/){
                    push @mTD, $1;
                }
                if ($inf =~ /INSREAD-(\d+)/){
                    $insread += $1;
                }
                if ($inf =~ /INSLEN-(\d+)/){
                    push @insduplen, $1;
                }
                if ($TAG !~ /intDUP/){
                    if ($inf =~ /DUPPOS-(\d+)/){
                        push @duppos, $1;
                    }
                    if ($inf =~ /DUPLEN-(\d+)/){
                        push @duplen, $1 if ($1 > 0);
                    }
                }
                if ($inf =~ /INSLEN-(\d+)/){
                    push @indel_len, $1 if ($1 > 0);
                }
                if ($inf =~ /VRR-([\d\.]+)/){
                    push @vrr, $1;
                }
                if ($inf =~ /GT-([^;]+)/){
                    $gt{$1} ++;
                }
                if ($inf =~ /SAR-([\d\.]+)/){
                    push @sar, $1;
                }
                if ($inf =~ /TR-([^;]+)/){
                    $strid = $1;
                }
                if ($inf =~ /MOTIF-([^;]+)/){
                    $motif_info = $1;
                }
                if (($TAG =~ /intDUP/) and ($inf =~ /DUPPOS-([^:]+:\d+)/)){
                    if ($iduppos eq ''){
                        $iduppos = $1;
                    }
                    if ($inf =~ /DUPLEN-(\d+)/){
                        if ($iduplen == 0){
                            $iduplen = $1;
                        }
                    }
                }
                if (($type eq 'INS') and ($inf =~ /DUPBPLEN-(\d+)/)){
                    $dupbplen = $1;
                }
                if (($type eq 'INS') and ($inf =~ /DELBPLEN-(\d+)/)){
                    $delbplen = $1;
                }
                if (($type eq 'INS') and ($inf =~ /BPLEN2-(\d+)/)){
                    $bplen2 = $1;
                }
                if ($inf =~ /MEI-(.+),MEILEN-(\d+),MEICN-(\d+)/){
                    $MEI{$1} = "MEILEN=$2;MEICN=$3";
                }
                if ($type eq 'INV'){
                    if ($inf =~ /RN12-(\d+)/){
                        $readnum2 += $1;
                    }
                }
            }
            my $dprate = 0;
            my $sum_dprate = 0;
            if (@dprate > 0){
                map{$sum_dprate += $_} @dprate;
                $dprate = int ($sum_dprate / @dprate * 100 + 0.5) / 100;
            }
            my $ave_mTD = 0;
            my $sum_mTD = 0;
            if (@mTD > 0){
                map{$sum_mTD += $_} @mTD;
                $ave_mTD = int ($sum_mTD / @mTD + 0.5);
            }
            my $dup_pos = 0;
            my $sum_duppos = 0;
            if (@duppos > 0){
                map{$sum_duppos += $_} @duppos;
                $dup_pos = int ($sum_duppos / @duppos + 0.5);
            }
            my $dup_len = 0;
            my $sum_duplen = 0;
            if (@duplen > 0){
                map{$sum_duplen += $_} @duplen;
                $dup_len = int ($sum_duplen / @duplen + 0.5);
            }
            my $indel_len = 0;
            my $sum_indel = 0;
            if (@indel_len > 0){
                map{$sum_indel += $_} @indel_len;
                $indel_len = int ($sum_indel / @indel_len + 0.5);
            }
            my $sum_insduplen = 0;
            if (@insduplen > 0){
                map{$sum_insduplen += $_} @insduplen;
                $insduplen = int ($sum_insduplen / @insduplen);
            }
            my $sum_vrr = 0;
            if (@vrr > 0){
                map{$sum_vrr += $_} @vrr;
            }
            my $sum_insbplen = 0;
            my $ave_insbplen = 0;
            if (@insbplen > 0){
                map{$sum_insbplen += $_} @insbplen;
                $ave_insbplen = int ($sum_insbplen / @insbplen + 0.5);
            }
            my $int_dup_info = '';
            my $MEI_info = '';
            if ($iduppos eq ''){
                $int_dup_info = "$iduppos-$iduplen";
            }
            if (scalar keys %MEI > 0){
                foreach my $mei (keys %MEI){
                    if ($MEI_info eq ''){
                        $MEI_info = "$mei;$MEI{$mei}";
                    }
                    else{
                        my $pre_meilen = $1 if ($MEI_info =~ /MEILEN=(\d+)/);
                        my $meilen = $1 if ($MEI{$mei} =~ /MEILEN=(\d+)/);
                        if ($meilen > $pre_meilen){
                            $MEI_info = "$mei;$MEI{$mei}";
                        }
                    }
                }
            }
            
            my $top_gt = '';
            my $sec_gt = '';
            my $top_num = 0;
            my $sec_num = 0;
            my $select_gt = 'NA';
            foreach my $gt (sort {$gt{$b} <=> $gt{$a}} keys %gt){
                if ($top_gt eq ''){
                    $top_gt = $gt;
                    $top_num = $gt{$gt};
                }
                elsif ($sec_gt eq ''){
                    $sec_gt = $gt;
                    $sec_num = $gt{$gt};
                    last;
                }
            }
            if ($top_num == $sec_num){
                if ($top_gt eq 'NA'){
                    $select_gt = $sec_gt;
                }
                if ($sec_gt eq 'NA'){
                    $select_gt = $top_gt;
                }
            }
            else{
                $select_gt = $top_gt;
            }
            my $sum_sar = 0;
            my $ave_sar = 0;
            map{$sum_sar += $_} @sar;
            $ave_sar = int ($sum_sar / @sar * 100 + 0.5) / 100 if (@sar > 0);
            my $line = '';
            if ($type eq 'INS'){
                my $ins_info = '';
                my $type_info = '';
                if ($dupbplen > 0){
                    $ins_info .= "BPLEN=$dupbplen;"; 
                    $type_info = "BP:DUP:";
                    $end += $dupbplen;
                }
                elsif ($delbplen > 0){
                    $ins_info .= "BPLEN=$delbplen;";
                    $type_info = "BP:DEL:";
                }
                elsif ($bplen2 > 0){
                    $ins_info .= "BPLEN=$bplen2;";
                }
                elsif ($ave_insbplen > 0){
                    $ins_info .= "BPLEN=$ave_insbplen;" if ($insbplen_tag eq '');
                    $ins_info .= "BPLEN=$ave_insbplen($insbplen_tag);" if ($insbplen_tag ne '');
                }
                elsif ($ave_len == 0){
                    $ins_info .= "BPLEN=0;";
                    $type_info = "BP:";
                }
                if ($iduplen > 0){
                    $ins_info .= "DUPPOS=$iduppos;DUPLEN=$iduplen;";
                    $type_info .= "$TAG:";
                }
                elsif ($dup_len > 0){
                    $ins_info .= "DUPPOS=$dup_pos;DUPLEN=$dup_len;" if ($dup_pos > 0);
                    $ins_info .= "DUPLEN=$dup_len;" if ($dup_pos == 0);
                    $type_info .= "$TAG:";
                }
                if ($MEI_info ne ''){
                    $ins_info .= "MEI=$MEI_info;";
                    my $mei_str = '';
                    foreach my $mei (keys %MEI){
                        $mei_str .= "$mei.";
                    }
                    $mei_str =~ s/\.$//;
                    $type_info .= "ME:$mei_str:";
                }
                if (($ins_info =~ /DUPLEN=/) and ($ins_info !~ /DUPPOS=/) and ($ins_info !~ /DUPBPLEN=/)){
                    $ins_info =~ s/DUPLEN=\d+;//;
                }
                if ($dprate > 0){
                    $ins_info .= "DPR=$dprate;";
                }
                if ($strid ne ''){
                    $ins_info .= "TRID=$strid;"
                }
                if ($motif_info ne ''){
                    $ins_info .= "TRUNIT=$motif_info;"
                }
                $ins_info =~ s/;$// if ($ins_info =~ /;$/);
                $type_info = "INS:$type_info";
                $type_info =~ s/:*$//;
                if ($type_info =~ /INS:IT/){
                    $type_info =~ s/INS:IT/INS/;
                    $type_info .= ':DUP' if ($ins_info =~ /DUP/) and ($type_info !~ /DUP/);
                }
                $type_info =~ s/DUP:intDUP/intDUP/ if ($type_info =~ /DUP:intDUP/);
                $type_info =~ s/^INS:/INS:BP:/ if ($ave_len == 0) and ($type_info !~ /:BP/);
                if ($ins_info eq ''){
                    $line = "$chr2\t$pos2\t.\t.\t<$type_info>\t.\tPASS\tSVTYPE=$type;SVLEN=$ave_len;READS=$readnum;VRR=$sum_vrr;SAR=$ave_sar;GT=$select_gt;END=$end" if ($BP == 0);
                    $line = "$chr2\t$pos2\t.\t.\t<$type_info>\t.\tPASS\tSVTYPE=$type;SVLEN=$ave_len;READS=$readnum;BP=$BP;VRR=$sum_vrr;SAR=$ave_sar;GT=$select_gt;END=$end" if ($BP > 0);
                }
                else{
                    $line = "$chr2\t$pos2\t.\t.\t<$type_info>\t.\tPASS\tSVTYPE=$type;SVLEN=$ave_len;READS=$readnum;$ins_info;VRR=$sum_vrr;SAR=$ave_sar;GT=$select_gt;END=$end" if ($BP == 0);
                    $line = "$chr2\t$pos2\t.\t.\t<$type_info>\t.\tPASS\tSVTYPE=$type;SVLEN=$ave_len;READS=$readnum;BP=$BP;$ins_info;VRR=$sum_vrr;SAR=$ave_sar;GT=$select_gt;END=$end" if ($BP > 0);
                }
            }
            elsif ($type eq 'INV'){
                $readnum += $inv_align;
                $readnum += $readnum2;
                $line = "$chr2\t$pos2\t.\t.\t<$type>\t.\tPASS\tSVTYPE=$type;SVLEN=$ave_len;READS=$readnum;INVALIN=$inv_align;VRR=$sum_vrr;SAR=$ave_sar;GT=$select_gt;END=$end";
            }
            elsif ($type eq 'DUP'){
                my $dup_info = '';
                my $type_info = '';
                my $CN = 2;
                if (($select_gt eq 'HM') and ($sum_vrr <= 0.65)){
                   $select_gt = 'HT';
                }
                if ($MEI_info ne ''){
                    $dup_info .= "MEI=$MEI_info;";
                    my $mei_str = '';
                    foreach my $mei (keys %MEI){
                        $mei_str .= "$mei.";
                    }
                    $mei_str =~ s/\.$//;
                    $type_info .= "MEI:$mei_str:";
                }
                if (@mTD > 0){
                    if ($select_gt eq 'HT'){
                        $CN = $ave_mTD + 1;
                    }
                    else{
                        $CN += int ($dprate + 0.5);
                    }
                }
                else{
                    if ($select_gt eq 'HT'){
                        my $add_cn = 1;
                        $add_cn = 2 if ($dprate >= 1.9);
                    }
                    else{
                        my $add_cn = int ($dprate + 0.5);
                        $add_cn = 2 if ($add_cn < 2);
                        $CN += $add_cn;
                    }
                }
                $dup_info .= "CN=$CN;";
                if (@dprate > 0){
                    $dup_info .= "DPR=$dprate;";
                }
                if ($insread > 0){
                    $dup_info .= "INSREAD=$insread;"
                }
                if ($insduplen > 0){
                    $dup_info .= "INSLEN=$insduplen;"
                }
                $dup_info =~ s/;$//;
                $type_info =~ s/:$//;
                $type2 .= ":$type_info" if ($type_info ne '');
                $line = "$chr2\t$pos2\t.\t.\t<$type2>\t.\tPASS\tSVTYPE=$type;SVLEN=$ave_len;READS=$readnum;$dup_info;VRR=$sum_vrr;SAR=$ave_sar;GT=$select_gt;END=$end";
            }
            else{
                my $pos2 = $pos;
                $pos2 -- if (exists ${$STR{$chr}}{$pos});
                if ($BP == 0){
                    if (@dprate > 0){
                        $line = "$chr2\t$pos2\t.\t.\t<$type>\t.\tPASS\tSVTYPE=$type;SVLEN=$ave_len;READS=$readnum;DPR=$dprate;VRR=$sum_vrr;SAR=$ave_sar;GT=$select_gt;END=$end";
                    }
                    else{
                        $line = "$chr2\t$pos2\t.\t.\t<$type>\t.\tPASS\tSVTYPE=$type;SVLEN=$ave_len;READS=$readnum;VRR=$sum_vrr;SAR=$ave_sar;GT=$select_gt;END=$end";
                    }
                }
                else{
                    if (@dprate > 0){
                        $line = "$chr2\t$pos2\t.\t.\t<$type>\t.\tPASS\tSVTYPE=$type;SVLEN=$ave_len;READS=$readnum;BP=$BP;DPR=$dprate;VRR=$sum_vrr;SAR=$ave_sar;GT=$select_gt;END=$end" if ($indel_len == 0);
                        $line = "$chr2\t$pos2\t.\t.\t<$type>\t.\tPASS\tSVTYPE=$type;SVLEN=$ave_len;READS=$readnum;BP=$BP;DPR=$dprate;VRR=$sum_vrr;SAR=$ave_sar;GT=$select_gt;INSLEN=$indel_len;END=$end" if ($indel_len > 0);
                    }
                    else{
                        $line = "$chr2\t$pos2\t.\t.\t<$type>\t.\tPASS\tSVTYPE=$type;SVLEN=$ave_len;READS=$readnum;BP=$BP;VRR=$sum_vrr;SAR=$ave_sar;GT=$select_gt;END=$end" if ($indel_len == 0);
                        $line = "$chr2\t$pos2\t.\t.\t<$type>\t.\tPASS\tSVTYPE=$type;SVLEN=$ave_len;READS=$readnum;BP=$BP;VRR=$sum_vrr;SAR=$ave_sar;GT=$select_gt;INSLEN=$indel_len;END=$end" if ($indel_len > 0);
                    }
                }
            }
            ${${$sv2{$chr}}{$pos2}}{$type} = $line;
        }
    }
}

foreach my $chr (sort keys %sv){
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    my $pre_pos = 0;
    my $pre_end = 0;
    my $pre_read = 0;
    my $pre_vrr = 0;
    my $pre_len = 0;
    my $pre_type = '';
    foreach my $pos (sort {$a <=> $b} keys %{$sv{$chr}}){
        my $stype1 = '';
        my $stype2 = '';
        my $len1 = 0;
        my $len2 = 0;
        my $info1 = '';
        my $info2 = '';
        foreach my $type (keys %{${$sv{$chr}}{$pos}}){
            next if ($type !~ /TR/);
            my ($tag, $len, $info) = split (/=/, ${${$sv{$chr}}{$pos}}{$type});
            my $type2 = 'INS';
            $type2 = 'DEL' if ($type =~ /del/);
            if ($stype1 eq ''){
                $stype1 = $type2;
                $len1 = $len;
                $info1 = $info;
            }
            else{
                $stype2 = $type2;
                $len2 = $len;
                $info2 = $info;
            }
        }
        next if ($stype1 eq '');
        my $line = '';
        my $end = $pos + $len1 - 1;
        $end = $pos if ($stype1 eq 'INS');
        my $readnum = $1 if ($info1 =~ /RN-(\d+)/);
        my $readnum1 = $readnum;
        my $readnum2 = 0;
        $readnum2 = $1 if ($info2 =~ /RN-(\d+)/);
        my $read_cov = $1 if ($info1 =~ /RD-([\d\.]+)/);
        my $vrr = int ($readnum / $read_cov * 100 + 0.5) / 100;
        my $BP = 0;
        $BP = $1 if ($info2 =~ /BP-(\d+)/);
        $BP = $1 if ($info1 =~ /BP-(\d+)/);
        if (($stype1 eq 'DEL') and ($stype2 eq 'DEL')){
            if (($readnum1 > $readnum2) and ($readnum1 / $readnum2 >= 4)){
                $readnum += $readnum2;
                $stype2 = '';
            }
            elsif (($readnum2 > $readnum1) and ($readnum2 / $readnum1 >= 4)){
                $readnum2 += $readnum1;
                $vrr = 0;
                $stype1 = '';
            }
        }
        my $gt = $1 if ($info1 =~ /GT-(.+?);/);
        my $strid = $1 if ($info1 =~ /TR-([^;]+)/);
        my $sar = $1 if ($info1 =~ /SAR-([\d\.]+)/);
        if ($info2 =~ /SAR-([\d\.]+)/){
            my $sar2 = $1;
            $sar = int (($sar + $sar2) * 0.5 * 100 + 0.5) / 100;
        }
        my ($spos, $send, $unit_size) = split (/=/, $STR2{$strid});
        my $cn = int ($len1 / $unit_size * 10 + 0.5) / 10;
        if ($info1 =~ /TRCN-([\d\.]+)/){
            $cn = $1;
        }
        my $cn_str = "gain+$cn" if ($stype1 eq 'INS');
        $cn_str = "loss-$cn" if ($stype1 eq 'DEL');
        my $strdup = 0;
        $strdup = $1 if ($info1 =~ /TRDUP-([\d\.]+)/);
        $strdup = $1 if ($info2 =~ /TRDUP-([\d\.]+)/) and ($strdup == 0);
        my $bp_dup_info = '';
        if ($BP > 0){
            $bp_dup_info = "BP=$BP";
            $bp_dup_info .= ";TRDUP=$strdup" if ($strdup > 0);
        }
        elsif ($strdup > 0){
            $bp_dup_info = "TRDUP=$strdup"
        }
        if ($stype2 ne ''){
            my $end2 = $pos + $len2 - 1;
            $end2 = $pos if ($stype2 eq 'INS');
            my $vrr2 = int ($readnum2 / $read_cov * 100 + 0.5) / 100;
            my $strmatch2 = 0;
            $strmatch2 = $1 if ($info2 =~ /TRMATCH-(\d+)/);
            if ($vrr2 >= 0.1){
                my $cn2 = int ($len2 / $unit_size * 10 + 0.5) / 10;
                $cn2 = int ($strmatch2 / $unit_size * 10 + 0.5) / 10 if ($strmatch2 > 0);
                my $cn_str2 = "gain+$cn2" if ($stype2 eq 'INS');
                $cn_str2 = "loss-$cn2" if ($stype2 eq 'DEL');
                my $gt2 = $1 if ($info2 =~ /GT-(.+?);/);
                if ($stype1 eq ''){
                    $len1 = $len2;
                    $readnum = $readnum2;
                    $end = $end2;
                    $vrr = $vrr2;
                    $gt = $gt2;
                    $cn_str = $cn_str2;
                    $stype1 = $stype2;
                }
                else{
                    if ($readnum >= $readnum2){
                        $len1 .= ",$len2";
                        $readnum .= ",$readnum2";
                        $end = $end2 if ($end2 > $end);
                        $vrr .= ",$vrr2";
                        $gt = 'HT2';
                        $cn_str .= ",$cn_str2";
                        $stype1 .= ",$stype2" if ($stype1 ne $stype2);
                    }
                    else{
                        $len1 = "$len2,$len1";
                        $readnum = "$readnum2,$readnum";
                        $end = $end2 if ($end2 > $end);
                        $vrr = "$vrr2,$vrr";
                        $gt = 'HT2';
                        $cn_str = "$cn_str2,$cn_str";
                        $stype1 = "$stype2,$stype1" if ($stype1 ne $stype2);
                    }
                }
            }
        }
        next if ($stype1 eq '');
        $line = "$chr2\t$pos\t.\t.\t<TR:CNV>\t.\tPASS\tSVTYPE=$stype1;SVLEN=$len1;READS=$readnum;CN=$cn_str;VRR=$vrr;SAR=$sar;GT=$gt;END=$end;TRID=$strid;TREND=$send;TRULEN=$unit_size";
        $line = "$chr2\t$pos\t.\t.\t<TR:CNV>\t.\tPASS\tSVTYPE=$stype1;SVLEN=$len1;READS=$readnum;CN=$cn_str;VRR=$vrr;SAR=$sar;GT=$gt;$bp_dup_info;END=$end;TRID=$strid;TREND=$send;TRULEN=$unit_size" if ($bp_dup_info ne '');
        ${${$sv2{$chr}}{$pos}}{'TR'} = $line;
        $pre_pos = $pos;
        $pre_end = $send;
        $pre_read = $readnum;
        $pre_vrr = $vrr;
        $pre_len = $len1;
        $pre_type = $stype1;
    }
}

foreach my $chr (sort keys %sv2){
    my $pre_pos = 0;
    my $pre_end = 0;
    my $pre_read = 0;
    my $pre_svlen = 0;
    my $pre_vrr = 0;
    my $pre_type = '';
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$sv2{$chr}}){
        foreach my $type (keys %{${$sv2{$chr}}{$pos}}){
            next if ($type =~ /TR/);
            my @line = split (/\t/, ${${$sv2{$chr}}{$pos}}{$type});
            my $end = $1 if ($line[7] =~ /END=(\d+)/);
            my $read = $1 if ($line[7] =~ /READS=(\d+)/);
            my $svlen = $1 if ($line[7] =~ /SVLEN=(\d+)/);
            next if ($svlen < 20);
            my $vrr = $1 if ($line[7] =~ /VRR=([\d\.]+)/);
            ${$del_pos{$chr2}}{$pos} = $svlen if ($type eq 'DEL');
            if ((abs ($pre_end - $pos) < $bp_diff) or (abs ($pre_pos - $pos) < $bp_diff)){
                if (($type eq 'INV') and ($pre_type eq 'INS') and ($pre_svlen <= 150)){
                    delete ${${$sv2{$chr}}{$pre_pos}}{$pre_type};
                }
                elsif (($pre_type eq 'INV') and ($type eq 'INS') and ($svlen <= 150)){
                    delete ${${$sv2{$chr}}{$pos}}{$type};
                    next;
                }
            }
            if ((abs ($pre_end - $pos) < $ins_bp_diff) or (abs ($pre_pos - $pos) < $ins_bp_diff)){
                if (($type =~ /INS|DUP/) and ($pre_type =~ /DUP|INS/) and ($type ne $pre_type) and (abs ($svlen - $pre_svlen) <= $svlen * 0.5) and (abs ($svlen - $pre_svlen) <= $pre_svlen * 0.5)){
                    my $new_read = $read + $pre_read;
                    my $new_vrr = $vrr + $pre_vrr;
                    $new_vrr = 1 if ($new_vrr > 1);
                    if ($read >= $pre_read){
                        $line[7] =~ s/READS=\d+/READS=$new_read/;
                        $line[7] =~ s/VRR=[\d\.]+/VRR=$new_vrr/;
                        my $new_line = join ("\t", @line);
                        ${${$sv2{$chr}}{$pos}}{$type} = $new_line;
                        delete ${${$sv2{$chr}}{$pre_pos}}{$pre_type};
                    }
                    else{
                        my $pre_line = ${${$sv2{$chr}}{$pre_pos}}{$pre_type};
                        my @pre_line = split (/\t/, $pre_line);
                        $pre_line[7] =~ s/READS=\d+/READS=$new_read/;
                        $pre_line[7] =~ s/VRR=[\d\.]+/VRR=$new_vrr/;
                        my $new_line = join ("\t", @pre_line);
                        ${${$sv2{$chr}}{$pre_pos}}{$pre_type} = $new_line;
                        $pre_read = $new_read;
                        $pre_vrr = $new_vrr;
                        delete ${${$sv2{$chr}}{$pos}}{$type};
                        next;
                    }
                }
            }
            if (($pre_type eq 'DEL') and ($pre_end > $end) and (exists ${${$sv2{$chr}}{$pre_pos}}{$pre_type})){
                if (($pre_vrr >= 0.9) and ($vrr <= 0.2)){
                    delete ${${$sv2{$chr}}{$pos}}{$type};
                    next;
                }
                elsif (($type eq 'DEL') and ($pre_vrr <= 0.2) and ($vrr >= 0.9)){
                    delete ${${$sv2{$chr}}{$pre_pos}}{$pre_type};
                }
                else{
                    next;
                }
            }
            $pre_pos = $pos;
            $pre_type = $type;
            $pre_end = $end;
            $pre_read = $read;
            $pre_svlen = $svlen;
            $pre_vrr = $vrr;
        }
    }
}

foreach my $chr (sort keys %sv2){       # merge homozygous neigboring/overlapping DELs
    my $pre_pos = 0;
    my $pre_end = 0;
    my $pre_len = 0;
    my $pre_read = 0;
    my $pre_vrr = 0;
    my $pre_dpr = 0;
    my $pre_sar = 0;
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$sv2{$chr}}){
        foreach my $type (keys %{${$sv2{$chr}}{$pos}}){
            next if ($type !~ /DEL/);
            my @line = split (/\t/, ${${$sv2{$chr}}{$pos}}{$type});
            my $gt = $1 if ($line[7] =~ /GT=(.+?);/);
            next if ($gt ne 'HM');
            my $end = $1 if ($line[7] =~ /END=(\d+)/);
            my $len = $1 if ($line[7] =~ /SVLEN=(\d+)/);
            next if ($len < 3);
            my $vrr = $1 if ($line[7] =~ /VRR=([\d\.]+)/);
            my $read = $1 if ($line[7] =~ /READS=(\d+)/);
            my $dpr = 1;
            $dpr = $1 if ($line[7] =~ /DPR=([\d\.]+)/);
            my $sar = $1 if ($line[7] =~ /SAR=([\d\.]+)/);
            my $distance = $pos - $pre_end + 1;
            my $flag = 0;
            if (($pre_end > 0) and ($end > $pre_end) and (abs ($pre_end - $pos) < $bp_diff) and ($distance <= $len * 0.5) and ($distance <= $pre_len * 0.5) and ($pre_vrr >= 0.2) and ($vrr >= 0.2) and (($pre_vrr >= 0.8) or ($vrr >= 0.8))){
                $flag = 1;
            }
            elsif (($distance <= 500) and ($end > $pre_end) and ($distance <= $len * 0.2) and ($distance <= $pre_len * 0.2) and ($pre_vrr >= 0.2) and ($vrr >= 0.2) and (($pre_vrr >= 0.8) or ($vrr >= 0.8))){
                $flag = 1;
            }
            elsif (($pre_end > 0) and ($distance < 0) and ($pre_vrr >= 0.8) and ($vrr < 0.2)){
                $flag = 2;
            }
            elsif (($pre_end > 0) and ($distance < 0) and ($pre_vrr < 0.2) and ($vrr >= 0.8)){
                $flag = 3;
            }
            if ($flag == 1){
                my $new_len = $end - $pre_pos + 1;
                my $pre_line = ${${$sv2{$chr}}{$pre_pos}}{$type};
                $read = $pre_read if ($pre_read > $read);
                $vrr = $pre_vrr if ($pre_vrr > $vrr);
                $dpr = $pre_dpr if ($pre_len > $len);
                $sar = int (($pre_sar + $sar) / 2 * 100 + 0.5) / 100;
                $pre_line =~ s/SVLEN=\d+/SVLEN=$new_len/;
                $pre_line =~ s/READS=\d+/READS=$read/;
                $pre_line =~ s/VRR=[\d\.]+/VRR=$vrr/;
                $pre_line =~ s/DPR=[\d\.]+/DPR=$dpr/ if ($pre_line =~ /DPR=/) and ($line[7] =~ /DPR=/);
                $pre_line =~ s/END=\d+/END=$end/;
                $pre_line =~ s/SAR=[\d\.]+/SAR=$sar/;
                ${${$sv2{$chr}}{$pre_pos}}{$type} = $pre_line;
                delete ${${$sv2{$chr}}{$pos}}{$type};
                $pre_end = $end;
                $pre_len = $new_len;
                $pre_read = $read;
                $pre_vrr = $vrr;
                $pre_dpr = $dpr;
                $pre_sar = $sar;
                next;
            }
            elsif ($flag == 2){
                if (($len / $pre_len <= 1.5) and ($len / $pre_len >= 0.67)){
                    my $pre_line = ${${$sv2{$chr}}{$pre_pos}}{$type};
                    $pre_read += $read;
                    $pre_vrr += $vrr;
                    $pre_vrr = 1 if ($pre_vrr > 1);
                    $pre_line =~ s/READS=\d+/READS=$pre_read/;
                    $pre_line =~ s/VRR=[\d\.]+/VRR=$pre_vrr/;
                    ${${$sv2{$chr}}{$pre_pos}}{$type} = $pre_line;
                }
                delete ${${$sv2{$chr}}{$pos}}{$type};
                next;
            }
            elsif ($flag == 3){
                if (($len / $pre_len <= 1.5) and ($len / $pre_len >= 0.67)){
                    my $line = ${${$sv2{$chr}}{$pos}}{$type};
                    $read += $pre_read;
                    $vrr += $pre_vrr;
                    $vrr = 1 if ($vrr > 1);
                    $line =~ s/READS=\d+/READS=$read/;
                    $line =~ s/VRR=[\d\.]+/VRR=$vrr/;
                    ${${$sv2{$chr}}{$pos}}{$type} = $line;
                }
                delete ${${$sv2{$chr}}{$pre_pos}}{$type};
            }
            $pre_pos = $pos;
            $pre_end = $end;
            $pre_len = $len;
            $pre_read = $read;
            $pre_vrr = $vrr;
            $pre_dpr = $dpr;
            $pre_sar = $sar;
        }
    }
}

my %str_sv;
my %ins_sv;
my %dup_sv;

foreach my $chr (keys %sv2){
    foreach my $pos (keys %{$sv2{$chr}}){
        foreach my $type (keys %{${$sv2{$chr}}{$pos}}){
            if ($type eq 'TR'){
                my $str_end = $1 if (${${$sv2{$chr}}{$pos}}{$type} =~ /TREND=(\d+)/);
                ${$str_sv{$chr}}{$pos} = $str_end;
            }
            else{
                my $end = $1 if (${${$sv2{$chr}}{$pos}}{$type} =~ /;END=(\d+)/);
                if ($type eq 'DUP'){
                    ${$dup_sv{$chr}}{$pos} = $end;
                }
                elsif ($type eq 'INS'){
                    my $len = $1 if (${${$sv2{$chr}}{$pos}}{$type} =~ /;SVLEN=(\d+)/);
                    ${$ins_sv{$chr}}{$pos} = $len if ($len >= 20);
                }
            }
        }
    }
}

my $pre_ipos = 0;
my $pre_ilen = 0;
foreach my $chr (keys %ins_sv){ 
    foreach my $ipos (sort {$a <=> $b} keys %{$ins_sv{$chr}}){
        next if (!exists ${${$sv2{$chr}}{$ipos}}{'INS'});
        my $ilen = ${$ins_sv{$chr}}{$ipos};
        my $distance = $ipos - $pre_ipos;
        if (($distance <= $bp_diff3) and ($distance < $pre_ilen) and ($distance < $ilen) and ($pre_ilen > 0) and ($ilen > 0) and ($pre_ilen / $ilen <= 1.5) and ($pre_ilen / $ilen >= 0.66) and (exists ${${$sv2{$chr}}{$pre_ipos}}{'INS'})){
            my $pre_line = ${${$sv2{$chr}}{$pre_ipos}}{'INS'};
            my $line = ${${$sv2{$chr}}{$ipos}}{'INS'};
            my $pre_read = $1 if ($pre_line =~ /READS=(\d+)/);
            my $read = $1 if ($line =~ /READS=(\d+)/);
            my $pre_vrr = $1 if ($pre_line =~ /VRR=([\d\.]+)/);
            my $vrr = $1 if ($line =~ /VRR=([\d\.]+)/);
            my $new_read = $pre_read + $read;
            my $new_vrr = $pre_vrr + $vrr;
            if (($pre_read > 0) and (($read / $pre_read < 0.8) or ($read / $pre_read > 1.25))){
                if ($read >= $pre_read){
                    $line =~ s/READS=\d+/READS=$new_read/;
                    $line =~ s/VRR=[\d\.]+/VRR=$new_vrr/;
                    delete ${${$sv2{$chr}}{$pre_ipos}}{'INS'};
                    ${${$sv2{$chr}}{$ipos}}{'INS'} = $line;
                }
                else{
                    $pre_line =~ s/READS=\d+/READS=$new_read/;
                    $pre_line =~ s/VRR=[\d\.]+/VRR=$new_vrr/;
                    delete ${${$sv2{$chr}}{$ipos}}{'INS'};
                    ${${$sv2{$chr}}{$pre_ipos}}{'INS'} = $pre_line;
                    next;
                }
            }
        }
        $pre_ipos = $ipos;
        $pre_ilen = $ilen;
    }
}

foreach my $chr (sort keys %str_sv){            # merge and delete non-TR-INSs within TR-SVs +/- 100 bp
    foreach my $pos (sort {$a <=> $b} keys %{$str_sv{$chr}}){
        my $end = ${$str_sv{$chr}}{$pos};
        my $sline = ${${$sv2{$chr}}{$pos}}{'TR'};
        my $stype = $1 if ($sline =~ /SVTYPE=(.+?);/);
        my $stype1 = '';
        my $stype2 = '';
        my $read1 = 0;
        my $read2 = 0;
        my $vrr1 = 0;
        my $vrr2 = 0;
        my $len1 = 0;
        my $len2 = 0;
        my $cn1 = 0;
        $stype1 = $stype if ($stype !~ /,/);
        if ($stype =~ /(.+),(.+)/){
            $stype1 = $1;
            $stype2 = $2;
        }
        $read1 = $1 if ($sline =~ /READS=(\d+)/);
        $read2 = $1 if ($sline =~ /READS=\d+,(\d+)/);
        $vrr1 = $1 if ($sline =~ /VRR=([\d\.]+)/);
        $vrr2 = $1 if ($sline =~ /VRR=[\d\.]+,([\d\.]+)/);
        $len1 = $1 if ($sline =~ /SVLEN=(\d+)/);
        $len2 = $1 if ($sline =~ /SVLEN=\d+,(\d+)/);
        $cn1 = $1 if ($sline =~ /CN=gain\+([\d\.]+)/) and ($stype1 eq 'INS');
        $cn1 = $1 if ($sline =~ /CN=lossn-([\d\.]+)/) and ($stype1 eq 'DEL');
        my $strulen = 0;
        $strulen = $1 if ($sline =~ /TRULEN=(\d+)/);
        my $maxlen = $len1;
        $maxlen = $len2 if ($len2 > $len1);
        my $strlen = $end - $pos + 1;
        my $diff = 20;
        foreach my $ipos (sort {$a <=> $b} keys %{$ins_sv{$chr}}){
            last if ($ipos > $end + $diff);
            my $type = 'INS';
            my $ilen = ${$ins_sv{$chr}}{$ipos};
            next if ($ipos < $pos - $diff);
            next if (!exists ${${$sv2{$chr}}{$ipos}}{$type});
            my $iline = ${${$sv2{$chr}}{$ipos}}{$type};
            next if ($iline =~ /TRID=/);
            my $iread = $1 if ($iline =~ /READS=(\d+)/);
            my $ivrr = $1 if ($iline =~ /VRR=([\d\.]+)/);
            my $iend2 = $ipos;
            my $bp_distance = 0;
            $bp_distance = $1 if ($iline =~ /BPLEN=(\d+)/);
            $iend2 += $bp_distance;
            if (($stype1 eq 'INS') and ($stype2 eq '') and ($len2 > 0)){
                if (($bp_distance > 0) and (abs ($ipos - $pos) <= $diff) and (abs ($iend2 - $ipos) <= $diff)){
                    if ($maxlen < $bp_distance){
                        my $cn = int ($strlen / $strulen * 10 + 0.5) / 10;
                        if ($len1 >= $len2){
                            $read1 += $iread;
                            $sline =~ s/SVLEN=$len1/SVLEN=$strlen/;
                            $sline =~ s/VRR=$vrr1/VRR=1/;
                            $sline =~ s/READS=\d+/READS=$read1/;
                            $sline =~ s/CN=gain\+[\d\.]+/CN=gain+$cn/;
                        }
                        else{
                            $read2 += $iread;
                            $sline =~ s/SVLEN=$len1,$len2/SVLEN=$strlen,$len1/;
                            $sline =~ s/VRR=$vrr1,$vrr2/VRR=1,$vrr1/;
                            $sline =~ s/READS=$read1,\d+/READS=$read2,$read1/;
                            $sline =~ s/CN=gain\+[\d\.]+,gain\+[\d\.]+/CN=gain+$cn,gain+$cn1/;
                        }
                    }
                }
                elsif (($ilen > 0) and ($len1 / $ilen <= 1.5) and ($len1 / $ilen >= 0.67)){
                    $read1 += $iread;
                    $vrr1 += $ivrr;
                    $sline =~ s/READS=\d+/READS=$read1/;
                    $sline =~ s/VRR=[\d\.]+/VRR=$vrr1/;
                }
                elsif (($ilen > 0) and ($len2 / $ilen <= 1.5) and ($len2 / $ilen >= 0.67)){
                    $read2 += $iread;
                    $vrr2 += $ivrr;
                    if ($vrr1 >= $vrr2){
                        $sline =~ s/READS=\d+,\d+/READS=$read1,$read2/;
                        $sline =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$vrr1,$vrr2/;
                    }
                    else{
                        $sline =~ s/READS=\d+,\d+/READS=$read2,$read1/;
                        $sline =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$vrr2,$vrr1/;
                    }
                }
                else{
                    next;
                }
            }
            elsif (($stype1 eq 'INS') and ($len2 == 0)){
                if (($bp_distance > 0) and (abs ($ipos - $pos) <= $diff) and (abs ($iend2 - $ipos) <= $diff)){
                    if ($len1 < $bp_distance){
                        $read1 += $iread;
                        my $cn = int ($strlen / $strulen * 10 + 0.5) / 10;
                        $sline =~ s/SVLEN=$len1/SVLEN=$strlen/;
                        $sline =~ s/VRR=$vrr1/VRR=1/;
                        $sline =~ s/READS=\d+/READS=$read1/;
                        $sline =~ s/CN=gain\+[\d\.]+/CN=gain+$cn/;
                        if (($sline =~ /GT=HT;/) and ($vrr1 >= 0.75)){
                            $sline =~ s/GT=HT/GT=HM/;
                        }
                    }
                }
                elsif (($ilen > 0) and ($len1 / $ilen <= 1.5) and ($len1 / $ilen >= 0.67)){
                    $read1 += $iread;
                    $vrr1 += $ivrr;
                    $sline =~ s/READS=\d+/READS=$read1/;
                    $sline =~ s/VRR=[\d\.]+/VRR=$vrr1/;
                    if (($sline =~ /GT=HT;/) and ($vrr1 >= 0.75)){
                        $sline =~ s/GT=HT/GT=HM/;
                    }
                }
                elsif ($ilen > 0){
                    my $cn = int ($ilen / $strulen * 10 + 0.5) / 10;
                    if ($vrr1 >= $ivrr){
                        $sline =~ s/SVLEN=\d+/SVLEN=$len1,$ilen/;
                        $sline =~ s/READS=\d+/READS=$read1,$iread/;
                        $sline =~ s/VRR=[\d\.]+/VRR=$vrr1,$ivrr/;
                        $sline =~ s/gain\+[\d\.]+/gain+$cn1,gain+$cn/;
                        $sline =~ s/GT=H[TM]/GT=HT2/;
                    }
                    else{
                        $sline =~ s/SVLEN=\d+/SVLEN=$ilen,$len1/;
                        $sline =~ s/READS=\d+/READS=$iread,$read1/;
                        $sline =~ s/VRR=[\d\.]+/VRR=$ivrr,$vrr1/;
                        $sline =~ s/gain\+[\d\.]+/gain+$cn,gain+$cn1/;
                        $sline =~ s/GT=H[TM]/GT=HT2/;
                    }
                }
                else{
                    next;
                }
            }
            elsif (($stype1 eq 'INS') and ($stype2 eq 'DEL')){
                if (($bp_distance > 0) and (abs ($ipos - $pos) <= $diff) and (abs ($iend2 - $ipos) <= $diff)){
                    if ($len1 < $bp_distance){
                        $read1 += $iread;
                        my $cn = int ($strlen / $strulen * 10 + 0.5) / 10;
                        $sline =~ s/SVLEN=\d+/SVLEN=$strlen/;
                        $sline =~ s/VRR=[\d\.]+/VRR=1/;
                        $sline =~ s/READS=\d+/READS=$read1/;
                        $sline =~ s/CN=gain\+[\d\.]+/CN=gain+$cn/;
                    }
                }
                elsif (($ilen > 0) and ($len1 / $ilen <= 1.5) and ($len1 / $ilen >= 0.65)){
                    $read1 += $iread;
                    $vrr1 += $ivrr;
                    $sline =~ s/READS=\d+/READS=$read1/;
                    $sline =~ s/VRR=[\d\.]+/VRR=$vrr1/;
                }
                else{
                    next;
                }
            }
            elsif (($stype1 eq 'DEL') and ($stype2 eq 'INS')){
                if (($bp_distance > 0) and (abs ($ipos - $pos) <= $diff) and (abs ($iend2 - $ipos) <= $diff)){
                    if ($len2 < $bp_distance){
                        $read2 += $iread;
                        my $cn = int ($strlen / $strulen * 10 + 0.5) / 10;
                        $sline =~ s/SVTYPE=$stype1,$stype2/SVTYPE=$stype2,$stype1/;
                        $sline =~ s/SVLEN=$len1,$len2/SVLEN=$strlen,$len1/;
                        $sline =~ s/VRR=$vrr1,$vrr2/VRR=1,$vrr1/;
                        $sline =~ s/READS=$read1,\d+/READS=$read2,$read1/;
                        my $cn1 = $1 if ($sline =~ /CN=loss-([\d\.]+)/);
                        $sline =~ s/CN=loss-[\d\.]+,gain\+[\d\.]+/CN=gain+$cn,loss-$cn1/;
                    }
                }
                elsif (($ilen > 0) and ($len2 / $ilen <= 1.5) and ($len2 / $ilen >= 0.65)){
                    $read2 += $iread;
                    $vrr2 += $ivrr;
                    $sline =~ s/READS=\d+,\d+/READS=$read1,$read2/;
                    $sline =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$vrr1,$vrr2/;
                }
                else{
                    next;
                }
            }
            elsif (($stype1 eq 'DEL') and ($len2 == 0)){
                my $cn = int ($ilen / $strulen * 10 + 0.5) / 10;
                if (($bp_distance > 0) and (abs ($ipos - $pos) <= $diff) and (abs ($iend2 - $ipos) <= $diff)){
                    if ($len2 < $bp_distance){
                        $sline =~ s/SVTYPE=$stype1/SVTYPE=INS,$stype1/;
                        $sline =~ s/SVLEN=$len1/SVLEN=$ilen,$len1/;
                        $sline =~ s/VRR=$vrr1/VRR=1,$vrr1/;
                        $sline =~ s/READS=$read1/READS=$iread,$read1/;
                        $sline =~ s/CN=loss-[\d\.]+/CN=gain+$cn,loss-$cn1/;
                        $sline =~ s/GT=H[TM]/GT=HT2/;
                    }
                }
                elsif ($ilen > 0){
                    if ($vrr1 >= $ivrr){
                        $sline =~ s/SVTYPE=$stype1/SVTYPE=$stype1,INS/;
                        $sline =~ s/SVLEN=$len1/SVLEN=$len1,$ilen/;
                        $sline =~ s/READS=\d+/READS=$read1,$iread/;
                        $sline =~ s/VRR=[\d\.]+/VRR=$vrr1,$ivrr/;
                        $sline =~ s/CN=loss-[\d\.]+/CN=loss-$cn1,gain+$cn/;
                    }
                    else{
                        $sline =~ s/SVTYPE=$stype1/SVTYPE=INS,$stype1/;
                        $sline =~ s/SVLEN=$len1/SVLEN=$ilen,$len1/;
                        $sline =~ s/READS=\d+/READS=$iread,$read1/;
                        $sline =~ s/VRR=[\d\.]+/VRR=$ivrr,$vrr1/;
                        $sline =~ s/CN=loss-[\d\.]+/CN=gain+$cn,loss-$cn1/;
                    }
                    $sline =~ s/GT=H[TM]/GT=HT2/;
                }
                else{
                    next;
                }
            }
            ${${$sv2{$chr}}{$pos}}{'TR'} = $sline;
            delete ${${$sv2{$chr}}{$ipos}}{$type};
            delete ${$ins_sv{$chr}}{$ipos};
        }
    }
}

foreach my $chr (keys %dup_sv){     # merge DUP overlapping TR with TR-CNV
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$dup_sv{$chr}}){
        my $end = ${$dup_sv{$chr}}{$pos};
        my $len = $end - $pos + 1;
        my $Mbin = int ($pos / $Mbin_size);
        next if (!exists ${$repeat{$chr2}}{$Mbin});
        foreach my $spos (sort {$a <=> $b} keys %{${$repeat{$chr2}}{$Mbin}}){
            last if ($spos > $end + $ins_bp_diff);
            my $send = ${${$repeat{$chr2}}{$Mbin}}{$spos};
            next if ($send < $pos - $ins_bp_diff);
            my $slen = $send - $spos + 1;
            my $flag = 0;
            my $overlap = 0;
            if ((abs ($pos - $spos) <= $ins_bp_diff) and (abs ($end - $send) <= $ins_bp_diff)){
                $flag = 1;
            }
            elsif (($spos <= $pos) and ($send >= $end)){
                $overlap = $slen;
            }
            elsif (($spos >= $pos) and ($spos <= $end)){
                $overlap = $end - $spos + 1;
                $overlap = $slen if ($send < $end);
            }
            elsif (($send >= $pos) and ($send <= $end)){
                $overlap = $send - $pos + 1;
                $overlap = $slen if ($spos > $pos);
            }
            if (($overlap >= $len * 0.8) and ($overlap >= $slen * 0.8)){
                $flag = 1;
            }
            if ($flag == 1){
                my $dupline = ${${$sv2{$chr}}{$pos}}{'DUP'};
                my $strid = ${$STR{$chr2}}{$spos};
                my ($spos2, $send2, $strulen) = split (/=/, $STR2{$strid});
                my $strlen = $send2 - $spos2 + 1;
                my $CN = $1 if ($dupline =~ /CN=([\d\.]+)/);
                my $expected_slen = $strlen * ($CN - 2);
                $expected_slen = $strlen if ($expected_slen == 0);
                my $cn = int ($expected_slen / $strulen * 10 + 0.5) / 10;
                my $read = $1 if ($dupline =~ /READS=(\d+)/);
                my $vrr = $1 if ($dupline =~ /VRR=([\d\.]+)/);
                my $sar = $1 if ($dupline =~ /SAR=([\d\.]+)/);
                my $gt = $1 if ($dupline =~ /GT=([^;]+)/);
                my $dpr = $1 if ($dupline =~ /DPR=([\d\.]+)/);
                if (!exists ${$str_sv{$chr}}{$spos}){
                    my $line = "$chr2\t$spos\t.\t.\t<TR:CNV>\t.\tPASS\tSVTYPE=INS;SVLEN=$len;READS=$read;CN=gain+$cn;VRR=$vrr;SAR=$sar;GT=$gt;BP=$read;TRDUP=$dpr;END=$spos;TRID=$strid;TREND=$send2;TRULEN=$strulen";
                    ${${$sv2{$chr}}{$spos}}{'TR'} = $line;
                }
                else{
                    my $sline = ${${$sv2{$chr}}{$spos}}{'TR'};
                    my $stype1 = '';
                    my $stype2 = '';
                    my $read1 = 0;
                    my $read2 = 0;
                    my $vrr1 = 0;
                    my $vrr2 = 0;
                    my $len1 = 0;
                    my $len2 = 0;
                    my $cn1 = 0;
                    my $cn2 = 0;
                    my $gt1 = '';
                    my $stype = $1 if ($sline =~ /SVTYPE=(.+?);/);
                    $stype1 = $stype if ($stype !~ /,/);
                    if ($stype =~ /(.+),(.+)/){
                        $stype1 = $1;
                        $stype2 = $2;
                    }
                    $read1 = $1 if ($sline =~ /READS=(\d+)/);
                    $read2 = $1 if ($sline =~ /READS=\d+,(\d+)/);
                    $vrr1 = $1 if ($sline =~ /VRR=([\d\.]+)/);
                    $vrr2 = $1 if ($sline =~ /VRR=[\d\.]+,([\d\.]+)/);
                    $len1 = $1 if ($sline =~ /SVLEN=(\d+)/);
                    $len2 = $1 if ($sline =~ /SVLEN=\d+,(\d+)/);
                    $cn1 = $1 if ($sline =~ /CN=gain\+([\d\.]+)/) and ($stype1 eq 'INS');
                    $cn1 = $1 if ($sline =~ /CN=loss-([\d\.]+)/) and ($stype1 eq 'DEL');
                    $cn2 = $1 if ($sline =~ /CN=gain\+[\d\.]+,gain\+([\d\.]+)/) and ($stype1 eq 'INS') and ($len2 > 0);
                    $cn2 = $1 if ($sline =~ /CN=loss-[\d\.]+,gain\+([\d\.]+)/) and ($stype1 eq 'DEL') and ($stype2 eq 'INS');
                    $gt1 = $1 if ($sline =~ /GT=(.+?);/);
                    my $flag2 = 0;
                    if ($stype1 eq 'INS'){
                        if (($len1 / $len >= 0.8) and ($len1 / $len <= 1.25)){
                            $read1 += $read;
                            $vrr1 += $vrr;
                            $vrr1 = 1 if ($vrr1 > 1);
                            $gt1 = 'HM' if ($gt1 eq 'HT') and ($vrr1 >= 0.8) and ($vrr2 == 0);
                            $flag2 = 1;
                            $sline =~ s/READS=\d+/READS=$read1/;
                            $sline =~ s/VRR=[\d\.]+/VRR=$vrr1/;
                            $sline =~ s/GT=HT/GT=$gt1/ if ($sline =~ /GT=HT/) and ($gt1 eq 'HM');
                            ${${$sv2{$chr}}{$spos}}{'TR'} = $sline;
                        }
                    }
                    if (($stype1 eq 'DEL') and ($stype2 eq 'INS') and ($flag2 == 0)){
                        if (($len2 / $len >= 0.8) and ($len2 / $len <= 1.25)){
                            $read2 += $read;
                            $vrr2 += $vrr;
                            $vrr2 = 1 if ($vrr2 > 1);
                            $flag2 = 1;
                            if ($vrr2 > $vrr1){
                                $sline =~ s/SVTYPE=DEL,INS/SVTYPE=INS,DEL/;
                                $sline =~ s/READS=\d+,\d+/READS=$read2,$read1/;
                                $sline =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$vrr2,$vrr1/;
                                $sline =~ s/SVLEN=\d+,\d+/SVLEN=$len2,$len1/;
                                $sline =~ s/CN=loss-[\d\.]+,gain\+[\d\.]+/CN=gain\+$cn2,loss-$cn1/;
                                ${${$sv2{$chr}}{$spos}}{'TR'} = $sline;
                            }
                            else{
                                $sline =~ s/READS=\d+,\d+/READS=$read1,$read2/;
                                $sline =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$vrr1,$vrr2/;
                                ${${$sv2{$chr}}{$spos}}{'TR'} = $sline;
                            }
                        }
                    }
                    elsif (($stype1 eq 'INS') and ($len2 > 0) and ($flag2 == 0)){
                        if (($len2 / $len >= 0.8) and ($len2 / $len <= 1.25)){
                            $read2 += $read;
                            $vrr2 += $vrr;
                            $vrr2 = 1 if ($vrr2 > 1);
                            $flag2 = 1;
                            if ($vrr2 > $vrr1){
                                $sline =~ s/READS=\d+,\d+/READS=$read2,$read1/;
                                $sline =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$vrr2,$vrr1/;
                                $sline =~ s/SVLEN=\d+,\d+/SVLEN=$len2,$len1/;
                                $sline =~ s/CN=gain\+[\d\.]+,gain\+[\d\.]+/CN=gain\+$cn2,gain\+$cn1/;
                                ${${$sv2{$chr}}{$spos}}{'TR'} = $sline;
                            }
                            else{
                                $sline =~ s/READS=\d+,\d+/READS=$read1,$read2/;
                                $sline =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$vrr1,$vrr2/;
                                ${${$sv2{$chr}}{$spos}}{'TR'} = $sline;
                            }
                        }
                    }
                }
                delete ${${$sv2{$chr}}{$pos}}{'DUP'};
                last;
            }
        }
    }
}

foreach my $chr (sort keys %str_sv){            # remove TR-CNV with < min_tr_vrr VRR
    foreach my $pos (sort {$a <=> $b} keys %{$str_sv{$chr}}){
        next if (!exists ${${$sv2{$chr}}{$pos}}{'TR'});
        my $sline = ${${$sv2{$chr}}{$pos}}{'TR'};
        my $vrr1 = $1 if ($sline =~ /VRR=([\d\.]+)/);
        my $len1 = $1 if ($sline =~ /SVLEN=(\d+)/);
        if (($len1 <= 3) and ($vrr1 < $min_str_vrr2)){
            delete ${${$sv2{$chr}}{$pos}}{'TR'};
        }
        elsif ($vrr1 < $min_str_vrr){
            delete ${${$sv2{$chr}}{$pos}}{'TR'};
        }
        elsif ($sline =~ /VRR=[\d\.]+,([\d\.]+)/){
            my $vrr2 = $1;
            my $len2 = $1 if ($sline =~ /SVLEN=\d+,(\d+)/);
            if ((($len2 > 3) and ($vrr2 < 0.1)) or (($len2 <= 3) and ($vrr2 <= $min_str_vrr2))){
                my $type = $1 if ($sline =~ /SVTYPE=(.+?);/);
                my $type1 = $type;
                if ($type =~ /,/){
                    ($type1) = split (/,/, $type);
                    $sline =~ s/SVTYPE=$type/SVTYPE=$type1/;
                }
                my $len1 = $1 if ($sline =~ /SVLEN=(\d+)/);
                my $read1 = $1 if ($sline =~ /READS=(\d+)/);
                my $cn1 = 0;
                $cn1 = $1 if ($sline =~ /CN=gain\+([\d\.]+)/) and ($type1 eq 'INS');
                $cn1 = $1 if ($sline =~ /CN=loss-([\d\.]+)/) and ($type1 eq 'DEL');
                if ($cn1 == 0){
                    my $ulen = $1 if ($sline =~ /TRULEN=(\d+)/);
                    $cn1 = int ($len1 / $ulen * 10 + 0.5) / 10;
                }
                $sline =~ s/SVLEN=\d+,\d+/SVLEN=$len1/;
                $sline =~ s/READS=\d+,\d+/READS=$read1/;
                $sline =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$vrr1/;
                $sline =~ s/CN=.+?;/CN=gain+$cn1;/ if ($type1 eq 'INS');
                $sline =~ s/CN=.+?;/CN=loss-$cn1;/ if ($type1 eq 'DEL');
                $sline =~ s/GT=HT2/GT=HT/ if ($vrr1 < 0.75) and ($sline =~ /GT=HT2/);
                $sline =~ s/GT=HT2/GT=HM/ if ($vrr1 >= 0.75) and ($sline =~ /GT=HT2/);
                ${${$sv2{$chr}}{$pos}}{'TR'} = $sline;
            }
        }
    }
}

my %gap;
my $pre_gchr = '';
my $pre_gend = 0;
my $pre_gpos = 0;

if (-f $gap_bed){
    open (FILE, $gap_bed) or die "$gap_bed is not found: $!\n" if ($gap_bed !~ /\.gz$/);
    open (FILE, "gzip -dc $gap_bed |") or die "$gap_bed is not found: $!\n" if ($gap_bed =~ /\.gz$/);
    while (my $line = <FILE>){
        chomp $line;
        my @line = split (/\s+/, $line);
        my $chr = $line[0];
        my $pos = $line[1];
        my $end = $line[2];
        if (($chr eq $pre_gchr) and ($pos - $pre_gend <= 1000)){
            ${$gap{$chr}}{$pre_gpos} = $end;
            $pre_gend = $end;
            next;
        }
        ${$gap{$chr}}{$pos} = $end;
        $pre_gchr = $chr;
        $pre_gend = $end;
        $pre_gpos = $pos;
    }
    close (FILE);
}

my $excluded_sv = 0;
my %lowconf_sv;
my %lowqual_sv;

foreach my $chr (sort keys %sv2){
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$sv2{$chr}}){
        foreach my $type (keys %{${$sv2{$chr}}{$pos}}){
            my $svline = ${${$sv2{$chr}}{$pos}}{$type};
            my $len = 0;
            $len = $1 if ($svline =~ /SVLEN=-*(\d+)/);
            my $end = $pos + $len - 1;
            $end = $pos if ($type eq 'INS');
            my $read = $1 if ($svline =~ /READS=([\d,]+)/);
            my $vrr = $1 if ($svline =~ /VRR=([\d\.,]+)/);
            if ($type eq 'TR'){
                if ($read !~ /,/){
                    $vrr = $1 if ($vrr =~ /(.+),/);
                    if (($read < $min_str_reads) or ($vrr < $min_str_vrr)){
                        delete ${${$sv2{$chr}}{$pos}}{$type};
                        delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                        next;
                    }
                    if (($read < $min_str_reads) or (($len <= 3) and ($vrr < $min_str_vrr2))){
                        delete ${${$sv2{$chr}}{$pos}}{$type};
                        delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                        next;
                    }
                }
                else{
                    my ($read1, $read2) = split (/,/, $read);
                    if (($read1 < $min_str_reads) and ($read2 < $min_str_reads)){
                        delete ${${$sv2{$chr}}{$pos}}{$type};
                        delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                        next;
                    }
                    if ($vrr =~ /,/){
                        my ($vrr1, $vrr2) = split (/,/, $vrr);
                        my ($len1, $len2) = split (/,/, $len);
                        if (($vrr1 < $min_str_vrr) and ($vrr2 < $min_str_vrr)){
                            delete ${${$sv2{$chr}}{$pos}}{$type};
                            delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                            next;
                        }
                        if (($len1 <= 3) and ($vrr1 < $min_str_vrr2) and ($len2 <= 3) and ($vrr2 < $min_str_vrr2)){
                            delete ${${$sv2{$chr}}{$pos}}{$type};
                            delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                            next;
                        }
                        if (($len1 <= 3) and ($vrr1 < $min_str_vrr2) and ($vrr2 < $min_str_vrr)){
                            delete ${${$sv2{$chr}}{$pos}}{$type};
                            delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                            next;
                        }
                        if (($vrr1 < $min_str_vrr) and ($len2 <= 3) and ($vrr2 < $min_str_vrr2)){
                            delete ${${$sv2{$chr}}{$pos}}{$type};
                            delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                            next;
                        }
                    }
                    else{
                        if ($vrr < $min_str_vrr){
                            delete ${${$sv2{$chr}}{$pos}}{$type};
                            delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                            next;
                        }
                        if (($vrr < $min_str_vrr2) and ($len !~ /,/) and ($len <= 3)){
                            delete ${${$sv2{$chr}}{$pos}}{$type};
                            delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                            next;
                        }
                    }
                }
            }
            else{
                if (($type eq 'INS') and ($read < $min_ins_reads)){
                    delete ${${$sv2{$chr}}{$pos}}{$type};
                    delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                    next;
                }
                elsif (($type eq 'INS') and ($len < 10) and ($vrr < $min_VRR2)){
                    delete ${${$sv2{$chr}}{$pos}}{$type};
                    delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                    next;
                }
                elsif (($type eq 'DEL') and ($read < $min_del_reads)){
                    delete ${${$sv2{$chr}}{$pos}}{$type};
                    delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                    next;
                }
                elsif (($type eq 'DEL') and ($len < 10) and ($vrr < $min_VRR2)){
                    delete ${${$sv2{$chr}}{$pos}}{$type};
                    delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                    next;
                }
                if (($len > 1) and ($len < $min_indel_size)){
                    delete ${${$sv2{$chr}}{$pos}}{$type};
                    delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                    next;
                }
                my $gap_overlap = 0;
                foreach my $gstart (sort {$a <=> $b} keys %{$gap{$chr}}){
                    my $gend = ${$gap{$chr}}{$gstart};
                    if (($pos >= $gstart) and ($pos <= $gend)){
                        if ($type eq 'INS'){
                            $gap_overlap = 1;
                            last;
                        }
                        else{
                            my $overlap = $gend - $pos + 1 if ($end >= $gend);
                            $overlap = $end - $pos + 1 if ($end < $gend);
                            if ($overlap >= $len * 0.5){
                                $gap_overlap = 1;
                                last;
                            }
                        }
                    }
                    elsif (($type ne 'INS') and ($gstart >= $pos) and ($gstart <= $end)){
                        my $overlap = $end - $gstart + 1 if ($end <= $gend);
                        $overlap = $gend - $gstart + 1 if ($end > $gend);
                        if ($overlap >= $len * 0.5){
                            $gap_overlap = 1;
                            last;
                        }
                    }
                    last if ($gstart > $end);
                }
                if ($gap_overlap == 1){
                    delete ${${$sv2{$chr}}{$pos}}{$type};
                    delete ${$sv2{$chr}}{$pos} if (scalar keys %{${$sv2{$chr}}{$pos}} == 0);
                    next;
                }
            }
            if (($type eq 'INS') and ($svline =~ /intDUP/)){
                my $intdup_chrpos = $1 if ($svline =~ /DUPPOS=(.+?);/);
                my $dup_len = $1 if ($svline =~ /DUPLEN=(\d+)/);
                my ($dup_chr, $dup_pos) = split (/:/, $intdup_chrpos);
                my $dup_end = $dup_pos + $dup_len - 1;
                if (exists $del_pos{$dup_chr}){
                    foreach my $dpos (sort {$a <=> $b} keys %{$del_pos{$dup_chr}}){
                        last if ($dpos > $dup_end);
                        my $dlen = ${$del_pos{$dup_chr}}{$dpos};
                        my $dend = $dpos + $dlen - 1;
                        next if ($dend < $dup_pos);
                        my $flag = 0;
                        if (($dpos <= $dup_pos) and ($dend >= $dup_end)){
                            if ($dup_len >= $dlen * $min_overlap_rate){
                                $flag = 1;
                            }
                        }
                        elsif (($dpos >= $dup_pos) and ($dpos <= $dup_end)){
                            my $overlap = $dup_end - $dpos + 1;
                            $overlap = $dlen if ($dend < $dup_end);
                            if (($overlap >= $dlen * $min_overlap_rate) and ($overlap >= $dup_len * $min_overlap_rate)){
                                $flag = 1;
                            }
                        }
                        elsif (($dend >= $dup_pos) and ($dend <= $dup_end)){
                            my $overlap = $dend - $dup_pos + 1;
                            $overlap = $dlen if ($dpos > $dup_pos);
                            if (($overlap >= $dlen * $min_overlap_rate) and ($overlap >= $dup_len * $min_overlap_rate)){
                                $flag = 1;
                            }
                        }
                        if ($flag == 1){
                            my @svline = split (/\t/, $svline);
                            $svline[4] = '<TRA:INS>';
                            $svline[7] =~ s/SVTYPE=$type/SVTYPE=TRA/;
                            $svline[7] =~ s/DUPPOS=/TRAPOS=/;
                            $svline[7] =~ s/DUPLEN=/TRALEN=/;
                            $svline = join ("\t", @svline);
                            $type = 'TRA';
                            ${${$sv2{$chr}}{$pos}}{$type} = $svline;
                            if (exists ${${$sv2{$dup_chr}}{$dpos}}{'DEL'}){
                                my $del_line = ${${$sv2{$dup_chr}}{$dpos}}{'DEL'};
                                my @del_line = split (/\t/, $del_line);
                                $del_line[4] = '<TRA:DEL>';
                                $del_line[7] =~ s/SVTYPE=DEL/SVTYPE=TRA/;
                                $del_line[7] .= ";TRAPOS=$chr:$pos;TRALEN=$len";
                                $del_line = join ("\t", @del_line);
                                delete ${${$sv2{$dup_chr}}{$dpos}}{'DEL'};
                                ${${$sv2{$dup_chr}}{$dpos}}{'TRA'} = $del_line;
                            }
                            last;
                        }
                    }
                }
            }
            if (($type ne 'TR') and (exists $lowconf_region{$chr2})){
                my $hit_flag = 0;
                foreach my $rpos (sort {$a <=> $b} keys %{$lowconf_region{$chr2}}){
                    last if ($rpos > $end);
                    my $rend = ${$lowconf_region{$chr2}}{$rpos};
                    next if ($rend < $pos);
                    my $rlen = $rend - $rpos + 1;
                    if ($type eq 'INS'){
                        $hit_flag = 1;
                        last;
                    }
                    else{
                        my $overlap = 0;
                        if (($rpos <= $pos) and ($rend >= $end)){
                            $hit_flag = 1;
                            last;
                        }
                        elsif (($rpos >= $pos) and ($rpos <= $end)){
                            $overlap = $end - $rpos + 1;
                            $overlap = $rlen if ($rend < $end);
                        }
                        elsif (($rend >= $pos) and ($rend <= $end)){
                            $overlap = $rend - $pos + 1;
                            $overlap = $rlen if ($rpos > $pos);
                        }
                        if ($overlap >= $len * $min_overlap_rate){
                            $hit_flag = 1;
                            last;
                        }
                    }
                }
                if ($hit_flag == 1){
                    my @svline = split (/\t/, $svline);
                    $svline[6] = 'LowConf';
                    $lowconf_sv{$type} ++;
                    $svline = join ("\t", @svline);
                    ${${$sv2{$chr}}{$pos}}{$type} = $svline;
                }
            }
        }
    }
}

foreach my $chr (sort keys %sv2){   # add uTR ID to the overlapping SV and mark LowConf for low-confident SVs within uTR region
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$sv2{$chr}}){
        foreach my $type (keys %{${$sv2{$chr}}{$pos}}){
            my $svline = ${${$sv2{$chr}}{$pos}}{$type};
            my @svline = split (/\t/, $svline);
#            next if ($svline =~ /LowConf/);
            my $len = 0;
            $len = $1 if ($svline =~ /SVLEN=-*(\d+)/);
            my $end = $pos + $len - 1;
            $end = $pos if ($type eq 'INS');
            my $sar = $1 if ($svline =~ /SAR=([\d\.]+)/);
            if (($type ne 'TR') and (exists $uSTR{$chr2})){
                my $hit_STR = '';
                my %hit_str;
                foreach my $rpos (sort {$a <=> $b} keys %{$uSTR{$chr2}}){
                    last if ($rpos > $end);
                    my $rend = ${$uSTR{$chr2}}{$rpos};
                    next if ($rend < $pos);
                    my $strid = ${$uSTR2{$chr2}}{$rpos};
                    my $rlen = $rend - $rpos + 1;
                    if ($type eq 'INS'){
                        $hit_STR = $strid;
                        last;
                    }
                    else{
                        my $overlap = 0;
                        if (($rpos <= $pos) and ($rend >= $end)){
                            $hit_STR = $strid;
                            last;
                        }
                        elsif (($rpos >= $pos) and ($rpos <= $end)){
                            $overlap = $end - $rpos + 1;
                            $overlap = $rlen if ($rend < $end);
                        }
                        elsif (($rend >= $pos) and ($rend <= $end)){
                            $overlap = $rend - $pos + 1;
                            $overlap = $rlen if ($rpos > $pos);
                        }
                        if ($overlap >= $len * $min_overlap_rate){
                            $hit_str{$strid} = $overlap;
                        }
                    }
                }
                if ((scalar keys %hit_str > 0) and ($hit_STR eq '')){
                    foreach my $strid (sort {$hit_str{$b} <=> $hit_str{$a}} keys %hit_str){
                        $hit_STR = $strid;
                        last;
                    }
                }
                if ($hit_STR ne ''){
                    $svline[7] .= ";TRID=$hit_STR";
                    if ((exists $lowconf_STR{$hit_STR}) and ($svline[6] ne 'LowConf')){
                        $svline[6] = 'LowConf';
                        $lowconf_sv{$type} ++;
                    }
                    $svline = join ("\t", @svline);
                    ${${$sv2{$chr}}{$pos}}{$type} = $svline;
                }
            }
            if ($sar > $max_SAR){
                if ($svline[6] eq 'PASS'){
                    $svline[6] = 'LowQual';
                    $lowqual_sv{$type} ++;
                    if ($type eq 'TR'){
                        $lowqual_sv{'TR2'} ++;
                        if ($svline[7] =~ /SVLEN=\d+,\d+/){
                            $lowqual_sv{'TR2'} ++;
                        }
                    }
                }
                else{
                    $svline[6] .= ',LowQual';
                }
                $svline = join ("\t", @svline);
                ${${$sv2{$chr}}{$pos}}{$type} = $svline;
            }
            if (($type eq 'INS') and ($len == 0) and ($svline[4] !~ /BP/)){
                $svline[4] =~ s/INS:/INS:BP:/;
                $svline = join ("\t", @svline);
                ${${$sv2{$chr}}{$pos}}{$type} = $svline;
            }
        }
    }
}

my $added_str_del = 0;
my $added_str_ins = 0;

foreach my $chr (sort keys %sv2){   # add TR-DEL located within large DELs and DUPs
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$sv2{$chr}}){
        foreach my $type (keys %{${$sv2{$chr}}{$pos}}){
            next if ($type eq 'INS') or ($type eq 'TR') or ($type eq 'INV');
            my $line = ${${$sv2{$chr}}{$pos}}{$type};
            my @line = split (/\t/, $line);
            next if ($line[6] ne 'PASS');
            next if ($line =~ /TRID=uTR/);
            my $len = $1 if ($line =~ /SVLEN=(\d+)/);
            my $end = $pos + $len - 1;
            my $vrr = $1 if ($line =~ /VRR=([\d\.]+)/);
            my $read = $1 if ($line =~ /READS=(\d+)/);
            my $gt = $1 if ($line =~ /GT=(.+?);/);
            my $sar = $1 if ($line =~ /SAR=([\d\.]+)/);
            if (($len >= 50) and ($vrr >= 0.15) and ($read >= 3) and ($sar < 0.3) and ($len <= 20000)){
                my $Mbin1 = int ($pos / $Mbin_size);
                my $Mbin2 = int ($end / $Mbin_size);
                my %hit_rpos;
                foreach my $rpos (sort {$a <=> $b} keys %{${$repeat{$chr2}}{$Mbin1}}){
                    last if ($rpos > $end);
                    my $rend = ${${$repeat{$chr2}}{$Mbin1}}{$rpos};
                    next if ($rend < $pos);
                    my $rlen = $rend - $rpos + 1;
                    if (($rpos >= $pos) and ($rend <= $end)){
                        $hit_rpos{$rpos} = $rlen;
                    }
                    elsif (($rpos >= $pos) and ($rpos <= $end)){
                        my $overlap = $end - $rpos + 1;
                        if ($overlap >= 50){
                            $hit_rpos{$rpos} = $overlap;
                        } 
                    }
                    elsif (($rend >= $pos) and ($rend <= $end)){
                        my $overlap = $rend - $pos + 1;
                        if ($overlap >= 50){
                            $hit_rpos{$rpos} = $overlap;
                        }
                    }
                }
                if ($Mbin2 > $Mbin1){
                    foreach my $rpos (sort {$a <=> $b} keys %{${$repeat{$chr2}}{$Mbin2}}){
                        last if ($rpos > $end);
                        my $rend = ${${$repeat{$chr2}}{$Mbin2}}{$rpos};
                        next if ($rend < $pos);
                        my $rlen = $rend - $rpos + 1;
                        if (($rpos >= $pos) and ($rend <= $end)){
                            $hit_rpos{$rpos} = 1;
                        }
                        elsif (($rpos >= $pos) and ($rpos <= $end)){
                            my $overlap = $end - $rpos + 1;
                            if ($overlap >= 50){
                                $hit_rpos{$rpos} = $overlap;
                            } 
                        }
                        elsif (($rend >= $pos) and ($rend <= $end)){
                            my $overlap = $rend - $pos + 1;
                            if ($overlap >= 50){
                                $hit_rpos{$rpos} = $overlap;
                            }
                        }
                    }
                }
                if (scalar keys %hit_rpos > 0){
                    if ($type eq 'DEL'){
                        foreach my $rpos (keys %hit_rpos){
                            my $del_len = $hit_rpos{$rpos};
                            if (exists ${${$sv2{$chr}}{$rpos}}{'TR'}){
                                my $str_line = ${${$sv2{$chr}}{$rpos}}{'TR'};
                                my $str_type = $1 if ($str_line =~ /SVTYPE=(.+?);/);
                                my $str_gt = $1 if ($str_line =~ /GT=(.+?);/);
                                next if ($str_gt eq 'HT2');
                                my @str_line = split (/\t/, $str_line);
                                my $str_ulen = $1 if ($str_line =~ /TRULEN=(\d+)/);
                                my $str_end = $1 if ($str_line =~ /TREND=(\d+)/);
                                my $del_cn = int ($del_len / $str_ulen * 10 + 0.5) / 10;
                                my $str_len = $1 if ($str_line =~ /SVLEN=(\d+)/);
                                my $str_vrr = $1 if ($str_line =~ /VRR=([\d\.]+)/);
                                my $str_read = $1 if ($str_line =~ /READS=(\d+)/);
                                my $str_cn = $1 if ($str_line =~ /CN=(.+?);/);
                                if (($str_type eq 'INS') and ($str_gt eq 'HT')){
                                    $str_line[7] =~ s/SVLEN=$str_len/SVLEN=$str_len,$del_len/;
                                    $str_line[7] =~ s/SVTYPE=INS/SVTYPE=INS,DEL/;
                                    $str_line[7] =~ s/CN=$str_cn/CN=$str_cn,loss-$del_cn/;
                                    $str_line[7] =~ s/VRR=$str_vrr/VRR=$str_vrr,$vrr/;
                                    $str_line[7] =~ s/READS=$str_read/READS=$str_read,$read/;
                                    $str_line[7] =~ s/GT=HT/GT=HT2/;
                                    $str_line = join ("\t", @str_line);
                                    $str_line .= ";ENCSV=$pos-$type-$len";
                                    ${${$sv2{$chr}}{$rpos}}{'TR'} = $str_line;
                                    $added_str_del ++;
                                }
                                elsif ($str_type eq 'DEL'){
                                    my $lenrate = int ($del_len / $str_len * 100 + 0.5) / 100;
                                    my $lenrate2 = int ($len / $str_len * 100 + 0.5) / 100;
                                    if (($lenrate2 >= 0.9) and ($lenrate2 <= 1.1)){
                                        $str_line[7] =~ s/SVLEN=$str_len/SVLEN=$del_len/;
                                        $str_line[7] =~ s/CN=$str_cn/CN=$del_cn/;
                                        $str_line = join ("\t", @str_line);
                                        $str_line .= ";ENCSV=$pos-$type-$len";
                                        ${${$sv2{$chr}}{$rpos}}{'TR'} = $str_line;
                                        $added_str_del ++;
                                    }
                                    elsif (($lenrate >= 2) and ($str_gt eq 'HT')){
                                        $str_line[7] =~ s/SVLEN=$str_len/SVLEN=$str_len,$del_len/;
                                        $str_line[7] =~ s/CN=$str_cn/CN=$str_cn,$del_cn/;
                                        $str_line[7] =~ s/VRR=$str_vrr/VRR=$str_vrr,$vrr/;
                                        $str_line[7] =~ s/READS=$str_read/READS=$str_read,$read/;
                                        $str_line[7] =~ s/GT=HT/GT=HT2/;
                                        $str_line = join ("\t", @str_line);
                                        $str_line .= ";ENCSV=$pos-$type-$len";
                                        ${${$sv2{$chr}}{$rpos}}{'TR'} = $str_line;
                                        $added_str_del ++;
                                    }
                                }
                            }
                            else{
                                my $strid = ${$STR{$chr2}}{$rpos};
                                my ($rpos2, $rend, $mlen) = split (/=/, $STR2{$strid});
                                my $cn = int ($del_len / $mlen * 10 + 0.5) / 10;
                                my $str_line = "$chr2\t$rpos\t.\t.\t<CNV:TR>\t.\tPASS\tSVTYPE=DEL;SVLEN=$del_len;READS=$read;CN=loss-$cn;VRR=$vrr;SAR=$sar;GT=$gt;END=$rend;TRID=$strid;TREND=$rend;TRULEN=$mlen;ENCSV=$pos-$type-$len";
                                ${${$sv2{$chr}}{$rpos}}{'TR'} = $str_line;
                                $added_str_del ++;
                            }
                        }
                    }
                    elsif ($type eq 'DUP'){
                        my $dpr = $1 if ($line =~ /DPR=([\d\.]+)/);
                        next if ($dpr < 1.2);
                        my $dup_cn = $1 if ($line =~ /CN=([\d\.]+)/);
                        $dup_cn -= 2;
                        $dup_cn = 1 if ($dup_cn < 1);
                        foreach my $rpos (keys %hit_rpos){
                            my $ovl_len = $hit_rpos{$rpos};
                            if (exists ${${$sv2{$chr}}{$rpos}}{'TR'}){
                                my $str_line = ${${$sv2{$chr}}{$rpos}}{'TR'};
                                my $str_type = $1 if ($str_line =~ /SVTYPE=(.+?);/);
                                my $str_gt = $1 if ($str_line =~ /GT=(.+?);/);
                                if ($str_gt eq 'HT'){
                                    my @str_line = split (/\t/, $str_line);
                                    my $str_ulen = $1 if ($str_line =~ /TRULEN=(\d+)/);
                                    my $str_end = $1 if ($str_line =~ /TREND=(\d+)/);
                                    my $ins_len = $ovl_len * $dup_cn;
                                    my $ins_cn = int ($ins_len / $str_ulen * 10 + 0.5) / 10 * $dup_cn;
                                    my $str_len = $1 if ($str_line =~ /SVLEN=(\d+)/);
                                    my $str_cn = $1 if ($str_line =~ /CN=(.+?);/);
                                    my $str_vrr = $1 if ($str_line =~ /VRR=([\d\.]+)/);
                                    my $str_read = $1 if ($str_line =~ /READS=(\d+)/);
                                    if ($str_type eq 'DEL'){
                                        $str_line[7] =~ s/SVTYPE=DEL/SVTYPE=DEL,INS/;
                                        $str_line[7] =~ s/SVLEN=$str_len/SVLEN=$str_len,$ins_len/;
                                        $str_line[7] =~ s/CN=$str_cn/CN=$str_cn,gain+$ins_cn/;
                                        $str_line[7] =~ s/VRR=$str_vrr/VRR=$str_vrr,$vrr/;
                                        $str_line[7] =~ s/READS=$str_read/READS=$str_read,$read/;
                                        $str_line[7] =~ s/GT=HT/GT=HT2/;
                                    }
                                    elsif ($str_type eq 'INS'){
                                        my $lenrate = int ($ins_len / $str_len * 100 + 0.5) / 100;
                                        if ($lenrate >= 2){
                                            $str_line[7] =~ s/SVLEN=$str_len/SVLEN=$str_len,$ins_len/;
                                            $str_line[7] =~ s/CN=$str_cn/CN=$str_cn,$ins_cn/;
                                            $str_line[7] =~ s/VRR=$str_vrr/VRR=$str_vrr,$vrr/;
                                            $str_line[7] =~ s/READS=$str_read/READS=$str_read,$read/;
                                            $str_line[7] =~ s/GT=HT/GT=HT2/;
                                        }
                                    }
                                    $str_line = join ("\t", @str_line);
                                    $str_line .= ";ENCSV=$pos-$type-$len";
                                    ${${$sv2{$chr}}{$rpos}}{'TR'} = $str_line;
                                    $added_str_ins ++;
                                }
                            }
                            else{
                                my $strid = ${$STR{$chr2}}{$rpos};
                                my ($rpos2, $rend, $mlen) = split (/=/, $STR2{$strid});
                                my $ins_len = $ovl_len * $dup_cn;
                                my $cn = int ($ins_len / $mlen * 10 + 0.5) / 10 * $dup_cn;
                                my $str_line = "$chr2\t$rpos\t.\t.\t<CNV:TR>\t.\tPASS\tSVTYPE=INS;SVLEN=$ins_len;READS=$read;CN=gain+$cn;VRR=$vrr;SAR=$sar;GT=$gt;END=$rend;TRID=$strid;TREND=$rend;TRULEN=$mlen;ENCSV=$pos-$type-$len";
                                ${${$sv2{$chr}}{$rpos}}{'TR'} = $str_line;
                                $added_str_ins ++;
                            }
                        }
                    }
                }
            }
        }
    }
}


my $out_vcf = "$out_prefix.discov.vcf";
open (OUT, "> $out_vcf");
print OUT "##INFO=<ID=SVTYPE,Number=.,Type=String,Description=\"Type of tandem repeat expansion/contraction (TR-CNV) and structural variation (SV)\">\n";
print OUT "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles (0 when undefined)\">\n";
print OUT "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variantdescribed in this record\">\n";
print OUT "##INFO=<ID=BPLEN,Number=1,Type=Integer,Description=\"Breakpoint distance of INS:BP\">\n";
print OUT "##INFO=<ID=DUPPOS,Number=1,Type=Integer,Description=\"Reference position of tandem duplication found between INS sequence and reference\">\n";
print OUT "##INFO=<ID=DUPLEN,Number=1,Type=Integer,Description=\"Length of tandem duplication between INS and reference\">\n";
print OUT "##INFO=<ID=TRAPOS,Number=1,Type=Integer,Description=\"Reference chr:pos of translocated segment\">\n";
print OUT "##INFO=<ID=TRALEN,Number=1,Type=Integer,Description=\"Length of translocated segment\">\n";
print OUT "##INFO=<ID=MEI,Number=1,Type=String,Description=\"Type of mobile element found in INS\">\n";
print OUT "##INFO=<ID=MEILEN,Number=1,Type=Integer,Description=\"Length of mobile element found in INS\">\n";
print OUT "##INFO=<ID=MEICN,Number=1,Type=Integer,Description=\"Extra copy numver of mobile element found in INS\">\n";
print OUT "##INFO=<ID=TRID,Number=1,Type=String,Description=\"TR ID containing TR-CNV\">\n";
print OUT "##INFO=<ID=TREND,Number=1,Type=Integer,Description=\"END position of TR-ID\">\n";
print OUT "##INFO=<ID=TRULEN,Number=1,Type=Integer,Description=\"Length of repeat unit of TR-ID\">\n";
print OUT "##INFO=<ID=TRDUP,Number=1,Type=Float,Description=\"Duplication of the entire TR region (value: DPR of DUP)\">\n";
print OUT "##INFO=<ID=TRUNIT,Number=1,Type=String,Description=\"motif:copy-number detected in INS sequence (only for SVTYPE=INS)\">\n";
print OUT "##INFO=<ID=READS,Number=1,Type=Integer,Description=\"Number of reads supporting the TR-CNV/SV\">\n";
print OUT "##INFO=<ID=BP,Number=1,Type=Integer,Description=\"Number of split-reads (breakpoints) supporting the TR-CNV/SV\">\n";
print OUT "##INFO=<ID=VRR,Number=1,Type=Float,Description=\"Ratio of TR-CNV/SV-supporting reads to read depth at the site\">\n";
print OUT "##INFO=<ID=DPR,Number=1,Type=Float,Description=\"Ratio of read depth in DEL/DUP region to that to the flanking regions (only non-TR-CNV)\">\n";
print OUT "##INFO=<ID=SAR,Number=1,Type=Float,Description=\"Ratio of TR-CNV/SV-supporting reads with mapping quality 0 to total supporting reads, including secondary alignments\">\n";
print OUT "##INFO=<ID=CN,Number=.,Type=Float,Description=\"Copy number of TR repeat unit in TR (gain/loss) or copy number in DUP\">\n";
print OUT "##INFO=<ID=ENCSV,Number=1,Type=String,Description=\"Position-Type-Size of DUP or DEL encompassing the TR region to be assigned as TR-INS or TR-DEL\">\n";
print OUT "##INFO=<ID=GT,Number=1,Type=String,Description=\"Genotype of SV (HT/HT2/HM/NA, HT2: multiallelic genotype)\">\n";

print OUT "##FILTER=<ID=LowConf,Description=\"Repeat region where TR-CNVs/SVs could be unreliably called\">\n";
print OUT "##FILTER=<ID=LowQual,Description=\"Low quality TR-CNVs/SVs with > $max_SAR SAR\">\n";

print OUT "##ALT=<ID=TR:CNV,Description=\"INS(repeat unit expansion) and/or DEL (repeat unit contraction) within Tandem Repeat (TR) regions\">\n";
print OUT "##ALT=<ID=DEL,Description=\"Deletion\">\n";
print OUT "##ALT=<ID=DUP,Description=\"Duplication\">\n";
print OUT "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">\n";
print OUT "##ALT=<ID=INS:DUP,Description=\"Insertion of tandem duplication\">\n";
print OUT "##ALT=<ID=INS:DUP:R,Description=\"Insertion of inverse tandem duplication\">\n";
print OUT "##ALT=<ID=INS:intDUP,Description=\"Insertion of interspersed tandem duplication\">\n";
print OUT "##ALT=<ID=INS:intDUP:R,Description=\"Insertion of interspersed inverse tandem duplication\">\n";
print OUT "##ALT=<ID=INS:BP,Description=\"Insertion with only breakpoints\">\n";
print OUT "##ALT=<ID=INS:BP:DUP,Description=\"Insertion with only breakpoints and with a potentially duplication segment\">\n";
print OUT "##ALT=<ID=INS:BP:DEL,Description=\"Insertion with only breakpoints and with a potentially deletion segment\">\n";

print OUT "##ALT=<ID=INS:ME:ALU,Description=\"Insertion of ALU element\">\n";
print OUT "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">\n";
print OUT "##ALT=<ID=INS:ME:SVA,Description=\"Insertion of SVA element\">\n";
print OUT "##ALT=<ID=INS:ME:HERVK,Description=\"Insertion of HERVK element\">\n";

print OUT "##ALT=<ID=INV:Description=\"Inversion\">\n";
print OUT "##ALT=<ID=TRA:INS:Description=\"Inserted translocation\">\n";
print OUT "##ALT=<ID=TRA:DEL:Description=\"Deleted translocation\">\n";
print OUT "##ALT=<ID=REP:Description=\"Replacement\">\n";
foreach my $chr (sort keys %sv2){
    foreach my $pos (sort {$a <=> $b} keys %{$sv2{$chr}}){
        foreach my $type (keys %{${$sv2{$chr}}{$pos}}){
            my $svline = ${${$sv2{$chr}}{$pos}}{$type};
            print OUT "$svline\n";
            $svtype{$type} ++;
            if ($type eq 'TR'){
                $svtype{'TR2'} ++;
                if ($svline =~ /SVLEN=\d+,\d+/){
                    $svtype{'TR2'} ++;
                }
            }
        }
    }
}
close (OUT);

my $out_ins = "$out_prefix.discov.INS.fa";
open (OUT, "> $out_ins");
foreach my $chr (sort keys %INS){
    foreach my $pos (sort {$a <=> $b} keys %{$INS{$chr}}){
        print OUT ${$INS{$chr}}{$pos};
    }
}
close (OUT);

print STDERR "#Newly added TR-INSs from DUPs: $added_str_ins\n";
print STDERR "#Newly added TR-DELs from DELs: $added_str_del\n";
print STDERR "#TR-CNVs/SVs removed from the excluded regions: $excluded_sv\n" if ($exclude_bed ne '');

my $total_calls = 0;
my $total_high_calls = 0;

print STDERR "\n#Summary of called SVs\n";
print STDERR "#Type\tNumber\tNumber(excluding LowConf&LowQual)\n";
foreach my $type (sort keys %svtype){
    my $type2 = $type;
    $type2 = 'CNV-TR-site' if ($type eq 'TR');
    $type2 = 'CNV-TR-allele' if ($type eq 'TR2');
    my $num = $svtype{$type};
    my $lowconf = 0;
    my $lowqual = 0;
    $lowconf = $lowconf_sv{$type} if (exists $lowconf_sv{$type});
    $lowqual = $lowqual_sv{$type} if (exists $lowqual_sv{$type});
    my $num2 = $num - $lowconf - $lowqual;
    $total_calls += $num if ($type ne 'TR');
    $total_high_calls += $num2 if ($type ne 'TR');
    print STDERR "$type2\t$num\t$num2\n";
}
print STDERR "Total\t$total_calls\t$total_high_calls\n";


sub run_step1{
    my ($config, $chr) = @_;
    system ("$Bin/TRsv_hifi_step1.pl $config $chr") if ($min_indel_size > 3);
    system ("$Bin/TRsv_hifi_ml1_step1.pl $config $chr") if ($min_indel_size <= 3);
    my $out_file = "$temp_dir/$out_prefix.chr$chr.discov.txt";
    $out_file = "$temp_dir/$out_prefix.$chr.discov.txt" if ($chr !~ /^[\dXY]+$/);
    threads->yield();
    sleep 1;

    return ($out_file, $chr);
}

sub run_step2{
    my ($config, $chr, $step1_out) = @_;

    system ("$Bin/TRsv_hifi_step2.pl $config $chr $step1_out");
    my $out_file = "$temp_dir/$out_prefix.chr$chr.discov.out";
    $out_file = "$temp_dir/$out_prefix.$chr.discov.out" if ($chr !~ /^[\dXY]+$/);

    threads->yield();
    sleep 1;

    return ($out_file, $chr);
}

sub run_step3{
    my ($config, $chr, $step2_out) = @_;

    system ("$Bin/TRsv_step3.pl $config $chr $step2_out");
    my $out_file = "$temp_dir/$out_prefix.chr$chr.discov.out";
    $out_file = "$temp_dir/$out_prefix.$chr.discov.out" if ($chr !~ /^[\dXY]+$/);

    threads->yield();
    sleep 1;

    return ($out_file, $chr);
}


