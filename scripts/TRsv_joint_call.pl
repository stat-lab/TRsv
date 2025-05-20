#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long; 
use Pod::Usage;
use FindBin qw($Bin);
use File::Basename;

my $vcf_list = '';

my $out_prefix = 'TRsv.joint';

my $out_dir = '.';

my $target_chr = 'ALL';

my $min_sv_len = 50;
my $min_str_len = 0;

my $ins_sd = 150;

my $min_overlap_ratio = 0.5;
my $min_overlap_ratio2 = 0.8;

my $max_SAR = 0.6;

my $pos_ave = 1;

my $non_human = 0;

my $build = 37;

my $ref_index = '';

my $data_dir = "$Bin/../Data";

my $gap_bed = '';

my $gender_list = '';

my $lowQ_sample_list = '';

my $help;

GetOptions(
    'vcf_list|v=s' => \$vcf_list,
    'prefix|p=s' => \$out_prefix,
    'out_dir|od=s' => \$out_dir,
    'target_chr|c=s' => \$target_chr,
    'gender|g=s' => \$gender_list,
    'gap_bed|gap=s' => \$gap_bed,
    'non_human|nh=i' => \$non_human,
    'build|b=i' => \$build,
    'ref|r=s' => \$ref_index,
    'lqs_list|lqs=s' => \$lowQ_sample_list,
    'min_len|ml=i' => \$min_sv_len,
    'min_slen|msl=i' => \$min_str_len,
    'max_sar|sar=f' => \$max_SAR,
    'ins_sd|is=i' => \$ins_sd,
    'overlap_rate|or=f' => \$min_overlap_ratio2,
    'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  TRsv_joint_call.pl -v <vcf list> -p <out_prefix> -od <output directory> -r <ref index> (-gap <gap bed -nh 1 if samples are non-human sepecies)

  Options:
   --vcf_list or -v <STR>   a table file of sample names and corresponding vcf files. Vcf files should be absolue path or relative path from the working directory [mandatory]
   --ref or -r <STR>        index file of reference fasta (ref.fa.fai generated with samtools faidx) [mandatory for non human sepecies]
							
   --prefix or -p <STR>     prefix name of an output joint called vcf file [default: LRsv.joint]
   --out_dir or -od <STR>   output directory [default: ./]
   --target_chr or -c <STR> target chromosome [default: ALL]
   --non_human or -nh <INT> samples are non-human species (0: human, 1: non-human) [default: 0]
   --build or -b <STR>      human reference build (GRCh37, GRCh38) number (37, 38, or T2T, only effective for human) [default: 37]
   --gap or -gap <STR>      gap bed file, indicating gap regions in reference genome. may be specified for non-human species [default for human: Data/gap.bed or gap.b38,bed]
   --gender or -g <STR>     sample name-gender table file to exclude chrY for female (sample_name and M/F, separated with tab in each line) [optional]
   --lqs_list or -lqs <STR> list file indicating sample names of low quality long read data, such as PacBio CLR and ONT (singleton DELs and INSs, with < 100 bp and < 0.2 VRR, of the corresponding sample are removed) [optional]
   --min_len or -ml <INT>   minimum len (bp) of SV, except for INS-BP [default: 50]
   --min_slen or -msl <INT> minimum len (bp) of TR-CNV [default: N spcified with -ml]
   --max_sar or -sar <FLOAT> maximum rate of SVs supported by alignments with maping quality 0, including secondary alignmnents (SVs exceeding this value are marked as 'LowQual' in the FILTER field) [default: 0.6]
   --ins_sd or -is <INT>    maximum distance (bp) between proximal INS breakpoints to be merged [default: 200]
   --overlap_rate or -or <FLOAT>  minimum reciprocal overlap rate for merging overlapped DELs, DUPs, or INVs [default: 0.8]
   --help or -h             output help message
   
=cut

die "sample-name/vcf list file is not specified: \n" if ($vcf_list eq '');

system ("mkdir $out_dir") if ($out_dir ne '.') and (!-d $out_dir);

my $temp_dir = "$out_dir/temp";
system ("mkdir $temp_dir") if (!-d $temp_dir);

my $var_sd = 150;
my $mei_sd = 150;

if ($non_human == 0){
	if ($gap_bed eq ''){
		$gap_bed = "$data_dir/gap.bed";
		$gap_bed = "$data_dir/gap.b38.bed" if ($build eq '38');
	}
	if ($ref_index eq ''){
		$ref_index = "$data_dir/hs37.fa.fai";
		$ref_index = "$data_dir/hs38.fa.fai" if ($build eq '38');
		$ref_index = "$data_dir/chm13v2.0.fa.fai" if ($build eq 'T2T');
	}
}
print STDERR "Gap-bed: $gap_bed\tRef-index: $ref_index\n";

$min_str_len = $min_sv_len if ($min_str_len == 0);

my %gap;

if ($gap_bed ne ''){
	my $pre_gchr = '';
	my $pre_gend = 0;
	my $pre_gpos = 0;
	open (FILE, $gap_bed) or die "$gap_bed is not found: $!\n";
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

my @ref;
my %ref_len;
open (FILE, $ref_index) or die "$ref_index is not found:\n";
while (my $line = <FILE>){
	chomp $line;
	my ($chr, $len) = split (/\t/, $line);
	next if ($chr eq 'M') or ($chr eq 'chrM') or ($chr eq 'Mit') or ($chr eq 'Mt');
	last if ($chr !~ /^c*h*r*[\dXY]+$/) and ($non_human == 0);
	push @ref, $chr;
	$ref_len{$chr} = $len; 
}
close (FILE);

my %sample_id;
my %Gid;
my %Gid_order;
my %gender;
my %lowQ_sample;
my @var_file;
my @sample_id;

if ($gender_list ne ''){
    open (FILE, $gender_list) or die "$gender_list is not found: $!\n";
    while (my $line = <FILE>){
        chomp $line;
        next if ($line =~ /^#|^$/);
        my ($ID, $sex) = split (/\s+/, $line);
        $sex = 'F' if ($sex eq 'Female') or ($sex eq 'female');
        $sex = 'M' if ($sex eq 'Male') or ($sex eq 'male');
        $gender{$ID} = $sex;
    }
    close (FILE);
}

if ($lowQ_sample_list ne ''){
	open (FILE, $lowQ_sample_list) or die "$lowQ_sample_list is not found: $!\n";
	while (my $line = <FILE>){
		chomp $line;
		next if ($line =~ /^#|^$/);
		my ($ID) = split (/\t/, $line);
		$lowQ_sample{$ID} = 1;
	}
	close (FILE);
}

open (FILE, $vcf_list) or die "$vcf_list is not found: $!\n";
while (my $line = <FILE>){
	chomp $line;
	next if ($line =~ /^#|^$/);
	my ($ID, $vcf) = split (/\s+/, $line);
	push @sample_id, $ID;
	push @var_file, $vcf;
	$sample_id{$ID} = $vcf;
	$Gid{$ID} = 1;
}
close (FILE);

my $sample_count = 9;
foreach my $id (sort keys %Gid){
	$sample_count ++;
	$Gid_order{$sample_count} = $id;
}

my $total_sample_num = scalar @sample_id;
print STDERR "Total samples: $total_sample_num\n";

my %call;
my %call_line;
my %cnv_line;
my %cnv_info;
my %call_cons;
my %call_cons_len;
my %used_pos;
my %call_ave_len;
my %call_diverged;
my %call_diverged_ave;
my %call_samples;
my %ins_bp;
my %ins_bp_line;
my $count = 0;
my @header;
my $total_str2 = 0;

foreach my $id (@sample_id){
    my $var_file = $sample_id{$id};
	$count ++;
    open (FILE, $var_file) or die "$var_file is not found: $!\n";
    while (my $line = <FILE>){
		chomp $line;
		next if ($line =~ /^#/);
		my @line = split (/\t/, $line);
		my $chr = $line[0];
		next if ($target_chr ne 'ALL') and ($chr ne $target_chr);
		next if (exists $gender{$id}) and ($gender{$id} eq 'F') and ($chr eq 'Y');
		my $pos = $line[1];
		my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
		if ($line[4] =~ /CNV/){
			my $len1 = 0;
			my $len2 = 0;
			$len1 = $1 if ($line[7] =~ /SVLEN=(\d+),\d+/);
			$len2 = $1 if ($line[7] =~ /SVLEN=\d+,(\d+)/);
			$len1 = $1 if ($line[7] =~ /SVLEN=(\d+);/);
			my $strid = $1 if ($line[7] =~ /TRID=(.+?);/);
			my $strend = $1 if ($line[7] =~ /TREND=(\d+)/);
			my $strulen = $1 if ($line[7] =~ /TRULEN=(\d+)/);
			if (($type eq 'DEL') or ($type eq 'INS')){
				if (($len1 < $min_str_len) and ($len2 < $min_str_len)){
					next;
				}
				elsif (($len1 >= $min_str_len) and ($len2 > 0) and ($len2 < $min_str_len)){
					$line[7] =~ s/SVLEN=$len1,$len2/SVLEN=$len1/;
					if ($line[7] =~ /READS=(\d+),\d+/){
						my $r = $1;
						$line[7] =~ s/READS=\d+,\d+/READS=$r/;
					}
					if ($line[7] =~ /VRR=([\d\.]+),[\d\.]+/){
						my $v = $1;
						$line[7] =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$v/;
					}
					if (($type eq 'DEL') and ($line[7] =~ /CN=loss\-([\d\.]+),loss\-[\d\.]+/)){
						my $c = $1;
						$line[7] =~ s/CN=loss\-[\d\.]+,loss\-[\d\.]+/CN=loss-$c/;
					}
					elsif (($type eq 'INS') and ($line[7] =~ /CN=gain\+([\d\.]+),gain\+[\d\.]+/)){
						my $c = $1;
						$line[7] =~ s/CN=gain\+[\d\.]+,gain\+[\d\.]+/CN=gain+$c/;
					}
				}
				elsif (($len1 < $min_str_len) and ($len2 >= $min_str_len)){
					$line[7] =~ s/SVLEN=$len1,$len2/SVLEN=$len2/;
					if ($line[7] =~ /READS=\d+,(\d+)/){
						my $r = $1;
						$line[7] =~ s/READS=\d+,\d+/READS=$r/;
					}
					if ($line[7] =~ /VRR=[\d\.]+,([\d\.]+)/){
						my $v = $1;
						$line[7] =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$v/;
					}
					if (($type eq 'DEL') and ($line[7] =~ /CN=loss\-[\d\.]+,loss\-([\d\.]+)/)){
						my $c = $1;
						$line[7] =~ s/CN=loss\-[\d\.]+,loss\-[\d\.]+/CN=loss-$c/;
					}
					elsif (($type eq 'INS') and ($line[7] =~ /CN=gain\+[\d\.]+,gain\+([\d\.]+)/)){
						my $c = $1;
						$line[7] =~ s/CN=gain\+[\d\.]+,gain\+[\d\.]+/CN=gain+$c/;
					}
				}
				if (($type eq 'DEL') and ($line[7] =~ /CN=loss\-([\d\.]+),([\d\.]+)/)){
					my $cn1 = $1;
					my $cn2 = $2;
					$line[7] =~ s/CN=loss\-[\d\.]+,[\d\.]+/CN=loss-$cn1,loss-$cn2/;
				}
				elsif (($type eq 'DEL') and ($line[7] =~ /CN=([\d\.]+),([\d\.]+)/)){
					my $cn1 = $1;
					my $cn2 = $2;
					$line[7] =~ s/CN=loss\-[\d\.]+,[\d\.]+/CN=loss-$cn1,loss-$cn2/;
				}
				elsif (($type eq 'INS') and ($line[7] =~ /CN=gain\+([\d\.]+),([\d\.]+)/)){
					my $cn1 = $1;
					my $cn2 = $2;
					$line[7] =~ s/CN=gain\+[\d\.]+,[\d\.]+/CN=gain+$cn1,gain+$cn2/;
				}
			}
			elsif ($type eq 'DEL,INS'){
				if (($len1 < $min_str_len) and ($len2 < $min_str_len)){
					next;
				}
				elsif (($len1 >= $min_str_len) and ($len2 > 0) and ($len2 < $min_str_len)){
					$line[7] =~ s/SVLEN=$len1,$len2/SVLEN=$len1/;
					$line[7] =~ s/SVTYPE=DEL,INS/SVTYPE=DEL/;
					if ($line[7] =~ /READS=(\d+),\d+/){
						my $r = $1;
						$line[7] =~ s/READS=\d+,\d+/READS=$r/;
					}
					if ($line[7] =~ /VRR=([\d\.]+),[\d\.]+/){
						my $v = $1;
						$line[7] =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$v/;
					}
					if (($line[7] =~ /CN=loss\-([\d\.]+),gain\+[\d\.]+/)){
						my $c = $1;
						$line[7] =~ s/CN=loss\-[\d\.]+,gain\+[\d\.]+/CN=loss-$c/;
					}
				}
				elsif (($len1 < $min_str_len) and ($len2 >= $min_str_len)){
					$line[7] =~ s/SVLEN=$len1,$len2/SVLEN=$len2/;
					$line[7] =~ s/SVTYPE=DEL,INS/SVTYPE=INS/;
					if ($line[7] =~ /READS=\d+,(\d+)/){
						my $r = $1;
						$line[7] =~ s/READS=\d+,\d+/READS=$r/;
					}
					if ($line[7] =~ /VRR=[\d\.]+,([\d\.]+)/){
						my $v = $1;
						$line[7] =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$v/;
					}
					if (($line[7] =~ /CN=loss\-[\d\.]+,gain\+([\d\.]+)/)){
						my $c = $1;
						$line[7] =~ s/CN=loss\-[\d\.]+,gain\+[\d\.]+/CN=gain+$c/;
					}
				}
				if ($line[7] =~ /CN=loss\-([\d\.]+);/){
					my $cn1 = $1;
					my $strulen = $1 if ($line[7] =~ /TRULEN=(\d+)/);
					my $cn2 = int ($len2 / $strulen * 10 + 0.5) / 10;
					$line[7] =~ s/CN=loss\-[\d\.]+;/CN=loss-$cn1,gain+$cn2;/;
				}
			}
			elsif ($type eq 'INS,DEL'){
				if (($len1 < $min_str_len) and ($len2 < $min_str_len)){
					next;
				}
				elsif (($len1 >= $min_str_len) and ($len2 > 0) and ($len2 < $min_str_len)){
					$line[7] =~ s/SVLEN=$len1,$len2/SVLEN=$len1/;
					$line[7] =~ s/SVTYPE=INS,DEL/SVTYPE=INS/;
					if ($line[7] =~ /READS=(\d+),\d+/){
						my $r = $1;
						$line[7] =~ s/READS=\d+,\d+/READS=$r/;
					}
					if ($line[7] =~ /VRR=([\d\.]+),[\d\.]+/){
						my $v = $1;
						$line[7] =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$v/;
					}
					if (($line[7] =~ /CN=gain\+([\d\.]+),loss\-[\d\.]+/)){
						my $c = $1;
						$line[7] =~ s/CN=gain\+[\d\.]+,loss\-[\d\.]+/CN=gain+$c/;
					}
				}
				elsif (($len1 < $min_str_len) and ($len2 >= $min_str_len)){
					$line[7] =~ s/SVLEN=$len1,$len2/SVLEN=$len2/;
					$line[7] =~ s/SVTYPE=INS,DEL/SVTYPE=DEL/;
					if ($line[7] =~ /READS=\d+,(\d+)/){
						my $r = $1;
						$line[7] =~ s/READS=\d+,\d+/READS=$r/;
					}
					if ($line[7] =~ /VRR=[\d\.]+,([\d\.]+)/){
						my $v = $1;
						$line[7] =~ s/VRR=[\d\.]+,[\d\.]+/VRR=$v/;
					}
					if (($line[7] =~ /CN=gain\+[\d\.]+,loss\-([\d\.]+)/)){
						my $c = $1;
						$line[7] =~ s/CN=gain\+[\d\.]+,loss\-[\d\.]+/CN=loss-$c/;
					}
				}
				if ($line[7] =~ /CN=gain\+([\d\.]+);/){
					my $cn1 = $1;
					my $strulen = $1 if ($line[7] =~ /TRULEN=(\d+)/);
					my $cn2 = int ($len2 / $strulen * 10 + 0.5) / 10;
					$line[7] =~ s/CN=gain\+[\d\.]+;/CN=gain+$cn1,loss-$cn2;/;
				}
			}
			$line[4] = '<TR:CNV>';
			$line = join ("\t", @line);
			${${$cnv_line{$chr}}{$pos}}{$id} = $line;
			$cnv_info{$strid} = "$chr=$pos=$strend=$strulen";
		}
		else{
			my $len = 0;
			$len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
			next if ($len < $min_sv_len) and ($len > 0);
			my $end = $1 if ($line[7] =~ /END=(\d+)/);
			my $bp_distance = 0;
			$bp_distance = $1 if ($line[7] =~ /BPLEN=(\d+)/);
			${${${$call{$type}}{$chr}}{$pos}}{$id} = $len if ($len > 0);
			${${${$call_line{$type}}{$chr}}{$pos}}{$id} = "$line[4]==$line[7]~$line[6]" if ($len > 0);
			${${$ins_bp{$chr}}{$pos}}{$id} = $bp_distance if ($len == 0) and ($type eq 'INS');
			${${$ins_bp_line{$chr}}{$pos}}{$id} = "$line[4]==$line[7]~$line[6]" if ($len == 0) and ($type eq 'INS');
		}
		$call_samples{$id} ++;
    }
    close (FILE);
}

my $sum_calls = 0;
my $ave_call = 0;
foreach my $id (keys %call_samples){
	$sum_calls += $call_samples{$id};
}
$ave_call = int ($sum_calls / $count);

foreach my $type (keys %call){			# clustering DEL/DUP/INV exhibiting > 2-fold larger or smaller size of the averaged length at the same called positions
    foreach my $chr (keys %{$call{$type}}){
		foreach my $pos (sort {$a <=> $b} keys %{${$call{$type}}{$chr}}){
		    my @len = ();
		    foreach my $id (keys %{${${$call{$type}}{$chr}}{$pos}}){
				my $len = ${${${$call{$type}}{$chr}}{$pos}}{$id};
				push @len, $len;
				push @{${${$call_cons{$type}}{$chr}}{$pos}}, "$id==$pos==$len==${${${$call_line{$type}}{$chr}}{$pos}}{$id}";
			}
			if (@len == 1){
				${${$call_ave_len{$type}}{$chr}}{$pos} = $len[0];
		    }
		    elsif (@len > 1){
				my %clust_len;
				my %clust_ave;
				my %clust_cont;
				my $clust_num = 1;
				my $clust_len = 0;
				my $count_len = 0;
				foreach (sort {$a <=> $b} @len){
				    if ($clust_len == 0){
						push @{$clust_len{$clust_num}}, $_;
						$clust_len = $_;
				    }
				    elsif (($clust_len > $_ * 0.556) or ($_ <= 3)){
						push @{$clust_len{$clust_num}}, $_;
				    }
				    else{
						$clust_num ++;
						$clust_len = $_;
						push @{$clust_len{$clust_num}}, $_;
				    }
				}
				foreach my $clustn (sort {$a <=> $b} keys %clust_len){
				    my $sum_len = 0;
				    my $clust_cont = scalar @{$clust_len{$clustn}};
				    foreach my $len1 (@{$clust_len{$clustn}}){
						$sum_len += $len1;
				    }
				    my $ave_len = int ($sum_len / $clust_cont);
				    $clust_ave{$clustn} = $ave_len;
				    $clust_cont{$clustn} = $clust_cont;
				}
				my $pre_ave = 0;
				my $pre_clustn = 0;
				foreach my $clustn (sort {$a <=> $b} keys %clust_ave){
				    my $avelen = $clust_ave{$clustn};
				    if ($pre_ave == 0){
						$pre_ave = $avelen;
						$pre_clustn = $clustn;
						next;
				    }
				    if (($avelen < $pre_ave * 1.8) or (($pre_ave == 1) and ($avelen <= 3))){
						foreach my $len1 (@{$clust_len{$clustn}}){
						    push @{$clust_len{$pre_clustn}}, $len1;
						}
						my $new_avelen = int (($pre_ave * $clust_cont{$pre_clustn} + $avelen * $clust_cont{$clustn}) / ($clust_cont{$pre_clustn} + $clust_cont{$clustn}));
						$clust_cont{$pre_clustn} += $clust_cont{$clustn};
						delete $clust_len{$clustn};
						delete $clust_cont{$clustn};
						$clust_ave{$pre_clustn} = $new_avelen;
						$pre_ave = $new_avelen;
				    }
				    else{
						$pre_ave = $avelen;
						$pre_clustn = $clustn;
				    }
				}
				my $top_clust = 0;
				my %top_clust;
				my %low_clust;
				foreach my $clustn (sort {$clust_cont{$b} <=> $clust_cont{$a}} keys %clust_cont){
				    $top_clust = $clustn;
				    foreach my $len1 (@{$clust_len{$clustn}}){
						$top_clust{$len1} = $clustn;
				    }
				    last;
				}
				my $new_clust = 0;
				foreach my $clustn (sort {$a <=> $b} keys %clust_len){
				    next if ($clustn == $top_clust);
				    $new_clust ++;
				    foreach my $len1 (@{$clust_len{$clustn}}){
						$low_clust{$len1} = $new_clust;
				    }
				    ${${${$call_diverged_ave{$type}}{$chr}}{$pos}}{$new_clust} = $clust_ave{$clustn};
				}
				if (scalar keys %clust_len > 1){
				    my @info = (@{${${$call_cons{$type}}{$chr}}{$pos}});
				    delete ${${$call_cons{$type}}{$chr}}{$pos};
				    foreach (@info){
						my ($id1, $pos1, $len1) = split (/==/, $_);
						if (exists $top_clust{$len1}){
						    push @{${${$call_cons{$type}}{$chr}}{$pos}}, $_;
						    ${${$call_ave_len{$type}}{$chr}}{$pos} = $clust_ave{$top_clust};
						}
						else{
						    my $clustn = $low_clust{$len1};
						    push @{${${${$call_diverged{$type}}{$chr}}{$pos}}{$clustn}}, $_;
						}
				    }
				}
				else{
				    my $sum_len = 0;
				    map{$sum_len += $_} @len;
				    my $ave_len2 = int ($sum_len / @len);
				    ${${$call_ave_len{$type}}{$chr}}{$pos} = $ave_len2;
				}
		    }
		}
    }
}
%call = ();

foreach my $type (keys %call_diverged){			# reassign the above extracted diverged sites to the most appropriate clusetr
    foreach my $chr (keys %{$call_diverged{$type}}){
		foreach my $pos (sort {$a <=> $b} keys %{${$call_diverged{$type}}{$chr}}){
		    foreach my $clustn (sort {$a <=> $b} keys %{${${$call_diverged{$type}}{$chr}}{$pos}}){
				my $avelen = ${${${$call_diverged_ave{$type}}{$chr}}{$pos}}{$clustn};
				my $end = $pos + $avelen - 1;
				my %len_diff;
				my $added_pos = 0;
				foreach my $pos2 (sort {$a <=> $b} keys %{${$call_ave_len{$type}}{$chr}}){
				    my $len2 = ${${$call_ave_len{$type}}{$chr}}{$pos2};
				    my $end2 = $pos2 + $len2 - 1;
				    next if ($end2 < $pos);
				    last if ($pos2 > $end);
				    my $diff = abs ($avelen - $len2);
				    if ($type eq 'INS'){
				    	if (((abs ($pos - $pos2) <= 50) and ($avelen >= 50)) or ((abs ($pos - $pos2) <= $avelen) and ($avelen < 50))){
					    	if (($len2 / $avelen > 0.556) and ($len2 / $avelen < 1.8)){
								$len_diff{$pos2} = $diff;
							}
						}
				    }
				    else{
					    if (($pos <= $pos2) and ($end >= $end2)){
							if ($len2 > $avelen * 0.556){
							    $len_diff{$pos2} = $diff;
							}
					    }
					    elsif (($pos >= $pos2) and ($end <= $end2)){
							if ($avelen > $len2 * 0.556){
							    $len_diff{$pos2} = $diff;
							}
					    }
					    elsif (($pos <= $end2) and ($end >= $end2)){
							if (($end2 - $pos > $avelen * 0.556) and ($end2 - $pos > $len2 * 0.556)){
							    $len_diff{$pos2} = $diff;
							}
					    }
					    elsif (($pos <= $pos2) and ($end >= $pos2)){
							if (($end - $pos2 > $avelen * 0.556) and ($end - $pos2 > $len2 * 0.556)){
							    $len_diff{$pos2} = $diff;
							}
					    }
					    elsif (($avelen <= 2) and ($len2 <= 2) and (($end == $pos2 - 1) or ($end2 == $pos - 1))){
					    	$len_diff{$pos2} = $diff;
					    }
					}
				}
				if (scalar keys %len_diff > 0){
				    foreach my $pos3 (sort {$len_diff{$a} <=> $len_diff{$b}} keys %len_diff){
						push @{${${$call_cons{$type}}{$chr}}{$pos3}}, @{${${${$call_diverged{$type}}{$chr}}{$pos}}{$clustn}};
						$added_pos = $pos3;
						last;
				    }
				}
				else{
				    my $new_pos = $pos;
				    my $flag = 0;
				    while ($flag == 0){
						if (!exists ${${$call_cons{$type}}{$chr}}{$new_pos}){
						    push @{${${$call_cons{$type}}{$chr}}{$new_pos}}, @{${${${$call_diverged{$type}}{$chr}}{$pos}}{$clustn}};
						    $added_pos = $new_pos;
						    $flag = 1;
						}
						else{
						    $new_pos ++;
						}
				    }
				}
				if ($added_pos != 0){
				    my $sum_pos = 0;
				    my $sum_len = 0;
				    my $item_num = 0;
				    my @info = (@{${${$call_cons{$type}}{$chr}}{$added_pos}});
				    foreach (@{${${$call_cons{$type}}{$chr}}{$added_pos}}){
						my ($id1, $pos1, $len1) = split (/==/, $_);
						$sum_pos += $pos1;
						$sum_len += $len1;
						$item_num ++;
				    }
				    my $ave_pos = int ($sum_pos / $item_num);
				    my $ave_len = int ($sum_len / $item_num);
				    delete ${${$call_cons{$type}}{$chr}}{$added_pos};
				    delete ${${$call_ave_len{$type}}{$chr}}{$added_pos};
				    if ((!exists ${${$call_cons{$type}}{$chr}}{$ave_pos}) or ($added_pos == $ave_pos)){
						push @{${${$call_cons{$type}}{$chr}}{$ave_pos}}, @info;
						${${$call_ave_len{$type}}{$chr}}{$ave_pos} = $ave_len;
				    }
				    else{
						my $new_pos = $ave_pos;
						my $flag = 0;
						while ($flag == 0){
						    if (!exists ${${$call_cons{$type}}{$chr}}{$new_pos}){
								push @{${${$call_cons{$type}}{$chr}}{$new_pos}}, @info;
								${${$call_ave_len{$type}}{$chr}}{$new_pos} = $ave_len;
								$flag = 1;
						    }
						    else{
								$new_pos ++;
						    }
						}
				    }
				}
		    }
		}
    }
}
%call_ave_len = ();
%call_diverged = ();
%call_diverged_ave = ();
print STDERR "1st step completed:\n";

foreach my $type (keys %call_cons){	# merge posiotions present within 30-bp of consensus pos
    foreach my $chr (keys %{$call_cons{$type}}){
		foreach my $pos (sort {$a <=> $b} keys %{${$call_cons{$type}}{$chr}}){
		    next if (!exists ${${$call_cons{$type}}{$chr}}{$pos});
		    next if (exists ${${$used_pos{$type}}{$chr}}{$pos});
		    my $id_num = scalar @{${${$call_cons{$type}}{$chr}}{$pos}};
		    my $sum_len = 0;
		    my $ave_len = 0;
			foreach my $info (@{${${$call_cons{$type}}{$chr}}{$pos}}){
			    my ($id, $pos1, $len1) = split (/==/, $info);
			    $sum_len += $len1;
			}
			$ave_len = int ($sum_len / $id_num);
			for (my $i = $pos + 1; $i <= $pos + 30; $i++){
				next if (exists ${${$used_pos{$type}}{$chr}}{$i});
				my $distance = $i - $pos;
				next if ($distance > $ave_len * 2);
				next if ($ave_len < 50) and ($distance > $ave_len);
				if (exists ${${$call_cons{$type}}{$chr}}{$i}){
					my $sum_leni = 0;
					foreach my $infoi (@{${${$call_cons{$type}}{$chr}}{$i}}){
					    my ($idi, $posi, $leni) = split (/==/, $infoi);
					    $sum_leni += $leni;
					}
					my $ave_leni = int ($sum_leni / scalar @{${${$call_cons{$type}}{$chr}}{$i}});
					next if ($ave_len > $ave_leni * 1.8) or ($ave_leni > $ave_len * 1.8);
					next if ($distance > $ave_leni * 2);
				    push @{${${$call_cons{$type}}{$chr}}{$pos}}, @{${${$call_cons{$type}}{$chr}}{$i}};
				    delete ${${$call_cons{$type}}{$chr}}{$i};
				    ${${$used_pos{$type}}{$chr}}{$i} = 1;
				}
		    }
		}
    }
}

foreach my $type (keys %call_cons){	# merge posiotions present within 100-bp of consensus pos
    foreach my $chr (keys %{$call_cons{$type}}){
		foreach my $pos (sort {$a <=> $b} keys %{${$call_cons{$type}}{$chr}}){
		    next if (!exists ${${$call_cons{$type}}{$chr}}{$pos});
		    next if (exists ${${$used_pos{$type}}{$chr}}{$pos});
		    my $id_num = scalar @{${${$call_cons{$type}}{$chr}}{$pos}};
		    my $sum_len = 0;
		    my $ave_len = 0;
			foreach my $info (@{${${$call_cons{$type}}{$chr}}{$pos}}){
			    my ($id, $pos1, $len1) = split (/==/, $info);
			    $sum_len += $len1;
			}
			$ave_len = int ($sum_len / $id_num);
		    my %ids;
		    foreach my $item (@{${${$call_cons{$type}}{$chr}}{$pos}}){
				my ($id) = split (/==/, $item);
				$ids{$id} = 1;
		    }
		    for (my $i = $pos + 1; $i <= $pos + 100; $i++){
				next if (exists ${${$used_pos{$type}}{$chr}}{$i});
				my $distance = $i - $pos;
				next if ($distance > $ave_len * 2);
				next if ($ave_len < 50) and ($distance > $ave_len);
				if (exists ${${$call_cons{$type}}{$chr}}{$i}){
				    my $flag = 0;
				    foreach my $item (@{${${$call_cons{$type}}{$chr}}{$i}}){		# next if samples at pos and pre_pos share a common sample id
						my ($id) = split (/==/, $item);
						if (exists $ids{$id}){
						    $flag = 1;
						    last;
						}
					}
					my $sum_leni = 0;
					foreach my $infoi (@{${${$call_cons{$type}}{$chr}}{$i}}){
					    my ($idi, $posi, $leni) = split (/==/, $infoi);
					    $sum_leni += $leni;
					}
					my $ave_leni = int ($sum_leni / scalar @{${${$call_cons{$type}}{$chr}}{$i}});
					next if ($ave_len > $ave_leni * 1.8) or ($ave_leni > $ave_len * 1.8);
					next if ($distance > $ave_leni * 2);
					if ($flag == 0){
						push @{${${$call_cons{$type}}{$chr}}{$pos}}, @{${${$call_cons{$type}}{$chr}}{$i}};
						delete ${${$call_cons{$type}}{$chr}}{$i};
						${${$used_pos{$type}}{$chr}}{$i} = 1;
					}
					else{
						foreach my $info (@{${${$call_cons{$type}}{$chr}}{$i}}){
			    			my ($id, $pos1, $len1) = split (/==/, $info);
			    			if (!exists $ids{$id}){
			    				push @{${${$call_cons{$type}}{$chr}}{$pos}}, $info;
			    			}
			    		}
			    		delete ${${$call_cons{$type}}{$chr}}{$i};
						${${$used_pos{$type}}{$chr}}{$i} = 1;
				    }
				}
		    }
		}
    }
}

foreach my $type (keys %call_cons){		# assign a median position from multiple sample-positions to the pos
    foreach my $chr (keys %{$call_cons{$type}}){
		foreach my $pos (sort {$a <=> $b} keys %{${$call_cons{$type}}{$chr}}){
		    my $sample_num = @{${${$call_cons{$type}}{$chr}}{$pos}};
		    my $median_pos = $pos;
		    my %id_pos;
		    foreach my $item (@{${${$call_cons{$type}}{$chr}}{$pos}}){
				my ($id, $pos2, $len2) = split (/==/, $item);
				my $id_pos = "$id=$pos2";
				$id_pos{$id_pos} = $pos2;
		    }
		    my $sample_num_2 = scalar keys %id_pos;
		    my $half_num = int ($sample_num_2 * 0.5 + 0.5);
		    my $count = 0;
		    if ($sample_num_2 > 2){
				foreach my $id_pos (sort {$id_pos{$a} <=> $id_pos{$b}} keys %id_pos){
				    $count ++;
				    if ($count == $half_num){
						my ($id, $pos2) = split (/=/, $id_pos);
						$median_pos = $pos2;
				    }
				}
		    }
		    if ($pos != $median_pos){
				push @{${${$call_cons{$type}}{$chr}}{$median_pos}}, @{${${$call_cons{$type}}{$chr}}{$pos}};
				delete ${${$call_cons{$type}}{$chr}}{$pos};
				${${$used_pos{$type}}{$chr}}{$pos} = 1;
				delete ${${$used_pos{$type}}{$chr}}{$median_pos} if (exists ${${$used_pos{$type}}{$chr}}{$median_pos});
		    }
		}
    }
}

foreach my $type (keys %call_cons){		# assign median SV length for each consensus pos
    foreach my $chr (keys %{$call_cons{$type}}){
		foreach my $pos (sort {$a <=> $b} keys %{${$call_cons{$type}}{$chr}}){
		    next if (exists ${${$used_pos{$type}}{$chr}}{$pos});
		    my @len;
		    my $sum_len = 0;
		    my $new_len = 0;
	        my $maxlen = 0;
		    my $num = scalar @{${${$call_cons{$type}}{$chr}}{$pos}};
		    foreach (@{${${$call_cons{$type}}{$chr}}{$pos}}){
				my ($id, $pos2, $len) = split (/==/, $_);
				if (($len == 0) or ($len == 1)){
				    $num --;
				    next;
				}
				$sum_len += $len;
				push @len, $len;
	            $maxlen = $len if ($maxlen < $len);
		    }
		    if ($type eq 'INS'){
                if ($maxlen >= $min_sv_len){
                    $new_len = $maxlen;
                }
                else{
                    $new_len = 0;
                }
            }
		    else{
				if ($num == 1){
				    $new_len = $sum_len;
				}
				elsif ($num == 2){
				    $new_len = int ($sum_len / $num + 0.5);
				}
				elsif ($num >= 3){
				    my $hn = int ($num * 0.5);
				    @len = sort {$a <=> $b} @len;
				    my $median = $len[$hn];
				    ($new_len) = &cons_len (\@len, $median);
				}
				if (($type eq 'TRA') and ($new_len == 0)){

				}
		    }
		    ${${$call_cons_len{$type}}{$chr}}{$pos} = $new_len;
		}
    }
}

foreach my $type (keys %call_cons){		# merge neighboring consensus pos with < 200 bp distance for >= 50 bp INS and reciprocal 50%-overlap for the other types of SVs
    foreach my $chr (keys %{$call_cons{$type}}){
		my %pre_info;
		foreach my $pos (sort {$a <=> $b} keys %{${$call_cons{$type}}{$chr}}){
		    next if (!exists ${${$call_cons{$type}}{$chr}}{$pos});
		    next if (exists ${${$used_pos{$type}}{$chr}}{$pos});
		    my $len = ${${$call_cons_len{$type}}{$chr}}{$pos};
		    next if ($len < 20);
		    my $end = $pos + $len - 1;
		    my %select_pos;
		    my @match;
		    foreach my $pos2 (sort {$b <=> $a} keys %pre_info){
				next if ($pos == $pos2);
				next if (exists ${${$used_pos{$type}}{$chr}}{$pos2});
				my $len2 = $pre_info{$pos2};
				next if ($len2 < 20);
				my $end2 = $pos2 + $len2 - 1;
				$end2 = $pos2 if ($type eq 'INS');
				last if ($pos - $pos2 > 20000000);
				next if ($end2 + $ins_sd < $pos);
				my $ovlrate = 0;
				if ($type eq 'INS'){
				    $end = $pos;
				    my $distance = abs ($pos - $pos2);
				    if (($len >= 50) and ($len2 >= 50)){
					    if (($distance <= $ins_sd) and ($distance < $len * 2) and ($distance < $len2 * 2)){
							$ovlrate = 1 / $distance if ($distance > 0);
							$ovlrate = 1 if ($distance == 0);
					    }
					}
					else{
						if (($distance <= $len) and ($distance <= $len2)){
							$ovlrate = 1 / $distance if ($distance > 0);
							$ovlrate = 1 if ($distance == 0);
					    }
					}
				}
				else{
				    my $overlap = 0;
					if (($pos2 <= $pos) and ($end2 >= $end)){
						$overlap = $len;
					}
					elsif (($pos2 >= $pos) and ($pos2 <= $end)){
						$overlap = $end - $pos2 + 1;
						$overlap = $len2 if ($end2 < $end);
					}
					elsif (($end2 >= $pos) and ($end2 <= $end)){
						$overlap = $end2 - $pos + 1;
						$overlap = $len2 if ($pos2 > $pos);
					}
					if (($overlap >= $len * $min_overlap_ratio2) and ($overlap >= $len2 * $min_overlap_ratio2)){
						if ($len >= $len2){
						    $ovlrate = $len2 / $len;
						}
						else{
						    $ovlrate = $len / $len2;
						}
					}
				}
				if ($ovlrate > 0){
				    push @match, "$pos2=$len2=$ovlrate";
				}
		    }
		    if (@match > 0){
				if (@match == 1){
				    my ($mpos, $mlen) = split (/=/, $match[0]);
				    $select_pos{$mpos} = $mlen if (exists ${${$call_cons{$type}}{$chr}}{$mpos});
				}
				else{
				    my $max_rate = 0;
				    my $best_pos = 0;
				    my $best_len = 0;
				    foreach (@match){
						my @info = split (/=/, $_);
						my $rate = $info[2];
						if ($rate > $max_rate){
						    $max_rate = $rate;
						    $best_pos = $info[0];
						    $best_len = $info[1];
						}
				    }
				    $select_pos{$best_pos} = $best_len if (exists ${${$call_cons{$type}}{$chr}}{$best_pos});
				}
		    }
		    if (scalar keys %select_pos > 0){		# check if samples at pos and pre_pos share a common sample id
				my @cur_info = ();
				my @pre_info = ();
				my @pre_pos = ();
				@cur_info = (@{${${$call_cons{$type}}{$chr}}{$pos}});
				foreach my $prepos (sort {$a <=> $b} keys %select_pos){
				    push @pre_info, @{${${$call_cons{$type}}{$chr}}{$prepos}};
				    push @pre_pos, $prepos;
				}
				my $overlap_flag = 0;
				my %ids;
				my %overlap_ids;
				my $overlap_num = 0;
				my $small_sample = 0;
				foreach my $item (@cur_info){
				    my ($id) = split (/==/, $item);
				    $ids{$id} = 1;
				}
				foreach my $item (@pre_info){
				    my ($pre_id) = split (/==/, $item);
				    if ((exists $ids{$pre_id}) and ($ids{$pre_id} == 1)){
						$ids{$pre_id} ++;
						$overlap_num ++;
						$overlap_ids{$pre_id} = 1;
				    }
				}
				if ($overlap_num <= 1){
				    $overlap_flag = 0;
				}
				elsif ($overlap_num > 1){
				    my @pos = ();
				    my @pre_pos2 = ();
				    my @overlap_pos = ();
				    my @overlap_pre_pos = ();
				    my @len = ();
				    my @pre_len = ();
				    my @overlap_len = ();
				    my @overlap_pre_len = ();
				    my @item;
				    my @pre_item;
				    my @overlap_item = ();
				    my @overlap_pre_item = ();
				    my $sum_pos = 0;
				    my $sum_pre_pos = 0;
				    my $sum_overlap_pos = 0;
				    my $sum_overlap_pre_pos = 0;
				    my $sum_len = 0;
				    my $sum_pre_len = 0;
				    my $sum_overlap_len = 0;
				    my $sum_overlap_pre_len = 0;
				    my %used_item;
				    my %used_pre_item;
				    foreach my $item (@cur_info){
						next if (exists $used_item{$item});
						my ($id, $pos1, $len1) = split (/==/, $item);
						if (exists $overlap_ids{$id}){
						    push @overlap_item, $item;
						    push @overlap_pos, $pos1;
						    push @overlap_len, $len1 if ($len1 >= $min_sv_len);
						}
						else{
						    push @item, $item;
						    push @pos, $pos1;
						    push @len, $len1 if ($len1 >= $min_sv_len);
						}
						$used_item{$item} = 1;
				    }
				    foreach my $item (@pre_info){
						next if (exists $used_pre_item{$item});
						my ($pre_id, $pre_pos1, $pre_len1) = split (/==/, $item);
						if (exists $ids{$pre_id}){
						    push @overlap_pre_item, $item;
						    push @overlap_pre_pos, $pre_pos1;
						    push @overlap_pre_len, $pre_len1 if ($pre_len1 >= $min_sv_len);
						}
						else{
						    push @pre_item, $item;
						    push @pre_pos2, $pre_pos1;
						    push @pre_len, $pre_len1 if ($pre_len1 >= $min_sv_len);
						}
						$used_pre_item{$item} = 1;
				    }
				    map{$sum_pos += $_} @pos;
				    map{$sum_pre_pos += $_} @pre_pos2;
				    map{$sum_overlap_pos += $_} @overlap_pos;
				    map{$sum_overlap_pre_pos += $_} @overlap_pre_pos;
				    map{$sum_len += $_} @len;
				    map{$sum_pre_len += $_} @pre_len;
				    map{$sum_overlap_len += $_} @overlap_len;
				    map{$sum_overlap_pre_len += $_} @overlap_pre_len;
				    my $ave_pos = 0;
				    my $ave_pre_pos = 0;
				    my $ave_overlap_pos = 0;
				    my $ave_overlap_pre_pos = 0;
				    my $ave_len = 0;
				    my $ave_pre_len = 0;
				    my $ave_overlap_len = 0;
				    my $ave_overlap_pre_len = 0;
				    $ave_pos = int ($sum_pos / @pos + 0.5) if (@pos > 0);
				    $ave_pre_pos = int ($sum_pre_pos / @pre_pos2 + 0.5) if (@pre_pos2 > 0);
				    $ave_overlap_pos = int ($sum_overlap_pos / @overlap_pos + 0.5) if (@overlap_pos > 0);
				    $ave_overlap_pre_pos = int ($sum_overlap_pre_pos / @overlap_pre_pos + 0.5) if (@overlap_pre_pos > 0);
				    $ave_len = int ($sum_len / @pos + 0.5) if (@len > 0);
				    $ave_pre_len = int ($sum_pre_len / @pre_len + 0.5) if (@pre_len > 0);
				    $ave_overlap_len = int ($sum_overlap_len / @overlap_len + 0.5) if (@overlap_len > 0);
				    $ave_overlap_pre_len = int ($sum_overlap_pre_len / @overlap_pre_len + 0.5) if (@overlap_pre_len > 0);
				    $ave_pos = $ave_overlap_pos if ($ave_pos == 0);
				    $ave_pre_pos = $ave_overlap_pre_pos if ($ave_pre_pos == 0);
				    $ave_len = $ave_overlap_len if ($ave_len == 0);
				    $ave_pre_len = $ave_overlap_pre_len if ($ave_pre_len == 0);
				    my $AB = abs ($ave_pos - $ave_pre_pos);
				    my $ovAovB = abs ($ave_overlap_pos - $ave_overlap_pre_pos);
				    my $AB_len = abs ($ave_len - $ave_pre_len);
				    my $ovAovB_len = abs ($ave_overlap_len - $ave_overlap_pre_len);
				    if ($type eq 'INS'){	# compare absolute distances between SV positions/lengths for overlapped and non-overlapped IDs
						if ($AB < $ovAovB){
						    $overlap_flag = 1;
						}
						else{
						    $overlap_flag = 0;
						}
				    }
				    else{
						if (($AB < $ovAovB) and ($AB_len < $ovAovB_len)){
						    $overlap_flag = 1;
						}
						else{
						    $overlap_flag = 0;
						}
				    }
				    if ($overlap_flag == 1){		# re-assign unoverlapped items of pre-pos and pos to overlapped pre-pos and overlapped pos
						foreach my $item (@item){
						    my ($id, $pos1, $len1) = split (/==/, $item);
						    if (abs ($pos1 - $ave_overlap_pos) <= abs ($pos1 - $ave_overlap_pre_pos)){
							push @overlap_item, $item;
						    }
						    else{
							push @overlap_pre_item, $item;
						    }
						}
						foreach my $item (@pre_item){
						    my ($id, $pos1, $len1) = split (/==/, $item);
						    if (abs ($pos1 - $ave_overlap_pos) <= abs ($pos1 - $ave_overlap_pre_pos)){
							push @overlap_item, $item;
						    }
						    else{
							push @overlap_pre_item, $item;
						    }
						}
						my $median_pos = 0;
						my $median_pre_pos = 0;
						my $median_len = 0;
						my $median_pre_len = 0;
						my @sum_pos = ();
						my @sum_pre_pos = ();
						my @sum_len = ();
						my @sum_pre_len = ();
						foreach my $item (@overlap_item){
						    my ($id, $pos1, $len1) = split (/==/, $item);
						    push @sum_pos, $pos1;
						    push @sum_len, $len1 if ($len1 >= $min_sv_len);
						}
						foreach my $item (@overlap_pre_item){
						    my ($id, $pos1, $len1) = split (/==/, $item);
						    push @sum_pre_pos, $pos1;
						    push @sum_pre_len, $len1 if ($len1 >= $min_sv_len);
						}
						my $half_num = int (@sum_pos * 0.5 + 0.5);
						my $half_pre_num = int (@sum_pre_pos * 0.5 + 0.5);
						my $count = 0;
						foreach my $pos1 (sort {$a <=> $b} @sum_pos){
						    $count ++;
						    if ($count == $half_num){
								$median_pos = $pos1;
								last;
						    }
						}
						$count = 0;
						foreach my $pos1 (sort {$a <=> $b} @sum_pre_pos){
						    $count ++;
						    if ($count == $half_pre_num){
								$median_pre_pos = $pos1;
								last;
						    }
						}
						$count = 0;
						foreach my $len1 (sort {$a <=> $b} @sum_len){
						    $count ++;
						    if ($count == $half_num){
								$median_len = $len1;
								last;
						    }
						}
						$count = 0;
						foreach my $len1 (sort {$a <=> $b} @sum_pre_len){
						    $count ++;
						    if ($count == $half_pre_num){
								$median_pre_len = $len1;
								last;
						    }
						}
						$median_pre_pos -- if ($median_pos == $median_pre_pos);
						delete ${${$call_cons{$type}}{$chr}}{$pos};
						map {delete ${${$call_cons{$type}}{$chr}}{$_}} @pre_pos;
						map {delete $pre_info{$_}} @pre_pos;
                        if (exists ${${$call_cons{$type}}{$chr}}{$median_pre_pos}){
                            while (1){
                                $median_pre_pos += 10;
                                last if (!exists ${${$call_cons{$type}}{$chr}}{$median_pre_pos});
                            }
                        }
                        if (exists ${${$call_cons{$type}}{$chr}}{$median_pos}){
                            while (1){
                                $median_pos += 10;
                                last if (!exists ${${$call_cons{$type}}{$chr}}{$median_pos});
                            }
                        }
						push @{${${$call_cons{$type}}{$chr}}{$median_pre_pos}}, @overlap_pre_item;
						${${$call_cons_len{$type}}{$chr}}{$median_pre_pos} = $median_pre_len;
						$pre_info{$median_pre_pos} = $median_pre_len;
						push @{${${$call_cons{$type}}{$chr}}{$median_pos}}, @overlap_item;
						${${$call_cons_len{$type}}{$chr}}{$median_pos} = $median_len;
						$pre_info{$median_pos} = $median_len;
						delete ${${$used_pos{$type}}{$chr}}{$median_pre_pos} if (exists ${${$used_pos{$type}}{$chr}}{$median_pre_pos});
						delete ${${$used_pos{$type}}{$chr}}{$median_pos} if (exists ${${$used_pos{$type}}{$chr}}{$median_pos});
						next;
				    }
				}
				if ($overlap_flag == 0){
				    my %pos;
				    my $sum_pos = 0;;
				    my $top_pos1 = 0;
				    my $top_pos2 = 0;
				    my $pos1_freq = 0;
				    my $pos2_freq = 0;
				    my $select_pos = 0;
				    my $sum_len2 = 0;
				    my $sum_len_num = 0;
				    my $ave_len2 = 0;
				    my $total_num = @cur_info + @pre_info;
				    delete ${${$call_cons{$type}}{$chr}}{$pos};
				    delete ${${$call_cons_len{$type}}{$chr}}{$pos};
				    map {delete ${${$call_cons{$type}}{$chr}}{$_}} @pre_pos;
				    map {delete ${${$call_cons_len{$type}}{$chr}}{$_}} @pre_pos;
				    map {delete $pre_info{$_}} @pre_pos;
					
				    foreach (@cur_info){
						my ($id, $bp, $len) = split (/==/, $_);
						$pos{$bp} ++;
						$sum_pos += $bp;
						$sum_len2 += $len if ($len >= $min_sv_len);
						$sum_len_num ++ if ($len >= $min_sv_len);
				    }
				    foreach (@pre_info){
						my ($id, $bp, $len) = split (/==/, $_);
						$pos{$bp} ++;
						$sum_pos += $bp;
						$sum_len2 += $len if ($len >= $min_sv_len);
						$sum_len_num ++ if ($len >= $min_sv_len);
				    }
				    $ave_len2 = int ($sum_len2 / $sum_len_num) if ($sum_len_num > 0);
				    
				    foreach my $pos2 (sort {$pos{$b} <=> $pos{$a}} keys %pos){
						$top_pos2 = $pos2 if ($top_pos1 > 0) and ($top_pos2 == 0);
						$pos2_freq = $pos{$pos2} if ($pos1_freq > 0) and ($pos2_freq == 0);
						$top_pos1 = $pos2 if ($top_pos1 == 0);
						$pos1_freq = $pos{$pos2} if ($pos1_freq == 0);
						last if ($top_pos2 > 0);
				    }
				    if ($total_num == 2){
						$select_pos = int (($top_pos1 + $top_pos2) / 2 + 0.5);
				    }
				    elsif ($total_num == 3){
						if ($pos1_freq >= 2){
						    $select_pos = $top_pos1;
						}
						else{
						    $select_pos = int ($sum_pos / 3 + 0.5);
						}
				    }
				    else{
						if (($pos1_freq >= 2) and ($pos1_freq >= $pos2_freq * 2)){
						    $select_pos = $top_pos1;
						}
						elsif ($pos2_freq >= 2){
						    $select_pos = int (($top_pos1 + $top_pos2) / 2 + 0.5);
						}
						else{
						    $select_pos = int ($sum_pos / $total_num + 0.5);
						}
				    }
                    if (exists ${${$call_cons{$type}}{$chr}}{$select_pos}){
                        while (1){
                            $select_pos += 10;
                            last if (!exists ${${$call_cons{$type}}{$chr}}{$select_pos});
                        }
                    }
				    push @{${${$call_cons{$type}}{$chr}}{$select_pos}}, @pre_info;
				    push @{${${$call_cons{$type}}{$chr}}{$select_pos}}, @cur_info;
				    ${${$call_cons_len{$type}}{$chr}}{$select_pos} = $ave_len2;
				    map {${${$used_pos{$type}}{$chr}}{$_} = 1} @pre_pos;
				    ${${$used_pos{$type}}{$chr}}{$pos} = 1;
				    delete ${${$used_pos{$type}}{$chr}}{$select_pos} if (exists ${${$used_pos{$type}}{$chr}}{$select_pos});
				    $pre_info{$select_pos} = $ave_len2;
				}
		    }
		    else{
				$pre_info{$pos} = $len;
		    }
		}
    }
}

print STDERR "2nd step completed:\n";

my %vcf_type;
my %vcf_cons;

foreach my $type (keys %call_cons){
    foreach my $chr (keys %{$call_cons{$type}}){
		my $chr02d = $chr;
		$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/) and ($chr !~ /^0/);
		foreach my $pos (sort {$a <=> $b} keys %{${$call_cons{$type}}{$chr}}){
		    next if (exists ${${$used_pos{$type}}{$chr}}{$pos});
		    next if (scalar @{${${$call_cons{$type}}{$chr}}{$pos}} == 0);
		    my %id_pos;
		    my %id_pos_2;
		    my $ID_pos = '';
		    my $subtype = '';
		    my $len = ${${$call_cons_len{$type}}{$chr}}{$pos};
		    my $end = $pos + $len - 1 if ($type ne 'INS');
		    $end = $pos if ($type eq 'INS');
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
		    next if ($gap_overlap == 1);
		    foreach my $item (@{${${$call_cons{$type}}{$chr}}{$pos}}){
				my ($id, $pos2, $len2, $alt, $info) = split (/==/, $item);
				${$id_pos{$id}}{$pos2} = "$len2==$alt==$info";
                $id_pos_2{$id} = $pos2;
		    }
            foreach my $id (sort keys %id_pos){
                foreach my $pos2 (sort {$a <=> $b} keys %{$id_pos{$id}}){
                    my $info2 = ${$id_pos{$id}}{$pos2};
                    my $id_pos = "$id==$pos2==$info2";
                    $ID_pos .= $id_pos . ',,';
                }
            }
            $ID_pos =~ s/,,$//;
            my $sample_num = scalar keys %id_pos;
            my $half_num = int ($sample_num * 0.5 + 0.5);
            my $count = 0;
            my $new_pos = $pos;
            if ($sample_num >= 3){
                foreach my $id (sort {$id_pos_2{$a} <=> $id_pos_2{$b}} keys %id_pos_2){		# determine median position of multiple sample-positions
                    $count ++;
                    if ($count == $half_num){
                        $new_pos = $id_pos_2{$id};
                        last;
                    }
                }
            }
		    ${${$vcf_type{$type}}{$chr02d}}{$new_pos} = "SVLEN2=$len;SAMPLES=$ID_pos";
		}
    }
}
%call_cons = ();
%call_cons_len = ();

my %ins_bp_cons;
my %ins_bp_len;
my %ins_bp_cons2;
my %ins_bp_convert;

foreach my $chr (keys %ins_bp){
	foreach my $pos (sort {$a <=> $b} keys %{$ins_bp{$chr}}){
		my @len;
		foreach my $id (sort keys %{${$ins_bp{$chr}}{$pos}}){
			my $len = ${${$ins_bp{$chr}}{$pos}}{$id};
			push @{${$ins_bp_cons{$chr}}{$pos}}, "$id==$pos==$len==${${$ins_bp_line{$chr}}{$pos}}{$id}";
			push @len, $len;
		}
		my $hn = int (@len * 0.5);
		@len = sort {$a <=> $b} @len;
		my $med_len = $len[$hn];
		${$ins_bp_len{$chr}}{$pos} = $med_len;
	}
}

foreach my $chr (keys %ins_bp_cons){
	my %delete;
	foreach my $pos1 (sort {$a <=> $b} keys %{$ins_bp_cons{$chr}}){
		next if (exists $delete{$pos1});
		my $len1 = ${$ins_bp_len{$chr}}{$pos1};
		my $end1 = $pos1 + $len1 - 1;
		$end1 = $pos1 if ($len1 == 0);
		my @pos2;
		foreach my $pos2 (sort {$a <=> $b} keys %{$ins_bp_cons{$chr}}){
			last if ($pos2 > $pos1 + $var_sd);
			next if ($pos2 <= $pos1);
			next if (exists $delete{$pos2});
			my $len2 = ${$ins_bp_len{$chr}}{$pos2};
			my $end2 = $pos2 + $len2 - 1;
			$end2 = $pos2 if ($len2 == 0);
			if (($pos2 - $pos1 <= $var_sd) and (abs ($end1 - $end2) <= $var_sd)){
				push @pos2, $pos2;
			}
		}
		if (@pos2 > 0){
			my $top_num = 0;
			my $top_pos = scalar @{${$ins_bp_cons{$chr}}{$pos1}};
			foreach my $pos2 (@pos2){
				my $num2 = scalar @{${$ins_bp_cons{$chr}}{$pos2}};
				if ($num2 > $top_num){
					$top_num = $num2;
					$top_pos = $pos2;
				}
			}
			push @pos2, $pos1;
			my @info;
			foreach my $pos2 (@pos2){
				next if ($pos2 == $top_pos);
				push @info,  @{${$ins_bp_cons{$chr}}{$pos2}};
				delete ${$ins_bp_cons{$chr}}{$pos2};
				$delete{$pos2} = 1;
			}
			push @{${$ins_bp_cons{$chr}}{$top_pos}}, @info;
		}
	}
}

foreach my $chr (keys %ins_bp_cons){
	foreach my $pos (sort {$a <=> $b} keys %{$ins_bp_cons{$chr}}){
		my %ids;
		my @info;
		foreach (@{${$ins_bp_cons{$chr}}{$pos}}){
			my ($id) = split (/==/, $_);
			if (!exists $ids{$id}){
				push @info, $_;
			}
			$ids{$id} = 1;
		}
		@{${$ins_bp_cons{$chr}}{$pos}} = (@info);
	}
}

my %ins_bplen;
my %ins_bplen2;
foreach my $chr (keys %{$vcf_type{'INS'}}){			# Merge > 200 bp INSs within INS with > 100 bp distance of breakpoints
	foreach my $pos1 (sort {$a <=> $b} keys %{${$vcf_type{'INS'}}{$chr}}){
		next if (!exists ${${$vcf_type{'INS'}}{$chr}}{$pos1});
		my $line1 = ${${$vcf_type{'INS'}}{$chr}}{$pos1};
		my $len1 = $1 if ($line1 =~ /SVLEN2=(\d+)/);
		next if ($len1 < 200);
		my $idpos1 = $1 if ($line1 =~ /SAMPLES=(.+)/);
		my @idpos1 = split (/,,/, $idpos1);
		my @bplen;
		foreach (@idpos1){
			my $BP = 0;
			my $bplen = 0;
			$BP = $1 if ($_ =~ /BP=(\d+)/);
			next if ($BP < 2);
			$bplen = $1 if ($_ =~ /BPLEN=(\d+)/);
			next if ($bplen < 100);
			push @bplen, $bplen;
		}
		next if (@bplen == 0);
		my $max_bplen = 0;
		my $min_bplen = 1000000000;
		foreach (@bplen){
			next if ($_ * 0.5 > $len1);
			if ($_ > $max_bplen){
				$max_bplen = $_;
			}
			if ($_ < $min_bplen){
				$min_bplen = $_;
			}
		}
		next if ($max_bplen == 0);
		${$ins_bplen{$chr}}{$pos1} = "$min_bplen=$max_bplen=$len1";
		push @{${$ins_bplen2{$chr}}{$pos1}}, @bplen;
		my @hit_pos;
		my @hit_len;
		foreach my $pos2 (sort {$a <=> $b} keys  %{${$vcf_type{'INS'}}{$chr}}){
			next if ($pos2 <= $pos1);
			last if ($pos2 > $pos1 + $max_bplen + 20);
			my $ilen = $1 if (${${$vcf_type{'INS'}}{$chr}}{$pos2} =~ /SVLEN2=(\d+)/);
			if (($ilen >= $len1 * 0.8) and ($ilen <= $len1 * 1.25)){
				push @hit_pos, $pos2;
				push @hit_len, $ilen;
			}
		}
		next if (@hit_pos == 0);
		my $info1 = $idpos1;
		my $ids1 = '';
		my $ids2 = '';
		my @info1 = split (/,,/, $info1);
		foreach (@info1){
			my ($id) = split (/==/, $_);
			$ids1 .= "$id,";
		}
		foreach my $ipos (@hit_pos){
			my $inf = $1 if (${${$vcf_type{'INS'}}{$chr}}{$ipos} =~ /SAMPLES=(.+)/);
			$info1 .= ",,$inf";
			delete ${${$vcf_type{'INS'}}{$chr}}{$ipos};
			my @inf = split (/,,/, $inf);
			foreach (@inf){
				my ($id) = split (/==/, $_);
				$ids2 .= "$id,";
			}
		}
		$line1 =~ s/SAMPLES=.+/SAMPLES=$info1/;
		${${$vcf_type{'INS'}}{$chr}}{$pos1} = $line1;
	}
}

my %ins_bp_assign;
foreach my $chr (keys %ins_bp_cons){			# Assign INS-BP to INS matching the breakpoints
	my $chr02d = $chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	foreach my $pos (sort {$a <=> $b} keys %{$ins_bp_cons{$chr}}){
		my $bplen = ${$ins_bp_len{$chr}}{$pos};
		my $end = $pos + $bplen - 1;
		my %hit_pos;
		my @hit_len;
		foreach my $ipos (sort {$a <=> $b} keys %{$ins_bplen{$chr02d}}){
			next if (!exists ${${$vcf_type{'INS'}}{$chr02d}}{$ipos});
			last if ($ipos > $pos + 30);
			my ($min_bplen, $max_bplen, $ilen) = split (/=/, ${$ins_bplen{$chr02d}}{$ipos});
			next if ($ipos < $pos - 30);
			if (($ipos + $min_bplen >= $end - 30) and ($ipos + $max_bplen <= $end + 30) and ($ilen >= 300) and ($ilen > $bplen * 0.5)){
				$hit_pos{$ipos} = $ilen;
				push @hit_len, $ilen;
			}
			elsif (($bplen < 100) and ($ilen >= 300)){
				$hit_pos{$ipos} = $ilen;
				push @hit_len, $ilen;
			}
		}
		next if (@hit_len == 0);
		@hit_len = sort {$a <=> $b} @hit_len;
		my $hn = int (@hit_len * 0.5);
		my $med_len = $hit_len[$hn];
		my $str = '';
		my $ids1 = '';
		my $ids2 = '';
		my $count = 0;
		foreach (@{${$ins_bp_cons{$chr}}{$pos}}){
			my ($id, $bpos, $blen, $balt, $bline) = split (/==/, $_);
			$ids1 .= "$id,";
			${${$ins_bp_cons{$chr}}{$pos}}[$count] = "$id==$bpos==$med_len==$balt==$bline";
			$count ++;
		}
		my $info = join (',,', @{${$ins_bp_cons{$chr}}{$pos}});
		foreach my $ipos (keys %hit_pos){
			my $ilen = $hit_pos{$ipos};
			if (($ilen / $med_len < 2) and ($ilen / $med_len > 0.5)){
				$str .= "$ipos-$ilen ";
				my $inf = $1 if (${${$vcf_type{'INS'}}{$chr02d}}{$ipos} =~ /SAMPLES=(.+)/);
				$info .= ",,$inf";
				delete ${${$vcf_type{'INS'}}{$chr02d}}{$ipos};
				my @inf = split (/,,/, $inf);
				foreach (@inf){
					my ($id) = split (/==/, $_);
					$ids2 .= "$id,";
				}
			}
		}
		${${$vcf_type{'INS'}}{$chr02d}}{$pos} = "SVLEN2=$med_len;SAMPLES=$info";
		${$ins_bplen{$chr02d}}{$pos} = "$bplen=$bplen=$med_len";
		push @{${$ins_bplen2{$chr02d}}{$pos}}, $bplen;
		${$ins_bp_assign{$chr02d}}{$pos} = "$bplen=$med_len";
		${$ins_bp_convert{$chr}}{$pos} = $bplen;
		delete ${$ins_bp_cons{$chr}}{$pos};
		delete ${$ins_bp_len{$chr}}{$pos};
	}
}

foreach my $chr (keys %ins_bp_assign){
	foreach my $pos1 (sort {$a <=> $b} keys %{$ins_bp_assign{$chr}}){
		my ($bplen, $len1) = split (/=/, ${$ins_bp_assign{$chr}}{$pos1});
		next if (!exists ${${$vcf_type{'INS'}}{$chr}}{$pos1});
		my $line1 = ${${$vcf_type{'INS'}}{$chr}}{$pos1};
		my @hit_pos;
		foreach my $pos2 (sort {$a <=> $b} keys %{${$vcf_type{'INS'}}{$chr}}){
			last if ($pos2 > $pos1 + $bplen + 20);
			next if ($pos2 < $pos1 - 20);
			next if ($pos2 == $pos1);
			my $ilen = $1 if (${${$vcf_type{'INS'}}{$chr}}{$pos2} =~ /SVLEN2=(\d+)/);
			if (($ilen >= $len1 * 0.8) and ($ilen <= $len1 * 1.25)){
				push @hit_pos, $pos2;
			}
		}
		next if (@hit_pos == 0);
		my $info1 = $1 if ($line1 =~ /SAMPLES=(.+)/);
		foreach my $ipos (@hit_pos){
			my $inf = $1 if (${${$vcf_type{'INS'}}{$chr}}{$ipos} =~ /SAMPLES=(.+)/);
			$info1 .= ",,$inf";
			delete ${${$vcf_type{'INS'}}{$chr}}{$ipos};
		}
		$line1 =~ s/SAMPLES=.+/SAMPLES=$info1/;
		${${$vcf_type{'INS'}}{$chr}}{$pos1} = $line1;
	}
}

foreach my $chr (keys %{$vcf_type{'DUP'}}){			# Assign INS-BP to DUP matching the breakpoints
	my $chr2 = $chr;
	$chr2 =~ s/^0*//;
	foreach my $pos (sort {$a <=> $b} keys %{${$vcf_type{'DUP'}}{$chr}}){
		my $dline = ${${$vcf_type{'DUP'}}{$chr} }{$pos};
		my $dlen = $1 if ($dline =~ /SVLEN2=(\d+)/);
		next if ($dlen > 20000);
		my $end = $pos + $dlen - 1;
		my @hit_pos;
		my @hit_len;
		foreach my $bpos (sort {$a <=> $b} keys %{$ins_bp_cons{$chr2}}){
			next if (!exists ${$ins_bp_cons{$chr2}}{$bpos}) or (!exists ${$ins_bp_len{$chr2}}{$bpos});
			my $bplen = ${$ins_bp_len{$chr2}}{$bpos};
			my $bend = $bpos + $bplen - 1;
			last if ($bpos > $pos + 50);
			next if ($bpos < $pos - 50);
			if (($bend >= $end - 50) and ($bend <= $end + 50)){
				push @hit_pos, $bpos;
				push @hit_len, $bplen;
			}
		}
		next if (@hit_pos == 0);
		my $info = $1 if ($dline =~ /SAMPLES=(.+)/);
		my $ids1 = '';
		my $ids2 = '';
		my @info = split (/,,/, $info);
		foreach (@info){
			my ($id) = split (/==/, $_);
			$ids1 .= "$id,";
		}
		foreach my $ipos (@hit_pos){
			my $count = 0;
			foreach (@{${$ins_bp_cons{$chr2}}{$ipos}}){
				my ($id, $bpos, $blen, $balt, $bline) = split (/==/, $_);
				${${$ins_bp_cons{$chr}}{$ipos}}[$count] = "$id==$bpos==$dlen==$balt==$bline";
				$ids2 .= "$id,";
				$count ++;
			}
			my $inf = join (',,', @{${$ins_bp_cons{$chr2}}{$ipos}});
			$info .= ",,$inf";
			delete ${$ins_bp_cons{$chr2}}{$ipos};
			delete ${$ins_bp_len{$chr2}}{$ipos};
		}
		$dline =~ s/SAMPLES=.+/SAMPLES=$info/;
		${${$vcf_type{'DUP'}}{$chr}}{$pos} = $dline;
	}
}

foreach my $chr (keys %{$vcf_type{'DUP'}}){			# Merge INSs within < 20 Kb DUP
	my $chr2 = $chr;
	$chr2 =~ s/^0*//;
	foreach my $pos (sort {$a <=> $b} keys %{${$vcf_type{'DUP'}}{$chr}}){
		my $dline = ${${$vcf_type{'DUP'}}{$chr} }{$pos};
		my $dlen = $1 if ($dline =~ /SVLEN2=(\d+)/);
		next if ($dlen > 20000);
		my $end = $pos + $dlen - 1;
		my @hit_pos;
		my @hit_len;
		foreach my $ipos (sort {$a <=> $b} keys %{${$vcf_type{'INS'}}{$chr}}){
			last if ($ipos > $end + 50);
			next if ($ipos < $pos - 50);
			next if (exists ${$ins_bp_convert{$chr2}}{$ipos});
			my $ilen = $1 if (${${$vcf_type{'INS'}}{$chr}}{$ipos} =~ /SVLEN2=(\d+)/);
			next if ($ilen < $dlen * 0.5);
			if (exists ${$ins_bplen{$chr}}{$ipos}){
				my ($min_bplen, $max_bplen, $ilen2) = split (/=/, ${$ins_bplen{$chr}}{$ipos});
				next if ($dlen < $min_bplen * 0.8) and ($dlen > $max_bplen * 1.2);
			}
			push @hit_pos, $ipos;
			push @hit_len, $ilen;
		}
		next if (@hit_pos == 0);
		my $info = $1 if ($dline =~ /SAMPLES=(.+)/);
		my $ids1 = '';
		my $ids2 = '';
		my @info = split (/,,/, $info);
		foreach (@info){
			my ($id) = split (/==/, $_);
			$ids1 .= "$id,";
		}
		foreach my $ipos (@hit_pos){
			my $inf = $1 if (${${$vcf_type{'INS'}}{$chr}}{$ipos} =~ /SAMPLES=(.+)/);
			$info .= ",,$inf";
			delete ${${$vcf_type{'INS'}}{$chr}}{$ipos};
			my @inf = split (/,,/, $inf);
			foreach (@inf){
				my ($id) = split (/==/, $_);
				$ids2 .= "$id,";
			}
		}
		$dline =~ s/SAMPLES=.+/SAMPLES=$info/;
		${${$vcf_type{'DUP'}}{$chr}}{$pos} = $dline;
	}
}

print STDERR "3rd step completed:\n";

foreach my $type (keys %vcf_type){              # adjust sample variants between neighboring positions
	next if ($type eq 'INS');
    foreach my $chr (keys %{$vcf_type{$type}}){
        my $pre_pos = 0;
        my $pre_end = 0;
        my $pre_len = 0;
        my $pre_line = '';
        my $chr2 = $chr;
        $chr2 =~ s/^0*//;
        foreach my $pos (sort {$a <=> $b} keys %{${$vcf_type{$type}}{$chr}}){
            my $line = ${${$vcf_type{$type}}{$chr}}{$pos};
            my $len = $1 if ($line =~ /SVLEN2=-*(\d+)/);
            my $end = $pos + $len - 1;
            if ($pre_pos == 0){
            	$pre_pos = $pos;
	            $pre_end = $end;
	            $pre_len = $len;
	            $pre_line = $line;
	            next;
            }
            if (($pos < $pre_end) or (($len <= 5) and ($pos < $pre_end + $len))){
            	my $ovl_flag = 0;
            	my $overlap = $pre_end - $pos + 1;
				$overlap = $len if ($end < $pre_end);
            	if (($overlap >= $pre_len * $min_overlap_ratio) and ($overlap >= $len * $min_overlap_ratio)){
            		$ovl_flag = 1;
            	}
            	elsif (($len <= 10) and ($overlap > 0) and ($len / $pre_len <= 2) and ($len / $pre_len >= 0.5)){
            		$ovl_flag = 1;
            	}
            	elsif (($len <= 5) and (abs ($pos - $pre_end) <= 5) and ($len / $pre_len <= 2) and ($len / $pre_len >= 0.5)){
            		$ovl_flag = 1;
            	}
				if ($ovl_flag == 1){
				    my $idpos = $1 if ($line =~ /SAMPLES=(.+)/);
				    my $pre_idpos = $1 if ($pre_line =~ /SAMPLES=(.+)/);
				    my @idpos = split (/,,/, $idpos);
				    my @pre_idpos = split (/,,/, $pre_idpos);
				    my $sn = scalar @idpos;
				    my $pre_sn = scalar @pre_idpos;
				    my %add_idpos;
				    my %add_preidpos;
				    my @add_preidpos;
				    my @add_idpos;
				    foreach (@pre_idpos){
						my ($gid, $pos2, $len2) = split (/==/, $_);
						my $diff_rate1 = abs ($pos2 - $pre_pos) / $pre_pos;
						$diff_rate1 = $diff_rate1 * (abs ($len2 - $pre_len) / $pre_len);
						my $diff_rate2 = abs ($pos2 - $pos) / $pos;
						$diff_rate2 = $diff_rate2 * (abs ($len2 - $len) / $len);
						if ($diff_rate1 > $diff_rate2){
						    $add_idpos{$_} = 1;
						    push @add_idpos, "$_|$diff_rate1==$diff_rate2";
						}
				    }
				    foreach (@idpos){
						my ($gid, $pos2, $len2) = split (/==/, $_);
						my $diff_rate1 = abs ($pos2 - $pre_pos) / $pre_pos;
						$diff_rate1 = $diff_rate1 * (abs ($len2 - $pre_len) / $pre_len);
						my $diff_rate2 = abs ($pos2 - $pos) / $pos;
						$diff_rate2 = $diff_rate2 * (abs ($len2 - $len) / $len);
						if ($diff_rate1 < $diff_rate2){
						    $add_preidpos{$_} = 1;
						    push @add_preidpos, "$_|$diff_rate1==$diff_rate2";
						}
				    }
				    foreach (@pre_idpos){
						$add_preidpos{$_} = 1 if (!exists $add_idpos{$_});
				    }
				    foreach (@idpos){
						$add_idpos{$_} = 1 if (!exists $add_preidpos{$_});
				    }
				    if (@pre_idpos != scalar keys %add_preidpos){
						my $new_sn = 0;
						my $new_presn = 0;
						my $new_pos = 0;
						my $new_prepos = 0;
						my $new_len = 0;
						my $new_prelen = 0;
						my $sum_pos = 0;
						my $sum_prepos = 0;
						my $sum_len = 0;
						my $sum_prelen = 0;
						my $new_idpos = '';
						my $new_preidpos = '';
						$new_presn = scalar keys %add_preidpos;
						$new_sn = scalar keys %add_idpos;
						foreach my $idpos2 (sort keys %add_preidpos){
						    my ($gid, $pos2, $len2) = split (/==/, $idpos2);
						    $sum_prepos += $pos2;
						    $sum_prelen += $len2;
						    $new_preidpos .= "$idpos2,,"
						}
						$new_prepos = int ($sum_prepos / $new_presn + 0.5) if ($new_presn > 0);
						$new_prelen = int ($sum_prelen / $new_presn + 0.5) if ($new_presn > 0);
						$new_preidpos =~ s/,,$//;
						foreach my $idpos2 (sort keys %add_idpos){
						    my ($gid, $pos2, $len2) = split (/==/, $idpos2);
						    $sum_pos += $pos2;
						    $sum_len += $len2;
						    $new_idpos .= "$idpos2,,"
						}
						$new_pos = int ($sum_pos / $new_sn + 0.5) if ($new_sn > 0);
						$new_len = int ($sum_len / $new_sn + 0.5) if ($new_sn > 0);
						$new_idpos =~ s/,,$//;
						delete ${${$vcf_type{$type}}{$chr}}{$pre_pos};
						delete ${${$vcf_type{$type}}{$chr}}{$pos};
						${${$vcf_type{$type}}{$chr}}{$new_prepos} = "SVLEN2=$new_prelen;SAMPLES=$new_preidpos" if ($new_presn > 0);
						${${$vcf_type{$type}}{$chr}}{$new_pos} = "SVLEN2=$new_len;SAMPLES=$new_idpos" if ($new_sn > 0);
						if ($new_sn > 0){
						    $pre_pos = $new_pos;
						    $pre_end = $new_pos + $new_len - 1;
						    $pre_len = $new_len;
						    $pre_line = ${${$vcf_type{$type}}{$chr}}{$new_pos};
						}
						else{
						    $pre_pos = $new_prepos;
						    $pre_end = $new_prepos + $new_prelen - 1;
						    $pre_len = $new_prelen;
						    $pre_line = ${${$vcf_type{$type}}{$chr}}{$new_prepos};
						}
						next;
				    }
				}
            }
            $pre_pos = $pos;
            $pre_end = $end;
            $pre_len = $len;
            $pre_line = $line;
        }
    }
}

if (($ave_call > 100000) or ($total_sample_num >= 1000)){     # merge overlapping variants
	my %vcf_type2;
	my $Mbin_size = 1000000;
	foreach my $type (keys %vcf_type){
	    foreach my $chr (keys %{$vcf_type{$type}}){
	    	foreach my $pos (sort {$a <=> $b} keys %{${$vcf_type{$type}}{$chr}}){
	    		my $Mbin1 = int ($pos / $Mbin_size);
	    		my $line = ${${$vcf_type{$type}}{$chr}}{$pos};
	            my $len = $1 if ($line =~ /SVLEN2=-*(\d+)/);
	            my $end = $pos + $len - 1;
	            my $Mbin2 = int ($end / $Mbin_size);
	            ${${${$vcf_type2{$type}}{$chr}}{$Mbin1}}{$pos} = $line;
	            if ($Mbin2 > $Mbin1){
	            	${${${$vcf_type2{$type}}{$chr}}{$Mbin2}}{$pos} = $line;
	            }
	        }
	    }
	}
	%vcf_type = ();
	foreach my $type (keys %vcf_type2){                  # merge overlapping variants
	    foreach my $chr (keys %{$vcf_type2{$type}}){
	        my $chr2 = $chr;
	        $chr2 =~ s/^0*//;
	        my %removed;
	        foreach my $bin1 (sort {$a <=> $b} keys %{${$vcf_type2{$type}}{$chr}}){
		        foreach my $pos1 (sort {$a <=> $b} keys %{${${$vcf_type2{$type}}{$chr}}{$bin1}}){
		        	next if (exists $removed{$pos1});
		            my $line1 = ${${${$vcf_type2{$type}}{$chr}}{$bin1}}{$pos1};
		            my $len1 = $1 if ($line1 =~ /SVLEN2=-*(\d+)/);
		            my $end1 = $pos1 + $len1 - 1;
		            $end1 = $pos1 if ($type eq 'INS');
		            my $idpos1 = $1 if ($line1 =~ /SAMPLES=(.+)/);
		            my @idpos1 = split (/,,/, $idpos1);
		            my $sn1 = scalar @idpos1;
		            foreach my $pos2 (sort {$a <=> $b} keys %{${${$vcf_type2{$type}}{$chr}}{$bin1}}){
		            	next if ($pos2 <= $pos1);
		            	next if (exists $removed{$pos2});
		            	last if ($pos2 > $end1 + 1000);
			            my $line2 = ${${${$vcf_type2{$type}}{$chr}}{$bin1}}{$pos2};
			            my $len2 = $1 if ($line2 =~ /SVLEN2=-*(\d+)/);
			            next if ($len2 == 0);
			            my $lenrate = int ($len1 / $len2 * 100 + 0.5) / 100;
			            next if ($lenrate > 2) or ($lenrate < 0.5);
			            my $end2 = $pos2 + $len2 - 1;
			            $end2 = $pos2 if ($type eq 'INS');
			            my $idpos2 = $1 if ($line2 =~ /SAMPLES=(.+)/);
			            my @idpos2 = split (/,,/, $idpos2);
		            	my $sn2 = scalar @idpos2;
			            my $distance = $pos2 - $end1;
			            if ($type eq 'INS'){
			            	my %ids;
		                    my $new_sn = 0;
		                    my $new_idpos = '';
		                    my $new_pos = 0;
		                    my $new_len = 0;
		                    my @len;
		                    my $flag = 0;
		                    my $ovl_flag = 0;
		                    if (($distance <= 150) and ($len1 * 0.5 >= $distance) and ($len2 * 0.5 >= $distance) and ($lenrate > 0.556) and ($lenrate < 1.8)){
		                    	$ovl_flag = 1;
		                    }
		                    elsif (($distance <= 300) and ($len1 >= $distance) and ($len2 >= $distance) and ($lenrate >= 0.8) and ($lenrate <= 1.25)){
		                    	$ovl_flag = 1;
		                    }
		                    elsif (($distance <= 1000) and ($len1 >= $distance) and ($len2 >= $distance) and ($lenrate >= 0.95) and ($lenrate <= 1.05)){
		                    	$ovl_flag = 1;
		                    }
		                    elsif (($distance <= 10) and ($len1 <= 10) and ($len1 > 3) and ($len2 <= 10) and ($len2 > 3) and ($lenrate >= 0.5) and ($lenrate <= 2)){
		                    	$ovl_flag = 1;
		                    }
		                    elsif (($distance <= 2) and ($len1 <= 3) and ($len2 <= 3) and ($lenrate >= 0.5) and ($lenrate <= 2)){
		                    	$ovl_flag = 1;
		                    }
		                    if ($ovl_flag == 1){
		                        if ($sn2 >= $sn1){
		                            foreach (@idpos1){
			                            my ($id, $spos) = split (/==/, $_);
			                            ${$ids{$id}}{$spos} = $_;
			                        }
			                        foreach (@idpos2){
			                            my ($id, $spos) = split (/==/, $_);
			                            ${$ids{$id}}{$spos} = $_;
			                        }
			                        $new_pos = $pos2;
			                        delete ${${${$vcf_type2{$type}}{$chr}}{$bin1}}{$pos1};
			                        $removed{$pos1} = 1;
			                        $flag = 2;
		                        }
		                        else{
		                            foreach (@idpos2){
			                            my ($id, $spos) = split (/==/, $_);
			                            ${$ids{$id}}{$spos} = $_;
			                        }
			                        foreach (@idpos1){
			                            my ($id, $spos) = split (/==/, $_);
			                            ${$ids{$id}}{$spos} = $_;
			                        }
			                        $new_pos = $pos1;
			                        delete ${${${$vcf_type2{$type}}{$chr}}{$bin1}}{$pos2};
			                        $removed{$pos2} = 1;
			                        $flag = 1;
		                        }
		                    }
		                    if ($flag >= 1){
			                    foreach my $gid (sort keys %ids){
			                    	foreach my $spos (sort {$a <=> $b} keys %{$ids{$gid}}){
				                        my ($gid2, $spos2, $len2) = split (/==/, ${$ids{$gid}}{$spos});
				                        $new_idpos .= "${$ids{$gid}}{$spos},,";
				                        push @len, $len2;
				                    }
			                    }
			                    my $hnum = int (@len * 0.5);
			                    $new_len = $len[$hnum];
			                    $new_idpos =~ s/,,$//;
			                    ${${${$vcf_type2{$type}}{$chr}}{$bin1}}{$new_pos} = "SVLEN2=$new_len;SAMPLES=$new_idpos";
			                    delete $removed{$new_pos} if (exists $removed{$new_pos});
			                    last if ($flag == 2);
		                    }
			            }
			            else{
			            	if ($distance <= 5){
		                        my $overlap = $end1 - $pos2 + 1;
		                        $overlap = $len2 if ($end2 < $end1);
		                        my $ovl_flag = 0;
		                        if (($overlap >= $len1 * $min_overlap_ratio2) and ($overlap >= $len2 * $min_overlap_ratio2)){
		                        	$ovl_flag = 1;
		                        }
		                        elsif (($overlap > 0) and ($len1 <= 10) and ($len2 <= 10) and ($lenrate <= 2) and ($lenrate >= 0.5)){
		                        	$ovl_flag = 1;
		                        }
		                        elsif (($len1 <= 5) and ($len2 <= 5) and ($lenrate <= 2) and ($lenrate >= 0.5)){
		                        	$ovl_flag = 1;
		                        }
		                        if ($ovl_flag == 1){
		                        	my %ids;
				                    my $new_sn = 0;
				                    my $new_idpos = '';
				                    my $new_pos = 0;
				                    my $new_len = 0;
				                    my @len;
				                    my $flag = 0;
				                    if ($sn2 >= $sn1){
				                        foreach (@idpos1){
				                            my ($id, $spos) = split (/==/, $_);
				                            ${$ids{$id}}{$spos} = $_;
				                        }
				                        foreach (@idpos2){
				                            my ($id, $spos) = split (/==/, $_);
				                            ${$ids{$id}}{$spos} = $_;
				                        }
				                        $new_pos = $pos2;
				                        delete ${${${$vcf_type2{$type}}{$chr}}{$bin1}}{$pos1};
				                        $removed{$pos1} = 1;
				                        $flag = 2;
				                    }
				                    else{
				                        foreach (@idpos2){
				                            my ($id, $spos) = split (/==/, $_);
				                            ${$ids{$id}}{$spos} = $_;
				                        }
				                        foreach (@idpos1){
				                            my ($id, $spos) = split (/==/, $_);
				                            ${$ids{$id}}{$spos} = $_;
				                        }
				                        $new_pos = $pos1;
				                        delete ${${${$vcf_type2{$type}}{$chr}}{$bin1}}{$pos2};
				                        $removed{$pos2} = 1;
				                        $flag = 1;
				                    }
				                    foreach my $gid (sort keys %ids){
				                    	foreach my $spos (sort {$a <=> $b} keys %{$ids{$gid}}){
					                        my ($gid2, $spos2, $len2) = split (/==/, ${$ids{$gid}}{$spos});
					                        $new_idpos .= "${$ids{$gid}}{$spos},,";
					                        push @len, $len2;
					                    }
				                    }
				                    my $hnum = int (@len * 0.5);
				                    $new_len = $len[$hnum];
				                    $new_idpos =~ s/,,$//;
				                    ${${${$vcf_type2{$type}}{$chr}}{$bin1}}{$new_pos} = "SVLEN2=$new_len;SAMPLES=$new_idpos";
				                    delete $removed{$new_pos} if (exists $removed{$new_pos});
				                    last if ($flag == 2);
		                        }
		                    }
			            }
			        }
			    }
			}
		}
	}
	foreach my $type (keys %vcf_type2){
	    foreach my $chr (keys %{$vcf_type2{$type}}){
	    	foreach my $bin (keys %{${$vcf_type2{$type}}{$chr}}){
		        foreach my $pos (keys %{${${$vcf_type2{$type}}{$chr}}{$bin}}){
		            ${${$vcf_cons{$chr}}{$pos}}{$type} = ${${${$vcf_type2{$type}}{$chr}}{$bin}}{$pos} if (!exists ${${$vcf_cons{$chr}}{$pos}}{$type});
		        }
		    }
	    }
	}
	%vcf_type2 = ();
}
else{
	foreach my $type (keys %vcf_type){
	    foreach my $chr (keys %{$vcf_type{$type}}){
	        my $chr2 = $chr;
	        $chr2 =~ s/^0*//;
	        my %removed;
	        foreach my $pos1 (sort {$a <=> $b} keys %{${$vcf_type{$type}}{$chr}}){
	        	next if (exists $removed{$pos1});
	            my $line1 = ${${$vcf_type{$type}}{$chr}}{$pos1};
	            my $len1 = $1 if ($line1 =~ /SVLEN2=-*(\d+)/);
	            my $end1 = $pos1 + $len1 - 1;
	            $end1 = $pos1 if ($type eq 'INS');
	            my $idpos1 = $1 if ($line1 =~ /SAMPLES=(.+)/);
	            my @idpos1 = split (/,,/, $idpos1);
	            my $sn1 = scalar @idpos1;
	            foreach my $pos2 (sort {$a <=> $b} keys %{${$vcf_type{$type}}{$chr}}){
	            	next if ($pos2 <= $pos1);
	            	next if (exists $removed{$pos2});
	            	last if ($pos2 > $end1 + 1000);
		            my $line2 = ${${$vcf_type{$type}}{$chr}}{$pos2};
		            my $len2 = $1 if ($line2 =~ /SVLEN2=-*(\d+)/);
		            next if ($len2 == 0);
		            my $lenrate = int ($len1 / $len2 * 100 + 0.5) / 100;
		            next if ($lenrate > 2) or ($lenrate < 0.5);
		            my $end2 = $pos2 + $len2 - 1;
		            $end2 = $pos2 if ($type eq 'INS');
		            my $idpos2 = $1 if ($line2 =~ /SAMPLES=(.+)/);
		            my @idpos2 = split (/,,/, $idpos2);
	            	my $sn2 = scalar @idpos2;
		            my $distance = $pos2 - $end1;
		            if ($type eq 'INS'){
		            	my %ids;
	                    my $new_sn = 0;
	                    my $new_idpos = '';
	                    my $new_pos = 0;
	                    my $new_len = 0;
	                    my @len;
	                    my $flag = 0;
	                    my $ovl_flag = 0;
	                    if (($distance <= 150) and ($len1 * 0.5 >= $distance) and ($len2 * 0.5 >= $distance) and ($lenrate > 0.556) and ($lenrate < 1.8)){
	                    	$ovl_flag = 1;
	                    }
	                    elsif (($distance <= 300) and ($len1 >= $distance) and ($len2 >= $distance) and ($lenrate >= 0.8) and ($lenrate <= 1.25)){
	                    	$ovl_flag = 1;
	                    }
	                    elsif (($distance <= 1000) and ($len1 >= $distance) and ($len2 >= $distance) and ($lenrate >= 0.95) and ($lenrate <= 1.05)){
	                    	$ovl_flag = 1;
	                    }
	                    elsif (($distance <= 10) and ($len1 <= 10) and ($len1 > 3) and ($len2 <= 10) and ($len2 > 3) and ($lenrate >= 0.5) and ($lenrate <= 2)){
	                    	$ovl_flag = 1;
	                    }
	                    elsif (($distance <= 2) and ($len1 <= 3) and ($len2 <= 3) and ($lenrate >= 0.5) and ($lenrate <= 2)){
	                    	$ovl_flag = 1;
	                    }
	                    if ($ovl_flag == 1){
	                        if ($sn2 >= $sn1){
	                            foreach (@idpos1){
		                            my ($id, $spos) = split (/==/, $_);
		                            ${$ids{$id}}{$spos} = $_;
		                        }
		                        foreach (@idpos2){
		                            my ($id, $spos) = split (/==/, $_);
		                            ${$ids{$id}}{$spos} = $_;
		                        }
		                        $new_pos = $pos2;
		                        delete ${${$vcf_type{$type}}{$chr}}{$pos1};
		                        $removed{$pos1} = 1;
		                        $flag = 2;
	                        }
	                        else{
	                            foreach (@idpos2){
		                            my ($id, $spos) = split (/==/, $_);
		                            ${$ids{$id}}{$spos} = $_;
		                        }
		                        foreach (@idpos1){
		                            my ($id, $spos) = split (/==/, $_);
		                            ${$ids{$id}}{$spos} = $_;
		                        }
		                        $new_pos = $pos1;
		                        delete ${${$vcf_type{$type}}{$chr}}{$pos2};
		                        $removed{$pos2} = 1;
		                        $flag = 1;
	                        }
	                    }
	                    if ($flag >= 1){
		                    foreach my $gid (sort keys %ids){
		                    	foreach my $spos (sort {$a <=> $b} keys %{$ids{$gid}}){
			                        my ($gid2, $spos2, $len2) = split (/==/, ${$ids{$gid}}{$spos});
			                        $new_idpos .= "${$ids{$gid}}{$spos},,";
			                        push @len, $len2;
			                    }
		                    }
		                    my $hnum = int (@len * 0.5);
		                    $new_len = $len[$hnum];
		                    $new_idpos =~ s/,,$//;
		                    ${${$vcf_type{$type}}{$chr}}{$new_pos} = "SVLEN2=$new_len;SAMPLES=$new_idpos";
		                    delete $removed{$new_pos} if (exists $removed{$new_pos});
		                    last if ($flag == 2);
	                    }
		            }
		            else{
		            	if ($distance <= 5){
	                        my $overlap = $end1 - $pos2 + 1;
	                        $overlap = $len2 if ($end2 < $end1);
	                        my $ovl_flag = 0;
	                        if (($overlap >= $len1 * $min_overlap_ratio2) and ($overlap >= $len2 * $min_overlap_ratio2)){
	                        	$ovl_flag = 1;
	                        }
	                        elsif (($overlap > 0) and ($len1 <= 10) and ($len2 <= 10) and ($lenrate <= 2) and ($lenrate >= 0.5)){
	                        	$ovl_flag = 1;
	                        }
	                        elsif (($len1 <= 5) and ($len2 <= 5) and ($lenrate <= 2) and ($lenrate >= 0.5)){
	                        	$ovl_flag = 1;
	                        }
	                        if ($ovl_flag == 1){
	                        	my %ids;
			                    my $new_sn = 0;
			                    my $new_idpos = '';
			                    my $new_pos = 0;
			                    my $new_len = 0;
			                    my @len;
			                    my $flag = 0;
			                    if ($sn2 >= $sn1){
			                        foreach (@idpos1){
			                            my ($id, $spos) = split (/==/, $_);
			                            ${$ids{$id}}{$spos} = $_;
			                        }
			                        foreach (@idpos2){
			                            my ($id, $spos) = split (/==/, $_);
			                            ${$ids{$id}}{$spos} = $_;
			                        }
			                        $new_pos = $pos2;
			                        delete ${${$vcf_type{$type}}{$chr}}{$pos1};
			                        $removed{$pos1} = 1;
			                        $flag = 2;
			                    }
			                    else{
			                        foreach (@idpos2){
			                            my ($id, $spos) = split (/==/, $_);
			                            ${$ids{$id}}{$spos} = $_;
			                        }
			                        foreach (@idpos1){
			                            my ($id, $spos) = split (/==/, $_);
			                            ${$ids{$id}}{$spos} = $_;
			                        }
			                        $new_pos = $pos1;
			                        delete ${${$vcf_type{$type}}{$chr}}{$pos2};
			                        $removed{$pos2} = 1;
			                        $flag = 1;
			                    }
			                    foreach my $gid (sort keys %ids){
			                    	foreach my $spos (sort {$a <=> $b} keys %{$ids{$gid}}){
				                        my ($gid2, $spos2, $len2) = split (/==/, ${$ids{$gid}}{$spos});
				                        $new_idpos .= "${$ids{$gid}}{$spos},,";
				                        push @len, $len2;
				                    }
			                    }
			                    my $hnum = int (@len * 0.5);
			                    $new_len = $len[$hnum];
			                    $new_idpos =~ s/,,$//;
			                    ${${$vcf_type{$type}}{$chr}}{$new_pos} = "SVLEN2=$new_len;SAMPLES=$new_idpos";
			                    delete $removed{$new_pos} if (exists $removed{$new_pos});
			                    last if ($flag == 2);
	                        }
	                    }
		            }
		        }
		    }
		}
	}
	foreach my $type (keys %vcf_type){
	    foreach my $chr (keys %{$vcf_type{$type}}){
	        foreach my $pos (keys %{${$vcf_type{$type}}{$chr}}){
	            ${${$vcf_cons{$chr}}{$pos}}{$type} = ${${$vcf_type{$type}}{$chr}}{$pos};
	        }
	    }
	}
	%vcf_type = ();
}

foreach my $chr (sort keys %vcf_cons){		# divide multi-allelic INSs
	my $chr2 = $chr;
	$chr2 =~ s/^0*//;
	foreach my $pos (sort {$a <=> $b} keys %{$vcf_cons{$chr}}){
		foreach my $type (keys %{${$vcf_cons{$chr}}{$pos}}){
			next if ($type ne 'INS');
			my $len = $1 if (${${$vcf_cons{$chr}}{$pos}}{$type} =~ /SVLEN2=-*(\d+)/);
			my $id_info = $1 if (${${$vcf_cons{$chr}}{$pos}}{$type} =~ /SAMPLES=(\S+)/);
			my @id_info = split (/,,/, $id_info);
			my $sample_num = scalar @id_info;
			next if ($sample_num == 1);
			my @len;
			my %id_pos_len;
			foreach (@id_info){
				my ($id, $pos2, $len2) = split (/==/, $_);
				if (($_ !~ /==/) or (!defined $pos2)){
					next;
				}
				$len2 = $len if ($len2 == 0);
				push @len, $len2;
				${$id_pos_len{$id}}{$pos2} = $len2;
			}
			my $max_len = 0;
			my $min_len = 10000000000000;
			foreach (@len){
				if ($_ > $max_len){
					$max_len = $_;
				}
				if ($_ < $min_len){
					$min_len = $_;
				}
			}
			if (($max_len / $min_len > 1.8) and ($max_len > 3)){
				my %max_id;
				my %min_id;
				my %med_id;
				my @max_pos;
				my @min_pos;
				my @med_pos;
				my @min_len;
				my @med_len;
				my @max_len;
				foreach my $id (keys %id_pos_len){
					foreach my $pos2 (keys %{$id_pos_len{$id}}){
						my $ilen = ${$id_pos_len{$id}}{$pos2};
						my $max_diff = $max_len - $ilen;
				        my $min_diff = $ilen - $min_len;
				        if (($max_diff >= 20) and ($min_diff >= 20) and ($max_diff / $ilen > 0.5) and ($min_diff / $ilen > 0.5) and ($max_len >= 50)){
				        	$med_id{$id} = $pos2;
				            push @med_pos, $pos2;
				            push @med_len, $ilen;
				        }
				        elsif ($max_diff <= $min_diff){
				            $max_id{$id} = $pos2;
				            push @max_pos, $pos2;
				            push @max_len, $ilen;
				        }
				        elsif ($max_diff > $min_diff){
				            $min_id{$id} = $pos2;
				            push @min_pos, $pos2;
				            push @min_len, $ilen;
				        }
				    }
				}
				if ((@med_len > 0) and (@min_len > 0) and (@max_len > 0)){
					my $sum_minlen = 0;
					my $sum_medlen = 0;
					my $sum_maxlen = 0;
					map{$sum_minlen += $_} @min_len;
					map{$sum_medlen += $_} @med_len;
					map{$sum_maxlen += $_} @max_len;
					my $ave_minlen = int ($sum_minlen / @min_len + 0.5);
					my $ave_medlen = int ($sum_medlen / @med_len + 0.5);
					my $ave_maxlen = int ($sum_maxlen / @max_len + 0.5);
					my $max_diff2 = $ave_maxlen - $ave_medlen;
					my $med_diff2 = $ave_medlen - $ave_minlen;
					if (($max_diff2 < 20) or ($med_diff2 < 20)){
						%max_id = ();
						%min_id = ();
						%med_id = ();
						@max_pos = ();
						@min_pos = ();
						@med_pos = ();
						foreach my $id (keys %id_pos_len){
							foreach my $pos2 (keys %{$id_pos_len{$id}}){
								my $ilen = ${$id_pos_len{$id}}{$pos2};
								my $max_diff = $max_len - $ilen;
						        my $min_diff = $ilen - $min_len;
						        if ($max_diff <= $min_diff){
						            $max_id{$id} = $pos2;
						            push @max_pos, $pos2;
						        }
						        elsif ($max_diff > $min_diff){
						            $min_id{$id} = $pos2;
						            push @min_pos, $pos2;
						        }
						    }
						}
					}
				}
				my $max_info = '';
				my $min_info = '';
				my $med_info = '';
				my @maxid;
				my @minid;
				my @medid;
				foreach (@id_info){
					my ($id, $pos2) = split (/==/, $_);
					if ((exists $max_id{$id}) and ($max_id{$id} == $pos2)){
						$max_info .= "$_,,";
						push @maxid, $id;
					}
					elsif ((exists $min_id{$id}) and ($min_id{$id} == $pos2)){
						$min_info .= "$_,,";
						push @minid, $id;
					}
					elsif ((exists $med_id{$id}) and ($med_id{$id} == $pos2)){
						$med_info .= "$_,,";
						push @medid, $id;
					}
				}
				$max_info =~ s/,,$//;
				$min_info =~ s/,,$//;
				$med_info =~ s/,,$// if ($med_info =~ /,,$/);
				my $max_pos = 0;
				my $min_pos = 0;
				my $med_pos = 0;
				if (@max_pos == 1){
					$max_pos = $max_pos[0];
				}
				elsif (@max_pos == 2){
					$max_pos = int (($max_pos[0] + $max_pos[1]) * 0.5 + 0.5);
				}
				elsif (@max_pos > 2){
					my $hmax = int (@max_pos * 0.5);
					@max_pos = sort {$a <=> $b} @max_pos;
					$max_pos = $max_pos[$hmax];
				}
				if (@min_pos == 1){
					$min_pos = $min_pos[0];
				}
				elsif (@min_pos == 2){
					$min_pos = int (($min_pos[0] + $min_pos[1]) * 0.5 + 0.5);
				}
				elsif (@min_pos > 2){
					my $hmin = int (@min_pos * 0.5);
					@min_pos = sort {$a <=> $b} @min_pos;
					$min_pos = $min_pos[$hmin];
				}
				if (@med_pos == 1){
					$med_pos = $med_pos[0];
				}
				elsif (@med_pos == 2){
					$med_pos = int (($med_pos[0] + $med_pos[1]) * 0.5 + 0.5);
				}
				elsif (@med_pos > 2){
					my $hmin = int (@med_pos * 0.5);
					@med_pos = sort {$a <=> $b} @med_pos;
					$med_pos = $med_pos[$hmin];
				}
				if (($max_pos > 0) and ($min_pos > 0) and ($med_pos > 0)){
					if (($max_pos == $min_pos) and ($max_pos == $med_pos)){
						while (1){
							$med_pos --;
							last if (!exists ${${$vcf_cons{$chr}}{$med_pos}}{$type});
						}
						$min_pos = $med_pos;
						while (1){
							$min_pos --;
							last if (!exists ${${$vcf_cons{$chr}}{$min_pos}}{$type});
						}
					}
					elsif ($max_pos == $min_pos){
						while (1){
							$min_pos --;
							last if (!exists ${${$vcf_cons{$chr}}{$min_pos}}{$type});
						}
					}
					elsif ($max_pos == $med_pos){
						while (1){
							$med_pos --;
							last if (!exists ${${$vcf_cons{$chr}}{$med_pos}}{$type});
						}
					}
					elsif ($min_pos == $med_pos){
						while (1){
							$min_pos --;
							last if (!exists ${${$vcf_cons{$chr}}{$min_pos}}{$type});
						}
					}
					my $arg1 = ${${$vcf_cons{$chr}}{$pos}}{$type};
					my $arg2 = $arg1;
					my $arg3 = $arg1;
					$arg1 =~ s/SAMPLES=\S+/SAMPLES=$max_info/;
					$arg2 =~ s/SAMPLES=\S+/SAMPLES=$min_info/;
					$arg3 =~ s/SAMPLES=\S+/SAMPLES=$med_info/;
					delete ${${$vcf_cons{$chr}}{$pos}}{$type};
					delete ${$vcf_cons{$chr}}{$pos} if (scalar keys %{${$vcf_cons{$chr}}{$pos}} == 0);
					${${$vcf_cons{$chr}}{$max_pos}}{$type} = $arg1;
					${${$vcf_cons{$chr}}{$min_pos}}{$type} = $arg2;
					${${$vcf_cons{$chr}}{$med_pos}}{$type} = $arg3;
					if (exists ${$ins_bp_convert{$chr2}}{$pos}){
						my $bplen = ${$ins_bp_convert{$chr2}}{$pos};
						${$ins_bp_convert{$chr2}}{$max_pos} = $bplen;
						${$ins_bp_convert{$chr2}}{$min_pos} = $bplen;
						${$ins_bp_convert{$chr2}}{$med_pos} = $bplen;
					}
				}
				elsif (($max_pos > 0) and ($min_pos > 0)){
					if ($max_pos == $min_pos){
						while (1){
							$min_pos --;
							last if (!exists ${${$vcf_cons{$chr}}{$min_pos}}{$type});
						}
					}
					my $arg1 = ${${$vcf_cons{$chr}}{$pos}}{$type};
					my $arg2 = $arg1;
					$arg1 =~ s/SAMPLES=\S+/SAMPLES=$max_info/;
					$arg2 =~ s/SAMPLES=\S+/SAMPLES=$min_info/;
					delete ${${$vcf_cons{$chr}}{$pos}}{$type};
					delete ${$vcf_cons{$chr}}{$pos} if (scalar keys %{${$vcf_cons{$chr}}{$pos}} == 0);
					${${$vcf_cons{$chr}}{$max_pos}}{$type} = $arg1;
					${${$vcf_cons{$chr}}{$min_pos}}{$type} = $arg2;
					if (exists ${$ins_bp_convert{$chr2}}{$pos}){
						my $bplen = ${$ins_bp_convert{$chr2}}{$pos};
						${$ins_bp_convert{$chr2}}{$max_pos} = $bplen;
						${$ins_bp_convert{$chr2}}{$min_pos} = $bplen;
					}
				}
				elsif (($med_pos > 0) and ($min_pos > 0)){
					if ($med_pos == $min_pos){
						while (1){
							$min_pos --;
							last if (!exists ${${$vcf_cons{$chr}}{$min_pos}}{$type});
						}
					}
					my $arg1 = ${${$vcf_cons{$chr}}{$pos}}{$type};
					my $arg2 = $arg1;
					$arg1 =~ s/SAMPLES=\S+/SAMPLES=$med_info/;
					$arg2 =~ s/SAMPLES=\S+/SAMPLES=$min_info/;
					delete ${${$vcf_cons{$chr}}{$pos}}{$type};
					delete ${$vcf_cons{$chr}}{$pos} if (scalar keys %{${$vcf_cons{$chr}}{$pos}} == 0);
					${${$vcf_cons{$chr}}{$med_pos}}{$type} = $arg1;
					${${$vcf_cons{$chr}}{$min_pos}}{$type} = $arg2;
					if (exists ${$ins_bp_convert{$chr2}}{$pos}){
						my $bplen = ${$ins_bp_convert{$chr2}}{$pos};
						${$ins_bp_convert{$chr2}}{$med_pos} = $bplen;
						${$ins_bp_convert{$chr2}}{$min_pos} = $bplen;
					}
				}
				elsif (($max_pos > 0) and ($med_pos > 0)){
					if ($max_pos == $med_pos){
						while (1){
							$med_pos --;
							last if (!exists ${${$vcf_cons{$chr}}{$med_pos}}{$type});
						}
					}
					my $arg1 = ${${$vcf_cons{$chr}}{$pos}}{$type};
					my $arg2 = $arg1;
					$arg1 =~ s/SAMPLES=\S+/SAMPLES=$max_info/;
					$arg2 =~ s/SAMPLES=\S+/SAMPLES=$med_info/;
					delete ${${$vcf_cons{$chr}}{$pos}}{$type};
					delete ${$vcf_cons{$chr}}{$pos} if (scalar keys %{${$vcf_cons{$chr}}{$pos}} == 0);
					${${$vcf_cons{$chr}}{$max_pos}}{$type} = $arg1;
					${${$vcf_cons{$chr}}{$med_pos}}{$type} = $arg2;
					if (exists ${$ins_bp_convert{$chr2}}{$pos}){
						my $bplen = ${$ins_bp_convert{$chr2}}{$pos};
						${$ins_bp_convert{$chr2}}{$max_pos} = $bplen;
						${$ins_bp_convert{$chr2}}{$med_pos} = $bplen;
					}
				}
			}
		}
	}
}

my %vcf_cons2;
my %vcf_type2;
my %ins_str;
my $ins2str = 0;
my $str2ins = 0;
my $delete_str = 0;
my $delete_ins = 0;

foreach my $chr (sort keys %vcf_cons){
	my $chr2 = $chr;
	$chr2 =~ s/^0*//;
	foreach my $pos (sort {$a <=> $b} keys %{$vcf_cons{$chr}}){
		foreach my $type (keys %{${$vcf_cons{$chr}}{$pos}}){
			my $len = $1 if (${${$vcf_cons{$chr}}{$pos}}{$type} =~ /SVLEN2=-*(\d+)/);
			next if ($len == 0);
			my $id_info = $1 if (${${$vcf_cons{$chr}}{$pos}}{$type} =~ /SAMPLES=(\S+)/);
			my @id_info = split (/,,/, $id_info);
			my %gt;
			my %vrr;
			my %pos;
			my %len;
			my %dpr;
			my %bp;
			my %read;
			my @duppos = ();
			my @duplen = ();
			my @trapos;
			my @tralen;
			my %mei;
			my @meilen;
			my @len;
			my %insdup;
			my %traindel;
			my %strid;
			my %strid2;
			my %lowconf;
			my %inslen;
			my %ins_tr;
			my %ins_strunit;
			my $ins_str_count = 0;
			my @sar;
			foreach (@id_info){
				my ($id, $pos2, $len2, $alt, $inf2) = split (/==/, $_);
				$len2 = $len if ($len2 == 0);
				if (exists $len{$id}){
					my $pre_len = $len{$id};
					my $diff = abs ($len - $len2);
					my $pre_diff = abs ($len - $pre_len);
					next if ($diff >= $pre_diff);
				}
				push @len, $len2 if ($len2 > 0);
				$pos{$id} = $pos2;
				$len{$id} = $len2;
				my $gt = $1 if ($inf2 =~ /GT=(.+?);/);
				$gt{$id} = '0/1';
				$gt{$id} = '1/1' if ($gt eq 'HM');
				$gt{$id} = './.' if ($gt eq 'NA');
				my $vrr = $1 if ($inf2 =~ /VRR=([\d\.]+)/);
				$vrr{$id} = $vrr;
				my $dpr = 1;
				$dpr = $1 if ($inf2 =~ /DPR=([\d\.]+)/);
				$dpr{$id} = $dpr;
				my $read = $1 if ($inf2 =~ /READS=(\d+)/);
				$read{$id} = $read;
				my $bp = 0;
				$bp = $1 if ($inf2 =~ /BP=(\d+)/);
				$bp{$id} = $bp;
				my $inslen = 0;
				$inslen = $1 if ($inf2 =~ /INSLEN=(\d+)/);
				$inslen{$id} = $inslen if ($inslen > 0);;
				if ($inf2 =~ /MEI=(.+?);/){
					$mei{$1} ++;
				}
				if ($inf2 =~ /MEILEN=(\d+)/){
					push @meilen, $1;
				}
				if ($inf2 =~ /DUPPOS=(.+?);/){
					push @duppos, $1;
				}
				if ($inf2 =~ /DUPLEN=(\d+)/){
					push @duplen, $1;
				}
				if ($inf2 =~ /TRAPOS=(.+?);/){
					push @trapos, $1;
				}
				if ($inf2 =~ /TRALEN=(\d+)/){
					push @tralen, $1;
				}
				if (($alt=~ /INS/) and ($alt =~ /DUP/)){
					if ($alt =~ /intDUP:R/){
						$insdup{'intDUP:R'} ++;
					}
					elsif ($alt =~ /intDUP/){
						$insdup{'intDUP'} ++;
					}
					elsif ($alt =~ /DUP:R/){
						$insdup{'DUP:R'} ++;
					}
					else{
						$insdup{'DUP'} ++;
					}
				}
				if ($alt =~ /TRA:DEL/){
					$traindel{'DEL'} ++;
				}
				elsif ($alt =~ /TRA:INS/){
					$traindel{'INS'} ++;
				}
				if ($inf2 =~ /TRID=([^;~]+)/){
					$strid{$id} = $1;
					$strid2{$1} ++;
				}
				if ($inf2 =~ /LowConf/){
					$lowconf{$id} = 1;
				}
				if ($inf2 =~ /SAR=([\d\.]+)/){
					push @sar, $1;
				}
				if (($type eq 'INS') and ($inf2 =~ /TRUNIT=([^;]+)/)){
					my $ins_unit = $1;
					my ($umotif, $ucn) = split (/:/, $ins_unit);
					my $ulen = length ($umotif);
					$ins_tr{$id} = "$ulen-$ucn";
					${$ins_strunit{$ulen}}{$ucn} ++;
					$ins_str_count ++;
				}
			}
			my $sample_num = scalar keys %pos;
			my $hnum = int ($sample_num * 0.5);
			@len = sort {$a <=> $b} @len;
			my $med_len = $len[$hnum] if (@len > 0);
			$med_len = $len if (@len == 0);
			if (@len == 2){
				$med_len = int (($len[0] + $len[1]) * 0.5 + 0.5);
			}
			my $end = $pos + $med_len - 1;
			$end = $pos if ($type eq 'INS');
			$med_len = 0 - $med_len if ($type eq 'DEL');
			next if ($med_len == 0);
			my $ave_sar = 0;
			if (@sar > 0){
				my $sum_sar = 0;
				map{$sum_sar += $_} @sar;
				$ave_sar = int ($sum_sar / @sar * 100 + 0.5) / 100;
			}
			my $format = 'GT:VP:VL:VN:VR:BP:DR:TR';

			my $format_info = '';
			my %gt_info;
			my $ac = 0;
			my $sc = 0;
			foreach my $id (keys %pos){
				my $pos2 = $pos{$id};
				my $len2 = $len{$id};
				$len2 = $inslen{$id} if ($type eq 'DUP') and (exists $inslen{$id});
				my $gt = $gt{$id};
				my $read = $read{$id};
				my $vrr = $vrr{$id};
				my $bp = $bp{$id};
				my $dpr = $dpr{$id};
				my $info = "$gt:$pos2:$len2:$read:$vrr:$bp:$dpr";
				if (exists $ins_tr{$id}){
					$info .= ":$ins_tr{$id}";
				}
				else{
					$info .= ":0";
				}
				$gt_info{$id} = $info;
				$ac ++;
				$ac ++ if ($gt eq '1/1');
				$sc ++;
			}
			foreach my $id (sort keys %Gid){
				if (exists $gt_info{$id}){
					$format_info .= "$gt_info{$id}\t";
				}
				else{
					$format_info .= "0/0:0:0:0:0:0:0:0\t";
				}
			}
			$format_info =~ s/\t$//;
			my $new_line = '';
			my $alt = $type;
			if ($type eq 'INS'){
				my $meitype = '';
				my $duppos = '';
				my $duplen = 0;
				my $meilen = 0;
				if (@duplen > 0){
					if (@duplen == 1){
						$duplen = $duplen[0];
						$duppos = $duppos[0] if (@duppos > 0);
					}
					else{
						my $dup_hnum = int (@duplen * 0.5);
						@duplen = sort {$a <=> $b} @duplen;
						$duplen = $duplen[$dup_hnum];
						my %dupchr;
						my $dupchr = '';
						my %duppos2;
						my @duppos2 = ();
						foreach (@duppos){
							if ($_ =~ /:/){
								my ($chr2, $pos2) = split (/:/, $_);
								$dupchr{$chr2} ++;
								${$duppos2{$chr2}}{$pos2} ++;
							}
							else{
								push @duppos2, $_;
							}
						}
						if (scalar keys %dupchr > 0){
							foreach my $chr2 (sort {$dupchr{$b} <=> $dupchr{$a}} keys %dupchr){
								$dupchr = $chr2;
								last;
							}
							foreach my $pos2 (sort {$a <=> $b} keys %{$duppos2{$dupchr}}){
								my $num = ${$duppos2{$dupchr}}{$pos2};
								while ($num > 0){
									push @duppos2, $pos2;
									$num --;
								}
							}
						}
						my $duphnum = int (@duppos2 * 0.5);
						$duppos = $duppos2[$duphnum] if (@duppos2 > 0);
						$duppos = $dupchr . ":$duppos" if ($dupchr ne '') and ($duppos ne '');
					}
				}
				if (scalar keys %insdup > 0){
					foreach my $tag (sort {$insdup{$b} <=> $insdup{$a}} keys %insdup){
						$alt .= ":$tag";
						last;
					}
				}
				if (@meilen > 0){
					if (@meilen == 1){
						$meilen = $meilen[0];
					}
					else{
						my $mei_hnum = int (@meilen * 0.5);
						@meilen = sort {$a <=> $b} @meilen;
						$meilen = $meilen[$mei_hnum];
					}
					foreach my $mei (sort {$mei{$b} <=> $mei{$a}} keys %mei){
						$meitype = $mei;
						last;
					}
					$alt .= ":ME:$meitype";
				}
				$new_line = "$chr2\t$pos\t.\t.\t<$alt>\t.\tPASS\tSVTYPE=$type;SVLEN=$med_len";
				if (exists ${$ins_bplen2{$chr}}{$pos}){
					my @bplen = (@{${$ins_bplen2{$chr}}{$pos}});
					my $hn = int (@bplen * 0.5);
					@bplen = sort {$a <=> $b} @bplen;
					my $med_bplen = $bplen[$hn];
					$new_line .= ";BPLEN=$med_bplen";
				}
				if ($duppos ne ''){
					$new_line .= ";DUPPOS=$duppos;DUPLEN=$duplen";
				}
				if ($meitype ne ''){
					$new_line .= ";MEILEN=$meilen";
				}
				my $strid_count = scalar keys %strid;
				if (($ins_str_count / $sc >= 0.5) and ($strid_count / $sc < 0.5)){
					my $top_ulen = 0; 
					my $min_cn = 1000000;
					my $max_cn = 0;
					foreach my $ulen (sort {scalar keys %{$ins_strunit{$b}} <=> scalar keys %{$ins_strunit{$a}}} keys %ins_strunit){
						$top_ulen = $ulen;
						foreach my $cn (sort {$a <=> $b} keys %{$ins_strunit{$ulen}}){
							if ($cn < $min_cn){
								$min_cn = $cn;
							}
							if ($cn > $max_cn){
								$max_cn = $cn;
							}
						}
					}
					my $unit_range = "$min_cn-$max_cn";
					$unit_range = $min_cn if ($max_cn == $min_cn);
					$new_line .= ";TRUNIT=$top_ulen:$unit_range";
				}
				$new_line .= ";SAR=$ave_sar" if (@sar > 0);
				$new_line .= ";END=$end";
			}
			elsif ($type eq 'TRA'){
				my $trapos = '';
				my $tralen = 0;
				my $alt = 'TRA';
				if (@tralen > 0){
					if (@tralen == 1){
						$tralen = $tralen[0];
						$trapos = $trapos[0];
					}
					else{
						my $tra_hnum = int (@tralen * 0.5);
						@tralen = sort {$a <=> $b} @tralen;
						$tralen = $tralen[$tra_hnum];
						my %trachr;
						my $trachr = '';
						my %trapos2;
						my @trapos2;
						foreach (@trapos){
							my ($chr2, $pos2) = split (/:/, $_);
							$trachr{$chr2} ++;
							${$trapos2{$chr2}}{$pos2} ++;
						}
						foreach my $chr2 (sort {$trachr{$b} <=> $trachr{$a}} keys %trachr){
							$trachr = $chr2;
							last;
						}
						foreach my $pos2 (sort {$a <=> $b} keys %{$trapos2{$trachr}}){
							my $num = ${$trapos2{$trachr}}{$pos2};
							while ($num > 0){
								push @trapos2, $pos2;
								$num --;
							}
						}
						my $trahnum = int (@trapos2 * 0.5);
						$trapos = "$trachr:$trapos2[$trahnum]";
					}
				}
				if (scalar keys %traindel > 0){
					foreach my $tag (sort {$traindel{$b} <=> $traindel{$a}} keys %traindel){
						$alt .= ":$tag";
						last;
					}
				}
				$new_line = "$chr2\t$pos\t.\t.\t<$alt>\t.\tPASS\tSVTYPE=$type;SVLEN=$med_len";
				if ($trapos ne ''){
					$new_line .= "TRAPOS=$trapos;TRALEN=$tralen";
				}
				$new_line .= ";SAR=$ave_sar" if (@sar > 0);
				$new_line .= ";END=$end";
			}
			elsif ($type eq 'DUP'){
				$end = $pos + $len - 1 if ($end < $pos);
				$new_line = "$chr2\t$pos\t.\t.\t<$alt>\t.\tPASS\tSVTYPE=$type;SVLEN=$len;END=$end";
				$new_line = "$chr2\t$pos\t.\t.\t<$alt>\t.\tPASS\tSVTYPE=$type;SVLEN=$len;SAR=$ave_sar;END=$end" if (@sar > 0);
			}
			else{
				$end = $pos + $med_len - 1 if ($end < $pos);
				$new_line = "$chr2\t$pos\t.\t.\t<$alt>\t.\tPASS\tSVTYPE=$type;SVLEN=$med_len;END=$end";
				$new_line = "$chr2\t$pos\t.\t.\t<$alt>\t.\tPASS\tSVTYPE=$type;SVLEN=$med_len;SAR=$ave_sar;END=$end" if (@sar > 0);
			}
			my $strid_count = scalar keys %strid;
			my $lowconf_count = scalar keys %lowconf;
			if (($strid_count > 0) and ($type eq 'INS') and (@meilen == 0)){
				my $top_strid = '';
				foreach my $strid (sort {$strid2{$b} <=> $strid2{$a}} keys %strid2){
					$top_strid = $strid;
					last;
				}
				if (exists $cnv_info{$top_strid}){
					my ($schr, $spos, $strend, $strulen) = split (/=/, $cnv_info{$top_strid});
					my $match = 0;
					my %match;
					foreach my $id (keys %{${$cnv_line{$schr}}{$spos}}){
						my $cnv_line = ${${$cnv_line{$schr}}{$spos}}{$id};
						my $stype = $1 if ($cnv_line =~ /SVTYPE=(.+?);/);
						my $stype1 = $stype;
						my $stype2 = '';
						($stype1, $stype2) = split (/,/, $stype) if ($stype =~ /,/);
						my $vlen = $1 if ($cnv_line =~ /SVLEN=(.+?);/);
						my $svrr = $1 if ($cnv_line =~ /VRR=(.+?);/);
						my $sread = $1 if ($cnv_line =~ /READS=(.+?);/);
						my $sgt = $1 if ($cnv_line =~ /GT=(.+?);/);
						$sgt = 'HT' if ($sgt eq 'HT2');
						if ($vlen =~ /,/){
							my ($vlen1, $vlen2) = split (/,/, $vlen);
							$vlen1 = 0 - $vlen1 if ($vlen1 < 0);
							$vlen2 = 0 - $vlen2 if ($vlen2 < 0);
							if (($stype1 eq 'INS') and ($vlen1 >= $len * 0.8) and ($vlen1 <= $len * 1.25)){
								$match ++;
								my ($svrr1, $svrr2) = split (/,/, $svrr);
								my ($sread1, $sread2) = split (/,/, $sread);
								$match{$id} = "$vlen1=$svrr1=$sread1=$sgt";
							}
							elsif (($stype1 eq 'INS') and ($stype !~ /,/) and ($vlen2 >= $len * 0.8) and ($vlen2 <= $len * 1.25)){
								$match ++;
								my ($svrr1, $svrr2) = split (/,/, $svrr);
								my ($sread1, $sread2) = split (/,/, $sread);
								$match{$id} = "$vlen2=$svrr2=$sread2=$sgt";
							}
							elsif (($stype2 eq 'INS') and ($stype =~ /,/) and ($vlen2 >= $len * 0.8) and ($vlen2 <= $len * 1.25)){
								$match ++;
								my ($svrr1, $svrr2) = split (/,/, $svrr);
								my ($sread1, $sread2) = split (/,/, $sread);
								$match{$id} = "$vlen2=$svrr2=$sread2=$sgt";
							}
						}
						elsif ($stype1 eq 'INS'){
							if (($vlen > 0) and ($vlen >= $len * 0.8) and ($vlen <= $len * 1.25)){
								$match ++;
								$match{$id} = "$vlen=$svrr=$sread=$sgt";
							}
						}
					}
					if ($match >= $sc){						# INS within STR region is merged to STR-INS
						foreach my $id (keys %pos){
							my $slen = $len{$id};
							my $svrr = $vrr{$id};
							my $sread = $read{$id};
							my $sgt = $gt{$id};
							my $cn = int ($slen / $strulen * 10 + 0.5) / 10;
							if (!exists ${${$cnv_line{$schr}}{$spos}}{$id}){
								my $new_line2 = "$chr2\t$spos\t.\t.\t<CNV:TR>\t.\tPASS\tSVTYPE=INS;SVLEN=$slen;READS=$sread;CN=gain+$cn;VRR=$svrr;SAR=0;GT=$sgt;END=$spos;TRID=$top_strid;TREND=$strend;TRULEN=$strulen";
								${${$cnv_line{$schr}}{$spos}}{$id} = $new_line2;
							}
							else{
								my $cnv_line = ${${$cnv_line{$schr}}{$spos}}{$id};
								my $vlen = $1 if ($cnv_line =~ /SVLEN=(.+?);/);
								if ($vlen !~ /,/){
									$vlen = 0 - $vlen if ($vlen < 0);
									my $vtype = $1 if ($cnv_line =~ /SVTYPE=(.+?);/);
									my $vvrr = $1 if ($cnv_line =~ /VRR=(.+?);/);
									my $vread = $1 if ($cnv_line =~ /READS=(.+?);/);
									my $vcn = int ($vlen / $strulen * 10 + 0.5) / 10;
									my $cn = int ($vlen / $strulen * 10 + 0.5) / 10;
									if (($vtype eq 'INS') and (($vlen / $slen < 0.8) or ($vlen / $slen > 1.25))){
										$cnv_line =~ s/SVLEN=$vlen/SVLEN=$vlen,$slen/;
										$cnv_line =~ s/VRR=$vvrr/VRR=$vvrr,$svrr/;
										$cnv_line =~ s/READS=$vread/READS=$vread,$sread/;
										$cnv_line =~ s/GT=.+?;/GT=HT2;/;
										$cnv_line =~ s/CN=.+?;/CN=gain+$vcn,gain+$cn;/;
										${${$cnv_line{$schr}}{$spos}}{$id} = $cnv_line;
									}
									elsif ($vtype eq 'DEL'){
										$cnv_line =~ s/SVLEN=$vlen/SVLEN=$vlen,$slen/;
										$cnv_line =~ s/VRR=$vvrr/VRR=$vvrr,$svrr/;
										$cnv_line =~ s/READS=$vread/READS=$vread,$sread/;
										$cnv_line =~ s/GT=.+?;/GT=HT2;/;
										$cnv_line =~ s/CN=.+?;/CN=loss-$vcn,gain+$cn;/;
										$cnv_line =~ s/SVTYPE=DEL/SVTYPE=DEL,INS/;
										${${$cnv_line{$schr}}{$spos}}{$id} = $cnv_line;
									}
								}
							}
							$ins2str ++;
						}
						$delete_ins ++;
#print STDERR "$chr:$pos INS-$len was incorporated into STR-CNV ($top_strid): $schr:$spos ($match >= $sc)\n";
						next;
					}
					elsif ($strid_count / $sc >= 0.5){		# STR-INS is merged to INS with no STR copies
						$new_line .= ";TRID=$top_strid";
#						if (($top_strid =~ /^TR/) and ($type eq 'INS')){
#							${$ins_str{$top_strid}}{"$chr=$pos=$med_len"} = $ac;
#						}
						if ($strid_count >= $match){
							$format_info = '';
							foreach my $id (sort keys %Gid){
								if (exists $gt_info{$id}){
									$format_info .= "$gt_info{$id}\t";
								}
								elsif (exists $match{$id}){
									my ($slen, $svrr, $sread, $sgt) = split (/=/, $match{$id});
									$format_info .= "$sgt:$pos:$slen:$sread:$svrr:0:1:0\t";
									$str2ins ++;
								}
								else{
									$format_info .= "0/0:0:0:0:0:0:0:0\t";
								}
							}
							$format_info =~ s/\t$//;
#print STDERR "$top_strid ($schr:$spos) INS was incorporated into $chr:$pos INS-$len ($match < $sc $strid_count)\n";
							foreach my $id (keys %{${$cnv_line{$schr}}{$spos}}){
								next if (!exists $match{$id});
								my $cnv_line = ${${$cnv_line{$schr}}{$spos}}{$id};
								my $stype = $1 if ($cnv_line =~ /SVTYPE=(.+?);/);
								my $stype1 = $stype;
								my $stype2 = '';
								($stype1, $stype2) = split (/,/, $stype) if ($stype =~ /,/);
								my $vlen = $1 if ($cnv_line =~ /SVLEN=(.+?);/);
								my $svrr = $1 if ($cnv_line =~ /VRR=(.+?);/);
								my $sread = $1 if ($cnv_line =~ /READS=(.+?);/);
								my $sgt = $1 if ($cnv_line =~ /GT=(.+?);/);
								my ($slen) = split (/=/, $match{$id});
								if ($vlen =~ /,/){
									my ($vlen1, $vlen2) = split (/,/, $vlen);
									$vlen1 = 0 - $vlen1 if ($vlen1 < 0);
									$vlen2 = 0 - $vlen2 if ($vlen2 < 0);
									if (($stype1 eq 'INS') and ($vlen1 == $slen)){
										my ($svrr1, $svrr2) = split (/,/, $svrr);
										my ($sread1, $sread2) = split (/,/, $sread);
										$cnv_line =~ s/SVLEN=.+?;/SVLEN=$vlen2;/;
										$cnv_line =~ s/VRR=.+?;/VRR=$svrr2;/;
										$cnv_line =~ s/READS=.+?;/READS=$sread2;/;
										$cnv_line =~ s/GT=.+?;/GT=HT;/;
										my $cn = int ($vlen2 / $strulen * 10 + 0.5) / 10;
										if ($vlen2 < 0){
											$cnv_line =~ s/SVTYPE=.+?;/SVTYPE=DEL;/;
											$cnv_line =~ s/CN=.+?;/CN=loss$cn;/;
										}
										else{
											$cnv_line =~ s/CN=.+?;/CN=gain+$cn;/;
										}
										${${$cnv_line{$schr}}{$spos}}{$id} = $cnv_line;
									}
									elsif ((($stype2 eq 'INS') or (($stype1 eq 'INS') and ($stype2 eq ''))) and ($vlen2 == $slen)){
										my ($svrr1, $svrr2) = split (/,/, $svrr);
										my ($sread1, $sread2) = split (/,/, $sread);
										$cnv_line =~ s/SVLEN=.+?;/SVLEN=$vlen1;/;
										$cnv_line =~ s/VRR=.+?;/VRR=$svrr1;/;
										$cnv_line =~ s/READS=.+?;/READS=$sread1;/;
										$cnv_line =~ s/GT=.+?;/GT=HT;/;
										my $cn = int ($vlen1 / $strulen * 10 + 0.5) / 10;
										if ($vlen1 < 0){
											$cnv_line =~ s/SVTYPE=.+?;/SVTYPE=DEL;/;
											$cnv_line =~ s/CN=.+?;/CN=loss$cn;/;
										}
										else{
											$cnv_line =~ s/CN=.+?;/CN=gain+$cn;/;
										}
										${${$cnv_line{$schr}}{$spos}}{$id} = $cnv_line;
									}
								}
								else{
									$vlen = 0 - $vlen if ($vlen < 0);
									if (($stype1 eq 'INS') and ($vlen == $slen)){
										delete ${${$cnv_line{$schr}}{$spos}}{$id};
									}
								}
								if (scalar keys %{${$cnv_line{$schr}}{$spos}} == 0){
									delete ${$cnv_line{$schr}}{$spos};
									delete $cnv_info{$top_strid};
									$delete_str ++;
print STDERR "$top_strid: $schr:$spos: was deleted\n";
								}
							}
						}
					}
				}
				else{
					if ($strid_count / $sc >= 0.5){
						$new_line .= ";TRID=$top_strid";
#						if (($top_strid =~ /^TR/) and ($type eq 'INS')){
#							${$ins_str{$top_strid}}{"$chr=$pos=$med_len"} = $ac;
#						}
					}
				}
			}
			if ($lowconf_count / $sc >= 0.5){
				my @new_line = split (/\t/, $new_line);
				$new_line[6] = 'LowConf';
				$new_line[6] = 'LowConf,LowQual' if ($ave_sar > $max_SAR);
				$new_line = join ("\t", @new_line);
			}
			elsif ($ave_sar > $max_SAR){
				my @new_line = split (/\t/, $new_line);
				$new_line[6] = 'LowQual';
				$new_line = join ("\t", @new_line);
			}
			$new_line .= "\t$format\t$format_info";
			${${$vcf_cons2{$chr}}{$pos}}{$type} = $new_line;
			${${$vcf_type2{$type}}{$chr}}{$pos} = $new_line;
		}
	}
}
%vcf_cons = ();

print STDERR "INS to STR converted: $ins2str\n";
print STDERR "STR to INS converted: $str2ins\n";
print STDERR "Deleted INS: $delete_ins\n";
print STDERR "Deleted STR-INS: $delete_str\n";


foreach my $type (keys %vcf_type2){			# merge or rearrange alleles of proximal sites depending on their lengths
	foreach my $chr (keys %{$vcf_type2{$type}}){
		my $pre_pos = 0;
		my $pre_len = 0;
		my @mpos;
		my %removed;
		foreach my $pos (sort {$a <=> $b} keys %{${$vcf_type2{$type}}{$chr}}){
			next if (exists $removed{$pos});
			my $len = $1 if (${${$vcf_type2{$type}}{$chr}}{$pos} =~ /SVLEN=-*(\d+)/);
			if (($pre_pos > 0) and ($pos - $pre_pos <= 100) and ($len > 0)){
				if ((($len >= 50) or ($pre_len >= 50)) and ($pre_len / $len < 1.8) and ($pre_len / $len > 0.556)){
					push @mpos, $pre_pos if (@mpos == 0);
					push @mpos, $pos;
					next;
				}
				if (($pos - $pre_pos <= $len) and ($pos - $pre_pos <= $pre_len) and ($pre_len / $len < 1.8) and ($pre_len / $len > 0.556)){
					push @mpos, $pre_pos if (@mpos == 0);
					push @mpos, $pos;
					next;
				}
			}
			elsif ($pre_pos > 0){
				if (@mpos >= 2){
					foreach my $pos1 (@mpos){
						next if (exists $removed{$pos1});
						my $line1 = ${${$vcf_type2{$type}}{$chr}}{$pos1};
						my @line1 = split (/\t/, $line1);
						my $len1 = $1 if ($line1[7] =~ /SVLEN=-*(\d+)/);
						next if ($len1 == 0);
						foreach my $pos2 (@mpos){
							next if (exists $removed{$pos2});
							next if ($pos2 <= $pos1);
							my $line2 = ${${$vcf_type2{$type}}{$chr}}{$pos2};
							my @line2 = split (/\t/, $line2);
							my $len2 = $1 if ($line2[7] =~ /SVLEN=-*(\d+)/);
							next if ($len2 == 0);
							if (($len1 / $len2 < 1.8) and ($len1 / $len2 > 0.556)){
								my %ids1;
								my %ids2;
								my %ids;
								my $sc1 = 0;
								my $sc2 = 0;
								my $count = 0;
								foreach (@line1){
									$count ++;
									next if ($count <= 9);
									my ($gt, $vp, $vl, $vn, $vr) = split (/:/, $_);
									next if ($_ =~ /^0\/0/);
									$ids{$count} = 1;
									$ids1{$count} = 1;
									$sc1 ++;
								}
								$count = 0;
								foreach (@line2){
									$count ++;
									next if ($count <= 9);
									my ($gt, $vp, $vl, $vn, $vr) = split (/:/, $_);
									next if ($_ =~ /^0\/0/);
									$ids{$count} = 1;
									$ids2{$count} = 1;
									$sc2 ++;
								}
								foreach my $cnt (keys %ids){
									if (($sc1 >= $sc2) and (!exists $ids1{$cnt})){
										$line1[$cnt - 1] = $line2[$cnt - 1];
									}
									elsif (($sc1 < $sc2) and (!exists $ids2{$cnt})){
										$line2[$cnt - 1] = $line1[$cnt - 1];
									}
									elsif ((exists $ids1{$cnt}) and (exists $ids2{$cnt})){
										my ($gt1, $vp1, $vl1, $vn1, $vr1, $bp1, $dr1, $tr1) = split (/:/, $line1[$cnt - 1]);
										my ($gt2, $vp2, $vl2, $vn2, $vr2, $bp2, $dr2, $tr2) = split (/:/, $line2[$cnt - 1]);
										if (($vp1 != $vp2) or ($vn1 != $vn2)){
											my $new_vn = $vn1 + $vn2;
											my $new_bp = $bp1 + $bp2;
											my $new_vl = int (($vl1 * $vn1 + $vl2 * $vn2) / ($vn1 + $vn2) + 0.5);
											my $new_vr = $vr1 + $vr2;
											$new_vr = 1 if ($new_vr > 1);
											my $new_gt = '0/1';
											$new_gt = '1/1' if ($new_vr >= 0.85);
											if ($sc1 >= $sc2){
												$line1[$cnt - 1] = "$new_gt:$vp1:$new_vl:$new_vn:$new_vr:$new_bp:$dr1:$tr1";
											}
											else{
												$line2[$cnt - 1] = "$new_gt:$vp2:$new_vl:$new_vn:$new_vr:$new_bp:$dr2:$tr2";
											}
										}
									}
								}
								if ($sc1 >= $sc2){
									my @len;
									my $count = 0;
									foreach (@line1){
										$count ++;
										next if ($count <= 9);
										next if ($_ =~ /^0\/0/);
										my ($gt, $vp, $vl) = split (/:/, $_);
										push @len, $vl;
									}
									my $hn = int (@len * 0.5);
									@len = sort {$a <=> $b} @len;
									my $medlen = $len[$hn];
									my $end = $pos2 + $medlen - 1;
									if ($medlen != $len1){
										$medlen = 0 - $medlen if ($type eq 'DEL');
										$line1[7] =~ s/SVLEN=-*\d+/SVLEN=$medlen/;
										$line1[7] =~ s/END=\d+/END=$end/ if ($type ne 'INS');
									}
									delete ${${$vcf_type2{$type}}{$chr}}{$pos2};
									delete ${${$vcf_cons2{$chr}}{$pos2}}{$type};
									$removed{$pos2} = 1;
									$line1 = join ("\t", @line1);
									${${$vcf_type2{$type}}{$chr}}{$pos1} = $line1;
									${${$vcf_cons2{$chr}}{$pos1}}{$type} = $line1;
								}
								else{
									my @len;
									my $count = 0;
									foreach (@line2){
										$count ++;
										next if ($count <= 9);
										next if ($_ =~ /^0\/0/);
										my ($gt, $vp, $vl) = split (/:/, $_);
										push @len, $vl;
									}
									my $hn = int (@len * 0.5);
									@len = sort {$a <=> $b} @len;
									my $medlen = $len[$hn];
									my $end = $pos2 + $medlen - 1;
									if ($medlen != $len2){
										$medlen = 0 - $medlen if ($type eq 'DEL');
										$line2[7] =~ s/SVLEN=-*\d+/SVLEN=$medlen/;
										$line2[7] =~ s/END=\d+/END=$end/ if ($type ne 'INS');
									}
									delete ${${$vcf_type2{$type}}{$chr}}{$pos1};
									delete ${${$vcf_cons2{$chr}}{$pos1}}{$type};
									$removed{$pos1} = 1;
									$line2 = join ("\t", @line2);
									${${$vcf_type2{$type}}{$chr}}{$pos2} = $line2;
									${${$vcf_cons2{$chr}}{$pos2}}{$type} = $line2;
									last;
								}
							}
							else{
								my %ids1;
								my %ids2;
								my %ids;
								my $sc1 = 0;
								my $sc2 = 0;
								my $count = 0;
								my $flag = 0;
								foreach (@line1){
									$count ++;
									next if ($count <= 9);
									my ($gt, $vp, $vl, $vn, $vr) = split (/:/, $_);
									next if ($_ =~ /^0\/0/);
									$ids1{$count} = 1;
									$ids{$count} = 1;
									$sc1 ++;
								}
								$count = 0;
								foreach (@line2){
									$count ++;
									next if ($count <= 9);
									my ($gt, $vp, $vl, $vn, $vr) = split (/:/, $_);
									next if ($_ =~ /^0\/0/);
									$ids2{$count} = 1;
									$ids{$count} = 1;
									$sc2 ++;
								}
								foreach my $cnt (keys %ids){
									my $ID = $Gid_order{$cnt};
									my $sex = '';
									$sex = $gender{$ID} if (exists $gender{$ID});
									if ((exists $ids1{$cnt}) and (exists $ids2{$cnt})){
										my ($gt1, $vp1, $vl1, $vn1, $vr1, $bp1, $dr1, $tr1) = split (/:/, $line1[$cnt - 1]);
										my ($gt2, $vp2, $vl2, $vn2, $vr2, $bp2, $dr2, $tr2) = split (/:/, $line2[$cnt - 1]);
										my $diff1 = int (abs ($vl1 - $len1) / $len1 * 100 + 0.5) / 100;
										my $diff2 = int (abs ($vl1 - $len2) / $len2 * 100 + 0.5) / 100;
										my $diff3 = int (abs ($vl2 - $len1) / $len1 * 100 + 0.5) / 100;
										my $diff4 = int (abs ($vl2 - $len2) / $len2 * 100 + 0.5) / 100;
										if (($diff1 > $diff2) and ($diff3 < $diff4)){
											my $info1 = $line1[$cnt - 1];
											my $info2 = $line2[$cnt - 1];
											$line1[$cnt - 1] = $info2;
											$line2[$cnt - 1] = $info1;
											$flag = 1;
										}
										elsif ($diff1 > $diff2){
											if (($vp1 != $vp2) or ($vn1 != $vn2)){
												my $new_vn = $vn1 + $vn2;
												my $new_bp = $bp1 + $bp2;
												my $new_vr = $vr1 + $vr2;
												$new_vr = 1 if ($new_vr > 1);
												my $new_vl = $vl2;
												$new_vl = $vl1 if ($vn1 > $vn2);
												my $new_vp = $vp2;
												$new_vp = $vp1 if ($vn1 > $vn2);
												my $new_gt = '0/1';
												$new_gt = '1/1' if ($new_vr >= 0.85);
												$line2[$cnt - 1] = "$new_gt:$new_vp:$new_vl:$new_vn:$new_vr:$new_bp:$dr2:$tr2";
											}
											$line1[$cnt - 1] = '0/0:0:0:0:0:0:0:0';
											$flag = 1;
										}
										elsif ($diff3 < $diff4){
											if (($vp1 != $vp2) or ($vn1 != $vn2)){
												my $new_vn = $vn1 + $vn2;
												my $new_bp = $bp1 + $bp2;
												my $new_vr = $vr1 + $vr2;
												$new_vr = 1 if ($new_vr > 1);
												my $new_vl = $vl1;
												$new_vl = $vl2 if ($vn1 < $vn2);
												my $new_vp = $vp1;
												$new_vp = $vp2 if ($vn1 < $vn2);
												my $new_gt = '0/1';
												$new_gt = '1/1' if ($new_vr >= 0.85);
												$line1[$cnt - 1] = "$new_gt:$new_vp:$new_vl:$new_vn:$new_vr:$new_bp:$dr1:$tr1";
											}
											$line2[$cnt - 1] = '0/0:0:0:0:0:0:0:0';
											$flag = 1;
										}
									}
									elsif (exists $ids1{$cnt}){
										my ($gt1, $vp1, $vl1) = split (/:/, $line1[$cnt - 1]);
										my $diff1 = int (abs ($vl1 - $len1) / $len1 * 100 + 0.5) / 100;
										my $diff2 = int (abs ($vl1 - $len2) / $len2 * 100 + 0.5) / 100;
										if ($diff1 > $diff2){
											$line2[$cnt - 1] = $line1[$cnt - 1];
											$line1[$cnt - 1] = '0/0:0:0:0:0:0:0:0';
											$flag = 1;
										}
									}
									elsif (exists $ids2{$cnt}){
										my ($gt2, $vp2, $vl2) = split (/:/, $line2[$cnt - 1]);
										my $diff1 = int (abs ($vl2 - $len1) / $len1 * 100 + 0.5) / 100;
										my $diff2 = int (abs ($vl2 - $len2) / $len2 * 100 + 0.5) / 100;
										if ($diff1 < $diff2){
											$line1[$cnt - 1] = $line2[$cnt - 1];
											$line2[$cnt - 1] = '0/0:0:0:0:0:0:0:0';
											$flag = 1;
										}
									}
								}
								if ($flag == 1){
									my $count = 0;
									my @len1;
									my @len2;
									foreach (@line1){
										$count ++;
										next if ($count <= 9);
										next if ($_ =~ /^0\/0/);
										my ($gt, $vp, $vl) = split (/:/, $_);
										push @len1, $vl;
									}
									$count = 0;
									foreach (@line2){
										$count ++;
										next if ($count <= 9);
										next if ($_ =~ /^0\/0/);
										my ($gt, $vp, $vl) = split (/:/, $_);
										push @len2, $vl;
									}
									my $hn1 = int (@len1 * 0.5);
									my $hn2 = int (@len2 * 0.5);
									if (@len1 == 0){
										delete ${${$vcf_type2{$type}}{$chr}}{$pos1};
										delete ${${$vcf_cons2{$chr}}{$pos1}}{$type};
										$removed{$pos1} = 1;
									}
									if (@len2 == 0){
										delete ${${$vcf_type2{$type}}{$chr}}{$pos2};
										delete ${${$vcf_cons2{$chr}}{$pos2}}{$type};
										$removed{$pos2} = 1;
									}
									if (@len1 > 0){
										my $medlen1 = $len1[$hn1];
										$medlen1 = 0 - $medlen1 if ($type eq 'DEL');
										my $end = $pos1 + $medlen1 - 1;
										if ($medlen1 != $len1){
											$line1[7] =~ s/SVLEN=-*\d+/SVLEN=$medlen1/;
											$line1[7] =~ s/END=\d+/END=$end/ if ($type ne 'INS');
										}
										$line1 = join ("\t", @line1);
										${${$vcf_type2{$type}}{$chr}}{$pos1} = $line1;
										${${$vcf_cons2{$chr}}{$pos1}}{$type} = $line1;
									}
									if (@len2 > 0){
										my $medlen2 = $len2[$hn2];
										$medlen2 = 0 - $medlen2 if ($type eq 'DEL');
										my $end = $pos2 + $medlen2 - 1;
										if ($medlen2 != $len2){
											$line2[7] =~ s/SVLEN=-*\d+/SVLEN=$medlen2/;
											$line2[7] =~ s/END=\d+/END=$end/ if ($type ne 'INS');
										}
										$line2 = join ("\t", @line2);
										${${$vcf_type2{$type}}{$chr}}{$pos2} = $line2;
										${${$vcf_cons2{$chr}}{$pos2}}{$type} = $line2;
									}
									last if (@len1 == 0);
								}
							}
						}
					}
				}
				@mpos = ();
			}
			$pre_pos = $pos;
			$pre_len = $len;
		}
	}
}
%vcf_type2 = (); 

foreach my $chr (keys %ins_bp_cons){
	my $chr02d = $chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	foreach my $pos (sort {$a <=> $b} keys %{$ins_bp_cons{$chr}}){
		my %gt;
		my %vrr;
		my %pos;
		my %bplen;
		my %bp;
		my %read;
		my @duppos = ();
		my @duplen = ();
		my @bplen;
		my %mei;
		my @meilen;
		my @len;
		my %insdup;
		my $sample_num = scalar @{${$ins_bp_cons{$chr}}{$pos}};
		my $hnum = int ($sample_num * 0.5);
		my %strid;
		my %strid2;
		my %lowconf;
		my @sar;
		foreach (@{${$ins_bp_cons{$chr}}{$pos}}){
			my ($id, $pos2, $len2, $alt, $inf2) = split (/==/, $_);
			$pos{$id} = $pos2;
			my $gt = $1 if ($inf2 =~ /GT=(.+?);/);
			$gt{$id} = '0/1';
			$gt{$id} = '1/1' if ($gt eq 'HM');
			$gt{$id} = './.' if ($gt eq 'NA');
			my $vrr = $1 if ($inf2 =~ /VRR=([\d\.]+)/);
			$vrr{$id} = $vrr;
			my $read = $1 if ($inf2 =~ /READS=(\d+)/);
			$read{$id} = $read;
			my $bp = 0;
			$bp = $1 if ($inf2 =~ /BP=(\d+)/);
			$bp{$id} = $bp;
			if ($inf2 =~ /MEI=(.+?);/){
				$mei{$1} ++;
			}
			if ($inf2 =~ /MEILEN=(\d+)/){
				push @meilen, $1;
			}
			if ($inf2 =~ /DUPPOS=(.+?);/){
				push @duppos, $1;
			}
			if ($inf2 =~ /DUPLEN=(\d+)/){
				push @duplen, $1;
			}
			if ($alt =~ /DUP/){
				if ($alt =~ /intDUP:R/){
					$insdup{'intDUP:R'} ++;
				}
				elsif ($alt =~ /intDUP/){
					$insdup{'intDUP'} ++;
				}
				elsif ($alt =~ /DUP:R/){
					$insdup{'DUP:R'} ++;
				}
				else{
					$insdup{'DUP'} ++;
				}
			}
			if ($inf2 =~ /TRID=([^;~]+)/){
				$strid{$id} = $1;
				$strid2{$1} ++;
			}
			if ($inf2 =~ /LowConf/){
				$lowconf{$id} = 1;
			}
			if ($inf2 =~ /BPLEN=(\d+)/){
				$bplen{$id} = $1;
				push @bplen, $1;
			}
			else{
				$bplen{$id} = 0;
				push @bplen, 0;
			}
			if ($inf2 =~ /SAR=([\d\.]+)/){
				push @sar, $1;
			}
		}
		my $med_len = $len[$hnum];
		if (@len == 2){
			$med_len = int (($len[0] + $len[1]) * 0.5 + 0.5);
		}
		my $ave_sar = 0;
		if (@sar > 0){
			my $sum_sar = 0;
			map{$sum_sar += $_} @sar;
			$ave_sar = int ($sum_sar / @sar * 100 + 0.5) / 100;
		}
		my $end = $pos;
		my $format = 'GT:VP:VL:VN:VR:BP:DR:TR';

		my $format_info = '';
		my %gt_info;
		foreach my $id (keys %pos){
			my $pos2 = $pos{$id};
			my $len2 = $bplen{$id};
			my $gt = $gt{$id};
			my $read = $read{$id};
			my $vrr = $vrr{$id};
			my $bp = $bp{$id};
			my $info = "$gt:$pos2:$len2:$read:$vrr:$bp:0:0";
			$gt_info{$id} = $info;
		}
		foreach my $id (sort keys %Gid){
			if (exists $gt_info{$id}){
				$format_info .= "$gt_info{$id}\t";
			}
			else{
				$format_info .= "0/0:0:0:0:0:0:0:0\t";
			}
		}
		$format_info =~ s/\t$//;
		my $new_line = '';
		my $alt = 'INS:BP';
		my $meitype = '';
		my $duppos = '';
		my $duplen = 0;
		my $meilen = 0;
		my $bplen = 0;
		if (@duplen > 0){
			if (@duplen == 1){
				$duplen = $duplen[0];
				$duppos = $duppos[0] if (@duppos > 0);
			}
			else{
				my $dup_hnum = int (@duplen * 0.5);
				@duplen = sort {$a <=> $b} @duplen;
				$duplen = $duplen[$dup_hnum];
				my %dupchr;
				my $dupchr = '';
				my %duppos2;
				my @duppos2 = ();
				foreach (@duppos){
					if ($_ =~ /:/){
						my ($chr2, $pos2) = split (/:/, $_);
						$dupchr{$chr2} ++;
						${$duppos2{$chr2}}{$pos2} ++;
					}
					else{
						push @duppos2, $_;
					}
				}
				if (scalar keys %dupchr > 0){
					foreach my $chr2 (sort {$dupchr{$b} <=> $dupchr{$a}} keys %dupchr){
						$dupchr = $chr2;
						last;
					}
					foreach my $pos2 (sort {$a <=> $b} keys %{$duppos2{$dupchr}}){
						my $num = ${$duppos2{$dupchr}}{$pos2};
						while ($num > 0){
							push @duppos2, $pos2;
							$num --;
						}
					}
				}
				my $duphnum = int (@duppos2 * 0.5);
				$duppos = $duppos2[$duphnum] if (@duppos2 > 0);
				$duppos = $dupchr . ":$duppos" if ($dupchr ne '') and ($duppos ne '');
			}
		}
		if (scalar keys %insdup > 0){
			foreach my $tag (sort {$insdup{$b} <=> $insdup{$a}} keys %insdup){
				$alt .= ":$tag";
				last;
			}
		}
		if (@meilen > 0){
			if (@meilen == 1){
				$meilen = $meilen[0];
			}
			else{
				my $mei_hnum = int (@meilen * 0.5);
				@meilen = sort {$a <=> $b} @meilen;
				$meilen = $meilen[$mei_hnum];
			}
			foreach my $mei (sort {$mei{$b} <=> $mei{$a}} keys %mei){
				$meitype = $mei;
				last;
			}
			$alt .= ":ME:$meitype";
		}
		@bplen = sort {$a <=> $b} @bplen;
		my $med_bplen = $bplen[$hnum];
		$new_line = "$chr\t$pos\t.\t.\t<$alt>\t.\tPASS\tSVTYPE=INS;SVLEN=0;BPLEN=$med_bplen";
		if ($duppos ne ''){
			$new_line .= ";DUPPOS=$duppos;DUPLEN=$duplen";
		}
		if ($meitype ne ''){
			$new_line .= ";MEILEN=$meilen";
		}
		$new_line .= ";SAR=$ave_sar" if (@sar > 0);
		$new_line .= ";END=$end";
		my $sc = scalar keys %gt_info;
		my $strid_count = scalar keys %strid;
		my $lowconf_count = scalar keys %lowconf;
		if ($strid_count / $sc >= 0.5){
			my $top_strid = '';
			foreach my $strid (sort {$strid2{$b} <=> $strid2{$a}} keys %strid2){
				$top_strid = $strid;
				last;
			}
			$new_line .= ";TRID=$top_strid"
		}
		if ($lowconf_count / $sc >= 0.5){
			my @new_line = split (/\t/, $new_line);
			$new_line[6] = 'LowConf';
			$new_line[6] = 'LowConf,LowQual' if ($ave_sar > $max_SAR);
			$new_line = join ("\t", @new_line);
		}
		elsif ($ave_sar > $max_SAR){
			my @new_line = split (/\t/, $new_line);
			$new_line[6] = 'LowQual';
			$new_line = join ("\t", @new_line);
		}
		$new_line .= "\t$format\t$format_info";
		${$ins_bp_cons2{$chr02d}}{$pos} = $new_line;
	}
}
%ins_bp_cons = ();

foreach my $chr (keys %ins_bp_cons2){
	foreach my $pos (keys %{$ins_bp_cons2{$chr}}){
		${${$vcf_cons2{$chr}}{$pos}}{'INS'} = ${$ins_bp_cons2{$chr}}{$pos};
	}
}

my $delete_str2 = 0;
my $delete_str3 = 0;
my $delete_ins2 = 0;
foreach my $chr (keys %cnv_line){
	my $chr02d = $chr;
	$chr02d = sprintf ("%02d", $chr) if ($chr =~ /^\d+$/);
	foreach my $pos (sort {$a <=> $b} keys %{$cnv_line{$chr}}){
		my %gt;
		my %vrr;
		my %len;
		my %bp;
		my %read;
		my %cn;
		my %ctype;
		my @inslen;
		my @dellen;
		my @sar;
		my $strid = '';
		my $strend = 0;
		my $strulen = 0;
		my $ins_ac = 0;
		foreach my $id (keys %{${$cnv_line{$chr}}{$pos}}){
			my $line = ${${$cnv_line{$chr}}{$pos}}{$id};
			my @line = split (/\t/, $line);
			my $type = $1 if ($line[7] =~ /SVTYPE=(.+?);/);
			my $len = $1 if ($line[7] =~ /SVLEN=-*([\d,]+)/);
			my $read = $1 if ($line[7] =~ /READS=([\d,]+)/);
			my $vrr = $1 if ($line[7] =~ /VRR=([\d\.,]+)/);
			if ($line[7] =~ /SAR=([\d\.]+)/){
				push @sar, $1;
			}
			my $bp = 0;
			$bp = $1 if ($line[7] =~ /BP=(\d+)/);
			my $gt = '';
			$gt = $1 if ($line[7] =~ /GT=(.+?);/);
			$gt = 'HT2' if (($len =~ /,/) or ($type =~ /,/));
			$strid = $1 if ($strid eq '') and ($line[7] =~ /TRID=(.+?);/);
			$strend = $1 if ($strend == 0) and ($line[7] =~ /TREND=(\d+)/);
			$strulen = $1 if ($strulen == 0) and ($line[7] =~ /TRULEN=(\d+)/);
			my $cn_gain = 0;
			my $cn_loss = 0;
			$cn_gain = $1 if ($line[7] =~ /gain\+([\d\.]+)/);
			$cn_gain = $1 if ($line[7] =~ /gain\+([A-Za-z]+\d+)/) and ($cn_gain == 0);
			$cn_loss = $1 if ($line[7] =~ /loss\-([\d\.]+)/);
			$gt{$id} = '0/1';
			$gt{$id} = '1/1' if ($gt eq 'HM');
			$gt{$id} = '2/2' if ($gt eq 'HT2');
			$gt{$id} = './.' if ($gt eq 'NA');
			$bp{$id} = $bp;
			$read{$id} = $read;
			$vrr{$id} = $vrr;
			if ($len =~ /,/){
				my ($len1, $len2) = split (/,/, $len);
				if ($type =~ /,/){
					my ($type1, $type2) = split (/,/, $type);
					$ctype{$type1} ++;
					$ctype{$type2} ++;
					if ($type1 eq 'DEL'){
						$cn{$id} = "-$cn_loss,$cn_gain";
						$len{$id} = "-$len1,$len2";
						push @dellen, $len1;
						push @inslen, $len2;
					}
					else{
						$cn{$id} = "$cn_gain,-$cn_loss";
						$len{$id} = "$len1,-$len2";
						push @inslen, $len1;
						push @dellen, $len2;
						$ins_ac ++;
					}
				}
				else{
					if ($type eq 'DEL'){
						if ($line[7] =~ /CN=loss\-([\d\.]+),loss\-([\d\.]+)/){
							$cn{$id} = "-$1,-$2";
							$len{$id} = "-$len1,-$len2";
						}
						elsif ($line[7] =~ /CN=loss\-([\d\.]+),([\d\.]+)/){
							$cn{$id} = "-$1,-$2";
							$len{$id} = "-$len1,-$len2";
						}
						push @dellen, $len1, $len2;
					}
					else{
						if ($line[7] =~ /CN=gain\+([\d\.]+),gain\+([\d\.]+)/){
							$cn{$id} = "$1,$2";
							$len{$id} = "$len1,$len2";
						}
						elsif ($line[7] =~ /CN=gain\+([A-Za-z]+\d+),gain\+([\d\.]+)/){
							$cn{$id} = "$1,$2";
							$len{$id} = "$len1,$len2";
						}
						elsif ($line[7] =~ /CN=gain\+([\d\.]+),gain\+([A-Za-z]+\d+)/){
							$cn{$id} = "$1,$2";
							$len{$id} = "$len1,$len2";
						}
						elsif (($line[7] =~ /CN=.+,gain\+[\d\.]+/) and ($line[7] =~ /MEILEN=/)){
							my $cn2 = $1 if ($line[7] =~ /CN=.+,gain\+([\d\.]+)/);
							my $mei = $1 if ($line[7] =~ /MEI=([^=]+);/);
							my $meilen = $1 if ($line[7] =~ /MEILEN=(\d+)/);
							$cn{$id} = "$mei$meilen,$cn2";
							$len{$id} = "$len1,$len2";
						}
						elsif (($line[7] =~ /CN=gain\+[\d\.]+,gain\+/) and ($line[7] =~ /MEILEN=/)){
							my $cn1 = $1 if ($line[7] =~ /CN=gain\+([\d\.]+),gain\+/);
							my $mei = $1 if ($line[7] =~ /MEI=([^=]+);/);
							my $meilen = $1 if ($line[7] =~ /MEILEN=(\d+)/);
							$cn{$id} = "$cn1,$mei$meilen";
							$len{$id} = "$len1,$len2";
						}
						else{
							$cn{$id} = $cn_gain;
							$len{$id} = "$len1,$len2";
						}
						push @inslen, $len1, $len2;
						$ins_ac += 2;
					}
					$ctype{$type} ++;
				}
			}
			else{
				if ($type eq 'DEL'){
					$cn{$id} = "-$cn_loss";
					$len{$id} = "-$len";
					push @dellen, $len;
				}
				else{
					$cn{$id} = $cn_gain;
					$len{$id} = $len;
					push @inslen, $len;
					$ins_ac ++;
					$ins_ac ++ if ($gt eq 'HM');
				}
				$ctype{$type} ++;
			}
		}
		my $ins_range = '';
		my $del_range = '';
		my $type2 = '';
		my $max_inslen = 0;
		my $min_inslen = 0;
		if (scalar keys %ctype == 2){
			$type2 = 'DEL,INS';
			next if (@dellen == 0);
			next if (@inslen == 0);
			@dellen = sort {$a <=> $b} @dellen;
			$del_range = "$dellen[0]-$dellen[-1]";
			$del_range = $dellen[0] if ($dellen[0] == $dellen[-1]);
			@inslen = sort {$a <=> $b} @inslen;
			$ins_range = "$inslen[0]-$inslen[-1]";
			$ins_range = $inslen[0] if ($inslen[0] == $inslen[-1]);
			$max_inslen = $inslen[-1];
			$min_inslen = $inslen[0];
		}
		else{
			foreach my $ctype (keys %ctype){
				$type2 = $ctype;
			}
			if ($type2 eq 'DEL'){
				next if (@dellen == 0);
				@dellen = sort {$a <=> $b} @dellen;
				$del_range = "$dellen[0]-$dellen[-1]";
				$del_range = $dellen[0] if ($dellen[0] == $dellen[-1]);
			}
			else{
				if ($ins_ac == 0){
					$delete_str3 ++;
					next;
				}
				@inslen = sort {$a <=> $b} @inslen;
				$ins_range = "$inslen[0]-$inslen[-1]";
				$ins_range = $inslen[0] if ($inslen[0] == $inslen[-1]);
				$max_inslen = $inslen[-1];
				$min_inslen = $inslen[0];
			}
		}
=pod
		if (exists $ins_str{$strid}){					# merge TR-INS and INS within the same strid when meeting the criteria
			foreach my $pos_info (sort {${$ins_str{$strid}}{$a} <=> ${$ins_str{$strid}}{$b}} keys %{$ins_str{$strid}}){
				my $ac2 = ${$ins_str{$strid}}{$pos_info};
				my ($chr2, $pos2, $ilen2) = split (/=/, $pos_info);
				if (($ins_ac >= $ac2 * 2) and ($max_inslen >= $ilen2)){
					if (exists ${${$vcf_cons2{$chr2}}{$pos2}}{'INS'}){
						my $iline = ${${$vcf_cons2{$chr2}}{$pos2}}{'INS'};
						my @iline = split (/\t/, $iline);
						my $count = 0;
						foreach (@iline){
							$count ++;
							next if ($count <= 9);
							next if ($_ =~ /^0\/0/);
							my $ID = $Gid_order{$count};
							my ($igt, $ipos, $ilen, $iread, $ivrr, $ibp, $idpr) = split (/:/, $_);
							my $icn = int ($ilen / $strulen * 0.5 * 10 + 0.5) / 10;
							if (!exists $gt{$ID}){
								$len{$ID} = $ilen;
								$gt{$ID} = $igt ;
								$read{$ID} = $iread;
								$vrr{$ID} = $ivrr;
								$bp{$ID} = $ibp;
								$cn{$ID} = $icn;
							}
						}
						delete ${${$vcf_cons2{$chr2}}{$pos2}}{'INS'};
						$delete_ins2 ++;
					}
				}
				elsif (($ins_ac > 0) and ($ins_ac * 4 < $ac2)){
					if (($ilen2 * 0.7 >= $min_inslen) and ($ilen2 * 1.4 <= $max_inslen)){
						if (exists ${${$vcf_cons2{$chr2}}{$pos2}}{'INS'}){
							my $iline = ${${$vcf_cons2{$chr2}}{$pos2}}{'INS'};
							my @iline = split (/\t/, $iline);
							my $count = 0;
							my $merge_flag = 0;
							foreach (@iline){
								$count ++;
								next if ($count <= 9);
								next if ($_ !~ /^0\/0/);
								my $ID = $Gid_order{$count};
								next if (!exists $gt{$ID}) or (!exists $read{$ID});
								my $str_len = $len{$ID};
								my $flag = 0;
								if ($str_len =~ /,/){
									my ($str_len1, $str_len2) = split (/,/, $str_len);
									$str_len = 0;
									if ($str_len1 > 0){
										if (($str_len1 / $ilen2 >= 0.5) and ($str_len1 / $ilen2 <= 2)){
											$str_len = $str_len1;
											$flag = 1;
										}
									}
									if (($str_len == 0) and ($str_len2 > 0)){
										if (($str_len2 / $ilen2 >= 0.5) and ($str_len2 / $ilen2 <= 2)){
											$str_len = $str_len2;
											$flag = 2;
										}
									}
								}
								else{
									if ($str_len > 0){
										if (($str_len / $ilen2 >= 0.5) and ($str_len / $ilen2 <= 2)){
											$flag = 1;
										}
									}
								}
								next if ($flag == 0);
								$merge_flag = 1;
								my $ilen = $str_len;
								my $igt = $gt{$ID};
								$igt = '0/1' if ($igt eq '2/2');
								my $iread = $read{$ID};
								if ($iread =~ /,/){
									my @iread = split (/,/, $iread);
									$iread = $iread[$flag - 1];
								}
								my $ivrr = $vrr{$ID};
								if ($ivrr =~ /,/){
									my @ivrr = split (/,/, $ivrr);
									$ivrr = $ivrr[$flag - 1];
								}
								my $ibp = $bp{$ID};
								my $ipos = $pos;
								my $dpr = 1;
								my $info = "$igt:$ipos:$ilen:$iread:$ivrr:$ibp:$dpr:0";
								$iline[$count - 1] = $info;
							}
							if ($merge_flag == 1){
								my $new_iline = join ("\t", @iline);
								${${$vcf_cons2{$chr2}}{$pos2}}{'INS'} = $new_iline;
								$delete_str2 ++;
								next;
							}
						}
					}
				}
			}
		}
=cut
		my $ave_sar = 0;
		if (@sar > 0){
			my $sum_sar = 0;
			map{$sum_sar += $_} @sar;
			$ave_sar = int ($sum_sar / @sar * 100 + 0.5) / 100;
		}
		my $format = 'GT:VL:VN:VR:CN:BP';

		my $format_info = '';
		my %gt_info;
		foreach my $id (keys %read){
			my $len2 = $len{$id};
			my $gt = $gt{$id};
			my $read = $read{$id};
			my $vrr = $vrr{$id};
			my $bp = $bp{$id};
			my $cn = $cn{$id};
			my $info = "$gt:$len2:$read:$vrr:$cn:$bp";
			$gt_info{$id} = $info;
		}
		foreach my $id (sort keys %Gid){
			if (exists $gt_info{$id}){
				$format_info .= "$gt_info{$id}\t";
			}
			else{
				$format_info .= "0/0:0:0:0:0:0\t";
			}
		}
		$format_info =~ s/\t$//;
		my $new_line = "$chr\t$pos\t.\t.\t<TR:CNV>\t.\tPASS\tSVTYPE=$type2";
		if ($del_range ne ''){
			$new_line .= ";DELLEN=$del_range";
		}
		if ($ins_range ne ''){
			$new_line .= ";INSLEN=$ins_range";
		}
		$new_line .= ";SAR=$ave_sar" if (@sar > 0);
		$new_line .= ";TRID=$strid;TREND=$strend;TRULEN=$strulen\t$format\t$format_info";
		${${$vcf_cons2{$chr02d}}{$pos}}{'CNV'} = $new_line;
		$total_str2 ++;
	}
}

print STDERR "Deleted INS2: $delete_ins2\n";
print STDERR "Deleted STR-INS2: $delete_str2 ($delete_str3)\n";
print STDERR "Initial STR-CNV: ", scalar keys %cnv_info, "\n";
print STDERR "Final STR-CNV: $total_str2\n";

print STDERR "4th step completed:\n";

my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
$year += 1900;
$mon = sprintf ("%02d", $mon);
$mday = sprintf ("%02d", $mday);
my $time = $year . $mon . $mday;

my $ID_str = '';

foreach my $gid (sort keys %Gid){
    $ID_str .= "$gid\t";
}
$ID_str =~ s/\t$//;

my $AN = $total_sample_num * 2;

system ("mkdir $out_dir") if ($out_dir ne '.') and (!-d $out_dir);

my $merge_out = "$out_dir/$out_prefix.All-samples.vcf";

open (OUT, "> $merge_out");

print OUT "##fileformat=VCFv4.0\n";
print OUT "##fileDate=$time\n";
print OUT "##source=MOPlineV1.0\n";
if ($non_human == 0){
	print OUT "##reference=GRCh37\n" if ($build eq '37');
	print OUT "##reference=GRCh38\n" if ($build eq '38');
	print OUT "##reference=CHM13\n" if ($build eq 'T2T');
}
else{
	my $ref = basename ($ref_index);
	$ref = $1 if ($ref =~ /(.+)\.fas*t*a*/);
	print OUT "##reference=$ref\n";
}
print OUT "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of tandem repeat expansion/contraction (TR-CNV) and structural variation (SV)\">\n";
print OUT "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles (0 when undefined)\">\n";
print OUT "##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Range of INS length observed in the TR (TR-CNV only)\">\n";
print OUT "##INFO=<ID=DELLEN,Number=1,Type=Integer,Description=\"Range of DEL length observed in the TR (TR-CNV only)\">\n";
print OUT "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variantdescribed in this record\">\n";
print OUT "##INFO=<ID=BPLEN,Number=1,Type=Integer,Description=\"Breakpoint distance of INS\">\n";
print OUT "##INFO=<ID=SAR,Number=1,Type=Float,Description=\"Mean ratio of SV-supporting reads with mapping quality 0 to total SV-supporting reads, including secondary alignments\">\n";
print OUT "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">\n";
print OUT "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count\">\n";
print OUT "##INFO=<ID=SC,Number=1,Type=Integer,Description=\"Number of samples carrying the SV allele\">\n";
print OUT "##INFO=<ID=DUPPOS,Number=1,Type=Integer,Description=\"Reference position of tandem duplication found between INS sequence and reference\">\n";
print OUT "##INFO=<ID=DUPLEN,Number=1,Type=Integer,Description=\"Length of tandem duplication between INS and reference\">\n";
print OUT "##INFO=<ID=TRAPOS,Number=1,Type=Integer,Description=\"Reference chr:pos of translocated segment\">\n";
print OUT "##INFO=<ID=TRALEN,Number=1,Type=Integer,Description=\"Length of translocated segment\">\n";
print OUT "##INFO=<ID=MEILEN,Number=1,Type=Integer,Description=\"Length of mobile element found in INS\">\n";
print OUT "##INFO=<ID=TRID,Number=1,Type=String,Description=\"TR ID containing TR-CNV\">\n";
print OUT "##INFO=<ID=TREND,Number=1,Type=Integer,Description=\"END position of TR-ID\">\n";
print OUT "##INFO=<ID=TRULEN,Number=1,Type=Integer,Description=\"Length of repeat unit of TR-ID\">\n";
print OUT "##INFO=<ID=TRUNIT,Number=1,Type=String,Description=\"motifsize:copy-number-range detected in INS sequence (only for SVTYPE=INS)\">\n";

print OUT "##FILTER=<ID=LowConf,Description=\"Repeat region where TR-CNVs/SVs could be unreliably called\">\n";
print OUT "##FILTER=<ID=LowQual,Description=\"Low quality TR-CNVs/SVs with > $max_SAR SAR\">\n";

print OUT "##ALT=<ID=TR:CNV,Description=\"INS(repeat expansion) and/or DEL (repeat contraction) within Tandem Repeat (TR) regions\">\n";
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

print OUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print OUT "##FORMAT=<ID=VP,Number=1,Type=Integer,Description=\"Position of SV for a sample (only non-TR-CNVs)\">\n";
print OUT "##FORMAT=<ID=VL,Number=.,Type=Integer,Description=\"Length of TR-CNV/SV for a sample (length of included INS for DUP)\">\n";
print OUT "##FORMAT=<ID=VN,Number=.,Type=Integer,Description=\"Number of reads supporting TR-CNV/SV allele for a sample\">\n";
print OUT "##FORMAT=<ID=VR,Number=.,Type=Float,Description=\"Ratio of TR-CNV/SV-supporting reads to read depth at the site for a sample\">\n";
print OUT "##FORMAT=<ID=CN,Number=.,Type=Float,Description=\"Copy number of CNV (increased or decreased number of a repeat unit) in TR for a sample (TR-CNV only)\">\n";
print OUT "##FORMAT=<ID=BP,Number=1,Type=Float,Description=\"Number of break ends (soft-clipped read ends) supporting TR-CNV/SV for a sample\">\n";
print OUT "##FORMAT=<ID=DR,Number=1,Type=Float,Description=\"Ratio of read depth in DEL/DUP region to that to the flanking regions for a sample (only non-TR-CNV)\">\n";
print OUT "##FORMAT=<ID=TR,Number=1,Type=String,Description=\"Tandem repeat (unit length-copy number) observed in non-TR INS sequence\">\n";

foreach my $chr (@ref){
	print OUT "##contig=<ID=$chr,length=$ref_len{$chr}>\n";
}

print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$ID_str\n";

foreach my $chr (sort keys %vcf_cons2){
    my $chr2 = $chr;
    $chr2 =~ s/^0*//;
    foreach my $pos (sort {$a <=> $b} keys %{$vcf_cons2{$chr}}){
        foreach my $type (sort keys %{${$vcf_cons2{$chr}}{$pos}}){
            my @line = split (/\t/, ${${$vcf_cons2{$chr}}{$pos}}{$type});
            my $ac = 0;
            my $sc = 0;
            my $count = 0;
            my $single_inf = '';
            foreach (@line){
            	$count ++;
            	next if ($count <= 9);
            	my ($gt) = split (/:/, $_);
            	next if ($gt eq '0/0');
            	$sc ++;
            	my $ID = $Gid_order{$count};
            	if (exists $lowQ_sample{$ID}){
            		$single_inf = $_;
            	}
            	next if ($gt eq './.');
            	if ($gt eq '0/1'){
            		$ac ++;
            	}
            	else{
            		$ac += 2;
            	}
            }
            my $af = 0;
            if ($AN < 100){
            	$af = int ($ac / $AN * 100 + 0.5) / 100;
            }
            elsif ($AN < 1000){
            	$af = int ($ac / $AN * 1000 + 0.5) / 1000;
            }
            elsif ($AN < 10000){
            	$af = int ($ac / $AN * 10000 + 0.5) / 10000;
            }
            else{
            	$af = int ($ac / $AN * 100000 + 0.5) / 100000;
            }
            if (($sc == 1) and ($single_inf ne '') and ($line[4] =~ /DEL|INS/)){
            	my ($gt, $vp, $vl, $vn, $vr) = split (/:/, $single_inf);
            	if (($vl < 200) and ($vr < 0.2)){
            		next;
            	}
            }
		    $line[7] .= ";AF=$af;AC=$ac;SC=$sc";
		    $line[2] =~ s/^0// if ($line[0] =~ /^0/) and ($non_human == 1);
		    if ($type eq 'DEL'){
		    	my $len = $1 if ($line[7] =~ /SVLEN=-*(\d+)/);
		    	my $end = $1 if ($line[7] =~ /END=(\d+)/);
		    	my $end2 = $pos + $len - 1;
		    	if ($end != $end2){
		    		$line[7] =~ s/END=\d+/END=$end2/;
		    	}
		    }
		    if ($line[4] eq '<INS:BP>'){
		    	my $count = 0;
		    	foreach (@line){
		    		$count ++;
		    		next if ($count <= 9);
		    		next if ($_ =~ /^0\/0/);
		    		my @info = split (/:/, $_);
		    		if ($info[2] != 0){
		    			$info[2] = 0;
		    			$line[$count - 1] = join (':', @info);
		    		}
		    	}
		    }
            my $line = join ("\t", @line);
            print OUT $line, "\n";
        }
    }
}
close (OUT);


sub cons_len{
    my ($ref_len, $len) = @_;
    my $num = scalar @{$ref_len};
    my %len;
    my $len_sd = 3;
    if ($len > 300){
		$len_sd = int ($len * 0.01);
    }
    map{$len{$_} ++} @{$ref_len};
    my $pre_len = 0;
    foreach my $slen (sort {$a <=> $b} keys %len){
		if (($pre_len > 0) and ($slen - $pre_len <= $len_sd)){
		    my $num1 = $len{$pre_len};
		    my $num2 = $len{$slen};
		    if ($num2 > $num1){
				$len{$slen} += $num1;
		    }
		    elsif ($num2 < $num1){
				$len{$pre_len} += $num2;
		    }
		    else{
				$len{$slen} += $num1;
				$len{$pre_len} += $num2;
		    }
		}
		$pre_len = $slen;
    }
    my $high_freq_len = 0;
    my $high_freq_len2 = 0;
    my $high_freq = 0;
    my $high_freq2 = 0;
    my $new_len = 0;
    foreach my $slen (sort {$len{$b} <=> $len{$a}} keys %len){
		next if ($len{$slen} == 1);
		if ($high_freq_len == 0){
		    $high_freq_len = $slen;
		    $high_freq = $len{$slen};
		    next;
		}
		if ($high_freq_len2 == 0){
		    $high_freq_len2 = $slen;
		    $high_freq2 = $len{$slen};
		    last;
		}
    }
    if (($high_freq > 1) and ($high_freq >= $num * 0.3)){
		if ($high_freq == $high_freq2){
		    $new_len = int (($high_freq_len + $high_freq_len2) * 0.5 + 0.5);
		}
		else{
		    $new_len = $high_freq_len;
		}
    }
    else{
		my $half = int ($num * 0.5);
		my $half2 = $half - 1;
		if ($num % 2 == 0){
		    $new_len = int ((${$ref_len}[$half2] + ${$ref_len}[$half]) * 0.5 + 0.5);
		}
		else{
		    $new_len = ${$ref_len}[$half];
		}
    }
    return ($new_len);
}
