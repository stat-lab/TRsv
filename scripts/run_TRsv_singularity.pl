#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;

my $sif_file = '';
my $TRsv_command = '';
my $no_home = 0;
my $help;

my @args;
my @argv2;
my $sif_flag = 0;
my $com_flag = 0;
foreach (@ARGV){
	if ($_ =~ /-sif/){
		$sif_flag = 1;
		push @argv2, $_;
		next;
	}
	elsif ($_ =~ /-com/){
		$com_flag = 1;
		push @argv2, $_;
		next;
	}
	elsif ($sif_flag == 1){
		$sif_flag = 0;
		push @argv2, $_;
		next;
	}
	elsif ($com_flag == 1){
		$com_flag = 0;
		push @argv2, $_;
		next;
	}
	else{
		push @args, $_;
	}
}
@ARGV = (@argv2);
	
GetOptions(
	'sif=s' => \$sif_file,
	'com=s' => \$TRsv_command,
  'no_home|noh' => \$no_home,
  'help' => \$help
) or pod2usage(-verbose => 0);
pod2usage(-verbose => 0) if $help;

=head1 SYNOPSIS

  run_TRsv_singularity.pl --sif <TRsv.sif> --com <command> [options for specified command: specified files should be absolute path]
  (example: run_TRsv_singularity.pl --sif /home/tools/TRsv_v1.1/TRsv.sif --com call -x hifi -b /home/TR_call/sample1/sample1.bam -r /home/Ref/GRCh38.fa --build 38 -ml 50 -n 6 -p sample1.hifi.ml50)

  Options:
   -sif <STR>      absolute path of TRsv.sif file [mandatory]
   -com <STR>		TRsv command (any of call, joint_call, or annotate) [mandatory]
   -no_home        <BOOLEAN>  do not add $HOME on the host to singularity container [default: false]
   -help or -h     output help message
   
=cut

die "sif file is not specified: or does not exist\n" if ($sif_file eq '') or (!-f $sif_file);
die "TRsv command (call|joint_call|annotate) is not specified:\n" if ($TRsv_command ne 'call') and ($TRsv_command ne 'joint_call') and ($TRsv_command ne 'annotate');

my $args = join (' ', @args);

print STDERR "$args\n";
	
my %path;

my $work_dir = `pwd`;
chomp $work_dir;

my $TRsv_data_dir = '/opt/local/tools/TRsv/Data';

#$path{$TRsv_data_dir} = 1;

my $out_prefix = '';
if ($args =~ /-prefix\s+(\S+)/){
	$out_prefix = $1;
}
elsif ($args =~ /-p\s+(\S+)/){
	$out_prefix = $1;
}
print STDERR "prefix: $out_prefix\n";

my $bind_dir = '';

if ($TRsv_command eq 'call'){
	my $temp_dir = "$work_dir/$out_prefix.temp";
	system ("mkdir $temp_dir") if (!-d $temp_dir);
	my $platform = '';
	my $bam_file = '';
	my $ref_file = '';
	my $conf_file = '';
	my $repeat_bed = '';
	my $repeat_u = '';
	my $lowconf_list = '';
	my $gap_bed = '';
	my $exclude_bed = '';
	my $te_fasta = '';
	if ($args =~ /-platform\s+(\S+)/){
		$platform = $1;
	}
	elsif ($args =~ /-x\s+(\S+)/){
		$platform = $1;
	}
	if ($args =~ /-bam_file\s+(\S+)/){
		$bam_file = $1;
	}
	elsif ($args =~ /-b\s+(\S+)/){
		$bam_file = $1;
	}
	if ($args =~ /-ref_file\s+(\S+)/){
		$ref_file = $1;
	}
	elsif ($args =~ /-r\s+(\S+)/){
		$ref_file = $1;
	}
	if ($args =~ /-conf\s+(\S+)/){
		$conf_file = $1;
	}
	elsif ($args =~ /-c\s+(\S+)/){
		$conf_file = $1;
	}
	if ($args =~ /-repeat_bed\s+(\S+)/){
		$repeat_bed = $1;
	}
	elsif ($args =~ /-rep\s+(\S+)/){
		$repeat_bed = $1;
	}
	if ($args =~ /-repeat_u\s+(\S+)/){
		$repeat_u = $1;
	}
	elsif ($args =~ /-reu\s+(\S+)/){
		$repeat_u = $1;
	}
	if ($args =~ /-lowconf_tr\s+(\S+)/){
		$lowconf_list = $1;
	}
	elsif ($args =~ /-lcs\s+(\S+)/){
		$lowconf_list = $1;
	}
	if ($args =~ /-exclude_bed\s+(\S+)/){
		$exclude_bed = $1;
	}
	elsif ($args =~ /-exc\s+(\S+)/){
		$exclude_bed = $1;
	}
	if ($args =~ /-gap_bed\s+(\S+)/){
		$gap_bed = $1;
	}
	elsif ($args =~ /-gb\s+(\S+)/){
		$gap_bed = $1;
	}
	if ($args =~ /-te_fasta\s+(\S+)/){
		$te_fasta = $1;
	}
	elsif ($args =~ /-tf\s+(\S+)/){
		$te_fasta = $1;
	}
	my $bam_dir = dirname ($bam_file);
	my $ref_dir = dirname ($ref_file);
	my $conf_dir = '';
	my $rep_dir = '';
	my $rep_u_dir = '';
	my $exclude_dir = '';
	my $lowconf_dir = '';
	my $gap_dir = '';
	my $te_dir = '';
	if ($conf_file ne ''){
		$conf_dir = dirname ($conf_file);
	}
	if ($repeat_bed ne ''){
		$rep_dir = dirname ($repeat_bed);
	}
	if ($repeat_u ne ''){
		$rep_u_dir = dirname ($repeat_u);
	}
	if ($exclude_bed ne ''){
		$exclude_dir = dirname ($exclude_bed);
	}
	if ($lowconf_list ne ''){
		$lowconf_dir = basenmae ($lowconf_list);
	}
	if ($gap_bed ne ''){
		$gap_dir = dirname ($gap_bed);
	}
	if ($te_fasta ne ''){
		$te_dir = dirname ($te_fasta);
	}
	$path{$temp_dir} = 1;
	$path{$bam_dir} = 1;
	$path{$ref_dir} = 1;
	$path{$conf_dir} = 1 if ($conf_dir ne '');
	$path{$rep_dir} = 1 if ($rep_dir ne '');
	$path{$rep_u_dir} = 1 if ($rep_u_dir ne '');
	$path{$exclude_dir} = 1 if ($exclude_dir ne '');
	$path{$lowconf_dir} = 1 if ($lowconf_dir  ne '');
	$path{$gap_dir} = 1 if ($gap_dir  ne '');
	$path{$te_dir} = 1 if ($te_dir  ne '');
	if ($platform eq 'ont'){
		my $ML_dir = "$TRsv_data_dir/ML_train/ONT";
		$path{$ML_dir} = 1;
	}
	elsif ($platform eq 'clr'){
		my $ML_dir = "$TRsv_data_dir/ML_train/CLR";
		$path{$ML_dir} = 1;
	}
}
elsif ($TRsv_command eq 'joint_call'){
	my $vcf_list = '';
	my $out_dir = '';
	my $gap_bed = '';
	my $ref_index = '';
	my $gender_list = '';
	my $lowQ_sample_list = '';
	if ($args =~ /-vcf_list\s+(\S+)/){
		$vcf_list = $1;
	}
	elsif ($args =~ /-v\s+(\S+)/){
		$vcf_list = $1;
	}
	if ($args =~ /-out_dir\s+(\S+)/){
		$out_dir = $1;
	}
	elsif ($args =~ /-od\s+(\S+)/){
		$out_dir = $1;
	}
	if ($args =~ /-gap_bed\s+(\S+)/){
		$gap_bed = $1;
	}
	elsif ($args =~ /-g\s+(\S+)/){
		$gap_bed = $1;
	}
	if ($args =~ /-ref\s+(\S+)/){
		$ref_index = $1;
	}
	elsif ($args =~ /-r\s+(\S+)/){
		$ref_index = $1;
	}
	if ($args =~ /-gender\s+(\S+)/){
		$gender_list = $1;
	}
	elsif ($args =~ /-g\s+(\S+)/){
		$gender_list = $1;
	}
	if ($args =~ /-lqs_list\s+(\S+)/){
		$lowQ_sample_list = $1;
	}
	elsif ($args =~ /-lqs\s+(\S+)/){
		$lowQ_sample_list = $1;
	}
	my $vcf_dir = dirname ($vcf_list);
	my $gap_dir = '';
	my $ref_dir = '';
	my $gender_dir = '';
	my $lowQ_dir = '';
	$out_dir = '' if ($out_dir eq '.');
	if ($gap_bed ne ''){
		$gap_dir = dirname ($gap_bed);
	}
	if ($ref_index ne ''){
		$ref_dir = dirname ($ref_index);
	}
	if ($gender_list ne ''){
		$gender_dir = basenmae ($gender_list);
	}
	if ($lowQ_sample_list ne ''){
		$lowQ_dir = dirname ($lowQ_sample_list);
	}
	$path{$vcf_dir} = 1;
	$path{$out_dir} = 1 if ($out_dir ne '');
	$path{$gap_dir} = 1 if ($gap_dir ne '');
	$path{$ref_dir} = 1 if ($ref_dir ne '');
	$path{$gender_dir} = 1 if ($gender_dir ne '');
	$path{$lowQ_dir} = 1 if ($lowQ_dir ne '');
}
elsif ($TRsv_command eq 'annotate'){
	my $vcf_file = '';
	my $ref_gff = '';
	if ($args =~ /-vcf\s+(\S+)/){
		$vcf_file = $1;
	}
	elsif ($args =~ /-v\s+(\S+)/){
		$vcf_file = $1;
	}
	if ($args =~ /-ref\s+(\S+)/){
		$ref_gff = $1;
	}
	elsif ($args =~ /-r\s+(\S+)/){
		$ref_gff = $1;
	}
	my $vcf_dir = dirname ($vcf_file);
	my $ref_dir = '';
	$ref_dir = dirname ($ref_gff) if ($ref_gff ne '');
	$path{$vcf_dir} = 1;
	$path{$ref_dir} = 1 if ($ref_dir ne '');
}

foreach my $dir (keys %path){
	$bind_dir .= "$dir,";
}
$bind_dir =~ s/,$//;

my $command = "singularity exec --bind $bind_dir $sif_file TRsv $TRsv_command $args";
$command = "singularity exec --bind $bind_dir --no-home $sif_file TRsv $TRsv_command $args" if ($no_home == 1);

my $command_log = "$out_prefix.command.log";
my $error_log = "$out_prefix.error.log";

open (OUT, "> $command_log");
print OUT "singularity TRsv command: $command\n";
close (OUT);

my $run = `$command 2>$error_log`;

