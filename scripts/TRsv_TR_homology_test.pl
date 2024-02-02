#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;

my $str_prox_file = shift @ARGV;

my $target_chr = shift @ARGV;

my $ref_fasta = shift @ARGV;

my $vcf_list = shift @ARGV;

my $out_prefix = shift @ARGV;

my $temp_dir = shift @ARGV;

my $out_file = shift @ARGV;

my $min_sv_len = shift @ARGV;

my %ref;
my %sample_id;
my $header = '';
my $seq = '';
my $Mb_count = 0;
my $Mbin_size = 1000000;

open (FILE, $ref_fasta) or die "$ref_fasta is not found: $!\n";
while (my $line = <FILE>){
    chomp $line;
    if ($line =~ /^>(\S+)/){
        if (($seq ne '') and ($header eq $target_chr)){
            while (length $seq > 0){
                if (length $seq >= $Mbin_size){
                    my $seq_Mb = substr ($seq, 0, $Mbin_size, '');
                    $ref{$Mb_count} = $seq_Mb;
                }
                else{
                    $ref{$Mb_count} = $seq;
                    $seq = '';
                }
                $Mb_count ++;
            }
        }
        $seq = '';
        $Mb_count = 0;
        $header = $1;
        last if ($header !~ /^[\dXY]+$/);
    }
    else{
        $seq .= uc $line;
    }
}
if (($seq ne '') and ($header eq $target_chr)){
    while (length $seq > 0){
        if (length $seq >= $Mbin_size){
            my $seq_Mb = substr ($seq, 0, $Mbin_size, '');
            $ref{$Mb_count} = $seq_Mb;
        }
        else{
            $ref{$Mb_count} = $seq;
            $seq = '';
        }
        $Mb_count ++;
    }
}
close (FILE);

open (FILE, $vcf_list) or die "$vcf_list is not found: $!\n";
while (my $line = <FILE>){
	chomp $line;
	next if ($line =~ /^#|^$/);
	my ($ID, $vcf) = split (/\t/, $line);
	$sample_id{$ID} = $vcf;
}
close (FILE);

open (OUT, "> $out_file");
open (FILE, $str_prox_file) or die "$str_prox_file is not found: $!\n";
while (my $line = <FILE>){
	chomp $line;
	my @line = split (/\t/, $line);
	my $chr = $line[0];
	next if ($chr ne $target_chr);
	my $pos = $line[1];
	my $type = $line[2];
	my $len = $line[3];
	my $hit_spos = $line[4];
	if ($type ne 'INS'){
		print OUT "$line\n";
	}
	else{
		if ($len < 50){
			print OUT "$line\n";
			next;
		}
		my $send = $line[5];
		my $ulen = $line[6];
		my $id_str = $line[7];
		my @ids = split (/,/, $id_str);
		my %ids;
		my @ins_line = ();
		my $ins_seq = '';
		my $str_seq = '';
		foreach my $idlen (@ids){
			my ($id, $vp, $vl) = split (/=/, $idlen);
			$ids{$id} = "$vp=$vl";
		}
		foreach my $id (keys %ids){
			my ($vp, $vl) = split (/=/, $ids{$id});
			next if ($vl / $len > 1.5) or ($vl / $len < 0.67);
			my $vcf = $sample_id{$id};
			my $vcf_base = $1 if ($vcf =~ /(.+)\.vcf$/);
			$vcf_base = $1 if ($vcf =~ /(.+)\.filt\.vcf$/);
			my $ins_fasta = "$vcf_base.INS.fa";
			next if (!-f $ins_fasta);
			@ins_line = `grep -A 1 $chr:$vp $ins_fasta`;
			if (@ins_line == 0){
				next;
			}
			else{
				last;
			}
		}
		if (@ins_line == 0){
			next;
		}
		foreach (@ins_line){
			chomp $_;
			next if ($_ =~ /^>/);
			$ins_seq = $_;
		}
print STDERR "No NS seq: $chr:$pos $hit_spos\n" if ($ins_seq eq '');
	    my $Mbin = int ($hit_spos / $Mbin_size);
	    my $Mbin_res = $hit_spos % $Mbin_size;
	    my $Mbin_end = int ($send / $Mbin_size);
	    my $str_len = $send - $hit_spos + 1;
	    if ($Mbin == $Mbin_end){
	        $str_seq = substr ($ref{$Mbin}, $Mbin_res - 1, $str_len);
	    }
	    else{
	        my $Mbin_endres = $send % $Mbin_size;
	        my $str_seq1 = substr ($ref{$Mbin}, $Mbin_res - 1);
	        my $str_seq2 = substr ($ref{$Mbin_end}, 0, $Mbin_endres);
	        $str_seq = $str_seq1 . $str_seq2;
	    }

		my ($identity, $ins_cov, $cov_len) = &blast_ins ($chr, $str_seq, $ins_seq);
		my $flag1 = 0;
		if (($identity >= 85) and ($ins_cov >= 50) and ($cov_len >= $ulen)){
			$flag1 = 1;
		}
		elsif (($identity >= 70) and ($ins_cov >= 70) and ($cov_len >= $ulen)){
			$flag1 = 1;
		}
		elsif (($len < 200) and ($identity + $ins_cov >= 130) and ($cov_len >= $ulen)){
            $flag1 = 1;
        }
		if ($flag1 == 1){
			print OUT "$line\n";
		}
		else{
			my ($identity2, $ins_cov2, $cov_len2) = &multalin_ins ($chr, $str_seq, $ins_seq, $pos);
			my $flag2 = 0;
			if (($identity2 >= 60) and ($ins_cov2 >= 70)){
				$flag2 = 1;
			}
			elsif (($identity2 >= 70) and ($ins_cov2 >= 60)){
				$flag2 = 1;
			}
			if (($identity2 >= 60) and ($ins_cov2 >= 70) and ($cov_len2 >= $ulen)){
                $flag2 = 1;
            }
            elsif (($identity2 >= 70) and ($ins_cov2 >= 40) and ($cov_len2 >= $min_sv_len) and ($cov_len2 >= $ulen)){
                $flag2 = 1;
            }
            elsif (($len < 200) and ($identity2 + $ins_cov2 >= 130) and ($cov_len2 >= $ulen)){
                $flag2 = 1;
            }
			if ($flag2 == 1){
				print OUT "$line\n";
			}
			else{
				my $rc_ins_seq = reverse $ins_seq;
				$rc_ins_seq =~ tr/ACGT/TGCA/;
				my ($identity3, $ins_cov3, $cov_len3) = &multalin_ins ($chr, $str_seq, $rc_ins_seq, $pos);
				my $flag2 = 0;
				if (($identity3 >= 60) and ($ins_cov3 >= 70)){
					$flag2 = 1;
				}
				elsif (($identity3 >= 70) and ($ins_cov3 >= 60)){
					$flag2 = 1;
				}
				if (($identity3 >= 60) and ($ins_cov3 >= 70) and ($cov_len3 >= $ulen)){
	                $flag2 = 1;
	            }
	            elsif (($identity3 >= 70) and ($ins_cov3 >= 40) and ($cov_len3 >= $min_sv_len) and ($cov_len3 >= $ulen)){
	                $flag2 = 1;
	            }
	            elsif (($len < 200) and ($identity3 + $ins_cov3 >= 130) and ($cov_len3 >= $ulen)){
	                $flag2 = 1;
	            }
				if ($flag2 == 1){
					print OUT "$line\n";
				}
				else{
					print STDERR "Low-homology $chr:$pos INS ($identity2% ident $ins_cov2% cov)\n" if ($identity2 + $ins_cov2 >= $identity3 + $ins_cov3);
					print STDERR "Low-homology $chr:$pos INS ($identity3% ident $ins_cov3% cov)\n" if ($identity2 + $ins_cov2 < $identity3 + $ins_cov3);
					next;
				}
			}
		}
	}
}
close (FILE);
close (OUT);


sub blast_ins{
	my ($chr, $str_seq, $ins_seq) = @_;
	$ins_seq =~ s/-//g if ($ins_seq =~ /-/);
	my $ins_len = length $ins_seq;
	my $str_len = length $str_seq;
    
    my $target_fasta = "$temp_dir/$out_prefix.chr$chr.target.fasta";
    my $query_fasta = "$temp_dir/$out_prefix.chr$chr.query.fasta";
    my $out_blast = "$temp_dir/$out_prefix.chr$chr.blast.out";
    open (OUT1, "> $target_fasta");
    print OUT1 ">STR\n";
    print OUT1 "$str_seq\n";
    close (OUT1);

    open (OUT2, "> $query_fasta");
    print OUT2 ">INS\n";
    print OUT2 "$ins_seq\n";
    close (OUT2);

    system ("makeblastdb -in $target_fasta -dbtype nucl -logfile $temp_dir/blastdb.chr$chr.log");
    system ("blastn -db $target_fasta -query $query_fasta -out $out_blast -task blastn-short -outfmt 6");

    if (-z $out_blast){
#    	print STDERR "No blast result: $chr: $str_seq\t$ins_seq\n";
    	system ("rm $target_fasta $query_fasta $out_blast");
        return (0, 0, 0);
    }

    my %align;
    my %ident;
    open (NUCM, $out_blast) or die "$out_blast is not found: $!\n";
    while (my $aline = <NUCM>){
        chomp $aline;
        my @aline = split (/\t/, $aline);
        my $qstart = $aline[6];
        my $qend = $aline[7];
        my $ident = $aline[2];
        if (!exists $align{$qstart}){
            $align{$qstart} = $qend;
            $ident{$qstart} = $ident;
        }
        else{
            if ($qend > $align{$qstart}){
                $align{$qstart} = $qend;
                $ident{$qstart} = $ident;
            }
        }
        
    }
    close (NUCM);

    my %align_ident;
    my $total_cov = 0;
    my $ave_ident = 0;
    my $pre_pos = 0;
    my $pre_end = 0;
    my $last_pos = 0;

    foreach my $pos (sort {$a <=> $b} keys %align){
        my $end = $align{$pos};
        my $ident = $ident{$pos};
        if ($pre_end == 0){
            $align_ident{"$pos-$end"} = $ident;
            $pre_pos = $pos;
            $pre_end = $end;
            $last_pos = $pos;
            next;
        }
        $last_pos = $pos;
    }
    my $flag = 0;
    while ($flag == 0){
        my %ovl_align;
        my $pos2 = 0;
        my $end2 = 0;
        $flag = 1 if (scalar keys %align == 0);
        foreach my $pos (sort {$a <=> $b} keys %align){
            $flag = 1 if ($pos == $last_pos);
            my $end = $align{$pos};
            next if ($end <= $pre_end);
            my $ident = $ident{$pos};
            if ($pos <= $pre_end){
                $ovl_align{$pos} = $end;
            }
            else{
                $pos2 = $pos;
                $end2 = $end;
                last;
            }
        }
        if (scalar keys %ovl_align > 0){
            my $top_pos = 0;
            my $top_end = 0;
            foreach my $apos (sort {$ovl_align{$b} <=> $ovl_align{$a}} keys %ovl_align){
                $top_pos = $apos;
                $top_end = $ovl_align{$apos};
                last;
            }
            my $pos2 = $pre_end + 1;
            my $ident2 = $ident{$top_pos};
            $align_ident{"$pos2-$top_end"} = $ident2;
            $pre_pos = $top_pos;
            $pre_end = $top_end;
        }
        else{
            my $ident = $ident{$pos2} if (exists $ident{$pos2});
			  $ident = $ident{$last_pos} if (!exists $ident{$pos2});
            $align_ident{"$pos2-$end2"} = $ident;
            $pre_pos = $pos2;
            $pre_end = $end2;
        }
    }
    foreach my $range (keys %align_ident){
        my ($pos, $end) = split (/-/, $range);
        $total_cov += $end - $pos + 1;
    }
    foreach my $range (keys %align_ident){
        my ($pos, $end) = split (/-/, $range);
        my $len = $end - $pos + 1;
        my $ident = $align_ident{$range};
        $ave_ident += int (($len / $total_cov) * $ident * 10 + 0.5) / 10;
    }
    my $cov_rate = int ($total_cov / $ins_len * 1000 + 0.5) / 10;
    system ("rm $target_fasta $query_fasta $out_blast");

    return ($ave_ident, $cov_rate, $total_cov);
}

sub multalin_ins{
	my ($chr, $str_seq, $ins_seq, $pos) = @_;
	$ins_seq =~ s/-//g if ($ins_seq =~ /-/);
	my $ins_len = length $ins_seq;
	my $str_len = length $str_seq;
    my $input_fasta = "$temp_dir/$out_prefix.chr$chr.fasta";
    my $input_base = "$temp_dir/$out_prefix.chr$chr";
    my $multalin_out = "$input_base.msf";
    my $multalin_out1 = "$input_base.clu";
    my $multalin_out2 = "$input_base.cl2";
    my @ins_seq;
    if ($ins_len > $str_len + 30){
    	if ($str_len >= 50){
    		my $bin = int ($ins_len / 50) + 1;
	    	my $res = $ins_len % $bin;
	    	while (1){
	    		if (length ($ins_seq) >= 50 + $res){
	    			my $subseq = substr ($ins_seq, 0, 50, '');
	    			push @ins_seq, $subseq;
	    		}
	    		else{
	    			push @ins_seq, $ins_seq;
	    			last;
	    		}
	    	}
    	}
    	else{
	    	my $bin = int ($ins_len / $str_len) + 1;
	    	my $sublen = int ($ins_len / $bin);
	    	my $res = $ins_len % $bin;
	    	while (1){
	    		if (length ($ins_seq) >= $sublen + $res){
	    			my $subseq = substr ($ins_seq, 0, $sublen, '');
	    			push @ins_seq, $subseq;
	    		}
	    		else{
	    			push @ins_seq, $ins_seq;
	    			last;
	    		}
	    	}
	    }
    }
    else{
    	push @ins_seq, $ins_seq;
    }
    my $match_len = 0;
    my $mismatch_len = 0;
    my $match_len2 = 0;
    my $mismatch_len2 = 0;
    my $cov_len = 0;
    foreach my $iseq (@ins_seq){
    	open (OUT1, "> $input_fasta");
	    print OUT1 ">STR\n";
	    print OUT1 "$str_seq\n";
	    print OUT1 ">INS\n";
	    print OUT1 "$iseq\n";
	    close (OUT1);
	    system ("multalin -q $input_fasta > $temp_dir/multalin.chr$chr.log");
	    
	    if (!-f $multalin_out){
#	        print STDERR "$chr:$pos INS multalin 2nd alignment failed:\n";
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
	            $qalign .= $subseq;
	        }
	        elsif ($line =~ /^\s+STR\S+\s+(.+)/){
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
	    system ("rm -f $input_fasta $multalin_out $multalin_out1 $multalin_out2");
	    $match_len += $match;
	    $mismatch_len += $mismatch;
	    my $match_rate = int ($match / ($match + $mismatch) * 1000 + 0.5) / 10;
	    if ($match_rate >= 60){
	        $cov_len += $cov;
	        $match_len2 += $match;
	        $mismatch_len2 += $mismatch;
	    }
	    elsif ($match_rate >= 40){
	    	$cov_len += $match;
	    }
    }
    my $identity = int ($match_len / ($match_len + $mismatch_len) * 1000 + 0.5) / 10;
    my $identity2 = 0;
	$identity2 = int ($match_len2 / ($match_len2 + $mismatch_len2) * 1000 + 0.5) / 10 if ($match_len2 > 0);
	$identity = $identity2 if ($identity < $identity2);
	my $coverage = int ($cov_len / $ins_len * 1000 + 0.5) / 10;
	
    return ($identity, $coverage, $cov_len);
}

