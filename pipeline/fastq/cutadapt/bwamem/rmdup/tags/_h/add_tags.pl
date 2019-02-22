#!/usr/bin/env perl
$^W=1;

use strict;
use Gapafim;
use Bio::DB::Sam;
use Getopt::Long;

my $fai;

# We add the custom flags below to the bam file:
# YA: score of semiglobal alignment between the first $prefix_length bp of read2 and $consensus.
# YG: YG: score of prefix vs flanking genome.
# YR: a 0 or 1 valued integer, indicating if the pair is a reference read pair or not. 1: the beginning of read 2 contains less than $soft_clip_length_threshold soft clip letters, and it is in proper pair with read 1. 0: otherwise.


sub max
{
    my ($a, $b)=@_;
    return $a > $b ? $a : $b;
}


sub min
{
    my ($a, $b)=@_;
    return $a < $b ? $a : $b;
}


sub rc
{
    # Reverse complement
    my $s=shift;
    $s=reverse $s;
    $s=~ tr/ACTGactg/TGACtgac/;
    return $s;
}


sub cigar2reflen
{
    #calend
    my $cigar=shift;
    my $r=0;
    while($cigar=~/(\d+)[MDN]/g)
    {
	#Match, Del, or N (reference skip, for spliced reads)
	$r+=$1;
    }
    return $r;
}


sub calculate_yr
{
    my ($r2, $soft_clip_length_threshold)= @_;
    
    return 0 if ($r2->[1] & 4);

    return 0 if ($r2->[1] & 2) == 0;
    
    return 0 if (($r2->[1] & 16) == 0) && ($r2->[5]=~/^(\d+)S/) && ($1>=$soft_clip_length_threshold);

    return 0 if ($r2->[1] & 16) && ($r2->[5]=~/(\d+)S$/) && ($1>=$soft_clip_length_threshold);

    return 1;
}

sub calculate_ys
{
    my ($r2)= @_;
    
    my $s=0;
    if ((($r2->[1] & 16) == 0) && ($r2->[5]=~/^(\d+)S/))
    {
	$s=$1;
    }
    elsif (($r2->[1] & 16) && ($r2->[5]=~/(\d+)S$/))
    {
	$s=$1;
    }
    return $s;
}

sub calculate_yg_from_r2
{
    # Calculates YG from a window of 2*$extra_flank centered on the beginning of R2 - prefix_length

    my ($r2, $prefix, $prefix_length, $extra_flank)=@_;

    my $flank;
    if($r2->[1] & 16)
    {
	my $e = $r2->[3] + cigar2reflen($r2->[5])-1+$prefix_length+$extra_flank;
	my $b = max(0, $e-2*$extra_flank-$prefix_length+1);
	$flank = rc($fai->fetch("$r2->[2]:$b-$e"));
    }
    else
    {
	my $b = max(0, $r2->[3]-$prefix_length-$extra_flank);
	my $e = $r2->[3]-1+$extra_flank;
	$flank = $fai->fetch("$r2->[2]:$b-$e");
    }
    
    my @b= Gapafim::sw($flank, lc($prefix));
    
    return $b[3];
}

sub calculate_yg_from_r1
{
    # Calculates YG from a window of $extra_flank from the beginning of R1

    my ($r1, $prefix, $prefix_length, $extra_flank)=@_;
    if($r1->[1] & 4)
    {
	#return YG=0 for unmapped reads
	return 0;
    }

    my $flank;

    if(!($r1->[1] & 16))
    {
	my $b = max(0, $r1->[3]);
	my $e = $r1->[3] + $extra_flank;
	
	$flank = rc($fai->fetch("$r1->[2]:$b-$e"));
    }
    else
    {
	my $b = max(0, $r1->[3] - $extra_flank);
	my $e = $r1->[3] + cigar2reflen($r1->[5]);
	
	$flank = $fai->fetch("$r1->[2]:$b-$e");
    }
    
    my @b= Gapafim::sw($flank, lc($prefix));
    
    return $b[3];
}


sub calculate_yg
{
    my ($r1, $r2, $prefix, $prefix_length, $r1_flank_length, $r2_flank_length)=@_;
    
    #If properly paired, use r2, else use r1

    if( (($r2->[1] & 4)==0) && ($r2->[1] & 2))
    {
	return calculate_yg_from_r2($r2, $prefix, $prefix_length, $r2_flank_length);
    }
    else
    {
	return calculate_yg_from_r1($r1, $prefix, $prefix_length, $r1_flank_length);
    }
}



sub read_name
{
    my $n=shift;
    $n=~/^([^\/]+)/;
    return $1;
}

my $prefix_length;
my $consensus;
my $r1_flank_length;
my $r2_flank_length;
my $soft_clip_length_threshold;
my $genome_fasta_file;


my $oo= GetOptions ("prefix_length=i" => \$prefix_length,
		    "consensus=s"=> \$consensus,
		    "r1_flank_length=i" => \$r1_flank_length,
		    "r2_flank_length=i" => \$r2_flank_length,
		    "soft_clip_length_threshold=i"=> \$soft_clip_length_threshold,
		    "genome_fasta_file=s"=> \$genome_fasta_file,
    );


$fai = Bio::DB::Sam::Fai->load($genome_fasta_file);


my @r1;
my @r2;

while(<>)
{
    if(/^\@/)
    {
	print;
	next;
    }

    my @a=split /\t/;
    chomp $a[-1];

    if(($a[1] & 256) || ($a[1] & 2048))
    {
	#skip secondary mappings
	print;
	next;
    }

    if($a[1] & 64)
    {
	@r1=@a;
    }
    elsif($a[1] & 128)
    {
	@r2=@a;
	
	die "read names don't match" if read_name($r2[0]) ne read_name($r1[0]);

	my $seq2=$r2[1] & 16 ? rc($r2[9]) : $r2[9];
	
	my $prefix=substr($seq2, 0, $prefix_length);
	
	{
	    push @r1, "Y2:Z:$seq2";
	}
	
	{
	    my @z= Gapafim::sw($consensus, $prefix);
	    #print join("\n",@z),"\n\n";
	    #print ">$prefix_length $prefix\n";
	    push @r1, "YA:i:$z[3]";
	    push @r2, "YA:i:$z[3]";
	}

	{
	    my $score=calculate_yg(\@r1, \@r2, $prefix, $prefix_length, $r1_flank_length, $r2_flank_length);
	    push @r1, "YG:i:$score";
	    push @r2, "YG:i:$score";
	    #print join("\n", $alu_aln, $aln, $prefix_aln, $score, $a[5], $a[1] & 16, "$a[2]:$a[3]-".($a[3]+100)),"\n";
	}

	{
	    my $yr=calculate_yr(\@r2, $soft_clip_length_threshold);
	    push @r1, "YR:i:$yr";
	    push @r2, "YR:i:$yr";
	}

	{
	    my $ys=calculate_ys(\@r2);
	    push @r1, "YS:i:$ys";
	    push @r2, "YS:i:$ys";
	}

	print join("\t", @r1),"\n";
	print join("\t", @r2),"\n";
    }
    
}
