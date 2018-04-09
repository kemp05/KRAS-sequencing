#!/usr/bin/perl -w

# Processes NUCLEOTIDE_COUNT_FOR_HASH.pl output
# Calculates SNV allele frequencies and calls mutations


# name file for results
print "\nEnter filename for results\n";
$results = <STDIN>;
chomp($results);
open(OUTPUT, ">$results");

# hash table population (using output from NUCLEOTIDE_COUNT_FOR_HASH.pl)
%TABLE=();
open (INPUT3, "amp1_SLX_14588_BH964_read_counts_hash.txt") or die "\nAmp file could not be found\n\n";
while (<INPUT3>) {
	chomp ($_);
	@AMP = split('\t',$_);
	$index = 1;
	while ($index <= scalar@AMP){
		$key = "$AMP[0]"."_$index";
		$TABLE{$key} = $AMP[$index];
		$index = $index +1;
	}
}

# make array of sample barcode IDs - allowing for differing leading 0s
@SAMPLES = ();
$sample = 1;

while ($sample <= 9){
	$sample = "00"."$sample";
	push(@SAMPLES,$sample);
	$sample = $sample +1;
}

while ($sample <= 99){
	$sample = "0"."$sample";
	push(@SAMPLES,$sample);
	$sample = $sample +1;
}

while ($sample <= 384){
	push(@SAMPLES,$sample);
	$sample = $sample +1;
}


# calculate mean and st.dev for every nucleotide-position/variant-nucleotide 
@BASES = ('A','C','G','T');
$nt_position = 28;
while ($nt_position <= 139){

	foreach $base (@BASES){

# collect relevant nucleotide-position/variant-nucleotide values.  Not enough reads = -1
		@NTs = ();
		$values = 0;
		foreach $val (@SAMPLES){
			$key = "FLD0$val"."_$base"."_$nt_position";
			$sum_key = "FLD0$val"."_Sum"."_$nt_position";
			if ($TABLE{$sum_key} >= 1000){
				$percent_nt = ($TABLE{$key} / $TABLE{$sum_key}) *100;
			}else{
				$percent_nt = -1;
			}
			push (@NTs,$percent_nt);
		}

# calculate statistics

		$array_sum = 0;
		$array_mean = 0;
		$array_stdev = 0;
		$array_variance = 0;
		$element_variance = 0;
		$array_sum_variance = 0;
		
		foreach $val (@NTs){
			if ($val >= 0){
				$array_sum = $array_sum + $val;
				$values = $values + 1;				
			}
		}

		$array_mean = $array_sum / $values;

		foreach $val (@NTs){
			if ($val >= 0){	
				$element_variance = ($val - $array_mean) * ($val - $array_mean);
				$array_sum_variance = $array_sum_variance + $element_variance;
			}
		}

		$array_variance = $array_sum_variance / ($values -1);
		$array_stdev = sqrt($array_variance);
		$cut_off = $array_mean + ($array_stdev * 3.290527);
	
# mutation calling: filter all individual sample values with cut-offs.
		$index =0;
		while ($index <= scalar@NTs){
			if (($NTs[$index] >= 0.1) && (($NTs[$index] >= $cut_off) || ($NTs[$index] >= $array_mean *4))){
				print OUTPUT "FLD0$SAMPLES[$index]\t$nt_position\t$base\t$NTs[$index]\t$array_mean\t$array_stdev\t$values\n";				
			}
			$index = $index +1;
		}
	}
	$nt_position = $nt_position +1;
}

exit;

