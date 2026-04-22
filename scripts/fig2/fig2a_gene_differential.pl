use strict;use warnings;
use Cwd;

use lib "./";
use ERR_CAT;

&ERR_CAT;

my $Ribo_SAM_List=shift;

open(F,$Ribo_SAM_List);
my $SAM_dir=<F>;chomp $SAM_dir;
my $Control=<F>;chomp $Control;
my $Treatment=<F>;chomp $Treatment;
close(F);

my @Controls=split(/,/,$Control);
my @Treatments=split(/,/,$Treatment);

shift @Controls;shift @Treatments;

my %data;my %rpkm_table;
foreach my $file(@Controls){
	&READ_FILE("$SAM_dir\\RPKM\\read_count_$file",\%data,"control",$file);
	&READ_FILE("$SAM_dir\\RPKM\\rpkm_$file",\%rpkm_table,"control",$file);
}


foreach my $file(@Treatments){
	&READ_FILE("$SAM_dir\\RPKM\\read_count_$file",\%data,"treatment",$file);
	&READ_FILE("$SAM_dir\\RPKM\\rpkm_$file",\%rpkm_table,"treatment",$file);
}


my $nControl=scalar(@Controls);my $nTreatment=scalar(@Treatments);

if($nControl<2 or $nTreatment<2){
	my $output_file="$SAM_dir\\RPKM\\rpkm_table.csv";
	&OUTPUT_DATA(\%rpkm_table,$output_file);
}
else{
	&OUTPUT_DATA(\%data,"read_count");
	system("Rscript gene_differential.R $nControl $nTreatment");
}


sub READ_FILE{
	my $file=shift;
	my $ref_hash=shift;
	
	my $flag=shift;
	my $sample=shift;
	
	open(F,$file);
	$sample=~s/\.bam//;
	while(my $line=<F>){
		my($name,$value)=$line=~/(\S+),(\S+)/;
		$$ref_hash{$name}{$flag}{$sample}=$value;
	}
	close(F);
}
	
sub OUTPUT_DATA{
	my $ref_hash=shift;
	my $output_file=shift;
	
	my $header='gene_name';my $header_line=0;
	
	open(F,">$output_file");
	foreach my $name(keys %{$ref_hash}){
		my $output=$name;
		foreach my $control_file(@Controls){
			$control_file=~s/\.bam//;
			$header.=",$control_file";
		
			unless(exists $$ref_hash{$name}{'control'}{$control_file}){
				$output.=",0";
			}
			else{
				$output.=",$$ref_hash{$name}{'control'}{$control_file}";
			}
		}
		
		foreach my $treatment_file(@Treatments){
			$treatment_file=~s/\.bam//;
			$header.=",$treatment_file";
			unless(exists $$ref_hash{$name}{'treatment'}{$treatment_file}){
				$output.=",0";
			}
			else{
				$output.=",$$ref_hash{$name}{'treatment'}{$treatment_file}";
			}
		}
	
		if($header_line==0){
			print F $header,"\n";
			$header_line=1;
		}
		print F $output,"\n";
	}
	close(F);
}
