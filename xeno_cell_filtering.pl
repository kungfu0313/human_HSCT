#!usr/bin/perl -w
use List::Util qw/sum max min/;
use Statistics::Descriptive;
use strict;

##--------- Step1: process sequencing data by mapping to both human and pig reference genome first 

##--------- Step2: Prepare input files by R based on cellranger output 
#Read in pig count
#data <- Read10X(data.dir = "filtered_feature_bc_matrix")
#scRNA = CreateSeuratObject(counts = data, min.cells = 0,  min.features = 200, project = "cellranger")
#count <- GetAssayData(scRNA,slot="data",assay="RNA")
#data <- t(count)
#write.csv(data, "./data_pig.csv", row.names = T,quote=F)
#Read in human count
#scRNA <- readRDS("./Thymus_human_20230301_cellcycled.rds")
#count <- GetAssayData(scRNA,slot="data",assay="RNA")
#data <- t(count)
#write.csv(data, "./data_human.csv", row.names = T,quote=F)


cell_HP_counting(); ##--------- Step3: Counting number of expressed human/pig genes in each cell
merge_HP_2023();    ##--------- Step4: put counts of expressed human genes and pig genes together, one cell per row
collect_pigmaybe(); ##--------- Step5: potential xeno(pig) cells


sub cell_HP_counting
{
		#consider expressed MT genes specifically    
    #open(IN, "/PATH/TO/FILES/human_MTglist.txt") or die ("Cannot open the file $_!");
    open(IN, "/PATH/TO/FILES/pig_MTglist.txt") or die ("Cannot open the file $_!");

    my %mt;
    while(<IN>)
    {
			chomp;
			%mt = (%mt, $_, 1);
    }
    close IN;
    #open(IN, "/PATH/TO/FILES/data_human.csv") or die ("Cannot open the file $_!");
    #open(OUT, ">/PATH/TO/FILES/counting_human.txt") or die ("Cannot open the file $_!");
    open(IN, "/PATH/TO/FILES/data_pig.csv") or die ("Cannot open the file $_!");
    open(OUT, ">/PATH/TO/FILES/counting_pig.txt") or die ("Cannot open the file $_!");

    my %mt_pos;
    my %count_all;
    my %count_mt;
    my %count_sum;
    my %mt_sum;
    my $i = 0;
    while(<IN>)
    {
	chomp;
	my @info = split(/\,/, $_);
	if($i eq 0) 
	{
	    for(my $j = 0; $j <$#info+1; $j ++)
	    {
		if(exists($mt{$info[$j]})) {%mt_pos = (%mt_pos, $j, $info[$j]);}
	    }
	    foreach my $key(keys %mt_pos)
	    {
		print OUT "$mt_pos{$key},";
	    }
	    print OUT "\ncell_ID\tcount_all\tcount_mt\tsum_all\tsum_mt\n";
	
	}
	else
	{
	    %count_all = (%count_all, $info[0], 0);
	    %count_mt = (%count_mt, $info[0], 0);
	    %count_sum = (%count_sum, $info[0], 0);
	    %mt_sum = (%mt_sum, $info[0], 0);

	    for(my $j = 1; $j <$#info+1; $j ++)
	    {
		if($info[$j] ne 0)
		{
		   $count_all{$info[0]} ++;
		   $count_sum{$info[0]}+= $info[$j];
		   if(exists($mt_pos{$j}))
		   {
		       $count_mt{$info[0]} ++;
		       $mt_sum{$info[0]}+= $info[$j];
		   }
		}
	    }
	    print OUT "$info[0]\t$count_all{$info[0]}\t$count_mt{$info[0]}\t$count_sum{$info[0]}\t$mt_sum{$info[0]}\n";
	}
	$i ++;
    }
    close IN;
    close OUT;
}

sub merge_HP_2023
{
	open(IN, "/PATH/TO/FILES/counting_pig.txt") or die ("Cannot open the file $_!");
	my %pcell;
	my $i = 0;
	while(<IN>)
	{
			chomp;
			if($i > 0)
			{
					my @info = split(/\t/, $_); 
					%pcell = (%pcell, $info[0], $_);
			}
			$i ++;
	}
	close IN;
	open(IN, "/PATH/TO/FILES/counting_human.txt") or die ("Cannot open the file $_!");
	open(OUT, ">/PATH/TO/FILES/counting_HP.txt") or die ("Cannot open the file $_!");
	$i = 0;
	while(<IN>)
	{
			chomp;
			if($i > 0)
			{
					my @info = split(/\t/, $_); 
					if(exists($pcell{$info[0]}))
					{
						print OUT "$_\t$pcell{$info[0]}\n";	
					}
					else
					{
							print OUT "$_\t$info[0]\t0\t0\t0\t0\n";
					}
			}
			$i ++;
	}
	close IN;
	close IN;
	close OUT;
}

sub collect_pigmaybe
{
	open(IN, "/PATH/TO/FILES/counting_HP.txt") or die ("Cannot open the file $_!");
	open(OUT, ">/PATH/TO/FILES/pigmaybe.txt") or die ("Cannot open the file $_!");
	my $i = 0;
	while(<IN>)
	{
		chomp;
		if($i > 0)
		{
			my @info = split(/\t/, $_);
			my $count1 = $info[6] - $info[1];
			my $count2 = $info[7] - $info[2];
			if($count1 >0 && $count2 >= 4) { print OUT "$_\tpig\n";}
			else {print OUT "$_\thuman\n";}
		}
		else
		{
			print OUT "$_\tpigmaybe\n";	
		}
		$i ++;
	}	
	close IN;
	close OUT;
}


##--------- Step6: update metadata and filter out potential xeno(pig) clusters(cells) by R 