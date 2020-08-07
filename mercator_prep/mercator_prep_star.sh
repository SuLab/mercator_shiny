#!/bin/bash

chrom_file="$(dirname $0)/chrom_sizes_GRCh38.p7.txt"
bed="$(dirname $0)/Gencode-v25.mod.bed"
out_flag=false
bw_flag=true

# if [[ "$#" -eq 0 ]]
# then
#     printf >&2 "No files specified, run with -h for help
# Usage: $(basename $0) [-h] [-b bedfile] [-c chrom_size_file] -o out_file in_files\n"
#     exit 1
# fi

function script_help { printf >&1 \
"Script for gene counting from wiggle files, in a manner similar to Recount2, for mercatorproject.com
Requires: UCSC tools wigToBigWig and bigWigAverageOverBed
Usage: $(basename $0) [-hw] [-b bedfile] [-c chrom_size_file] -o out_file in_files
       -h: print this help message then exit
       -w: keep bigwig files, default is to delete them in order to reduce footprint
       -b: location of modified bedfile; default is /path/of/script/Gencode-v25.mod.bed
       -c: location of chromosome sizes; default is /path/of/script/chrom_sizes_GRCh38.p7.txt
       -o: desired output file
       in_files: input wiggle files generated using STAR or Rails-RNA
Output: .tsv file with gene counts for samples (wiggle files) in each column\n"
		     }

while getopts "c:b:o:hw" opt
do
    case $opt in
	b) bed=$OPTARG;;
	c) chrom_file=$OPTARG;;
	h) script_help; exit;;
	w) bw_flag=false;;
	o) out_file=$OPTARG; out_flag=true;;
	[?])	printf >&2 "Usage: $(basename $0) [-hw] [-b bedfile] [-c chrom_size_file] -o out_file in_files\n"
		exit 1;; 
    esac
done	    

if ! [[ -f $bed ]]
then
    printf >&2 "Must specify valid bedfile, specified does not exist
Usage: $(basename $0) [-hw] [-b bedfile] [-c chrom_size_file] -o out_file in_files\n"
    exit 1
fi

if ! [[ -f $chrom_file ]]
then
    printf >&2 "Must specify valid chromosome size file, specified does not exist
Usage: $(basename $0) [-hw] [-b bedfile] [-c chrom_size_file] -o out_file in_files\n"
    exit 1
fi

if ! $out_flag || ! [[ -d $(dirname $out_file) ]] 
then
    printf >&2 "Output file was not specified or the directory specified with -o does not exist
Usage: $(basename $0) [-hw] [-b bedfile] [-c chrom_size_file] -o out_file in_files\n"
    exit 1
fi

shift $((OPTIND-1))

if [[ $# -eq 0 ]]
then
    printf >&2 "No input files specified, run with -h for help
Usage: $(basename $0) [-hw] [-b bedfile] [-c chrom_size_file] -o out_file in_files\n"
    exit 1
fi

wigToBigWig $1 $chrom_file ${1%wig}bw 
bigWigAverageOverBed ${1%wig}bw $bed stdout | \
    sed 's/xXx[[:digit:]]\+//' | \
    awk '{a[$1] += $4} END {for (i in a) {printf "%-15s\t%s\n",i,a[i];}}' | \
    sed "1 s@^@$1\n@" > $out_file

if $bw_flag
then
    rm ${1%wig}bw
fi


shift 1    

for file in $@
do

    wigToBigWig $file $chrom_file ${file%wig}bw
    bigWigAverageOverBed ${file%wig}bw $bed stdout | \
	sed 's/xXx[[:digit:]]\+//' | \
	awk '{a[$1] += $4} END {for (i in a) {printf "%s\n",a[i];}}' | \
	sed "1 s@^@$file\n@" |  \
	paste $out_file -  > ${out_file}.tmp

    mv ${out_file}.tmp $out_file

    if $bw_flag
    then
	rm ${file%wig}bw
    fi

    
done
