#!/bin/bash
# 
# Download Novaseq FASTQ data from Johns Hopkins
# Create user-only-readable password file in home directory, named accordingly
#
# Usage:  ./aspera-download.sh <dataset>
#
# Dataset names are typically something like 'hzhan_123456' per Hopkins email.
#

USAGE="usage: aspera-download.sh <dataset>"

ASP_URL="https://jhg-aspera2.jhgenomics.jhu.edu"
ASP_USER="hzhan@cshl.edu"
ASP_OUTDIR="/grid/mbseq/data_norepl/mapseq/raw_data/novaseq/"
ASP_PWFILE=~/hzhanpw.txt
ASP_PWS=`cat $ASP_PWFILE`
#echo "pw is $ASP_PWS"
NARGS=$#
#echo "arg number is $NARGS"

if [ "$NARGS" -ne 1 ] ; then 
	echo $USAGE
	exit 1
fi
DATASET=$1

echo ascli shares files download --url=$ASP_URL --username=$ASP_USER --password=$ASP_PWS /cidr-zhan/$DATASET --ts=@json:'{"target_rate_kbps":300000,"resume_policy":"sparse_csum"}' --to-folder=$ASP_OUTDIR
ascli shares files download --url=$ASP_URL --username=$ASP_USER --password=$ASP_PWS /cidr-zhan/$DATASET --ts=@json:'{"target_rate_kbps":300000,"resume_policy":"sparse_csum"}' --to-folder=$ASP_OUTDIR