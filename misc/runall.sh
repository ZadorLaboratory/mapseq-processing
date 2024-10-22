#!/bin/bash -x

SCRIPTDIR=~/git/mapseq-processing/scripts
EXP='M205.htna24'
EXPTAG='M205'
OUTROOT=~/project/mapseq/$EXP
SIF=$OUTROOT/${EXPTAG}_sampleinfo.xlsx
#DEBUG='$DEBUG'
DEBUG='-v'

echo "SIF=$SIF"

echo time $SCRIPTDIR/process_fastq_pairs.py $DEBUG -o $OUTROOT/fastq.out/$EXP.reads.tsv  $OUTROOT/fastq/$EXPTAG*
time $SCRIPTDIR/process_fastq_pairs.py $DEBUG -o $OUTROOT/fastq.out/$EXP.reads.tsv  $OUTROOT/fastq/$EXPTAG*


echo time $SCRIPTDIR/align_collapse.py $DEBUG -o $OUTROOT/collapse.out/$EXP.collapse.tsv  $OUTROOT/fastq.out/$EXP.reads.tsv
time $SCRIPTDIR/align_collapse.py $DEBUG -o $OUTROOT/collapse.out/$EXP.collapse.tsv  $OUTROOT/fastq.out/$EXP.reads.tsv

echo time $SCRIPTDIR/make_readtable.py $DEBUG -s $SIF  -o $OUTROOT/readtable.out/$EXP.readtable.tsv  $OUTROOT/collapse.out/$EXP.collapse.tsv
time $SCRIPTDIR/make_readtable.py $DEBUG -s $SIF -o $OUTROOT/readtable.out/$EXP.readtable.tsv  $OUTROOT/collapse.out/$EXP.collapse.tsv


echo time $SCRIPTDIR/make_vbctable.py $DEBUG  -o $OUTROOT/vbctable.out/$EXP.vbctable.tsv  $OUTROOT/readtable.out/$EXP.readtable.tsv
time $SCRIPTDIR/make_vbctable.py $DEBUG  -o $OUTROOT/vbctable.out/$EXP.vbctable.tsv  $OUTROOT/readtable.out/$EXP.readtable.tsv

echo time $SCRIPTDIR/make_matrices.py $DEBUG -O $OUTROOT/matrices.out/  $OUTROOT/vbctable.out/$EXP.vbctable.tsv
time $SCRIPTDIR/make_matrices.py $DEBUG -O $OUTROOT/matrices.out/  $OUTROOT/vbctable.out/$EXP.vbctable.tsv
