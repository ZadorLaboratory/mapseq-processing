#!/bin/bash -x
#
# retrieve all relevant inputs for completing project locally. 
#
DIRLIST="reads.out aggregated.out filtered.out readtable.out collapsed.out vbctable.out plots.out vbcfiltered.out matrices.out"
USERHOST="hover@bamdev2.cshl.edu"
PROJECTROOT="/grid/zador/home/hover/project/mapseq"
SEQTECH="novaseq"
#PROJECT="M295_207690"
#PROJECT_SHORT="M295"

if [ "$#" -ne 3 ]; then
	echo "usage: get_elzar_project.sh <PROJECT_ID> <RUN_SUBDIR> <SEQ_TECH> "
	echo "    SEQ_TECH = nextseq OR novaseq "
	exit 1
fi
	
PROJECT=$1
RUNDIR=$2
SEQTECH=$3
echo "project is $PROJECT seqtech is $SEQTECH"

mkdir $DIRLIST

for DIR in $DIRLIST ; do
 scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/$DIR/stats* ./$DIR/
 scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/$DIR/*.xlsx ./$DIR/
 scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/$DIR/*.pdf ./$DIR/
 scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/$DIR/*.compsize.tsv ./$DIR/
 scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/$DIR/brain*/*component_info.tsv ./$DIR/
done

scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/readtable.out/sampleinfo.tsv ./readtable.out/
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/collapsed.out/*.collapsed.tsv ./collapsed.out/
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/collapsed.out/*.collapsed.parquet ./collapsed.out/
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/vbctable.out/*.vbctable.* ./vbctable.out/
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/vbctable.out/target*.tsv ./vbctable.out/
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/vbctable.out/injection*.tsv ./vbctable.out/
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/vbctable.out/*.controls.tsv ./vbctable.out/
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/*.mapseq.conf ./
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/*.sampleinfo*.xlsx ./

scp -r $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/vbcfiltered.out ./
scp -r $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$RUNDIR/matrices.out ./

