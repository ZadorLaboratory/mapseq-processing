#!/bin/bash -x
#
# retrieve all relevant inputs for completing project locally. 
#
DIRLIST="reads.out aggregated.out filtered.out readtable.out collapsed.out vbctable.out plots.out"
USERHOST="hover@zelk.cshl.edu"
PROJECTROOT="/home/hover/project/mapseq"
SEQTECH="novaseq"
#PROJECT="M295_207690"
#PROJECT_SHORT="M295"

if [ "$#" -ne 2 ]; then
	echo "usage: get_elzar_project.sh <PROJECT_ID> <SEQ_TECH> "
	echo "    SEQ_TECH = nextseq OR novaseq "
	exit 1
fi
	
PROJECT=$1
SEQTECH=$2
echo "project is $PROJECT seqtech is $SEQTECH"

mkdir $DIRLIST

for DIR in $DIRLIST ; do
 scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$DIR/stats* ./$DIR/
 scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$DIR/*.xlsx ./$DIR/
 scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$DIR/*.pdf ./$DIR/
 scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$DIR/*.compsize.tsv ./$DIR/
 scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/$DIR/brain*/*component_info.tsv ./$DIR/
done

scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/readtable.out/sampleinfo.tsv ./readtable.out/
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/vbctable.out/*.vbctable.* ./vbctable.out/
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/vbctable.out/target*.tsv ./vbctable.out/
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/vbctable.out/injection*.tsv ./vbctable.out/
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/vbctable.out/*.controls.tsv ./vbctable.out/
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/*.mapseq.conf ./
scp $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/*.sampleinfo.xlsx ./

scp -r $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/vbcfiltered.out ./
scp -r $USERHOST:$PROJECTROOT/$SEQTECH/$PROJECT/matrices.out ./

