#!/bin/bash -x
#
# retrieve all relevant inputs for completing project locally. 
#
DIRLIST="reads.out aggregated.out filtered.out readtable.out collapsed.out vbctable.out"
USERHOST="hover@bamdev2.cshl.edu"
PROJECTROOT="/grid/zador/home/hover/project/mapseq"
PROJECT="M295_207690"
PROJECT_SHORT="M295"

mkdir $DIRLIST

for DIR in $DIRLIST ; do
 scp $USERHOST:$PROJECTROOT/$PROJECT/$DIR/stats* ./$DIR/
 scp $USERHOST:$PROJECTROOT/$PROJECT/$DIR/*.xlsx ./$DIR/
 scp $USERHOST:$PROJECTROOT/$PROJECT/$DIR/*.pdf ./$DIR/
done

scp $USERHOST:$PROJECTROOT/$PROJECT/readtable.out/sampleinfo.tsv ./readtable.out/
scp $USERHOST:$PROJECTROOT/$PROJECT/vbctable.out/*.vbctable.* ./vbctable.out/
scp $USERHOST:$PROJECTROOT/$PROJECT/$PROJECT_SHORT.mapseq.conf ./
