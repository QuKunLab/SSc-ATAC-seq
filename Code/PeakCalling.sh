#!/bin/sh

#a sample of peakcalling for all normal samples
cd ./ATAC-pipe-master
python atac-pipe.py --PeakCalling --bed ./BedFiles -r hg19 -o ./outputFiles --group ./SourceData/NormalSamples_Groups.txt --p1 3 --q1 5 --pipeup 10 -u
