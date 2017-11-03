#!/bin/bash

baseDir="$1"
codeDir="$baseDir/package"
outDir="$baseDir/output/fig2/"
tempDir="$baseDir/temp/fig2/"
pList="exdScr exdUbxIVa Max Hnf4a Ftz-f1"
fList="MGW.5mer ProT.5mer Roll.4mer HelT.4mer"
nRandom=1000

mkdir -p $outDir $tempDir

for p in $pList; do
	echo ">  Starting $p."
	modelFile="$baseDir/input/bindingModels/$p.mono.tsv"
	seqFile="$baseDir/input/highAffinitySequences/$p.txt"
	sampledSeqFile="$tempDir/$p.mono.sampled.fa"

	echo ">> Samples random sequences and bins by energy."
	$codeDir/sampleSequences.py $modelFile -n $nRandom >  $sampledSeqFile

	echo ">> Computes average shape profiles"

	for f in $fList; do

		#Creates input fhape files
		shapeFile="$tempDir/$f.csv"
		if [ "$(echo "$f" | cut -d '.' -f2)" = "5mer" ]; then
			cp $baseDir/input/shapeTables/$f.csv $shapeFile
		else
			fIn="$(echo $f | sed 's/4mer/5mer/g')"
			cat $baseDir/input/shapeTables/$fIn.csv | awk -F ',' '{s[substr($1,1,4)]+=$2;s[substr($1,2,4)]+=$4}  END {for(i in s)print i","s[i]/8}'| sort > $shapeFile
		fi

		#Output shape profile files
		profileFile="$outDir/$p.$f.highAffinity.lst"; >$profileFile
		sampledProfileFile="$outDir/$p.$f.sampled.lst"; >$sampledProfileFile

		#Number of "N" added to the sampled files
		if [ $(echo "$f" | cut -d '.' -f2) = "4mer" ]; then
			nPad=1
		else
			nPad=2
		fi

		#Mean shape of binnned sampled sequences
		seq 10 | while read iE; do 
			profile="$(cat $sampledSeqFile | grep -A 1 "iE=$iE," | grep -v "^--$" | grep -v ">" | awk -v nPad=$nPad '{s=$0;for(i=0;i<nPad;i++)s="N"s"N";print s}' | $codeDir/computeMeanProfile.py --scoreN - $shapeFile | tr '\n' ','| sed 's/,$//g')"
			echo "iE_$iE,$profile" >> $sampledProfileFile
		done

		#mean shape of high-affninity probes
		$codeDir/computeMeanProfile.py --scoreN $seqFile $shapeFile | tr '\n' ','| sed 's/,$//g' > $profileFile
	done
done
