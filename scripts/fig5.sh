#!/bin/bash

baseDir="$1"
codeDir="$baseDir/package"
outDir="$baseDir/output/fig5"
tempDir="$baseDir/temp/fig5"
modelDir="$baseDir/input/bindingModels"
pList="exdScr exdUbxIVa Max CEBPb"
fList="MGW.5mer ProT.5mer Roll.4mer HelT.4mer"

nSample=10

mkdir -p $outDir $tempDir


echo ">  Performs shape projection for all TF-shape combinations:"
for f in $fList; do

	#Creates input fhape files
	if [ "$(echo "$f" | cut -d '.' -f2)" = "5mer" ]; then
		cp $baseDir/input/shapeTables/$f.csv $tempDir/
	else
		fIn="$(echo $f | sed 's/4mer/5mer/g')"
		cat $baseDir/input/shapeTables/$fIn.csv | awk -F ',' '{s[substr($1,1,4)]+=$2;s[substr($1,2,4)]+=$4}  END {for(i in s)print i","s[i]/8}'| sort > $tempDir/$f.csv
	fi
	shapeFile="$tempDir/$f.csv"

	echo ">> Creates mono+di shape model file for $f" 
	$codeDir/kMerLinearRegression.py $shapeFile --di   --betas > $tempDir/$f.di.betas.tsv

	for p in $pList; do
		echo ">> Running $f projection for $p"
		$codeDir/shapeProjection.py $modelDir/$p.monoDi.tsv $tempDir/$f.di.betas.tsv > $outDir/$p.$f.projectionBetas.tsv
	done
done

echo ">  Performs shape projection with L1 or L2 penalty terms and with different values of the lambda scale paramters (Supplemental figure)"
f="MGW.5mer"
lambdaList="0,0.125,0.25,0.5,1,2,4,8"
for p in exdScr exdUbxIVa; do 
	for LMono in L1Mono L2Mono; do
		for LShape in L1Shape L2Shape; do
			echo ">> Running $f projection for $p with -$LMono -$LShape"
			$codeDir/shapeProjection.py $modelDir/$p.monoDi.tsv $tempDir/$f.di.betas.tsv -$LMono $lambdaList -$LShape $lambdaList --header > $outDir/$p.$f.$LMono.$LShape.projectionBetas.tsv
		done
	done
done



echo "> Runs shape projection for random tables."
echo "! This scripts samples nSample=$nSample random tables. Figure 5 is based on 5,000 random tables."
seq $nSample | while read s; do 
	echo ">> Starting randomization with seed=${s}."
	for f in $fList; do
		randomShapeDir="$tempDir/randomTables/$f"
		mkdir -p $randomShapeDir

		echo ">>>Generating random table: $f, $s."
		$codeDir/randomKmerTable.py $shapeFile -s $s --symmetric | tee $randomShapeDir/$f.s$s.csv | $codeDir/kMerLinearRegression.py - --di --betas | grep -v "^intercept" > $randomShapeDir/$f.s$s.di.betas.tsv 

		echo ">>>Running shape projection: $f, $s."
		for p in $pList; do
			randomProjectionDir="$outDir/randomTableProjection/$p.$f/"
			mkdir -p $randomProjectionDir
			$codeDir/shapeProjection.py $modelDir/$p.monoDi.tsv $randomShapeDir/$f.s$s.di.betas.tsv > $randomProjectionDir/$p.$f.s$s.projectionBetas.tsv
		done
	done
done

