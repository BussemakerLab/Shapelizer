#!/bin/bash

baseDir="$1"
codeDir="$baseDir/package"
outDir="$baseDir/output/fig3/"
tempDir="$baseDir/temp/fig3/"
fList="MGW.5mer ProT.5mer Roll.4mer HelT.4mer"

mkdir -p $outDir $tempDir

for f in $fList; do

	#Creates input shape files
	if [ "$(echo "$f" | cut -d '.' -f2)" = "5mer" ]; then
		cp $baseDir/input/shapeTables/$f.csv $tempDir/
	else
		fIn="$(echo $f | sed 's/4mer/5mer/g')"
		cat $baseDir/input/shapeTables/$fIn.csv | awk -F ',' '{s[substr($1,1,4)]+=$2;s[substr($1,2,4)]+=$4}  END {for(i in s)print i","s[i]/8}'| sort > $tempDir/$f.csv
	fi
	shapeFile="$tempDir/$f.csv"

	echo ">  Starting $f"
	echo ">> $f, Mono,   cross-validated R2."
	r2Mono="$r2Mono,$(      $codeDir/kMerLinearRegression.py $shapeFile --mono         --crossR2)"
	echo ">> $f, Di,     cross-validated R2."
	r2Di="$r2Di,$(          $codeDir/kMerLinearRegression.py $shapeFile --di           --crossR2)"
	echo ">> $f, All,    cross-validated R2."
	r2All="$r2All,$(        $codeDir/kMerLinearRegression.py $shapeFile --all          --crossR2)"
	echo ">> $f, Central,cross-validated R2."
	r2Central="$r2Central,$($codeDir/kMerLinearRegression.py $shapeFile --centralMono  --crossR2)"

	echo $r2Central
	echo ">> $f, Computes regression coefficients."
	$codeDir/kMerLinearRegression.py $shapeFile --mono --betas > $outDir/$f.mono.betas.tsv
	$codeDir/kMerLinearRegression.py $shapeFile --di   --betas > $outDir/$f.di.betas.tsv
	$codeDir/kMerLinearRegression.py $shapeFile --all  --betas > $outDir/$f.all.betas.tsv

	echo ">> $f, Computes predicted values."
	$codeDir/kMerLinearRegression.py $shapeFile --mono --yHat  > $outDir/$f.mono.yHat.csv
	$codeDir/kMerLinearRegression.py $shapeFile --di   --yHat  > $outDir/$f.di.yHat.csv
	$codeDir/kMerLinearRegression.py $shapeFile --all  --yHat  > $outDir/$f.all.yHat.csv

	
	echo ">> $f, Creating permuted shape table."
	$codeDir/randomKmerTable.py      $shapeFile --permute -s 0 > $tempDir/$f.permuted.csv

	echo ">> $f, Mono, cross-validated R2, permuted."
	r2MonoPerm="$r2MonoPerm,$($codeDir/kMerLinearRegression.py $tempDir/$f.permuted.csv --mono --crossR2)"
	echo ">> $f, Di,   cross-validated R2, permuted."
	r2DiPerm="$r2DiPerm,$(    $codeDir/kMerLinearRegression.py $tempDir/$f.permuted.csv --di   --crossR2)"
	echo ">> $f, All,  cross-validated R2, permuted."
	r2AllPerm="$r2AllPerm,$(  $codeDir/kMerLinearRegression.py $tempDir/$f.permuted.csv --all  --crossR2)"
done

echo " $fList" | tr ' ' ','     >  $outDir/r2.permuted.csv
echo "mono$r2MonoPerm"          >> $outDir/r2.permuted.csv
echo "di$r2DiPerm"              >> $outDir/r2.permuted.csv
echo "all$r2AllPerm"            >> $outDir/r2.permuted.csv

echo " $fList" | tr ' ' ',' >  $outDir/r2.csv
echo "mono$r2Mono"          >> $outDir/r2.csv
echo "di$r2Di"              >> $outDir/r2.csv
echo "all$r2All"            >> $outDir/r2.csv
echo "central$r2Central"    >> $outDir/r2.csv

echo "Cross-validated R2:"
cat $outDir/r2.csv
echo "Cross-validated R2, permuted table"
cat $outDir/r2.permuted.csv

