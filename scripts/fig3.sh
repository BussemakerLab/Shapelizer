#!/bin/bash

baseDir="$1"
codeDir="$baseDir/package"
outDir="$baseDir/output/fig3/"
tempDir="$baseDir/temp/fig3/"
tableDir="$baseDir/input/shapeTables"
fList="MGW.5mer ProT.5mer Roll.4mer HelT.4mer"

mkdir -p $outDir $tempDir

for f in $fList; do
	echo ">  Starting $f"
	echo ">> $f, Mono, cross-validated R2."
	r2Mono="$r2Mono,$($codeDir/kMerLinearRegression.py $tableDir/$f.csv --mono --crossR2)"
	echo ">> $f, Di,   cross-validated R2."
	r2Di="$r2Di,$(    $codeDir/kMerLinearRegression.py $tableDir/$f.csv --di   --crossR2)"
	echo ">> $f, All,  cross-validated R2."
	r2All="$r2All,$(  $codeDir/kMerLinearRegression.py $tableDir/$f.csv --all  --crossR2)"

	echo ">> $f, Computes regression coefficients."
	$codeDir/kMerLinearRegression.py $tableDir/$f.csv --mono --betas > $outDir/$f.mono.betas.tsv
	$codeDir/kMerLinearRegression.py $tableDir/$f.csv --di   --betas > $outDir/$f.di.betas.tsv
	$codeDir/kMerLinearRegression.py $tableDir/$f.csv --all  --betas > $outDir/$f.all.betas.tsv

	echo ">> $f, Computes predicted values."
	$codeDir/kMerLinearRegression.py $tableDir/$f.csv --mono --yHat  > $outDir/$f.mono.yHat.csv
	$codeDir/kMerLinearRegression.py $tableDir/$f.csv --di   --yHat  > $outDir/$f.di.yHat.csv
	$codeDir/kMerLinearRegression.py $tableDir/$f.csv --all  --yHat  > $outDir/$f.all.yHat.csv

	
	echo ">> $f, Creating permuted shape table."
	$codeDir/randomKmerTable.py      $tableDir/$f.csv --permute -s 0 > $tempDir/$f.permuted.csv

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

echo "Cross-validated R2:"
cat $outDir/r2.csv
echo "Cross-validated R2, permuted table"
cat $outDir/r2.permuted.csv

