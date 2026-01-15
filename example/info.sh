################################# parameter settings below #################################
inputFolder=input
outputFolder=output
genoFile=./${inputFolder}/geno
abdFile=./${inputFolder}/abd.txt
outputSnpFile=./${outputFolder}/info_snp.txt
outputFeatureFile=./${outputFolder}/info_feature.txt

################################# workflow below (do not modify) #################################
mkdir -p "${outputFolder}"
echo "Generating SNP information from genotype data..."
pixi run --manifest-path=../pixi.toml Rscript ../extdata/step0_snpInfo.R \
    --genoFile=${genoFile} \
    --outputFile=${outputSnpFile}

echo "Generating feature information from phenotype data..."
pixi run --manifest-path=../pixi.toml Rscript ../extdata/step0_featureInfo.R \
    --abdFile=${abdFile} \
    --outputFile=${outputFeatureFile}