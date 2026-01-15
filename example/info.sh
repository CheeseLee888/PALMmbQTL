################################# parameter settings below #################################
inputFolder=input
outputFolder=output
genoFile=./${inputFolder}/geno
outputFile=./${outputFolder}/snp_info.txt

################################# workflow below (do not modify) #################################
mkdir -p "${outputFolder}"
echo "Generating SNP information from genotype data..."
pixi run --manifest-path=../pixi.toml Rscript ../extdata/step0_snpInfo.R \
    --genoFile=${genoFile} \
    --outputFile=${outputFile}