################################# parameter settings below #################################
######### required below #########
inputFolder=input
outputFolder=output

abdFile=abd.tsv
covFile=cov.tsv
sampleIDCol=sample_id

covariateInterest=sex
covariateAdjust="age,PC1,PC2,batch"

######### optional below #########
outPrefix=PALM_results






###################### change to absolute paths in docker (do not modify) ########################
workDir="$(pwd)"
inputFolder="${workDir}/${inputFolder}"
outputFolder="${workDir}/${outputFolder}"
abdFile="${inputFolder}/${abdFile}"
covFile="${inputFolder}/${covFile}"
outPrefix="${outputFolder}/${outPrefix}"

################################# workflow below (do not modify) #################################
mkdir -p "${outputFolder}"

echo "[INFO] Running PALM on example data ..."
run_PALM.R \
  --abdFile "${abdFile}" \
  --covFile "${covFile}" \
  --sampleIDCol "${sampleIDCol}" \
  --covariateInterest "${covariateInterest}" \
  --covariateAdjust "${covariateAdjust}" \
  --outPrefix "${outPrefix}"

echo "[INFO] Done. Results in ${outPrefix}.rds and ${outPrefix}.tsv"
