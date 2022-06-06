#!/bin/sh

# These paths should be relative to the root folder of this repository.
dataDirPath="data"
outputFileName="results"
sopangExe="sopang"
bndmExe="../eds_search-master/eds_search"

minChromosomeId=25
maxChromosomeId=28

cd ..

if [ ! -f ${sopangExe} ];
then
    echo "Executable does not exist: ${sopangExe}"
    exit 1
fi

# Running ED-text matching with sources.
#for i in $(seq ${minChromosomeId} ${maxChromosomeId});
#do
#    for sourceCount in 32 64 128 256 512
#    do
        #./${sopangExe} --ceds ${dataDirPath}/chr${i}.eds ${dataDirPath}/patterns8.txt -d -o ./results/${outputFileName}_8.txt --in-sources-file ${dataDirPath}/chr${i}_${sourceCount}.edss > ./results/${outputFileName}_${i}_${sourceCount}.txt
        #./${sopangExe} ${dataDirPath}/chr${i}.eds ${dataDirPath}/patterns8.txt -d -o ./results/${outputFileName}_8sop.txt --in-sources-file ${dataDirPath}/chr${i}_${sourceCount}.edss --full-sources-output > ./results/${outputFileName}_${i}_${sourceCount}_sopang.txt
        #./${sopangExe} --ceds ${dataDirPath}/chr${i}.eds ${dataDirPath}/patterns8.txt -d -o ${outputFileName}_8_sources.txt --in-sources-file ${dataDirPath}/chr${i}_${sourceCount}.edss --full-sources-output

        #./${sopangExe} --ceds ${dataDirPath}/chr${i}.eds ${dataDirPath}/patterns16.txt -d -o ${outputFileName}_16_sources.txt --in-sources-file ${dataDirPath}/chr${i}_${sourceCount}.edss
        #./${sopangExe} --ceds ${dataDirPath}/chr${i}.eds ${dataDirPath}/patterns16.txt -d -o ${outputFileName}_16_sources.txt --in-sources-file ${dataDirPath}/chr${i}_${sourceCount}.edss --full-sources-output

        #./${sopangExe} --ceds ${dataDirPath}/chr${i}.eds ${dataDirPath}/patterns32.txt -d -o ${outputFileName}_32_sources.txt --in-sources-file ${dataDirPath}/chr${i}_${sourceCount}.edss
        #./${sopangExe} --ceds ${dataDirPath}/chr${i}.eds ${dataDirPath}/patterns32.txt -d -o ${outputFileName}_32_sources.txt --in-sources-file ${dataDirPath}/chr${i}_${sourceCount}.edss --full-sources-output

        #./${sopangExe} --ceds ${dataDirPath}/chr${i}.eds ${dataDirPath}/patterns64.txt -d -o ${outputFileName}_64_sources.txt --in-sources-file ${dataDirPath}/chr${i}_${sourceCount}.edss
        #./${sopangExe} --ceds ${dataDirPath}/chr${i}.eds ${dataDirPath}/patterns64.txt -d -o ${outputFileName}_64_sources.txt --in-sources-file ${dataDirPath}/chr${i}_${sourceCount}.edss --full-sources-output
#    done
#done

# Running ED-text matching with sources.


for i in $(seq ${minChromosomeId} ${maxChromosomeId});
do
    for sourceCount in 32 64 128 256 512
    do
        #./${bndmExe} ${dataDirPath}/chr${i}.eds 1000 16 cedsm  ${dataDirPath}/chr${i}_${sourceCount}.edss > ./results/${outputFileName}_${i}_${sourceCount}_bndm16.txt
        ./${bndmExe} ${dataDirPath}/chr${i}.eds 1000 32 cedsm ATCTAATAATCTAATAATCTAATAATCTAATA ${dataDirPath}/chr${i}_${sourceCount}.edss ./results/${outputFileName}_bndm32.txt
        ./${sopangExe} --ceds ${dataDirPath}/chr${i}.eds sample/patterns32.txt -d -o ./results/${outputFileName}_32.txt --in-sources-file ${dataDirPath}/chr${i}_${sourceCount}.edss
    done
done