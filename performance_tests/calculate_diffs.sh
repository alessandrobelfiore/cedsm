#!/bin/sh

resultsDirPath="results"
scriptsDirPath="scripts"

cd ..

for sourceCount in 32 64 128 256 512
do
    python3 ./${scriptsDirPath}/count_results.py ${resultsDirPath}/results_25_${sourceCount}.txt ${resultsDirPath}/results_25_${sourceCount}_sopang.txt > perc_25_${sourceCount}.txt
done

for sourceCount in 32 64 128 256 512
do
    python3 ./${scriptsDirPath}/count_results.py ${resultsDirPath}/results_26_${sourceCount}.txt ${resultsDirPath}/results_26_${sourceCount}_sopang.txt > perc_26_${sourceCount}.txt
done

for sourceCount in 32 64 128 256 512
do
    python3 ./${scriptsDirPath}/count_results.py ${resultsDirPath}/results_27_${sourceCount}.txt ${resultsDirPath}/results_27_${sourceCount}_sopang.txt > perc_27_${sourceCount}.txt
done

for sourceCount in 32 64 128 256 512
do
    python3 ./${scriptsDirPath}/count_results.py ${resultsDirPath}/results_28_${sourceCount}.txt ${resultsDirPath}/results_28_${sourceCount}_sopang.txt > perc_28_${sourceCount}.txt
done