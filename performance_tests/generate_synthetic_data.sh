#!/bin/sh

# These paths should be relative to the root folder of this repository.
dataDirPath="data"
scriptsDirPath="scripts"

#additionalArgs="--alphabet ACGTN --degenerate-positions-factor 0.1 --max-segment-variants 10 --max-variant-size 10"
additionalArgs="--alphabet ACGTN --segment-count 500000 --max-segment-variants 10 --max-variant-size 5"
additionalArgs1="--alphabet ACGTN --segment-count 500000 --degenerate-positions-factor 0.1 --max-variant-size 10"
additionalArgs2="--alphabet ACGTN --segment-count 500000 --degenerate-positions-factor 0.04 --max-variant-size 5"
additionalArgs3="--alphabet ACGTN --segment-count 500000 --degenerate-positions-factor 0.01 --max-variant-size 5"

cd ..

mkdir -p ${dataDirPath}

#python3 ./${scriptsDirPath}/generate_synth_ed_text.py ${dataDirPath}/chr25.eds --degenerate-positions-factor 0.04 ${additionalArgs}
#python3 ./${scriptsDirPath}/generate_synth_ed_text.py ${dataDirPath}/chr26.eds --degenerate-positions-factor 0.06 ${additionalArgs}
#python3 ./${scriptsDirPath}/generate_synth_ed_text.py ${dataDirPath}/chr27.eds --degenerate-positions-factor 0.08 ${additionalArgs}
#python3 ./${scriptsDirPath}/generate_synth_ed_text.py ${dataDirPath}/chr28.eds --degenerate-positions-factor 0.10 ${additionalArgs}

python3 ./${scriptsDirPath}/generate_synth_ed_text.py ${dataDirPath}/chr25.eds --max-segment-variants 2 ${additionalArgs3}
python3 ./${scriptsDirPath}/generate_synth_ed_text.py ${dataDirPath}/chr26.eds --max-segment-variants 4 ${additionalArgs3}
python3 ./${scriptsDirPath}/generate_synth_ed_text.py ${dataDirPath}/chr27.eds --max-segment-variants 8 ${additionalArgs3}
python3 ./${scriptsDirPath}/generate_synth_ed_text.py ${dataDirPath}/chr28.eds --max-segment-variants 16 ${additionalArgs3}

for sourceCount in 32 64 128 256 512
do
    python3 ./${scriptsDirPath}/generate_synth_sources.py ${dataDirPath}/chr25.eds ${dataDirPath}/chr25_${sourceCount}.edss --source-count ${sourceCount}

    python3 ./${scriptsDirPath}/generate_synth_sources.py ${dataDirPath}/chr26.eds ${dataDirPath}/chr26_${sourceCount}.edss --source-count ${sourceCount}

    python3 ./${scriptsDirPath}/generate_synth_sources.py ${dataDirPath}/chr27.eds ${dataDirPath}/chr27_${sourceCount}.edss --source-count ${sourceCount}

    python3 ./${scriptsDirPath}/generate_synth_sources.py ${dataDirPath}/chr28.eds ${dataDirPath}/chr28_${sourceCount}.edss --source-count ${sourceCount}
done

#for sourceCount in 32 64 128 256 512
#do
#    python3 ./${scriptsDirPath}/generate_synth_sources.py ${dataDirPath}/chr25x.eds ${dataDirPath}/chr25x_${sourceCount}.edss --source-count ${sourceCount}
#    python3 ./${scriptsDirPath}/generate_synth_sources.py ${dataDirPath}/chr26x.eds ${dataDirPath}/chr26x_${sourceCount}.edss --source-count ${sourceCount}
#    python3 ./${scriptsDirPath}/generate_synth_sources.py ${dataDirPath}/chr27x.eds ${dataDirPath}/chr27x_${sourceCount}.edss --source-count ${sourceCount}
#    python3 ./${scriptsDirPath}/generate_synth_sources.py ${dataDirPath}/chr28x.eds ${dataDirPath}/chr28x_${sourceCount}.edss --source-count ${sourceCount}
#done

#for sourceCount in 32 64 128 256 512
#do
#    python3 ./${scriptsDirPath}/generate_synth_sources.py ${dataDirPath}/chr25y.eds ${dataDirPath}/chr25y_${sourceCount}.edss --source-count ${sourceCount}
#    python3 ./${scriptsDirPath}/generate_synth_sources.py ${dataDirPath}/chr26y.eds ${dataDirPath}/chr26y_${sourceCount}.edss --source-count ${sourceCount}
#    python3 ./${scriptsDirPath}/generate_synth_sources.py ${dataDirPath}/chr27y.eds ${dataDirPath}/chr27y_${sourceCount}.edss --source-count ${sourceCount}
#    python3 ./${scriptsDirPath}/generate_synth_sources.py ${dataDirPath}/chr28y.eds ${dataDirPath}/chr28y_${sourceCount}.edss --source-count ${sourceCount}
#done