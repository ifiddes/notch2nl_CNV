#!/usr/bin/bash
batchSystem="singleMachine"
jobTree=".jobTree"
maxThreads="10"
log="test.log"
queries="data/queries.pickle"
genomes="GRCh37-lite"
tissue_types="10"
debug_cutoff="2"

export PATH=./sonLib/bin:./jobTree/bin:${PATH}
export PYTHONPATH=./:${PYTHONPATH}

if [ -d ${jobTree} ]; then
    rm -rf ${jobTree}
fi

if [ ! -e ${queries} ]; then
    python src/cgqueryHandler.py --genomes ${genomes} --tissue_types ${tissue_types} --debug_cutoff ${debug_cutoff}
fi

python src/main.py --key_file test.log --cgquery_file ${queries} --jobTree ${jobTree} --maxThreads ${maxThreads} --logLevel DEBUG &> ${log}
