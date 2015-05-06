export PATH=./sonLib/bin:./jobTree/bin:${PATH}
export PYTHONPATH=./:${PYTHONPATH}

python scripts/overlay_raw_data.py --individual_graph output/mike_sny/mike_sny_individual_graph.pickle \
--mole_graph output/mike_sny/mike_sny_graph.pickle --sun_results output/mike_sny/sun_results.pickle \
--uuid mike_sny --out_dir output/ --individual_raw_data output/mike_sny/mike_sny.Individual.RawData.pickle \
--mole_raw_data output/mike_sny/mike_sny.RawData.pickle 

python scripts/overlay_raw_data_v2.py --individual_graph output/mike_sny/mike_sny_individual_graph.pickle \
--mole_graph output/mike_sny/mike_sny_graph.pickle --sun_results output/mike_sny/sun_results.pickle \
--uuid mike_sny --out_dir output/ --individual_raw_data output/mike_sny/mike_sny.Individual.RawData.pickle \
--mole_raw_data output/mike_sny/mike_sny.RawData.pickle 

