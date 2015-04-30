export PATH=./sonLib/bin:./jobTree/bin:${PATH}
export PYTHONPATH=./:${PYTHONPATH}

python scripts/overlay_raw_data.py --individual_graph output/b88da23c/b88da23c_individual_graph.pickle --mole_graph output/b88da23c/b88da23c_graph.pickle --sun_results output/b88da23c/sun_results.pickle --uuid b88da23c --out_dir output/b88da23c/ --individual_raw_data output/b88da23c/b88da23c.Individual.RawData.pickle  --mole_raw_data output/b88da23c/b88da23c.RawData.pickle 
