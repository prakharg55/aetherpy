#!/bin/sh

../bin/run_plot_model_results.py -var=3 3DALL_t021221_000000.bin -alt=200
../bin/run_plot_block_model_results.py -var=Temperature -alt=30 3DALL_20110320_010000.nc
