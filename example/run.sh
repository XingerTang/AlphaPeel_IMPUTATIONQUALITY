#!/bin/sh

python ../converter.py -input_file .geno_prob.txt \
                        -output_file converted_geno_prob.txt \
                        -params_file params.txt \
                        -program ../ImputeAccure.py
