#!/bin/sh

source ../.env

if [ -f "accuracy_list.txt" ]; then 
    rm accuracy_list.txt
fi

if [ -f "params.txt" ]; then 
    rm params.txt
fi

while read geno_prob; do

    echo "-i converted_$geno_prob.txt" >> params.txt
    echo "-l 3" >> params.txt
    echo "-c 1" >> params.txt
    echo "-n SNP_no,SSSSS,PPPPPPP" >> params.txt

    python ../converter.py -input_file $DATA_DIR/$geno_prob/.geno_prob.txt \
                            -output_file converted_$geno_prob.txt \
                            -params_file params.txt \
                            -program ../ImputeAccure.py
    
    echo "converted_$geno_prob.accuracy" >> accuracy_list.txt

    rm params.txt
    
done <geno_prob_list.txt