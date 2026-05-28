#!/bin/sh

source ../.env

# if [ -f "accuracy_list.txt" ]; then 
#     rm accuracy_list.txt
# fi

if [ -f "params.txt" ]; then 
    rm params.txt
fi

while read geno_prob; do

    echo "Processing $geno_prob ..."

    echo "-i converted_wider_$geno_prob.txt" >> params.txt
    echo "-l 3" >> params.txt
    echo "-c 0" >> params.txt
    echo "-n SNP_no,SSSSS,PPPPPPP" >> params.txt

    python ../converter.py -input_file /Users/xtang3/Downloads/genoProbsforXinger/$geno_prob.geno_prob.txt \
                            -output_file converted_wider_$geno_prob.txt \
                            -params_file params.txt \
                            -program ../ImputeAccure.py
    
    echo "converted_wider_$geno_prob.accuracy" >> accuracy_list.txt

    rm params.txt
    
done <wider_geno_prob_list.txt

# while read true_geno_prob; do

#     echo "Processing $true_geno_prob ..."

#     echo "-i converted_$true_geno_prob.txt" >> params.txt
#     echo "-l 3" >> params.txt
#     echo "-c 0" >> params.txt
#     echo "-n SNP_no,SSSSS,PPPPPPP" >> params.txt

#     python ../converter.py -input_file $TRUE_DATA_DIR/$true_geno_prob.txt \
#                             -output_file converted_$true_geno_prob.txt \
#                             -params_file params.txt \
#                             -program ../ImputeAccure.py
    
#     # echo "converted_$true_geno_prob.accuracy" >> accuracy_list.txt

#     rm params.txt
    
# done <true_geno_prob_list.txt

