#!/bin/sh

source ../.env

if [ -f "accuracy_list.txt" ]; then 
    rm accuracy_list.txt
fi

if [ -f "params.txt" ]; then 
    rm params.txt
fi

while read geno_prob; do

    echo "Processing $geno_prob ..."

    echo "-i converted_$geno_prob.txt" >> params.txt
    echo "-l 3" >> params.txt
    echo "-c 0" >> params.txt
    echo "-n SNP_no,SSSSS,PPPPPPP" >> params.txt

    python ../converter.py -input_file $DATA_DIR/$geno_prob/.geno_prob.txt \
                            -output_file converted_$geno_prob.txt \
                            -params_file params.txt \
                            -program ../ImputeAccure.py
    
    echo "converted_$geno_prob.accuracy" >> accuracy_list.txt

    rm params.txt
    
done <geno_prob_list.txt

while read true_geno_prob; do

    echo "Processing $true_geno_prob ..."

    echo "-i converted_$true_geno_prob.txt" >> params.txt
    echo "-l 3" >> params.txt
    echo "-c 0" >> params.txt
    echo "-n SNP_no,SSSSS,PPPPPPP" >> params.txt

    python ../converter.py -input_file $TRUE_DATA_DIR/$true_geno_prob.txt \
                            -output_file converted_$true_geno_prob.txt \
                            -params_file params.txt \
                            -program ../ImputeAccure.py
    
    echo "converted_$true_geno_prob.accuracy" >> accuracy_list.txt

    rm params.txt
    
done <true_geno_prob_list.txt


while read sex_geno_prob; do

    echo "Processing $sex_geno_prob female only..."

    echo "-i converted_${sex_geno_prob}_female.txt" >> params.txt
    echo "-l 3" >> params.txt
    echo "-c 0" >> params.txt
    echo "-n SNP_no,SSSSS,PPPPPPP" >> params.txt

    python converter_female_only.py -input_file $DATA_DIR/$sex_geno_prob/.geno_prob.txt \
                            -output_file converted_${sex_geno_prob}_female.txt \
                            -ped_file ${TRUE_DATA_DIR}/X_chr_ped_file.txt \
                            -params_file params.txt \
                            -program ../ImputeAccure.py
    
    echo "converted_${sex_geno_prob}_female.accuracy" >> accuracy_list.txt

    rm params.txt
    
done <sex_geno_prob_list.txt


while read sex_geno_prob; do

    echo "Processing $sex_geno_prob male only..."

    echo "-i converted_${sex_geno_prob}_male_switched.txt" >> params.txt
    echo "-l 3" >> params.txt
    echo "-c 0" >> params.txt
    echo "-n SNP_no,SSSSS,PPPPPPP" >> params.txt

    python converter_male_only.py -input_file $DATA_DIR/$sex_geno_prob/.geno_prob.txt \
                            -output_file converted_${sex_geno_prob}_male_switched.txt \
                            -ped_file ${TRUE_DATA_DIR}/X_chr_ped_file.txt \
                            -params_file params.txt \
                            -program ../ImputeAccure.py \
                            -switched 1
    
    echo "converted_${sex_geno_prob}_male_switched.accuracy" >> accuracy_list.txt

    rm params.txt
    
done <sex_geno_prob_list.txt
