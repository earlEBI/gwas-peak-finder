#!/usr/bin/env python

import pandas as pd
import argparse
import os
import math
from tqdm import tqdm

# we don't want to get warnings. It works.
pd.options.mode.chained_assignment = None 

# Function to calculate -log10(x) without handling x as a float:
def get_log(pv):
    try:
        pv_str = str(pv).strip()
        if not pv_str:
            return float('nan')

        if 'E' in pv_str.upper():
            mantissa, exponent = pv_str.upper().split('E', 1)
            mantissa_float = float(mantissa)
            if mantissa_float == 0:
                return float('-inf')
            return math.log10(mantissa_float) + float(exponent)

        pv_float = float(pv_str)
        if pv_float == 0:
            return float('-inf')
        return math.log10(pv_float)
    except (TypeError, ValueError):
        return float('nan')

# Function to find and annotate top association:
def find_top_association(trait, chrom):
    global input_df
    global window
    global threshold

    test_df = input_df.loc[(input_df['trait'] == trait) & (input_df['chromosome'] == chrom)].copy()
    if test_df.empty:
        return

    # test_df.loc[:,'pvalue'] = test_df['pvalue'].astype(float) # changing type of p-value column <- too dangerous. Underflow hazard.
    test_df['mLogPv'] = test_df['pvalue'].apply(get_log)
    test_df['bp_location_num'] = pd.to_numeric(test_df['bp_location'], errors='coerce')

    # sorting rows by p-value:
    for index in test_df.sort_values(['mLogPv'], na_position='last').index:
        
        # Check if the index has 'false' or 'review' flag:
        if test_df.loc[index]['isTopAssociation'] == 'false' or \
            test_df.loc[index]['isTopAssociation'] == 'AUTO-REVIEWED-FALSE':
            continue

        # Excluding sub significant variants.
        elif pd.isna(test_df.loc[index]['mLogPv']) or test_df.loc[index]['mLogPv'] > threshold:
            test_df.loc[index, 'isTopAssociation'] = 'false'
            continue

        # We have found a top snp!
        else:
            pos = test_df.loc[index]['bp_location_num']
            if pd.isna(pos):
                test_df.loc[index, 'isTopAssociation'] = 'false'
                continue

            in_window = abs(test_df['bp_location_num'] - pos) <= window
            # Setting false flag for ALL variants within the window.
            test_df.loc[in_window, 'isTopAssociation'] = 'false'
            
            # If there are multiple variants with the same p-value, request review, othervise it's a true peak.
            same_pvalue_in_window = in_window & (test_df['mLogPv'] == test_df.loc[index, 'mLogPv'])
            if same_pvalue_in_window.sum() > 1:
                test_df.loc[same_pvalue_in_window, 'isTopAssociation'] = 'AUTO-REVIEWED-FALSE'
                similar = test_df.loc[same_pvalue_in_window]
                
                # Sort similar entries by bp location
                similar = similar.sort_values(['bp_location_num'], ascending=True)
                
                num_snps = len(similar)
                middle_index = (num_snps - 1) // 2
                actual_middle_index = similar.iloc[middle_index].name
                similar.loc[actual_middle_index, 'isTopAssociation'] = 'true'  # Set the middle or first middle to true
                test_df.update(similar)
            else:
                test_df.loc[index, 'isTopAssociation'] = 'true'


    # Modifying the original dataframe: 
    input_df.loc[ test_df.index, 'isTopAssociation'] = test_df.isTopAssociation

    

if __name__ == '__main__':

    # Parsing commandline arguments
    parser = argparse.ArgumentParser(description='This script finds the most significant association within a defined range (100kbp by default).')

    parser.add_argument('-f', '--input', help='Input file name with table of associations.')
    parser.add_argument('-o', '--output', help='Output file name.')
    parser.add_argument('-w', '--window', default=100000, help='Window size.', type = int)
    parser.add_argument('-t', '--threshold', default=1e-5, help='p-value threshold.', type = float)
    parser.add_argument('-p', '--prune', default=False, help='Prune out sub significant associations from the output.', action='store_true')

    args = parser.parse_args()

    # inputFile = args.input
    outputFile = args.output
    window = args.window
    if args.threshold <= 0 or args.threshold >= 1:
        raise ValueError("[Error] --threshold must be between 0 and 1 (exclusive).")
    threshold = math.log10(args.threshold)
    prune = args.prune

    if not outputFile:
        raise(Exception("[Error] A output file needs to be specified! Exiting."))

    # If inputfile is not specified or not submitted exiting:
    if args.input and os.path.isfile(args.input):
        inputFile = args.input
    else:
        raise(Exception("[Error] A valid input file is required! Exiting."))

    # Reading input file into pandas dataframe:
    try:
        input_df = pd.read_csv(inputFile, sep="\t", dtype=str)
    except UnicodeDecodeError:
        input_df = pd.read_csv(inputFile, sep="\t", dtype=str, encoding="latin1")
    except:
        raise Exception("[Error] Input file could not be read. Please check a .txt file is given or try pasting values without formatting to a fresh template. Exiting.")

    # Checking header (setting lowercase):
    input_df.columns = map(str.lower, input_df.columns)
    if not pd.Series(['trait', 'rs_id', 'pvalue', 'chromosome', 'bp_location']).isin(input_df.columns).all():
        raise Exception('[Error] Not all required columns were found in the file header. Required columns: "trait", "rs_id", "pvalue", "chromosome", "bp_location"')

    if prune:
        input_df['_mLogPv'] = input_df['pvalue'].apply(get_log)
        input_df = input_df[input_df['_mLogPv'] <= threshold].copy()
        input_df.drop(columns=['_mLogPv'], inplace=True)

    input_df['isTopAssociation'] = ''

    unique_pairs = input_df[['trait', 'chromosome']].drop_duplicates()
    for (trait, chromosome) in tqdm(unique_pairs.itertuples(index=False), total=unique_pairs.shape[0]):
        find_top_association(trait, chromosome)
        
    
    # Saving the modified table into a tab separated file:
    # input_df.drop(["mLogPv"], inplace=True)
    input_df.to_csv(outputFile, sep="\t", index= False, na_rep = 'NA')
