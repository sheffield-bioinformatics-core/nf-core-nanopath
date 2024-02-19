#!/usr/bin/env python

import pandas as pd
from functools import reduce
import requests
import json
import numpy as np
import argparse
import logging

logger = logging.getLogger()

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--infile", help="Classification results file", type=str)
    parser.add_argument("--prefix", help="ID of the sample, usually barcode.", type=str)
    parser.add_argument("--outfile", help="Output file name.", type=str, default="rel_abundance")

    return parser.parse_args()

def get_taxname(tax_id,tax_level):
    tags = {"S": "species_name","G": "genus_name","F": "family_name","O":'order_name', "C": "class_name"}
    tax_level_tag = tags[tax_level]
    #Avoids pipeline crash due to "nan" classification output. Thanks to Qi-Maria from Github
    if str(tax_id) == "nan":
        tax_id = 1
    
    path = 'http://api.unipept.ugent.be/api/v1/taxonomy.json?input[]=' + str(int(tax_id)) + '&extra=true&names=true'
    complete_tax = requests.get(path).text

    #Checks for API correct response (field containing the tax name). Thanks to devinbrown from Github
    try:
        name = json.loads(complete_tax)[0][tax_level_tag]
        if name == "":
            name = json.loads(complete_tax)[0]["taxon_name"]
    except:
        logger.error("Error retrieving taxonomic name for tax_id {tax_id} at level {tax_level}.".format(tax_id=tax_id, tax_level=tax_level))
        name = str(int(tax_id))

    return name


def get_taxname_from_dmp(data, tax_id, tax_level):
    tags = {"S": "species","G": "genus","F": "family", "O": "order"}
    tax_level_tag = tags[tax_level]

    if str(tax_id) == "nan":
        name = 'unclassified'
    else:
        name = data.loc[data['taxid'] == tax_id, tax_level_tag].iloc[0]
        if type(name) != str:
            name = data.loc[data['taxid'] == tax_id, "name"].iloc[0]
            if type(name) != str:
                name = data.loc[data['taxid'] == tax_id, "sciname"].iloc[0]

    return name


def get_abundance_values(names,paths):
    dfs = []
    for name,path in zip(names,paths):
        data1 = pd.read_csv(path, index_col=False, sep=';').iloc[:,1:]

        total = sum(data1['reads_in_cluster'])
        rel_abundance=[]

        data=choose_classification(data1)

        for index,row in data.iterrows():
            rel_abundance.append(row['reads_in_cluster'] / total * 100)
            
        data['rel_abundance'] = rel_abundance
        dfs.append(pd.DataFrame({'taxid': data['taxid'], 'rel_abundance': rel_abundance, 'reads': data['reads_in_cluster']}))
        data.to_csv("" + name + "_nanoclust_out.txt")

    return dfs, data

def choose_classification(dataframe):
    print(dataframe)
    if len(dataframe.columns)>13:
        chosen_frame=[]
        classification_score={}
        for index, row in dataframe.iterrows():
            print(row['class_level'])
            if row['class_level']=="S":
                chosen_frame.append(row.iloc[:12].tolist())
            else:
                classification_score["kraken2"]=sum(row.notna()[8:12])
                classification_score["blast"]=sum(row.notna()[24:])
                classification_score["seqmatch"]=sum(row.notna()[16:20])
                choice=max(classification_score, key=classification_score.get)
                
                if choice == "kraken2":
                    chosen_frame.append(row.iloc[:12].tolist())
                elif choice == "seqmatch":
                    chosen_frame.append(row.iloc[np.r_[0:4,13:20]].tolist())
                else:
                    chosen_frame.append(row.iloc[np.r_[0:4,20:28]].tolist())
            
        logger.info("Choosing classification")
        logger.debug(chosen_frame)

        chosen_df=pd.DataFrame(chosen_frame, columns=['reads_in_cluster', 'used_for_consensus', 'reads_after_corr', 'draft_id', 'classifier_name', 'taxid', 'stat', 'name', 'species', 'genus', 'family', 'order'])
        logger.debug(len(chosen_df))
        logger.debug(chosen_df)

        return chosen_df
    else:
        return dataframe


def merge_abundance(dfs, data, tax_level):
    df_final = reduce(lambda left,right: pd.merge(left,right,on='taxid',how='outer').fillna(0), dfs)
    all_tax=[]
    for index, row in df_final.iterrows():
        try:
            all_tax.append(get_taxname_from_dmp(data, row["taxid"], tax_level))
        except:
            logger.error("Error getting taxonomic name for tax_id {} in merge_abundance.".format(row["taxid"]))
            all_tax.append(get_taxname(row["taxid"], tax_level))
    df_final["taxid"] = all_tax
    df_final_grp = df_final.groupby(["taxid"], as_index=False).sum()
    df_final_sorted = df_final_grp.sort_values(by='rel_abundance', ascending=False)
    return df_final_sorted


def get_abundance(names,paths,tax_level, outfile):
    if(not isinstance(paths, list)):
        paths = [paths]
        names = [names]

    dfs, data = get_abundance_values(names,paths)
    df_final_grp = merge_abundance(dfs, data, tax_level)
    df_final_grp.to_csv(outfile + "_"+ names[0] + "_" + tax_level + ".csv", index = False)


def main(args):

    for level in ["G", "S", "O", "F"]:
        get_abundance(args.prefix, args.infile, level, args.outfile)

if __name__=="__main__":
    args = parse_args()

    main(args)