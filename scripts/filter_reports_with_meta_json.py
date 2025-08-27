import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import json
import math
from Bio.Seq import Seq
from collections import defaultdict, Counter

def parseReport(report):

    prefix = report.split("/")[-1].split(".")[0]

    seqtups = []
    fragtups = []
    with open(report, 'r') as fin:
        lines = fin.read().split("Sequence variants:")[-1].split("\nLength variants:")[0]
        for line in lines.split("\n"):
            if line == "":
                continue
            seqtup, count = line.split("\t")
            _, seq, _, family, _ = seqtup.split("'")
            count = int(count)

            fullseq = str(Seq(seq).reverse_complement())
            
            subtype = inferGP60Subtype(fullseq, family)

            seqtups.append((fullseq, subtype, count))
    
    return sorted(seqtups, key=lambda x: x[1], reverse=True)


def filterAlleles(sid, seqtups, min_cov, min_af=0.05, stutter_filter=True, l3=0.16, l6=0.08):
    # Handle empty input
    if not seqtups:
        return []

    filtered_seqtups = []
    # Sort alleles by count (coverage) in descending order
    s_alleles = sorted(seqtups, key=lambda x: x[-1], reverse=True)
    major_allele = s_alleles[0]

    # Calculate total coverage to determine allele frequencies
    total_cov = sum(count for _, _, count in seqtups)
    if total_cov == 0:  # Avoid division by zero if counts are all zero
        return []

    primary_cov = major_allele[-1]
    lmaj = len(major_allele[0])
    
    maj_fam = major_allele[1][-3:]

    # The major allele is always kept, assuming it represents the primary infection
    filtered_seqtups.append(major_allele)

    # If minor alleles exist, filter them
    if len(s_alleles) > 1:
        minor_alleles = s_alleles[1:]

        for seq, subtype, count in minor_alleles:
            # --- General Filters for all minor alleles ---
            
            ## handle polyfamily
            if subtype[-3:] != maj_fam and count >= min_cov:
                # print(f"{sid} poly family: ({maj_fam}) {subtype[-3:]} {count}")
                filtered_seqtups.append((seq, subtype, count))
                continue

            # 1. Filter by minimum coverage
            if count < min_cov:
                continue

            # 2. Filter by minimum allele frequency
            af = count / total_cov
            if af < min_af:
                continue

            # --- Stutter-specific Filtering ---
            if stutter_filter:
                lmin = len(seq)
                is_stutter_3 = (lmin == lmaj - 3)
                is_stutter_6 = (lmin == lmaj - 6)

                if is_stutter_3 or is_stutter_6:
                    # Calculate frequency relative to the major allele for stutter check
                    primary_af = count / primary_cov

                    if is_stutter_3 and primary_af < l3:
                        continue  # Filtered as stutter artifact
                    if is_stutter_6 and primary_af < l6:
                        continue  # Filtered as stutter artifact

            # If the allele has passed all applicable filters, add it to the list
            filtered_seqtups.append((seq, subtype, count))

    return filtered_seqtups


def inferGP60Subtype(seq, family):

    # family = "IIa"
    subtype=f""

    if family == "IId":
        seq = seq[:-1]  ## adjust for janky OBO probe selection

        r_index = seq.find("ACATCG")
        seq = seq[:r_index]
        
        numTCA = seq.count("TCA")
        numTCG = seq.count("TCG")
        numTCT = seq.count("TCT")
        
        if (numTCA!=0):
            subtype += "A" + str(numTCA)

        if (numTCG!=0):
            subtype += "G" + str(numTCG)

        if (numTCT!=0):
            subtype += "T" + str(numTCT)

        return subtype + "-" + family

    elif family == "IIa" or family == "IIb" or family == "IIc" or family == "IIr":

        r_index = seq.find("ACATCA")

        pre_r_seq = seq[:r_index]

        numAcatca = seq.count("ACATCA")
        numTCA = pre_r_seq.count("TCA")
        numTCG = pre_r_seq.count("TCG")
        numTCT = pre_r_seq.count("TCT")

        if (numTCA!=0):
            subtype += "A" + str(numTCA)

        if (numTCG!=0):
            subtype += "G" + str(numTCG)

        if (numAcatca!=0):
            subtype += "R" + str(numAcatca)

        return subtype + "-" + family


def filter_only_mixed(moi_dict):

    filtered_dict = defaultdict(dict)

    for sample, allele_dict in moi_dict.items():
        if len(allele_dict) > 1:
            filtered_dict[sample] = allele_dict
    
    return filtered_dict


def read_meta_csv():

    meta_csv = "/home/arthur/BioInf/Crypto_popgen/data/crypto_meta_new2025.csv"

    meta_df = pd.read_csv(meta_csv)

    return {row['name'] : row['Host'] for index, row in meta_df.iterrows()}, {row['name'] : row['continent'] for index, row in meta_df.iterrows()}, {row['name'] : row['country'] for index, row in meta_df.iterrows()}


def run(poly_json, af, mincov, outdir="./", l3=0.15, l6=0.75):

    print(f"AF: {af}, cov: {mincov}, l3: {l3}, l6: {l6}")

    host_dict, continent_dict, country_dict = read_meta_csv()

    fails = []

    with open(poly_json, 'r') as fin:
        jdata = json.load(fin)
    
    moi_dict = defaultdict(dict)
    sample_ids = set([x.split("/")[-1].split(".")[0] for x in jdata.keys()])

    cov_dict = defaultdict(dict)

    for sid, seqtups in jdata.items():
        if len(seqtups) == 0:
            fails.append(sid)
            continue
        
        # cov_dict[sid][locus_id] = int(np.sum([ count for _, _, count in seqtups ]))
        # print(prefix)
        seqtups = filterAlleles(sid, seqtups, min_af=af, min_cov=mincov, stutter_filter=True, l3=l3, l6=l6)

        if len(seqtups) == 0:
            continue

        moi_dict[sid] = seqtups

    with open(f"{outdir}/filtered_gp60_MOI_{l3}_{l6}_AF{round(af,2)}_cov{mincov}.json", 'w') as f:
        json.dump(moi_dict, f, indent=4)

    pd.DataFrame(cov_dict).transpose().to_csv("BM_cov.csv")

    filtered_dict = filter_only_mixed(moi_dict)

    with open(f"{outdir}/filtered_gp60_MOI_{l3}_{l6}_AF{round(af,2)}_cov{mincov}_polyclonal.json", 'w') as f:
        json.dump(filtered_dict, f, indent=4)

    # unique_alleles = set()
    # for key, value in filtered_dict.items():
    #     for seq, subtype, count in value:
    #         unique_alleles.add(subtype)

    si_data = []
    for sid, alleles in moi_dict.items():
        total = sum([c for _, _, c in alleles])
        allele_ids = [aid for _, aid, g in sorted(alleles, key=lambda x: x[2], reverse=True)]
        major_allele = allele_ids[0]

        if len(allele_ids) > 1:
            minor_alleles = "/".join(allele_ids[1:])
        else:
            minor_alleles = None

        H = -(sum([(a/total)*math.log(a/total) for _, _, a in alleles]))

        if host_dict[sid] == "Homo sapiens":
            host_cat1 = "human"
        elif host_dict[sid] == "Bos taurus" or host_dict[sid] == "Ovis aries":
            host_cat1 = "livestock"
        else:
            host_cat1 = None

        if host_dict[sid] == "Homo sapiens":
            host_cat2 = "human"
        elif host_dict[sid] == "Bos taurus":
            host_cat2 = "cattle"
        else:
            host_cat2 = None

        si_data.append([sid, len(alleles), H, np.exp(H), host_dict[sid], continent_dict[sid], country_dict[sid], host_cat1, host_cat2, major_allele, minor_alleles])

    df = pd.DataFrame(si_data, columns=["Sample_ID", "num_alleles", "Shannon_index", "eMOI", "host", "continent", "country", "host_cat1", "host_cat2", "major_allele", "minor_alleles"])
    df.to_csv(f"{outdir}/gp60_MOI_{l3}_{l6}_AF{round(af,2)}_cov{mincov}.csv")


def main(poly_json):
    
    # l3 = 0.15
    # l6 = l3/2
    # af = 0.05
    # mincov = 3

    # run(rdir, af, mincov, outdir="./", l3=0, l6=0)


    # for a in range(1, 31, 1):
    #     ## AF 0.01-0.3
    #     af = a/100

    #     for l in range(1, 31):
    #         ## stutter filter l-3 0.01-0.3
    #         l3 = l/100
    #         l6 = l3/2

    #         run(rdir, af, 3, outdir="./stutter_filter_range", l3=l3, l6=l6)

    # for a in range(1, 31, 1):
    #     ## AF 0.01-0.3
    #     af = a/100

    #     for l in range(1, 31):
    #         ## mincov 1-30
    #         mincov = l

    #         run(rdir, af, mincov, outdir="./p_range_jsons", l3=0, l6=0)

    for l in range(1, 31, 1):
        ## stutter filter l-3 0.01-0.3
        l3 = l/100
        l6 = l3/2
        af = 0.05

        for c in range(1, 31):
            ## cov 1-30
            run(poly_json, af, c, outdir="./stutter_filter_c_vs_af/", l3=l3, l6=l6)

if __name__=="__main__":
    main(sys.argv[1])