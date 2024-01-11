#!/bin/env python

def parseMakeDB(makedb_file):
    file_dict = dict()
    input_file = open(makedb_file, "r")
    header = input_file.readline().rstrip()
    header = header.split("\t")
    for row in input_file:
        line = row.rstrip("\n").split("\t")
        key = line[0]
        row_dict = dict(zip(header, line))
        file_dict[key] = row_dict

    return(header, file_dict)


def correct_sequences(file_dict):
    modified_dict = file_dict.copy()
    for key in modified_dict.keys():
        if 'SEQUENCE_IMGT' in modified_dict.keys():
            if (  (re.match("^IG[HK]", modified_dict[key]['V_CALL'] )) and  ( modified_dict[key]['J_CALL'] != ""  ) ):
                modified_dict[key]['FWR1_IMGT'] = modified_dict[key]['SEQUENCE_IMGT'][:78]
                modified_dict[key]['CDR1_IMGT'] = modified_dict[key]['SEQUENCE_IMGT'][78:114]
                modified_dict[key]['FWR2_IMGT'] = modified_dict[key]['SEQUENCE_IMGT'][114:165]
                modified_dict[key]['CDR2_IMGT'] = modified_dict[key]['SEQUENCE_IMGT'][165:195]
                modified_dict[key]['FWR3_IMGT'] = modified_dict[key]['SEQUENCE_IMGT'][195:312]
                #modified_dict[key]['CDR3_IMGT'] = modified_dict[key]['SEQUENCE_IMGT'][312:351]
                modified_dict[key]['CDR3_IMGT'] = modified_dict[key]['CDR3_IGBLAST']
            if (  (re.match("^IGL", modified_dict[key]['V_CALL'] )) and  ( modified_dict[key]['J_CALL'] != ""  ) ):
                modified_dict[key]['FWR1_IMGT'] = modified_dict[key]['SEQUENCE_IMGT'][:81]
                modified_dict[key]['CDR1_IMGT'] = modified_dict[key]['SEQUENCE_IMGT'][81:117]
                modified_dict[key]['FWR2_IMGT'] = modified_dict[key]['SEQUENCE_IMGT'][117:174]
                modified_dict[key]['CDR2_IMGT'] = modified_dict[key]['SEQUENCE_IMGT'][174:204]
                modified_dict[key]['FWR3_IMGT'] = modified_dict[key]['SEQUENCE_IMGT'][204:321]
                #modified_dict[key]['CDR3_IMGT'] = modified_dict[key]['SEQUENCE_IMGT'][321:357]
                modified_dict[key]['CDR3_IMGT'] = modified_dict[key]['CDR3_IGBLAST']
        elif 'sequence_alignment' in modified_dict.keys():
            if (  (re.match("^IG[HK]", modified_dict[key]['V_CALL'] )) and  ( modified_dict[key]['J_CALL'] != ""  ) ):
                modified_dict[key]['fwr1'] = modified_dict[key]['sequence_alignment'][:78]
                modified_dict[key]['cdr1'] = modified_dict[key]['sequence_alignment'][78:114]
                modified_dict[key]['fwr2'] = modified_dict[key]['sequence_alignment'][114:165]
                modified_dict[key]['cdr2'] = modified_dict[key]['sequence_alignment'][165:195]
                modified_dict[key]['fwr3'] = modified_dict[key]['sequence_alignment'][195:312]
                #modified_dict[key]['cdr3'] = modified_dict[key]['sequence_alignment'][312:351]
                modified_dict[key]['cdr3'] = modified_dict[key]['cdr3_IGBLAST']
            if (  (re.match("^IGL", modified_dict[key]['V_CALL'] )) and  ( modified_dict[key]['J_CALL'] != ""  ) ):
                modified_dict[key]['fwr1'] = modified_dict[key]['sequence_alignment'][:81]
                modified_dict[key]['cdr1'] = modified_dict[key]['sequence_alignment'][81:117]
                modified_dict[key]['fwr2'] = modified_dict[key]['sequence_alignment'][117:174]
                modified_dict[key]['cdr2'] = modified_dict[key]['sequence_alignment'][174:204]
                modified_dict[key]['fwr3'] = modified_dict[key]['sequence_alignment'][204:321]
                #modified_dict[key]['cdr3'] = modified_dict[key]['sequence_alignment'][321:357]
                modified_dict[key]['cdr3'] = modified_dict[key]['cdr3']
    return(modified_dict)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism", help = "Organism")
    parser.add_argument("--makedb", help = "MakeDB file")
    parser.add_argument("--outdir", help = "output dir")
    parser.add_argument("--analysisname", help = "analysis name")

    args = parser.parse_args()

    organism = args.organism
    makedb_file = args.makedb
    outdir = args.outdir
    analysisname = args.analysisname

    output_file = outdir + (analysisname+"_db-pass.tab")

    header, file_dict = parseMakeDB(makedb_file)

    export_table = open( output_file, "w" )
    export_table.write("\t".join(header) + "\n" )
    if organism != "rhesus_monkey":
        for key in file_dict.keys():
            col_list = []
            for col in header:
                col_list.append( file_dict[key][col] )

            export_table.write( "\t".join(col_list) + "\n" )
    else:
        corrected_dict = correct_sequences(file_dict)
        for key in corrected_dict.keys():
            col_list = []
            for col in header:
                col_list.append( corrected_dict[key][col] )

            export_table.write( "\t".join(col_list) + "\n" )
    export_table.close()


if __name__ == '__main__':
    import argparse
    import pandas as pd
    import re

    main()
