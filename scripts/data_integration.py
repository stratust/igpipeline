#!/usr/local/bin/python


def pprint_color(obj, flat=False):
    jsonpickle.set_preferred_backend('json')
    jsonpickle.set_encoder_options('json', sort_keys=True, indent=4)
    if flat is True:
        parsed = jsonpickle.encode(obj, unpicklable=False)
    else:
        parsed = jsonpickle.encode(obj, make_refs=True)

    print(
        highlight(
            parsed,
            JsonLexer(),
            Terminal256Formatter(style='rainbow_dash')
        )
    )


def get_dict_from_table(table, clones_dict, check_dict):

    table_file = open(table, "r")
    table_file_dict = dict()

    header = []
    for row in table_file:
        if re.match('^SEQUENCE_ID', row, re.IGNORECASE):
            header = row.rstrip().split("\t")
            continue

        if not header:
            print(header)
            print("No header in the file")
            sys.exit()

        row_list = row.rstrip().split("\t")
        row_dict = dict(zip(header, row_list))

        if check_dict:
            if row_list[0] in clones_dict:
                table_file_dict[row_list[0]] = row_dict
        else:
            table_file_dict[row_list[0]] = row_dict

    return(table_file_dict)


def get_sequences(igblast_airr_dict, v_germline_sequences, organism, hv_primer, kv_primer, lv_primer, corrected_regions_file_dict):

    header = [
        'full_input',
        'corrected_input',
        'full_input_from_start',
        'corrected_input_from_start'
    ]
    sequences_dict = dict()
    aux_dict = dict()

    for key in igblast_airr_dict.keys():

        full_input_from_start = ""
        corrected_input = ""
        corrected_input_from_start = ""

        full_input = igblast_airr_dict[key]['sequence']
        vdj_sequences = corrected_regions_file_dict[key]['SEQUENCE_VDJ']
        vdj_sequences = re.sub("-", "", vdj_sequences)

        full_input = re.match(r'(^\S*' + vdj_sequences + ')', full_input).group(1)

        fwr1_start = int(igblast_airr_dict[key]['v_sequence_start']) - 1
        v_germline_start = int(igblast_airr_dict[key]['v_germline_start']) - 1
        v_germline_id = igblast_airr_dict[key]['v_call'].split(",")[0]

        if re.search(r"IGH", v_germline_id):
            correction_length = int(hv_primer)

        elif re.search(r"IGK", v_germline_id):
            correction_length = int(kv_primer)

        elif re.search(r"IGL", v_germline_id):
            correction_length = int(lv_primer)


        v_germ_sequence = v_germline_sequences[v_germline_id].seq

        if fwr1_start <= v_germline_start:

            if v_germline_start > correction_length:
                from_start_nth_nt_germ_seq = v_germ_sequence[:v_germline_start]

                corrected_input_from_start = from_start_nth_nt_germ_seq + full_input[fwr1_start:]
                corrected_input = full_input
                full_input_from_start = corrected_input_from_start

            else:
                from_start_nth_nt_germ_seq = v_germ_sequence[:correction_length]

                full_input_end = (correction_length - v_germline_start) + fwr1_start
                relative_germline_start = correction_length - full_input_end
                germline_overlap_seq = from_start_nth_nt_germ_seq[relative_germline_start:]
                corrected_input = germline_overlap_seq + full_input[full_input_end :]
                corrected_input_from_start = from_start_nth_nt_germ_seq + full_input[full_input_end:]
                full_input_from_start = from_start_nth_nt_germ_seq[:relative_germline_start] + full_input

        elif fwr1_start > v_germline_start:

            if v_germline_start > correction_length:
                from_start_nth_nt_germ_seq = v_germ_sequence[:v_germline_start]

                corrected_input_from_start = from_start_nth_nt_germ_seq + full_input[fwr1_start : ]
                corrected_input = full_input[:fwr1_start - v_germline_start] + from_start_nth_nt_germ_seq[:v_germline_start] + full_input[fwr1_start: ]
                full_input_from_start = corrected_input

            else:
                from_start_nth_nt_germ_seq = v_germ_sequence[:correction_length]

                full_input_end = (correction_length - v_germline_start) + fwr1_start
                corrected_input_from_start = from_start_nth_nt_germ_seq + full_input[full_input_end :]
                corrected_input = full_input[: fwr1_start - v_germline_start ] + corrected_input_from_start
                full_input_from_start = full_input[: fwr1_start - v_germline_start ] + from_start_nth_nt_germ_seq[:v_germline_start] + full_input[fwr1_start:]

        sequences_list = [str(full_input), str(corrected_input), str(full_input_from_start), str(corrected_input_from_start)]
        aux_dict = dict(zip(header, sequences_list))
        sequences_dict[key] = aux_dict

    return(sequences_dict)


def check_dict_keys(igblast_dict):

    keys_to_check = ['CDR3-IMGT (germline)_from', 'CDR3-IMGT (germline)_to', 'CDR3-IMGT (germline)_length', 'CDR3-IMGT (germline)_matches', 'CDR3-IMGT (germline)_mismatches', 'CDR3-IMGT (germline)_gaps',
                    'FR1-IMGT_from', 'FR1-IMGT_to', 'FR1-IMGT_length', 'FR1-IMGT_matches', 'FR1-IMGT_mismatches', 'FR1-IMGT_gaps',
                    'CDR1-IMGT_from', 'CDR1-IMGT_to', 'CDR1-IMGT_length', 'CDR1-IMGT_matches', 'CDR1-IMGT_mismatches', 'CDR1-IMGT_gaps',
                    'FR2-IMGT_from', 'FR2-IMGT_to', 'FR2-IMGT_length', 'FR2-IMGT_matches', 'FR2-IMGT_mismatches', 'FR2-IMGT_gaps',
                    'CDR2-IMGT_from', 'CDR2-IMGT_to', 'CDR2-IMGT_length', 'CDR2-IMGT_matches', 'CDR2-IMGT_mismatches', 'CDR2-IMGT_gaps',
                    'FR3-IMGT_from', 'FR3-IMGT_to', 'FR3-IMGT_length', 'FR3-IMGT_matches', 'FR3-IMGT_mismatches', 'FR3-IMGT_gaps']

    for seq in igblast_dict:
        for key in keys_to_check:
            if key not in igblast_dict[seq]:
                igblast_dict[seq][key] = np.nan

    return(igblast_dict)


def get_dict_from_igblast_fmt7(clones_dict, igblast_fmt7):

    igblast_file = open(igblast_fmt7, "r")
    igblast_file_dict = dict()
    information_dict = dict()

    key = None
    header = []
    header_list = []
    information_all_regions = []

    information_flag = False
    for row in igblast_file:
        if re.match(".*Query: ", row):
            key = row.split(" ")[2].rstrip()
            continue

        if re.match(".*Alignment summary", row):
            header = re.search(r'\(.*\)', row).group(0)
            header = header.split(",")
            header = [element.strip() for element in header]
            header[0] = header[0].replace("(", "")
            header[-1] = header[-1].replace(")", "")
            header_aux = header
            information_flag = True
            continue

        if (re.match("^(?!Total)", row)) and (information_flag):
            information_list = row.rstrip().split("\t")
            region = information_list[0]
            header = [region + "_" + element for element in header]
            header_list.append(header)
            information_all_regions.append(information_list[1:])
            header = header_aux
            continue

        elif re.match("^Total\t", row):
            information_flag = False

            flat_header_list = [
                item for sublist in header_list for item in sublist
            ]

            flat_information_list = [
                item for sublist in information_all_regions for item in sublist
            ]

            information_dict = dict(
                zip(flat_header_list, flat_information_list)
            )

            header_list = []
            information_all_regions = []

        if key is not None and key in clones_dict:
            igblast_file_dict[key] = information_dict

    igblast_file_dict_corrected = check_dict_keys(igblast_file_dict)
    print("Correction:")
    print(igblast_file_dict_corrected)
    return(igblast_file_dict_corrected)


def hamming_distance(chaine1, chaine2):
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))


def aminoacids_mismatches(aminoacids_sequences_table):

    mismatches_list = []
    for i in range(0, aminoacids_sequences_table.shape[0]):
        v_germ_seq = str(
            aminoacids_sequences_table.iloc[i]['v_germline_alignment_aa']
        )
        v_seq_aa = str(
            aminoacids_sequences_table.iloc[i]['v_sequence_alignment_aa'])

        if len(v_germ_seq) > len(v_seq_aa):
            v_germ_seq_subset = v_germ_seq[:len(v_seq_aa)]
            mismatches_list.append(
                hamming_distance(
                    v_germ_seq_subset,
                    v_seq_aa))

        elif len(v_germ_seq) < len(v_seq_aa):
            v_seq_aa_subset = v_seq_aa[:len(v_germ_seq)]
            mismatches_list.append(
                hamming_distance(
                    v_germ_seq,
                    v_seq_aa_subset))

        elif len(v_germ_seq) == len(v_seq_aa):
            mismatches_list.append(hamming_distance(v_germ_seq, v_seq_aa))

    return(mismatches_list)


def select_information(define_clones_dict, igblast_airr_dict, igblast_fmt7_dict, corrected_sequences_dict, correction):

    define_clones_pd = pd.DataFrame(define_clones_dict).T

    igblast_airr_pd = pd.DataFrame(igblast_airr_dict).T
    igblast_airr_pd.insert(0, 'SEQUENCE_ID', list(igblast_airr_pd.index))

    igblast_fmt7_pd = pd.DataFrame(igblast_fmt7_dict).T
    igblast_fmt7_pd.insert(0, 'SEQUENCE_ID', list(igblast_fmt7_pd.index))

    corrected_sequences_pd = pd.DataFrame(corrected_sequences_dict).T
    corrected_sequences_pd.insert(0, 'SEQUENCE_ID', list(corrected_sequences_pd.index))

    merge_1 = pd.merge(left = define_clones_pd, right = igblast_airr_pd, left_on = 'SEQUENCE_ID', right_on = 'SEQUENCE_ID')
    merge_2 = pd.merge(left = igblast_fmt7_pd, right = corrected_sequences_pd, left_on = 'SEQUENCE_ID', right_on = 'SEQUENCE_ID')
    table_all_columns = pd.merge(left = merge_1, right = merge_2, left_on = 'SEQUENCE_ID', right_on='SEQUENCE_ID')
    table_all_columns.fillna(0)


    ##################################################################################################################################
    ################################################ Generating Informations #########################################################
    ##################################################################################################################################

    ############################################### Aminoacids count #################################################################

    get_columns = ['cdr3_aa','fwr1_aa','cdr1_aa','fwr2_aa','cdr2_aa','fwr3_aa']
    result_columns = ['cdr3_aa_length','fwr1_aa_length','cdr1_aa_length','fwr2_aa_length','cdr2_aa_length','fwr3_aa_length']
    length_list = []
    for i in range(0, len(get_columns)):
        current_column = get_columns[i]
        aa_list = table_all_columns[ current_column ]

        for aa_seq in aa_list:
            length_list.append( len(aa_seq) )

        table_all_columns.insert(0, result_columns[i] , length_list)
        length_list = []

    ################################################ Aminoacids mismatches calculation ##############################################
    aminoacids_sequences_table = table_all_columns[['v_germline_alignment_aa', 'v_sequence_alignment_aa']]
    aminoacids_list_mismatches = aminoacids_mismatches(aminoacids_sequences_table)
    table_all_columns.insert(0, 'v_region_aa_mismatches', aminoacids_list_mismatches)

    ############################################### V insertions and deletions ######################################################
    v_cigar = table_all_columns['v_cigar']
    v_insertions = []
    v_deletions = []

    for desc in v_cigar:
        insertions_match = re.search(r"\d+I", desc)
        deletions_match = re.search(r"\d+D", desc)

        if insertions_match:
            insertions_desc = insertions_match.group(0)
            inserstions_number = int(re.search("\d+", insertions_desc).group(0))
            v_insertions.append(inserstions_number)
        else:
            v_insertions.append(0)

        if deletions_match:
            deletions_desc = deletions_match.group(0)
            deletions_number = int(re.search(r"\d+", deletions_desc).group(0))
            v_deletions.append(deletions_number)
        else:
            v_deletions.append(0)

    table_all_columns.insert(0, 'v_insertions', v_insertions)
    table_all_columns.insert(0, 'v_deletions', v_deletions)


    ############################################### V region mismatches ######################################################
    v_composition = ['FR1-IMGT_mismatches','CDR1-IMGT_mismatches','FR2-IMGT_mismatches','CDR2-IMGT_mismatches','FR3-IMGT_mismatches']
    nt_mismatches_V_region = []
    for i in range(0, table_all_columns.shape[0] ):
        v_mismatches_list = list(table_all_columns.iloc[i][v_composition])
        print(v_mismatches_list)

        v_mismatches_list_int = [int(item) for item in v_mismatches_list if str(item) != 'nan']
        print(v_mismatches_list_int)
        #v_mismatches_list_int = [ int(mismatch) for mismatch in v_mismatches_list_filtered ]

        nt_mismatches_V_region.append( sum(v_mismatches_list_int) )
    table_all_columns.insert(0, 'nt_mismatches_V_region', nt_mismatches_V_region)


    ############################################### Information extraction #################################################
    if correction == "False":
        table_all_columns_filtered = table_all_columns[['CLONE', 'SEQUENCE_ID', 'FUNCTIONAL', 'IN_FRAME', 'STOP' ,'V_CALL', 'D_CALL', 'J_CALL',
                                                        'full_input','full_input_from_start', 'v_germline_start',
                                                        'fwr1_start', 'fwr1',
                                                        'v_insertions', 'v_deletions', 'nt_mismatches_V_region',
                                                        'v_region_aa_mismatches',
                                                        'CDR3_IMGT','CDR3-IMGT (germline)_length','CDR3-IMGT (germline)_matches','CDR3-IMGT (germline)_mismatches', 'CDR3-IMGT (germline)_gaps', 'cdr3_aa', 'cdr3_aa_length',
                                                        'FWR1_IMGT','FR1-IMGT_length','FR1-IMGT_matches','FR1-IMGT_mismatches', 'FR1-IMGT_gaps', 'fwr1_aa', 'fwr1_aa_length',
                                                        'CDR1_IMGT','CDR1-IMGT_length','CDR1-IMGT_matches','CDR1-IMGT_mismatches', 'CDR1-IMGT_gaps', 'cdr1_aa', 'cdr1_aa_length',
                                                        'FWR2_IMGT','FR2-IMGT_length','FR2-IMGT_matches','FR2-IMGT_mismatches', 'FR2-IMGT_gaps', 'fwr2_aa', 'fwr2_aa_length',
                                                        'CDR2_IMGT','CDR2-IMGT_length','CDR2-IMGT_matches','CDR2-IMGT_mismatches', 'CDR2-IMGT_gaps', 'cdr2_aa', 'cdr2_aa_length',
                                                        'FWR3_IMGT','FR3-IMGT_length','FR3-IMGT_matches','FR3-IMGT_mismatches', 'FR3-IMGT_gaps', 'fwr3_aa', 'fwr3_aa_length',
                                                        'JUNCTION',
                                                        'JUNCTION_LENGTH',
                                                        'germline_alignment'
                                                        ]]

    else:
        table_all_columns_filtered = table_all_columns[['CLONE', 'SEQUENCE_ID', 'FUNCTIONAL', 'IN_FRAME', 'STOP' ,'V_CALL', 'D_CALL', 'J_CALL',
                                                        'corrected_input','corrected_input_from_start','full_input','full_input_from_start', 'v_germline_start',
                                                        'fwr1_start', 'fwr1',
                                                        'v_insertions', 'v_deletions', 'nt_mismatches_V_region',
                                                        'v_region_aa_mismatches',
                                                        'CDR3_IMGT','CDR3-IMGT (germline)_length','CDR3-IMGT (germline)_matches','CDR3-IMGT (germline)_mismatches', 'CDR3-IMGT (germline)_gaps', 'cdr3_aa', 'cdr3_aa_length',
                                                        'FWR1_IMGT','FR1-IMGT_length','FR1-IMGT_matches','FR1-IMGT_mismatches', 'FR1-IMGT_gaps', 'fwr1_aa', 'fwr1_aa_length',
                                                        'CDR1_IMGT','CDR1-IMGT_length','CDR1-IMGT_matches','CDR1-IMGT_mismatches', 'CDR1-IMGT_gaps', 'cdr1_aa', 'cdr1_aa_length',
                                                        'FWR2_IMGT','FR2-IMGT_length','FR2-IMGT_matches','FR2-IMGT_mismatches', 'FR2-IMGT_gaps', 'fwr2_aa', 'fwr2_aa_length',
                                                        'CDR2_IMGT','CDR2-IMGT_length','CDR2-IMGT_matches','CDR2-IMGT_mismatches', 'CDR2-IMGT_gaps', 'cdr2_aa', 'cdr2_aa_length',
                                                        'FWR3_IMGT','FR3-IMGT_length','FR3-IMGT_matches','FR3-IMGT_mismatches', 'FR3-IMGT_gaps', 'fwr3_aa', 'fwr3_aa_length',
                                                        'JUNCTION',
                                                        'JUNCTION_LENGTH',
                                                        'germline_alignment'
                                                        ]]


    split_regions = ['V_CALL', 'D_CALL', 'J_CALL']
    for i in range(0, len(split_regions)):
        region_call = split_regions[i]
        region_call_all_columns_filtered = table_all_columns_filtered[region_call]
        region_call_all_columns_filtered = [ element.split(",")[0]  for element in region_call_all_columns_filtered ]
        table_all_columns_filtered[region_call] = region_call_all_columns_filtered

    return(table_all_columns_filtered)


def corrected_sequences_to_fasta(corrected_sequences_dict, output):

    seq_record_list = []
    for key in corrected_sequences_dict:
        record = SeqRecord( Seq( corrected_sequences_dict[key]["corrected_input"] ), id = key, description = ""  )
        seq_record_list.append(record)
    SeqIO.write(seq_record_list, os.path.dirname(output)+"/corrected_sequences.fasta", "fasta")



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--igblast_airr", help="Igblast output file | --outfmt19")
    parser.add_argument("--igblast_fmt7", help="Igblast output file | --outfmt7")
    parser.add_argument("--define_clones", help="DefineClones.py output")
    parser.add_argument("--corrected_regions_file", help = "Corrected output")
    parser.add_argument("--germline_ig_v_seq", help="Germline V sequences database")
    parser.add_argument("--organism", help="Organism")
    parser.add_argument("--HV_primer", help="Heavy (V) primer length")
    parser.add_argument("--KV_primer", help="Kappa (V) primer length")
    parser.add_argument("--LV_primer", help="Lambda (V) primer length")
    parser.add_argument("--correction", help="Boolean value indicating whether the sequences will be corrected")
    parser.add_argument("--output", help="Output directory/file")
    args = parser.parse_args()

    igblast_airr = args.igblast_airr
    igblast_fmt7 = args.igblast_fmt7
    define_clones = args.define_clones
    corrected_regions_file = args.corrected_regions_file
    germline_ig_v_seq = args.germline_ig_v_seq
    organism = args.organism
    hv_primer = args.HV_primer
    kv_primer = args.KV_primer
    lv_primer = args.LV_primer
    correction = args.correction
    output = args.output

    v_germline_sequences = SeqIO.to_dict(
        SeqIO.parse(open(germline_ig_v_seq), 'fasta')
    )

    define_clones_dict = get_dict_from_table(define_clones,{}, False)

    corrected_regions_file_dict = get_dict_from_table(corrected_regions_file, {}, False)

    igblast_airr_dict = get_dict_from_table(igblast_airr, define_clones_dict, True)

    igblast_fmt7_dict = get_dict_from_igblast_fmt7(define_clones_dict, igblast_fmt7)

    corrected_sequences_dict = get_sequences(igblast_airr_dict, v_germline_sequences, organism, hv_primer, kv_primer, lv_primer, corrected_regions_file_dict)
    corrected_sequences_to_fasta(corrected_sequences_dict, output)

    final_pd = select_information(
        define_clones_dict,
        igblast_airr_dict,
        igblast_fmt7_dict,
        corrected_sequences_dict,
        correction)
    final_pd.to_csv(output, sep="\t", header=True, index=False)

if __name__ == '__main__':
    import argparse
    import sys
    import re
    import os
    import pandas as pd
    import numpy as np
    from pygments import highlight
    from pygments.lexers import JsonLexer
    from pygments.formatters import Terminal256Formatter
    import jsonpickle
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    main()
