#!/usr/local/bin/python

def Insert_row_(row_number, df, row_value):

    df1 = df[0:row_number]
    df2 = df[row_number:]

    df1.loc[row_number]=row_value

    df_result = pd.concat([df1, df2])
    df_result.index = [*range(df_result.shape[0])]

    return df_result


def skip_rows(changeo_df_ordered):

    clone_ids = list(np.unique(changeo_df_ordered['clone_id']))

    for i in range(0, len(clone_ids)):
        clones_df = changeo_df_ordered[changeo_df_ordered['clone_id'] == clone_ids[i]]
        last_row = clones_df.index[len(clones_df) - 1] + 1
        changeo_df_ordered = Insert_row_(int(last_row), changeo_df_ordered, '')

    return changeo_df_ordered


def order_clones_by_count(changeo_df):

    count = [1] * len(changeo_df)
    aux = pd.concat([changeo_df['CLONE'], pd.DataFrame(count)], axis = 1)
    aux.columns = ['CLONE', 'count']
    counts_per_clone = aux.groupby('CLONE').sum()


    count_list = []
    for i in range(0, len(counts_per_clone)):
        count_list.append([counts_per_clone.iloc[i, 0] ] * counts_per_clone.iloc[i, 0] )

    count_list = list(itertools.chain(*count_list))

    changeo_df = changeo_df.join(pd.DataFrame(count_list))
    changeo_df.rename(columns={0: 'count'}, inplace=True)
    changeo_df = changeo_df.sort_values(['count', 'CLONE'], ascending = False)
    changeo_df.insert(0, 'clone_id', changeo_df['CLONE'])
    changeo_df.insert(1, 'num_seqs', changeo_df['count'])
    changeo_df = changeo_df.drop(columns=['count','CLONE'])
    changeo_df.index = range(len(changeo_df))

    return(changeo_df)


def correct_df_functional(changeo_df_functional_only):

    clone_ids = np.unique(changeo_df_functional_only['clone_id'])

    for clone_id in clone_ids:
        changeo_df_functional_only_subset = changeo_df_functional_only[changeo_df_functional_only['clone_id']  == clone_id]
        current_count = changeo_df_functional_only_subset.iloc[0, ]['num_seqs']
        count_by_clone_id = len(changeo_df_functional_only_subset)

        modified_clones = []
        if current_count != count_by_clone_id:
            indexes = changeo_df_functional_only_subset.index
            modified_clones.append(clone_id)

            for index in indexes:
                changeo_df_functional_only.iloc[index, changeo_df_functional_only.columns.get_loc('num_seqs')] = count_by_clone_id

    #changeo_df_functional_only.sort_values(by = ['num_seqs','clone_id'], ascending=False )
    changeo_df_functional_only = changeo_df_functional_only.sort_values(by = ['num_seqs','clone_id'], ascending=False )
    changeo_df_functional_only.index = [*range(changeo_df_functional_only.shape[0])]
    #print(changeo_df_functional_only[['clone_id','num_seqs']])
    return(changeo_df_functional_only, modified_clones)


def get_dict(changeo_df, workbook):
    changeo_dict = dict()
    clones_dict = dict()

    CLONE_COLORS = ['#FF0000','#FF6666','#FF6600','#FFCC99','#FFFF00','#FFFFCC','#008000','#00FF00', '#CCFFCC', '#00FFFF','#99CCFF','#0000FF', '#800080','#660099','#FF00FF','#CC99FF','#FF99CC','#E0E0E0']
    clone_ids = changeo_df['clone_id']
    indexes = np.unique(clone_ids, return_index=True)[1]
    clone_ids = [clone_ids[index] for index in sorted(indexes)]

    current_color = 0

    for i in clone_ids:
        size = len(changeo_df[changeo_df['clone_id'] == i])

        if size == 1:
            color = workbook.add_format({})
            merge_format = color
            clones_dict = {"color": color, "merge": merge_format}
        else:
            color = workbook.add_format({'bg_color': CLONE_COLORS[current_color]})
            chart_color = CLONE_COLORS[current_color]
            merge_format = workbook.add_format({
                'bold': 1,
                'border': 1,
                'align': 'center',
                'valign': 'vcenter',
                'fg_color': CLONE_COLORS[current_color]
            })
            current_color = current_color + 1
            clones_dict = {"color": color, "merge": merge_format, "chart_color": chart_color}

        changeo_dict[i] = clones_dict

        if current_color == len(CLONE_COLORS) - 1:
            current_color = 0

    return(changeo_dict)


def color_and_merge(workbook, worksheet, changeo_df, clone_ids, changeo_dict):

    header_format = workbook.add_format({'bold': 1})
    for head_col in range(0,len(changeo_df.columns)):
        worksheet.write(0, head_col ,changeo_df.columns[head_col], header_format)

    to_merge = dict()
    for j in range(0, changeo_df.shape[0]):
        current_id = changeo_df.iloc[j, ]['clone_id']
        if current_id in to_merge:
            to_merge[current_id].append(j + 1)
        else:
            to_merge[current_id] = [j + 1]

        if (current_id == ''):
            continue
        for k in range(0, changeo_df.shape[1]):
            value = changeo_df.iloc[j, k]
            if (pd.isnull(value)):
                worksheet.write(j + 1, k, "", changeo_dict[current_id]["color"])
            else:
                worksheet.write(j + 1, k, value , changeo_dict[current_id]["color"])

    for clone_id in clone_ids:
        clone_size = len(to_merge[clone_id])
        if clone_size > 1:
            worksheet.merge_range(to_merge[clone_id][0],0, to_merge[clone_id][-1],0, clone_id, changeo_dict[clone_id]["merge"])
            worksheet.merge_range(to_merge[clone_id][0],1, to_merge[clone_id][-1],1, clone_size, changeo_dict[clone_id]["merge"])

    return(workbook)


def separate_sequences_by_chain(chains_dict, changeo_df, workbook):


    for key in chains_dict.keys():
        changeo_df_subset = None

        worksheet_name = key
        filtered_worksheet_name = key + " - filtered"

        if key == "heavy":
            if chains_dict[key]:
                changeo_df_subset = changeo_df[changeo_df['V_CALL'].str.contains('IGH.*')]

        if key == "kappa":
            if chains_dict[key]:
                changeo_df_subset = changeo_df[changeo_df['V_CALL'].str.contains('IGK.*')]

        if key == "lambda":
            if chains_dict[key]:
                changeo_df_subset = changeo_df[changeo_df['V_CALL'].str.contains('IGL.*')]

        if changeo_df_subset is not None:
            changeo_df_subset.index = [*range(changeo_df_subset.shape[0])]
            worksheet = workbook.add_worksheet(worksheet_name)

            changeo_df_subset_ordered = order_clones_by_count(changeo_df_subset)
            changeo_dict = get_dict(changeo_df_subset_ordered, workbook)
            clone_ids = list(np.unique(changeo_df_subset_ordered['clone_id']))
            changeo_df_subset_ordered_sep = skip_rows(changeo_df_subset_ordered)

            workbook = color_and_merge(workbook, worksheet, changeo_df_subset_ordered_sep, clone_ids, changeo_dict)

            # merge and coloring functional sequences only
            print(changeo_df_subset_ordered_sep)
            changeo_df_functional_only = changeo_df_subset_ordered[changeo_df_subset_ordered['FUNCTIONAL'] == 'T' ]
            changeo_df_functional_only.index = [*range(changeo_df_functional_only.shape[0])]

            if changeo_df_functional_only.shape[0] > 0:
                changeo_df_functional_only, modified_clones_list = correct_df_functional(changeo_df_functional_only)
                worksheet_functional_only = workbook.add_worksheet(filtered_worksheet_name)
                #changeo_functional_dict = get_dict(changeo_df_functional_only, workbook)
                clone_ids_functional = np.unique(changeo_df_functional_only['clone_id'])
                changeo_df_functional_only = skip_rows(changeo_df_functional_only)
                workbook = color_and_merge(workbook, worksheet_functional_only, changeo_df_functional_only, clone_ids_functional, changeo_dict)


    return(workbook)


def generating_plot(changeo_df_subset_ordered, workbook, worksheet_name, key, changeo_dict):


    worksheet = workbook.add_worksheet(worksheet_name)
    #changeo_dict = get_dict(changeo_df_subset_ordered, workbook)

    singles = len(np.where(changeo_df_subset_ordered['num_seqs'] == 1)[0])
    clones = list(np.where(changeo_df_subset_ordered['num_seqs'] != 1)[0])

    changeo_df_ordered_clones = changeo_df_subset_ordered.iloc[clones]
    clone_ids = changeo_df_ordered_clones['clone_id']
    indexes = np.unique(clone_ids, return_index=True)[1]
    clone_ids = [clone_ids[index] for index in sorted(indexes)]

    worksheet.write(0, 0, "clone_id")
    worksheet.write(0, 1, "num_seqs")
    clones_size = dict()
    color_points = []
    i = 1
    total = 0
    for clone_id in clone_ids:
        size = changeo_df_ordered_clones[changeo_df_ordered_clones['clone_id']  == clone_id].shape[0]
        clones_size[clone_id] = size
        total = total + size
        color_points.append({'fill': { 'color': changeo_dict[clone_id]["chart_color"]}, 'border': {'color': '#606060'}})
        worksheet.write(i, 0, clone_id)
        worksheet.write(i, 1, changeo_df_ordered_clones[changeo_df_ordered_clones['clone_id'] == clone_id].shape[0])
        i = i + 1
    worksheet.write(i, 0, "singles")
    worksheet.write(i, 1, singles)

    color_points.append({'fill': {'color': 'white'}, 'border': {'color': '#606060'}})

    chart1 = workbook.add_chart({'type': 'doughnut'})
    chart1.add_series({
        'categories': [worksheet_name, 1, 0, i, 0],
        'values': [worksheet_name, 1, 1, i, 1],
        'points': color_points
    })
    chart1.set_title({'name': key + ' chain clones'})
    chart1.set_style(10)

    worksheet.insert_chart('C2', chart1, {'x_offset': 15, 'y_offset': 7})
    worksheet.insert_textbox('F10',
                             str(total + singles),
                            {'font': {'size': 12},
                                'width': 60,
                                'height': 30,
                                'border': {'none': True},
                                'align': {'vertical': 'middle',
                                        'horizontal': 'center'}})
    return(workbook)


def plotting_by_chain(chains_dict, changeo_df, workbook):

    for key in chains_dict.keys():

        changeo_df_subset = None
        if key == "heavy":
            if chains_dict[key]:
                changeo_df_subset = changeo_df[changeo_df['V_CALL'].str.contains('IGH.*')]

        elif key == "kappa":
            if chains_dict[key]:
                changeo_df_subset = changeo_df[changeo_df['V_CALL'].str.contains('IGK.*')]

        elif key == "lambda":
            if chains_dict[key]:
                changeo_df_subset = changeo_df[changeo_df['V_CALL'].str.contains('IGL.*')]


        if changeo_df_subset is not None:
            worksheet_name = key + " - chart"
            filtered_worksheet_name = key + " filtered - chart"

            changeo_df_subset.index = [*range(changeo_df_subset.shape[0])]
            changeo_df_subset_ordered = order_clones_by_count(changeo_df_subset)
            changeo_dict = get_dict(changeo_df_subset_ordered, workbook)

            workbook = generating_plot(changeo_df_subset_ordered, workbook, worksheet_name, key, changeo_dict)

            changeo_df_functional_only = changeo_df_subset_ordered[changeo_df_subset_ordered['FUNCTIONAL'] == 'T' ]
            changeo_df_functional_only.index = [*range(changeo_df_functional_only.shape[0])]

            if changeo_df_functional_only.shape[0] > 0:
                changeo_df_functional_only, modified_clones_list = correct_df_functional(changeo_df_functional_only)
                workbook = generating_plot(changeo_df_functional_only, workbook, filtered_worksheet_name, key, changeo_dict)

    return(workbook)



def from_tab_to_xlsx(changeo_integrated_file, output):

    print(changeo_integrated_file + "\n")
    #output_file = re.sub('integration_output', 'xlsx_output', changeo_output)
    #output_file = re.sub('.tsv', '_summary.xlsx', output_file)

    #writer = pd.ExcelWriter(output, engine='xlsxwriter')
    changeo_df = pd.read_csv(changeo_integrated_file, sep='\t', header=0)
    workbook = xlsxwriter.Workbook(output)

    chains_dict = {
        "heavy": re.search('IGH.*', changeo_df['V_CALL'].to_string() ),
        "kappa": re.search('IGK.*', changeo_df['V_CALL'].to_string() ),
        "lambda": re.search('IGL.*', changeo_df['V_CALL'].to_string() )
    }


    workbook = separate_sequences_by_chain(chains_dict, changeo_df, workbook)
    workbook = plotting_by_chain(chains_dict, changeo_df, workbook)
    workbook.worksheets_objs.sort(key=lambda x: x.name)

    workbook.close()

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--changeo_integrated_file", help = "Changeo integrated file")
    parser.add_argument("--output", help="output directory/file")
    args = parser.parse_args()

    changeo_integrated_file = args.changeo_integrated_file
    output = args.output

    from_tab_to_xlsx(changeo_integrated_file, output)


if __name__ == '__main__':
    import sys
    import argparse
    import pandas as pd
    import re
    import itertools
    import numpy as np
    import xlsxwriter

    main()
