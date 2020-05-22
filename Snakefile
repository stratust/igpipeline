'''
Antibody Analysis Snakefile
'''
import os, re, glob

CWD=os.getcwd()

configfile: CWD+'/config.yml'

cmd_result = os.popen("  for i in $( find " + CWD + "/data -iname \*.ab1 );do test=$(dirname $i); echo $test;done | uniq | perl -plane 's/data/results/g; s/^(.*\/(.*))$/$1\/igpairs\/$2\_strict\.xlsx/g' ").read().rstrip()
cmd_result_splitted = cmd_result.split("\n")

RESULT_FILE = cmd_result_splitted

correct = config['CORRECTION']
bcell_database = config['BCR_DB']


rule all:
    input: [RESULT_FILE, CWD + "/results/circos_plots/plots/", CWD + "/results/SHM_analysis/", CWD + "/results/hydrophobicity_analysis/"]


rule chromatogram_to_fasta:
    input: CWD + "/data/{directory}"
    output: CWD + "/results/{directory}/fasta_output/fasta_conversion.fasta"
    log: CWD + "/results/{directory}/{directory}.log"
    params:
        result_dir = CWD + "/results/{directory}/fasta_output"
    shell:
        """
        scripts/prepare_data.sh {input} {params.result_dir}/fasta_conversion.fasta
        """


rule fix_fasta:
    input: CWD + "/results/{directory}/fasta_output/fasta_conversion.fasta"
    output: CWD + "/results/{directory}/fasta_output/{directory}.fasta"
    log: CWD + "/results/{directory}/{directory}.log"
    params:
        organism = config['organism'],
        result_dir = CWD + "/results/{directory}/fasta_output"
    shell:
        """
        if [ '{params.organism}' = 'rhesus_monkey' ]; then
            scripts/monkey_ig_tookit.pl fix_heavy_chain --input_file {params.result_dir}/fasta_conversion.fasta --output_file {output}
        else
            mv {input} {output}
        fi
        """


rule igblast_one:
    input: CWD + "/results/{directory}/fasta_output/{directory}.fasta"
    output: CWD + "/results/{directory}/run1/igblast_output/{directory}.fmt7"
    log: CWD + "/results/{directory}/run1/igblast_output/{directory}.log"
    params:
        organism = config['organism'],
        germline_db_V = config['igblast_germline_db_V'],
        germline_db_D = config['igblast_germline_db_D'],
        germline_db_J = config['igblast_germline_db_J'],
        auxiliary_data = config['igblast_auxiliary_data']
    shell: 
        """
        igblastn -germline_db_V {params.germline_db_V} \
            -germline_db_D {params.germline_db_D} \
            -germline_db_J {params.germline_db_J} \
            -auxiliary_data {params.auxiliary_data} \
            -domain_system imgt \
            -ig_seqtype Ig \
            -organism {params.organism} \
            -outfmt '7 std qseq sseq btop' \
            -query {input} \
            -out {output} &> {log}
        """


rule igblast_airr_one:
    input: CWD + "/results/{directory}/fasta_output/{directory}.fasta"
    output: CWD + "/results/{directory}/run1/igblast_output/{directory}.fmt19"
    log: CWD + "/results/{directory}/run1/igblast_output/{directory}.log"
    params:
        organism = config['organism'],
        germline_db_V = config['igblast_germline_db_V'],
        germline_db_D = config['igblast_germline_db_D'],
        germline_db_J = config['igblast_germline_db_J'],
        auxiliary_data = config['igblast_auxiliary_data']
    shell:
        """
        igblastn -germline_db_V {params.germline_db_V} \
            -germline_db_D {params.germline_db_D} \
            -germline_db_J {params.germline_db_J} \
            -auxiliary_data {params.auxiliary_data} \
            -domain_system imgt \
            -ig_seqtype Ig \
            -organism {params.organism} \
            -outfmt 19 \
            -query {input} \
            -out {output} &> {log}
        """


rule makedb_one:
    input:
        igblast_output = CWD + "/results/{directory}/run1/igblast_output/{directory}.fmt7",
        fasta_file = CWD + "/results/{directory}/fasta_output/{directory}.fasta"
    output: CWD + "/results/{directory}/run1/makedb_output/{directory}_db-pass.tab"
    log: CWD + "/results/{directory}/run1/makedb_output/{directory}_db-pass.log"
    params:
        imgt_ig_v = config['imgt_germline_db_V'],
        imgt_ig_d = config['imgt_germline_db_D'],
        imgt_ig_j = config['imgt_germline_db_J'],
        outdir = CWD + "/results/{directory}/run1/makedb_output"
    shell:
        """
        MakeDb.py igblast -i {input.igblast_output} -s {input.fasta_file} -r {params.imgt_ig_v} {params.imgt_ig_d} {params.imgt_ig_j} --outdir {params.outdir} --regions --scores --failed --partial --cdr3 &> {log}
        """


rule correct_regions_one:
    input: CWD + "/results/{directory}/run1/makedb_output/{directory}_db-pass.tab"
    output: CWD + "/results/{directory}/run1/corrected_output/{directory}_db-pass.tab"
    log: CWD + "/results/{directory}/run1/corrected_output/{directory}_db-pass.log"
    params:
        organism = config['organism'],
        outdir = CWD + "/results/{directory}/run1/corrected_output/",
        analysisname = "{directory}"
    shell: "scripts/correct_regions.py --organism {params.organism} --makedb {input} --outdir {params.outdir} --analysisname {params.analysisname}"


rule define_clones_one:
    input: CWD + "/results/{directory}/run1/corrected_output/{directory}_db-pass.tab"
    output: CWD + "/results/{directory}/run1/clones_output/{directory}_db-pass_clone-pass.tab"
    log: CWD + "/results/{directory}/run1/clones_output/{directory}_db-pass-clones.log"
    params:
        outdir = CWD + "/results/{directory}/run1/clones_output"
    shell: "DefineClones.py -d {input} --outdir {params.outdir} --failed  --act set --model ham --norm len --dist 0.15 &> {log}"


rule data_integration_one:
    input:
        igblast_airr_out = CWD+"/results/{directory}/run1/igblast_output/{directory}.fmt19",
        igblast_fmt7_out = CWD+"/results/{directory}/run1/igblast_output/{directory}.fmt7",
        define_clones_out = CWD+"/results/{directory}/run1/clones_output/{directory}_db-pass_clone-pass.tab",
        corrected_out = CWD+"/results/{directory}/run1/corrected_output/{directory}_db-pass.tab"
    output: 
        tsv = CWD+"/results/{directory}/run1/integration_output/{directory}.tsv",
        fasta =  CWD+"/results/{directory}/run1/integration_output/corrected_sequences.fasta"
    log: CWD+"/results/{directory}/run1/integration_output/{directory}.log"
    params:
        organism = config['organism'],
        correct = config['CORRECTION'],
        HV_primer = config['HV_primer'],
        KV_primer = config['KV_primer'],
        LV_primer = config['LV_primer'],
        #imgt_ig_v = "/data04-scratch/toliveira/vramos/changeo_config_files/igblast/fasta/imgt_mouse_ig_v.fasta "
        imgt_ig_v = config['imgt_germline_v_db_fasta']
    shell: " scripts/data_integration.py --igblast_airr {input.igblast_airr_out} --igblast_fmt7 {input.igblast_fmt7_out}  --define_clones {input.define_clones_out} --corrected_regions_file {input.corrected_out} --germline_ig_v_seq {params.imgt_ig_v} --organism {params.organism} --HV_primer {params.HV_primer} --KV_primer {params.KV_primer} --LV_primer {params.LV_primer} --correction {params.correct} --output {output.tsv} &> {log} "


rule igblast_two:
    input: CWD+"/results/{directory}/run1/integration_output/corrected_sequences.fasta",
    output: CWD+"/results/{directory}/run2/igblast_output/{directory}.fmt7"
    log: CWD+"/results/{directory}/run2/igblast_output/{directory}.log"
    params:
        organism = config['organism'],
        germline_db_V = config['igblast_germline_db_V'],
        germline_db_D = config['igblast_germline_db_D'],
        germline_db_J = config['igblast_germline_db_J'],
        auxiliary_data = config['igblast_auxiliary_data']
    shell: " igblastn -germline_db_V {params.germline_db_V}  -germline_db_D {params.germline_db_D} -germline_db_J  {params.germline_db_J} -auxiliary_data {params.auxiliary_data} -domain_system imgt -ig_seqtype Ig -organism {params.organism} -outfmt '7 std qseq sseq btop' -query {input} -out {output} &> {log} "


rule igblast_airr_two:
    input: CWD+"/results/{directory}/run1/integration_output/corrected_sequences.fasta",
    output: CWD+"/results/{directory}/run2/igblast_output/{directory}.fmt19"
    log: CWD+"/results/{directory}/run2/igblast_output/{directory}.log"
    params:
        organism = config['organism'],
        germline_db_V = config['igblast_germline_db_V'],
        germline_db_D = config['igblast_germline_db_D'],
        germline_db_J = config['igblast_germline_db_J'],
        auxiliary_data = config['igblast_auxiliary_data']
    shell: " igblastn -germline_db_V {params.germline_db_V}  -germline_db_D {params.germline_db_D} -germline_db_J {params.germline_db_J} -auxiliary_data {params.auxiliary_data} -domain_system imgt -ig_seqtype Ig -organism {params.organism} -outfmt 19 -query {input} -out {output} &> {log} "


rule makedb_two:
    input:
        igblast_output = CWD+"/results/{directory}/run2/igblast_output/{directory}.fmt7",
        fasta_file = CWD+"/results/{directory}/run1/integration_output/corrected_sequences.fasta" 
    output: CWD+"/results/{directory}/run2/makedb_output/{directory}_db-pass.tab"
    log: CWD+"/results/{directory}/run2/makedb_output/{directory}_db-pass.log"
    params:
        imgt_ig_v = config['imgt_germline_db_V'],
        imgt_ig_d = config['imgt_germline_db_D'],
        imgt_ig_j = config['imgt_germline_db_J'],
        outdir    = CWD+"/results/{directory}/run2/makedb_output"
    shell: " MakeDb.py igblast -i {input.igblast_output} -s {input.fasta_file} -r {params.imgt_ig_v} {params.imgt_ig_d} {params.imgt_ig_j} --outdir {params.outdir} --regions --scores --failed --partial --cdr3 &> {log} "


rule correct_regions_two:
    input: CWD+"/results/{directory}/run2/makedb_output/{directory}_db-pass.tab"
    output: CWD+"/results/{directory}/run2/corrected_output/{directory}_db-pass.tab"
    log: CWD+"/results/{directory}/run2/corrected_output/{directory}_db-pass.log"
    params:
        organism = config['organism'],
        outdir = CWD+"/results/{directory}/run2/corrected_output/",
        analysisname = "{directory}" 
    shell: " scripts/correct_regions.py --organism {params.organism} --makedb {input} --outdir {params.outdir} --analysisname {params.analysisname} "


rule define_clones_two:
    input: CWD+"/results/{directory}/run2/corrected_output/{directory}_db-pass.tab"
    output: CWD+"/results/{directory}/run2/clones_output/{directory}_db-pass_clone-pass.tab"
    log: CWD+"/results/{directory}/run2/clones_output/{directory}_db-pass-clones.log"
    params: 
        outdir = CWD+"/results/{directory}/run2/clones_output"
    shell: " DefineClones.py -d {input} --outdir {params.outdir} --failed --act set --model ham --norm len --dist 0.15 &> {log} "


rule data_integration_two:
    input:
        igblast_airr_out = CWD+"/results/{directory}/run2/igblast_output/{directory}.fmt19",
        igblast_fmt7_out = CWD+"/results/{directory}/run2/igblast_output/{directory}.fmt7",
        define_clones_out = CWD+"/results/{directory}/run2/clones_output/{directory}_db-pass_clone-pass.tab",
        corrected_out = CWD+"/results/{directory}/run2/makedb_output/{directory}_db-pass.tab"
    output: CWD+"/results/{directory}/run2/integration_output/{directory}.tsv",
    log: CWD+"/results/{directory}/run2/integration_output/{directory}.log"
    params:
        organism = config['organism'],
        correct = config['CORRECTION'],
        HV_primer = config['HV_primer'],
        KV_primer = config['KV_primer'],
        LV_primer = config['LV_primer'],
        #imgt_ig_v = "/data04-scratch/toliveira/vramos/changeo_config_files/igblast/fasta/imgt_mouse_ig_v.fasta "
        imgt_ig_v = config['imgt_germline_v_db_fasta']
    shell: " scripts/data_integration.py --igblast_airr {input.igblast_airr_out} --igblast_fmt7 {input.igblast_fmt7_out}  --define_clones {input.define_clones_out} --corrected_regions_file {input.corrected_out} --germline_ig_v_seq {params.imgt_ig_v} --organism {params.organism} --HV_primer {params.HV_primer} --KV_primer {params.KV_primer} --LV_primer {params.LV_primer} --correction {params.correct} --output {output} &> {log} "


rule xlsx_file:
    input: CWD+"/results/{directory}/run2/integration_output/{directory}.tsv" if correct == True else  CWD+"/results/{directory}/run1/integration_output/{directory}.tsv"
    output: CWD+"/results/{directory}/xlsx_output/{directory}_summary.xlsx"
    params:
        analysisname = config['ANALYSISNAME'],
        clones_one_output = CWD+"/results/{directory}/run1/clones_output/{directory}_db-fail.tab",
    log: CWD+"/results/{directory}/xlsx_output/{directory}_summary.log"
    shell: 
        """
        scripts/generate_xlsx.py --changeo_integrated_file {input} --output {output} &> {log}
        """

rule create_ig_pairs:
    input: CWD+"/results/{directory}/xlsx_output/{directory}_summary.xlsx"
    output:
        all = CWD+"/results/{directory}/igpairs/{directory}_strict.xlsx",
        selected = CWD+"/results/{directory}/igpairs/{directory}_strict_selected_columns.xlsx",
        all_combo = CWD + "/results/circos_plots/data/{directory}_strict.xlsx",
        selected_combo = CWD + "/results/circos_plots/data/{directory}_strict_selected_columns.xlsx",
    log: CWD+"/results/{directory}/igpairs/{directory}.log",
    params:
        output_prefix=CWD+"/results/{directory}/igpairs/{directory}",
        pairs_file=CWD+"/data/pairs.xlsx"
    shell:
        """
        scripts/IgPairs.pl parse_excel --input_file {input} \
            --pairs_file {params.pairs_file} \
            --output_prefix {params.output_prefix} \
            &> {log}
            cp {output.all} {output.all_combo}
            cp {output.selected} {output.selected_combo}
        """


rule create_circos_plot:
    input:
        xlsx_files = RESULT_FILE
    output: directory(CWD + "/results/circos_plots/plots/")
    params:
        data = CWD + "/results/circos_plots/data/",
        sample_order = config['SAMPLE_ORDER']
    log: CWD + "/results/circos_plots/plots.log"
    shell:
        """
        Rscript scripts/ParseToCircos2.R -i {params.data} -o {output} --sample_order {params.sample_order}  &> {log}
        """


rule SHM:
    input:
        xlsx_files = RESULT_FILE,
        plots  = CWD + "/results/circos_plots/plots/"
    output: directory(CWD + "/results/SHM_analysis/")
    params:
        data = CWD + "/results/circos_plots/data/"

    log: CWD + "/results/SHM_analysis/SHM_analysis.log"
    shell:
        """
        Rscript scripts/SHM_analysis.R -i {params.data} -o {output} &> {log}
        """


rule hydrophobicity:
    input:
        xlsx_files = RESULT_FILE,
        shm = CWD + "/results/SHM_analysis/"
    output: directory(CWD + "/results/hydrophobicity_analysis/")
    params:
        data = CWD + "/results/circos_plots/data/",
        bcelldbrds = bcell_database
    log: CWD + "/results/hydrophobicity_analysis/hydrophobicity_analysis.log"
    shell:
        """
        Rscript scripts/hydrophobicity_analysis.R -i {params.data} -o {output} --rdsbcelldb {params.bcelldbrds} &> {log}
        """
