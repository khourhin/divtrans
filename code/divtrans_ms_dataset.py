# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.3'
#       jupytext_version: 1.0.5
#   kernelspec:
#     display_name: Python [conda env:classic]
#     language: python
#     name: conda-env-classic-py
# ---

import pandas as pd
from pynextgen.divtrans import DivTransFromBam

# # Application to MS dataset
# ## Input

# +
ms_meta = pd.read_csv("../../ms_enh/data/metadata.csv")

# Get Healthy bams
healthy_bams = [
    "../../ms_enh/results/snakemake/star/bams/" + x + ".bam"
    for x in ms_meta.loc[ms_meta.Main.isin(["ia", "Healthy_control"])].loc[
        :, "Petient ID"
    ]
]

# Get ms bams
rrms_bams = [
    "../../ms_enh/results/snakemake/star/bams/" + x + ".bam"
    for x in ms_meta.loc[ms_meta.Diagnose == "RRMS"].loc[:, "Petient ID"]
]
# -

# Detect divergent transcription in MS dataset
healthy_divtrans = [
    DivTransFromBam(
        bam,
        distance=100,
        max_template_length=500,
        no_overlap=True,
        remove_GL_contigs=True,
    )
    for bam in healthy_bams
]

# Run the analysis
[divtrans.run() for divtrans in healthy_divtrans]

healthy_div_trans_beds = [get_divergent_transcription_beds(bam) for bam in healthy_bams]
healthy_div_trans_merged = [
    bed.merge(outfolder="../results/from_bam") for bed in healthy_div_trans_beds
]

ms_div_trans_beds = [get_divergent_transcription_beds(bam) for bam in ms_bams]
ms_div_trans_merged = [
    bed.merge(outfolder="../results/from_bam") for bed in ms_div_trans_beds
]

# ## Select intervals overlapping with CAGE dataset divergent transcription events

div_trans_beds_merged_filter = [
    bed.intersect(cage_minus_bed, supp_args="-wa").intersect(
        cage_plus_bed, supp_args="-wa"
    )
    for bed in div_trans_beds_merged
]
