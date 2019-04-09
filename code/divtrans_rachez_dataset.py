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

# # Divergent transcription validation
#
# We want to validate the divergent identification dataset by comparing RNA-seq and K27Ac Chip-seq data from Rachez's dataset.

# +
from pynextgen import divtrans
from pynextgen.bed import Bed
from pynextgen.basics_bam import Bam
import glob
import pandas as pd
from matplotlib import rcParams
import seaborn as sns

rcParams["figure.figsize"] = (16.0, 8.0)
# -

# %load_ext autoreload
# %autoreload 2

# ## Input
# Get mm9 genome size for fisher test:

# + {"language": "bash"}
#
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from mm9.chromInfo" \
#     | sort -k1,1 -k2,2n > ../data/rachez_test_data/mm9.genome
# -

chro_sizes = "../data/rachez_test_data/mm9.genome"

# Bam files:

bams = glob.glob("../data/rachez_test_data/*.bam")
bams

# Subset test bam:

subset_bam = Bam(bams[0]).multi_fetch(["chr1:0-249250621"], bam_out="test.bam")


subset_bam

# ## Generate bigwigs

# + {"language": "bash"}
#
# bamCoverage --bam test.bam  --filterRNAstrand forward -o test_forward.bw --region chr1:0:249250621
# bamCoverage --bam test.bam  --filterRNAstrand reverse -o test_reverse.bw --region chr1:0:249250621
# -

# ## Divergent transcription detection (test dataset)

test_analysis = divtrans.DivTransFromBam(
    subset_bam.path,
    distance=100,
    max_template_length=500,
    no_overlap=True,
    remove_GL_contigs=True,
)

test_analysis.run()

# ### Filtering

bedfilter = divtrans.BedFilter(test_analysis.bed, subset_bam.path, threads=4)

bedfilter.multiple_run(
    chro_sizes,
    flanks=[100, 250, 500],
    count_thress=[2, 5, 10],
    count_ratio_thress=[1, 2, 4],
)

# ## Compare with K27Ac peaks
#
# Using bed file with quality filtered peaks (select)
#
# Rep2 was empty after quality filtering (apparently all peaks called before were noise)

ref_bed = Bed(
    "../data/rachez_test_data/H3K27ac_vs_INPUT-Mouse_Rep1_select_sort.narrowPeak"
).filter_by_chromosome("chr1")

# ### Statistics for parameter tuning:
#
# #### True positive rate:
#
# $TPR = TP/P$
#
# With:
# - TP: number of true positives
# - P: number of real positive cases in the data
#
# #### False positive rate:
#
# $FPR = FP/N$
#
# With:
# - FP: number of false positives
# - N: number of real negative cases in the data
#
# Source:
# https://en.wikipedia.org/wiki/Receiver_operating_characteristic

test_analysis.bed

divtrans_beds = list(bedfilter.filtered_beds.values()).append(test_analysis.bed)

# +
jaccard_df = pd.concat(
    [bed.jaccard(ref_bed) for bed in bedfilter.filtered_beds.values()]
)
fisher_df = pd.concat(
    [bed.fisher(ref_bed, chro_sizes) for bed in bedfilter.filtered_beds.values()]
)

# To Compute the true positive ratio
len_df = pd.DataFrame(
    {bed.name: len(bed) for bed in bedfilter.filtered_beds.values()},
    index=["bed_length"],
).T

# To compute false positives
fp_df = pd.DataFrame(
    {
        bed.name: len(bed.intersect(ref_bed, supp_args="-v"))
        for bed in bedfilter.filtered_beds.values()
    },
    index=["false_positives"],
).T

# Total negative number is considered to be the number of intervals when complementing the K27Ac peaks
total_negatives = len(ref_bed.complement(chro_sizes))

# Join results
sum_df = pd.concat(
    [fisher_df.set_index("bed1"), jaccard_df.set_index("bed1"), len_df, fp_df],
    axis=1,
    sort=False,
).sort_values("jaccard", ascending=False)

# n_intersection from jaccard test can be assimilated to True Positives:
# ie (number of divergent transcription interval overlapping a K27Ac peak)
sum_df.rename(columns={"n_intersections": "true_positives"}, inplace=True)

# Compute true/false positive rate
sum_df["true_positive_rate"] = sum_df["true_positives"] / len(ref_bed)
sum_df["false_positive_rate"] = sum_df["false_positives"] / total_negatives

# Cannot compute the false positive rate.

# Table cleanup
sum_df.index = sum_df.index.str.replace("test_divtrans_", "")
sum_df.index = sum_df.index.str.replace("_counts_filtered", "")
# -

sum_df.sort_values('true_positive_rate', ascending=False).head()

sns.pairplot(
    x_vars="false_positive_rate",
    y_vars="true_positive_rate",
    data=sum_df.reset_index(),
    hue="bed1",
)

# Plot
sum_df["true_positive_rate"].sort_values().plot(kind="bar", color="salmon")

# ## Divergent transcription detection

# +

analysis = [
    divtrans.DivTransFromBam(
        bam,
        distance=100,
        max_template_length=500,
        no_overlap=True,
        remove_GL_contigs=True,
    )
    for bam in bams
]
# -

[an.run() for an in analysis]

# ## Divergent transcription bed filtering

divtrans.BedFilter(analysis[0].bed, bams[0], threads=4).run(
    chro_sizes, flank=500, count_thres=5, count_ratio_thres=2
)


# + {"toc-hr-collapsed": true, "cell_type": "markdown"}
# ## Compare with K27Ac Chip peaks
# -

rep1_k27ac = Bed(
    "../data/rachez_test_data/H3K27ac_Rep1_vs_INPUT-Mouse_Rep1_peaks.narrowPeak"
)
rep2_k27ac = Bed(
    "../data/rachez_test_data/H3K27ac_Rep1_vs_INPUT-Mouse_Rep1_peaks.narrowPeak"
)
k27ac_vs_input = Bed(
    "../data/rachez_test_data/H3K27ac_vs_INPUT-Mouse_Rep1_select.narrowPeak"
).sort()

analysis[0].bed.sort().fisher(rep2_k27ac, "../data/rachez_test_data/mm9.genome")
