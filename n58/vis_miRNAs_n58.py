import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import matplotlib.figure
import matplotlib.pyplot as plt
import seaborn as sns

import os
import sys

# sys.argv[1] = "ref/mirna_pirnas_seq.tsv"
# sys.argv[2] = "counts/N58_miRNA_counts.xlsx"

miRNA_ref = pd.read_table(sys.argv[1], usecols=["name", "seq"])
miRNA_ref.set_index("name", inplace=True)
miRNA_ref = miRNA_ref[~miRNA_ref.index.duplicated(keep="first")]


df = pd.read_excel(sys.argv[2], sheet_name="N2 N58", index_col="miRNA")
df = df[~df.index.duplicated(keep="first")]
df["miRNA"] = df.index
df["seq"] = df["miRNA"].map(lambda x: miRNA_ref["seq"].get(x, ""))

fc_col = "FC_12"
log2fc_col = "log2_fc"
rel_log10expr_col = "log10_rel_expr_ctrl(EV)"
log2fc_text = "log<sub>2</sub> Fold Change (<i>nol-58</i> RNAi vs EV)"
rel_log10expr_text = "log<sub>10</sub> Normalized Expression in EV (Ctrl.)"
ctrl_reads_col, expt_reads_col = "N2_N58_EV", "N2_N58_RNAi"
ctrl_rpm_col, expt_rpm_col = "N_N2_N58_EV", "N_N2_N58_RNAi"
normalization_spike = "hsa-miR-122-5p"
change_col = "Change"

df[ctrl_rpm_col] = df[ctrl_reads_col] * 1e6 / df[ctrl_reads_col].sum()
df[expt_rpm_col] = df[expt_reads_col] * 1e6 / df[expt_reads_col].sum()

ctrl_norm_count = df[ctrl_rpm_col].loc[normalization_spike]
expt_norm_count = df[expt_rpm_col].loc[normalization_spike]

df[ctrl_rpm_col] = df[ctrl_rpm_col] * 1e4 / ctrl_norm_count
df[expt_rpm_col] = df[expt_rpm_col] * 1e4 / expt_norm_count
# df.drop(normalization_spike, inplace=True)

fc_thresholds = [1.1, 1.5, 2.0, 3.0, 4.0]
baseline = fc_thresholds[0]
comparisons = fc_thresholds[1:]

detection_threshold = 0.1  # RPM in Ctrl sample
df = df[(df[ctrl_rpm_col] >= detection_threshold)]
df[fc_col] = df[expt_rpm_col] / df[ctrl_rpm_col]
df[log2fc_col] = np.log2(df[fc_col])
print(len(df))
detected_miRNAs = set(df.index)

max_seq_len = max(df["seq"].map(len))


def pos_base_abundances(seq_df):
    abundances = pd.DataFrame(
        data={
            "pos": list(range(1, max_seq_len + 1)),
            "A": [0] * max_seq_len,
            "T": [0] * max_seq_len,
            "G": [0] * max_seq_len,
            "C": [0] * max_seq_len,
        }
    )
    abundances.set_index("pos", inplace=True)
    for entry in seq_df["seq"]:
        for i, c in enumerate(entry):
            abundances[c.upper()].iloc[i] += 1

    abundances = abundances.dropna().apply(
        lambda row: row * 100 / row.sum(), axis="columns"
    )
    return abundances


baseline_abundances = pos_base_abundances(df[df[fc_col] <= baseline])
fig: matplotlib.figure.Figure = plt.figure(figsize=(8, 7 * len(comparisons)), dpi=300)

for i, comparison in enumerate(comparisons):
    ax = fig.add_subplot(len(comparisons), 1, i + 1)
    comparison_abundances = pos_base_abundances(df[df[fc_col] >= comparison])
    ratio = comparison_abundances / baseline_abundances
    ratio2 = pd.melt(
        ratio.dropna(how="all").reset_index(),
        id_vars="pos",
        var_name="base",
        value_name="abundance",
    )
    sns.barplot(data=ratio2, x="pos", y="abundance", hue="base", ax=ax)
    ax.set_title(
        rf"$\mathrm{{FC }}\geq {comparison}\:\mathrm{{ vs. FC }}\leq {baseline}$"
    )
    ax.set_xticklabels(
        [
            f"{int(pos.get_text())}\n{ratio.loc[int(pos.get_text())].idxmax()}\n{baseline_abundances.loc[int(pos.get_text())].idxmax()}"
            for pos in ax.get_xticklabels()
        ]
    )
    ax.set_xlabel(
        "Position (top), Most Enriched Base in Comparison (middle) vs. Baseline (bottom)"
    )
    ax.set_ylabel("Relative Abundance")
    ax.axhline(y=1, linewidth=1, linestyle="--", color="black")
    ax.set_ylim(0, 3)

fig.subplots_adjust(hspace=0.3)
fig.savefig("extra/base_abundance_comparison.svg")

up_cutoff, down_cutoff = 1.5, 0.7

changes_colors = {
    "Unchanged": "grey",
    "Up": "blue",
    "Down": "red",
    "Undetected": "green",
}
changes = list(changes_colors.keys())


df[change_col] = df[["miRNA", fc_col]].apply(
    lambda x: (
        changes[1]
        if x[fc_col] >= up_cutoff
        else (
            changes[2]
            if x[fc_col] <= down_cutoff
            else (changes[-1] if pd.isna(x[fc_col]) else changes[0])
        )
    )
    if x["miRNA"] in detected_miRNAs
    else changes[-1],
    axis="columns",
)
change_counts_dict = {
    change: (len(df[df[change_col] == change]), len(detected_miRNAs), color)
    for change, color in changes_colors.items()
}
change_counts_dict = {
    change: (f"{change} - {val[0]} / {val[1]} ({val[0]/val[1]*100:.1f}%)", val[-1])
    for change, val in change_counts_dict.items()
}

changes_colors_dict = {val[0]: val[1] for val in change_counts_dict.values()}

df[change_col] = df[change_col].map(lambda x: change_counts_dict[x][0])
df[rel_log10expr_col] = np.log10(df[ctrl_rpm_col])
ifig = px.scatter(
    df,
    x=rel_log10expr_col,
    labels={rel_log10expr_col: rel_log10expr_text, log2fc_col: log2fc_text},
    y=log2fc_col,
    color=change_col,
    color_discrete_map=changes_colors_dict,
    hover_name="miRNA",
    hover_data={col: False for col in df.columns},
    template="simple_white",
)

ifig.add_hline(
    y=0, line_width=0.5, line_color="black", opacity=1,
)
ifig.add_hline(
    y=np.log2(up_cutoff),
    line_width=0.5,
    line_color="grey",
    opacity=1,
    line_dash="dash",
)
ifig.add_hline(
    y=np.log2(down_cutoff),
    line_width=0.5,
    line_color="grey",
    opacity=1,
    line_dash="dash",
)
ifig.add_vline(
    x=np.log10(1),
    line_width=0.5,
    line_color="grey",
    opacity=1,
    line_dash="dash",
)

fw = go.FigureWidget(ifig)
content = fw.data[0]
content.on_click(
    lambda trace, points, selector: [
        fw.add_annotation(point, text=str(point)) for point in points.point_inds
    ]
)

with open("onclick_annotate.js", "r") as extra_js:
    fw.write_html(
        f"{os.path.basename(__file__)}.html",
        include_plotlyjs="cdn",
        post_script=extra_js.read(),
    )
    fw.write_image(f"{os.path.basename(__file__)}.svg", engine="orca")
