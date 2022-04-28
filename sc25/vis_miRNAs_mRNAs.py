#%%
import pickle
from itertools import chain
from pathlib import Path
from typing import Callable, Dict, List, Optional, Set, Tuple

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from sklearn.preprocessing import minmax_scale
from scipy.stats import gmean

low_cutoff, hi_cutoff = 0.7, 1.3
pval_cutoff, padj_cutoff = 0.1, 0.1
expr_cutoff = 0.1

# species = "hsa"
species = "cel"

supported_species_ref_cols = {
    "cel": ("Wormbase Gene ID", "miRNA family"),
    "hsa": ("Gene Symbol", "miRNA family"),
}

supported_species = set(supported_species_ref_cols.keys())

if species not in supported_species:
    raise NotImplementedError(
        f"Only supported species are {', '.join(supported_species)}"
    )

mrna_file: Path = Path(f"counts/{species}/{species}_mRNAs_deseq2.xlsx")
mirna_file: Path = Path(f"counts/{species}/{species}_miRNAs_deseq2.xlsx")
sirna_file: Optional[Path] = Path(f"counts/{species}/{species}_22g_siRNAs_deseq2.xlsx")
sirna_file = sirna_file if (sirna_file and sirna_file.exists()) else None

# mirnas_of_interest = None
mirnas_of_interest: Set[str] = (
    {
        "hsa-let-7a-5p",
        "hsa-let-7b-5p",
        "hsa-let-7c-5p",
        "hsa-let-7d-5p",
        "hsa-let-7e-5p",
        "hsa-let-7f-5p",
        "hsa-let-7g-5p",
        "hsa-let-7i-5p",
        "hsa-miR-98-5p",
    }
    if species == "hsa"
    else set()
)
mirna_group_of_interest = "let-7 family miRNAs"

mirna_targeting_file = Path(f"ref/{species}/mirna_targeting.csv.gz")
transcript_id_mapping_file = Path(f"ref/{species}/transcript_to_gene_id.csv")
mirna_mature_sequences_file = Path("ref/misc/mature.fa")

highlight_mRNAs_files: List[Tuple[Path, str]] = [
    (
        Path(f"counts/{species}/DAF-12_targets_sc25_dauer.xlsx"),
        "detected DAF-12 target transcripts",
    ),
    (
        Path(f"counts/{species}/DAF-16_targets_sc25_dauer.xlsx"),
        "detected DAF-16 target transcripts",
    ),
] if species == "cel" else []

log2fc_col, expr_col = "log2FoldChange", "baseMean"
pval_col, padj_col = "pvalue", "padj"
mrna_col, mirna_col, sirna_col = "geneID", "miRNA", "transcriptID"

mrna_ref_col, mirna_seed_ref_col = supported_species_ref_cols[species]

mirna_targeting_cache, transcript_mapping_cache = (
    Path(f"ref/{species}/mirna_targeting.pkl"),
    Path(f"ref/{species}/transcript_mapping.pkl"),
)


def load_genes_target_mirnas_mapping(
    mirna_targeting_file: Path,
    mirna_mature_sequences_file: Path,
    mirna_targeting_cache: Path,
) -> Dict[str, Set[str]]:
    mirna_targeting = pd.read_csv(
        mirna_targeting_file,
        usecols=(mrna_ref_col, mirna_seed_ref_col, "Representative miRNA"),
    ).drop_duplicates()

    mirna_seq_map = {}
    mirna_mature_sequences: SeqIO.FastaIO.FastaIterator = SeqIO.parse(  # type: ignore
        mirna_mature_sequences_file, format="fasta"
    )
    for record in mirna_mature_sequences:
        if species in record.id:
            mirna_seq_map[str(record.seq)] = record.id
    mirna_seq_map = list(mirna_seq_map.items())

    # drop any families which don't have a representative miRNA for the species being considered
    mirna_targeting: pd.DataFrame = mirna_targeting[
        mirna_targeting["Representative miRNA"].str.startswith(species, na=False)
    ]

    genes_target_seeds_map = mirna_targeting.groupby(mrna_ref_col, sort=False)
    genes_target_mirnas_map: Dict[str, Set[str]] = {}
    total_genes = len(genes_target_seeds_map)
    for i, (gene, gene_info) in enumerate(genes_target_seeds_map, start=1):
        try:
            genes_target_mirnas_map[gene] = set(
                chain.from_iterable(
                    (
                        mirna
                        for seq, mirna in mirna_seq_map
                        if seed
                        in seq[
                            0 : len(seq) // 2
                        ]  # seed should be in the 5' half of the miRNA sequence
                    )
                    for seed in gene_info[mirna_seed_ref_col]
                )
            )
            representative_targeting_mirnas: set[str] = set(
                mirna[
                    : (dotidx if (dotidx := mirna.find(".")) >= 0 else len(mirna))
                ]  # strip weird disambiguation info not present in mirbase files
                for mirna in gene_info["Representative miRNA"]
            )
            if not genes_target_mirnas_map[gene].issuperset(
                representative_targeting_mirnas
            ):
                genes_target_mirnas_map[gene].update(representative_targeting_mirnas)
            if i % 1000 == 0:
                print(f"Processed target miRNAs for {i} of {total_genes} genes")
        except KeyError as key_error:
            missing_key = key_error.args[0]
            source = genes_target_seeds_map.get_group(gene)[
                [mirna_seed_ref_col, "Representative miRNA"]
            ]
            print(source[source[mirna_seed_ref_col] == missing_key])
    pickle.dump(
        genes_target_mirnas_map,
        mirna_targeting_cache.open("wb"),
        pickle.HIGHEST_PROTOCOL,
    )
    return genes_target_mirnas_map


genes_target_mirnas_map = None
if mirna_targeting_cache.exists():
    try:
        genes_target_mirnas_map = pickle.load(mirna_targeting_cache.open("rb"))
        print("Loaded miRNA targeting info from cache")
    except Exception:
        genes_target_mirnas_map = load_genes_target_mirnas_mapping(
            mirna_targeting_file, mirna_mature_sequences_file, mirna_targeting_cache
        )
else:
    genes_target_mirnas_map = load_genes_target_mirnas_mapping(
        mirna_targeting_file, mirna_mature_sequences_file, mirna_targeting_cache
    )

# import json
# json.dump({k: list(vs) for k, vs in genes_target_mirnas_map.items()}, Path(f"{mirna_targeting_cache}.json").open('w'), indent=2)


def load_transcript_id_mapping(
    transcript_id_mapping_file: Path, transcript_mapping_cache: Path
) -> Optional[Dict[str, Set[str]]]:
    transcript_id_mapping = (
        pd.read_csv(
            transcript_id_mapping_file,
            usecols=("initial_alias", "converted_alias", "name"),
        ).drop_duplicates()
        if Path(transcript_id_mapping_file).exists()
        else None
    )
    if transcript_id_mapping is not None:
        transcript_id_mapping.rename(
            {"initial_alias": "transcript_id", "converted_alias": "gene_id",},
            axis="columns",
            inplace=True,
        )
        transcript_id_mapping = transcript_id_mapping.groupby("gene_id", sort=False)
        transcript_id_mapping = {
            gene_id: set(info["transcript_id"].unique())
            for gene_id, info in transcript_id_mapping
        }
    pickle.dump(
        transcript_id_mapping,
        transcript_mapping_cache.open("wb"),
        pickle.HIGHEST_PROTOCOL,
    )
    return transcript_id_mapping


transcript_id_mapping = None
if transcript_mapping_cache.exists():
    try:
        transcript_id_mapping = pickle.load(transcript_mapping_cache.open("rb"))
        print("Loaded transcript gene id mapping from cache")
    except Exception:
        transcript_id_mapping = load_transcript_id_mapping(
            transcript_id_mapping_file, transcript_mapping_cache
        )
else:
    transcript_id_mapping = load_transcript_id_mapping(
        transcript_id_mapping_file, transcript_mapping_cache
    )


def classify_fc(
    mrna_log2fc: float,
    sRNA_log2fc: float,
    sRNA_type: str = "miRNA",
    split_not_significant: bool = True,
) -> str:
    mrna_fc, sRNA_fc = 2 ** mrna_log2fc, 2 ** sRNA_log2fc
    if pd.isna(mrna_fc) and pd.isna(sRNA_fc) and split_not_significant:
        return f"mRNA undetected, {sRNA_type} undetected/unknown"
    elif pd.isna(mrna_fc) and split_not_significant:
        if sRNA_fc <= low_cutoff:
            return f"{sRNA_type} down, mRNA undetected"
        elif sRNA_fc >= hi_cutoff:
            return f"{sRNA_type} up, mRNA undetected"
        else:
            return "neither"
    elif pd.isna(sRNA_fc) and split_not_significant:
        if mrna_fc <= low_cutoff:
            return f"mRNA down, {sRNA_type} undetected/unknown"
        elif mrna_fc >= hi_cutoff:
            return f"mRNA up, {sRNA_type} undetected/unknown"
        else:
            return "neither"
    elif (pd.isna(mrna_fc) or pd.isna(sRNA_fc)) and not split_not_significant:
        return "neither"
    elif mrna_fc <= low_cutoff and sRNA_fc >= hi_cutoff:
        return f"mRNA down, {sRNA_type} up"
    elif sRNA_fc <= low_cutoff and mrna_fc >= hi_cutoff:
        return f"mRNA up, {sRNA_type} down"
    elif mrna_fc <= low_cutoff and split_not_significant:
        return "only mRNA down"
    elif mrna_fc >= hi_cutoff and split_not_significant:
        return "only mRNA up"
    else:
        return "neither"


def savefig(
    fig: matplotlib.figure.Figure, filename: Optional[Path]
) -> matplotlib.figure.Figure:
    if filename:
        fig.savefig(filename, bbox_inches="tight")
    plt.close(fig=fig)
    return fig


def expr_plot(
    df: pd.DataFrame,
    pval_col: str = padj_col,
    fc_col: str = log2fc_col,
    expr_col: str = expr_col,
    label_col: str = "label",
    figsize: Tuple[float, float] = (10, 8),
    filename: Optional[Path] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    convert_log2fc: bool = False,
    convert_log10expr: bool = False,
    mairplot: bool = False,
    colors: Tuple[str, str, str] = ("red", "gray", "blue"),
    pval_cutoff: Tuple[float, float] = (0.01, 0.01),
    fc_cutoff: Tuple[float, float] = (-0.05, 0.05),
    expr_cutoff: float = 0.1,
    clip_limits: Tuple[float, float] = (1e-10, 10),
    lines: bool = False,
    labels: Optional[List[str]] = None,
    label_neg_outliers: int = 0,
    label_pos_outliers: int = 0,
    plot_postprocess: Optional[
        Callable[
            [matplotlib.figure.Figure, matplotlib.axes.Axes],
            Tuple[matplotlib.figure.Figure, matplotlib.axes.Axes],
        ]
    ] = None,
    title_item: str = "",
) -> Tuple[pd.DataFrame, matplotlib.axes.Axes]:
    df = df.copy()

    def mark_change(fc, expr, p):
        return (
            "-ve"
            if fc <= fc_cutoff[0] and p <= pval_cutoff[0] and expr >= expr_cutoff
            else "+ve"
            if fc >= fc_cutoff[1] and p <= pval_cutoff[0] and expr >= expr_cutoff
            else "/"
        )

    df["change"] = df.apply(
        lambda x: mark_change(x[fc_col], x[expr_col], x[pval_col]), axis=1
    )

    df["mlog10pval"] = -np.log10(np.clip(df[pval_col], clip_limits[0], clip_limits[1]))

    if convert_log2fc:
        df["log2fc"] = np.log2(df[fc_col])
        fc_col = "log2fc"
    if convert_log10expr:
        df["log10expr"] = np.log10(df[expr_col])
        expr_col = "log10expr"

    y_col, x_col = "mlog10pval", fc_col
    if mairplot:
        y_col, x_col = fc_col, expr_col

    dist_scale_factor = np.mean(df[y_col]) / np.mean(abs(df[x_col]))
    df["distance"] = df.apply(
        lambda x: np.hypot(dist_scale_factor * abs(x[y_col]), x[x_col]), axis=1,
    )

    fig = plt.figure(dpi=300, figsize=figsize)
    ax = fig.gca()
    ax = sns.scatterplot(
        x=x_col,
        y=y_col,
        hue="change",
        data=df,
        palette={"-ve": colors[0], "/": colors[1], "+ve": colors[2]},
        linewidth=0.2,
        edgecolor="black",
        s=24,
    )

    grouped_by_dir = df.groupby(by="change")
    neg = (
        grouped_by_dir.get_group("-ve").sort_values(by="distance", ascending=False)
        if "-ve" in grouped_by_dir.groups
        else None
    )
    pos = (
        grouped_by_dir.get_group("+ve").sort_values(by="distance", ascending=False)
        if "+ve" in grouped_by_dir.groups
        else None
    )

    current_handles, current_labels = ax.get_legend_handles_labels()
    new_labels = {current_label: "" for current_label in current_labels}

    def get_count(df, col, val=None):
        if val is None:
            return df[col].count()
        return df.loc[df[col] == val, col].count()

    total_cnt = get_count(df, "change")

    if "-ve" in new_labels:
        new_labels[
            "-ve"
        ] = f"Downregulated vs. ctrl. ({get_count(df, 'change', '-ve')})"
    if "/" in new_labels:
        new_labels["/"] = f"Unaltered ({get_count(df, 'change', '/')})"
    if "+ve" in new_labels:
        new_labels["+ve"] = f"Upregulated vs. ctrl. ({get_count(df, 'change', '+ve')})"

    new_labels_list = list(new_labels.values())
    ax.legend(
        current_handles,
        new_labels_list,
        title=f"Of {total_cnt} {title_item}:",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.2),
        ncol=len(new_labels_list),
        handletextpad=0.1,
        columnspacing=0.3,
    )

    if lines:
        if mairplot:
            for yc in fc_cutoff:
                ax.axhline(
                    y=np.log2(max(0, yc)) if convert_log2fc else yc,
                    linestyle="--",
                    color="black",
                    linewidth=0.5,
                )
            ax.axvline(
                x=np.log10(max(0, expr_cutoff)) if convert_log10expr else expr_cutoff,
                linestyle="--",
                color="black",
                linewidth=0.5,
            )
            ax.axhline(y=0, color="black", linewidth=0.5)
        else:
            for xc in fc_cutoff:
                ax.axvline(
                    x=np.log2(max(0, xc)) if convert_log2fc else xc,
                    linestyle="--",
                    color="black",
                    linewidth=0.5,
                )
            for yc in pval_cutoff:
                ax.axhline(
                    y=-np.log10(yc), linestyle="--", color="black", linewidth=0.5
                )
            ax.axvline(x=0, color="black", linewidth=0.5)

    if labels == None:
        labels = []

    if label_neg_outliers > 0 and (neg is not None):
        print(neg.head(label_neg_outliers))
        labels += neg.head(label_neg_outliers)[label_col].tolist()
    if label_pos_outliers > 0 and (pos is not None):
        print(neg.head(label_pos_outliers))
        labels += pos.head(label_pos_outliers)[label_col].tolist()

    all_labels = set(df[label_col].to_list())
    if labels:
        for label in labels:
            if label not in all_labels:
                continue
            text_x = df.loc[df[label_col] == label, x_col].values[0]
            text_y = df.loc[df[label_col] == label, y_col].values[0]
            ax.text(
                x=text_x, y=text_y, s=label, fontsize="xx-small", color="black",
            )

    ax.grid(False)
    ax.yaxis.set_major_locator(mticker.FixedLocator(ax.get_yticks().tolist()))
    ax.set_yticklabels([f"{ts:.1f}" for ts in ax.get_yticks()])

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    if plot_postprocess:
        fig, ax = plot_postprocess(fig, ax)

    fig = savefig(fig, filename)
    return df, ax


def gse_plot(
    gse_file: Path,
    size_func: Callable[[List[str]], float],
    size_col_name: Optional[str] = None,
) -> Tuple[pd.DataFrame, matplotlib.axes.Axes]:
    gse_data = pd.read_csv(gse_file, index_col="term_id")
    gse_data["percent_mapped"] = 100 * (
        gse_data["intersection_size"] / gse_data["query_size"]
    )
    gse_data["point_size"] = [
        size_func(intersection.split(",")) for intersection in gse_data["intersections"]
    ]
    new_term_col = "Enriched GO Terms"
    new_pt_size_col = "point_size" if size_col_name is None else size_col_name
    gse_data.rename(
        columns={"term_name": new_term_col, "point_size": new_pt_size_col}, inplace=True
    )
    gse_data.drop_duplicates(subset=[new_term_col], inplace=True)
    gse_data[new_term_col] = gse_data[new_term_col].apply(
        lambda term: term[0].upper() + term[1:]
    )
    sns.set_theme(style="ticks", font_scale=2.5)
    fig = plt.figure(dpi=300, figsize=(10, 10))
    default_markersize = plt.rcParams["lines.markersize"]
    plt.rcParams["lines.markersize"] = default_markersize * 3
    ax = fig.gca()
    ax = sns.scatterplot(
        data=gse_data,
        x="percent_mapped",
        y="negative_log10_of_adjusted_p_value",
        hue=new_term_col,
        size=new_pt_size_col,
        ax=ax,
        palette=sns.color_palette("husl", len(gse_data)),
        linewidth=1,
        edgecolor="black",
    )
    plt.rcParams["lines.markersize"] = default_markersize
    handles, labels = ax.get_legend_handles_labels()

    last_term_legend_idx = len(gse_data[new_term_col]) + 1  # +1 for title
    title_1_mm = labels[0].replace(" ", "\\ ")
    labels[0] = f"$\\hspace{{{128/len(title_1_mm)}}}\\bf{{\\mathrm{{{title_1_mm}}}}}$"
    title_2_mm = labels[last_term_legend_idx].replace(" ", "\\ ")
    labels[
        last_term_legend_idx
    ] = f"$\\hspace{{{128/len(title_2_mm)}}}\\bf{{\\mathrm{{{title_2_mm}}}}}$"
    # labels = labels[:last_term_legend_idx+1] + [f"{float(val):0.2f}" for val in labels[last_term_legend_idx+1:]]
    ax.legend(
        handles, labels, bbox_to_anchor=(1.025, 1), loc="upper left", labelspacing=0.5
    )

    ax.xaxis.set_major_locator(mticker.FixedLocator(ax.get_xticks().tolist()))
    ax.set_xticklabels([f"{xt}%" for xt in ax.get_xticks()])

    ax.set_xlabel("Target genes mapped to GO term")
    ax.set_ylabel("$-\\log_{10}\\left(p_{adj}\\right)$")
    sns.set_theme(style="ticks", font_scale=1)
    savefig(fig, f"figures/{species}/{Path(gse_file).stem}.svg")
    return gse_data, ax


def make_relplot(
    df: pd.DataFrame,
    x: str,
    y: str,
    category: str,
    xlabel: str,
    ylabel: str,
    colormap: Optional[Dict[str, str]] = None,
) -> matplotlib.figure.Figure:
    data_fill = df.copy().fillna(0)
    fig1 = plt.figure(dpi=300, figsize=(10, 8))
    ax1 = fig1.gca()

    ax1 = sns.scatterplot(
        data=data_fill,
        y=y,
        x=x,
        hue=category,
        s=24,
        linewidth=0.4,
        edgecolor="black",
        ax=ax1,
        palette=colormap,
    )
    handles, labels = ax1.get_legend_handles_labels()
    labels = {
        label: data_fill[data_fill[category] == label].count().min() for label in labels
    }
    total = sum(labels.values())
    labels = [
        f"{label} = {count} ({count/total*100:.02f}%)"
        for label, count in labels.items()
    ]
    if len(labels) > 0:
        handles, labels = zip(*sorted(zip(handles, labels), key=lambda x: x[1]))
        ax1.legend(handles, labels, bbox_to_anchor=(1.025, 0.75))
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.axvline(x=0, linewidth=1, color="black")
    ax1.axvline(x=np.log2(low_cutoff), linewidth=1, color="grey")
    ax1.axvline(x=np.log2(hi_cutoff), linewidth=1, color="grey")
    ax1.axhline(y=0, linewidth=1, color="black")
    ax1.axhline(y=np.log2(low_cutoff), linewidth=1, color="grey")
    ax1.axhline(y=np.log2(hi_cutoff), linewidth=1, color="grey")
    return fig1


def make_plots(
    prefix: str, mrna: pd.DataFrame, mirna: pd.DataFrame, sirna: Optional[pd.DataFrame]
) -> pd.DataFrame:
    mrna_2 = mrna[mrna[expr_col] > 0]

    expr_plot(
        mrna_2,
        label_col=mrna_col,
        xlabel="$\\log_2$ fold change of transcript levels (expt. vs. ctrl.)",
        ylabel="$-\\log_{10}$ Bonferroni-corrected p-value",
        pval_cutoff=(pval_cutoff, pval_cutoff),
        fc_cutoff=(np.log2(low_cutoff), np.log2(hi_cutoff)),
        expr_cutoff=expr_cutoff,
        title_item="detected polyA+ transcripts",
        lines=True,
        filename=Path(f"figures/{species}/mRNA_volcano_{prefix}.svg"),
    )
    expr_plot(
        mrna_2,
        label_col=mrna_col,
        ylabel="$\\log_2$ fold change of transcript levels (expt. vs. ctrl.)",
        xlabel="$-\\log_{10}$ normalized mean expression",
        pval_cutoff=(pval_cutoff, pval_cutoff),
        fc_cutoff=(np.log2(low_cutoff), np.log2(hi_cutoff)),
        expr_cutoff=expr_cutoff,
        mairplot=True,
        convert_log10expr=True,
        title_item="detected polyA+ transcripts",
        lines=True,
        filename=Path(f"figures/{species}/mRNA_expr_{prefix}.svg"),
    )

    for (highlight_mrna_file, title) in highlight_mRNAs_files:
        highlight_mRNAs: Set[str] = set(pd.read_excel(highlight_mrna_file)[mrna_col])
        mrna_df_highlighted = mrna_2[mrna_2[mrna_col].isin(highlight_mRNAs)]
        expr_plot(
            mrna_df_highlighted,
            label_col=mrna_col,
            xlabel="$\\log_2$ fold change of transcript levels (expt. vs. ctrl.)",
            ylabel="$-\\log_{10}$ Bonferroni-corrected p-value",
            pval_cutoff=(pval_cutoff, pval_cutoff),
            fc_cutoff=(np.log2(low_cutoff), np.log2(hi_cutoff)),
            expr_cutoff=expr_cutoff,
            title_item=title,
            lines=True,
            filename=Path(f"figures/{species}/mRNA_volcano_{prefix}_{title}.svg"),
        )
        expr_plot(
            mrna_df_highlighted,
            label_col=mrna_col,
            ylabel="$\\log_2$ fold change of transcript levels (expt. vs. ctrl.)",
            xlabel="$-\\log_{10}$ normalized mean expression",
            pval_cutoff=(pval_cutoff, pval_cutoff),
            fc_cutoff=(np.log2(low_cutoff), np.log2(hi_cutoff)),
            expr_cutoff=expr_cutoff,
            mairplot=True,
            convert_log10expr=True,
            title_item=title,
            lines=True,
            filename=Path(f"figures/{species}/mRNA_expr_{prefix}_{title}.svg"),
        )

    mirna = mirna[mirna[expr_col] > 0]

    expr_plot(
        mirna,
        label_col=mirna_col,
        xlabel="$\\log_2$ fold change of miRNA levels (expt. vs. ctrl.)",
        ylabel="$-\\log_{10}$ Bonferroni-corrected p-value",
        pval_cutoff=(pval_cutoff, pval_cutoff),
        fc_cutoff=(np.log2(low_cutoff), np.log2(hi_cutoff)),
        expr_cutoff=expr_cutoff,
        title_item="detected miRNAs",
        lines=True,
        filename=Path(f"figures/{species}/miRNA_volcano_{prefix}.svg"),
    )
    expr_plot(
        mirna,
        label_col=mirna_col,
        ylabel="$\\log_2$ fold change of miRNA levels (expt. vs. ctrl.)",
        xlabel="$-\\log_{10}$ normalized mean expression",
        pval_cutoff=(pval_cutoff, pval_cutoff),
        fc_cutoff=(np.log2(low_cutoff), np.log2(hi_cutoff)),
        expr_cutoff=expr_cutoff,
        mairplot=True,
        convert_log10expr=True,
        title_item="detected miRNAs",
        lines=True,
        filename=Path(f"figures/{species}/miRNA_expr_{prefix}.svg"),
    )

    mirna = mirna[mirna[padj_col] < padj_cutoff]

    mrna["targeting_miRNAs"] = mrna[mrna_col].apply(
        lambda geneid: genes_target_mirnas_map.get(geneid, set())
    )

    if mirnas_of_interest:
        mrna["is_mRNA_of_interest"] = mrna["targeting_miRNAs"].apply(
            lambda mirnas: len(mirnas_of_interest & mirnas) > 0
        )
        mrnas_of_interest = mrna[mrna["is_mRNA_of_interest"] & (mrna[expr_col] > 0)]
        expr_plot(
            mrnas_of_interest,
            label_col=mrna_col,
            xlabel="$\\log_2$ fold change of transcript levels (expt. vs. ctrl.)",
            ylabel="$-\\log_{10}$ Bonferroni-corrected p-value",
            pval_cutoff=(pval_cutoff, pval_cutoff),
            fc_cutoff=(np.log2(low_cutoff), np.log2(hi_cutoff)),
            expr_cutoff=expr_cutoff,
            title_item=f" detected polyA+ transcripts targeted by {mirna_group_of_interest}",
            lines=True,
            filename=Path(f"figures/{species}/mRNAs_of_interest_volcano_{prefix}.svg"),
        )
        expr_plot(
            mrnas_of_interest,
            label_col=mrna_col,
            ylabel="$\\log_2$ fold change of transcript levels (expt. vs. ctrl.)",
            xlabel="$-\\log_{10}$ normalized mean expression",
            pval_cutoff=(pval_cutoff, pval_cutoff),
            fc_cutoff=(np.log2(low_cutoff), np.log2(hi_cutoff)),
            expr_cutoff=expr_cutoff,
            mairplot=True,
            convert_log10expr=True,
            title_item=f" detected polyA+ transcripts targeted by {mirna_group_of_interest}",
            lines=True,
            filename=Path(f"figures/{species}/mRNAs_of_interest_expr_{prefix}.svg"),
        )
        mrnas_of_interest = mrnas_of_interest.copy()
        mrnas_of_interest["targeting_mirnas_log2fc"] = mrnas_of_interest[
            "targeting_miRNAs"
        ].apply(
            lambda targeting_miRNAs: np.array(
                mirna[
                    mirna[mirna_col].isin(targeting_miRNAs)
                    & mirna[mirna_col].isin(mirnas_of_interest)
                ][log2fc_col]
                if targeting_miRNAs
                else []
            )
        )
        mrnas_of_interest["gmean_mirna_log2fc"] = mrnas_of_interest[
            "targeting_mirnas_log2fc"
        ].apply(
            lambda targeting_mirnas_log2fc: np.log2(gmean(2 ** targeting_mirnas_log2fc))
            if targeting_mirnas_log2fc.size > 0
            else pd.NA
        )
        mrnas_of_interest["mirna_rel_type"] = mrnas_of_interest[
            [log2fc_col, "gmean_mirna_log2fc"]
        ].apply(
            lambda row: classify_fc(row[log2fc_col], row["gmean_mirna_log2fc"]),
            axis="columns",
        )
        mrnas_of_interest["mirna_rel_type_2"] = mrnas_of_interest[
            [log2fc_col, "gmean_mirna_log2fc"]
        ].apply(
            lambda row: classify_fc(
                row[log2fc_col], row["gmean_mirna_log2fc"], split_not_significant=False
            ),
            axis="columns",
        )
        mrna_mirna_of_interest_known = mrnas_of_interest[
            (mrnas_of_interest[expr_col] > 0)
            & (mrnas_of_interest[padj_col] < padj_cutoff)
        ].dropna(axis="index", subset=["gmean_mirna_log2fc"])

        mirna_fig = make_relplot(
            mrna_mirna_of_interest_known,
            y=log2fc_col,
            x="gmean_mirna_log2fc",
            category="mirna_rel_type",
            xlabel="$\\log_2$ geometric mean of fold change of targeting miRNAs of interest (expt. vs. ctrl.)",
            ylabel="$\\log_2$  fold change of transcript levels of interest (expt. vs. ctrl.)",
            colormap={
                "mRNA down, miRNA up": "blue",
                "mRNA up, miRNA down": "red",
                "only mRNA down": "yellow",
                "only mRNA up": "green",
            },
        )
        savefig(
            mirna_fig, Path(f"figures/{species}/miRNA_relplot_of_interest_{prefix}.svg")
        )

        mirna_fig_2 = make_relplot(
            mrna_mirna_of_interest_known,
            y=log2fc_col,
            x="gmean_mirna_log2fc",
            category="mirna_rel_type_2",
            xlabel="$\\log_2$ geometric mean of fold change of targeting miRNAs of interest (expt. vs. ctrl.)",
            ylabel="$\\log_2$  fold change of transcript levels of interest (expt. vs. ctrl.)",
            colormap={
                "mRNA down, miRNA up": "blue",
                "mRNA up, miRNA down": "red",
                "neither": "grey",
            },
        )
        savefig(
            mirna_fig_2,
            Path(f"figures/{species}/miRNA_relplot_simple_of_interest_{prefix}.svg"),
        )

    mrna["targeting_mirnas_log2fc"] = mrna["targeting_miRNAs"].apply(
        lambda targeting_miRNAs: np.array(
            mirna[mirna[mirna_col].isin(targeting_miRNAs)][log2fc_col]
            if targeting_miRNAs
            else []
        )
    )
    mrna["gmean_mirna_log2fc"] = mrna["targeting_mirnas_log2fc"].apply(
        lambda targeting_mirnas_log2fc: np.log2(gmean(2 ** targeting_mirnas_log2fc))
        if targeting_mirnas_log2fc.size > 0
        else pd.NA
    )
    mrna["mirna_rel_type"] = mrna[[log2fc_col, "gmean_mirna_log2fc"]].apply(
        lambda row: classify_fc(row[log2fc_col], row["gmean_mirna_log2fc"]),
        axis="columns",
    )
    mrna["mirna_rel_type_2"] = mrna[[log2fc_col, "gmean_mirna_log2fc"]].apply(
        lambda row: classify_fc(
            row[log2fc_col], row["gmean_mirna_log2fc"], split_not_significant=False
        ),
        axis="columns",
    )
    mrna_mirna_known = mrna[
        (mrna[expr_col] > 0) & (mrna[padj_col] < padj_cutoff)
    ].dropna(axis="index", subset=["gmean_mirna_log2fc"])

    mirna_fig = make_relplot(
        mrna_mirna_known,
        y=log2fc_col,
        x="gmean_mirna_log2fc",
        category="mirna_rel_type",
        xlabel="$\\log_2$ geometric mean of fold change of targeting miRNAs (expt. vs. ctrl.)",
        ylabel="$\\log_2$  fold change of transcript levels (expt. vs. ctrl.)",
        colormap={
            "mRNA down, miRNA up": "blue",
            "mRNA up, miRNA down": "red",
            "only mRNA down": "yellow",
            "only mRNA up": "green",
        },
    )
    savefig(mirna_fig, Path(f"figures/{species}/miRNA_relplot_{prefix}.svg"))

    mirna_fig_2 = make_relplot(
        mrna_mirna_known,
        y=log2fc_col,
        x="gmean_mirna_log2fc",
        category="mirna_rel_type_2",
        xlabel="$\\log_2$ geometric mean of fold change of targeting miRNAs (expt. vs. ctrl.)",
        ylabel="$\\log_2$  fold change of transcript levels (expt. vs. ctrl.)",
        colormap={
            "mRNA down, miRNA up": "blue",
            "mRNA up, miRNA down": "red",
            "neither": "grey",
        },
    )
    savefig(mirna_fig_2, Path(f"figures/{species}/miRNA_relplot_simple_{prefix}.svg"))

    if sirna is not None:
        sirna = sirna[sirna[expr_col] > 0]

        expr_plot(
            sirna,
            label_col=sirna_col,
            xlabel="$\\log_2$ fold change of 22G siRNA vs. transcript levels (expt. vs. ctrl.)",
            ylabel="$-\\log_{10}$ Bonferroni-corrected p-value",
            pval_cutoff=(pval_cutoff, pval_cutoff),
            fc_cutoff=(np.log2(low_cutoff), np.log2(hi_cutoff)),
            expr_cutoff=0.01,
            title_item="detected transcripts with 22G siRNAs",
            lines=True,
            filename=Path(f"figures/{species}/siRNA_volcano_{prefix}.svg"),
        )
        expr_plot(
            sirna,
            label_col=sirna_col,
            ylabel="$\\log_2$ fold change of 22G siRNA vs. transcript levels (expt. vs. ctrl.)",
            xlabel="$-\\log_{10}$ normalized mean expression",
            pval_cutoff=(pval_cutoff, pval_cutoff),
            fc_cutoff=(np.log2(low_cutoff), np.log2(hi_cutoff)),
            expr_cutoff=0.01,
            mairplot=True,
            convert_log10expr=True,
            title_item="detected transcripts with 22G siRNAs",
            lines=True,
            filename=Path(f"figures/{species}/siRNA_expr_{prefix}.svg"),
        )

        sirna = sirna[sirna[padj_col] < 0.2]

        mrna["transcript_ids"] = mrna[mrna_col].apply(
            lambda geneid: transcript_id_mapping.get(geneid, [])
        )
        mrna["targeting_sirnas_log2fc"] = mrna["transcript_ids"].apply(
            lambda transcript_ids: np.array(
                sirna[sirna[sirna_col].isin(transcript_ids)][log2fc_col]
                if transcript_ids
                else []
            )
        )
        mrna["gmean_sirna_log2fc"] = mrna["targeting_sirnas_log2fc"].apply(
            lambda targeting_sirnas_log2fc: np.log2(gmean(2 ** targeting_sirnas_log2fc))
            if targeting_sirnas_log2fc.size > 0
            else pd.NA
        )
        mrna["sirna_rel_type"] = mrna[[log2fc_col, "gmean_sirna_log2fc"]].apply(
            lambda row: classify_fc(
                row[log2fc_col],
                row["gmean_sirna_log2fc"],
                sRNA_type="siRNA",
                split_not_significant=False,
            ),
            axis="columns",
        )

        mrna_sirna_known = mrna[
            (mrna[expr_col] > 0) & (mrna[padj_col] < padj_cutoff)
        ].dropna(axis="index", subset=["gmean_sirna_log2fc"])
        sirna_fig = make_relplot(
            mrna_sirna_known,
            y=log2fc_col,
            x="gmean_sirna_log2fc",
            category="sirna_rel_type",
            xlabel="$\\log_2$ geometric mean of fold change of targeting siRNAs (expt. vs. ctrl.)",
            ylabel="$\\log_2$ FC of gene (expt. vs. ctrl.)",
            colormap={
                "mRNA down, siRNA up": "blue",
                "mRNA up, siRNA down": "red",
                "neither": "grey",
            },
        )
        savefig(sirna_fig, Path(f"figures/{species}/siRNA_relplot_{prefix}.svg"))
    return mrna


mrna: Dict[str, pd.DataFrame] = pd.read_excel(str(mrna_file), sheet_name=None)
mirna: Dict[str, pd.DataFrame] = pd.read_excel(str(mirna_file), sheet_name=None)
sirna: Optional[Dict[str, pd.DataFrame]] = (
    pd.read_excel(str(sirna_file), sheet_name=None)
    if sirna_file and sirna_file.exists()
    else None
)
mod_mrna: Dict[str, pd.DataFrame] = {}
for prefix, df in mrna.items():
    mod_mrna[prefix] = (
        make_plots(prefix, df, mirna[prefix], sirna[prefix])
        if sirna
        else make_plots(prefix, df, mirna[prefix], None)
    )
    mod_mrna[prefix].to_excel(
        Path(f"{Path(mrna_file).stem}_with_miRNA_{prefix}.xlsx"), index=False
    )


def gse_marker_size(expr_df: pd.DataFrame) -> Callable[[List[str],], float]:
    def impl(geneids: List[str]) -> float:
        expr_df["normalized_geneid"] = expr_df[mrna_col].apply(lambda x: x.upper())
        relevant_mrnas = expr_df[expr_df["normalized_geneid"].isin(geneids)]
        return round(np.log2(gmean(2 ** relevant_mrnas[log2fc_col])), ndigits=2)

    return impl


if species == "cel":
    cel_mrna_df = list(mrna.values())[0]
    gse_plot(
        Path("gene_ontology/cel/cel_mrna_down_mirna_up_go_filtered.csv"),
        gse_marker_size(cel_mrna_df),
        "log2 geometric mean of FC of involved genes",
    )
    gse_plot(
        Path("gene_ontology/cel/cel_mrna_up_mirna_down_go_filtered.csv"),
        gse_marker_size(cel_mrna_df),
        "log2 geometric mean of FC of involved genes",
    )
    gse_plot(
        Path("gene_ontology/cel/daf12_targets_go_filtered.csv"),
        gse_marker_size(cel_mrna_df),
        "log2 geometric mean of FC of involved genes",
    )
    gse_plot(
        Path("gene_ontology/cel/daf16_targets_go_filtered.csv"),
        gse_marker_size(cel_mrna_df),
        "log2 geometric mean of FC of involved genes",
    )
elif species == "hsa":
    prefixes_to_use = [
        "A549",
    ]
    for prefix_to_use in prefixes_to_use:
        gse_plot(
            Path(
                f"gene_ontology/{species}/{prefix_to_use}_upregulated_genes_targeted_by_let7_filtered.csv"
            ),
            gse_marker_size(mod_mrna[prefix_to_use]),
            "log2 geometric mean of FC of involved genes",
        )
        heatmap_genes: pd.DataFrame = pd.read_excel(
            str(Path(f"raw_data/{species}/target_heatmap_{prefix_to_use}.xlsx"))
        )
        raw_expr_data: pd.DataFrame = pd.read_excel(
            str(Path(f"counts/{species}/{species}_mRNAs.xlsx")),
            sheet_name=prefix_to_use,
        )
        raw_expr_data.set_index(mrna_col, inplace=True)
        heatmap_genes.set_index(mrna_col, inplace=True)
        ctrl_cols = [col for col in raw_expr_data.columns if col.startswith("CTRL")]
        expt_cols = [col for col in raw_expr_data.columns if col.startswith("EXPT")]
        expr_cols = ctrl_cols + expt_cols
        # begin converting fractional counts to TPM
        for expr_col in expr_cols:
            print(f"Converting counts to TPM for column {expr_col}")
            raw_expr_data[expr_col] = raw_expr_data[expr_col] / (
                raw_expr_data["Length"] / 1e3  # convert length from bases to kilobases
            )
            rpk_scale_factor = raw_expr_data[expr_col].sum() / 1e6
            raw_expr_data[expr_col] = raw_expr_data[expr_col] / rpk_scale_factor
            heatmap_genes[expr_col] = raw_expr_data[expr_col]
        # end converting fractional counts to TPM
        heatmap_genes[expr_cols] = minmax_scale(heatmap_genes[expr_cols], axis=1)
        prefix_mod_mrna_df = mod_mrna[prefix_to_use].set_index(mrna_col)
        heatmap_genes["is_direct_target"] = heatmap_genes.index.map(
            lambda geneid: len(
                prefix_mod_mrna_df.loc[geneid, "targeting_miRNAs"] & mirnas_of_interest
            )
            > 0
        )
        if "Type" not in heatmap_genes.columns:
            heatmap_genes["Type"] = heatmap_genes["is_direct_target"].apply(
                lambda is_direct_target: "direct" if is_direct_target else "indirect"
            )

        heatmap_genes_by_type = heatmap_genes.groupby("Type")
        for type, heatmap_df in heatmap_genes_by_type:
            fig = plt.figure(figsize=(10, 10), dpi=300)
            ax = fig.gca()
            ax = sns.heatmap(
                heatmap_df[expr_cols], square=True, robust=True, cmap="viridis", ax=ax
            )
            savefig(
                fig,
                Path(f"figures/{species}/target_heatmap_{prefix_to_use}_{type}.svg"),
            )

        heatmap_genes.to_excel(
            str(
                Path(
                    f"raw_data/{species}/target_heatmap_{prefix_to_use}_with_expr.xlsx"
                )
            )
        )

# %%
