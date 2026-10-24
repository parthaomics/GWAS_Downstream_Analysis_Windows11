#-------------------------------------------------------------------------------
# Name:        plink_funcs.py
# Purpose:     helper functions for PLINK utility GUI
# Author:      Michal Martin & Partha Talukdar
# Created:     09/06/2026
# Version:     1.1.0
#-------------------------------------------------------------------------------

import csv
import math
import os
import re
import webbrowser
from collections import defaultdict

import matplotlib.pyplot as plt
import pandas as pd
from openpyxl import Workbook

try:
    from scipy.stats import f as scipy_f_distribution
except Exception:
    scipy_f_distribution = None

HAPLOTYPE_MAX_DISTANCE_BP = 1000
HAPLOTYPE_LD_THRESHOLD = 0.6
VALID_HOMOZYGOTES = {"AA", "TT", "GG", "CC"}

R2_THRESHOLD = 0.6
SIGNIFICANCE_THRESHOLD = 0.05
MIN_GROUP_SIZE = 3
HOMOZYGOUS_CALL_RATE_THRESHOLD = 0.90


def normalize_homozygous_genotype(allele_1, allele_2):
    """Return a homozygous diploid code; treat missing/heterozygous calls as unavailable."""
    allele_1 = str(allele_1).strip().upper()
    allele_2 = str(allele_2).strip().upper()
    if allele_1 in {"", "0"} or allele_2 in {"", "0"}:
        return None
    if allele_1 != allele_2:
        return None
    return allele_1 + allele_2


def add_multiple_testing_results(results_df, p_value_column):
    """Sort valid tests by P value and calculate Bonferroni and BH results."""
    if results_df.empty:
        return results_df.copy()

    adjusted = results_df.copy()
    adjusted[p_value_column] = pd.to_numeric(adjusted[p_value_column], errors="coerce")
    adjusted = adjusted.dropna(subset=[p_value_column]).sort_values(
        p_value_column, kind="stable"
    ).reset_index(drop=True)

    number_of_tests = len(adjusted)
    adjusted["Rank"] = range(1, number_of_tests + 1)
    adjusted["Bonferroni_threshold"] = SIGNIFICANCE_THRESHOLD / number_of_tests
    adjusted["Bonferroni_significant"] = adjusted[p_value_column].apply(
        lambda value: "YES" if value < SIGNIFICANCE_THRESHOLD / number_of_tests else "NO"
    )
    adjusted["BH_q_raw"] = (
        adjusted[p_value_column] * number_of_tests / adjusted["Rank"]
    ).clip(upper=1.0)
    adjusted["BH_q_adjusted"] = adjusted["BH_q_raw"][::-1].cummin()[::-1].clip(upper=1.0)
    adjusted["BH_significant"] = adjusted["BH_q_adjusted"].apply(
        lambda value: "YES" if value < SIGNIFICANCE_THRESHOLD else "NO"
    )
    return adjusted


def validate_inputs(form):
    proceed_flag = True

    bp_lower = form.w_bp_lower.text()
    try:
        ibp_lower = int(bp_lower)
    except ValueError:
        print("Lower base-pair position " + bp_lower + " must be an integer")
        proceed_flag = False

    bp_upper = form.w_bp_upper.text()
    try:
        ibp_upper = int(bp_upper)
    except ValueError:
        print("Upper base-pair position " + bp_upper + " must be an integer")
        proceed_flag = False

    if proceed_flag and ibp_lower >= ibp_upper:
        print("Upper base-pair position " + bp_upper + " must exceed lower base-pair position " + bp_lower)
        proceed_flag = False

    out_dir = form.w_lbl_outdir.text()
    if not os.path.isdir(out_dir):
        print("Output directory " + out_dir + " does not exist")
        proceed_flag = False

    if " " in out_dir:
        print("Blanks not permitted in output directory " + out_dir)
        proceed_flag = False

    root_name = form.w_root_name.text()
    if " " in root_name:
        print("Blanks not permitted in root filename " + root_name)
        proceed_flag = False

    if proceed_flag:
        print("Input fields are valid")
    else:
        print("Errors found - cannot proceed")

    return proceed_flag


def combine_map_ped(form):
    """
    Combine PED, MAP, and phenotype (Fe_Phene.csv) into one table.
    Samples are rows, SNPs are columns, and phenotype is merged by ID.
    """
    out_dir = form.w_lbl_outdir.text()
    root_name = form.w_root_name.text() or "plink"

    ped_file = os.path.join(out_dir, root_name + ".ped")
    map_file = os.path.join(out_dir, root_name + ".map")
    pheno_file = os.path.join(out_dir, "Fe_Phene.csv")
    if not os.path.isfile(pheno_file):
        pheno_file = os.path.join(os.path.dirname(out_dir), "Fe_Phene.csv")

    output_file = os.path.join(out_dir, "map_ped_phene_combined_sample_table.csv")

    if not os.path.isfile(ped_file):
        print("PED file not found: " + ped_file)
        print("Set Root filename to the PED/MAP prefix you want to combine.")
        return

    if not os.path.isfile(map_file):
        print("MAP file not found: " + map_file)
        print("Set Root filename to the PED/MAP prefix you want to combine.")
        return

    if not os.path.isfile(pheno_file):
        print("Phenotype file not found: " + pheno_file)
        print("Expected Fe_Phene.csv in the output directory.")
        return

    print("Reading MAP file...")
    snp_names = []
    with open(map_file, "r", encoding="utf-8") as f:
        for line_num, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
            if not line:
                continue
            cols = line.split()
            if len(cols) < 4:
                raise ValueError(f"MAP line {line_num} has fewer than 4 columns.")
            snp_names.append(cols[1])
    print(f"SNPs in MAP: {len(snp_names)}")

    expected_genotype_columns = len(snp_names) * 2
    rows = []

    print("Reading PED file...")
    with open(ped_file, "r", encoding="utf-8") as f:
        for line_num, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
            if not line:
                continue

            cols = line.split()
            if len(cols) < 6:
                raise ValueError(f"PED line {line_num} has fewer than 6 mandatory columns.")

            genotype_cols = cols[6:]
            if len(genotype_cols) < expected_genotype_columns:
                raise ValueError(
                    f"PED line {line_num} has {len(genotype_cols)} genotype columns; "
                    f"expected {expected_genotype_columns}."
                )

            row = {
                "FID": cols[0],
                "IID": cols[1],
                "PID": cols[2],
                "MID": cols[3],
                "SEX": cols[4],
                "PHENO": cols[5],
            }

            trimmed = genotype_cols[:expected_genotype_columns]
            for i, snp_name in enumerate(snp_names):
                a1 = trimmed[i * 2]
                a2 = trimmed[i * 2 + 1]
                genotype = normalize_homozygous_genotype(a1, a2)
                row[snp_name] = genotype if genotype is not None else "NaN"

            rows.append(row)

    ped_df = pd.DataFrame(rows)
    print(f"Samples in PED: {len(ped_df)}")

    print("Reading phenotype file...")
    pheno_df = pd.read_csv(pheno_file, encoding="utf-8-sig")
    if "ID" not in pheno_df.columns:
        raise ValueError("Phenotype file must contain 'ID'")

    ped_df["FID"] = ped_df["FID"].astype(str).str.strip()
    pheno_df["ID"] = pheno_df["ID"].astype(str).str.strip()

    merged = ped_df.merge(
        pheno_df,
        left_on="FID",
        right_on="ID",
        how="inner"
    )

    if merged.empty:
        raise ValueError("No matching samples between PED FID and phenotype ID.")

    print(f"Merged samples: {len(merged)}")
    print("Writing PED + MAP + phenotype table...")

    with open(output_file, "w", encoding="utf-8", newline="") as f:
        merged.to_csv(f, index=False)

    print(f"Combined file written to: {output_file}")


def extract_markers(form):
    plink_chr = form.w_combo01.currentText()
    bp_lower = form.w_bp_lower.text()
    bp_upper = form.w_bp_upper.text()
    out_dir = form.w_lbl_outdir.text()
    root_name = form.w_root_name.text() or "plink"

    prev_dir = os.getcwd()
    os.chdir(out_dir)

    cmd = (
        form.plink_exe
        + " --pedmap " + form.map_ped_file
        + " --export ped"
        + " --chr " + plink_chr
        + " --from-bp " + bp_lower
        + " --to-bp " + bp_upper
        + " --out " + root_name
    )

    print("Extracting markers in current working directory: " + out_dir)
    print("Using command: " + cmd)
    os.system(cmd)

    os.chdir(prev_dir)
    print("Extraction complete")


def ld_decay(form):
    plink_chr = form.w_combo01.currentText()
    bp_lower = form.w_bp_lower.text()
    bp_upper = form.w_bp_upper.text()
    out_dir = form.w_lbl_outdir.text()
    root_name = form.w_root_name.text() or "plink"

    prev_dir = os.getcwd()
    os.chdir(out_dir)

    ld_root = root_name + "_ld_decay"

    cmd = (
        form.plink_exe
        + " --pedmap " + form.map_ped_file
        + " --chr " + plink_chr
        + " --from-bp " + bp_lower
        + " --to-bp " + bp_upper
        + " --r2-unphased"
        + " --ld-window 99999"
        + " --ld-window-kb 1000"
        + " --ld-window-r2 0"
        + " --out " + ld_root
    )

    print("Running LD decay in current working directory: " + out_dir)
    print("Using command: " + cmd)
    os.system(cmd)

    os.chdir(prev_dir)
    print("LD decay analysis complete")


def ld_decay_graph(form):
    out_dir = form.w_lbl_outdir.text()
    root_name = form.w_root_name.text() or "plink"

    vcor_file = os.path.join(out_dir, root_name + "_ld_decay.vcor")
    png_file = os.path.join(out_dir, root_name + "_ld_decay_plot.png")

    if not os.path.isfile(vcor_file):
        print("VCOR file not found: " + vcor_file)
        print("Run LD decay first.")
        return

    distances = []
    r2_values = []

    with open(vcor_file, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            pos_a = int(row["POS_A"])
            pos_b = int(row["POS_B"])
            r2 = float(row["UNPHASED_R2"])
            distances.append(abs(pos_b - pos_a))
            r2_values.append(r2)

    if not distances:
        print("No LD pairs found in the .vcor file.")
        return

    bin_size = 100
    binned = defaultdict(list)
    for dist, r2 in zip(distances, r2_values):
        binned[(dist // bin_size) * bin_size].append(r2)

    bin_x = sorted(binned.keys())
    bin_y = [sum(binned[x]) / len(binned[x]) for x in bin_x]

    if not bin_x:
        print("No binned LD values available.")
        return
    
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(bin_x, bin_y, color="red", linewidth=2)

    ax.set_xlabel("Physical distance (bp)")
    ax.set_ylabel("LD ($r^2$)")

    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0, top=1)

    ax.spines["left"].set_position(("data", 0))
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.tight_layout()
    plt.savefig(png_file, dpi=200)
    plt.show()
    
    print("LD decay graph saved to: " + png_file)
    webbrowser.open("file://" + os.path.abspath(png_file))


def potential_haplotype(form):
    out_dir = form.w_lbl_outdir.text()
    root_name = form.w_root_name.text() or "plink"

    vcor_file = os.path.join(out_dir, root_name + "_ld_decay.vcor")
    xlsx_file = os.path.join(out_dir, root_name + "_potential_haplotype.xlsx")

    if not os.path.isfile(vcor_file):
        print("VCOR file not found: " + vcor_file)
        print("Run LD decay first.")
        return

    pairs = []
    with open(vcor_file, "r", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            r2 = float(row["UNPHASED_R2"])
            pos_a = int(row["POS_A"])
            pos_b = int(row["POS_B"])
            distance_bp = abs(pos_b - pos_a)

            if distance_bp <= HAPLOTYPE_MAX_DISTANCE_BP:
                pairs.append({
                    "chrom_a": row["#CHROM_A"],
                    "pos_a": pos_a,
                    "id_a": row["ID_A"],
                    "chrom_b": row["CHROM_B"],
                    "pos_b": pos_b,
                    "id_b": row["ID_B"],
                    "r2": r2,
                    "distance_bp": distance_bp,
                })

    if not pairs:
        print(f"No SNP pairs found within {HAPLOTYPE_MAX_DISTANCE_BP} bp in the .vcor file.")
        wb = Workbook()
        ws = wb.active
        ws.title = "candidate_haplotypes"
        ws.append(["message"])
        ws.append([f"No SNP pairs found within {HAPLOTYPE_MAX_DISTANCE_BP} bp"])
        wb.save(xlsx_file)
        print("Potential haplotype Excel saved to: " + xlsx_file)
        webbrowser.open("file://" + os.path.abspath(xlsx_file))
        return

    strong_pairs = [p for p in pairs if p["r2"] >= HAPLOTYPE_LD_THRESHOLD]
    strong_pairs.sort(key=lambda x: (x["pos_a"], x["pos_b"]))

    wb = Workbook()

    ws0 = wb.active
    ws0.title = "settings"
    ws0.append(["parameter", "value"])
    ws0.append(["max_distance_bp", HAPLOTYPE_MAX_DISTANCE_BP])
    ws0.append(["ld_threshold", HAPLOTYPE_LD_THRESHOLD])
    ws0.append(["selected_bp_lower", form.w_bp_lower.text()])
    ws0.append(["selected_bp_upper", form.w_bp_upper.text()])

    ws1 = wb.create_sheet("candidate_haplotypes")
    ws1.append([
        "block_id",
        "chromosome",
        "start_bp",
        "end_bp",
        "block_size_bp",
        "num_snps",
        "mean_r2",
        "min_r2",
        "max_r2",
        "snp_ids",
    ])

    if strong_pairs:
        blocks = []
        current = None

        for pair in strong_pairs:
            block_start = min(pair["pos_a"], pair["pos_b"])
            block_end = max(pair["pos_a"], pair["pos_b"])

            if current is None:
                current = {
                    "chromosome": pair["chrom_a"],
                    "start_bp": block_start,
                    "end_bp": block_end,
                    "snps": set([pair["id_a"], pair["id_b"]]),
                    "r2_values": [pair["r2"]],
                }
                continue

            if block_start <= current["end_bp"] + HAPLOTYPE_MAX_DISTANCE_BP:
                current["end_bp"] = max(current["end_bp"], block_end)
                current["snps"].update([pair["id_a"], pair["id_b"]])
                current["r2_values"].append(pair["r2"])
            else:
                blocks.append(current)
                current = {
                    "chromosome": pair["chrom_a"],
                    "start_bp": block_start,
                    "end_bp": block_end,
                    "snps": set([pair["id_a"], pair["id_b"]]),
                    "r2_values": [pair["r2"]],
                }

        if current is not None:
            blocks.append(current)

        for idx, block in enumerate(blocks, start=1):
            mean_r2 = sum(block["r2_values"]) / len(block["r2_values"])
            ws1.append([
                idx,
                block["chromosome"],
                block["start_bp"],
                block["end_bp"],
                block["end_bp"] - block["start_bp"],
                len(block["snps"]),
                round(mean_r2, 4),
                round(min(block["r2_values"]), 4),
                round(max(block["r2_values"]), 4),
                ", ".join(sorted(block["snps"])),
            ])
    else:
        ws1.append([
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            f"No candidate haplotype blocks found at r2 >= {HAPLOTYPE_LD_THRESHOLD}",
        ])

    ws2 = wb.create_sheet("pairs_within_1000bp")
    ws2.append([
        "chrom_a",
        "pos_a",
        "id_a",
        "chrom_b",
        "pos_b",
        "id_b",
        "distance_bp",
        "r2",
        "passes_ld_threshold",
    ])

    for pair in pairs:
        ws2.append([
            pair["chrom_a"],
            pair["pos_a"],
            pair["id_a"],
            pair["chrom_b"],
            pair["pos_b"],
            pair["id_b"],
            pair["distance_bp"],
            round(pair["r2"], 4),
            "YES" if pair["r2"] >= HAPLOTYPE_LD_THRESHOLD else "NO",
        ])

    wb.save(xlsx_file)
    print("Potential haplotype Excel saved to: " + xlsx_file)
    webbrowser.open("file://" + os.path.abspath(xlsx_file))


def read_candidate_haplotypes(xlsx_path):
    df = pd.read_excel(xlsx_path, sheet_name="candidate_haplotypes")
    df = df.dropna(subset=["snp_ids"])
    df = df[df["snp_ids"].astype(str).str.strip() != ""].copy()

    blocks = []
    for _, row in df.iterrows():
        snp_ids = [x.strip() for x in str(row["snp_ids"]).split(",") if x.strip()]
        if not snp_ids:
            continue
        blocks.append({
            "block_id": row["block_id"],
            "chromosome": row.get("chromosome", ""),
            "start_bp": row.get("start_bp", ""),
            "end_bp": row.get("end_bp", ""),
            "block_size_bp": row.get("block_size_bp", ""),
            "num_snps": row.get("num_snps", len(snp_ids)),
            "mean_r2": row.get("mean_r2", ""),
            "min_r2": row.get("min_r2", ""),
            "max_r2": row.get("max_r2", ""),
            "snp_ids": snp_ids,
        })
    return blocks


def read_pairs(xlsx_path):
    df = pd.read_excel(xlsx_path, sheet_name="pairs_within_1000bp")
    needed = {"chrom_a", "pos_a", "id_a", "chrom_b", "pos_b", "id_b", "r2"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in pairs_within_1000bp: {sorted(missing)}")

    df["r2"] = pd.to_numeric(df["r2"], errors="coerce")
    df["pos_a"] = pd.to_numeric(df["pos_a"], errors="coerce")
    df["pos_b"] = pd.to_numeric(df["pos_b"], errors="coerce")
    if "distance_bp" in df.columns:
        df["distance_bp"] = pd.to_numeric(df["distance_bp"], errors="coerce")
    else:
        df["distance_bp"] = (df["pos_b"] - df["pos_a"]).abs()

    df = df.dropna(subset=["r2", "pos_a", "pos_b", "distance_bp"]).copy()
    df = df[
        (df["r2"] >= R2_THRESHOLD) &
        (df["distance_bp"] <= HAPLOTYPE_MAX_DISTANCE_BP)
    ].copy()

    df["pair_key"] = df.apply(
        lambda r: f"{r['chrom_a']}:{int(r['pos_a'])}__{r['chrom_b']}:{int(r['pos_b'])}",
        axis=1
    )
    df = df.drop_duplicates(subset=["pair_key"]).copy()
    return df


def read_map(map_path):
    return pd.read_csv(map_path, sep=r"\s+", header=None, names=["chr", "snp_id", "cm", "bp"])


def build_bp_to_index(map_df):
    bp_to_info = {}
    for idx, row in map_df.iterrows():
        bp_to_info[int(row["bp"])] = {
            "index": idx,
            "snp_id": row["snp_id"],
            "chr": row["chr"],
        }
    return bp_to_info


def read_ped_subset_raw(ped_path, snp_index):
    sorted_snps = sorted(snp_index.items(), key=lambda x: x[1])
    snp_names = [x[0] for x in sorted_snps]
    snp_positions = [x[1] for x in sorted_snps]

    subset_rows = []
    with open(ped_path, "r", encoding="utf-8") as f:
        for line_num, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
            if not line:
                continue

            cols = line.split()
            if len(cols) < 6:
                raise ValueError(f"PED line {line_num} has fewer than 6 mandatory columns.")

            genotype_cols = cols[6:]
            row = {
                "FID": str(cols[0]).strip(),
                "IID": str(cols[1]).strip(),
            }

            for snp_name, snp_pos in zip(snp_names, snp_positions):
                a1_col = snp_pos * 2
                a2_col = a1_col + 1
                if a2_col >= len(genotype_cols):
                    row[snp_name] = None
                    continue

                a1 = genotype_cols[a1_col]
                a2 = genotype_cols[a2_col]
                row[snp_name] = normalize_homozygous_genotype(a1, a2)

            subset_rows.append(row)

    return pd.DataFrame(subset_rows)


def build_haplotype_code(row, snp_ids):
    vals = []
    for snp in snp_ids:
        gt = row.get(snp)
        if pd.isna(gt) or gt is None:
            return None
        vals.append(str(gt))
    return ",".join(vals)


def summarize_group(values):
    n = len(values)
    mean_value = sum(values) / n
    if n == 1:
        std_dev = 0.0
    else:
        variance = sum((value - mean_value) ** 2 for value in values) / (n - 1)
        std_dev = math.sqrt(variance)
    return n, mean_value, std_dev


def one_way_anova(groups):
    groups = [g for g in groups if len(g) > 0]
    k = len(groups)
    total_n = sum(len(g) for g in groups)

    if k < 2 or total_n <= k:
        return "", "", "", "", ""

    group_means = [sum(g) / len(g) for g in groups]
    grand_mean = sum(sum(g) for g in groups) / total_n

    ss_between = sum(len(g) * (gm - grand_mean) ** 2 for g, gm in zip(groups, group_means))
    ss_within = 0.0
    for g, gm in zip(groups, group_means):
        for x in g:
            ss_within += (x - gm) ** 2
    ss_total = ss_between + ss_within

    df_between = k - 1
    df_within = total_n - k
    if df_within <= 0:
        return "", "", "", "", ""

    ms_between = ss_between / df_between
    ms_within = ss_within / df_within
    anova_r2 = ""
    if ss_total > 0:
        anova_r2 = ss_between / ss_total
    if ms_within == 0:
        return "", df_between, df_within, "", anova_r2

    f_stat = ms_between / ms_within
    p_value = ""
    if scipy_f_distribution is not None:
        p_value = scipy_f_distribution.sf(f_stat, df_between, df_within)

    return f_stat, df_between, df_within, p_value, anova_r2


def _load_phenotype_and_merge(ped_df, pheno_file, trait_column):
    pheno_df = pd.read_csv(pheno_file, encoding="utf-8-sig")
    if "ID" not in pheno_df.columns:
        raise ValueError("Phenotype file must contain 'ID'")
    if trait_column not in pheno_df.columns:
        raise ValueError("Phenotype file must contain '{}'".format(trait_column))

    ped_df["FID"] = ped_df["FID"].astype(str).str.strip()
    pheno_df["ID"] = pheno_df["ID"].astype(str).str.strip()
    pheno_df[trait_column] = pd.to_numeric(pheno_df[trait_column], errors="coerce")

    merged = ped_df.merge(
        pheno_df[["ID", trait_column]],
        left_on="FID",
        right_on="ID",
        how="inner"
    ).copy()

    return merged


def formal_haplotype_block_detection(form):
    """
    Formal LD block detection from the strong pair table using connected components.
    Output:
    - settings
    - formal_haplotype_blocks
    - pairs_used
    """
    out_dir = form.w_lbl_outdir.text()
    root_name = form.w_root_name.text() or "plink"

    hap_xlsx = os.path.join(out_dir, root_name + "_potential_haplotype.xlsx")
    output_xlsx = os.path.join(out_dir, root_name + "_FHBD_results.xlsx")

    if not os.path.isfile(hap_xlsx):
        print("Potential haplotype Excel not found: " + hap_xlsx)
        print("Run Potential Haplotype first.")
        return

    print("Reading pairwise LD data for formal block detection...")
    pair_df = pd.read_excel(hap_xlsx, sheet_name="pairs_within_1000bp")
    needed = {"chrom_a", "pos_a", "id_a", "chrom_b", "pos_b", "id_b", "r2"}
    missing = needed - set(pair_df.columns)
    if missing:
        raise ValueError(f"Missing columns in pairs_within_1000bp: {sorted(missing)}")

    pair_df["r2"] = pd.to_numeric(pair_df["r2"], errors="coerce")
    pair_df["pos_a"] = pd.to_numeric(pair_df["pos_a"], errors="coerce")
    pair_df["pos_b"] = pd.to_numeric(pair_df["pos_b"], errors="coerce")
    if "distance_bp" in pair_df.columns:
        pair_df["distance_bp"] = pd.to_numeric(pair_df["distance_bp"], errors="coerce")
    else:
        pair_df["distance_bp"] = (pair_df["pos_b"] - pair_df["pos_a"]).abs()

    pair_df = pair_df.dropna(subset=["r2", "pos_a", "pos_b", "distance_bp"]).copy()
    pair_df["chrom_a"] = pair_df["chrom_a"].astype(str).str.strip()
    pair_df["chrom_b"] = pair_df["chrom_b"].astype(str).str.strip()
    pair_df["id_a"] = pair_df["id_a"].astype(str).str.strip()
    pair_df["id_b"] = pair_df["id_b"].astype(str).str.strip()

    strong_df = pair_df[
        (pair_df["r2"] >= HAPLOTYPE_LD_THRESHOLD) &
        (pair_df["distance_bp"] <= HAPLOTYPE_MAX_DISTANCE_BP) &
        (pair_df["chrom_a"] == pair_df["chrom_b"])
    ].copy()

    wb = Workbook()

    ws0 = wb.active
    ws0.title = "settings"
    ws0.append(["parameter", "value"])
    ws0.append(["method", "formal_connected_components"])
    ws0.append(["max_distance_bp", HAPLOTYPE_MAX_DISTANCE_BP])
    ws0.append(["ld_threshold", HAPLOTYPE_LD_THRESHOLD])
    ws0.append(["selected_bp_lower", form.w_bp_lower.text()])
    ws0.append(["selected_bp_upper", form.w_bp_upper.text()])

    ws1 = wb.create_sheet("formal_haplotype_blocks")
    ws1.append([
        "block_id",
        "chromosome",
        "start_bp",
        "end_bp",
        "block_size_bp",
        "num_snps",
        "num_pairs",
        "mean_r2",
        "min_r2",
        "max_r2",
        "snp_ids",
    ])

    ws2 = wb.create_sheet("pairs_used")
    ws2.append([
        "block_id",
        "chromosome",
        "id_a",
        "pos_a",
        "id_b",
        "pos_b",
        "distance_bp",
        "r2",
    ])

    if strong_df.empty:
        ws1.append(["message"])
        ws1.append([f"No strong LD pairs found at r2 >= {HAPLOTYPE_LD_THRESHOLD}"])
        ws2.append(["message"])
        ws2.append([f"No strong LD pairs found at r2 >= {HAPLOTYPE_LD_THRESHOLD}"])
        wb.save(output_xlsx)
        print("Formal haplotype block Excel saved to: " + output_xlsx)
        webbrowser.open("file://" + os.path.abspath(output_xlsx))
        return

    snp_meta = {}
    all_snps = set()

    for _, row in strong_df.iterrows():
        for side in ("a", "b"):
            snp = row[f"id_{side}"]
            pos = int(row[f"pos_{side}"])
            chrom = str(row[f"chrom_{side}"])
            all_snps.add(snp)
            if snp not in snp_meta:
                snp_meta[snp] = {"chrom": chrom, "pos": pos}

    parent = {snp: snp for snp in all_snps}
    rank = {snp: 0 for snp in all_snps}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra = find(a)
        rb = find(b)
        if ra == rb:
            return
        if rank[ra] < rank[rb]:
            parent[ra] = rb
        elif rank[ra] > rank[rb]:
            parent[rb] = ra
        else:
            parent[rb] = ra
            rank[ra] += 1

    for _, row in strong_df.iterrows():
        union(row["id_a"], row["id_b"])

    components = defaultdict(list)
    for snp in all_snps:
        components[find(snp)].append(snp)

    component_list = []
    for _, comp_snps in components.items():
        comp_snps_sorted = sorted(comp_snps, key=lambda s: snp_meta[s]["pos"])
        first_snp = comp_snps_sorted[0]
        component_list.append((
            snp_meta[first_snp]["chrom"],
            snp_meta[first_snp]["pos"],
            comp_snps_sorted,
        ))

    def _chrom_sort_key(value):
        text = str(value).strip()
        try:
            return (0, int(text))
        except ValueError:
            return (1, text)

    component_list.sort(key=lambda x: (_chrom_sort_key(x[0]), int(x[1])))

    block_rows = []
    edge_rows = []

    block_id = 0
    for chrom, _, comp_snps_sorted in component_list:
        if len(comp_snps_sorted) < 2:
            continue

        comp_set = set(comp_snps_sorted)
        comp_edges = strong_df[
            strong_df["id_a"].isin(comp_set) &
            strong_df["id_b"].isin(comp_set)
        ].copy()

        if comp_edges.empty:
            continue

        block_id += 1
        start_bp = min(snp_meta[s]["pos"] for s in comp_snps_sorted)
        end_bp = max(snp_meta[s]["pos"] for s in comp_snps_sorted)
        r2_values = comp_edges["r2"].astype(float).tolist()

        block_rows.append({
            "block_id": block_id,
            "chromosome": chrom,
            "start_bp": start_bp,
            "end_bp": end_bp,
            "block_size_bp": end_bp - start_bp,
            "num_snps": len(comp_snps_sorted),
            "num_pairs": len(comp_edges),
            "mean_r2": round(sum(r2_values) / len(r2_values), 4),
            "min_r2": round(min(r2_values), 4),
            "max_r2": round(max(r2_values), 4),
            "snp_ids": ", ".join(comp_snps_sorted),
        })

        for _, row in comp_edges.iterrows():
            edge_rows.append({
                "block_id": block_id,
                "chromosome": chrom,
                "id_a": row["id_a"],
                "pos_a": int(row["pos_a"]),
                "id_b": row["id_b"],
                "pos_b": int(row["pos_b"]),
                "distance_bp": round(float(row["distance_bp"]), 4),
                "r2": round(float(row["r2"]), 4),
            })

    if block_rows:
        for row in block_rows:
            ws1.append([
                row["block_id"],
                row["chromosome"],
                row["start_bp"],
                row["end_bp"],
                row["block_size_bp"],
                row["num_snps"],
                row["num_pairs"],
                row["mean_r2"],
                row["min_r2"],
                row["max_r2"],
                row["snp_ids"],
            ])
    else:
        ws1.append(["message"])
        ws1.append(["No formal haplotype blocks were detected."])

    if edge_rows:
        for row in edge_rows:
            ws2.append([
                row["block_id"],
                row["chromosome"],
                row["id_a"],
                row["pos_a"],
                row["id_b"],
                row["pos_b"],
                row["distance_bp"],
                row["r2"],
            ])
    else:
        ws2.append(["message"])
        ws2.append(["No pairs were assigned to any formal block."])

    wb.save(output_xlsx)
    print("Formal haplotype block Excel saved to: " + output_xlsx)
    webbrowser.open("file://" + os.path.abspath(output_xlsx))


def two_snp_pair_anova(form):
    """Run homozygote-only two-locus genotype-combination ANOVA."""
    out_dir = form.w_lbl_outdir.text()
    root_name = form.w_root_name.text() or "plink"

    ped_file = os.path.join(out_dir, root_name + ".ped")
    map_file = os.path.join(out_dir, root_name + ".map")
    pairs_xlsx = os.path.join(out_dir, root_name + "_potential_haplotype.xlsx")
    output_xlsx = os.path.join(out_dir, root_name + "_two_snp_pair_anova_results.xlsx")
    trait_column = "AWD_Fe57_grain_Y1"

    for path, description in [
        (ped_file, "PED file"),
        (map_file, "MAP file"),
        (pairs_xlsx, "Candidate LD-region Excel file"),
    ]:
        if not os.path.isfile(path):
            print(description + " not found: " + path)
            return

    if scipy_f_distribution is None:
        raise ImportError("SciPy is required to calculate ANOVA P values.")

    print("Reading qualifying SNP pairs...")
    pair_df = read_pairs(pairs_xlsx)
    if pair_df.empty:
        raise ValueError(f"No SNP pairs found with r2 >= {R2_THRESHOLD}")
    print(f"Pairs with r2 >= {R2_THRESHOLD}: {len(pair_df)}")

    map_df = read_map(map_file)
    bp_to_info = build_bp_to_index(map_df)
    needed_bps = set(pair_df["pos_a"].astype(int).tolist() + pair_df["pos_b"].astype(int).tolist())
    missing = [bp for bp in needed_bps if bp not in bp_to_info]
    if missing:
        raise ValueError(f"These BP positions were not found in MAP: {missing[:20]}")

    sorted_needed = sorted(
        [(bp, bp_to_info[bp]["snp_id"], bp_to_info[bp]["index"]) for bp in needed_bps],
        key=lambda value: value[2],
    )
    subset_rows = []
    with open(ped_file, "r", encoding="utf-8") as ped_handle:
        for line_num, raw_line in enumerate(ped_handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            cols = line.split()
            if len(cols) < 6:
                raise ValueError(f"PED line {line_num} has fewer than 6 mandatory columns.")
            genotype_cols = cols[6:]
            row = {"FID": str(cols[0]).strip(), "IID": str(cols[1]).strip()}
            for _, snp_id, snp_idx in sorted_needed:
                allele_1_index = snp_idx * 2
                allele_2_index = allele_1_index + 1
                if allele_2_index >= len(genotype_cols):
                    row[snp_id] = None
                else:
                    row[snp_id] = normalize_homozygous_genotype(
                        genotype_cols[allele_1_index], genotype_cols[allele_2_index]
                    )
            subset_rows.append(row)

    ped_df = pd.DataFrame(subset_rows)
    phenotype_file = os.path.join(out_dir, "Fe_Phene.csv")
    if not os.path.isfile(phenotype_file):
        phenotype_file = os.path.join(os.path.dirname(out_dir), "Fe_Phene.csv")
    if not os.path.isfile(phenotype_file):
        print("Phenotype file not found: " + phenotype_file)
        return

    merged = _load_phenotype_and_merge(ped_df, phenotype_file, trait_column)
    merged = merged.dropna(subset=[trait_column]).copy()
    phenotype_eligible_samples = len(merged)
    if phenotype_eligible_samples == 0:
        raise ValueError("No matched accessions have a non-missing phenotype value.")
    required_homozygous_samples = math.ceil(
        HOMOZYGOUS_CALL_RATE_THRESHOLD * phenotype_eligible_samples
    )
    print(f"Phenotype-eligible samples: {phenotype_eligible_samples}")
    print(
        "Minimum homozygous calls required per locus and pair: "
        f"{required_homozygous_samples} ({HOMOZYGOUS_CALL_RATE_THRESHOLD:.0%})"
    )

    locus_call_counts = {
        snp_id: int(merged[snp_id].notna().sum()) for _, snp_id, _ in sorted_needed
    }
    result_rows = []
    excluded_locus_threshold = 0
    excluded_joint_threshold = 0
    excluded_group_threshold = 0

    total_pairs = len(pair_df)
    for pair_number, (_, pair_row) in enumerate(pair_df.iterrows(), start=1):
        if pair_number == 1 or pair_number % 100 == 0:
            print(f"Analyzing pair {pair_number}/{total_pairs}")

        id_a = str(pair_row["id_a"])
        id_b = str(pair_row["id_b"])
        if (
            locus_call_counts[id_a] < required_homozygous_samples
            or locus_call_counts[id_b] < required_homozygous_samples
        ):
            excluded_locus_threshold += 1
            continue

        pair_data = merged[["FID", "IID", trait_column, id_a, id_b]].dropna(
            subset=[id_a, id_b]
        ).copy()
        joint_homozygous_samples = len(pair_data)
        if joint_homozygous_samples < required_homozygous_samples:
            excluded_joint_threshold += 1
            continue

        pair_data["genotype_combination"] = (
            pair_data[id_a].astype(str) + "," + pair_data[id_b].astype(str)
        )
        group_counts = pair_data["genotype_combination"].value_counts()
        valid_groups = group_counts[group_counts >= MIN_GROUP_SIZE].index.tolist()
        pair_data = pair_data[
            pair_data["genotype_combination"].isin(valid_groups)
        ].copy()
        if pair_data["genotype_combination"].nunique() < 2:
            excluded_group_threshold += 1
            continue

        genotype_groups = {
            genotype: group[trait_column].astype(float).tolist()
            for genotype, group in pair_data.groupby("genotype_combination")
            if len(group) >= MIN_GROUP_SIZE
        }
        if len(genotype_groups) < 2:
            excluded_group_threshold += 1
            continue

        ordered_groups = sorted(genotype_groups)
        f_stat, df_between, df_within, p_value, anova_r2 = one_way_anova(
            [genotype_groups[group] for group in ordered_groups]
        )
        if p_value == "":
            continue

        group_summary = []
        for genotype in ordered_groups:
            sample_count, mean_value, sd_value = summarize_group(genotype_groups[genotype])
            group_summary.append(
                f"{genotype} (n={sample_count}, mean={mean_value:.4f}, sd={sd_value:.4f})"
            )

        pos_a = int(pair_row["pos_a"])
        pos_b = int(pair_row["pos_b"])
        result_rows.append({
            "chrom_a": str(pair_row["chrom_a"]),
            "pos_a": pos_a,
            "id_a": id_a,
            "chrom_b": str(pair_row["chrom_b"]),
            "pos_b": pos_b,
            "id_b": id_b,
            "start_bp": min(pos_a, pos_b),
            "end_bp": max(pos_a, pos_b),
            "distance_bp": int(pair_row["distance_bp"]),
            "LD_r2": float(pair_row["r2"]),
            "trait": trait_column,
            "joint_homozygous_call_count": joint_homozygous_samples,
            "joint_homozygous_call_rate": joint_homozygous_samples / phenotype_eligible_samples,
            "non_missing_samples_after_group_filter": len(pair_data),
            "genotype_combination_groups": len(genotype_groups),
            "observed_genotype_combinations": "; ".join(ordered_groups),
            "group_summary": "; ".join(group_summary),
            "f_statistic": round(f_stat, 6),
            "ANOVA_R2": round(anova_r2, 6),
            "df_between": df_between,
            "df_within": df_within,
            "p_value": p_value,
            "nominal_significant": "YES" if p_value < SIGNIFICANCE_THRESHOLD else "NO",
        })

    results_df = pd.DataFrame(result_rows)
    adjusted_df = add_multiple_testing_results(results_df, "p_value")
    if adjusted_df.empty:
        raise ValueError("No SNP pairs remained eligible for analysis.")

    base_columns = list(results_df.columns)
    bonferroni_df = adjusted_df[
        base_columns + ["Bonferroni_threshold", "Bonferroni_significant"]
    ].copy()
    bh_df = adjusted_df[
        base_columns + ["Rank", "BH_q_raw", "BH_q_adjusted", "BH_significant"]
    ].copy()
    results_sorted = adjusted_df[base_columns].copy()

    settings_df = pd.DataFrame([
        ["heterozygous_call_handling", "treated as missing"],
        ["homozygous_call_rate_threshold", HOMOZYGOUS_CALL_RATE_THRESHOLD],
        ["phenotype_eligible_samples", phenotype_eligible_samples],
        ["minimum_homozygous_samples", required_homozygous_samples],
        ["minimum_genotype_group_size", MIN_GROUP_SIZE],
        ["qualifying_LD_pairs_input", total_pairs],
        ["pairs_excluded_by_individual_locus_threshold", excluded_locus_threshold],
        ["pairs_excluded_by_joint_threshold", excluded_joint_threshold],
        ["pairs_excluded_by_group_requirement", excluded_group_threshold],
        ["pairs_tested", len(results_sorted)],
    ], columns=["parameter", "value"])

    with pd.ExcelWriter(output_xlsx, engine="openpyxl") as writer:
        settings_df.to_excel(writer, sheet_name="settings", index=False)
        results_sorted.to_excel(
            writer, sheet_name="pairwise_genotype_combinations", index=False
        )
        bonferroni_df.to_excel(writer, sheet_name="Bonferroni correction", index=False)
        bh_df.to_excel(writer, sheet_name="Benjamini-Hochberg", index=False)

    print(f"Pairs tested: {len(results_sorted)}")
    print(f"ANOVA results written to: {output_xlsx}")

def snp_trait_analysis(form):
    """
    Merge PED, MAP, and Fe_Phene.csv from the selected output directory, then
    run SNP-by-trait ANOVA for the phenotype column AWD_Fe57_grain_Y1.

    Outputs:
    - <root>_genotype_phenotype_merged.csv
    - <root>_snp_trait_anova_results.csv
    - <root>_anova_significant_boxplots/
    """
    out_dir = form.w_lbl_outdir.text()
    root_name = form.w_root_name.text() or "plink"

    ped_file = os.path.join(out_dir, root_name + ".ped")
    map_file = os.path.join(out_dir, root_name + ".map")
    phenotype_file = os.path.join(out_dir, "Fe_Phene.csv")
    if not os.path.isfile(phenotype_file):
        phenotype_file = os.path.join(os.path.dirname(out_dir), "Fe_Phene.csv")

    merged_file = os.path.join(out_dir, root_name + "_genotype_phenotype_merged.csv")
    results_file = os.path.join(out_dir, root_name + "_snp_trait_anova_results.xlsx")
    boxplot_dir = os.path.join(out_dir, root_name + "_anova_significant_boxplots")

    trait_column = "AWD_Fe57_grain_Y1"
    missing_genotypes = {"", "NaN", "NA", "nan", "None"}
    significance_threshold = 0.05

    if not os.path.isfile(ped_file):
        print("PED file not found: " + ped_file)
        return
    if not os.path.isfile(map_file):
        print("MAP file not found: " + map_file)
        return
    if not os.path.isfile(phenotype_file):
        print("Phenotype file not found: " + phenotype_file)
        print("Expected Fe_Phene.csv in the output directory or its parent folder.")
        return

    def read_map_file(map_path):
        snp_names = []
        with open(map_path, "r", encoding="utf-8") as map_handle:
            for line_number, raw_line in enumerate(map_handle, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                cols = line.split()
                if len(cols) < 2:
                    raise ValueError("MAP file line {} does not contain at least 2 columns.".format(line_number))
                snp_names.append(cols[1])
        return snp_names

    def build_genotype_rows(ped_path, snp_names):
        genotype_rows = []
        expected_genotype_columns = len(snp_names) * 2
        with open(ped_path, "r", encoding="utf-8") as ped_handle:
            for line_number, raw_line in enumerate(ped_handle, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                cols = line.split()
                if len(cols) < 6:
                    raise ValueError("PED file line {} has fewer than 6 required PED columns.".format(line_number))

                genotype_cols = cols[6:]
                if len(genotype_cols) < expected_genotype_columns:
                    raise ValueError(
                        "PED file line {} has {} genotype columns, but the MAP file requires {}."
                        .format(line_number, len(genotype_cols), expected_genotype_columns)
                    )

                row = {
                    "FID": cols[0],
                    "IID": cols[1],
                    "PID": cols[2],
                    "MID": cols[3],
                    "SEX": cols[4],
                    "PHENO": cols[5],
                }

                trimmed_genotype_cols = genotype_cols[:expected_genotype_columns]
                for index, snp_name in enumerate(snp_names):
                    genotype = normalize_homozygous_genotype(
                        trimmed_genotype_cols[index * 2],
                        trimmed_genotype_cols[index * 2 + 1],
                    )
                    row[snp_name] = genotype if genotype is not None else "NaN"
                genotype_rows.append(row)

        return ["FID", "IID", "PID", "MID", "SEX", "PHENO"] + snp_names, genotype_rows

    def read_phenotype_rows(phenotype_path):
        with open(phenotype_path, "r", encoding="utf-8-sig", newline="") as phenotype_handle:
            reader = csv.DictReader(phenotype_handle)
            if reader.fieldnames is None:
                raise ValueError("Phenotype file does not contain a header row.")
            return list(reader.fieldnames), list(reader)

    def merge_rows(genotype_headers, genotype_rows, phenotype_headers, phenotype_rows):
        if "ID" not in phenotype_headers:
            raise ValueError("Fe_Phene.csv must contain a column named 'ID'.")
        if trait_column not in phenotype_headers:
            raise ValueError("Fe_Phene.csv must contain a column named '{}'.".format(trait_column))

        phenotype_by_id = {}
        for row in phenotype_rows:
            row_id = row.get("ID", "").strip()
            if row_id:
                phenotype_by_id[row_id] = row

        phenotype_headers_to_add = [header for header in phenotype_headers if header != "ID"]
        merged_rows = []
        for genotype_row in genotype_rows:
            phenotype_row = phenotype_by_id.get(genotype_row["FID"].strip(), {})
            merged_row = dict(genotype_row)
            for header in phenotype_headers_to_add:
                merged_row[header] = phenotype_row.get(header, "")
            merged_rows.append(merged_row)

        return genotype_headers + phenotype_headers_to_add, merged_rows

    def write_csv_dicts(output_path, headers, rows):
        with open(output_path, "w", encoding="utf-8", newline="") as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=headers)
            writer.writeheader()
            writer.writerows(rows)

    def detect_snp_columns(headers):
        excluded = {"FID", "IID", "PID", "MID", "SEX", "PHENO", "ID", trait_column}
        return [header for header in headers if header not in excluded]

    def summarize_group(values):
        n = len(values)
        mean_value = sum(values) / n
        if n == 1:
            std_dev = 0.0
        else:
            variance = sum((value - mean_value) ** 2 for value in values) / (n - 1)
            std_dev = math.sqrt(variance)
        return n, mean_value, std_dev

    def one_way_anova(groups):
        filtered_groups = [group for group in groups if len(group) > 0]
        k = len(filtered_groups)
        total_n = sum(len(group) for group in filtered_groups)
        if k < 2 or total_n <= k:
            return "", "", "", "", ""

        group_means = [sum(group) / len(group) for group in filtered_groups]
        grand_mean = sum(sum(group) for group in filtered_groups) / total_n
        ss_between = sum(len(group) * (group_mean - grand_mean) ** 2 for group, group_mean in zip(filtered_groups, group_means))
        ss_within = 0.0
        for group, group_mean in zip(filtered_groups, group_means):
            for value in group:
                ss_within += (value - group_mean) ** 2
        ss_total = ss_between + ss_within

        df_between = k - 1
        df_within = total_n - k
        if df_within <= 0:
            return "", "", "", "", ""

        ms_between = ss_between / df_between
        ms_within = ss_within / df_within
        r2 = ""
        if ss_total > 0:
            r2 = ss_between / ss_total
        if ms_within == 0:
            return "", df_between, df_within, "", r2

        f_stat = ms_between / ms_within
        p_value = ""
        if scipy_f_distribution is not None:
            p_value = scipy_f_distribution.sf(f_stat, df_between, df_within)
        return f_stat, df_between, df_within, p_value, r2

    def analyze_snps(rows, snp_columns, trait_col):
        phenotype_eligible_rows = []
        for row in rows:
            try:
                phenotype_value = float(str(row.get(trait_col, "")).strip())
            except Exception:
                continue
            phenotype_eligible_rows.append((row, phenotype_value))

        phenotype_eligible_samples = len(phenotype_eligible_rows)
        if phenotype_eligible_samples == 0:
            raise ValueError("No accessions have a matched, non-missing phenotype value.")
        required_homozygous_samples = math.ceil(
            HOMOZYGOUS_CALL_RATE_THRESHOLD * phenotype_eligible_samples
        )
        results = []
        for snp in snp_columns:
            genotype_groups = {}
            homozygous_call_count = 0
            for row, phenotype_value in phenotype_eligible_rows:
                genotype = str(row.get(snp, "")).strip()
                if genotype in missing_genotypes:
                    continue
                homozygous_call_count += 1
                genotype_groups.setdefault(genotype, []).append(phenotype_value)

            if homozygous_call_count < required_homozygous_samples:
                continue

            genotype_groups = {
                genotype: values
                for genotype, values in genotype_groups.items()
                if len(values) >= MIN_GROUP_SIZE
            }
            if len(genotype_groups) < 2:
                continue

            ordered_genotypes = sorted(genotype_groups)
            group_summary = []
            for genotype in ordered_genotypes:
                count, mean_value, std_dev = summarize_group(genotype_groups[genotype])
                group_summary.append(
                    "{}(n={},mean={:.4f},sd={:.4f})".format(
                        genotype, count, mean_value, std_dev
                    )
                )

            f_stat, df_between, df_within, p_value, r2 = one_way_anova(
                [genotype_groups[genotype] for genotype in ordered_genotypes]
            )
            if p_value == "":
                continue
            results.append({
                "SNP": snp,
                "Trait": trait_col,
                "HomozygousCallCount": homozygous_call_count,
                "HomozygousCallRate": homozygous_call_count / phenotype_eligible_samples,
                "NonMissingSamplesAfterGroupFilter": sum(
                    len(values) for values in genotype_groups.values()
                ),
                "GenotypeGroups": len(genotype_groups),
                "ObservedGenotypes": "; ".join(ordered_genotypes),
                "GroupSummary": "; ".join(group_summary),
                "FStatistic": f_stat,
                "ANOVA_R2": r2,
                "DF_Between": df_between,
                "DF_Within": df_within,
                "PValue": p_value,
                "NominalSignificant": "YES" if p_value < significance_threshold else "NO",
            })
        return results, phenotype_eligible_samples, required_homozygous_samples

    def safe_file_name(text):
        return re.sub(r"[^A-Za-z0-9._-]+", "_", text)

    def create_significant_boxplots(rows, results, trait_col):
        os.makedirs(boxplot_dir, exist_ok=True)
        plot_count = 0
        for result in results:
            p_value = result.get("PValue", "").strip()
            if not p_value:
                continue
            if float(p_value) >= significance_threshold:
                continue

            snp_name = result["SNP"]
            genotype_groups = {}
            for row in rows:
                genotype = str(row.get(snp_name, "")).strip()
                trait_value = str(row.get(trait_col, "")).strip()
                if genotype in missing_genotypes:
                    continue
                try:
                    genotype_groups.setdefault(genotype, []).append(float(trait_value))
                except Exception:
                    continue

            genotype_groups = {
                genotype: values
                for genotype, values in genotype_groups.items()
                if len(values) >= MIN_GROUP_SIZE
            }

            if len(genotype_groups) < 2:
                continue

            ordered_genotypes = sorted(genotype_groups.keys())
            plt.figure(figsize=(8, 6))
            plt.boxplot([genotype_groups[g] for g in ordered_genotypes], labels=ordered_genotypes)
            plt.title("ANOVA: {} vs {}".format(trait_col, snp_name))
            plt.xlabel("Genotype")
            plt.ylabel(trait_col)
            plt.tight_layout()
            plt.savefig(os.path.join(boxplot_dir, "{}_boxplot.png".format(safe_file_name(snp_name))), dpi=300)
            plt.close()
            plot_count += 1
        return plot_count

    print("Reading MAP file...")
    snp_names = read_map_file(map_file)
    genotype_headers, genotype_rows = build_genotype_rows(ped_file, snp_names)
    phenotype_headers, phenotype_rows = read_phenotype_rows(phenotype_file)
    merged_headers, merged_rows = merge_rows(genotype_headers, genotype_rows, phenotype_headers, phenotype_rows)

    write_csv_dicts(merged_file, merged_headers, merged_rows)

    if scipy_f_distribution is None:
        raise ImportError("SciPy is required to calculate ANOVA P values.")

    snp_columns = detect_snp_columns(merged_headers)
    results, phenotype_eligible_samples, required_homozygous_samples = analyze_snps(
        merged_rows, snp_columns, trait_column
    )
    results_df = pd.DataFrame(results)
    adjusted_df = add_multiple_testing_results(results_df, "PValue")
    if adjusted_df.empty:
        raise ValueError("No SNPs remained eligible for analysis.")

    base_columns = list(results_df.columns)
    results_sorted = adjusted_df[base_columns].copy()
    bonferroni_df = adjusted_df[
        base_columns + ["Bonferroni_threshold", "Bonferroni_significant"]
    ].copy()
    bh_df = adjusted_df[
        base_columns + ["Rank", "BH_q_raw", "BH_q_adjusted", "BH_significant"]
    ].copy()
    settings_df = pd.DataFrame([
        ["heterozygous_call_handling", "treated as missing"],
        ["homozygous_call_rate_threshold", HOMOZYGOUS_CALL_RATE_THRESHOLD],
        ["phenotype_eligible_samples", phenotype_eligible_samples],
        ["minimum_homozygous_samples", required_homozygous_samples],
        ["minimum_genotype_group_size", MIN_GROUP_SIZE],
        ["SNPs_tested", len(results_sorted)],
    ], columns=["parameter", "value"])

    with pd.ExcelWriter(results_file, engine="openpyxl") as writer:
        settings_df.to_excel(writer, sheet_name="settings", index=False)
        results_sorted.to_excel(writer, sheet_name="Exploratory_SNP_trait", index=False)
        bonferroni_df.to_excel(writer, sheet_name="Bonferroni correction", index=False)
        bh_df.to_excel(writer, sheet_name="Benjamini-Hochberg", index=False)

    plot_count = 0
    print("Merged file: {}".format(merged_file))
    print("ANOVA results: {}".format(results_file))
    print("Heterozygous calls were treated as missing.")
    print(
        "Minimum homozygous calls required: {} ({:.0%} of {} phenotype-eligible accessions)"
        .format(
            required_homozygous_samples,
            HOMOZYGOUS_CALL_RATE_THRESHOLD,
            phenotype_eligible_samples,
        )
    )
    print("SNPs analyzed: {}".format(len(results_sorted)))
    print("Automatic SNP box plots disabled.")

