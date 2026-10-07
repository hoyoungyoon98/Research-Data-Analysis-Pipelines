# ============================================================
# PLAQUE ASSAY AUTOMATED ANALYSIS PIPELINE
# ============================================================

from pathlib import Path
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import dunnett
from statsmodels.formula.api import ols


# ============================================================
# STEP 0. USER CONFIGURATION
# ============================================================

# Base folder directory
BASE_FOLDER = Path(__file__).resolve().parent / "example_data"

# Excel sheet to analyze.
TARGET_SHEET = "DATASET1"

# Sample order
SAMPLE_ORDER = [
    "Control",
    "Sample_A",
    "Sample_B"
]

# Statistics
CONTROL = "Control"

# Plaque assay setting
PLATING_VOLUME_ML = 0.5

PLAQUE_ASSAY_LOD = 2.0

# Graph
GRAPH_TITLE = "Title"

XLABEL = "X-label"

YLABEL = "Y-label"

SHOW_SIGNIFICANCE = True

SHOW_LOD = True

SHOW_LEGEND = True


# ============================================================
# COMMON FUNCTIONS
# ============================================================

# Significance
def significance_star(p_value):

    if pd.isna(p_value):
        return "NA"

    if p_value < 0.001:
        return "***"

    elif p_value < 0.01:
        return "**"

    elif p_value < 0.05:
        return "*"

    else:
        return "ns"


# ============================================================
# Data parsing helpers
# ============================================================
# Dilution factor
def parse_dilution_factor(value):

    if pd.isna(value):
        return np.nan

    value = str(value).strip()

    if value == "":
        return np.nan

    # Neat 
    if value.lower() == "neat":
        return 1.0

    # Fraction 
    fraction_match = re.fullmatch(r"1\s*/\s*(\d+)", value)

    if fraction_match:

        denominator = float(fraction_match.group(1))

        if denominator == 0:
            return np.nan

        return 1.0 / denominator

    # Exponent 
    exponent_match = re.fullmatch(r"-(\d+)", value)

    if exponent_match:

        exponent = float(exponent_match.group(1))

        return 10 ** (-exponent)

    # Unknown
    return np.nan

# hpi
def extract_hpi(series):

    hpi = (
        series
        .astype(str)
        .str.extract(r"(\d+)")[0]
        .astype(float)
    )

    return hpi


# ============================================================
# EXPERIMENT MANAGEMENT FUNCTIONS
# ============================================================

# ============================================================
# EXPERIMENT FOLDER SEARCH
# ============================================================

def find_experiment_folders(base_folder):

    if not base_folder.exists():

        raise FileNotFoundError(
            f"Data folder not found: {base_folder}"
        )

    experiment_folders = sorted(
        [
            folder 
            for folder in base_folder.iterdir()
            if folder.is_dir()
            and folder.name.startswith("Experiment_")
        ]
    )

    if not experiment_folders:
    
            raise FileNotFoundError(
                "Experiment_ 폴더를 찾을 수 없습니다."
            )

    return experiment_folders


# ============================================================
# EXPERIMENT STATUS CHECK
# ============================================================

def classify_experiment_folders(experiment_folders):

    completed_experiments = []
    new_experiment_folders = []

    completion_marker = "Plaque_Assay_Analysis_Results.xlsx"

    for folder in experiment_folders:

        marker_file = folder / "Final_Results" / completion_marker

        if marker_file.exists():

            completed_experiments.append(folder)

        else:

            new_experiment_folders.append(folder)

    return (completed_experiments, new_experiment_folders)


# ============================================================
# FIND EXCEL FILE
# ============================================================
def find_excel_files(experiment_folder):

    excel_files = sorted(
        [
            file
            for file in experiment_folder.iterdir()
            if file.is_file()
            and file.suffix.lower() in [".xlsx", ".xls"]
            and not file.name.startswith("~$")
        ]
    )

    if not excel_files:

        raise ValueError(
            f"{experiment_folder.name}에 Excel 파일이 없습니다."
        )

    return excel_files



# ============================================================
# ANALYSIS FUNCTIONS
# ============================================================
# Load plaque assay data
def load_plaque_data(excel_file, target_sheet):
    
    # --------------------------------------------------------
    # Check available sheets
    # --------------------------------------------------------

    excel_reader = pd.ExcelFile(excel_file, engine="openpyxl")

    print("Sheet 목록:")

    for sheet in excel_reader.sheet_names:
        print(" -", sheet)

    # --------------------------------------------------------
    # Check target sheet
    # --------------------------------------------------------

    if target_sheet not in excel_reader.sheet_names:

        raise ValueError(
            f"'{target_sheet}' sheet를 찾을 수 없습니다."
        )

    # --------------------------------------------------------
    # Load data
    # --------------------------------------------------------

    df = pd.read_excel(
        excel_file,
        sheet_name=target_sheet,
        header=2,
        engine="openpyxl"
    )


    # --------------------------------------------------------
    # Check required columns
    # --------------------------------------------------------

    required_columns = [
        "Time point",
        "Sample",
        "plaque count",
        "Dilution factor"
    ]

    missing_columns = [
        column
        for column in required_columns
        if column not in df.columns
    ]

    if missing_columns:

        raise ValueError(
            "필수 column이 없습니다: "
            + ", ".join(missing_columns)
        )

    df = df[required_columns].copy()

    return df


# Process sample information
def process_sample_information(df, sample_order):

    df = df.copy()

    # --------------------------------------------------------
    # Clean sample names
    # --------------------------------------------------------

    df["Sample"] = df["Sample"].astype(str).str.strip()

    # --------------------------------------------------------
    # Extract replicate number
    # --------------------------------------------------------

    df["Replicate"] = df["Sample"].str.extract(r"-(\d+)$")[0]

    df["Replicate"] = pd.to_numeric(
        df["Replicate"],
        errors="coerce"
    )

    # --------------------------------------------------------
    # Remove replicate number from sample name
    # --------------------------------------------------------

    df["Name"] = (
        df["Sample"]
        .str.replace(r"-\d+$", "", regex=True)
        .str.strip()
    )

    # --------------------------------------------------------
    # Detect samples
    # --------------------------------------------------------

    detected_samples = (
        df["Name"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    # --------------------------------------------------------
    # Determine sample order
    # --------------------------------------------------------

    current_sample_order = [
        sample 
        for sample in sample_order
        if sample in detected_samples
    ]

    additional_samples = [
        sample
        for sample in detected_samples
        if sample not in current_sample_order
    ]

    current_sample_order.extend(sorted(additional_samples))

    return df, detected_samples, current_sample_order


# Calculate pfu
def calculate_pfu(df, plating_volume_ml):

    df = df.copy()

    # --------------------------------------------------------
    # Parse dilution factor
    # --------------------------------------------------------

    df["Dilution"] = df["Dilution factor"].apply(parse_dilution_factor)

    # --------------------------------------------------------
    # Convert plaque count to numeric
    # --------------------------------------------------------

    df["plaque count"] = pd.to_numeric(
        df["plaque count"],
        errors="coerce"
    )

    # --------------------------------------------------------
    # Calculate PFU/ml
    # --------------------------------------------------------

    df["PFU/ml"] = (
        df["plaque count"]
        /
        (df["Dilution"] * plating_volume_ml)
    )

    # --------------------------------------------------------
    # QC
    # --------------------------------------------------------

    invalid_pfu = (
        df["plaque count"].isna()
        | df["Dilution"].isna()
        | (df["Dilution"] <= 0)
    )

    df.loc[invalid_pfu,"PFU/ml"] = np.nan

    return df


# Pfu summary
def create_pfu_summary(df, sample_order):

    # --------------------------------------------------------
    # Create summary
    # --------------------------------------------------------

    summary_df = (
        df
        .groupby(
            ["Time point", "Name"],
            as_index=False
        )
        .agg(
            N=("PFU/ml", "count"),
            Mean=("PFU/ml", "mean"),
            SD=("PFU/ml", "std")
        )
    )

    # --------------------------------------------------------
    # Extract hpi
    # --------------------------------------------------------

    summary_df["hpi"] = extract_hpi(summary_df["Time point"])

    # --------------------------------------------------------
    # Apply sample order
    # --------------------------------------------------------

    summary_df["Name"] = pd.Categorical(
        summary_df["Name"],
        categories=sample_order,
        ordered=True
    )

    # --------------------------------------------------------
    # Sort summary
    # --------------------------------------------------------

    summary_df = (
        summary_df
        .sort_values(
            ["hpi", "Name"]
        )
        .reset_index(drop=True)
    )

    return summary_df


# Prepare plaque assay data for statistical analysis
def prepare_plaque_statistics_data(df, sample_order):

    # ------------------------------------------------------
    # Select samples used for statistical analysis
    # ------------------------------------------------------

    stats_df = df[df["Name"].isin(sample_order)].copy()

    # ------------------------------------------------------
    # Remove invalid rows
    # ------------------------------------------------------

    stats_df = stats_df.dropna(
        subset=[
            "Time point",
            "Name",
            "PFU/ml"
        ]
    ).copy()

    # ------------------------------------------------------
    # Remove zero or negative PFU/ml
    # ------------------------------------------------------

    stats_df = stats_df[stats_df["PFU/ml"] > 0].copy()

    # ------------------------------------------------------
    # Extract time point
    # ------------------------------------------------------

    stats_df["hpi"] = extract_hpi(stats_df["Time point"])

    # ------------------------------------------------------
    # Calculate log10 PFU/ml
    # ------------------------------------------------------

    stats_df["log10_PFU_ml"] = np.log10(stats_df["PFU/ml"])

    # ------------------------------------------------------
    # Apply sample order
    # ------------------------------------------------------

    stats_df["Name"] = pd.Categorical(
        stats_df["Name"],
        categories=sample_order,
        ordered=True
    )

    # ------------------------------------------------------
    # Sort data
    # ------------------------------------------------------

    stats_df = (
        stats_df
        .sort_values(
            ["hpi", "Name", "Replicate"]
        )
        .reset_index(drop=True)
    )

    return stats_df

# Two-Way ANOVA
def perform_plaque_two_way_anova(stats_df):

    # ------------------------------------------------------
    # Fit Two-way ANOVA model
    # ------------------------------------------------------

    anova_model = ols(
        "log10_PFU_ml ~ C(hpi) * C(Name)",
        data=stats_df
    ).fit()

    # ------------------------------------------------------
    # Perform Type II ANOVA
    # ------------------------------------------------------

    anova_table = sm.stats.anova_lm(
        anova_model,
        typ=2
    )

    # ------------------------------------------------------
    # Convert results to DataFrame
    # ------------------------------------------------------

    anova_results_df = (
        anova_table
        .reset_index()
        .rename(
            columns={
                "index": "Effect",
                "PR(>F)": "p_value"
            }
        )
    )

    return anova_results_df


# Dunnett test
def perform_plaque_dunnett_test(stats_df, sample_order, control):

    # ------------------------------------------------------
    # Initialize results
    # ------------------------------------------------------

    dunnett_results = []

    # ------------------------------------------------------
    # Detect time points
    # ------------------------------------------------------

    time_points = sorted(
        stats_df["hpi"]
        .dropna()
        .unique()
    )

    # ------------------------------------------------------
    # Perform Dunnett test at each time point
    # ------------------------------------------------------

    for time_point in time_points:

        time_df = stats_df[stats_df["hpi"] == time_point].copy()

        # ---------------------------------------------------
        # Control data
        # ---------------------------------------------------

        control_data = (
            time_df[
                time_df["Name"] == control
            ]["log10_PFU_ml"]
            .dropna()
            .to_numpy()
        )

        # ---------------------------------------------------
        # Treatment data
        # ---------------------------------------------------

        treatment_arrays = []
        treatment_names = []

        for sample_name in sample_order:

            if sample_name == control:
                continue

            treatment_data = (
                time_df[
                    time_df["Name"] == sample_name
                ]["log10_PFU_ml"]
                .dropna()
                .to_numpy()
            )

            if (
                len(control_data) >= 2
                and len(treatment_data) >= 2
            ):

                treatment_arrays.append(treatment_data)
                treatment_names.append(sample_name)

        # ---------------------------------------------------
        # Skip if no valid comparison exists
        # ---------------------------------------------------

        if len(treatment_arrays) == 0:
            continue

        # ---------------------------------------------------
        # Dunnett post-hoc test
        # ---------------------------------------------------

        result = dunnett(
            *treatment_arrays,
            control=control_data
        )

        # ---------------------------------------------------
        # Save results
        # ---------------------------------------------------

        for (
            sample_name,
            statistic,
            p_value
        ) in zip(
            treatment_names,
            np.atleast_1d(result.statistic),
            np.atleast_1d(result.pvalue)
        ):

            significance = significance_star(p_value)

            dunnett_results.append(
                {
                    "hpi": time_point,
                    "Control": control,
                    "Comparison": sample_name,
                    "Control_N": len(control_data),
                    "Comparison_N": len(
                        time_df[
                            time_df["Name"] == sample_name
                        ]["log10_PFU_ml"]
                        .dropna()
                    ),
                    "Statistic": statistic,
                    "p_value": p_value,
                    "Significance": significance
                }
            )

    # ------------------------------------------------------
    # Convert results to DataFrame
    # ------------------------------------------------------

    dunnett_results_df = pd.DataFrame(dunnett_results)

    return dunnett_results_df


# Growth curve
def create_growth_curve(
    experiment_summary,
    dunnett_results_df,
    sample_order,
    control,
    lod,
    output_path
):

    experiment_summary = experiment_summary.copy()

    # --------------------------------------------------------
    # Apply LOD
    # --------------------------------------------------------

    experiment_summary["Growth_PFU"] = experiment_summary["Mean"].copy()

    experiment_summary.loc[
        experiment_summary["Growth_PFU"] < lod,
        "Growth_PFU"
    ] = lod

    # --------------------------------------------------------
    # Create figure
    # --------------------------------------------------------

    plt.figure(figsize=(8, 6))

    sample_colors = {}

    # --------------------------------------------------------
    # Plot growth curves
    # --------------------------------------------------------

    for sample in sample_order:

        sample_df = (
            experiment_summary[
                experiment_summary["Name"] == sample
            ]
            .sort_values("hpi")
        )

        if sample_df.empty:
            continue

        line = plt.errorbar(
            sample_df["hpi"],
            sample_df["Growth_PFU"],
            yerr=sample_df["SD"],
            marker="o",
            markersize=6,
            linewidth=2,
            capsize=4,
            label=sample
        )

        sample_colors[sample] = line[0].get_color()

    # --------------------------------------------------------
    # Display LOD
    # --------------------------------------------------------

    if SHOW_LOD:
        plt.axhline(
            y=lod,
            linestyle="--",
            linewidth=1.5,
            label=f"LOD ({lod:g} PFU/ml)"
        )

    # --------------------------------------------------------
    # Display statistical significance
    # --------------------------------------------------------

    if (
        SHOW_SIGNIFICANCE
        and not dunnett_results_df.empty
    ):

        for hpi in sorted(
            dunnett_results_df["hpi"]
            .dropna()
            .unique()
        ):

            control_data = experiment_summary[
                (experiment_summary["Name"] == control)
                &
                (experiment_summary["hpi"] == hpi)
            ]

            if control_data.empty:
                continue

            control_y = (control_data["Growth_PFU"].iloc[0])

            stat_at_time = (
                dunnett_results_df[
                    dunnett_results_df["hpi"] == hpi
                ]
                .copy()
            )

            sample_rank = {
                sample: i
                for i, sample in enumerate(sample_order)
            }

            stat_at_time["Sample_Order"] = (
                stat_at_time["Comparison"]
                .map(sample_rank)
            )

            stat_at_time = (
                stat_at_time
                .sort_values("Sample_Order")
            )

            significant_results = (
                stat_at_time[
                    stat_at_time["Significance"]
                    .isin(
                        [
                            "*",
                            "**",
                            "***"
                        ]
                    )
                ]
                .copy()
            )

            if significant_results.empty:
                continue

            n_significant = len(significant_results)

            for star_level, (_, row) in enumerate(
                significant_results.iterrows()
            ):

                comparison = row["Comparison"]
                significance = row["Significance"]

                if comparison not in sample_colors:
                    continue

                star_color = sample_colors[comparison]

                level_from_bottom = n_significant - star_level
                

                star_y = (
                    control_y
                    * (1.5 ** level_from_bottom)
                )

                plt.text(
                    hpi,
                    star_y,
                    significance,
                    color=star_color,
                    ha="center",
                    va="bottom",
                    fontsize=14,
                    fontweight="bold",
                    zorder=10
                )

    # --------------------------------------------------------
    # Graph settings
    # --------------------------------------------------------

    plt.xlabel(XLABEL)
    plt.ylabel(YLABEL)
    plt.title(GRAPH_TITLE)

    plt.xticks(
        sorted(
            experiment_summary["hpi"]
            .dropna()
            .unique()
        )
    )

    plt.yscale("log")

    # --------------------------------------------------------
    # Automatic Y-axis range
    # --------------------------------------------------------

    y_max = experiment_summary["Growth_PFU"].max()

    max_star_y = y_max

    if (
        SHOW_SIGNIFICANCE
        and not dunnett_results_df.empty
    ):

        significant_counts = (
            dunnett_results_df[
                dunnett_results_df["Significance"]
                .isin(
                    [
                        "*",
                        "**",
                        "***"
                    ]
                )
            ]
            .groupby("hpi")
            .size()
        )

        for hpi, n_stars in significant_counts.items():

            control_data = experiment_summary[
                (experiment_summary["Name"] == control)
                &
                (experiment_summary["hpi"] == hpi)
            ]

            if control_data.empty:
                continue

            control_y = control_data["Growth_PFU"].iloc[0]

            star_y = (
                control_y
                * (1.5 ** n_stars)
            )

            max_star_y = max(
                max_star_y,
                star_y
            )

    y_top = max_star_y * 2.5

    plt.ylim(
        bottom=lod * 0.8,
        top=y_top
    )

    if SHOW_LEGEND:

        plt.legend()

    plt.tight_layout()

    # --------------------------------------------------------
    # Save figure
    # --------------------------------------------------------

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()
    plt.close()

    return experiment_summary


# Save analysis results

def save_analysis_results(
    output_path,
    df,
    experiment_summary,
    anova_results_df,
    dunnett_results_df
):

    with pd.ExcelWriter(
        output_path,
        engine="openpyxl"
    ) as writer:

        # --------------------------------------------------------
        # Raw data
        # --------------------------------------------------------

        df.to_excel(
            writer,
            sheet_name="Raw_Data",
            index=False
        )

        # --------------------------------------------------------
        # PFU summary
        # --------------------------------------------------------

        experiment_summary.to_excel(
            writer,
            sheet_name="PFU_Summary",
            index=False
        )

        # --------------------------------------------------------
        # Two-way ANOVA
        # --------------------------------------------------------

        anova_results_df.to_excel(
            writer,
            sheet_name="Two_Way_Anova",
            index=False
        )

        # --------------------------------------------------------
        # Dunnett post-hoc test
        # --------------------------------------------------------

        dunnett_results_df.to_excel(
            writer,
            sheet_name="Dunnett_Results",
            index=False
        )


# ============================================================
# STEP 3. EXPERIMENT ANALYSIS
# ============================================================

def analyze_experiment(experiment_folder):

    print("\n" + "=" * 60)
    print(experiment_folder.name)
    print("=" * 60)

    # --------------------------------------------------------
    # Result paths
    # --------------------------------------------------------

    result_folder = experiment_folder / "Final_Results"

    growth_curve_path = result_folder / "Growth_Curve.png"

    result_excel_path = result_folder / "Plaque_Assay_Analysis_Results.xlsx"


    # --------------------------------------------------------
    # Create result folder
    # --------------------------------------------------------

    result_folder.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # INPUT FILE SEARCH
    # --------------------------------------------------------

    excel_files = find_excel_files(experiment_folder)


    # --------------------------------------------------------
    # Analyze Excel files
    # --------------------------------------------------------

    for excel_file in excel_files:

        print("Excel:", excel_file.name)

        # ----------------------------------------------------
        # STEP 3-1. DATA LOADING 
        # ----------------------------------------------------

        df = load_plaque_data(excel_file, TARGET_SHEET)


        # ----------------------------------------------------
        # STEP 3-2. SAMPLE INFORMATION PROCESSING
        # ----------------------------------------------------

        df, detected_samples, sample_order = process_sample_information(
            df,
            sample_order=SAMPLE_ORDER
        )

        # ----------------------------------------------------
        # Check control sample
        # ----------------------------------------------------

        if CONTROL not in detected_samples:

            print(
                f"WARNING: Control "
                f"'{CONTROL}'이 데이터에 없습니다."
            )

            print(
                "Dunnett 분석은 수행되지 않을 수 있습니다."
            )

        print("자동 인식된 Sample:")

        for sample in sample_order:
            print(" -", sample)

        # ----------------------------------------------------
        # STEP 3-3. PFU CALCULATION
        # ----------------------------------------------------

        df = calculate_pfu(df, PLATING_VOLUME_ML)

        # ----------------------------------------------------
        # STEP 3-4. FINAL SUMMARY
        # ----------------------------------------------------

        experiment_summary = (
            create_pfu_summary(
                df,
                sample_order
            )
        )


        # ----------------------------------------------------
        # STEP 3-5. STATISTICAL DATA PREPARATION
        # ----------------------------------------------------

        stats_df = prepare_plaque_statistics_data(
            df,
            sample_order
        )

        # ----------------------------------------------------
        # STEP 3-6. TWO-WAY ANOVA
        # ----------------------------------------------------

        print("\n" + "=" * 60)
        print("TWO-WAY ANOVA")
        print("=" * 60)


        anova_results_df = perform_plaque_two_way_anova(
            stats_df
        )

        print("\n" + "=" * 60)
        print("Two-Way ANOVA Result")
        print()
        print(anova_results_df)

        # ----------------------------------------------------
        # STEP 3-7. DUNNETT POST-HOC TEST
        # ----------------------------------------------------

        print("\n" + "=" * 60)
        print("DUNNETT POST-HOC")
        print("=" * 60)

        dunnett_results_df = (
            perform_plaque_dunnett_test(
                stats_df,
                sample_order,
                CONTROL
            )
        )

        # ----------------------------------------------------
        # Display Dunnett results
        # ----------------------------------------------------

        print("\n" + "=" * 60)
        print("Dunnett Post-hoc Result")
        print("=" * 60)

        if not dunnett_results_df.empty:

            print(
                dunnett_results_df[
                    [
                        "hpi",
                        "Control",
                        "Comparison",
                        "Control_N",
                        "Comparison_N",
                        "p_value",
                        "Significance"
                    ]
                ].to_string(
                    index=False
                )
            )

        else:

            print("Dunnett 결과 없음")

        # ----------------------------------------------------
        # STEP 3-8. GROWTH CURVE
        # ----------------------------------------------------

        experiment_summary = (
            create_growth_curve(
                experiment_summary,
                dunnett_results_df,
                sample_order,
                CONTROL,
                PLAQUE_ASSAY_LOD,
                growth_curve_path
            )
        )

        # ----------------------------------------------------
        # STEP 3-9. RESULT SAVE
        # ----------------------------------------------------

        save_analysis_results(
            result_excel_path,
            df,
            experiment_summary,
            anova_results_df,
            dunnett_results_df
        )

        # ----------------------------------------------------
        # Analysis complete
        # ----------------------------------------------------


        print("\n" + "=" * 60)
        print("결과 저장")
        print("=" * 60)

        print()
        print("Growth curve saved:")
        print(growth_curve_path)

        print()
        print("Analysis results saved:")
        print(result_excel_path)

        print("\n" + "=" * 60)
        print(f"{experiment_folder.name} 분석 완료")
        print("=" * 60)
        print()


# ============================================================
# MAIN 
# ============================================================

def main():

    # ========================================================
    # STEP 1. EXPERIMENT FOLDER SEARCH
    # ========================================================

    print("\n" + "=" * 60)
    print("STEP 1. EXPERIMENT FOLDER SEARCH")
    print("=" * 60)

    experiment_folders = find_experiment_folders(BASE_FOLDER)

    print(
        f"\n찾은 Experiment 수: {len(experiment_folders)}"
    )

    for folder in experiment_folders:
        
        print(" -", folder.name)

    # ========================================================
    # STEP 2. EXPERIMENT STATUS CHECK
    # ========================================================

    print("\n" + "=" * 60)
    print("STEP 2. EXPERIMENT STATUS CHECK")
    print("=" * 60)

    completed_experiments, new_experiment_folders = (
        classify_experiment_folders(experiment_folders)
    )

    print(
        f"\n완료된 Experiment : {len(completed_experiments)}"
    )

    for folder in completed_experiments:

        print(f"  [SKIP] {folder.name}")

    print(
        f"\n새로운 Experiment : {len(new_experiment_folders)}"
    )

    for folder in new_experiment_folders:

        print(f"  [NEW] {folder.name}")


    # ========================================================
    # NO NEW EXPERIMENT CHECK
    # ========================================================

    if not new_experiment_folders:

        print("\n새로운 Experiment이 없습니다.")

        return


    # ========================================================
    # MAIN EXPERIMENT LOOP
    # ========================================================

    for experiment_folder in new_experiment_folders:

        try:

            analyze_experiment(experiment_folder)

        except Exception as error:

            print()
            print("!" * 60)

            print(
                f"{experiment_folder.name} "
                f"분석 중 오류 발생"
            )

            print(f"Error: {error}")

            print(
                "이 Experiment는 건너뛰고 다음 Experiment를 분석합니다."
            )

            print("!" * 60)
            print()

    # ============================================================
    # ALL NEW EXPERIMENT COMPLETE
    # ============================================================
    
    print("\n\n" + "=" * 60)
    print("ALL NEW EXPERIMENT ANALYSIS FINISHED")
    print("=" * 60)


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    main()


