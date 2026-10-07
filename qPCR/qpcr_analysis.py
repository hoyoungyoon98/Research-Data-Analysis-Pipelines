# ============================================================
# qPCR AUTOMATED ANALYSIS PIPELINE
# ============================================================

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import dunnett
import statsmodels.api as sm
from statsmodels.formula.api import ols
from matplotlib.ticker import LogFormatterMathtext


# ============================================================
# STEP 0. USER CONFIGURATION
# ============================================================

# ANALYSIS TYPE
# "relative" / "absolute"
ANALYSIS_TYPE = "relative"

# Project directory
PROJECT_DIR = Path(__file__).resolve().parent

# Example data directory
BASE_FOLDER = PROJECT_DIR / "example_data" / ANALYSIS_TYPE


# ============================================================
# COMMON CONFIGURATION
# ============================================================
# Statistics
CONTROL = "Control"

# Samples excluded from analysis
EXCLUDED_SAMPLES = [
    "Cell only",
    "EV"
]

# Graph
GRAPH_TITLE = "Gene Expression"

XLABEL = "X-label"

YLABEL = "Y-label"

SHOW_SIGNIFICANCE = True


# ============================================================
# RELATIVE qPCR CONFIGURATION
# ============================================================
# Sample order
RELATIVE_SAMPLE_ORDER = [
    "Control",
    "Sample_A",
    "Sample_B",
    "Sample_C",
    "Sample_D"
]

RELATIVE_TARGET_SHEET = "Target"

RELATIVE_REFERENCE_SHEET = "Housekeeping"


# ============================================================
# ABSOLUTE qPCR CONFIGURATION
# ============================================================
# Sample order
ABSOLUTE_SAMPLE_ORDER = [
    "Control",
    "Sample_A",
    "Sample_B"
]

ABSOLUTE_TARGET_SHEET = "Experiment"

TEMPLATE_VOLUME_UL = 4.0

# Standard curve QC
MIN_R_SQUARED = 0.99

MIN_EFFICIENCY = 90.0
MAX_EFFICIENCY = 110.0

MIN_STANDARD_POINTS = 5

# Allowed deviation between observed and predicted Cq
MAX_CQ_RESIDUAL = 1.0

# Allowed deviation of ΔCq between adjacent standards
MAX_DELTA_CQ_DEVIATION = 1.0


# ============================================================
# BASIC VALIDATION
# ============================================================

if ANALYSIS_TYPE not in ["relative", "absolute"]:

    raise ValueError(
        "ANALYSIS_TYPE은 'relative' 또는 'absolute'이어야 합니다."
    )


# ============================================================
# COMMON FUNCTIONS
# ============================================================

# Significance
def significance_star(p_value):

    if pd.isna(p_value):
        return "NA"

    elif p_value < 0.001:
        return "***"

    elif p_value < 0.01:
        return "**"

    elif p_value < 0.05:
        return "*"

    else:
        return "ns"


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

    if ANALYSIS_TYPE == "relative":

        completion_marker = "Relative_qPCR_Analysis_Results.xlsx"

    elif ANALYSIS_TYPE == "absolute":

        completion_marker = "Absolute_qPCR_Analysis_Results.xlsx"

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
# RELATIVE qPCR FUNCTIONS
# ============================================================
# Load relative qPCR data 
def load_relative_qpcr_data(
    excel_files, 
    target_sheet_name, 
    reference_sheet_name
):

    required_sheets = [target_sheet_name, reference_sheet_name]

    target_dataframes = []
    reference_dataframes = []

    for replicate_number, excel_file in enumerate(excel_files, start=1):

        # ----------------------------------------------------
        # Check required sheets
        # ----------------------------------------------------

        excel_data = pd.ExcelFile(excel_file)

        missing_sheets = [
            sheet
            for sheet in required_sheets
            if sheet not in excel_data.sheet_names
        ]

        if missing_sheets:

            raise ValueError(
                f"{excel_file.name}에 필요한 sheet가 없습니다:"
                + ", ".join(missing_sheets)
            )

        # ----------------------------------------------------
        # Load Target / Reference data
        # ----------------------------------------------------

        print(
            f"\nReading Replicate {replicate_number}: {excel_file.name}"
        )

        target_df = pd.read_excel(
            excel_file,
            sheet_name=target_sheet_name
        )

        reference_df = pd.read_excel(
            excel_file,
            sheet_name=reference_sheet_name
        )

        target_df["Biological_Replicate"] = replicate_number
        reference_df["Biological_Replicate"] = replicate_number

        target_dataframes.append(target_df)
        reference_dataframes.append(reference_df)

    # -------------------------------------------------------
    # Combine biological replicates
    # -------------------------------------------------------

    target_all = pd.concat(target_dataframes, ignore_index=True)

    reference_all = pd.concat(reference_dataframes, ignore_index=True)


    return target_all, reference_all


# Prepare relative qPCR data
def prepare_relative_qpcr_data(
    target_all,
    reference_all,
    experiment_name
):

    # ------------------------------------------------------
    # Add experiment name 
    # ------------------------------------------------------

    target_all = target_all.copy()
    reference_all = reference_all.copy()

    target_all["Experiment"] = experiment_name
    reference_all["Experiment"] = experiment_name


    # ------------------------------------------------------
    # Fill sample names
    # ------------------------------------------------------

    target_all["name"] = (
        target_all["name"]
        .ffill()
        .astype(str)
        .str.strip()
    )

    reference_all["name"] = (
        reference_all["name"]
        .ffill()
        .astype(str)
        .str.strip()
    )


    # ------------------------------------------------------
    # Convert Cq to numeric
    # ------------------------------------------------------

    target_all["Cq"] = pd.to_numeric(
        target_all["Cq"],
        errors="coerce"
    )

    reference_all["Cq"] = pd.to_numeric(
        reference_all["Cq"],
        errors="coerce"
    )


    # ------------------------------------------------------
    # Keep required columns
    # ------------------------------------------------------

    target_df = target_all[
        [
            "Experiment",
            "Biological_Replicate",
            "name",
            "Cq"
        ]
    ].copy()

    reference_df = reference_all[
        [
            "Experiment",
            "Biological_Replicate",
            "name",
            "Cq"
        ]
    ].copy()

    target_df = target_df.rename(
        columns={"Cq": "Target_Cq"}
    )

    reference_df = reference_df.rename(
        columns={"Cq": "Reference_Cq"}
    )


    # ------------------------------------------------------
    # Check row counts 
    # ------------------------------------------------------

    if len(target_df) != len(reference_df):

        raise ValueError(
            "Target와 Reference의 전체 row 수가 다릅니다."
        )


    # ------------------------------------------------------
    # Create row numbers within biological replicates
    # ------------------------------------------------------

    target_df["row_number"] = (
        target_df
        .groupby(
            [
                "Experiment",
                "Biological_Replicate"
            ]
        )
        .cumcount()
    )

    reference_df["row_number"] = (
        reference_df
        .groupby(
            [
                "Experiment",
                "Biological_Replicate"
            ]
        )
        .cumcount()
    )


    # ------------------------------------------------------
    # Match Target / Reference 
    # ------------------------------------------------------

    merged = pd.merge(
        target_df,
        reference_df,
        on=[
            "Experiment",
            "Biological_Replicate",
            "row_number"
        ],
        suffixes=("_Target", "_Reference"),
        how="inner"
    )


    # ------------------------------------------------------
    # Check sample matching 
    # ------------------------------------------------------

    name_match = merged["name_Target"] == merged["name_Reference"]
    

    if not name_match.all():

        mismatch_df = merged.loc[
            ~name_match,
            [
                "Biological_Replicate",
                "row_number",
                "name_Target",
                "name_Reference"
            ]
        ]

        print("\nTarget / Reference mismatch:")

        print(mismatch_df.to_string(index=False))

        raise ValueError(
            "Target와 Reference의 Sample 순서가 일치하지 않습니다."
        )

    merged["name"] = merged["name_Target"]


    return merged



# Calculate relative expression
def calculate_relative_expression(
        merged,
        control,
        excluded_samples,
        sample_order
):

    df = merged.copy()

    # ---------------------------------------------------
    # Calculate ΔCt
    # ---------------------------------------------------

    df["Delta_Ct"] = df["Target_Cq"] - df["Reference_Cq"]


    # ---------------------------------------------------
    # Calculate technical replicate mean
    # ---------------------------------------------------

    technical_mean_df = (
        df
        .groupby(
            [
                "Experiment",
                "Biological_Replicate",
                "name"
            ],
            as_index=False
        )
        .agg(
            Mean_Delta_Ct=("Delta_Ct", "mean"),
            Technical_N=("Delta_Ct", "count")
        )
    )


    # ---------------------------------------------------
    # Detect samples
    # ---------------------------------------------------

    detected_samples = (
        technical_mean_df["name"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )


    # ---------------------------------------------------
    # Exclude samples not used for analysis
    # ---------------------------------------------------

    analysis_samples = [
        sample
        for sample in detected_samples
        if sample not in excluded_samples
    ]


    # ---------------------------------------------------
    # Determine sample order 
    # ---------------------------------------------------

    current_sample_order = [
        sample 
        for sample in sample_order
        if sample in analysis_samples
    ]

    additional_samples = [
        sample
        for sample in analysis_samples
        if sample not in current_sample_order
    ]

    current_sample_order.extend(sorted(additional_samples))


    # ---------------------------------------------------
    # Check control
    # ---------------------------------------------------

    if control not in analysis_samples:

        raise ValueError(
            f"CONTROL '{control}'을 sample 데이터에서 찾을 수 없습니다."
        )


    # ---------------------------------------------------
    # Apply sample order 
    # ---------------------------------------------------

    technical_mean_df["name"] = pd.Categorical(
        technical_mean_df["name"],
        categories=(excluded_samples + current_sample_order),
        ordered=True
    )

    technical_mean_df = (
        technical_mean_df
        .sort_values(
            [
                "Biological_Replicate",
                "name"
            ]
        )
        .reset_index(drop=True)
    )


    # ---------------------------------------------------
    # Calculate control ΔCt
    # ---------------------------------------------------

    control_df = (
        technical_mean_df[
            technical_mean_df["name"] == control
        ][
            [
                "Experiment",
                "Biological_Replicate",
                "Mean_Delta_Ct"
            ]
        ]
        .rename(
            columns={
                "Mean_Delta_Ct": "Control_Mean_Delta_Ct"
            }
        )
    )


    # ---------------------------------------------------
    # Calculate ΔΔCt
    # ---------------------------------------------------

    technical_mean_df = pd.merge(
        technical_mean_df,
        control_df,
        on=[
            "Experiment",
            "Biological_Replicate"
        ],
        how="left"
    )

    technical_mean_df["DeltaDelta_Ct"] = (
        technical_mean_df["Mean_Delta_Ct"]
        - technical_mean_df["Control_Mean_Delta_Ct"]
    )

    # ---------------------------------------------------
    # Calculate relative expression
    # ---------------------------------------------------

    technical_mean_df["Relative_Expression"] = (
        2 ** (-technical_mean_df["DeltaDelta_Ct"])
    )


    return (
        df,
        technical_mean_df,
        detected_samples,
        current_sample_order
    )


# Relative qPCR Statistics
def analyze_relative_statistics(
        technical_mean_df,
        sample_order,
        control
):

    # ---------------------------------------------------
    # Prepare statistical analysis data
    # ---------------------------------------------------

    stats_df = (
        technical_mean_df[
            technical_mean_df["name"].isin(sample_order)
        ].copy()
    )

    stats_df = stats_df.dropna(
        subset=[
            "Biological_Replicate",
            "name",
            "Mean_Delta_Ct",
            "Relative_Expression"
        ]
    ).copy()

    stats_df = (
        stats_df
        .sort_values(
            [
                "name",
                "Biological_Replicate"
            ]
        )
        .reset_index(drop=True)
    )


    # ---------------------------------------------------
    # Paired t-test
    # ---------------------------------------------------

    stats_results = []

    comparison_samples = [
        sample
        for sample in sample_order
        if sample != control
    ]

    for sample_name in comparison_samples:

        sample_data = (
            stats_df[
                stats_df["name"] == sample_name
            ][
                [
                    "Biological_Replicate",
                    "Mean_Delta_Ct"
                ]
            ].copy()
        )

        control_data = (
            stats_df[
                stats_df["name"] == control
            ][
                [
                    "Biological_Replicate",
                    "Mean_Delta_Ct"
                ]
            ].copy()
        )

        # Pair by biological replicate
        comparison_df = pd.merge(
            sample_data,
            control_data,
            on="Biological_Replicate",
            suffixes=("_Sample", "_Control"),
            how="inner"
        )

        # Paired t-test
        n = len(comparison_df)

        if n >= 2:

            t_statistic, p_value = stats.ttest_rel(
                comparison_df["Mean_Delta_Ct_Sample"],
                comparison_df["Mean_Delta_Ct_Control"]
            )

        else:

            t_statistic = np.nan
            p_value = np.nan
            significance = "NA"

        significance = significance_star(p_value)

        stats_results.append(
            {
                "Control": control,
                "Comparison": sample_name,
                "N": n,
                "t_statistic": t_statistic,
                "p_value": p_value,
                "Significance": significance
            }
        )

    stats_results_df = pd.DataFrame(stats_results)

    return stats_df, stats_results_df


# Create relative qPCR summary
def create_relative_summary(
    technical_mean_df,
    stats_results_df,
    sample_order,
    control
):

    # ---------------------------------------------------
    # Calculate Mean, SD, N
    # ---------------------------------------------------

    summary_df = (
        technical_mean_df[
            technical_mean_df["name"].isin(sample_order)
        ]
        .groupby("name", observed=True)
        .agg(
            Mean_Relative_Expression=("Relative_Expression", "mean"),
            SD_Relative_Expression=("Relative_Expression", "std"),
            N=("Relative_Expression", "count")
        )
        .reset_index()
    )


    # ---------------------------------------------------
    # Apply sample order 
    # ---------------------------------------------------

    summary_df["name"] = pd.Categorical(
        summary_df["name"],
        categories=sample_order,
        ordered=True
    )

    summary_df = (
        summary_df
        .sort_values("name")
        .reset_index(drop=True)
    )

    # ---------------------------------------------------
    # Merge statistical results
    # ---------------------------------------------------

    summary_df = pd.merge(
        summary_df,
        stats_results_df[
            [
                "Comparison",
                "p_value",
                "Significance"
            ]
        ],
        left_on="name",
        right_on="Comparison",
        how="left"
    )

    summary_df = summary_df.drop(columns=["Comparison"])


    # --------------------------------------------------
    # Control 
    # --------------------------------------------------

    summary_df.loc[
        summary_df["name"] == control,
        "Significance"
    ] = "-"

    summary_df.loc[
        summary_df["name"] == control,
        "p_value"
    ] = np.nan

    summary_df["Significance"] = summary_df["Significance"].fillna("NA")


    return summary_df


# Create relative qPCR graph
def create_relative_qpcr_graph(
    summary_df,
    result_folder,
    graph_title,
    xlabel,
    ylabel,
    show_significance
):

    # --------------------------------------------------
    # Figure size 
    # --------------------------------------------------

    figure_width = max(8, len(summary_df) * 0.8)

    fig, ax = plt.subplots(
        figsize=(figure_width, 6)
    )


    # --------------------------------------------------
    # Bar graph
    # --------------------------------------------------

    bars = ax.bar(
        summary_df["name"],
        summary_df["Mean_Relative_Expression"],
        yerr=summary_df["SD_Relative_Expression"],
        capsize=5
    )


    # --------------------------------------------------
    # Axis labels
    # --------------------------------------------------

    ax.set_title(graph_title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        axis="both",
        direction="out"
    )


    # --------------------------------------------------
    # Determine data range
    # --------------------------------------------------

    y_max = (
        summary_df["Mean_Relative_Expression"]
        + summary_df["SD_Relative_Expression"].fillna(0)
    ).max()

    max_annotation_y = y_max


    # --------------------------------------------------
    # Add significance
    # --------------------------------------------------

    if show_significance:

        for (
            bar,
            mean,
            sd,
            significance
        ) in zip(
            bars,
            summary_df["Mean_Relative_Expression"],
            summary_df["SD_Relative_Expression"].fillna(0),
            summary_df["Significance"]
        ):

            if significance in ["-", "NA", "ns"]:
                continue

            y_position = mean + sd

            ax.text(
                bar.get_x() + bar.get_width() / 2,
                y_position * 1.03,
                significance,
                ha="center",
                va="bottom",
                fontsize=12
            )

            max_annotation_y = max(
                max_annotation_y,
                y_position * 1.12
            )

    # --------------------------------------------------
    # Automatic Y-axis range
    # --------------------------------------------------

    final_y_max = max(
        y_max,
        max_annotation_y
    )

    ax.set_ylim(
        0,
        final_y_max * 1.15
    )


    # --------------------------------------------------
    # Save graph
    # --------------------------------------------------

    fig.tight_layout()

    result_folder.mkdir(
        parents=True,
        exist_ok=True
    )

    graph_path = (
        result_folder
        / "Relative_qPCR_Graph.png"
    )

    fig.savefig(
        graph_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()

    plt.close(fig)


    return graph_path


# Save relative qPCR analysis results 
def save_relative_qpcr_results(
    result_file,
    merged,
    technical_mean_df,
    stats_df,
    stats_results_df,
    summary_df
):

    # --------------------------------------------------
    # Create result folder 
    # --------------------------------------------------

    result_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------
    # Save all results to one excel file
    # --------------------------------------------------

    with pd.ExcelWriter(
        result_file,
        engine="openpyxl"
    ) as writer:

        # Matched Target / Reference data
        merged.to_excel(
            writer,
            sheet_name="Matched_Data",
            index=False
        )

        # Technical replicate mean
        technical_mean_df.to_excel(
            writer,
            sheet_name="Technical_Mean",
            index=False
        )

        # Statistical analysis data
        stats_df.to_excel(
            writer,
            sheet_name="Statistical_Data",
            index=False
        )

        # Paired t-test results
        stats_results_df.to_excel(
            writer,
            sheet_name="Statistics",
            index=False
        )

        # Final summary
        summary_df.to_excel(
            writer,
            sheet_name="Summary",
            index=False
        )


    return result_file


# ============================================================
# ABSOLUTE qPCR FUNCTIONS
# ============================================================
# Load absolute qPCR data
def load_absolute_qpcr_data(
    excel_file,
    target_sheet_name
):

    # ------------------------------------------------------
    # Open Excel file
    # ------------------------------------------------------

    excel_reader = pd.ExcelFile(
        excel_file,
        engine="openpyxl"
    )


    # ------------------------------------------------------
    # Check sheet names
    # ------------------------------------------------------

    print("\nSheet 목록:")

    for sheet_name in excel_reader.sheet_names:
        print(" -", sheet_name)


    # ------------------------------------------------------
    # Check required sheet
    # ------------------------------------------------------

    if target_sheet_name not in excel_reader.sheet_names:

        raise ValueError(
            f"'{target_sheet_name}' sheet를 찾을 수 없습니다."
        )


    # ------------------------------------------------------
    # Load data
    # ------------------------------------------------------

    df = pd.read_excel(
        excel_file,
        sheet_name=target_sheet_name,
        engine="openpyxl"
    )


    # ------------------------------------------------------
    # Check required columns
    # ------------------------------------------------------

    required_columns = [
        "Well",
        "Content",
        "Sample",
        "Cq",
        "Starting Quantity (SQ)"
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


    return df


# Prepare absolute qPCR sample data
def prepare_absolute_sample_data(
    df,
    sample_order,
    control,
    excluded_samples
):

    # ------------------------------------------------------
    # Extract unknown samples 
    # ------------------------------------------------------

    sample_df = df[
        df["Content"]
        .astype(str)
        .str.strip()
        .str.lower()
        == "unkn"
    ].copy()


    # ------------------------------------------------------
    # Clean sample names 
    # ------------------------------------------------------

    sample_df["Sample"] = sample_df["Sample"].astype(str).str.strip()


    # ------------------------------------------------------
    # Exclude MOCK samples
    # ------------------------------------------------------

    sample_df = sample_df[
        ~sample_df["Sample"].str.contains(
            r"(?i)\bMOCK\b",
            regex=True,
            na=False
        )
    ].copy()


    # ------------------------------------------------------
    # Extract biological replicate
    # ------------------------------------------------------

    sample_df["Biological_Replicate"] = (
        sample_df["Sample"]
        .str.extract(r"-(\d+)\s*$")[0]
    )

    sample_df["Biological_Replicate"] = pd.to_numeric(
        sample_df["Biological_Replicate"],
        errors="coerce"
    )


    # ------------------------------------------------------
    # Extract time point 
    # ------------------------------------------------------

    sample_df["Time_h"] = (
        sample_df["Sample"]
        .str.extract(r"(?i)(\d+)\s*h")[0]
    )

    sample_df["Time_h"] = pd.to_numeric(
        sample_df["Time_h"],
        errors="coerce"
    )


    # ------------------------------------------------------
    # Extract sample name 
    # ------------------------------------------------------

    sample_df["Name"] = (
        sample_df["Sample"]
        .str.extract(r"(?i)\d+\s*h\s+(.+?)\s*-\d+\s*$")[0]
    )

    sample_df["Name"] = (
        sample_df["Name"]
        .astype("string")
        .str.strip()
    )


    # ------------------------------------------------------
    # Convert Cq to numeric
    # ------------------------------------------------------

    sample_df["Cq"] = pd.to_numeric(
        sample_df["Cq"],
        errors="coerce"
    )


    # ------------------------------------------------------
    # Detect invalid rows
    # ------------------------------------------------------

    invalid_sample_rows = sample_df[
        sample_df[
            [
                "Time_h",
                "Name",
                "Biological_Replicate",
                "Cq"
            ]
        ].isna().any(axis=1)
    ]

    if not invalid_sample_rows.empty:

        print(
            "\n경고: Sample 정보 해석에 실패한 row가 있습니다."
        )

        print(
            invalid_sample_rows[
                [
                    "Well",
                    "Sample",
                    "Cq"
                ]
            ].to_string(index=False)
        )


    # ------------------------------------------------------
    # Remove invalid rows
    # ------------------------------------------------------

    sample_df = sample_df.dropna(
        subset=[
            "Time_h",
            "Name",
            "Biological_Replicate",
            "Cq"
        ]
    ).copy()


    sample_df = sample_df[
        [
            "Well",
            "Sample",
            "Time_h",
            "Name",
            "Biological_Replicate",
            "Cq"
        ]
    ].copy()


    # ------------------------------------------------------
    # Detect samples 
    # ------------------------------------------------------

    detected_samples = (
        sample_df["Name"]
        .dropna()
        .astype(str)
        .unique()
        .tolist()
    )

    # ------------------------------------------------------
    # Exclude samples not used for analysis
    # ------------------------------------------------------

    analysis_samples = [
        sample
        for sample in detected_samples
        if sample not in excluded_samples
    ]


    # ------------------------------------------------------
    # Determine sample order 
    # ------------------------------------------------------

    current_sample_order = [
        sample
        for sample in sample_order
        if sample in analysis_samples
    ]

    additional_samples = [
        sample
        for sample in analysis_samples
        if sample not in current_sample_order
    ]

    current_sample_order.extend(sorted(additional_samples))


    # ------------------------------------------------------
    # Check control
    # ------------------------------------------------------

    if control not in analysis_samples:

        raise ValueError(
            f"CONTROL '{control}'을 sample 데이터에서 찾을 수 없습니다."
        )


    return (
        sample_df,
        detected_samples,
        current_sample_order
    )


# Prepare standard curve data
def prepare_standard_curve_data(
    df,
    min_standard_points
):

    # ------------------------------------------------------
    # Extract standard data
    # ------------------------------------------------------

    standard_df = df[
        df["Content"]
        .astype(str)
        .str.strip()
        .str.lower()
        == "std"
    ].copy()


    # ------------------------------------------------------
    # Keep required columns 
    # ------------------------------------------------------

    standard_df = standard_df[
        [
            "Well",
            "Cq",
            "Starting Quantity (SQ)"
        ]
    ].copy()


    # ------------------------------------------------------
    # Convert values to numeric 
    # ------------------------------------------------------

    standard_df["Cq"] = pd.to_numeric(
        standard_df["Cq"],
        errors="coerce"
    )

    standard_df["Starting Quantity (SQ)"] = pd.to_numeric(
        standard_df["Starting Quantity (SQ)"],
        errors="coerce"
    )


    # ------------------------------------------------------
    # Remove invalid rows
    # ------------------------------------------------------

    standard_df = standard_df.dropna(
        subset=[
            "Cq",
            "Starting Quantity (SQ)"
        ]
    ).copy()


    # ------------------------------------------------------
    # Remove zero or negative quantity
    # ------------------------------------------------------

    standard_df = standard_df[
        standard_df["Starting Quantity (SQ)"] > 0
    ].copy()


    # ------------------------------------------------------
    # Check number of standard points
    # ------------------------------------------------------

    if len(standard_df) < min_standard_points:

        raise ValueError(
            f"Standard point 수가 "
            f"{min_standard_points}개보다 적습니다."
        )

    # ------------------------------------------------------
    # Calculate log10 Quantity
    # ------------------------------------------------------

    standard_df["log10_Quantity"] = np.log10(
        standard_df["Starting Quantity (SQ)"]
    )


    return standard_df


# Perform automatic standard curve QC
def perform_standard_curve_qc(
    standard_df,
    min_standard_points,
    min_r_squared,
    min_efficiency,
    max_efficiency,
    max_cq_residual,
    max_delta_cq_deviation
):

    # ------------------------------------------------------
    # Prepare standard data for QC 
    # ------------------------------------------------------

    standard_qc_df = (
        standard_df
        .sort_values("Starting Quantity (SQ)", ascending=False)
        .reset_index(drop=True)
        .copy()
    )


    # ------------------------------------------------------
    # Initialize QC results
    # ------------------------------------------------------

    qc_history = []

    final_standard_df = None

    final_slope = None
    final_intercept = None
    final_r_squared = None
    final_efficiency = None


    # ------------------------------------------------------
    # Automatic Standard Curve QC
    # ------------------------------------------------------

    while len(standard_qc_df) >= min_standard_points:

        # --------------------------------------------------
        # Linear regression
        # --------------------------------------------------

        slope, intercept = np.polyfit(
            standard_qc_df["log10_Quantity"],
            standard_qc_df["Cq"],
            1
        )


        # --------------------------------------------------
        # Predicted Cq
        # --------------------------------------------------

        standard_qc_df["Predicted_Cq"] = (
            slope
            * standard_qc_df["log10_Quantity"]
            + intercept
        )


        # --------------------------------------------------
        # Residual
        # --------------------------------------------------

        standard_qc_df["Residual_Cq"] = (
            standard_qc_df["Cq"]
            - standard_qc_df["Predicted_Cq"]
        )


        # --------------------------------------------------
        # R-squared 
        # --------------------------------------------------

        ss_res = np.sum(standard_qc_df["Residual_Cq"] ** 2)

        ss_tot = np.sum(
            (standard_qc_df["Cq"] - standard_qc_df["Cq"].mean()) ** 2
        )

        r_squared = 1 - (ss_res / ss_tot)


        # --------------------------------------------------
        # Efficiency
        # --------------------------------------------------

        efficiency = (10 ** (-1 / slope) - 1) * 100


        # --------------------------------------------------
        # Dilution-step QC
        # --------------------------------------------------

        standard_qc_df["Log_Quantity_Difference"] = (
            standard_qc_df["log10_Quantity"]
            .diff()
            .abs()
        )

        standard_qc_df["Observed_Delta_Cq"] = (
            standard_qc_df["Cq"].diff()
        )

        standard_qc_df["Expected_Delta_Cq"] = (
            abs(slope)
            * standard_qc_df["Log_Quantity_Difference"]
        )

        standard_qc_df["Delta_Cq_Deviation"] = (
            standard_qc_df["Observed_Delta_Cq"]
            - standard_qc_df["Expected_Delta_Cq"]
        )


        # --------------------------------------------------
        # Lowest standard 
        # --------------------------------------------------

        lowest_index = standard_qc_df.index[-1]

        lowest_quantity = standard_qc_df.loc[
            lowest_index,
            "Starting Quantity (SQ)"
        ]

        lowest_residual = standard_qc_df.loc[
            lowest_index,
            "Residual_Cq"
        ]

        lowest_delta_deviation = standard_qc_df.loc[
            lowest_index,
            "Delta_Cq_Deviation"
        ]


        # --------------------------------------------------
        # QC criteria
        # --------------------------------------------------

        regression_pass = (
            r_squared >= min_r_squared
            and
            min_efficiency
            <= efficiency
            <= max_efficiency
        )

        residual_pass = (
            abs(lowest_residual)
            <= max_cq_residual
        )

        delta_cq_pass = (
            abs(lowest_delta_deviation)
            <= max_delta_cq_deviation
        )

        qc_pass = (
            regression_pass
            and residual_pass
            and delta_cq_pass
        )


        # --------------------------------------------------
        # Save QC history
        # --------------------------------------------------

        qc_history.append(
            {
                "Number_of_Standards": len(standard_qc_df),
                "Minimum_Quantity": lowest_quantity,
                "Slope": slope,
                "Y_intercept": intercept,
                "R_squared": r_squared,
                "Efficiency_percent": efficiency,
                "Lowest_Point_Residual": lowest_residual,
                "Lowest_Delta_Cq_Deviation": lowest_delta_deviation,
                "Regression_Pass": regression_pass,
                "Residual_Pass": residual_pass,
                "Delta_Cq_Pass": delta_cq_pass,
                "QC_Pass": qc_pass
            }
        )


        # --------------------------------------------------
        # QC passed 
        # --------------------------------------------------

        if qc_pass:

            final_standard_df = standard_qc_df.copy()

            final_slope = slope
            final_intercept = intercept
            final_r_squared = r_squared
            final_efficiency = efficiency

            break


        # --------------------------------------------------
        # QC failed -> remove lowest concentration
        # --------------------------------------------------

        print(
            f"QC FAIL -> Lowest standard removed: "
            f"{lowest_quantity:.3g}"
        )

        standard_qc_df = (
            standard_qc_df
            .iloc[:-1]
            .copy()
            .reset_index(drop=True)
        )


    # ------------------------------------------------------
    # Create QC history DataFrame
    # ------------------------------------------------------

    standard_qc_results_df = pd.DataFrame(qc_history)


    # ------------------------------------------------------
    # Check final result 
    # ------------------------------------------------------

    if final_standard_df is None:

        raise ValueError(
            "Standard curve QC failed. "
            "No valid standard range was found"
        )


    return (
        final_standard_df,
        final_slope,
        final_intercept,
        final_r_squared,
        final_efficiency,
        standard_qc_results_df
    )


# Calculate absolute qPCR sample quantity
def calculate_absolute_quantity(
    sample_df,
    slope,
    intercept,
    template_volume_ul
):

    # ------------------------------------------------------
    # Calculate log10 quantity
    # ------------------------------------------------------

    result_df = sample_df.copy()

    result_df["log10_Quantity"] = (result_df["Cq"] - intercept) / slope


    # ------------------------------------------------------
    # Calculate copies per reaction
    # ------------------------------------------------------

    result_df["Copies_per_reaction"] = 10 ** result_df["log10_Quantity"]
    

    # ------------------------------------------------------
    # Calculate copies per ml
    # ------------------------------------------------------

    result_df["Copies_per_ml"] = (
        result_df["Copies_per_reaction"] / template_volume_ul * 1000
    )


    # ------------------------------------------------------
    # Calculate log10 copies per ml
    # ------------------------------------------------------

    result_df["log10_Copies_per_ml"] = np.log10(result_df["Copies_per_ml"])


    return result_df


# Create absolute qPCR summary 
def create_absolute_summary(
    sample_df,
    sample_order
):

    # ------------------------------------------------------
    # Create summary 
    # ------------------------------------------------------

    summary_df = (
        sample_df
        .groupby(
            ["Time_h", "Name"],
            as_index=False
        )
        .agg(
            Mean_Quantity=("Copies_per_ml", "mean"),
            SD_Quantity=("Copies_per_ml", "std"),
            N=("Copies_per_ml", "count")
        )
    )


    # ------------------------------------------------------
    # Apply sample order 
    # ------------------------------------------------------

    summary_df["Name"] = pd.Categorical(
        summary_df["Name"],
        categories=sample_order,
        ordered=True
    )


    # ------------------------------------------------------
    # Sort summary 
    # ------------------------------------------------------

    summary_df = (
        summary_df
        .sort_values(["Time_h", "Name"])
        .reset_index(drop=True)
    )


    return summary_df


# Prepare absolute qPCR data for statistical analysis
def prepare_absolute_statistics_data(
    sample_df,
    sample_order
):

    # ------------------------------------------------------
    # Select samples used for statistical analysis
    # ------------------------------------------------------

    stats_df = sample_df[sample_df["Name"].isin(sample_order)].copy()


    # ------------------------------------------------------
    # Remove invalid rows 
    # ------------------------------------------------------

    stats_df = stats_df.dropna(
        subset=[
            "Time_h",
            "Name",
            "log10_Copies_per_ml"
        ]
    ).copy()


    # ------------------------------------------------------
    # Determine sample order 
    # ------------------------------------------------------

    stats_sample_order = [
        sample
        for sample in sample_order
        if sample in stats_df["Name"].astype(str).unique()
    ]


    # ------------------------------------------------------
    # Apply sample order 
    # ------------------------------------------------------

    stats_df["Name"] = pd.Categorical(
        stats_df["Name"],
        categories=stats_sample_order,
        ordered=True
    )


    # ------------------------------------------------------
    # Sort data 
    # ------------------------------------------------------

    stats_df = (
        stats_df
        .sort_values(
            [
                "Time_h",
                "Name",
                "Biological_Replicate"
            ]
        )
        .reset_index(drop=True)
    )


    return stats_df, stats_sample_order


# Perform Two-way ANOVA for absolute qPCR 
def perform_absolute_two_way_anova(stats_df):

    # ------------------------------------------------------
    # Fit Two-way ANOVA model
    # ------------------------------------------------------

    anova_model = ols(
        "log10_Copies_per_ml ~ C(Time_h) * C(Name)",
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


# Perform Dunnett post-hoc test for absolute qPCR
def perform_absolute_dunnett_test(
    stats_df,
    sample_order,
    control
):

    # ------------------------------------------------------
    # Initialize results
    # ------------------------------------------------------

    dunnett_results = []


    # ------------------------------------------------------
    # Detect time points
    # ------------------------------------------------------

    time_points = sorted(stats_df["Time_h"].dropna().unique())

    # ------------------------------------------------------
    # Perform Dunnett test at each time point 
    # ------------------------------------------------------

    for time_point in time_points:

        time_df = stats_df[stats_df["Time_h"] == time_point].copy()


        # ---------------------------------------------------
        # Control data
        # ---------------------------------------------------

        control_data = (
            time_df[time_df["Name"] == control]["log10_Copies_per_ml"]
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
                time_df[time_df["Name"] == sample_name]["log10_Copies_per_ml"]
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
                    "Time_h": time_point,
                    "Control": control,
                    "Comparison": sample_name,
                    "Control_N": len(control_data),
                    "Comparison_N": len(
                        time_df[
                            time_df["Name"] == sample_name
                        ]["log10_Copies_per_ml"]
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


# Create absolute qPCR graph
def create_absolute_qpcr_graph(
    summary_df,
    dunnett_results_df,
    sample_order,
    control,
    result_folder,
    graph_title,
    xlabel,
    ylabel,
    show_significance
):

    # ------------------------------------------------------
    # Prepare graph data
    # ------------------------------------------------------

    time_points = sorted(summary_df["Time_h"].dropna().unique())

    plot_sample_order = [
        sample
        for sample in sample_order
        if sample in summary_df["Name"].astype(str).unique()
    ]

    number_of_samples = len(plot_sample_order)


    # ------------------------------------------------------
    # Automatic bar width
    # -----------------------------------------------------
    
    bar_width = min(0.25, 0.80 / max(number_of_samples, 1))


    # ------------------------------------------------------
    # X positions
    # ------------------------------------------------------

    x = np.arange(len(time_points)) * 1.5


    # ------------------------------------------------------
    # Automatic figure width
    # ------------------------------------------------------

    figure_width = max(9, 2.2 * len(time_points) + 3)

    fig, ax = plt.subplots(figsize=(figure_width, 6))


    # ------------------------------------------------------
    # Draw bars 
    # ------------------------------------------------------

    for i, sample_name in enumerate(plot_sample_order):

        sample_summary = (
            summary_df[summary_df["Name"] == sample_name]
            .set_index("Time_h")
            .reindex(time_points)
        )

        means = (
            sample_summary["Mean_Quantity"]
            .to_numpy(dtype=float)
        )

        sds = (
            sample_summary["SD_Quantity"]
            .to_numpy(dtype=float)
        )

        positions = (
            x
            + (i - (len(plot_sample_order) - 1) / 2) * bar_width
        )

        ax.bar(
            positions,
            means,
            width=bar_width,
            yerr=sds,
            capsize=4,
            label=sample_name
        )


    # ------------------------------------------------------
    # X axis
    # ------------------------------------------------------

    ax.set_xticks(x)

    ax.set_xticklabels(
        [
            f"{time:g}"
            for time in time_points
        ]
    )

    ax.set_xlabel(xlabel)


    # ------------------------------------------------------
    # Y axis data range 
    # ------------------------------------------------------

    lower_values = (
        summary_df["Mean_Quantity"]
        - summary_df["SD_Quantity"].fillna(0)
    )

    lower_values = lower_values[lower_values > 0]

    upper_values = (
        summary_df["Mean_Quantity"]
        + summary_df["SD_Quantity"].fillna(0)
    )


    if lower_values.empty:

        raise ValueError(
            "Y-axis에 사용할 양수 값을 찾을 수 없습니다."
        )


    y_min_data = lower_values.min()

    y_max_data = upper_values.max()


    min_exponent = int(np.floor(np.log10(y_min_data)))

    ax.set_ylabel(ylabel)

    ax.set_yscale("log")


    # ------------------------------------------------------
    # Legend
    # ------------------------------------------------------

    ax.legend(loc="upper left", frameon=False)


    # ------------------------------------------------------
    # Graph appearance
    # ------------------------------------------------------

    ax.set_title(graph_title)

    ax.spines["top"].set_visible(False)

    ax.spines["right"].set_visible(False)

    ax.tick_params(axis="both", direction="out")


    # ------------------------------------------------------
    # Add dunnett significance brackets
    # ------------------------------------------------------

    max_annotation_y = y_max_data


    if (
        show_significance
        and not dunnett_results_df.empty
    ):

        comparison_order = [
            sample
            for sample in plot_sample_order
            if sample != control
        ]


        sample_position_offset = {}


        for i, sample_name in enumerate(plot_sample_order):

            sample_position_offset[sample_name] = (
                i
                - (len(plot_sample_order) - 1) / 2
                ) * bar_width


        bracket_height_ratio = 1.08

        bracket_spacing_ratio = 1.35


        for time_index, time_point in enumerate(time_points):

            time_results = (
                dunnett_results_df[dunnett_results_df["Time_h"] == time_point
                ].copy()
            )

            time_summary = (
                summary_df[summary_df["Time_h"] == time_point
                ].copy()
            )


            if time_summary.empty:
                continue


            y_base = (
                time_summary["Mean_Quantity"]
                + time_summary["SD_Quantity"].fillna(0)
            ).max()


            level = 0


            if control not in sample_position_offset:
                continue


            x_control = (
                x[time_index]
                + sample_position_offset[control]
            )


            for comparison_name in comparison_order:

                comparison_result = (
                    time_results[
                        time_results["Comparison"]
                        == comparison_name
                    ]
                )


                if comparison_result.empty:
                    continue


                significance = comparison_result["Significance"].iloc[0]


                # Do not display ns
                if significance == "ns":
                    continue


                if (
                    comparison_name
                    not in sample_position_offset
                ):
                    continue


                x_comparison = (
                    x[time_index]
                    + sample_position_offset[comparison_name]
                )


                y_bracket = (
                    y_base
                    * (bracket_spacing_ratio ** level)
                    * bracket_height_ratio
                )


                y_tick_bottom = (y_bracket / 1.06)


                ax.plot(
                    [
                        x_control,
                        x_control,
                        x_comparison,
                        x_comparison
                    ],
                    [
                        y_tick_bottom,
                        y_bracket,
                        y_bracket,
                        y_tick_bottom
                    ],
                    color="black",
                    linewidth=1.2
                )


                ax.text(
                    (x_control+ x_comparison) / 2,
                    y_bracket * 1.03,
                    significance,
                    ha="center",
                    va="bottom",
                    fontsize=11,
                    color="black"
                )


                max_annotation_y = max(
                    max_annotation_y,
                    y_bracket * 1.12
                )


                level += 1


    # ------------------------------------------------------
    # Final automatic Y-axis range
    # ------------------------------------------------------

    final_y_max = max(y_max_data, max_annotation_y)


    final_max_exponent = int(np.ceil(np.log10(final_y_max)))


    ax.set_ylim(
        10 ** min_exponent,
        10 ** final_max_exponent
    )


    y_ticks = [
        10 ** exponent
        for exponent in range(
            min_exponent,
            final_max_exponent + 1
        )
    ]


    ax.set_yticks(y_ticks)

    ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))


    # ------------------------------------------------------
    # Save graph
    # ------------------------------------------------------

    fig.tight_layout()


    result_folder.mkdir(parents=True, exist_ok=True)


    graph_path = result_folder / "Absolute_qPCR_Graph.png"


    fig.savefig(
        graph_path,
        dpi=300,
        bbox_inches="tight"
    )


    plt.show()

    plt.close(fig)


    return graph_path


# Create standard curve graph
def create_standard_curve_graph(
    standard_df,
    final_standard_df,
    final_slope,
    final_intercept,
    final_r_squared,
    final_efficiency,
    result_folder
):

    # ------------------------------------------------------
    # Prepare standard curve data 
    # ------------------------------------------------------

    used_quantities = set(final_standard_df["Starting Quantity (SQ)"])


    standard_plot_df = standard_df.copy()


    standard_plot_df["QC_Status"] = np.where(
        standard_plot_df["Starting Quantity (SQ)"].isin(used_quantities),
        "Used",
        "Excluded"
    )


    # ------------------------------------------------------
    # Create figure 
    # ------------------------------------------------------

    fig_std, ax_std = plt.subplots(figsize=(7, 6))


    # ------------------------------------------------------
    # Used standards 
    # ------------------------------------------------------

    used_df = standard_plot_df[standard_plot_df["QC_Status"] == "Used"]


    ax_std.scatter(
        used_df["log10_Quantity"],
        used_df["Cq"],
        s=55,
        label="Used standards"
    )


    # ------------------------------------------------------
    # Excluded standards
    # ------------------------------------------------------

    excluded_df = standard_plot_df[standard_plot_df["QC_Status"] == "Excluded"]


    if not excluded_df.empty:

        ax_std.scatter(
            excluded_df["log10_Quantity"],
            excluded_df["Cq"],
            s=55,
            marker="x",
            label="Excluded by QC"
        )


    # ------------------------------------------------------
    # Regression line 
    # ------------------------------------------------------

    x_line = np.linspace(
        final_standard_df["log10_Quantity"].min(),
        final_standard_df["log10_Quantity"].max(),
        200
    )


    y_line = (final_slope * x_line + final_intercept)


    ax_std.plot(
        x_line,
        y_line,
        linewidth=1.5,
        label="Linear regression"
    )


    # ------------------------------------------------------
    # Equation
    # ------------------------------------------------------

    equation_text = (
        f"Cq = {final_slope:.4f} × log10(Q) "
        f"+ {final_intercept:.4f}\n"
        f"$R^2$ = {final_r_squared:.4f}\n"
        f"Efficiency = {final_efficiency:.2f}%"
    )


    ax_std.text(
        0.05,
        0.05,
        equation_text,
        transform=ax_std.transAxes,
        fontsize=10,
        verticalalignment="bottom"
    )


    # ------------------------------------------------------
    # Graph appearance
    # ------------------------------------------------------

    ax_std.set_xlabel("log10(Starting Quantity)")

    ax_std.set_ylabel("Cq")

    ax_std.set_title("qPCR Standard Curve")


    ax_std.spines["top"].set_visible(False)

    ax_std.spines["right"].set_visible(False)


    ax_std.tick_params(
        axis="both",
        direction="out"
    )


    ax_std.legend(frameon=False)


    # ------------------------------------------------------
    # Save graph
    # ------------------------------------------------------

    fig_std.tight_layout()


    result_folder.mkdir(parents=True, exist_ok=True)


    standard_curve_path = result_folder / "Standard_Curve.png"


    fig_std.savefig(
        standard_curve_path,
        dpi=300,
        bbox_inches="tight"
    )


    plt.show()

    plt.close(fig_std)


    return standard_curve_path


# Create final standard curve summary
def create_standard_curve_summary(
    final_standard_df,
    final_slope,
    final_intercept,
    final_r_squared,
    final_efficiency
):

    standard_curve_summary_df = pd.DataFrame(
        {
            "Parameter": [
                "Number_of_Standards",
                "Minimum_Quantity",
                "Maximum_Quantity",
                "Slope",
                "Y_intercept",
                "R_squared",
                "Efficiency_percent"
            ],

            "Value": [
                len(final_standard_df),

                final_standard_df["Starting Quantity (SQ)"].min(),

                final_standard_df["Starting Quantity (SQ)"].max(),

                final_slope,

                final_intercept,

                final_r_squared,

                final_efficiency
            ]
        }
    )


    return standard_curve_summary_df


# Save absolute qPCR analysis results
def save_absolute_qpcr_results(
    result_file,
    sample_df,
    summary_df,
    final_standard_df,
    standard_qc_results_df,
    standard_curve_summary_df,
    anova_results_df,
    dunnett_results_df
):

    with pd.ExcelWriter(
        result_file,
        engine="openpyxl"
    ) as writer:


        # Calculated sample data
        sample_df.to_excel(
            writer,
            sheet_name="Calculated_Samples",
            index=False
        )


        # Summary data
        summary_df.to_excel(
            writer,
            sheet_name="Summary",
            index=False
        )


        # Standards used for final curve
        final_standard_df.to_excel(
            writer,
            sheet_name="Final_Standards",
            index=False
        )


        # Standard curve QC
        standard_qc_results_df.to_excel(
            writer,
            sheet_name="Standard_QC",
            index=False
        )


        # Standard curve summary
        standard_curve_summary_df.to_excel(
            writer,
            sheet_name="Standard_Curve_Summary",
            index=False
        )


        # Two-way ANOVA
        anova_results_df.to_excel(
            writer,
            sheet_name="Two_Way_ANOVA",
            index=False
        )


        # Dunnett results
        dunnett_results_df.to_excel(
            writer,
            sheet_name="Dunnett",
            index=False
        )


    return result_file


# ============================================================
# STEP 3. EXPERIMENT ANALYSIS
# ============================================================

def analyze_experiment(experiment_folder):

    print("\n" + "=" * 60)
    print(experiment_folder.name)
    print("=" * 60)

    result_folder = experiment_folder / "Final_Results"


    try:

        # ========================================================
        # STEP 3-1. Find Excel files
        # ========================================================

        print("\n" + "=" * 60)
        print("STEP 3-1. EXCEL FILE SEARCH")
        print("=" * 60)

        excel_files = find_excel_files(experiment_folder)

        print("\n찾은 Excel 파일:")

        for file in excel_files:

            print(f" - {file.name}")


        # ========================================================
        # STEP 3-2. Select qPCR analysis type
        # ========================================================

        if ANALYSIS_TYPE == "relative":

            print("\nRelative qPCR analysis selected")

            if len(excel_files) < 2:

                print(
                    f"{experiment_folder.name}: "
                    "Relative qPCR 입력 조건과 맞지 않아 SKIP"
                )

                return


            # --------------------------------------------------------
            # Result file
            # --------------------------------------------------------

            result_file = (
                result_folder / "Relative_qPCR_Analysis_Results.xlsx"
            )


            # ========================================================
            # RELATIVE STEP 1. Load qPCR data
            # ========================================================

            print("\nReading Relative qPCR data")

            target_all, reference_all = load_relative_qpcr_data(
                excel_files=excel_files,
                target_sheet_name=RELATIVE_TARGET_SHEET,
                reference_sheet_name=RELATIVE_REFERENCE_SHEET
            )


            # ========================================================
            # RELATIVE STEP 2. Prepare qPCR data
            # ========================================================

            print("\nPreparing Relative qPCR data")

            merged = prepare_relative_qpcr_data(
                target_all=target_all,
                reference_all=reference_all,
                experiment_name=experiment_folder.name
            )

            print("\nTarget / Reference Sample Matching 정상")


            # ========================================================
            # RELATIVE STEP 3. Calculate Relative Expression
            # ========================================================

            print("\nCalculating Relative Expression")

            merged, technical_mean_df, detected_samples, sample_order = (
                calculate_relative_expression(
                    merged=merged,
                    control=CONTROL,
                    excluded_samples=EXCLUDED_SAMPLES,
                    sample_order=RELATIVE_SAMPLE_ORDER
                )
            )


            print("\n자동 인식된 Sample:")

            for sample in detected_samples:

                print(" -", sample)


            print("\nCalculated Relative Expression")

            print(
                technical_mean_df[
                    [
                        "Biological_Replicate",
                        "name",
                        "Mean_Delta_Ct",
                        "Control_Mean_Delta_Ct",
                        "DeltaDelta_Ct",
                        "Relative_Expression"
                    ]
                ].to_string(index=False)
            )


            # ========================================================
            # RELATIVE STEP 4. Statistical Analysis
            # ========================================================

            print("\n" + "=" * 60)
            print("Statistical Analysis")
            print("=" * 60)

            stats_df, stats_results_df = analyze_relative_statistics(
                technical_mean_df=technical_mean_df,
                sample_order=sample_order,
                control=CONTROL
            )


            # --------------------------------------------------------
            # Statistical analysis data
            # --------------------------------------------------------

            print("\nStatistical Analysis Data")

            print(
                stats_df[
                    [
                        "Biological_Replicate",
                        "name",
                        "Mean_Delta_Ct",
                        "DeltaDelta_Ct",
                        "Relative_Expression"
                    ]
                ].to_string(index=False)
            )


            # --------------------------------------------------------
            # Paired t-test results
            # --------------------------------------------------------

            print("\nPaired t-test")

            print(stats_results_df.to_string(index=False))


            # ========================================================
            # RELATIVE STEP 5. Create Summary
            # ========================================================

            print("\n" + "=" * 60)
            print("Summary Data")
            print("=" * 60)

            summary_df = create_relative_summary(
                technical_mean_df=technical_mean_df,
                stats_results_df=stats_results_df,
                sample_order=sample_order,
                control=CONTROL
            )

            print(summary_df.to_string(index=False))


            # ========================================================
            # RELATIVE STEP 6. Create Relative qPCR Graph
            # ========================================================

            print("\n" + "=" * 60)
            print("Creating Relative qPCR Graph")
            print("=" * 60)

            graph_path = create_relative_qpcr_graph(
                summary_df=summary_df,
                result_folder=result_folder,
                graph_title=GRAPH_TITLE,
                xlabel=XLABEL,
                ylabel=YLABEL,
                show_significance=SHOW_SIGNIFICANCE
            )

            print("\nRelative qPCR graph 저장 완료:")
            print(graph_path)


            # ========================================================
            # RELATIVE STEP 7. Save Analysis Results
            # ========================================================

            print("\nSaving Relative qPCR Results")

            result_file = save_relative_qpcr_results(
                result_file=result_file,
                merged=merged,
                technical_mean_df=technical_mean_df,
                stats_df=stats_df,
                stats_results_df=stats_results_df,
                summary_df=summary_df
            )

            print("\nExcel 결과 저장 완료:")
            print(result_file)


            print("\n" + "=" * 60)
            print(f"{experiment_folder.name} 분석 완료")
            print("=" * 60)


        elif ANALYSIS_TYPE == "absolute":

            print("\nAbsolute qPCR analysis selected")

            if len(excel_files) != 1:

                print(
                    f"{experiment_folder.name}: "
                    "Absolute qPCR 입력 조건과 맞지 않아 SKIP"
                )

                return


            # --------------------------------------------------------
            # Result file
            # --------------------------------------------------------

            result_file = (
                result_folder / "Absolute_qPCR_Analysis_Results.xlsx"
            )


            # ========================================================
            # ABSOLUTE STEP 1. Load qPCR data
            # ========================================================

            print("\nReading Absolute qPCR data")

            input_excel_file = excel_files[0]

            df = load_absolute_qpcr_data(
                excel_file=input_excel_file,
                target_sheet_name=ABSOLUTE_TARGET_SHEET
            )


            print("\nContent 종류:")
            print(df["Content"].value_counts(dropna=False))

            print("\nSample 종류:")
            print(df["Sample"].value_counts(dropna=False))


            # ========================================================
            # ABSOLUTE STEP 2. Prepare Sample Data
            # ========================================================

            print("\nPreparing Absolute qPCR sample data")

            sample_df, detected_samples, sample_order = (
                prepare_absolute_sample_data(
                    df=df,
                    sample_order=ABSOLUTE_SAMPLE_ORDER,
                    control=CONTROL,
                    excluded_samples=EXCLUDED_SAMPLES
                )
            )


            print("\n정리된 Sample 데이터:")
            print(sample_df.to_string(index=False))


            print("\n자동 인식된 Sample:")

            for sample in detected_samples:

                print(" -", sample)


            print("\n최종 Sample 순서:")

            for sample in sample_order:

                print(" -", sample)


            # ========================================================
            # ABSOLUTE STEP 3. Prepare Standard Curve Data
            # ========================================================

            print("\nPreparing Standard Curve data")

            standard_df = prepare_standard_curve_data(
                df=df,
                min_standard_points=MIN_STANDARD_POINTS
            )

            print("\nStandard 데이터:")
            print(standard_df.to_string(index=False))


            # ========================================================
            # ABSOLUTE STEP 4. Standard Curve QC
            # ========================================================

            print("\nPerforming Standard Curve QC")

            (
                final_standard_df,
                final_slope,
                final_intercept,
                final_r_squared,
                final_efficiency,
                standard_qc_results_df
            ) = perform_standard_curve_qc(
                standard_df=standard_df,
                min_standard_points=MIN_STANDARD_POINTS,
                min_r_squared=MIN_R_SQUARED,
                min_efficiency=MIN_EFFICIENCY,
                max_efficiency=MAX_EFFICIENCY,
                max_cq_residual=MAX_CQ_RESIDUAL,
                max_delta_cq_deviation=MAX_DELTA_CQ_DEVIATION
            )


            print("\n" + "=" * 60)
            print("Standard Curve QC History")
            print("=" * 60)

            print(standard_qc_results_df.to_string(index=False))


            # ========================================================
            # ABSOLUTE STEP 5. Final Standard Curve Result
            # ========================================================

            print("\n" + "=" * 60)
            print("Final Standard Curve Result")
            print("=" * 60)

            print(
                f"Number of standards : {len(final_standard_df)}"
            )

            print(
                f"Minimum quantity    : "
                f"{final_standard_df['Starting Quantity (SQ)'].min():.3g}"
            )

            print(
                f"Slope               : {final_slope:.4f}"
            )

            print(
                f"Y-intercept         : {final_intercept:.4f}"
            )

            print(
                f"R²                  : {final_r_squared:.4f}"
            )

            print(
                f"Efficiency          : {final_efficiency:.2f} %"
            )


            # ========================================================
            # ABSOLUTE STEP 6. Calculate Absolute Quantity
            # ========================================================

            print("\nCalculating sample quantity")

            sample_df = calculate_absolute_quantity(
                sample_df=sample_df,
                slope=final_slope,
                intercept=final_intercept,
                template_volume_ul=TEMPLATE_VOLUME_UL
            )


            print("\nCalculated Sample Data:")

            print(
                sample_df[
                    [
                        "Time_h",
                        "Name",
                        "Biological_Replicate",
                        "Cq",
                        "Copies_per_reaction",
                        "Copies_per_ml",
                        "log10_Copies_per_ml"
                    ]
                ].to_string(index=False)
            )


            # ========================================================
            # ABSOLUTE STEP 7. Create Summary
            # ========================================================

            print("\nCreating Absolute qPCR summary")

            summary_df = create_absolute_summary(
                sample_df=sample_df,
                sample_order=sample_order
            )


            print("\n" + "=" * 60)
            print("Summary Data")
            print("=" * 60)

            print(summary_df.to_string(index=False))


            # ========================================================
            # ABSOLUTE STEP 8. Prepare Statistical Data
            # ========================================================

            print("\nPreparing data for statistical analysis")

            stats_df, stats_sample_order = (
                prepare_absolute_statistics_data(
                    sample_df=sample_df,
                    sample_order=sample_order
                )
            )


            print("\n" + "=" * 60)
            print("Statistical Analysis Data")
            print("=" * 60)

            print(
                stats_df[
                    [
                        "Time_h",
                        "Name",
                        "Biological_Replicate",
                        "Copies_per_ml",
                        "log10_Copies_per_ml"
                    ]
                ].to_string(index=False)
            )


            # ========================================================
            # ABSOLUTE STEP 9. Two-Way ANOVA
            # ========================================================

            print("\n" + "=" * 60)
            print("Two-way ANOVA")
            print("=" * 60)

            anova_results_df = perform_absolute_two_way_anova(
                stats_df=stats_df
            )

            print(anova_results_df.to_string(index=False))


            # ========================================================
            # ABSOLUTE STEP 10. Dunnett Post-hoc Test
            # ========================================================

            print("\n" + "=" * 60)
            print("Dunnett Post-hoc Test")
            print("=" * 60)

            dunnett_results_df = perform_absolute_dunnett_test(
                stats_df=stats_df,
                sample_order=stats_sample_order,
                control=CONTROL
            )

            print(dunnett_results_df.to_string(index=False))


            # ========================================================
            # ABSOLUTE STEP 11. Create Result Folder
            # ========================================================

            result_folder.mkdir(
                parents=True,
                exist_ok=True
            )

            print("\n결과 저장 폴더:")
            print(result_folder)


            # ========================================================
            # ABSOLUTE STEP 12. Create Absolute qPCR Graph
            # ========================================================

            print("\n" + "=" * 60)
            print("Creating Absolute qPCR Graph")
            print("=" * 60)

            absolute_graph_path = create_absolute_qpcr_graph(
                summary_df=summary_df,
                dunnett_results_df=dunnett_results_df,
                sample_order=sample_order,
                control=CONTROL,
                result_folder=result_folder,
                graph_title=GRAPH_TITLE,
                xlabel=XLABEL,
                ylabel=YLABEL,
                show_significance=SHOW_SIGNIFICANCE
            )

            print("\nAbsolute qPCR graph 저장 완료:")
            print(absolute_graph_path)


            # ========================================================
            # ABSOLUTE STEP 13. Create Standard Curve Graph
            # ========================================================

            print("\n" + "=" * 60)
            print("Creating Standard Curve Graph")
            print("=" * 60)

            standard_curve_path = create_standard_curve_graph(
                standard_df=standard_df,
                final_standard_df=final_standard_df,
                final_slope=final_slope,
                final_intercept=final_intercept,
                final_r_squared=final_r_squared,
                final_efficiency=final_efficiency,
                result_folder=result_folder
            )

            print("\nStandard curve graph 저장 완료:")
            print(standard_curve_path)


            # ========================================================
            # ABSOLUTE STEP 14. Create Standard Curve Summary
            # ========================================================

            standard_curve_summary_df = (
                create_standard_curve_summary(
                    final_standard_df=final_standard_df,
                    final_slope=final_slope,
                    final_intercept=final_intercept,
                    final_r_squared=final_r_squared,
                    final_efficiency=final_efficiency
                )
            )


            # ========================================================
            # ABSOLUTE STEP 15. Save Analysis Results
            # ========================================================

            saved_result_file = save_absolute_qpcr_results(
                result_file=result_file,
                sample_df=sample_df,
                summary_df=summary_df,
                final_standard_df=final_standard_df,
                standard_qc_results_df=standard_qc_results_df,
                standard_curve_summary_df=standard_curve_summary_df,
                anova_results_df=anova_results_df,
                dunnett_results_df=dunnett_results_df
            )


            print("\nExcel 결과 저장 완료:")
            print(saved_result_file)


            print("\n" + "=" * 60)
            print(f"{experiment_folder.name} 분석 완료")
            print("=" * 60)


        else:

            raise ValueError(
                "ANALYSIS_TYPE은 'relative' 또는 'absolute' 이어야 합니다."
            )


    except Exception as error:

        print("\n" + "!" * 60)
        print(f"{experiment_folder.name} 분석 중 오류 발생")
        print(f"Error: {error}")
        print("이 Experiment는 건너뛰고 다음 Experiment를 분석합니다.")
        print("!" * 60)

        plt.close("all")

        return

# ============================================================
# MAIN 
# ============================================================

def main(): 

    # ============================================================
    # STEP 1. EXPERIMENT SEARCH
    # ============================================================

    print("\n" + "=" * 60)
    print("STEP 1. EXPERIMENT FOLDER SEARCH")
    print("=" * 60)
    
    all_experiment_folders = find_experiment_folders(BASE_FOLDER)
    
    
    print(
        f"\n찾은 Experiment 수: {len(all_experiment_folders)}"
    )
    
    for folder in all_experiment_folders:
    
        print(f" - {folder.name}")


    # ============================================================
    # STEP 2. EXPERIMENT STATUS CHECK
    # ============================================================

    print("\n" + "=" * 60)
    print("STEP 2. EXPERIMENT STATUS CHECK")
    print("=" * 60)

    (
        completed_experiments,
        new_experiment_folders
    ) = classify_experiment_folders(all_experiment_folders)

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


    # ============================================================
    # NO NEW EXPERIMENT CHECK
    # ============================================================

    if not new_experiment_folders:

        print("\n새로운 Experiment이 없습니다.")

        return


    # ============================================================
    # MAIN EXPERIMENT LOOP
    # ============================================================

    for experiment_folder in new_experiment_folders:

        analyze_experiment(experiment_folder)


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