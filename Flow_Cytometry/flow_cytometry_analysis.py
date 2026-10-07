# ============================================================
# FLOW CYTOMETRY AUTOMATED ANALYSIS PIPELINE
# ============================================================

from pathlib import Path
import pandas as pd
import numpy as np
import re
from striprtf.striprtf import rtf_to_text
from scipy import stats
import matplotlib.pyplot as plt


# ============================================================
# STEP 0. USER CONFIGURATION
# ============================================================

# ANALYSIS TYPE
# "standard" / "timecourse"
ANALYSIS_TYPE = "standard"

# Project directory
PROJECT_DIR = Path(__file__).resolve().parent

# Example data directory
BASE_FOLDER = PROJECT_DIR / "example_data" / ANALYSIS_TYPE

# Flow cytometry value column
VALUE_COLUMN = "FITC-A Geo Mean"


# ============================================================
# COMMON CONFIGURATION
# ============================================================
# Normalization
NORMALIZATION = False

NORMALIZATION_CONTROL = "Control"

# Statistics
CONTROL = "Control"

# Samples excluded from analysis
EXCLUDED_SAMPLES = [
    "Cell only",
    "EV"
]

# Graph
GRAPH_TITLE = "Gene expression"

XLABEL = "X-label"

YLABEL = "Y-label"

SHOW_SIGNIFICANCE = True


# ============================================================
# STANDARD ANALYSIS CONFIGURATION
# ============================================================
# Sample order 
SAMPLE_ORDER = [
    "EV",
    "Control",
    "Sample_A",
    "Sample_B",
    "Sample_C"
]


# ============================================================
# TIMECOURSE ANALYSIS CONFIGURATION
# ============================================================
# Timepoint order 
TIMEPOINT_ORDER = [
    6,
    12,
    24,
    48
]

# ============================================================
# BASIC VALIDATION
# ============================================================

if ANALYSIS_TYPE not in ["standard", "timecourse"]:

    raise ValueError(
        "ANALYSIS_TYPE은 'standard' 또는 'timecourse'이어야 합니다."
    )


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
# "STANDARD" NAME PARSER
# ============================================================

STANDARD_REPLICATE_PATTERN = r"^(.*?)\s*-\s*(\d+)"

def parse_sample_name(name):

    name = str(name).strip()

    # Pattern matching
    replicate_match = re.match(STANDARD_REPLICATE_PATTERN, name)

    if replicate_match is None:

        return (name, np.nan)

    # Sample
    sample = replicate_match.group(1).strip()

    # Technical replicate
    technical_replicate = int(replicate_match.group(2))


    return (sample, technical_replicate)

# ============================================================
# "TIMECOURSE" NAME PARSER
# ============================================================

TIMECOURSE_REPLICATE_PATTERN = r"^(.*?)\s*-\s*(\d+)\s*\("

TIMECOURSE_TIMEPOINT_PATTERN = r"\((\d+)\s*h\)"


def parse_timecourse_name(name):

    name = str(name).strip()

    # Pattern matching
    replicate_match = re.search(TIMECOURSE_REPLICATE_PATTERN, name)

    timepoint_match = re.search(TIMECOURSE_TIMEPOINT_PATTERN, name)

    if (replicate_match is None or timepoint_match is None):

        return (None, None, None)

    # Sample
    sample = replicate_match.group(1).strip()

    # Timepoint 
    timepoint = int(timepoint_match.group(1))

    # Biological replicate
    biological_replicate = int(replicate_match.group(2))


    return (sample, timepoint, biological_replicate)


# ============================================================
# RTF ENCODING READER
# ============================================================
def read_rtf_file(rtf_file):

    encodings = [
        "cp949",
        "utf-8",
        "utf-8-sig",
        "euc-kr"
    ]

    for encoding in encodings:

        try:

            with open(rtf_file, "r", encoding=encoding) as f:

                return f.read()

        except UnicodeDecodeError:

            continue

    raise ValueError(
        f"RTF 파일 encoding을 확인할 수 없습니다: {rtf_file.name}"
    )


# RTF TABLE PARSER
def extract_gmfi_from_rtf(rtf_file, experiment,rtf_index):

    # Read RTF
    rtf_content = read_rtf_file(rtf_file)

    text_content = rtf_to_text(rtf_content)

    lines = text_content.splitlines()

    # Header search
    header_index = None

    for i, line in enumerate(lines):

        if VALUE_COLUMN in line:

            header_index = i
            break


    if header_index is None:

        raise ValueError(
            f"{VALUE_COLUMN} header를 찾지 못했습니다: {rtf_file.name}"
        )

    header_line = lines[header_index]

    # Header
    header = [
        column.strip()
        for column in header_line.split("|")
    ]

    # gMFI column
    try:

        gMFI_index = header.index(VALUE_COLUMN)

    except ValueError:

        raise ValueError(
            f"{VALUE_COLUMN} column을 찾지 못했습니다: {rtf_file.name}"
        )


    # Data rows
    records = []

    for line in lines[header_index + 1:]:

        if not line.strip():

            continue

        columns = [
            column.strip()
            for column in line.split("|")
        ]

        # column 수 검사
        if len(columns) <= gMFI_index:

            continue

        name = columns[0].strip()

        gmfi_value = columns[gMFI_index].replace(",", "").strip()

        try:

            gmfi = float(gmfi_value)

        except ValueError:

            continue


        # ====================================================
        # STANDARD
        # ====================================================

        if ANALYSIS_TYPE == "standard":

            sample, technical_replicate = parse_sample_name(name)

            records.append(
                {
                    "Experiment": experiment,

                    "RTF_File": rtf_file.name,

                    "Name": name,

                    "Sample": sample,

                    "Technical_Replicate": technical_replicate,

                    "RTF_Index": rtf_index,

                    "gMFI": gmfi
                }
            )


        # ====================================================
        # TIMECOURSE
        # ====================================================

        elif ANALYSIS_TYPE == "timecourse":

            sample, timepoint, biological_replicate = (
                parse_timecourse_name(name)
            )


            if sample is None:

                print(
                    f"  [WARNING] "
                    f"Time-course parsing 실패: {name}"
                )

                continue


            records.append(
                {
                    "Experiment": experiment,

                    "RTF_File": rtf_file.name,

                    "Name": name,

                    "Sample": sample,

                    "Timepoint": timepoint,

                    "Biological_Replicate": biological_replicate,

                    "gMFI": gmfi
                }
            )


    return records


# ============================================================
# EXPERIMENT MANAGEMENT FUNCTIONS
# ============================================================

# ============================================================
# EXPERIMENT FOLDER SEARCH
# ============================================================

def find_experiment_folders(base_folder):

    if not base_folder.exists():

        raise FileNotFoundError(
            f"\nBASE_FOLDER를 찾을 수 없습니다: {base_folder}"
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

    completion_marker = "Flow_Cytometry_Analysis_Results.xlsx"

    for folder in experiment_folders:

        marker_file = folder / "Final_Results" / completion_marker

        if marker_file.exists():

            completed_experiments.append(folder)

        else:

            new_experiment_folders.append(folder)

    return (completed_experiments, new_experiment_folders)


# ============================================================
# RTF FILE SEARCH
# ============================================================
def find_rtf_files(experiment_folder):

    rtf_files = sorted(
        [
            file
            for file in experiment_folder.iterdir()
            if file.is_file()
            and file.suffix.lower() == ".rtf"
            and not file.name.startswith("~$")
        ]
    )

    if not rtf_files:

        raise ValueError(
            f"{experiment_folder.name}에 RTF 파일이 없습니다."
        )

    return rtf_files


# ============================================================
# EXPERIMENT DATA EXTRACTION
# ============================================================
def extract_experiment_data(rtf_files, experiment_name):

    records = []

    for rtf_index, rtf_file in enumerate(rtf_files, start=1):

        print(
            f"\n[RTF {rtf_index}] "
            f"{rtf_file.name}"
        )

        extracted = extract_gmfi_from_rtf(
            rtf_file, experiment_name, rtf_index
        )

        print(
            f"  추출된 row: {len(extracted)}"
        )

        records.extend(extracted)

    # --------------------------------------------------------
    # Extraction validation
    # --------------------------------------------------------
    if not records:

        raise ValueError(
            f"{experiment_name}에서 gMFI 데이터가 추출되지 않았습니다."
        )

    df = pd.DataFrame(records)

    return df


# ============================================================
# STANDARD ANALYSIS FUNCTIONS
# ============================================================
# STANDARD Data cleaning
def clean_standard_data(df):

    df = df.copy()

    df["gMFI"] = pd.to_numeric(df["gMFI"], errors="coerce")

    df["Sample"] = df["Sample"].astype(str).str.strip()

    df_valid = df.dropna(subset=["gMFI", "Sample"]).copy()

    return df, df_valid


# STANDARD Sample order
def get_standard_sample_order(df_valid):

    found_samples = sorted(df_valid["Sample"].unique())

    current_sample_order = [
        sample
        for sample in SAMPLE_ORDER
        if sample in found_samples
    ]


    return (found_samples, current_sample_order)


# STANDARD Technical replicate summary
def create_technical_summary(df_valid, sample_order):

    technical_summary = (
        df_valid
        .groupby(
            [
                "Experiment",
                "RTF_Index",
                "RTF_File",
                "Sample"
            ],
            as_index=False
        )
        .agg(
            Technical_N=("gMFI", "count"),
            Mean_gMFI=("gMFI", "mean"),
            SD_gMFI=("gMFI", "std")
        )
    )

    technical_summary["Sample"] = pd.Categorical(
        technical_summary["Sample"],
        categories=sample_order,
        ordered=True
    )

    technical_summary = (
        technical_summary
        .sort_values(
            [
                "RTF_Index",
                "Sample"
            ]
        )
        .reset_index(drop=True)
    )

    return technical_summary


# STANDARD Normalization
def normalize_standard_data(technical_summary):

    if NORMALIZATION:

        control = (
            technical_summary[
                technical_summary["Sample"] == NORMALIZATION_CONTROL
            ][
                [
                    "RTF_Index",
                    "Mean_gMFI"
                ]
            ]
            .rename(
                columns={
                    "Mean_gMFI": "Control_gMFI"
                }
            )
        )

        if control.empty:

            raise ValueError(
                f"Normalization control "
                f"'{NORMALIZATION_CONTROL}'을 찾지 못했습니다."
            )

        if (control["Control_gMFI"] <= 0).any():

            raise ValueError(
                "Control gMFI가 0 이하인 데이터가 존재합니다."
            )

        normalized_summary = (
            technical_summary
            .merge(
                control,
                on="RTF_Index",
                how="left"
            )
        )

        normalized_summary["Analysis_Value"] = (
            normalized_summary["Mean_gMFI"]
            /
            normalized_summary["Control_gMFI"]
        )

    else:

        normalized_summary = technical_summary.copy()

        normalized_summary["Control_gMFI"] = np.nan

        normalized_summary["Analysis_Value"] = normalized_summary["Mean_gMFI"]

    return normalized_summary


# STANDARD Final summary
def create_standard_final_summary(normalized_summary, sample_order):

    final_summary = (
        normalized_summary
        .groupby("Sample", observed=True)["Analysis_Value"]
        .agg(
            N="count",
            Mean="mean",
            SD="std"
        )
        .reset_index()
    )

    final_summary["Sample"] = pd.Categorical(
        final_summary["Sample"],
        categories=sample_order,
        ordered=True
    )

    final_summary = final_summary.sort_values("Sample").reset_index(drop=True)

    return final_summary


# STANDARD Statistics
def run_standard_statistics(normalized_summary, sample_order):

    stat_results = []

    control_values = (
        normalized_summary[
            normalized_summary["Sample"] == CONTROL
        ][
            [
                "RTF_Index",
                "Analysis_Value"
            ]
        ]
        .rename(
            columns={
                "Analysis_Value": "Control"
            }
        )
    )

    for sample in sample_order:

        if sample == CONTROL:

            continue

        if sample in EXCLUDED_SAMPLES:

            continue

        treatment = (
            normalized_summary[
                normalized_summary["Sample"] == sample
            ][
                [
                    "RTF_Index",
                    "Analysis_Value"
                ]
            ]
            .rename(
                columns={
                    "Analysis_Value": "Treatment"
                }
            )
        )

        paired = pd.merge(
            control_values,
            treatment,
            on="RTF_Index",
            how="inner"
        )

        paired = paired.dropna(
            subset=[
                "Control",
                "Treatment"
            ]
        )

        n = len(paired)

        if n >= 2:

            t_stat, p_value = stats.ttest_rel(
                paired["Treatment"],
                paired["Control"]
            )

            star = significance_star(p_value)

        else:

            t_stat = np.nan
            p_value = np.nan
            star = "NA"

        stat_results.append(
            {
                "Sample": sample,
                "N": n,
                "t_stat": t_stat,
                "p_value": p_value,
                "Significance": star
            }
        )

    stat_results = pd.DataFrame(stat_results)

    return stat_results


# STANDARD Graph
def plot_standard_results(final_summary, stat_results, result_folder):

    plot_summary = final_summary.copy()

    means = plot_summary["Mean"].astype(float).to_numpy()

    sds = plot_summary["SD"].astype(float).fillna(0).to_numpy()


    # --------------------------------------------------------
    # Remove EV error bar 
    # --------------------------------------------------------

    for i, sample in enumerate(plot_summary["Sample"]):

        if str(sample) == "EV":

            sds[i] = 0


    # --------------------------------------------------------
    # Bar graph
    # --------------------------------------------------------
    x = np.arange(len(plot_summary))

    fig, ax = plt.subplots(figsize=(10, 8))

    ax.bar(
        x,
        means,
        yerr=sds,
        capsize=5
    )

    ax.set_xticks(x)

    ax.set_xticklabels(
        plot_summary["Sample"]
        .astype(str),
        rotation=45,
        ha="right"
    )

    ax.set_xlabel(XLABEL)

    ax.set_ylabel(YLABEL)

    ax.set_title(GRAPH_TITLE)


    # --------------------------------------------------------
    # Significance
    # --------------------------------------------------------
    if SHOW_SIGNIFICANCE:

        for i, sample in enumerate(plot_summary["Sample"]):

            sample = str(sample)

            if sample == CONTROL:

                continue

            if sample in EXCLUDED_SAMPLES:

                continue

            result = stat_results[stat_results["Sample"] == sample]

            if result.empty:

                continue

            star = result["Significance"].iloc[0]

            if star == "NA":

                continue

            mean_value = means[i]

            sd_value = sds[i]


            y_position = mean_value + sd_value

            if y_position <= 0:

                y_position = mean_value + 0.1

            significance_y = (
                y_position * (1 + 0.08)
            )

            ax.text(
                i,
                significance_y,
                star,
                ha="center",
                va="bottom",
                fontsize=14
            )


    # --------------------------------------------------------
    # Automatic Y-axis margin
    # --------------------------------------------------------
    data_top = (plot_summary["Mean"] + plot_summary["SD"].fillna(0)).max()

    annotation_positions = [
        text_obj.get_position()[1]
        for text_obj in ax.texts
        if text_obj.get_text()
        in {
            "ns",
            "*",
            "**",
            "***",
            "****"
        }
    ]

    if annotation_positions:

        required_top = max(annotation_positions)

    else:

        required_top = data_top

    if (
        not np.isfinite(required_top)
        or
        required_top <= 0
    ):

        required_top = 1

    ax.set_ylim(
        bottom=0,
        top=required_top * 1.18
    )


    # --------------------------------------------------------
    # Graph path 
    # --------------------------------------------------------

    if NORMALIZATION:

        graph_path = result_folder /"Normalized_gMFI_Mean_SD.png"
        

    else:

        graph_path = result_folder /"Raw_gMFI_Mean_SD.png"


    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------
    plt.tight_layout()

    plt.savefig(
        graph_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()

    plt.close()

    return graph_path


# STANDARD Result save
def save_standard_results(
    df,
    technical_summary,
    normalized_summary,
    final_summary,
    stat_results,
    result_folder
):

    result_excel_path = result_folder / "Flow_Cytometry_Analysis_Results.xlsx"
    

    with pd.ExcelWriter(
        result_excel_path,
        engine="openpyxl"
    ) as writer:

        df.to_excel(
            writer,
            sheet_name="Raw_Data",
            index=False
        )

        technical_summary.to_excel(
            writer,
            sheet_name="Technical_Summary",
            index=False
        )

        normalized_summary.to_excel(
            writer,
            sheet_name="Analysis_Data",
            index=False
        )

        final_summary.to_excel(
            writer,
            sheet_name="Final_Summary",
            index=False
        )

        stat_results.to_excel(
            writer,
            sheet_name="Statistics",
            index=False
        )

    return result_excel_path


# STANDARD Final validation
def validate_standard_results(
    df_valid,
    normalized_summary,
    result_excel_path,
    graph_path
):

    validation_passed = True

    # --------------------------------------------------------
    # Valid data
    # --------------------------------------------------------

    if df_valid.empty:

        validation_passed = False

    # --------------------------------------------------------
    # Sample
    # -------------------------------------------------------- 

    found_samples = set(df_valid["Sample"].unique())

    unknown_samples = sorted(found_samples - set(SAMPLE_ORDER))

    if unknown_samples:

        print(
            f"Warning: SAMPLE_ORDER에 없는 sample이 있습니다: "
            f"{unknown_samples}"
        )


    # --------------------------------------------------------
    # Statistical control
    # --------------------------------------------------------

    if CONTROL not in found_samples:

        print(
            f"Warning: Statistical control "
            f"'{CONTROL}'을 찾지 못했습니다."
        )

        validation_passed = False


    # --------------------------------------------------------
    # Analysis value
    # --------------------------------------------------------

    normalization_failed = (
        normalized_summary["Analysis_Value"].isna()
        |
        ~np.isfinite(
            normalized_summary["Analysis_Value"]
        )
    )

    if normalization_failed.any():

        validation_passed = False


    # --------------------------------------------------------
    # Saved files 
    # --------------------------------------------------------

    saved_files = [result_excel_path, graph_path]

    for file in saved_files:

        if not file.exists():

            validation_passed = False

    return validation_passed


# ============================================================
# STANDARD EXPERIMENT ANALYSIS
# ============================================================
def analyze_standard_experiment(df, result_folder):

    print("\n" + "=" * 60)
    print("STANDARD ANALYSIS")
    print("=" * 60)

    # ========================================================
    # STEP 3-3. RAW DATA CLEANING
    # ========================================================

    print("\n" + "=" * 60)
    print("STANDARD STEP 3-3. RAW DATA CLEANING")
    print("=" * 60)

    df, df_valid = clean_standard_data(df)

    print(
        f"전체 raw rows: {len(df)}"
    )

    print(
        f"유효 gMFI rows: {len(df_valid)}"
    )

    print(
        f"Missing gMFI: {df['gMFI'].isna().sum()}"
    )

    # ========================================================
    # STEP 3-4. SAMPLE VALIDATION
    # ========================================================

    print("\n" + "=" * 60)
    print("STANDARD STEP 3-4. SAMPLE ORDER")
    print("=" * 60)

    (
        found_samples,
        current_sample_order
    ) = get_standard_sample_order(df_valid)

    print("\n찾은 Sample:")

    for sample in found_samples:

        print(f"  - {sample}")

    print("\n현재 Sample 순서:")

    print(current_sample_order)

    # ========================================================
    # STEP 3-5. TECHNICAL SUMMARY
    # ========================================================

    technical_summary = create_technical_summary(
        df_valid,
        current_sample_order
    )

    # ========================================================
    # STEP 3-6. NORMALIZATION
    # ========================================================

    normalized_summary = normalize_standard_data(technical_summary)

    # ========================================================
    # STEP 3-7. FINAL SUMMARY
    # ========================================================

    final_summary = create_standard_final_summary(
        normalized_summary,
        current_sample_order
    )

    # ========================================================
    # STEP 3-8. STATISTICS
    # ========================================================

    stat_results = run_standard_statistics(
        normalized_summary,
        current_sample_order
    )

    # ========================================================
    # STEP 3-9. GRAPH
    # ========================================================

    graph_path = plot_standard_results(
        final_summary,
        stat_results,
        result_folder
    )

    # ========================================================
    # STEP 3-10. RESULT SAVE
    # ========================================================

    result_excel_path = save_standard_results(
        df,
        technical_summary,
        normalized_summary,
        final_summary,
        stat_results,
        result_folder
    )

    # ========================================================
    # STEP 3-11. FINAL VALIDATION
    # ========================================================

    validation_passed = validate_standard_results(
        df_valid,
        normalized_summary,
        result_excel_path,
        graph_path
    )

    return validation_passed


# ============================================================
# TIMECOURSE ANALYSIS FUNCTIONS
# ============================================================
# TIMECOURSE Data cleaning
def clean_timecourse_data(df):

    df = df.copy()

    df["gMFI"] = pd.to_numeric(df["gMFI"], errors="coerce")

    df["Sample"] = (
        df["Sample"]
        .astype(str)
        .str.strip()
    )

    df["Timepoint"] = pd.to_numeric(df["Timepoint"], errors="coerce")

    df["Biological_Replicate"] = pd.to_numeric(
        df["Biological_Replicate"],
        errors="coerce"
    )

    df_valid = df.dropna(
        subset=[
            "gMFI",
            "Sample",
            "Timepoint",
            "Biological_Replicate"
        ]
    ).copy()

    return df, df_valid


# TIMECOURSE Timepoint order
def get_timepoint_order(df_valid):

    found_timepoints = sorted(df_valid["Timepoint"].unique())

    unexpected_timepoints = [
        timepoint
        for timepoint in found_timepoints
        if timepoint not in TIMEPOINT_ORDER
    ]

    current_timepoint_order = [
        timepoint
        for timepoint in TIMEPOINT_ORDER
        if timepoint in found_timepoints
    ]


    return (
        found_timepoints,
        unexpected_timepoints,
        current_timepoint_order
    )


# TIMECOURSE Sample validation
def get_timecourse_samples(df_valid):

    found_samples = sorted(df_valid["Sample"].unique())

    return found_samples


# TIMECOURSE Biological replicate validation
def validate_biological_replicates(df_valid):

    found_replicates = sorted(
        df_valid["Biological_Replicate"]
        .astype(int)
        .unique()
    )

    replicate_check = (
        df_valid
        .groupby(
            [
                "Sample",
                "Timepoint"
            ]
        )["Biological_Replicate"]
        .nunique()
        .reset_index(name="Biological_N")
    )

    return (found_replicates, replicate_check)


# TIMECOURSE Biological replicate summary
def create_biological_summary(df_valid, timepoint_order):

    biological_summary = (
        df_valid
        .groupby(
            [
                "Experiment",
                "Sample",
                "Timepoint",
                "Biological_Replicate"
            ],
            as_index=False
        )
        .agg(
            N_gMFI=("gMFI", "count"),
            Mean_gMFI=("gMFI", "mean"),
            SD_gMFI=("gMFI", "std")
        )
    )

    biological_summary["Timepoint"] = pd.Categorical(
        biological_summary["Timepoint"],
        categories=timepoint_order,
        ordered=True
    )

    biological_summary = (
        biological_summary
        .sort_values(
            [
                "Sample",
                "Timepoint",
                "Biological_Replicate"
            ]
        )
        .reset_index(drop=True)
    )

    return biological_summary


# TIMECOURSE Normalization
def normalize_timecourse_data(biological_summary):

    biological_summary = biological_summary.copy()

    if NORMALIZATION:

        # ----------------------------------------------------
        # Extract normalization control
        # ----------------------------------------------------

        timecourse_control = (
            biological_summary[
                biological_summary["Sample"]
                == NORMALIZATION_CONTROL
            ][
                [
                    "Timepoint",
                    "Biological_Replicate",
                    "Mean_gMFI"
                ]
            ]
            .rename(
                columns={
                    "Mean_gMFI": "Control_gMFI"
                }
            )
        )

        if timecourse_control.empty:

            raise ValueError(
                f"TIMECOURSE normalization control "
                f"'{NORMALIZATION_CONTROL}'을 찾지 못했습니다."
            )

        if (timecourse_control["Control_gMFI"] <= 0).any():

            raise ValueError(
                "TIMECOURSE Control gMFI가 0 이하인 데이터가 존재합니다."
            )

        biological_summary = (
            biological_summary
            .merge(
                timecourse_control,
                on=[
                    "Timepoint",
                    "Biological_Replicate"
                ],
                how="left"
            )
        )

        if biological_summary["Control_gMFI"].isna().any():

            raise ValueError(
                "일부 Timepoint / Biological Replicate에서 "
                "normalization control을 찾지 못했습니다."
            )

        biological_summary["Analysis_Value"] = (
            biological_summary["Mean_gMFI"]
            /
            biological_summary["Control_gMFI"]
        )

    else:

        biological_summary["Control_gMFI"] = np.nan

        biological_summary["Analysis_Value"] = (
            biological_summary["Mean_gMFI"]
        )

    return biological_summary


# TIMECOURSE Final summary
def create_timecourse_summary(biological_summary, timepoint_order):

    timecourse_summary = (
        biological_summary
        .groupby(
            [
                "Sample",
                "Timepoint"
            ],
            observed=True
        )["Analysis_Value"]
        .agg(
            N="count",
            Mean="mean",
            SD="std"
        )
        .reset_index()
    )

    timecourse_summary["Timepoint"] = pd.Categorical(
        timecourse_summary["Timepoint"],
        categories=timepoint_order,
        ordered=True
    )

    timecourse_summary = (
        timecourse_summary
        .sort_values(
            [
                "Sample",
                "Timepoint"
            ]
        )
        .reset_index(drop=True)
    )

    return timecourse_summary


# TIMECOURSE Statistics
def run_timecourse_statistics(
    biological_summary,
    found_samples,
    timepoint_order
):

    timecourse_stat_results = []

    statistical_control = CONTROL

    # -------------------------------------------------------
    # Analyze each timepoint 
    # -------------------------------------------------------

    for timepoint in timepoint_order:

        # ---------------------------------------------------
        # Control at current timepoint 
        # ---------------------------------------------------

        control = (
            biological_summary[
                (biological_summary["Sample"] == statistical_control
                )
                &
                (biological_summary["Timepoint"] == timepoint
                )
            ][
                [
                    "Biological_Replicate",
                    "Analysis_Value"
                ]
            ]
            .rename(
                columns={
                    "Analysis_Value": "Control"
                }
            )
        )


        # ---------------------------------------------------
        # Compare each sample with control
        # ---------------------------------------------------

        for sample in found_samples:

            if sample == statistical_control:

                continue

            if sample in EXCLUDED_SAMPLES:

                continue

            treatment = (
                biological_summary[
                    (biological_summary["Sample"] == sample)
                    &
                    (biological_summary["Timepoint"] == timepoint)
                ][
                    [
                        "Biological_Replicate",
                        "Analysis_Value"
                    ]
                ]
                .rename(
                    columns={
                        "Analysis_Value": "Treatment"
                    }
                )
            )


            # ---------------------------------------------------
            # Pair biological replicates
            # ---------------------------------------------------

            paired = pd.merge(
                control,
                treatment,
                on="Biological_Replicate",
                how="inner"
            )

            paired = paired.dropna(
                subset=[
                    "Control",
                    "Treatment"
                ]
            )

            n = len(paired)


            # ---------------------------------------------------
            # Paired t-test
            # ---------------------------------------------------

            if n >= 2:

                t_stat, p_value = stats.ttest_rel(
                    paired["Treatment"],
                    paired["Control"]
                )

                star = significance_star(p_value)

            else:

                t_stat = np.nan
                p_value = np.nan
                star = "NA"


            # ---------------------------------------------------
            # Save result
            # ---------------------------------------------------

            timecourse_stat_results.append(
                {
                    "Timepoint": int(timepoint),
                    "Control": statistical_control,
                    "Sample": sample,
                    "N": n,
                    "t_stat": t_stat,
                    "p_value": p_value,
                    "Significance": star
                }
            )

    return pd.DataFrame(timecourse_stat_results)


# TIMECOURSE Graph
def plot_timecourse_results(
    timecourse_summary,
    timecourse_stat_results,
    found_samples,
    timepoint_order,
    result_folder
):

    fig, ax = plt.subplots(figsize=(12, 8))

    # -------------------------------------------------------
    # Sample plot 
    # -------------------------------------------------------

    for sample in found_samples:

        plot_data = (
            timecourse_summary[
                timecourse_summary["Sample"] == sample
            ]
            .sort_values("Timepoint")
            .copy()
        )

        if plot_data.empty:

            continue

        x_values = plot_data["Timepoint"].astype(int).to_numpy()

        y_values = plot_data["Mean"].astype(float).to_numpy()

        error_values = (
            plot_data["SD"]
            .astype(float)
            .fillna(0)
            .to_numpy()
        )


        # ---------------------------------------------------
        # Line + error bar
        # ---------------------------------------------------

        ax.errorbar(
            x_values,
            y_values,
            yerr=error_values,
            marker="o",
            capsize=5,
            linewidth=2,
            label=sample
        )


    # -------------------------------------------------------
    # Maximum value at each timepoint
    # -------------------------------------------------------

    timepoint_max = {}

    for timepoint in timepoint_order:

        current_data = timecourse_summary[
            timecourse_summary["Timepoint"] == timepoint
        ]

        if current_data.empty:

            continue

        current_max = (
            current_data["Mean"]
            +
            current_data["SD"].fillna(0)
        ).max()

        timepoint_max[timepoint] = current_max


    # -------------------------------------------------------
    # Significance annotation
    # -------------------------------------------------------

    if (
        SHOW_SIGNIFICANCE
        and not timecourse_stat_results.empty
    ):

        for _, row in (timecourse_stat_results.iterrows()):

            star = row["Significance"]

            if star in ["NA", "ns"]:

                continue

            sample = row["Sample"]

            timepoint = int(row["Timepoint"])

            sample_data = (
                timecourse_summary[
                    (timecourse_summary["Sample"] == sample)
                    &
                    (timecourse_summary["Timepoint"] == timepoint)
                ]
            )

            if sample_data.empty:

                continue

            mean_value = float(sample_data["Mean"].iloc[0])

            sd_value = sample_data["SD"].iloc[0]

            if pd.isna(sd_value):

                sd_value = 0


            # -----------------------------------------------
            # Horizontal offset
            # -----------------------------------------------

            sample_index = found_samples.index(sample)

            x_offset = (
                sample_index
                - (len(found_samples) - 1) / 2
            ) * 0.08


            # -----------------------------------------------
            # Vertical position
            # -----------------------------------------------

            base_y = mean_value + sd_value
            

            global_max = (
                timepoint_max.get(
                    timepoint,
                    base_y
                )
            )

            y_position = max(base_y, global_max)

            if y_position <= 0:

                y_position = 1

            y_position = (
                y_position
                *
                (1 + 0.08)
            )

            ax.text(
                timepoint + x_offset,
                y_position,
                star,
                ha="center",
                va="bottom",
                fontsize=12,
                fontweight="bold"
            )

    # -------------------------------------------------------
    # Graph format 
    # -------------------------------------------------------

    ax.set_xticks(timepoint_order)

    ax.set_xticklabels(
        [
            f"{int(x)} h"
            for x in timepoint_order
        ]
    )

    ax.set_title(GRAPH_TITLE)

    ax.set_xlabel(XLABEL)

    ax.set_ylabel(YLABEL)

    ax.set_ylim(bottom=0)

    ax.legend(
        title="Sample",
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )


    # -------------------------------------------------------
    # Automatic Y-axis margin
    # -------------------------------------------------------

    current_ylim = ax.get_ylim()

    upper_limit = current_ylim[1]

    if (
        SHOW_SIGNIFICANCE
        and timepoint_max
    ):

        max_annotation = (
            max(timepoint_max.values())
            * (1 + 0.08)
            * 1.15
        )

        if max_annotation > upper_limit:

            upper_limit = max_annotation

    if upper_limit <= 0:

        upper_limit = 1

    ax.set_ylim(
        0,
        upper_limit
    )


    # -------------------------------------------------------
    # Graph path
    # -------------------------------------------------------

    if NORMALIZATION:

        graph_path = result_folder / "Timecourse_Normalized_gMFI_Mean_SD.png"

    else:

        graph_path = result_folder / "Timecourse_gMFI_Mean_SD.png"


    # -------------------------------------------------------
    # Save graph
    # -------------------------------------------------------

    plt.tight_layout()

    plt.savefig(
        graph_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()

    plt.close()

    return graph_path


# TIMECOURSE Result save
def save_timecourse_results(
    df,
    biological_summary,
    timecourse_summary,
    timecourse_stat_results,
    result_folder
):

    result_excel_path = result_folder / "Flow_Cytometry_Analysis_Results.xlsx"

    with pd.ExcelWriter(
        result_excel_path,
        engine="openpyxl"
    ) as writer:

        df.to_excel(
            writer,
            sheet_name="Raw_Data",
            index=False
        )

        biological_summary.to_excel(
            writer,
            sheet_name="Biological_Summary",
            index=False
        )

        timecourse_summary.to_excel(
            writer,
            sheet_name="Timecourse_Summary",
            index=False
        )

        timecourse_stat_results.to_excel(
            writer,
            sheet_name="Statistics",
            index=False
        )

    return result_excel_path


# TIMECOURSE Final validation
def validate_timecourse_results(
    df_valid,
    biological_summary,
    found_timepoints,
    found_samples,
    result_excel_path,
    graph_path
):

    validation_passed = True

    # --------------------------------------------------------
    # Valid data
    # --------------------------------------------------------

    if df_valid.empty:

        validation_passed = False


    # --------------------------------------------------------
    # Timepoint
    # --------------------------------------------------------

    missing_timepoints = [
        timepoint
        for timepoint in TIMEPOINT_ORDER
        if timepoint not in found_timepoints
    ]

    unexpected_timepoints = [
        timepoint
        for timepoint in found_timepoints
        if timepoint not in TIMEPOINT_ORDER
    ]

    if missing_timepoints:

        print(
            f"Warning: Missing timepoints: "
            f"{missing_timepoints}"
        )

    if unexpected_timepoints:

        print(
            f"Warning: Unexpected timepoints: "
            f"{unexpected_timepoints}"
        )


    # --------------------------------------------------------
    # Statistical control
    # --------------------------------------------------------

    if CONTROL not in found_samples:

        print(
            f"Warning: Statistical control "
            f"'{CONTROL}'을 찾지 못했습니다."
        )

        validation_passed = False


    # --------------------------------------------------------
    # Analysis value
    # --------------------------------------------------------

    analysis_failed = (
        biological_summary["Analysis_Value"].isna()
        |
        ~np.isfinite(biological_summary["Analysis_Value"])
    )

    if analysis_failed.any():

        validation_passed = False


    # --------------------------------------------------------
    # Saved files
    # --------------------------------------------------------

    saved_files = [result_excel_path, graph_path]

    for file in saved_files:

        if not file.exists():

            validation_passed = False

    return validation_passed


# ============================================================
# TIMECOURSE EXPERIMENT ANALYSIS
# ============================================================

def analyze_timecourse_experiment(
    df,
    result_folder
):

    print("\n" + "=" * 60)
    print("TIMECOURSE ANALYSIS")
    print("=" * 60)

    # ========================================================
    # STEP 3-3. RAW DATA CLEANING
    # ========================================================

    print("\n" + "=" * 60)
    print("TIMECOURSE STEP 3-3. RAW DATA CLEANING")
    print("=" * 60)

    df, df_valid = clean_timecourse_data(df)

    print(f"전체 raw rows: {len(df)}")

    print(f"유효 rows: {len(df_valid)}")

    print(f"Missing gMFI: {df['gMFI'].isna().sum()}")

    # ========================================================
    # STEP 3-4. TIMEPOINT VALIDATION
    # ========================================================

    print("\n" + "=" * 60)
    print("TIMECOURSE STEP 3-4. TIMEPOINT VALIDATION")
    print("=" * 60)

    (
        found_timepoints,
        unexpected_timepoints,
        current_timepoint_order
    ) = get_timepoint_order(df_valid)

    print(f"\n찾은 Timepoint: {found_timepoints}")

    if len(unexpected_timepoints) > 0:

        print(f"\n예상하지 않은 Timepoint: {unexpected_timepoints}")

    print(f"\n현재 Timepoint 순서: {current_timepoint_order}")

    # ========================================================
    # STEP 3-5. SAMPLE / REPLICATE VALIDATION
    # ========================================================

    print("\n" + "=" * 60)
    print("TIMECOURSE STEP 3-5. SAMPLE / REPLICATE VALIDATION")
    print("=" * 60)

    found_samples = get_timecourse_samples(df_valid)

    print("\n찾은 Sample:")

    for sample in found_samples:

        print(f"  - {sample}")

    (
        found_replicates,
        replicate_check
    ) = validate_biological_replicates(df_valid)

    print("\n찾은 Biological Replicate:")

    for rep in found_replicates:

        print(f"  - BR{rep}")

    print("\n===== REPLICATE CHECK =====")

    print(replicate_check.to_string(index=False))

    # ========================================================
    # STEP 3-6. BIOLOGICAL REPLICATE SUMMARY
    # ========================================================

    biological_summary = create_biological_summary(
        df_valid,
        current_timepoint_order
    )


    # ========================================================
    # STEP 3-7. NORMALIZATION
    # ========================================================

    biological_summary = normalize_timecourse_data(biological_summary)

    # ========================================================
    # STEP 3-8. FINAL SUMMARY
    # ========================================================

    timecourse_summary = create_timecourse_summary(
        biological_summary,
        current_timepoint_order
    )


    # ========================================================
    # STEP 3-9. STATISTICS
    # ========================================================

    timecourse_stat_results = run_timecourse_statistics(
        biological_summary,
        found_samples,
        current_timepoint_order
    )


    # ========================================================
    # STEP 3-10. GRAPH
    # ========================================================

    graph_path = plot_timecourse_results(
        timecourse_summary,
        timecourse_stat_results,
        found_samples,
        current_timepoint_order,
        result_folder
    )


    # ========================================================
    # STEP 3-11. RESULT SAVE
    # ========================================================

    result_excel_path = save_timecourse_results(
        df,
        biological_summary,
        timecourse_summary,
        timecourse_stat_results,
        result_folder
    )


    # ========================================================
    # STEP 3-12. FINAL VALIDATION
    # ========================================================

    validation_passed = validate_timecourse_results(
        df_valid,
        biological_summary,
        found_timepoints,
        found_samples,
        result_excel_path,
        graph_path
    )

    return validation_passed


# ============================================================
# STEP 3. EXPERIMENT ANALYSIS
# ============================================================

# ============================================================
# SINGLE EXPERIMENT ANALYSIS
# ============================================================

def analyze_experiment(experiment_folder):

    print("\n\n" + "#" * 60)

    print(f"# START ANALYSIS: {experiment_folder.name}")

    print("#" * 60)

    # ========================================================
    # RESULT FOLDER
    # ========================================================

    result_folder = experiment_folder / "Final_Results"

    result_folder.mkdir(
        parents=True,
        exist_ok=True
    )

    print(
        f"\n결과 저장 폴더:\n"
        f"{result_folder}"
    )

    # ========================================================
    # STEP 3-1. RTF FILE SEARCH
    # ========================================================

    print("\n" + "=" * 60)
    print("STEP 3-1. RTF FILE SEARCH")
    print("=" * 60)

    rtf_files = find_rtf_files(experiment_folder)

    print(f"\n{experiment_folder.name}")

    print(f"RTF 파일 수: {len(rtf_files)}")

    for i, file in enumerate(rtf_files, start=1):

        print(f"  RTF{i}: {file.name}")

    # ========================================================
    # STEP 3-2. RTF → gMFI EXTRACTION
    # ========================================================

    print("\n" + "=" * 60)
    print("STEP 3-2. RTF → gMFI EXTRACTION")
    print("=" * 60)

    df = extract_experiment_data(
        rtf_files,
        experiment_folder.name
    )

    # ========================================================
    # ANALYSIS
    # ========================================================

    if ANALYSIS_TYPE == "standard":

        validation_passed = analyze_standard_experiment(
            df,
            result_folder
        )

    elif ANALYSIS_TYPE == "timecourse":

        validation_passed = analyze_timecourse_experiment(
            df,
            result_folder
        )

    # ========================================================
    # FINAL ANALYSIS STATUS
    # ========================================================

    print("\n" + "=" * 60)

    if validation_passed:

        print("FINAL VALIDATION: PASS")

        print(
            f"{experiment_folder.name} ANALYSIS COMPLETE"
        )

    else:

        print("FINAL VALIDATION: FAIL")

        print(
            f"{experiment_folder.name} ANALYSIS FAILED"
        )

    print("=" * 60)

    return validation_passed


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
        f"\n찾은 Experiment 수: "
        f"{len(all_experiment_folders)}"
    )

    for folder in all_experiment_folders:

        print(f" - {folder.name}")

    
    # ============================================================
    # STEP 2. COMPLETED / NEW EXPERIMENT CHECK
    # ============================================================

    print("\n" + "=" * 60)
    print("STEP 2. ANALYSIS STATUS")
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

        print("\n새롭게 분석할 Experiment가 없습니다.")

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