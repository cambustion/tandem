# tandem_cpma_utils.py

import csv
from pathlib import Path
from typing import Tuple, Optional
import pandas as pd
import numpy as np
import logging


# -----------------------------------------------------------------------------
# LOGGING SETUP
# -----------------------------------------------------------------------------
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Only add handler if not already added (avoids duplicate logs on import)
if not logger.hasHandlers():
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)s - %(name)s - %(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

# -----------------------------------------------------------------------------
# CONSTANTS
# -----------------------------------------------------------------------------
LPM_TO_M3S = 1.0 / 60_000.0  # L/min -> m^3/s
P_REF = 101_325.0             # Reference pressure (Pa)
T_K_OFFSET = 273.15           # °C -> K


# -----------------------------------------------------------------------------
# FUNCTION DEFINITIONS
# -----------------------------------------------------------------------------
def read_raw_txt(txt_path: Path) -> pd.DataFrame:
    """
    Read the Cambustion tandem TXT file into a raw DataFrame.

    Parameters
    ----------
    txt_path : Path
        Path to the TXT file.

    Returns
    -------
    pd.DataFrame
        Raw DataFrame as read from the TXT file.
    """
    logger.info("Reading TXT file: %s", txt_path)

    with open(
        txt_path,
        "r",
        newline="",
        encoding="utf-8-sig"
    ) as f:
        content = f.read()

    if not content.strip():
        raise ValueError(f"TXT file is empty: {txt_path}")

    # -------------------------------------------------------------------------
    # Detect delimiter
    # -------------------------------------------------------------------------
    try:
        dialect = csv.Sniffer().sniff(
            content[:10000],
            delimiters=",\t"
        )
        delimiter = dialect.delimiter

    except csv.Error:
        # Fallback: compare the first several lines
        sample_lines = content.splitlines()[:10]

        comma_count = sum(
            line.count(",")
            for line in sample_lines
        )

        tab_count = sum(
            line.count("\t")
            for line in sample_lines
        )

        delimiter = "\t" if tab_count > comma_count else ","

    logger.info(
        "Detected TXT delimiter: %s",
        "TAB" if delimiter == "\t" else "COMMA"
    )

    # -------------------------------------------------------------------------
    # Parse CSV/TXT records
    # -------------------------------------------------------------------------
    reader = csv.reader(
        content.splitlines(),
        delimiter=delimiter,
        quotechar='"'
    )

    rows = list(reader)

    if not rows:
        raise ValueError(f"No data found in TXT file: {txt_path}")

    # -------------------------------------------------------------------------
    # Pad rows to equal length
    # -------------------------------------------------------------------------
    n_columns = max(len(row) for row in rows)

    rows_padded = [
        row + [None] * (n_columns - len(row))
        for row in rows
    ]

    return pd.DataFrame(rows_padded)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------
def _clean_value(value) -> str:
    """Convert a raw field to a clean string."""

    if value is None or pd.isna(value):
        return ""

    return (
        str(value)
        .replace("\ufeff", "")
        .replace("\r", "")
        .replace("\n", " ")
        .strip()
    )


def _find_row(df_raw: pd.DataFrame, first_cell: str) -> Optional[int]:
    """
    Find the first row whose first non-empty cell matches `first_cell`.

    Matching is case-insensitive and whitespace/BOM tolerant.
    """

    target = (
        str(first_cell)
        .replace("\ufeff", "")
        .strip()
        .lower()
    )

    for i in range(len(df_raw)):

        for j in range(df_raw.shape[1]):

            value = (
                _clean_value(df_raw.iloc[i, j])
                .replace("\ufeff", "")
                .strip()
                .lower()
            )

            if value == target:
                return i

    return None


def _find_row_containing(
    df_raw: pd.DataFrame,
    text: str
) -> Optional[int]:
    """
    Find the first row containing a specific text value.
    """

    target = text.strip().lower()

    for i in range(len(df_raw)):
        for j in range(df_raw.shape[1]):
            value = _clean_value(df_raw.iloc[i, j])

            if target in value.lower():
                return i

    return None


def _row_to_dict(header_row, data_row):
    """
    Convert a header row and data row into a dictionary.

    Only fields with non-empty header names are included.
    """

    result = {}

    for header, value in zip(header_row, data_row):

        header = _clean_value(header)

        if not header:
            continue

        result[header] = _clean_value(value)

    return result


def _to_numeric(value):
    """Safely convert a value to numeric."""

    return pd.to_numeric(value, errors="coerce")


def parse_metadata(
    df_raw: pd.DataFrame
) -> Tuple[dict, pd.DataFrame]:
    """
    Extract metadata from top rows and return metadata dict and the cleaned DataFrame.
    
    Parameters
    ----------
    df_raw : pd.DataFrame
        Raw DataFrame read from TXT.
    meta_rows : int
        Number of metadata rows at the top of TXT.
        
    Returns
    -------
    metadata : dict
        Extracted metadata (classifier names, flows, etc.).
    df : pd.DataFrame
        Cleaned DataFrame starting from measurement rows.
    """
    logger.info("Parsing TXT metadata")

    metadata = {}

    # -------------------------------------------------------------------------
    # 1. General file information
    # -------------------------------------------------------------------------
    if len(df_raw) > 0:

        general_row = [
            _clean_value(v)
            for v in df_raw.iloc[0].tolist()
        ]

        if len(general_row) >= 1:
            metadata["manufacturer"] = general_row[0]

        if len(general_row) >= 2:
            metadata["measurement_type"] = general_row[1]

        if len(general_row) >= 3:
            metadata["software_version"] = general_row[2]

        if len(general_row) >= 4:
            metadata["date"] = general_row[3]

        if len(general_row) >= 5:
            metadata["time"] = general_row[4]

        # The last two fields in the current format are typically:
        # "Bypass scans:" and True/False
        if len(general_row) >= 7:
            metadata["bypass_scans"] = (
                general_row[6].lower() == "true"
            )

    # -------------------------------------------------------------------------
    # 2. Locate Classifier 1 header
    # -------------------------------------------------------------------------
    cls1_header_idx = _find_row(
        df_raw,
        "Classifier 1"
    )

    if cls1_header_idx is None:
        raise ValueError(
            "Could not find 'Classifier 1' section in TXT file."
        )

    cls1_data_idx = cls1_header_idx + 1

    cls1_headers = [
        _clean_value(v)
        for v in df_raw.iloc[cls1_header_idx].tolist()
    ]

    cls1_values = [
        _clean_value(v)
        for v in df_raw.iloc[cls1_data_idx].tolist()
    ]

    cls1_meta = _row_to_dict(
        cls1_headers,
        cls1_values
    )

    metadata["classifier_1_info"] = cls1_meta

    # Classifier name
    metadata["classifier_1"] = cls1_meta.get(
        "Classifier 1",
        ""
    ).strip().lower()

    # Sample flow
    metadata["qa_1"] = (
        _to_numeric(
            cls1_meta.get("Sample flow (lpm)")
        ) * LPM_TO_M3S
    )

    # -------------------------------------------------------------------------
    # 3. Locate Classifier 2 header
    # -------------------------------------------------------------------------
    cls2_header_idx = _find_row(
        df_raw,
        "Classifier 2"
    )

    if cls2_header_idx is None:
        raise ValueError(
            "Could not find 'Classifier 2' section in TXT file."
        )

    cls2_data_idx = cls2_header_idx + 1

    cls2_headers = [
        _clean_value(v)
        for v in df_raw.iloc[cls2_header_idx].tolist()
    ]

    cls2_values = [
        _clean_value(v)
        for v in df_raw.iloc[cls2_data_idx].tolist()
    ]

    cls2_meta = _row_to_dict(
        cls2_headers,
        cls2_values
    )

    metadata["classifier_2_info"] = cls2_meta

    metadata["classifier_2"] = cls2_meta.get(
        "Classifier 2",
        ""
    ).strip().lower()

    # Sample flow is common to AAC/DMA
    metadata["qa_2"] = (
        _to_numeric(
            cls2_meta.get("Sample flow (lpm)")
        ) * LPM_TO_M3S
    )

    # Sheath flow may appear in different classifier outputs
    sheath_flow = cls2_meta.get(
        "Sheath flow SP (lpm)",
        cls2_meta.get("Sheath flow (lpm)")
    )

    if sheath_flow is not None:
        metadata["qsh_2"] = (
            _to_numeric(sheath_flow) * LPM_TO_M3S
        )
    else:
        metadata["qsh_2"] = np.nan

    # Other metadata if available
    metadata["serial_number_2"] = cls2_meta.get(
        "Serial number"
    )

    metadata["start_2"] = _to_numeric(
        cls2_meta.get("Start (nm)")
    )

    metadata["end_2"] = _to_numeric(
        cls2_meta.get("End (nm)")
    )

    # -------------------------------------------------------------------------
    # 4. Locate measurement header dynamically
    # -------------------------------------------------------------------------
    measurement_header_idx = None

    measurement_start_columns = {
        "m* (fg)",
        "dm* (nm)",
        "da* (nm)",
    }

    for i in range(len(df_raw)):

        first_value = _clean_value(df_raw.iloc[i, 0]).lower()

        if first_value in measurement_start_columns:
            measurement_header_idx = i
            break

    if measurement_header_idx is None:
        raise ValueError(
            "Could not find the measurement header. "
            "Expected the first column to be one of: "
            "'m* (fg)', 'dm* (nm)', or 'da* (nm)'."
        )

    logger.info(
        "Measurement header found at row %d: %s",
        measurement_header_idx,
        df_raw.iloc[measurement_header_idx, 0]
    )

    # -------------------------------------------------------------------------
    # 5. Build measurement DataFrame
    # -------------------------------------------------------------------------
    measurement_headers = [
        _clean_value(v)
        for v in df_raw.iloc[measurement_header_idx].tolist()
    ]

    measurement_headers = [
        col if col else f"Unnamed_{i}"
        for i, col in enumerate(measurement_headers)
    ]

    measurement_rows = df_raw.iloc[
        measurement_header_idx + 1:
    ].copy()

    # Keep only columns represented in measurement header
    measurement_rows = measurement_rows.iloc[
        :, :len(measurement_headers)
    ]

    measurement_rows.columns = measurement_headers

    # Remove completely empty rows
    measurement_rows = measurement_rows.dropna(
        how="all"
    ).reset_index(drop=True)


    return metadata, measurement_rows


def split_bypassed(df: pd.DataFrame, first_col: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split DataFrame into bypassed and normal rows based on the first column being "Bypassed".
    
    Parameters
    ----------
    df : pd.DataFrame
        Measurement DataFrame.
    first_col : str, optional
        Name of the column to check for "Bypassed" (default is first column).
        
    Returns
    -------
    df_bypassed : pd.DataFrame
        Rows where the first column == "Bypassed".
    df_normal : pd.DataFrame
        All other rows.
    """
    if first_col is None:
        first_col = df.columns[0] if len(df.columns) > 0 else None

    if first_col is None:
        logger.warning("DataFrame has no columns to check for 'Bypassed'.")
        return pd.DataFrame(columns=df.columns), df.copy()

    mask_bypassed = df[first_col].astype(str).str.strip() == "Bypassed"
    df_bypassed = df.loc[mask_bypassed].reset_index(drop=True)
    df_normal = df.loc[~mask_bypassed].reset_index(drop=True)
    logger.info("Split data: %d bypassed rows, %d tandem rows", len(df_bypassed), len(df_normal))

    return df_bypassed, df_normal


def coerce_numeric(df: pd.DataFrame, exclude: Optional[list] = None) -> pd.DataFrame:
    """
    Convert all columns to numeric where possible, optionally excluding some columns.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to convert.
    exclude : list, optional
        List of columns to exclude from numeric conversion.
        
    Returns
    -------
    pd.DataFrame
        DataFrame with numeric columns converted.
    """
    exclude = exclude or []
    cols_to_convert = [c for c in df.columns if c not in exclude]
    df_numeric = df.copy()
    df_numeric[cols_to_convert] = df_numeric[cols_to_convert].apply(pd.to_numeric, errors="coerce")
    return df_numeric


def extract_instrument_data(df: pd.DataFrame, classifier_name: str, is_classifier_2: bool = False) -> dict:
    """
    Dynamically maps columns to standard keys based on the instrument type.
    
        Parameters
    ----------
    df : pd.DataFrame
        DataFrame to convert.
    classifier_name : str
        Name of instrument used.
    is_classifier_2: boolian
        
    Returns
    -------
    data
        Dictionary consists of columns
    """
    data = {}
    if not classifier_name or df.empty:
        return data
        
    c_name = classifier_name.lower()
    suffix = "2" if is_classifier_2 else ""
    
    # --- CPMA Logic ---
    if "cpma" in c_name:
        data["size"] = df.get(f"m* (fg){suffix}").to_numpy(dtype=float) if f"m* (fg){suffix}" in df else np.array([])
        data["speed"] = df.get(f"Speed (rad/s){suffix}").to_numpy(dtype=float) if f"Speed (rad/s){suffix}" in df else np.array([])
        if f"Voltage (V){suffix}" in df:
            data["voltage"] = df.get(f"Voltage (V){suffix}").to_numpy(dtype=float)
        if f"Pressure (Pa){suffix}" in df:
            data["pressure"] = df.get(f"Pressure (Pa){suffix}").to_numpy(dtype=float) / P_REF
        if f"Temperature (C){suffix}" in df:
            data["temperature"] = df.get(f"Temperature (C){suffix}").to_numpy(dtype=float) + T_K_OFFSET
            
    # --- AAC Logic ---
    elif "aac" in c_name:
        data["size"] = df.get(f"da* (nm){suffix}").to_numpy(dtype=float) if f"da* (nm){suffix}" in df else np.array([])
        data["speed"] = df.get(f"Speed (rad/s){suffix}").to_numpy(dtype=float) if f"Speed (rad/s){suffix}" in df else np.array([])
        if f"Sheath flow (lpm){suffix}" in df:
            data["sheath_flow"] = df.get(f"Sheath flow (lpm){suffix}").to_numpy(dtype=float) * LPM_TO_M3S
        if f"Pressure (Pa){suffix}" in df:
            data["pressure"] = df.get(f"Pressure (Pa){suffix}").to_numpy(dtype=float) / P_REF
        if f"Temperature (C){suffix}" in df:
            data["temperature"] = df.get(f"Temperature (C){suffix}").to_numpy(dtype=float) + T_K_OFFSET
            
    # --- DMA Logic ---
    elif "dma" in c_name:
        data["size"] = df.get(f"dm* (nm){suffix}").to_numpy(dtype=float) if f"dm* (nm){suffix}" in df else np.array([])
        
        # DMA pressure is usually reported in kPa. We convert to Pa before dividing by P_REF.
        if f"Pressure (kPa){suffix}" in df:
            data["pressure"] = (df.get(f"Pressure (kPa){suffix}").to_numpy(dtype=float) * 1000.0) / P_REF
        elif f"Pressure (Pa){suffix}" in df:
            data["pressure"] = df.get(f"Pressure (Pa){suffix}").to_numpy(dtype=float) / P_REF
            
        if f"Temperature (C){suffix}" in df:
            data["temperature"] = df.get(f"Temperature (C){suffix}").to_numpy(dtype=float) + T_K_OFFSET

    return data

def main(txt_path_str: str) -> dict:
    """
    Main function to read CSV, extract metadata, clean data, split bypassed rows, and return arrays.
    
    Parameters
    ----------
    txt_path_str : str
        
    Returns
    -------
    dict
        Dictionary containing metadata, DataFrames, setpoints, and numpy arrays.
    """

    txt_path = Path(txt_path_str)

    # -------------------------------------------------------------------------
    # Read and parse file
    # -------------------------------------------------------------------------
    df_raw = read_raw_txt(txt_path)

    metadata, df_clean = parse_metadata(df_raw)

    # -------------------------------------------------------------------------
    # Split bypassed and tandem measurements
    # -------------------------------------------------------------------------
    df_bypass, df = split_bypassed(df_clean)

    # -------------------------------------------------------------------------
    # Convert numeric columns
    # -------------------------------------------------------------------------
    df = coerce_numeric(df)

    if not df_bypass.empty:
        df_bypass = coerce_numeric(df_bypass)

    # -------------------------------------------------------------------------
    # Extract classifier data
    # -------------------------------------------------------------------------
    c1_name = metadata.get("classifier_1", "")
    c2_name = metadata.get("classifier_2", "")

    arrays_c1 = extract_instrument_data(
        df,
        c1_name,
        is_classifier_2=False
    )

    arrays_c2 = extract_instrument_data(
        df,
        c2_name,
        is_classifier_2=True
    )

    arrays_bypass_c1 = extract_instrument_data(
        df_bypass,
        c1_name,
        is_classifier_2=False
    )

    arrays_bypass_c2 = extract_instrument_data(
        df_bypass,
        c2_name,
        is_classifier_2=True
    )

    # -------------------------------------------------------------------------
    # Concentrations
    # -------------------------------------------------------------------------
    conc = (
        df["Conc"].to_numpy(dtype=float)
        if "Conc" in df
        else np.array([])
    )

    bypass_conc = (
        df_bypass["Conc"].to_numpy(dtype=float)
        if "Conc" in df_bypass
        else np.array([])
    )

    # -------------------------------------------------------------------------
    # Setpoints
    # -------------------------------------------------------------------------
    def get_setpoints(data_dict):
        if "size" not in data_dict:
            return []

        size = data_dict["size"]

        if len(size) == 0:
            return []

        valid_sizes = size[~np.isnan(size)]

        return sorted(np.unique(valid_sizes))

    setpnt_1 = get_setpoints(arrays_c1)
    setpnt_2 = get_setpoints(arrays_c2)
    setpnt_bypass = get_setpoints(arrays_bypass_c2)

    logger.info(
        "Main function completed: %d tandem rows, %d bypassed rows",
        len(df),
        len(df_bypass)
    )

    return {
        "metadata": metadata,
        "df_raw": df_raw,
        "df_clean": df_clean,
        "df_bypass": df_bypass,
        "df": df,
        "setpoints": {
            "bypass": setpnt_bypass,
            "1": setpnt_1,
            "2": setpnt_2
        },
        "arrays": {
            "cls1": arrays_c1,
            "cls2": arrays_c2,
            "bypass_cls1": arrays_bypass_c1,
            "bypass_cls2": arrays_bypass_c2,
            "conc": conc,
            "bypass_conc": bypass_conc,
        },
    }



# -----------------------------------------------------------------------------
# Optional test when running this file directly
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    test_txt = Path(
        r"G:\My Drive\PhD Thesis\Tandem_CPMA_SASS_inversion"
        r"\Experimental Data\CPMA_SASS_DOS\your_file.txt"
    )

    try:
        results = main(test_txt)

        logger.info(
            "Test run complete. Extracted Classifiers: %s -> %s",
            results["metadata"]["classifier_1"],
            results["metadata"]["classifier_2"]
        )

    except FileNotFoundError:
        logger.error(
            "Could not find test TXT file. Please verify the path."
        )
