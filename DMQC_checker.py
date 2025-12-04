import os
import sys
import re
from urllib.parse import urljoin
from urllib.request import urlopen, urlretrieve
import glob
import shutil
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
# Argo QC flag mapping – extend as needed
QC_FLAGS = {
    "good": "1",
    "probably_good": "2",
    "probably_bad": "3",
    "bad": "4",
    "changed": "6",
    "missing": "9",
}


def download_r_files(float_id, out_dir,
                     dac="coriolis",
                     base_root="https://data-argo.ifremer.fr/dac"):
    """
    Downloads all R<float_id>_*.nc files for the given float
    into out_dir.

    Example:
      download_r_files("2904026", "./data/")
    """

    os.makedirs(out_dir, exist_ok=True)

    base_url = f"{base_root}/{dac}/{float_id}/profiles/"
    print(f"\nChecking: {base_url}")

    # Fetch HTML directory listing
    with urlopen(base_url) as resp:
        html = resp.read().decode("utf-8", errors="ignore")

    # extract all hrefs
    hrefs = re.findall(r'href="([^"]+)"', html)

    # is file name like R2904026_004.nc
    pattern = f"R{float_id}_"
    rfiles = sorted(
        fname for fname in hrefs
        if fname.startswith(pattern) and fname.endswith(".nc")
    )

    if not rfiles:
        print("No R-files found — maybe float ID incorrect?")
        return

    print(f"Found {len(rfiles)} R-files to download.\n")

    for fname in rfiles:
        fpath = os.path.join(out_dir, fname)
        if os.path.exists(fpath):
            print(f"File exists — skipping: {fname}")
            continue

        file_url = urljoin(base_url, fname)
        print(f"Downloading: {fname}")
        urlretrieve(file_url, fpath)

    print(f"\nDone. Files stored in: {out_dir}\n")


def fill_qc_var(ds, varname, entry):
    """
    Fill a QC variable in ds (if present) with the qc_flag in entry.
    Supports char ('S1') and numeric QC vars.

    entry must include: entry["qc_flag"] = "1" / "4" / ...
    """
    if varname not in ds.variables:
        return

    qc_flag_char = entry["qc_flag"]        # e.g. "1"
    v = ds.variables[varname]
    shape = v.shape
    kind = v.dtype.kind                    # 'S', 'i', 'u', etc.

    if kind == "S":
        arr = np.full(shape, qc_flag_char.encode("ascii"), dtype="S1")
    elif kind in ("i", "u"):
        arr = np.full(shape, int(qc_flag_char), dtype=v.dtype)
    else:
        print(f"[warn] QC var {varname} has unexpected dtype {v.dtype}, skipping.")
        return

    v[...] = arr


def flag_qc_by_user(float_id, out_dir, default_qc_input="good"):
    """
    Scan out_dir for R<float_id>_*.nc files and return a list
    of dicts, one per profile, each with a profile-level QC flag.

    For now, every profile gets the same QC flag (default_qc).
    Later you can replace default_qc with results from your GUI/logic.
    """
    default_qc = QC_FLAGS[default_qc_input]
    pattern = os.path.join(out_dir, f"R{float_id}_*.nc")
    rfiles = sorted(glob.glob(pattern))

    profile_qc_list = []
    for fpath in rfiles:
        fname = os.path.basename(fpath)
        # e.g. R6903708_012.nc -> "012"
        cycle = fname.split("_")[1].split(".")[0]

        profile_qc_list.append({
            "file": fpath,
            "cycle": cycle,
            "qc_flag": default_qc,  # e.g. "1" = good
        })

    return profile_qc_list

def write_d_files_from_qc_list(profile_qc_list):
    """
    For each entry in profile_qc_list (with keys: 'file', 'cycle', 'qc_flag'),
    create a D-file from the corresponding R-file and apply the profile-level QC
    to PRES/TEMP/PSAL QC variables and DATA_MODE.

    This version DOES NOT change any data values, only QC flags + DATA_MODE + DATE_UPDATE.
    """
    d_files = []
    for entry in profile_qc_list:
        r_file = entry["file"]
        qc_flag = entry["qc_flag"]  # already a single char like '1', '4', etc.

        dirname, basename = os.path.split(r_file)
        # R6903708_001.nc -> D6903708_001.nc
        d_basename = "D" + basename[1:]
        d_file = os.path.join(dirname, d_basename)
        d_files.append(d_file)
        
        # Copy R -> D
        shutil.copy(r_file, d_file)
        print(f"Created {d_file} from {basename}")

        # Open D-file for in-place editing
        ds = Dataset(d_file, "r+")
        try:
            # 1) Set DATA_MODE for the profile to 'D'
            if "DATA_MODE" in ds.variables:
                dm = ds.variables["DATA_MODE"]
                dm[0] = b"D" # Only adjust the main profile.


            # 2) Apply QC flag to raw QC vars (whole profile)
            for vname in ["PRES_QC", "TEMP_QC", "PSAL_QC"]:
                fill_qc_var(ds,vname,entry)

            # 3) If adjusted QC fields exist, also set them
            for vname in ["PRES_ADJUSTED_QC", "TEMP_ADJUSTED_QC", "PSAL_ADJUSTED_QC"]:
                fill_qc_var(ds,vname,entry)

            # 4) Update global DATE_UPDATE
            ds.setncattr("DATE_UPDATE", datetime.utcnow().strftime("%Y%m%d%H%M%S"))

        finally:
            ds.close()

        print(f"  -> updated QC to '{qc_flag}' and set DATA_MODE='D'")
    return d_files

def print_qc_summary(nc_path):
    """
    Print a compact QC summary for PRES/TEMP/PSAL
    and a short per-profile summary (type, levels, depth, time).
    """

    print(f"\nQC SUMMARY FOR: {os.path.basename(nc_path)}")
    ds = Dataset(nc_path)

    def show(var):
        if var in ds.variables:
            data = ds.variables[var][:]
            uniq = np.unique(data)
            print(f"  {var:<16} dtype={ds.variables[var].dtype}  unique={uniq}")
        else:
            print(f"  {var:<16} MISSING")

    # QC fields
    for v in ["PRES_QC", "TEMP_QC", "PSAL_QC",
              "PRES_ADJUSTED_QC", "TEMP_ADJUSTED_QC", "PSAL_ADJUSTED_QC"]:
        show(v)

    # DATA_MODE
    if "DATA_MODE" in ds.variables:
        print(f"  DATA_MODE        {ds.variables['DATA_MODE'][:]}")

    # Per-profile summary
    if "N_PROF" in ds.dimensions:
        n_prof = ds.dimensions["N_PROF"].size
        print(f"\n  N_PROF = {n_prof}")

        # Try to get VERTICAL_SAMPLING_SCHEME & JULD
        vss_var = ds.variables.get("VERTICAL_SAMPLING_SCHEME")
        juld_var = ds.variables.get("JULD")

        # Decode JULD if possible
        juld_dates = None
        if juld_var is not None:
            try:
                from netCDF4 import num2date
                juld_dates = num2date(juld_var[:], juld_var.units)
            except Exception:
                juld_dates = None

        for ip in range(n_prof):
            pres = ds.variables["PRES"][ip, :]
            nlev = int(np.sum(~np.isnan(pres)))
            maxz = float(np.nanmax(pres)) if nlev > 0 else np.nan

            # sampling scheme string
            vss_str = ""
            if vss_var is not None:
                try:
                    vss_bytes = vss_var[ip, :].tobytes()
                    vss_str = vss_bytes.decode("ascii", errors="ignore").strip()
                except Exception:
                    vss_str = ""

            # time string
            t_str = ""
            if juld_dates is not None:
                try:
                    t_str = juld_dates[ip].isoformat()
                except Exception:
                    t_str = ""

            print(f"  PROF {ip:2d}: levels={nlev:3d}, maxP={maxz:6.1f} dbar,"
                  f" VSS='{vss_str}', time={t_str}")

    ds.close()
    print()


if __name__ == "__main__":
    DO_DOWNLOAD = False  #testing purposes, if the download has been made and no need to redownload
    float_id ='6903708'
    out_dir_main = r'C:\Data\ARGO_Dataa\DMQCprocessing\\'
    out_dir = os.path.join(out_dir_main, float_id)
    #First download all files for the float in question:
    if( DO_DOWNLOAD):
        download_r_files(float_id, out_dir)
    #next let's go visual comparison, to mark wether they are okay:
    profile_qc_list = flag_qc_by_user(float_id, out_dir)        
    #All quality control steps palnend are made. Let's write the D-files
    d_files = write_d_files_from_qc_list(profile_qc_list)
    
    # Print QC summaries for the first few D-files
    for dfile in d_files[:3]:  # show first 3
        if os.path.exists(dfile):
            print_qc_summary(dfile)

