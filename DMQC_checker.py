import os
import re
import glob
import shutil
from urllib.parse import urljoin
from urllib.request import urlopen, urlretrieve
from datetime import datetime

import numpy as np
from netCDF4 import Dataset
import yaml


QC_FLAGS = {
    "good": "1",
    "probably_good": "2",
    "probably_bad": "3",
    "bad": "4",
    "changed": "6",
    "missing": "9",
}

DMQC_META_DEFAULT = {
    "institution": "IF",  # FMI
    "operator": "unknown",
    "software": "argo_helper",
    "software_release": "0.1",
    "history_action": "QC",
    "history_step": "DMQC",
    "comment": "DMQC performed; no value adjustment applied; QC/ERROR updated.",
    "parameter_data_mode": {"PRES": "D", "TEMP": "D", "PSAL": "D"},
    "calib_equation": "none",
    "calib_coefficient": "none",
}


# ---------- small helpers ----------

def _put_char_row(v, irow, text):
    """Write padded ASCII string into v[irow, :] where v is (N_ROW, N_CHAR) of S1."""
    nchar = v.shape[-1]
    s = (text + " " * nchar)[:nchar].encode("ascii", "ignore")
    v[irow, :] = np.frombuffer(s, dtype="S1")


def _put_char_3d(v, iprof, irow, text):
    """Write padded ASCII string into v[iprof, irow, :] where v is (N_PROF, N_ROW, N_CHAR) of S1."""
    nchar = v.shape[-1]
    s = (text + " " * nchar)[:nchar].encode("ascii", "ignore")
    v[iprof, irow, :] = np.frombuffer(s, dtype="S1")


def _get_profile_data(var, iprof):
    """Return (data_1d, ndim) for profile iprof from var; support 2D (N_PROF, N_LEVELS) and 1D (N_LEVELS)."""
    if var.ndim == 2:
        if var.shape[0] <= iprof:
            return None, var.ndim
        return var[iprof, :], var.ndim
    return var[:], var.ndim


def _set_profile_data(var, iprof, data):
    """Assign 1D array back into var at profile iprof (2D) or whole var (1D)."""
    if var.ndim == 2:
        var[iprof, :] = data
    else:
        var[:] = data


def get_meta(entry):
    """
    Return a safe meta dict for this entry:
    - copy outer dict
    - copy nested parameter_data_mode to avoid shared-mutation surprises
    """
    base = entry.get("meta", DMQC_META_DEFAULT)
    meta = dict(base)
    if "parameter_data_mode" in meta and isinstance(meta["parameter_data_mode"], dict):
        meta["parameter_data_mode"] = dict(meta["parameter_data_mode"])
    return meta

def ensure_float_layout(float_dir):
    r_dir = os.path.join(float_dir, "R")
    d_dir = os.path.join(float_dir, "D")
    cycles_dir = os.path.join(float_dir, "cycles")
    os.makedirs(r_dir, exist_ok=True)
    os.makedirs(d_dir, exist_ok=True)
    os.makedirs(cycles_dir, exist_ok=True)
    return r_dir, d_dir, cycles_dir


def load_or_init_meta(meta_path, float_id, dac="coriolis"):
    """Load meta.yaml; if missing, create from DMQC_META_DEFAULT."""
    if os.path.exists(meta_path):
        with open(meta_path, "r", encoding="utf-8") as f:
            meta = yaml.safe_load(f) or {}
    else:
        meta = dict(DMQC_META_DEFAULT)
        meta["parameter_data_mode"] = dict(meta["parameter_data_mode"])
        meta["float_id"] = str(float_id)
        meta["dac"] = str(dac)
        with open(meta_path, "w", encoding="utf-8") as f:
            yaml.safe_dump(meta, f, sort_keys=False, allow_unicode=True)
    return meta


def cycle_yaml_path(cycles_dir, cycle):
    return os.path.join(cycles_dir, f"{str(cycle).zfill(3)}.yaml")


def load_or_init_cycle_decision(cycles_dir, cycle, qc_default="1"):
    """
    Load cycles/NNN.yaml. If missing, create a stub.
    Return dict with at least keys: cycle, qc_flag, note (optional).
    """
    path = cycle_yaml_path(cycles_dir, cycle)
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as f:
            d = yaml.safe_load(f) or {}
    else:
        d = {"cycle": str(cycle).zfill(3), "qc_flag": str(qc_default), "note": ""}
        with open(path, "w", encoding="utf-8") as f:
            yaml.safe_dump(d, f, sort_keys=False, allow_unicode=True)
    return d


def build_profile_list_from_rfiles(float_id, r_dir, meta, cycles_dir):
    """
    Build your profile_qc_list from local R-files + cycle YAML.
    Each entry: {file, cycle, qc_flag, meta}
    """
    pattern = os.path.join(r_dir, f"R{float_id}_*.nc")
    rfiles = sorted(glob.glob(pattern))

    qc_default = str(meta.get("qc_flag_default", QC_FLAGS["good"]))
    out = []
    for fpath in rfiles:
        fname = os.path.basename(fpath)
        cycle = fname.split("_")[1].split(".")[0]  # "003"
        cyc_dec = load_or_init_cycle_decision(cycles_dir, cycle, qc_default=qc_default)

        entry = {
            "file": fpath,
            "cycle": cycle,
            "qc_flag": str(cyc_dec.get("qc_flag", qc_default)),
            "meta": meta,  # your writer will copy safely via get_meta(entry)
        }
        # (optional) carry note forward for later logging
        if "note" in cyc_dec:
            entry["note"] = cyc_dec["note"]
        out.append(entry)

    return out



# ---------- downloading ----------

def download_r_files(float_id, r_dir, dac="coriolis", base_root="https://data-argo.ifremer.fr/dac"):
    os.makedirs(r_dir, exist_ok=True)
    base_url = f"{base_root}/{dac}/{float_id}/profiles/"
    print(f"\nChecking: {base_url}")

    with urlopen(base_url) as resp:
        html = resp.read().decode("utf-8", errors="ignore")

    hrefs = re.findall(r'href="([^"]+)"', html)
    pattern = f"R{float_id}_"
    rfiles = sorted(fname for fname in hrefs if fname.startswith(pattern) and fname.endswith(".nc"))

    if not rfiles:
        print("No R-files found — maybe float ID incorrect?")
        return

    print(f"Found {len(rfiles)} R-files to download.\n")

    for fname in rfiles:
        fpath = os.path.join(r_dir, fname)
        if os.path.exists(fpath):
            print(f"File exists — skipping: {fname}")
            continue
        urlretrieve(urljoin(base_url, fname), fpath)
        print(f"Downloaded: {fname}")

    print(f"\nDone. Files stored in: {r_dir}\n")


# ---------- “user flags” placeholder ----------

def flag_qc_by_user(float_id, out_dir, default_qc_input="good", meta=None):
    """
    For now: every profile gets same qc flag + same meta.
    Later: replace with GUI results and per-level structures.
    """
    default_qc = QC_FLAGS[default_qc_input]
    pattern = os.path.join(out_dir, f"R{float_id}_*.nc")
    rfiles = sorted(glob.glob(pattern))

    if meta is None:
        meta = DMQC_META_DEFAULT

    out = []
    for fpath in rfiles:
        fname = os.path.basename(fpath)
        cycle = fname.split("_")[1].split(".")[0]
        out.append({"file": fpath, "cycle": cycle, "qc_flag": default_qc, "meta": meta})
    return out


# ---------- paperwork writers ----------

def set_parameter_data_mode(ds, iprof, pdm_map):
    """Update PARAMETER_DATA_MODE for selected parameters at profile iprof."""
    if "PARAMETER" not in ds.variables or "PARAMETER_DATA_MODE" not in ds.variables:
        return

    v_par = ds.variables["PARAMETER"]            # usually (N_PROF, N_PARAM, N_CHAR)
    v_pdm = ds.variables["PARAMETER_DATA_MODE"]  # usually (N_PROF, N_PARAM) S1
    if v_par.ndim != 3 or v_pdm.ndim != 2 or iprof >= v_par.shape[0]:
        return

    for j in range(v_par.shape[1]):
        name = v_par[iprof, j, :].tobytes().decode("ascii", "ignore").strip()
        if name in pdm_map:
            v_pdm[iprof, j] = pdm_map[name].encode("ascii")


def append_history(ds, meta, action=None, step=None, parameter=None, comment=None):
    """
    Append one HISTORY row if HISTORY vars exist.
    Uses meta defaults if action/step/parameter/comment not provided.
    """
    if "N_HISTORY" not in ds.dimensions:
        return False

    v_date = ds.variables.get("HISTORY_DATE")
    if v_date is None or v_date.ndim != 2:
        return False

    i_free = None
    for i in range(v_date.shape[0]):
        if v_date[i, :].tobytes().strip() == b"":
            i_free = i
            break
    if i_free is None:
        return False

    if action is None:
        action = meta.get("history_action", "QC")
    if step is None:
        step = meta.get("history_step", "DMQC")
    if parameter is None:
        pdm = meta.get("parameter_data_mode", {})
        parameter = " ".join(pdm.keys()) if pdm else "PRES TEMP PSAL"
    if comment is None:
        comment = meta.get("comment", "")

    now = datetime.utcnow().strftime("%Y%m%d%H%M%S")

    fields = {
        "HISTORY_DATE": now,
        "HISTORY_ACTION": action,
        "HISTORY_STEP": step,
        "HISTORY_PARAMETER": parameter,
        "HISTORY_INSTITUTION": meta.get("institution", ""),
        "HISTORY_SOFTWARE": meta.get("software", ""),
        "HISTORY_SOFTWARE_RELEASE": meta.get("software_release", ""),
        "HISTORY_COMMENT": comment,  # optional in some templates
    }

    for varname, text in fields.items():
        v = ds.variables.get(varname)
        if v is None or v.ndim != 2:
            continue
        _put_char_row(v, i_free, text)

    return True


def append_scientific_calib(ds, meta, iprof=0, params=("PRES", "TEMP", "PSAL")):
    """
    Append SCIENTIFIC_CALIB rows (one per param) if SCIENTIFIC_CALIB vars exist.
    Minimal “no adjustment” fill: equation/coeff = meta values, comment = meta comment.
    """
    if "N_CALIB" not in ds.dimensions:
        return

    needed = [
        "SCIENTIFIC_CALIB_DATE", "SCIENTIFIC_CALIB_PARAMETER",
        "SCIENTIFIC_CALIB_EQUATION", "SCIENTIFIC_CALIB_COEFFICIENT", "SCIENTIFIC_CALIB_COMMENT",
    ]
    if not all(v in ds.variables for v in needed):
        return

    v_date = ds.variables["SCIENTIFIC_CALIB_DATE"]  # (N_PROF, N_CALIB, N_CHAR)
    if v_date.ndim != 3 or iprof >= v_date.shape[0]:
        return

    # first empty row for this profile
    i_free = None
    for i in range(v_date.shape[1]):
        if v_date[iprof, i, :].tobytes().strip() == b"":
            i_free = i
            break
    if i_free is None:
        return

    now = datetime.utcnow().strftime("%Y%m%d%H%M%S")
    eq = meta.get("calib_equation", "none")
    co = meta.get("calib_coefficient", "none")
    cm = meta.get("comment", "DMQC performed.")

    for p in params:
        if i_free >= v_date.shape[1]:
            break
        _put_char_3d(ds.variables["SCIENTIFIC_CALIB_DATE"], iprof, i_free, now)
        _put_char_3d(ds.variables["SCIENTIFIC_CALIB_PARAMETER"], iprof, i_free, p)
        _put_char_3d(ds.variables["SCIENTIFIC_CALIB_EQUATION"], iprof, i_free, eq)
        _put_char_3d(ds.variables["SCIENTIFIC_CALIB_COEFFICIENT"], iprof, i_free, co)
        _put_char_3d(ds.variables["SCIENTIFIC_CALIB_COMMENT"], iprof, i_free, cm)
        i_free += 1


# ---------- D-file writer ----------

def write_d_files_from_qc_list(profile_qc_list, d_dir, default_errors=None):
    # Reasonable placeholders (tune later)
    if default_errors is None:
        default_errors = {"PRES": 0.1, "TEMP": 0.002, "PSAL": 0.01}

    d_files = []
    iprof_main = 0

    for entry in profile_qc_list:
        meta = get_meta(entry)
        r_file = entry["file"]
        qc_flag = entry["qc_flag"]

        basename = os.path.basename(r_file)
        d_file = os.path.join(d_dir, "D" + basename[1:])
        d_files.append(d_file)

        shutil.copy(r_file, d_file)
        print(f"Created {d_file} from {basename}")

        ds = Dataset(d_file, "r+")
        try:
            # 1) main profile DATA_MODE -> D
            if "DATA_MODE" in ds.variables:
                dm = ds.variables["DATA_MODE"]
                if dm.size > iprof_main:
                    dm[iprof_main] = b"D"

            # 2) main profile: set *_ADJUSTED_QC + *_ADJUSTED_ERROR
            for param in ["PRES", "TEMP", "PSAL"]:
                adj_name = f"{param}_ADJUSTED"
                adj_qc_name = f"{param}_ADJUSTED_QC"
                adj_err_name = f"{param}_ADJUSTED_ERROR"
                if adj_name not in ds.variables:
                    continue

                v_adj = ds.variables[adj_name]
                adj_data, _ = _get_profile_data(v_adj, iprof_main)
                if adj_data is None:
                    continue

                fill = getattr(v_adj, "_FillValue", None)
                mask_non_missing = (adj_data != fill) if fill is not None else ~np.isnan(adj_data)

                if adj_qc_name in ds.variables:
                    v_qc = ds.variables[adj_qc_name]
                    qc_arr, _ = _get_profile_data(v_qc, iprof_main)
                    qc_arr[:] = b"9"
                    qc_arr[mask_non_missing] = qc_flag.encode("ascii")
                    _set_profile_data(v_qc, iprof_main, qc_arr)

                if adj_err_name in ds.variables:
                    v_err = ds.variables[adj_err_name]
                    err_arr, _ = _get_profile_data(v_err, iprof_main)
                    err_arr[:] = fill if fill is not None else np.nan
                    err_arr[mask_non_missing] = default_errors.get(param, 0.0)
                    _set_profile_data(v_err, iprof_main, err_arr)

            # 3) fix blank *_ADJUSTED_QC in A-mode profiles
            if "DATA_MODE" in ds.variables:
                dm = ds.variables["DATA_MODE"]
                for ip in range(dm.size):
                    if dm[ip] != b"A":
                        continue
                    for param in ["PRES", "TEMP", "PSAL"]:
                        adj_name = f"{param}_ADJUSTED"
                        adj_qc_name = f"{param}_ADJUSTED_QC"
                        if adj_name not in ds.variables or adj_qc_name not in ds.variables:
                            continue

                        v_adj = ds.variables[adj_name]
                        v_qc = ds.variables[adj_qc_name]
                        adj_data, _ = _get_profile_data(v_adj, ip)
                        qc_arr, _ = _get_profile_data(v_qc, ip)
                        if adj_data is None or qc_arr is None:
                            continue

                        fill = getattr(v_adj, "_FillValue", None)
                        mask_non_missing = (adj_data != fill) if fill is not None else ~np.isnan(adj_data)
                        fix_mask = mask_non_missing & (qc_arr == b" ")
                        qc_arr[fix_mask] = b"1"
                        _set_profile_data(v_qc, ip, qc_arr)

            # 4) update DATE_UPDATE + DATA_STATE_INDICATOR (attr + first entry in var)
            ds.setncattr("DATE_UPDATE", datetime.utcnow().strftime("%Y%m%d%H%M%S"))
            ds.setncattr("DATA_STATE_INDICATOR", "2C")

            if "DATA_STATE_INDICATOR" in ds.variables:
                v_dsi = ds.variables["DATA_STATE_INDICATOR"]
                if v_dsi.dtype.kind == "S":
                    nchar = v_dsi.shape[-1] if len(v_dsi.shape) > 0 else 0
                    if nchar > 0:
                        row0 = np.frombuffer(((b"2C" + b" " * nchar)[:nchar]), dtype="S1")
                        if v_dsi.ndim == 1:
                            v_dsi[:] = row0
                        elif v_dsi.ndim == 2 and v_dsi.shape[0] >= 1:
                            arr = v_dsi[:]
                            arr[0, :] = row0
                            v_dsi[:] = arr

            # 5) DMQC paperwork (only if vars exist)
            set_parameter_data_mode(ds, iprof_main, meta.get("parameter_data_mode", {}))
            append_history(ds, meta)  # uses meta defaults
            append_scientific_calib(
                ds, meta, iprof=iprof_main,
                params=tuple(meta.get("parameter_data_mode", {}).keys()) or ("PRES", "TEMP", "PSAL")
            )

        finally:
            ds.close()

        print(f"  -> updated QC '{qc_flag}', DATA_MODE[0]='D', DATA_STATE_INDICATOR='2C' (+ paperwork if possible)")

    return d_files


# ---------- quick check ----------

def print_qc_summary(nc_path):
    print(f"\nQC SUMMARY FOR: {os.path.basename(nc_path)}")
    ds = Dataset(nc_path)

    def show(var):
        if var in ds.variables:
            uniq = np.unique(ds.variables[var][:])
            print(f"  {var:<16} dtype={ds.variables[var].dtype}  unique={uniq}")
        else:
            print(f"  {var:<16} MISSING")

    for v in ["PRES_ADJUSTED_QC", "TEMP_ADJUSTED_QC", "PSAL_ADJUSTED_QC",
              "PRES_ADJUSTED_ERROR", "TEMP_ADJUSTED_ERROR", "PSAL_ADJUSTED_ERROR"]:
        show(v)

    if "DATA_MODE" in ds.variables:
        print(f"  DATA_MODE        {ds.variables['DATA_MODE'][:]}")

    ds.close()
    print()


if __name__ == "__main__":
    DO_DOWNLOAD = True
    float_id = "6903708"
    out_dir_main = r"C:\Data\ARGO_Dataa\DMQCprocessing"
    float_dir = os.path.join(out_dir_main, float_id)

    r_dir, d_dir, cycles_dir = ensure_float_layout(float_dir)

    # meta.yaml
    meta = load_or_init_meta(os.path.join(float_dir, "meta.yaml"), float_id, dac="coriolis")

    if DO_DOWNLOAD:
        download_r_files(float_id, r_dir, dac=meta.get("dac", "coriolis"))

    # Build entries from R-files + per-cycle yaml
    profile_qc_list = build_profile_list_from_rfiles(float_id, r_dir, meta, cycles_dir)

    # Write D-files into D/
    d_files = write_d_files_from_qc_list(profile_qc_list, d_dir=d_dir)

    for dfile in d_files[:3]:
        if os.path.exists(dfile):
            print_qc_summary(dfile)
