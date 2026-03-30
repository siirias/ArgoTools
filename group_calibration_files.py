# -*- coding: utf-8 -*-

import os
import re
import shutil
import pdfplumber


def extract_text_from_pdf(pdf_path):
    lines = []
    with pdfplumber.open(pdf_path) as pdf:
        for page in pdf.pages:
            words = page.extract_words(x_tolerance=1, y_tolerance=1)
            rows = []
            for w in words:
                x, y, text = w["x0"], w["top"], w["text"]
                for row in rows:
                    if abs(row["y"] - y) <= 2:
                        row["items"].append((x, text))
                        break
                else:
                    rows.append({"y": y, "items": [(x, text)]})

            rows.sort(key=lambda r: r["y"])
            for row in rows:
                row["items"].sort(key=lambda t: t[0])
                lines.append(" ".join(text for x, text in row["items"]))
    return "\n".join(lines)


def extract_sensor_number(pdf_path):
    text = extract_text_from_pdf(pdf_path)

    parts = text.split("SENSOR SERIAL NUMBER")
    block = None
    for part in parts:
        if "CONDUCTIVITY CALIBRATION DATA" in part:
            block = "SENSOR SERIAL NUMBER" + part
            break

    if block is None:
        return None

    end_marker = "Date, Slope Correction"
    if end_marker in block:
        block = block[:block.rfind(end_marker)]

    m = re.search(r"SERIAL NUMBER\s*:?\s*(\d+)", block)
    return m.group(1) if m else None


def group_pdfs_by_sensor(directory):
    sensors = {}
    failed = []

    for fname in os.listdir(directory):
        if not fname.lower().endswith(".pdf"):
            continue

        path = os.path.join(directory, fname)
        try:
            sensor = extract_sensor_number(path)
            if not sensor:
                failed.append(fname)
                continue

            sensors.setdefault(sensor, []).append(fname)

            outdir = os.path.join(directory, f"sensor_{sensor}")
            os.makedirs(outdir, exist_ok=True)
            shutil.copy2(path, os.path.join(outdir, fname))

        except Exception as e:
            failed.append(fname)
            print(f"Failed: {fname} ({e})")

    print("Sensors found:")
    for sensor in sorted(sensors, key=int):
        print(f"{sensor}: {len(sensors[sensor])} file(s)")
        for fname in sensors[sensor]:
            print(f"  {fname}")

    if failed:
        print("\nFailed files:")
        for fname in failed:
            print(f"  {fname}")

    return sensors


if __name__ == "__main__":
    group_pdfs_by_sensor(r"C:/data/argodata/calib_sheet_test/morepdf/")