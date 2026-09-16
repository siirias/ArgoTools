# Argo D-file workflow

`DMQC_checker.py` preserves the historical implementation. Use `dmqc_process.py`
for the new framework described below.

`dmqc_process.py` is the entry point. Scientific checkers run independently and
supply decisions; this tool downloads R-files, combines decisions, and writes
D-files. The initial implementation only rejects whole core profiles. It does
not run checks, infer that profiles are good, or calculate corrections.

## Quick start

```bash
python dmqc_process.py
```

No arguments retains float **6903708** and the historical directory:

- Windows: `C:\Data\ARGO_Dataa\DMQCprocessing`
- WSL/Linux default: `/mnt/c/Data/ARGO_Dataa/DMQCprocessing`

The default `all` stage downloads missing R-files, reads available instructions,
and writes files for explicitly rejected profiles. Existing R-files are retained.
With no bad-profile decisions, no D-files are created. Existing D-files require
`--overwrite`; running the command again does not silently replace them.

```bash
python dmqc_process.py download
python dmqc_process.py combine
python dmqc_process.py write --dry-run
python dmqc_process.py write
python dmqc_process.py write --cycles 1,3-5 --overwrite
python dmqc_process.py all --float 6903708 --work-dir /path/to/work --dac coriolis
```

`--cycles` selects cycle numbers, including descending files of those cycles.
`--work-dir` is the parent of the float directory. `--instructions-dir` overrides
the partner instruction directory; legacy `cycles/` decisions are also read.
`--dry-run` writes nothing. A download dry-run still reads the remote listing;
`write --dry-run` is entirely local. An `all --dry-run` can validate only sources
already downloaded; missing source files referenced by instructions are errors.

The `combine` command writes a JSON review plan to stdout and status messages to
stderr, so it can be saved with `> combined.json`. This plan is an audit preview,
not another checker input. `write` recombines the current checker inputs, pins
source hashes, validates them, and writes outputs in the same invocation.

## Components

- `dmqc/download.py`: download sources atomically without overwriting them.
- `dmqc/instructions.py`: internal target/operation objects and provisional YAML adapters.
- `dmqc/combine.py`: merge suggestions, with bad winning over no finding.
- `dmqc/writer.py`: validate, stage, verify, and publish D-files and audit reports.

Python 3.10+ with `numpy`, `netCDF4`, and `PyYAML` is required. No plotting libraries
are needed by this workflow.

## Work directory

```text
DMQCprocessing/
  6903708/
    R/                          original downloaded sources
    D/                          generated D-files and *.report.json audit records
    instructions/               partner YAML files, searched recursively
      partner_a/
      partner_b/
    cycles/                     optional existing per-cycle YAML decisions
    meta.yaml                   optional institution/operator/software information
```

Do not put generated combined plans in `instructions/`. The tool does not create
automatic good decisions or launch partner checkers. Run automatic checkers and
complete human reviews before the writing stage. A present instruction with
`status: pending` blocks processing. This first version has no required-checker
manifest: an absent checker report is not detected as unfinished work.

## Provisional instruction adapter

The final interchange format is still open for discussion. YAML parsing is
isolated from the internal target and operation model so it can be replaced.
The following adapter lets partners exercise the first working path now:

```yaml
schema_version: 1
checker: partner_a.visual_review
instructions:
  - target:
      source: R6903708_001.nc
      profile_index: 0
      selection: whole_profile
    action: flag
    flag: '4'
    status: ready
    reason: 'Entire profile rejected after visual review: explain the finding here.'
    # Optional: SHA-256 of the exact source file used by the checker.
    # source_sha256: '<64 lowercase hexadecimal characters>'
```

`profile_index` is **zero-based**, including in MATLAB-generated YAML. It refers
to one `N_PROF` entry, not necessarily every profile in a file. The writer checks
the filename, `PLATFORM_NUMBER`, and `CYCLE_NUMBER`. Targeting another profile
requires another instruction. Both ascending and descending filenames must be
specified exactly, avoiding ambiguity from a cycle number alone.

Each whole-profile instruction covers `PRES`, `TEMP`, and `PSAL`. An optional
`parameters: [PRES, TEMP, PSAL]` explicitly states the same scope. Other parameter
sets, depth selections, flags, and correction actions are rejected for now.

A checker that finds nothing wrong reports `action: no_finding`, omits `flag`,
and provides a reason. This is **not** a QC upgrade or approval of the profile.
Bad wins independently of checker order; every contributing suggestion remains
in the combined decision and audit report. Missing-data flags remain missing.

Include `source_sha256` when possible. The combiner checks it against the local
source. Without it, the tool cannot tell whether a checker used an older source;
it still binds the combined decision to the current file hash before writing.

### Existing cycle YAML files

Existing `cycles/001.yaml` entries remain readable:

```yaml
cycle: '001'
qc_flag: '4'
note: 'Reason for rejecting the profile'
```

They target profile index 0. Legacy flag `1` means no finding and leaves QC
unchanged; flag `4` requires a nonempty reason. Other flags are rejected instead
of silently interpreted. Missing cycle YAML files are never assigned good QC.
Use new partner documents to address additional profiles or multiple checkers.

### Metadata

Only `institution`, `operator`, `software`, and `software_release` are read from
`meta.yaml`. Defaults are `IF`, `unknown`, `ADMW`, and `0.1`. Set the operator and
software information for real processing. Old mode maps, default errors, and
QC defaults are ignored. NetCDF text fields have fixed widths; full reasons and
metadata are retained in the JSON audit report even when NetCDF text is truncated.

## Writer behavior

For each targeted profile:

- Keep raw measurements and raw QC unchanged.
- Set adjusted QC to `4` for present, non-missing samples, and `9` for missing samples.
- Set rejected adjusted values and adjusted errors to their own fill values.
  Raw data identify rejected samples even if adjusted arrays were never populated.
- Update profile QC summaries and delayed-mode metadata for the targeted profiles.
- Append scientific calibration and per-parameter history records using the actual
  NetCDF dimensions. Existing records and untargeted profiles remain unchanged.
- Update the `DATE_UPDATE` variable rather than inventing a substitute attribute.

Filling rejected adjusted values/errors is the Argo bad-data representation,
not a numerical correction. See the [DMQC cookbook, “How to fill D files”](https://argo.ogs.it/pub/82152.pdf)
and the [Argo documentation catalogue](https://www.argodatamgt.org/Documentation).

The writer requires the standard core profile variables and metadata dimensions;
incomplete or unsupported templates stop with an error. This version does not
synthesize an entirely missing metadata schema. It supports fixed and unlimited
history dimensions and expands fixed calibration dimensions when needed.

All files are preflighted and staged before publication. Each D-file is reopened
and checked for expected flags/modes, fill values, and preservation of original
content outside intended edits. Each audit report records source/output hashes,
checker decisions, reasons, and sample counts. Publication is atomic per file,
not a transaction across the whole batch or the D-file/report pair. These are
local structural checks, not the complete GDAC submission validator.

## Later extensions

Keep operation and target separate. The target model reserves distinct sample
indices and pressure-range selectors. Later operations can carry corrected
adjusted values or uncertainties. They must be explicitly implemented and
validated; this version will not silently accept them. Corrections need conflict
rules distinct from the current bad/no-finding merge policy. The final YAML
schema, required-checker list, and review workflow can be agreed with partners.

## Tests

```bash
python -m unittest test_dmqc -v
```

Tests use synthetic NetCDFs, temporary directories, and mocked downloads. They
exercise merging, masks, profile targeting, metadata dimensions, stale sources,
overwrite protection, dry-runs, and the command-line defaults.
