# Argo D-file workflow

`DMQC_checker.py` preserves the historical implementation. Use `dmqc_process.py`
for the new framework described below.

`dmqc_process.py` is the entry point. Scientific checkers run independently and
supply decisions; this tool downloads R-files, combines decisions, and writes
D-files for every selected R-file and every profile index. Rejection instructions
flag whole core profiles; other profiles retain their existing QC. The tool does
not run checks, upgrade all values to good, or calculate corrections.

## Quick start

```bash
python dmqc_process.py
```

No arguments retains float **6903708** and the historical directory:

- Windows: `C:\Data\ARGO_Dataa\DMQCprocessing`
- WSL/Linux default: `/mnt/c/Data/ARGO_Dataa/DMQCprocessing`

The default `all` stage downloads missing R-files, reads available instructions,
and writes one D-file for every selected R-file, including those without rejection
instructions. Every profile index is processed. Existing R-files are retained.
An empty rejection list still exports all selected files with their existing QC.
Existing D-files require
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
the partner instruction directory. All decisions come from this directory.
`--dry-run` writes nothing. A download dry-run still reads the remote listing;
`write --dry-run` is entirely local. An `all --dry-run` can validate only sources
already downloaded; missing source files referenced by instructions are errors.

The `combine` command writes a YAML review plan to stdout and status messages to
stderr, so it can be saved with `> combined.yaml`. This plan is an audit preview,
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
    D/                          generated D-files
    reports/                    *.report.yaml audit records
    instructions/               partner YAML files, searched recursively
      partner_a/
      partner_b/
    meta.yaml                   optional institution/operator/software information
```

Do not put generated combined plans in `instructions/`. The tool does not create
automatic good decisions or launch partner checkers. Run automatic checkers and
complete human reviews before the writing stage. A present instruction with
`status: pending` blocks processing. This first version has no required-checker
manifest: an absent checker report is not detected as unfinished work.

## Inspecting and saving flags

```bash
python DMQC_inspector.py
```

- **Up/Down** selects a cycle; **P** changes profile index.
- **Space** toggles all indices in the selected cycle to the displayed index's new state.
- **F** toggles only the displayed index.
- **Enter** writes the instructions and keeps the figure open.
- **Q** writes the instructions and closes the figure only if saving succeeds.

The inspector saves to `<float>/instructions/visual_inspector.yaml`. Override
that path with `--instructions /path/to/report.yaml`. Each flagged index gets
one instruction. Saving replaces only this inspector's report; saving no flags
writes `instructions: []`, removing its previous rejections. Other checkers'
reports are untouched. Reopening restores saved flags, including indices not
currently displayed. The figure shows save results and unsaved changes.
Closing the window normally does not save: use Enter first, or Q to save and quit.
`--save` still exports a PNG and is separate from instruction saving.

The inspector records source hashes from the start of inspection and rejects
changed sources when loading saved flags or saving new ones. It also rejects
saved targets that are unavailable in the opened set of files. Saving errors
leave the window open and flags in memory. The inspector saves suggestions only;
run `python dmqc_process.py write --dry-run` and then `write` to create D-files.

## Provisional instruction adapter

A minimal file requires only the version, checker name, and explicit operations:

```yaml
schema_version: 1
checker: my_check
instructions:
  - target:
      source: R6903708_001.nc
      profile_index: 0
      selection: whole_profile
    action: flag
    flag: '4'
```

`reason`, `status`, and `source_sha256` are optional. An omitted reason becomes
empty text internally; omitted status means `ready`; omitted hash means the
checker did not bind its suggestion to a specific source version. An explicit
`status: pending` still blocks D-file processing.

Optional file-level `metadata` accepts `checker_version`, `operator`, and
`created_utc`, each as text. Metadata is retained with each parsed instruction
for the combined plan and D-file audit report. An example is
`metadata: {checker_version: '0.1', operator: 'Example reviewer'}`.

Examples are also available in `examples/dmqc/minimal.yaml` and
`examples/dmqc/with_metadata.yaml`. These are format examples, not actual QC
findings to copy into a production instruction directory.

### Python writer API

Partners can write YAML in any language. Python partners can use:

```python
from dmqc.instructions import write_profile_flags

write_profile_flags(
    'instructions/my_check.yaml',
    {('R6903708_001.nc', 0), ('R6903708_001.nc', 1)},
    checker='my_check',
)
```

`write_instructions(path, checker, instructions, metadata=None)` accepts a list
of plain instruction dictionaries for more general use. Both functions validate
the document before atomically replacing the destination. Passing an empty list
or empty set clears that checker's previous instructions.

Validate documents without writing D-files using `python dmqc_process.py combine`
(with `--instructions-dir` if needed). This also checks targets against local
R-files. To validate only the YAML structure, Python callers can use
`parse_document(read_yaml(path), path)` from `dmqc.instructions`.

### Example with optional details

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
and may provide a reason. This is **not** a QC upgrade or approval of the profile.
Bad wins independently of checker order; every contributing suggestion remains
in the combined decision and audit report. Missing-data flags remain missing.
The output plan also includes profiles without suggestions, with outcome `retain`.
No-finding and absent rejection instructions preserve existing QC, including
existing bad/questionable flags; they do not upgrade those flags to good.

Include `source_sha256` when possible. The combiner checks it against the local
source. Without it, the tool cannot tell whether a checker used an older source;
it still binds the combined decision to the current file hash before writing.

### Metadata

Only `institution`, `operator`, `software`, and `software_release` are read from
`meta.yaml`. Defaults are `IF`, `unknown`, `ADMW`, and `0.1`. Set the operator and
software information for real processing. Old mode maps, default errors, and
QC defaults are ignored. NetCDF text fields have fixed widths; full reasons and
metadata are retained in the YAML audit report even when NetCDF text is truncated.

## Writer behavior

For every profile in each selected R-file:

- Keep raw measurements and raw QC unchanged.
- If rejected, set adjusted QC to `4` for present samples and `9` for missing samples,
  and fill adjusted values/errors. Raw data identify samples even if adjusted
  arrays were never populated.
- Otherwise retain existing adjusted values and QC. Where adjusted values are
  absent and their QC is blank/unset, copy raw values and raw QC into the adjusted
  fields. Existing QC 3/4/9 is never upgraded. Values/errors with QC 4/9 are filled.
- Preserve existing uncertainty estimates for retained values. Missing estimates
  remain missing; newly copied raw values receive no invented estimate. Reports
  count retained samples without uncertainties, and the CLI prints their total.
- Recompute profile QC summaries and set delayed-mode metadata for all profiles.
- Append scientific calibration and per-parameter history records using the actual
  NetCDF dimensions. Existing calibration/history records remain unchanged.
- Update the `DATE_UPDATE` variable rather than inventing a substitute attribute.

`--cycles` restricts the source/output set, including when no instructions exist.
Use `write --overwrite` after changing the rejection list: outputs are always
rebuilt from the original R-files, so removing a rejection restores the source
values/QC instead of retaining the previous D-file's rejection.

This is a flagging-only export workflow. In particular, absent uncertainty
estimates still need to be supplied by a later checker before the outputs can
serve as fully assessed delayed-mode data. Creating a D-file does not estimate
sensor drift or certify its uncertainty.

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
checker decisions, reasons, and sample counts. Reports are saved in the sibling
`reports/` directory, referenced from each D-file as `../reports/<name>.report.yaml`.
Publication is atomic per file,
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
python -m unittest test_dmqc test_dmqc_inspector -v
```

Tests use synthetic NetCDFs, temporary directories, and mocked downloads. They
exercise merging, masks, profile targeting, metadata dimensions, stale sources,
overwrite protection, dry-runs, and the command-line defaults.
