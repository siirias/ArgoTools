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

All active DMQC tools share `config/local.yaml` for their working directory and
default float. On a new checkout, copy `config/local.example.yaml` to
`config/local.yaml` and edit the two settings:

```yaml
# Parent directory containing the individual float folders.
work_directory: '/mnt/c/Data/ARGO_Dataa/DMQCprocessing'
default_float: '6903708'
```

YAML supports comments starting with `#`. The example file includes Windows,
Linux and WSL paths. Windows accepts `'C:/Data/ARGO_Dataa/DMQCprocessing'` or a
single-quoted backslash path. Linux uses paths such as `'/home/name/argo'`; WSL
uses `'/mnt/c/...'` for files on Windows C:. Use the path style of the Python
environment running the script. Windows paths are rejected under Linux/WSL to
avoid accidentally creating directories literally named `C:\...`.

`~` expands to the current user's home. Relative **configured** paths resolve
against the settings file's directory. Relative **command-line** paths retain
their normal meaning relative to the shell's current directory. The default
settings file is found beside the scripts, regardless of the shell directory.

`local.yaml` is ignored by Git; only the example is shared with partners. Your
current local settings retain the existing test float and work directory. All
five tools (`dmqc_process.py`, `dmqc_inspector.py`, `assign_uncertainties.py`,
`check_surface_salinity.py`, `verify_dfiles.py`) also accept
`--settings /path/to/another.yaml`. This is separate from the salinity checker's
`--config`, which specifies scientific thresholds.

Explicit float/R/D directory arguments override the defaults. For the processor,
`--float` and `--work-dir` independently override their configured counterparts.
If both are supplied, no settings file is needed. Other tools likewise need no
settings when their directory is supplied. `--help` works without local settings;
running without required defaults gives a short setup message. Settings are
read on each invocation, never fixed when a module is imported.

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

- `dmqc/settings.py`: shared local path/float defaults and command-line overrides.
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
python dmqc_inspector.py
```

- **Up/Down** selects a cycle; **P** changes profile index.
- **Space** toggles all indices in the selected cycle to the displayed index's new state.
- **F** toggles only the displayed index.
- **Enter** writes the instructions and keeps the figure open.
- **Q** writes the instructions and closes the figure only if saving succeeds.

- **r** clears the visual decision for this index; **Shift+R** clears the cycle.

The inspector loads other checker YAML files from `<float>/instructions/`
(or `--instructions-dir`) and displays their combined accept/reject result.
Rejected profiles are red; the selected profile shows checker names, priorities,
and reasons. It also shows the current effective priority and any visual override.

Space/F toggle the current effective outcome and create explicit `flag` or
`accept` instructions at priority **100**. Merely viewing a profile does not
create an acceptance. Clearing an override restores the outcome from the other
instructions. A checker with priority above 100 can still defeat a visual
acceptance; rejection wins ties at 100.

The inspector saves only its own decisions to
`<float>/instructions/visual_inspector.yaml` (`--instructions` overrides this
path). Other checker reports remain untouched. Saving an empty override set
clears its own file. Existing explicit priorities are preserved on reopening;
old instructions without priority still mean 0 until changed with Space/F.
Reopening restores saved acceptances as well as rejections, including indices
not currently displayed. Put separate manual instructions in another file.
Closing normally does not save: use Enter first, or Q to save and quit.
`--save` exports a PNG and is separate from instruction saving.

Source hashes guard the plotted observations and saved decisions. Changed
source files or unavailable saved targets stop loading/saving. If another
checker report is added, removed or edited while the inspector is open, saving
asks you to reopen so the display reflects the new inputs. This check is for
that inspection session; previously saved acceptances remain applicable to later
checker results according to priority. Saving errors leave the window open and
unsaved decisions in memory. Run `dmqc_process.py write` separately to apply them.

### Inspector navigation performance

Arrow navigation reuses the rendered profile cloud and cached decisions. Only
the two selected-profile overlays and status text are repainted, using blitting
when the plotting backend supports it. Header/footer space is fixed so changing
reasons does not trigger layout calculation or move the plot axes. All reasons
remain available; unusually long text is fitted into the footer with smaller type.

Flag changes and profile-index changes rebuild the cloud background. Zoom/pan,
resize and export invalidate the cached image so later navigation cannot restore
an outdated view. PNG and toolbar exports include the selected curves and labels.
Backends without blitting use a normal full redraw. Switching to an uncached
profile index with P still reads its data once; Up/Down performs no file reads.

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

`reason`, `status`, `priority`, and `source_sha256` are optional. An omitted reason becomes
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

Each whole-profile instruction covers the declared members of `PRES`, `TEMP`, and `PSAL`. An optional
`parameters: [PRES, TEMP, PSAL]` explicitly states the same scope. Other parameter
sets, depth selections, flags, and correction actions are rejected for now.

A checker that finds nothing wrong reports `action: no_finding`, omits `flag`,
and may provide a reason. This is **not** a QC upgrade or approval of the profile.
The highest integer `priority` among `flag` and `accept` instructions wins;
omitting it means 0. Rejection wins ties, independently of checker/file order.
`no_finding` is neutral even with a higher priority. Every contributing suggestion
remains in the combined decision and audit report, which also identifies the
winning priority and decisions. Missing-data flags remain missing.
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
- Apply explicit uncertainty instructions when supplied; otherwise preserve
  existing estimates for retained values. Missing estimates remain missing;
  newly copied raw values receive no invented estimate. Reports
  count retained samples without uncertainties, and the CLI prints their total.
- Recompute profile QC summaries and set delayed-mode metadata for all profiles.
- Append scientific calibration and per-parameter history records using the actual
  NetCDF dimensions. Existing calibration/history records remain unchanged.
- Update the `DATE_UPDATE` variable rather than inventing a substitute attribute.

`--cycles` restricts the source/output set, including when no instructions exist.
Use `write --overwrite` after changing the rejection list: outputs are always
rebuilt from the original R-files, so removing a rejection restores the source
values/QC instead of retaining the previous D-file's rejection.

This workflow flags data and applies supplied uncertainty estimates. Missing
estimates need to be supplied by an instruction generator before the outputs can
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
adjusted values or per-sample uncertainties. They must be explicitly implemented and
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

## Independent D-file verification

```bash
python verify_dfiles.py
python verify_dfiles.py /path/to/6903708
python verify_dfiles.py /path/to/6903708/D --verbose --no-save
```

With no arguments, uses the work directory and float from local settings.
Accepts either the float directory or its `D/` directory and scans all `D*.nc`
files there. It opens NetCDFs read-only and continues if a file is unreadable.

The console lists each file, its profile count, parameters, and result. The full
YAML inventory is saved to `reports/verification.yaml`, alongside the existing
writer reports. It includes dimensions, variable names/types/shapes, profile
sampling schemes, raw/adjusted/error sample counts, QC counts, and findings with
profile indices and sample locations. Each rerun replaces this verification
report. Use `--output path.yaml` to choose another destination, `--no-save` to
write nothing, or `--verbose` to print individual findings. Findings are grouped:
a missing-error finding may represent many measurements, recorded in `count`.

Exit codes: **0** means no errors, **1** means validation errors, and **2** means
an input/report-writing problem. `--strict` also returns 1 for warnings.
If sibling `R/` files lack corresponding D-files, the report includes a warning.

Implemented checks cover core profile format 3.1: required variables, dimensions
and types; filename/profile identity; dates and position/time QC; data modes;
adjusted values, uncertainties and QC consistency; profile QC summaries;
pressure-QC dependencies; and calibration/history records. NaN/Inf and NULL
character padding are reported. Rejected samples (QC 4/9) require filled adjusted
values/errors; other assessed adjusted samples require values and uncertainties.
The writer does not invent missing uncertainties, so otherwise
well-formed retained profiles can fail this check.

This is a local subset, **not full GDAC certification or scientific QC**. Extra
parameters and unsupported format versions produce coverage warnings. It does
not validate all attributes, reference-table vocabularies, metadata-file links,
land masks, scientific calibration correctness, BGC or trajectory rules. Rules
are based on the published [Argo User's Manual v3.3a](https://cdn.ioos.noaa.gov/media/2020/03/argo_user_manual_v3.3a.pdf)
and [CTD QC Manual v3.3](https://cdn.ioos.noaa.gov/media/2020/03/Argo-QC-for-CTD-and-Trajectory-Data.pdf),
particularly sections 3.6 and 4.7 of the latter. For submission validation, use
the [official OneArgo format checker](https://github.com/OneArgo/ArgoFormatChecker)
with the applicable specifications as well.

Run verifier tests with `python -m unittest test_verify_dfiles -v`.


## Assigning default uncertainties

```bash
python assign_uncertainties.py
python assign_uncertainties.py /path/to/6903708
python assign_uncertainties.py /path/to/6903708 --inspect
```

The default uses local settings, as in the processing script. This simple
instruction generator asks for one positive uncertainty per core parameter:
PRES in decibar, TEMP in degrees Celsius, and PSAL in psu. It applies that
constant to every cycle and profile index; it does not infer sensor accuracy,
perform scientific assessment, or offer visual/per-depth editing.

Prompts suggest a previously saved default, otherwise a uniform positive
existing `*_ADJUSTED_ERROR` value on assessed source samples if one exists.
Varying source estimates are listed as a range and are not reduced to a single
default. Measurement `resolution` is not used as uncertainty. There are no
hard-coded scientific defaults. Source units must match the supported Argo
units. An optional source/justification is recorded with each estimate.

Enter accepts a displayed suggestion; with no suggestion, Enter skips the
parameter. `skip` omits it even if previously saved. Ctrl-C cancels without
saving. `--inspect` only lists source estimates. `--output path.yaml` selects
another report location.

The script replaces its own `instructions/uncertainties.yaml`, leaving visual
rejections alone. Parameters with the same justification share one short block:

```yaml
schema_version: 1
checker: default_uncertainties
instructions:
  - target:
      float: '6903708'
      selection: all_profiles
    action: set_uncertainty
    values:
      PRES: 2.4
      TEMP: 0.002
      PSAL: 0.01
    reason: "Example only: use estimates justified for this float"
```

`value` uses the parameter's native units; it changes uncertainty only, not
measurements or QC. Float-wide defaults intentionally cover future cycles and
updated R-files, so there are no source hashes in this short instruction file.
The combiner expands defaults onto the selected files/profiles and records their
current hashes in the processing reports. Cycle filters still apply. Applying
the instructions is a separate step:

```bash
python dmqc_process.py write --overwrite
python verify_dfiles.py
```

The writer assigns these values on present adjusted samples with QC 1/2/3/5/8,
replacing existing error estimates there. QC 3 remains QC 3. Rejected profiles
and QC 4/9 samples retain fill values. Omitted parameters retain the previous
writer behaviour. Explicit profile estimates override float-wide defaults for
their named parameters, regardless of instruction order. Conflicting values at
the same scope stop combining, even if overridden or rejected; identical
suggestions are compatible. Values must be positive, finite, and representable in the output
variable without colliding with its fill value. Audit YAML preserves full
instructions and the scientific-calibration comment records the assignment.

Run `python -m unittest test_uncertainties -v` for generator/writer tests.


For an exception, put a separate instruction in e.g.
`instructions/uncertainty_overrides.yaml` so rerunning the default generator
will not replace it. The existing explicit format is still supported:

```yaml
schema_version: 1
checker: uncertainty_review
instructions:
  - target:
      source: R6903708_001.nc
      profile_index: 1
      selection: whole_profile
      parameters: [TEMP]
    action: set_uncertainty
    value: 0.005
    reason: "Example profile-specific estimate"
```

This changes only the TEMP uncertainty for that profile index. Float-wide PRES
and PSAL defaults still apply. Explicit instructions can still pin a
`source_sha256`. Existing expanded default files remain readable; rerunning the
generator rewrites them in compact form. Per-date/per-depth selectors remain
future extensions. `all_profiles` currently supports only `set_uncertainty`,
with a nonempty `values` mapping and a float identifier; source/profile selectors
and source hashes cannot be mixed into this float-wide target.


## Automatic regional surface-salinity check

```bash
python check_surface_salinity.py --dry-run
python check_surface_salinity.py
python check_surface_salinity.py /path/to/6903698 --config config/surface_salinity.yaml
```

No arguments uses the float from local settings. The supplied
`config/surface_salinity.yaml` defines each area independently using a name,
latitude range, longitude range, and maximum surface salinity. Both ranges are
required. There is no shared bounding box or inherited range.

Each rectangle includes its lower boundaries and excludes its upper boundaries
(`minimum <= coordinate < maximum`). This avoids overlapping tests at shared
edges. Positions outside every rectangle are reported as not evaluated. For
intentionally overlapping areas, all matching rules apply; exceeding any
threshold rejects. The rectangles are not coastline masks.

The initial surface definition is median raw PSAL over raw PRES 0–5 dbar
(including both pressure endpoints), requiring at least two samples. These
settings are editable, not derived from the MATLAB threshold-selection code.
`statistic` accepts `median` or `maximum`. QC filters reproduce the selection
in the supplied MATLAB plotting example: PRES_QC 1 and PSAL_QC 1 or 4. They are
explicit lists of quoted codes in YAML. Fill values and NaN/Inf are excluded.
The check does not use adjusted fields. Salinity limits use the native PSAL
practical-salinity units; pressure ranges are in dbar.

Each profile index is evaluated independently. A statistic strictly greater
than the regional limit emits the existing whole-profile core rejection;
equality passes. Missing/invalid positions, no matching area, and insufficient
eligible surface samples are reported as `not_evaluated`, never as passes.
A malformed source or invalid configuration aborts the run before saving.

Normal runs replace this checker's `instructions/surface_salinity.yaml` and
write `reports/surface_salinity.yaml`. Other checker files remain unchanged.
The instruction reasons identify the failed area, statistic, limit, sample
count and pressure range. Instructions pin source hashes. The report includes
every profile, skipped-evaluation reasons, all matched tests, source hashes,
and a snapshot of the configuration. The two output files are each replaced
atomically, but publication is not a transaction across both files.

`--dry-run` prints counts without saving. Exit code 0 means the checker completed
(including scientific rejections); 2 means an input or processing failure. Run
`dmqc_process.py write --overwrite` separately to apply the generated rejections.
Run `python -m unittest test_surface_salinity -v` for focused tests.


### Explicit acceptance and priority

```yaml
schema_version: 1
checker: visual_inspector
instructions:
  - target:
      source: R6903708_120.nc
      profile_index: 1
      selection: whole_profile
    action: accept
    priority: 100
    reason: "Expert reviewed the automatic rejection."
```

`accept` cancels lower-priority whole-profile rejection. It does not set all
sample QC to 1; the writer retains source values and QC, including existing bad
or questionable samples. `accept` must not include `flag`. Negative integer
priorities are allowed; strings, booleans and fractional values are rejected.
Priority belongs to an individual instruction, not a checker filename.
Uncertainty assignments retain their existing float-default/profile-override
resolution; nonzero priority on those assignments is rejected for now.

Tests: `python -m unittest test_priorities test_dmqc_inspector -v`.


## Profiles with different parameter inventories

`STATION_PARAMETERS` is read separately for every profile index using
`dmqc/profiles.py`. A profile can contain only PRES (and auxiliary fields such
as MTIME), even when another profile in the same file contains PRES/TEMP/PSAL.
The meaning of an index is not assumed to be constant across cycles.

- The salinity checker reports `parameter_not_available` when a profile does not
  declare both PRES and PSAL. It emits no rejection for that profile. A declared
  but missing/malformed variable is still an error. Lack of sufficient eligible
  samples for a declared parameter remains `insufficient_surface_samples`.
- The writer applies whole-profile decisions to the supported core parameters
  actually declared. It retains pressure-only profiles, updates their pressure
  values/QC/uncertainties/calibration/history, and leaves absent parameter slices
  and auxiliary data unchanged. Duplicate inventories and data for undeclared
  core parameters are rejected. Profiles with no supported core parameters are
  outside this writer's scope and stop processing with a clear error.
- Float-wide uncertainty defaults expand only onto present parameters. An
  explicit profile uncertainty targeting an absent parameter is an error, so
  a mistaken profile index is not silently ignored. The uncertainty generator
  inspects only declared parameter slices and prompts only for parameters found.
- Verification checks declared parameters without requiring absent TEMP/PSAL
  fields to be populated. It still detects undeclared measurements and missing
  declared variables. Wholly empty calibration slots inherited from R-files do
  not count as calibration records; partially filled records are checked, and
  each declared core parameter still requires a populated calibration record.
  Auxiliary/noncore parameters are inventoried with a coverage warning, because
  their scientific content/calibration is outside this core verifier's scope.

Absent fields retain their original fill values and blank flags; the writer
never fabricates salinity or temperature for a pressure-only profile. Source
hashes and per-parameter changes remain in the audit reports. Regression tests:
`python -m unittest test_partial_profiles -v`.
