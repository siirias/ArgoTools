# DMQC instruction YAML — short guide

Four sections to reuse as slides, followed by a small reference for implementers.
This describes the **implemented instruction format, schema version 1**.

## 1. One common format between checks and the writer

```text
R-files → independent checks / human review → instruction YAMLs
                                               ↓
                                      combine → D-files
```

- Each checker can use Python, MATLAB, or another language; it writes YAML.
- Put checker outputs in the float's `instructions/` directory. Subdirectories
  and both `.yaml` and `.yml` extensions are supported.
- Each instruction says **what data it targets** and **what action to apply**.
- The framework combines instructions and writes every selected R-file's
  D-file, including profiles with no rejection instructions.
- Configuration files and output audit reports also use YAML, but have different
  structures. Only checker instructions belong in `instructions/`.

## 2. Minimal example: reject one profile

```yaml
schema_version: 1
checker: surface_salinity
instructions:
  - target:
      source: R6903708_001.nc
      profile_index: 0
      selection: whole_profile
    action: flag
    flag: '4'
```

These are all the required fields for rejection. Add more list entries for more
profiles. `source` is the exact R-file basename, without a directory.
`profile_index` is the **zero-based index inside that file**, not the cycle
number. To reject several indices in a file, give each its own instruction.

Currently rejection covers the whole profile's available core parameters:
**PRES, TEMP and PSAL**. It sets adjusted QC to bad, preserves missing-data flags,
and fills rejected adjusted values/errors. Raw measurements and raw QC remain
unchanged. Extra sensors are not flagged by this action.

## 3. Combining automatic checks and human review

The reviewer can write the same format with a higher priority:

```yaml
schema_version: 1
checker: visual_inspector
instructions:
  - target:
      source: R6903708_001.nc
      profile_index: 0
      selection: whole_profile
    action: accept
    priority: 100
    reason: "Reviewed: the automatic rejection is not justified."
```

| Action | Meaning |
| --- | --- |
| `flag` + `flag: '4'` | Reject the whole core profile. |
| `accept` | Override lower-priority rejection instructions. |
| `no_finding` | Neutral: this checker makes no rejection/acceptance decision. |

- Priority defaults to **0**; visual review normally uses **100**.
- Highest priority wins between `accept` and `flag`; **rejection wins ties**.
- `no_finding` never overrides a rejection, regardless of priority.
- `accept` does **not** reset existing sample QC to good. With no winning
  rejection, the writer retains existing QC, including bad/questionable flags.
- File order does not determine the result; audit reports retain all suggestions.

## 4. Compact uncertainty defaults for a whole float

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
    reason: "Illustrative values only; supply justified estimates."
```

Values must be positive finite numbers in the parameter's native units. These
set adjusted **uncertainties**, not measurement corrections or QC flags.
Defaults cover present core parameters across selected profiles, including
future cycles when processed. Bad/missing samples keep filled errors.

A profile-specific estimate overrides a float default. Conflicting estimates
at the same scope stop processing; identical estimates are compatible.
Uncertainties do not use the acceptance/rejection priority rules.

## Small reference for checker authors

**Document:** `schema_version`, `checker`, and `instructions` are required.
`instructions: []` is valid: the checker has no suggestions. Replace the
checker's previous file when rerunning it so obsolete instructions do not remain.

**Optional fields on individual profile instructions:**

| Field | Meaning / restriction |
| --- | --- |
| `reason` | Human-readable explanation. |
| `priority` | Integer, default `0`; nonzero only for `flag`, `accept`, `no_finding`. |
| `status` | Defaults to `ready`; any other value stops processing. |
| `source_sha256` | Source file's 64-character lowercase SHA-256; detects changed inputs. |

Optional document-level `metadata` accepts only `checker_version`, `operator`,
and `created_utc`, all strings. Quote timestamps. There is no implemented `tags`
field. Unknown fields/actions are rejected.

For a **profile-specific uncertainty**, use this entry under `instructions`:

```yaml
- target:
    source: R6903708_001.nc
    profile_index: 0
    selection: whole_profile
    parameters: [TEMP]
  action: set_uncertainty
  value: 0.005
```

`parameters` is required here and may contain distinct core parameters present
in that profile; the same `value` applies to each. For `flag`, `accept`, and
`no_finding`, omit `parameters`: these actions always target the core set together.

Float-wide instructions accept only `target: {float, selection: all_profiles}`,
`action: set_uncertainty`, `values`, and optional `reason`/`status`.
Do not attach profile indices, source hashes, or priorities to these defaults.

YAML supports `# comments`; use spaces for indentation and quote flag codes.
**Depth/sample-specific flags, numerical corrections, and QC of extra sensors
are future extensions, not accepted instructions today.**

For commands and writer details, see [the workflow guide](DMQC_WORKFLOW.md).
