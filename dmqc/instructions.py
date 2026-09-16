"""Provisional YAML adapter and format-independent instruction objects.

Only whole-profile rejection and no-finding reports are executable in v1.
Target selection and operation are separate so later adapters can describe
parameter/level flags, replacement values, and uncertainties without changing
how files, profiles, and provenance are identified.
"""
from dataclasses import asdict, dataclass
import hashlib
from pathlib import Path

import yaml

from .download import file_identity


CORE_PARAMETERS = ('PRES', 'TEMP', 'PSAL')


@dataclass(frozen=True)
class Target:
    source: str
    profile_index: int
    parameters: tuple = CORE_PARAMETERS
    selection: str = 'whole_profile'
    # Reserved for later versions: distinguish indices from pressure ranges.
    sample_indices: tuple | None = None
    pressure_range: tuple | None = None


@dataclass(frozen=True)
class Instruction:
    target: Target
    action: str
    checker: str
    reason: str
    origin: str
    flag: str | None = None
    source_sha256: str | None = None


@dataclass(frozen=True)
class Decision:
    target: Target
    source_sha256: str
    suggestions: tuple

    def as_dict(self):
        return asdict(self)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read_yaml(path):
    with open(path, encoding='utf-8') as stream:
        try:
            result = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            raise ValueError(f'{path}: invalid YAML: {exc}') from exc
    if not isinstance(result, dict):
        raise ValueError(f'{path}: expected a YAML mapping')
    return result


def _keys(mapping, allowed, context):
    if not isinstance(mapping, dict):
        raise ValueError(f'{context}: expected a mapping')
    unexpected = set(mapping) - set(allowed)
    if unexpected:
        raise ValueError(f'{context}: unsupported fields: {sorted(unexpected)}')


def _text(value, field):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{field} must be nonempty text')
    return value.strip()


def parse_document(document, origin):
    """Strict provisional adapter: unsupported actions/selections fail closed."""
    _keys(document, ('schema_version', 'checker', 'instructions'), origin)
    if type(document.get('schema_version')) is not int or document['schema_version'] != 1:
        raise ValueError(f'{origin}: expected schema_version: 1')
    checker = _text(document.get('checker'), 'checker')
    entries = document.get('instructions')
    if not isinstance(entries, list):
        raise ValueError(f'{origin}: instructions must be a list')
    result = []
    for entry in entries:
        _keys(entry, ('target', 'action', 'flag', 'reason', 'status', 'source_sha256'), origin)
        target = entry.get('target')
        _keys(target, ('source', 'profile_index', 'selection', 'parameters'), origin)
        source = _text(target.get('source'), 'target.source')
        file_identity(source)  # Exact basename prevents ambiguous cycles and path traversal.
        index = target.get('profile_index')
        if type(index) is not int or index < 0:
            raise ValueError(f'{origin}: profile_index must be a zero-based nonnegative integer')
        if target.get('selection') != 'whole_profile':
            raise ValueError(f'{origin}: only selection: whole_profile is implemented')
        if target.get('parameters', list(CORE_PARAMETERS)) != list(CORE_PARAMETERS):
            raise ValueError(f'{origin}: v1 flags PRES, TEMP and PSAL together')
        status = entry.get('status', 'ready')
        if status != 'ready':
            raise ValueError(f'{origin}: decision is {status!r}; finish human review before writing')
        action = entry.get('action')
        flag = entry.get('flag')
        if action == 'flag':
            if str(flag) != '4':
                raise ValueError(f'{origin}: v1 only supports flag 4 (bad)')
            flag = '4'
        elif action == 'no_finding':
            if flag is not None:
                raise ValueError(f'{origin}: no_finding must not set a QC flag')
        else:
            raise ValueError(f'{origin}: unsupported action {action!r}')
        reason = _text(entry.get('reason'), 'reason')
        checksum = entry.get('source_sha256')
        if checksum is not None:
            if not isinstance(checksum, str) or len(checksum) != 64 or any(c not in '0123456789abcdef' for c in checksum):
                raise ValueError(f'{origin}: invalid source_sha256')
        result.append(Instruction(Target(source, index), action, checker, reason, str(origin), flag, checksum))
    return result


def load_instructions(directory, r_dir, float_id, cycles=None, legacy_dir=None):
    """Read partner files recursively, plus the existing cycles/NNN.yaml format.

    Legacy flag 1 means no finding, never an instruction to upgrade QC.
    Legacy files address profile 0. New files must name their profile explicitly.
    """
    result = []
    directory = Path(directory)
    if directory.exists():
        for path in sorted(set(directory.rglob('*.yaml')) | set(directory.rglob('*.yml'))):
            result.extend(parse_document(read_yaml(path), path))
    if legacy_dir is not None and Path(legacy_dir).exists():
        for path in sorted(Path(legacy_dir).glob('*.yaml')):
            entry = read_yaml(path)
            _keys(entry, ('cycle', 'qc_flag', 'note'), path)
            cycle_text = str(entry.get('cycle', ''))
            # A trailing D represents a descending cycle in legacy filenames.
            source = f'R{float_id}_{cycle_text.zfill(3)}.nc'
            _, cycle = file_identity(source)
            if cycles is not None and cycle not in cycles:
                continue
            flag = str(entry.get('qc_flag', ''))
            if flag not in ('1', '4'):
                raise ValueError(f'{path}: legacy adapter supports only flags 1 and 4')
            reason = entry.get('note') or ('Legacy good decision; no QC upgrade' if flag == '1' else '')
            result.append(Instruction(Target(source, 0), 'flag' if flag == '4' else 'no_finding',
                                      'legacy_cycle', _text(reason, f'{path}: note'), str(path),
                                      '4' if flag == '4' else None))
    selected = []
    for item in result:
        source_float, cycle = file_identity(item.target.source)
        if source_float != str(float_id):
            raise ValueError(f'{item.origin}: source belongs to float {source_float}, not {float_id}')
        if cycles is None or cycle in cycles:
            if not (Path(r_dir) / item.target.source).is_file():
                raise ValueError(f'{item.origin}: missing source {item.target.source}')
            selected.append(item)
    return selected


