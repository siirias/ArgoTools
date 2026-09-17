"""Provisional YAML adapter and format-independent instruction objects.

Whole-profile rejection, no-finding reports and scalar uncertainties are executable.
Target selection and operation are separate so later adapters can describe
parameter/level flags, replacement values, and uncertainties without changing
how files, profiles, and provenance are identified.
"""
from dataclasses import asdict, dataclass, replace
import hashlib
import math
import os
import tempfile
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
class FloatTarget:
    float_id: str
    parameters: tuple
    selection: str = 'all_profiles'


@dataclass(frozen=True)
class Instruction:
    target: Target | FloatTarget
    action: str
    checker: str
    reason: str
    origin: str
    flag: str | None = None
    source_sha256: str | None = None
    metadata: dict | None = None
    value: float | None = None
    scope: str = 'profile'


@dataclass(frozen=True)
class Decision:
    target: Target
    source_sha256: str
    suggestions: tuple

    @property
    def rejected(self):
        return any(item.action == 'flag' for item in self.suggestions)

    def uncertainty(self, parameter):
        by_scope = {}
        for scope in ('float', 'profile'):
            values = {s.value for s in self.suggestions if s.action == 'set_uncertainty'
                      and parameter in s.target.parameters and s.scope == scope}
            if len(values) > 1:
                raise ValueError(f'{self.target.source} profile {self.target.profile_index}: conflicting {parameter} uncertainties at {scope} scope')
            if values:
                by_scope[scope] = next(iter(values))
        return by_scope.get('profile', by_scope.get('float'))

    def as_dict(self):
        return {**asdict(self), 'outcome': 'reject' if self.rejected else 'retain',
                'uncertainties': {p: self.uncertainty(p) for p in CORE_PARAMETERS
                                  if self.uncertainty(p) is not None}}


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
    _keys(document, ('schema_version', 'checker', 'instructions', 'metadata'), origin)
    if type(document.get('schema_version')) is not int or document['schema_version'] != 1:
        raise ValueError(f'{origin}: expected schema_version: 1')
    checker = _text(document.get('checker'), 'checker')
    metadata = document.get('metadata')
    if metadata is not None and not isinstance(metadata, dict):
        raise ValueError(f'{origin}: metadata must be a mapping')
    if metadata is not None:
        _keys(metadata, ('checker_version', 'operator', 'created_utc'), origin)
        for field, value in metadata.items():
            _text(value, f'metadata.{field}')
    entries = document.get('instructions')
    if not isinstance(entries, list):
        raise ValueError(f'{origin}: instructions must be a list')
    result = []
    for entry in entries:
        if isinstance(entry, dict) and isinstance(entry.get('target'), dict) and entry['target'].get('selection') == 'all_profiles':
            result.extend(_parse_float_default(entry, checker, metadata, origin))
            continue
        _keys(entry, ('target', 'action', 'flag', 'reason', 'status', 'source_sha256', 'value'), origin)
        target = entry.get('target')
        _keys(target, ('source', 'profile_index', 'selection', 'parameters'), origin)
        source = _text(target.get('source'), 'target.source')
        file_identity(source)  # Exact basename prevents ambiguous cycles and path traversal.
        index = target.get('profile_index')
        if type(index) is not int or index < 0:
            raise ValueError(f'{origin}: profile_index must be a zero-based nonnegative integer')
        if target.get('selection') != 'whole_profile':
            raise ValueError(f'{origin}: only selection: whole_profile is implemented')
        if entry.get('action') == 'set_uncertainty' and 'parameters' not in target:
            raise ValueError(f'{origin}: set_uncertainty requires explicit target.parameters')
        parameters = target.get('parameters', list(CORE_PARAMETERS))
        if not isinstance(parameters, list) or not parameters or any(p not in CORE_PARAMETERS for p in parameters) or len(set(parameters)) != len(parameters):
            raise ValueError(f'{origin}: parameters must be a nonempty list of distinct core parameters')
        if entry.get('action') != 'set_uncertainty' and parameters != list(CORE_PARAMETERS):
            raise ValueError(f'{origin}: v1 flags PRES, TEMP and PSAL together')
        status = entry.get('status', 'ready')
        if status != 'ready':
            raise ValueError(f'{origin}: decision is {status!r}; finish human review before writing')
        action = entry.get('action')
        flag = entry.get('flag')
        value = entry.get('value')
        if action != 'set_uncertainty' and value is not None:
            raise ValueError(f'{origin}: value is only supported for set_uncertainty')
        if action == 'set_uncertainty':
            if flag is not None:
                raise ValueError(f'{origin}: uncertainty must not set a QC flag')
            value = uncertainty_value(value)
        elif action == 'flag':
            if str(flag) != '4':
                raise ValueError(f'{origin}: v1 only supports flag 4 (bad)')
            flag = '4'
        elif action == 'no_finding':
            if flag is not None:
                raise ValueError(f'{origin}: no_finding must not set a QC flag')
        else:
            raise ValueError(f'{origin}: unsupported action {action!r}')
        reason = entry.get('reason', '')
        if not isinstance(reason, str):
            raise ValueError(f'{origin}: reason must be text when supplied')
        reason = reason.strip()
        checksum = entry.get('source_sha256')
        if checksum is not None:
            if not isinstance(checksum, str) or len(checksum) != 64 or any(c not in '0123456789abcdef' for c in checksum):
                raise ValueError(f'{origin}: invalid source_sha256')
        result.append(Instruction(Target(source, index, tuple(parameters)), action, checker, reason, str(origin), flag, checksum, metadata, value))
    return result


def _parse_float_default(entry, checker, metadata, origin):
    """Normalize the compact mapping to one internal instruction per parameter."""
    _keys(entry, ('target', 'action', 'values', 'reason', 'status'), origin)
    target = entry['target']
    _keys(target, ('float', 'selection'), origin)
    float_id = target.get('float')
    if type(float_id) not in (str, int):
        raise ValueError(f'{origin}: target.float must be a WMO identifier')
    float_id = str(float_id)
    source = f'R{float_id}_000.nc'
    file_identity(source)
    if entry.get('action') != 'set_uncertainty':
        raise ValueError(f'{origin}: all_profiles only supports set_uncertainty')
    values = entry.get('values')
    _keys(values, CORE_PARAMETERS, origin)
    if not values:
        raise ValueError(f'{origin}: values must contain at least one parameter')
    result = []
    for parameter, value in values.items():
        explicit = {k: v for k, v in entry.items() if k not in ('target', 'values')}
        explicit.update(target={'source': source, 'profile_index': 0,
                                'selection': 'whole_profile', 'parameters': [parameter]}, value=value)
        normalized = parse_document({'schema_version': 1, 'checker': checker,
                                     'metadata': metadata, 'instructions': [explicit]}, origin)[0]
        result.append(replace(normalized, target=FloatTarget(float_id, (parameter,)), scope='float'))
    return result


def uncertainty_value(value):
    """Require an explicit positive finite estimate, not a flag or missing value."""
    if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
        raise ValueError('Uncertainty must be a positive finite number')
    return float(value)


def load_instructions(directory, r_dir, float_id, cycles=None):
    """Read checker YAML reports recursively from the instructions directory."""
    result = []
    directory = Path(directory)
    if directory.exists():
        for path in sorted(set(directory.rglob('*.yaml')) | set(directory.rglob('*.yml'))):
            result.extend(parse_document(read_yaml(path), path))
    selected = []
    for item in result:
        if isinstance(item.target, FloatTarget):
            if item.target.float_id != str(float_id):
                raise ValueError(f'{item.origin}: instruction belongs to float {item.target.float_id}, not {float_id}')
            selected.append(item)
            continue
        source_float, cycle = file_identity(item.target.source)
        if source_float != str(float_id):
            raise ValueError(f'{item.origin}: source belongs to float {source_float}, not {float_id}')
        if cycles is None or cycle in cycles:
            if not (Path(r_dir) / item.target.source).is_file():
                raise ValueError(f'{item.origin}: missing source {item.target.source}')
            selected.append(item)
    return selected


def write_instructions(path, checker, instructions, metadata=None):
    """Validate and atomically replace one checker's complete YAML report.

    Instructions are plain dictionaries in the documented interchange format.
    An empty list clears previous suggestions. Optional metadata, reasons,
    status and source hashes are not needed for a minimal checker.
    """
    path = Path(path)
    document = {'schema_version': 1, 'checker': checker, 'instructions': list(instructions)}
    if metadata is not None:
        document['metadata'] = metadata
    parse_document(document, path)
    content = yaml.safe_dump(document, sort_keys=False, allow_unicode=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=path.parent, prefix='.' + path.name, suffix='.tmp')
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)
    return path


def write_profile_flags(path, flagged, checker='visual_inspector', source_hashes=None,
                        reason=None, metadata=None):
    """Write (source filename, zero-based profile index) pairs as bad flags."""
    entries = []
    for source, index in sorted(set(flagged)):
        entry = {
            'target': {'source': source, 'profile_index': index, 'selection': 'whole_profile'},
            'action': 'flag', 'flag': '4',
        }
        if reason is not None:
            entry['reason'] = reason
        if source_hashes is not None:
            entry['source_sha256'] = source_hashes[source]
        entries.append(entry)
    return write_instructions(path, checker, entries, metadata)
