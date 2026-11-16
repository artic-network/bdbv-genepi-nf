#!/usr/bin/env python3
"""Check that nextflow_schema.json is valid and consistent.

Catches the problems that stop a schema loading or silently misbehave in a parameter form:
  - invalid JSON (for example a trailing comma)
  - the same key twice in one object
  - the same parameter defined in two sections
  - a section listed in allOf that does not exist
  - a parameter with no type, or a default of the wrong type

Usage: python .github/scripts/check_schema.py [path/to/nextflow_schema.json]
"""
import json
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "nextflow_schema.json"
errors = []


def error(message):
    errors.append(message)
    # shows as an annotation on the file in a GitHub pull request
    print(f"::error file={path}::{message}")


def no_duplicate_keys(pairs):
    seen = {}
    for key, value in pairs:
        if key in seen:
            error(f"duplicate key '{key}' in the same object")
        seen[key] = value
    return seen


try:
    with open(path) as f:
        schema = json.load(f, object_pairs_hook=no_duplicate_keys)
except json.JSONDecodeError as e:
    error(f"not valid JSON: {e.msg} at line {e.lineno}, column {e.colno}")
    sys.exit(1)

for key in ("$schema", "title", "description", "type", "definitions", "allOf"):
    if key not in schema:
        error(f"missing top-level key '{key}'")
if errors:
    sys.exit(1)

sections = schema["definitions"]

for item in schema["allOf"]:
    ref = item.get("$ref", "")
    name = ref.split("/")[-1]
    if not ref.startswith("#/definitions/") or name not in sections:
        error(f"allOf refers to '{ref}', which is not a section in definitions")

TYPES = {
    "string": lambda v: isinstance(v, str),
    "integer": lambda v: isinstance(v, int) and not isinstance(v, bool),
    "number": lambda v: isinstance(v, (int, float)) and not isinstance(v, bool),
    "boolean": lambda v: isinstance(v, bool),
}

defined_in = {}
all_properties = [(name, sec.get("properties", {})) for name, sec in sections.items()]
all_properties.append(("(top level)", schema.get("properties", {})))

for section, props in all_properties:
    for param, spec in props.items():
        defined_in.setdefault(param, []).append(section)
        if not isinstance(spec, dict) or "type" not in spec:
            error(f"parameter '{param}' in {section} has no type")
            continue
        kind = spec["type"]
        if kind not in TYPES:
            error(f"parameter '{param}' in {section} has unknown type '{kind}'")
        elif "default" in spec and spec["default"] is not None and not TYPES[kind](spec["default"]):
            error(f"parameter '{param}' in {section}: default {spec['default']!r} is not a {kind}")

for param, where in defined_in.items():
    if len(where) > 1:
        error(f"parameter '{param}' is defined in more than one place: {', '.join(where)}")

if errors:
    print(f"\n{len(errors)} problem(s) found in {path}")
    sys.exit(1)

n = sum(len(p) for _, p in all_properties)
print(f"{path}: OK ({len(sections)} sections, {n} parameters)")
