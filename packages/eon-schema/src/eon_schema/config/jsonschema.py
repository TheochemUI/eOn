"""JSON Schema for eOn ``config.ini`` sections.

One document per :data:`eon_schema.config.ini.MODEL_INI_SECTION` entry, plus
``index.json``. Property names are config.ini option names: pydantic aliases,
then :data:`eon_schema.config.ini.INI_FIELD_ALIASES` (``dimer_improved`` is
``improved``). Defaults are the values a default instance carries, except
fields with a ``default_factory`` (``random_seed``), which have no fixed
default. Non-finite floats are the INI spellings ``inf`` and ``-inf`` because
JSON numbers cannot represent them.

The dimer user-guide block sets ``rotation_backend`` and ``lor_residual_tol``.
Those names are not fields on ``DimerConfig``, so that block is not copied
into ``examples``.
"""

from __future__ import annotations

import argparse
import ast
import copy
import enum
import inspect
import json
import math
import re
import textwrap
import warnings
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from pydantic import BaseModel
from pydantic.fields import FieldInfo

from eon_schema import __version__
from eon_schema.config.ini import INI_FIELD_ALIASES, MODEL_INI_SECTION

DEFAULT_SITE = "https://eondocs.org"
JSON_SCHEMA_DIALECT = "https://json-schema.org/draft/2020-12/schema"

_SECTION_HEADER = re.compile(r"^\[([^\]]+)\]\s*$")
_PYDANTIC_MODEL = re.compile(r"autopydantic_model::\s+eon\.schema\.(\w+)")
_TOCTREE_NAME = re.compile(r"^[A-Za-z0-9_]+$")
_META_DESCRIPTION = re.compile(r'"description"\s*:\s*"((?:\\.|[^"\\])*)"')


def section_slug(section: str) -> str:
    """URL slug for an INI section title (``Saddle Search`` → ``saddle-search``)."""
    return re.sub(r"[^a-z0-9]+", "-", section.lower()).strip("-")


def section_models() -> list[tuple[str, str, type[BaseModel]]]:
    """``(class name, INI section, model)`` in ``MODEL_INI_SECTION`` order."""
    import eon_schema.config.models as model_mod

    out: list[tuple[str, str, type[BaseModel]]] = []
    for cls_name, section in MODEL_INI_SECTION.items():
        cls = getattr(model_mod, cls_name)
        if not isinstance(cls, type) or not issubclass(cls, BaseModel):
            raise TypeError(f"{cls_name} is not a pydantic model")
        out.append((cls_name, section, cls))
    return out


def _explicit_alias(field: FieldInfo) -> Optional[str]:
    alias = field.alias
    if isinstance(alias, str) and alias:
        return alias
    return None


def ini_key(section: str, field_name: str, field: FieldInfo) -> str:
    """config.ini option name for one model field."""
    renamed = INI_FIELD_ALIASES.get((section, field_name))
    if renamed:
        return renamed
    return _explicit_alias(field) or field_name


def _schema_property_name(field_name: str, field: FieldInfo) -> str:
    """Name ``model_json_schema`` uses (alias when set, else the field)."""
    return _explicit_alias(field) or field_name


def _attribute_docstrings(model: type[BaseModel]) -> dict[str, str]:
    """Bare string literals that follow a field assignment in the class body."""
    try:
        source = inspect.getsource(model)
    except (OSError, TypeError):
        return {}
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", SyntaxWarning)
            tree = ast.parse(textwrap.dedent(source))
    except SyntaxError:
        return {}
    if not tree.body or not isinstance(tree.body[0], ast.ClassDef):
        return {}
    docs: dict[str, str] = {}
    body = tree.body[0].body
    for index, stmt in enumerate(body[:-1]):
        if not isinstance(stmt, ast.AnnAssign) or not isinstance(stmt.target, ast.Name):
            continue
        nxt = body[index + 1]
        if not isinstance(nxt, ast.Expr) or not isinstance(nxt.value, ast.Constant):
            continue
        text = nxt.value.value
        if not isinstance(text, str):
            continue
        cleaned = textwrap.dedent(text).strip()
        if cleaned:
            docs[stmt.target.id] = cleaned
    return docs


def _json_ready(value: Any) -> Any:
    """Value suitable for ``json.dumps`` without changing finite numbers."""
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return value
    if isinstance(value, float):
        if math.isinf(value):
            return "inf" if value > 0 else "-inf"
        if math.isnan(value):
            return "nan"
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, Mapping):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, BaseModel):
        return _json_ready(value.model_dump())
    if isinstance(value, enum.Enum):
        return _json_ready(value.value)
    raise TypeError(f"cannot encode schema default of type {type(value).__name__}")


def _instance_values(model: type[BaseModel]) -> dict[str, Any]:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="Pydantic serializer warnings",
            category=UserWarning,
        )
        inst = model()
    return {name: getattr(inst, name) for name in model.model_fields}


def section_json_schema(
    model: type[BaseModel],
    section: str,
    *,
    site: str = DEFAULT_SITE,
    filename: Optional[str] = None,
) -> dict[str, Any]:
    """JSON Schema for one config.ini section."""
    filename = filename or f"{section_slug(section)}.json"
    raw = model.model_json_schema()
    raw_props: dict[str, Any] = raw.get("properties", {})
    values = _instance_values(model)
    attr_docs = _attribute_docstrings(model)
    properties: dict[str, Any] = {}
    rename: dict[str, str] = {}
    for field_name, field in model.model_fields.items():
        src = _schema_property_name(field_name, field)
        if src not in raw_props:
            raise KeyError(
                f"{model.__name__}.{field_name} missing from JSON Schema as {src!r}"
            )
        dst = ini_key(section, field_name, field)
        if dst in properties:
            raise ValueError(f"duplicate config.ini key {dst!r} on [{section}]")
        prop = copy.deepcopy(raw_props[src])
        prop.pop("title", None)
        extra = attr_docs.get(field_name)
        if extra:
            current = prop.get("description") or ""
            if extra not in current:
                prop["description"] = (
                    f"{current}\n\n{extra}".strip() if current else extra
                )
        if field.default_factory is not None:
            prop.pop("default", None)
        else:
            prop["default"] = _json_ready(copy.deepcopy(values[field_name]))
        properties[dst] = prop
        rename[src] = dst

    cls_doc = (model.__doc__ or "").strip()
    description = cls_doc or (
        f"config.ini [{section}] section ({model.__name__}). "
        "Keys, defaults, and allowed values."
    )
    document: dict[str, Any] = {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": f"{site.rstrip('/')}/schema/{filename}",
        "title": section,
        "description": description,
        "type": "object",
        "additionalProperties": False,
        "properties": properties,
    }
    required = [rename[name] for name in raw.get("required", []) if name in rename]
    if required:
        document["required"] = required
    if "$defs" in raw:
        document["$defs"] = raw["$defs"]
    return document


def index_document(
    schemas: Sequence[tuple[str, str, str, dict[str, Any]]],
    *,
    site: str = DEFAULT_SITE,
) -> dict[str, Any]:
    """Catalog of section schema files.

    ``schemas`` is ``(class name, section, filename, schema document)``.
    """
    sections: dict[str, Any] = {}
    for cls_name, section, filename, schema in schemas:
        sections[section] = {
            "model": cls_name,
            "file": filename,
            "id": f"{site.rstrip('/')}/schema/{filename}",
            "keys": list(schema["properties"]),
        }
    return {
        "$id": f"{site.rstrip('/')}/schema/index.json",
        "title": "eOn config.ini",
        "description": (
            "Index of JSON Schema documents for config.ini. "
            "Each file is the copy of record for that section's keys, "
            "defaults, and allowed values."
        ),
        "eon_schema": __version__,
        "sections": sections,
    }


def dump_json(document: Mapping[str, Any]) -> str:
    """Pretty JSON with a trailing newline. Finite defaults stay numbers."""
    return json.dumps(document, indent=2, ensure_ascii=False) + "\n"


def write_config_schemas(out_dir: Path, *, site: str = DEFAULT_SITE) -> dict[str, Any]:
    """Write one schema file per section and ``index.json``. Return the index."""
    out_dir.mkdir(parents=True, exist_ok=True)
    built: list[tuple[str, str, str, dict[str, Any]]] = []
    used: dict[str, str] = {}
    written: set[str] = set()
    for cls_name, section, model in section_models():
        filename = f"{section_slug(section)}.json"
        if filename in used:
            raise ValueError(
                f"schema filename {filename} for [{section}] collides with [{used[filename]}]"
            )
        used[filename] = section
        schema = section_json_schema(model, section, site=site, filename=filename)
        (out_dir / filename).write_text(dump_json(schema), encoding="utf-8")
        written.add(filename)
        built.append((cls_name, section, filename, schema))
    index = index_document(built, site=site)
    (out_dir / "index.json").write_text(dump_json(index), encoding="utf-8")
    written.add("index.json")
    for stale in out_dir.glob("*.json"):
        if stale.name not in written:
            stale.unlink()
    return index


def _frontmatter(text: str) -> str:
    if not text.startswith("---"):
        return ""
    end = text.find("\n---", 3)
    if end < 0:
        return ""
    return text[3:end]


def _page_description(text: str) -> Optional[str]:
    match = _META_DESCRIPTION.search(_frontmatter(text))
    if not match:
        return None
    description = json.loads(f'"{match.group(1)}"')
    description = re.sub(r"\s+", " ", description).strip()
    return description or None


def _page_title(text: str) -> str:
    body = text
    if text.startswith("---"):
        end = text.find("\n---", 3)
        if end >= 0:
            body = text[end + 4 :]
    in_fence = False
    for line in body.splitlines():
        if line.startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        if line.startswith("# "):
            return line[2:].strip()
    raise ValueError("user guide page has no title")


def _toctree_names(index_text: str) -> list[str]:
    names: list[str] = []
    in_tree = False
    for line in index_text.splitlines():
        if line.startswith("```{toctree}"):
            in_tree = True
            continue
        if in_tree and line.startswith("```"):
            in_tree = False
            continue
        if in_tree and _TOCTREE_NAME.match(line.strip()):
            names.append(line.strip())
    return names


def user_guide_pages(guide_dir: Path) -> list[Path]:
    """Markdown pages under ``guide_dir``, index first, then its toctree order."""
    pages = {
        path.stem: path for path in sorted(guide_dir.glob("*.md")) if path.is_file()
    }
    order: list[str] = []
    index = guide_dir / "index.md"
    if index.is_file():
        order.append("index")
        for name in _toctree_names(index.read_text(encoding="utf-8")):
            if name in pages and name not in order:
                order.append(name)
    for name in sorted(pages):
        if name not in order:
            order.append(name)
    return [pages[name] for name in order]


def _sections_on_page(text: str) -> list[str]:
    by_model = MODEL_INI_SECTION
    by_section = {section: name for name, section in MODEL_INI_SECTION.items()}
    found: list[str] = []
    for line in text.splitlines():
        header = _SECTION_HEADER.match(line.strip())
        if header and header.group(1) in by_section:
            section = header.group(1)
            if section not in found:
                found.append(section)
            continue
        named = _PYDANTIC_MODEL.search(line)
        if named and named.group(1) in by_model:
            section = by_model[named.group(1)]
            if section not in found:
                found.append(section)
    return found


def render_llms_txt(
    guide_dir: Path,
    *,
    site: str = DEFAULT_SITE,
    sections: Optional[Sequence[tuple[str, str]]] = None,
) -> str:
    """llms.txt map of existing user-guide pages and config.ini schema files.

    ``sections`` is ``(INI section, filename)``. Defaults to every mapped model.
    """
    if sections is None:
        sections = [
            (section, f"{section_slug(section)}.json")
            for _name, section, _model in section_models()
        ]
    root = site.rstrip("/")
    file_for = {section: filename for section, filename in sections}
    lines = [
        "# eOn",
        "",
        "> User guide for the eOn algorithms and job types, and the config.ini section schemas.",
        "",
        "The pages below are the user guide. Each JSON Schema file is the copy of record",
        "for one config.ini section: keys, defaults, and allowed values.",
        "",
        "## User guide",
        "",
    ]
    for path in user_guide_pages(guide_dir):
        text = path.read_text(encoding="utf-8")
        title = _page_title(text)
        url = f"{root}/user_guide/{path.stem}.html"
        description = _page_description(text)
        schema_bits = []
        for section in _sections_on_page(text):
            filename = file_for.get(section)
            if not filename:
                continue
            schema_bits.append(f"[{section}]({root}/schema/{filename})")
        note = description or ""
        if schema_bits:
            schema_note = "Schema: " + ", ".join(schema_bits) + "."
            note = f"{note} {schema_note}".strip()
        if note:
            lines.append(f"- [{title}]({url}): {note}")
        else:
            lines.append(f"- [{title}]({url})")
    lines.extend(
        [
            "",
            "## config.ini schema",
            "",
            (
                f"- [Section index]({root}/schema/index.json): "
                "one JSON Schema file per config.ini section"
            ),
        ]
    )
    by_section = {section: name for name, section in MODEL_INI_SECTION.items()}
    for section, filename in sections:
        model_name = by_section.get(section, section)
        lines.append(f"- [{section}]({root}/schema/{filename}): {model_name}")
    lines.append("")
    return "\n".join(lines)


def write_llms_txt(
    path: Path,
    guide_dir: Path,
    *,
    site: str = DEFAULT_SITE,
    sections: Optional[Sequence[tuple[str, str]]] = None,
) -> Path:
    """Write ``llms.txt`` for the built site root."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        render_llms_txt(guide_dir, site=site, sections=sections),
        encoding="utf-8",
    )
    return path


def repo_root() -> Optional[Path]:
    """Monorepo root when this module sits in the eOn tree, else None."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "docs" / "source" / "conf.py").is_file() and (
            parent / "packages" / "eon-schema"
        ).is_dir():
            return parent
    return None


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Write config.ini JSON Schema files and the site llms.txt map.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        help="Directory for one JSON Schema per section and index.json",
    )
    parser.add_argument("--llms", type=Path, help="Path of llms.txt")
    parser.add_argument(
        "--user-guide",
        type=Path,
        help="Directory of user_guide markdown pages",
    )
    parser.add_argument(
        "--site", default=DEFAULT_SITE, help="Site origin, no trailing path"
    )
    args = parser.parse_args(list(argv) if argv is not None else None)
    root = repo_root()
    out = args.out
    llms = args.llms
    guide = args.user_guide
    if out is None and llms is None:
        if root is None:
            parser.error("pass --out and/or --llms (no monorepo docs tree found)")
        out = root / "docs" / "source" / "_extra" / "schema"
        llms = root / "docs" / "source" / "_extra" / "llms.txt"
        guide = guide or (root / "docs" / "source" / "user_guide")
    if out is not None:
        write_config_schemas(out, site=args.site)
    if llms is not None:
        if guide is None:
            if root is None:
                parser.error("--llms needs --user-guide")
            guide = root / "docs" / "source" / "user_guide"
        write_llms_txt(llms, guide, site=args.site)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
