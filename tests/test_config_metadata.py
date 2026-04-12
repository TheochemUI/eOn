"""Guard against drift between eon.schema and eon/config.yaml."""

from __future__ import annotations

import ast
import importlib.util
from pathlib import Path
import unittest
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
SCHEMA_PATH = REPO_ROOT / "eon/schema.py"


def _load_module(module_name: str, relative_path: str):
    module_path = REPO_ROOT / relative_path
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


config_module = _load_module("eon_config_for_tests", "eon/config.py")
ConfigClass = config_module.ConfigClass


SECTION_CLASS_NAMES = {
    "Saddle Search": "SaddleSearchConfig",
    "ARTn": "ARTnConfig",
    "IRA": "IRAConfig",
    "BGSD": "BGSDConfig",
}


def _config_sections() -> dict[str, dict[str, dict[str, Any]]]:
    config = ConfigClass()
    sections: dict[str, dict[str, dict[str, Any]]] = {}
    for section in config.format:
        options: dict[str, dict[str, Any]] = {}
        for key in section.keys:
            options[key.name] = {
                "default": key.default,
                "values": sorted(key.values),
            }
        sections[section.name] = options
    return sections


def _normalize_default(value: Any) -> str:
    if isinstance(value, str):
        lowered = value.lower()
        if lowered in {"true", "false"}:
            return lowered
        try:
            return format(float(value), ".15g")
        except ValueError:
            return value
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int | float):
        return format(value, ".15g")
    return str(value)


def _extract_schema_models() -> dict[str, dict[str, dict[str, Any]]]:
    tree = ast.parse(SCHEMA_PATH.read_text())
    classes = {
        node.name: node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name in SECTION_CLASS_NAMES.values()
    }
    return {
        section_name: _extract_fields(classes[class_name])
        for section_name, class_name in SECTION_CLASS_NAMES.items()
    }


def _extract_fields(class_node: ast.ClassDef) -> dict[str, dict[str, Any]]:
    fields: dict[str, dict[str, Any]] = {}
    for node in class_node.body:
        if not isinstance(node, ast.AnnAssign) or not isinstance(node.target, ast.Name):
            continue
        if not isinstance(node.value, ast.Call) or not isinstance(node.value.func, ast.Name):
            continue
        if node.value.func.id != "Field":
            continue

        default = None
        for keyword in node.value.keywords:
            if keyword.arg == "default":
                default = ast.literal_eval(keyword.value)
                break

        fields[node.target.id] = {
            "default": default,
            "values": _extract_literal_values(node.annotation),
        }
    return fields


def _extract_literal_values(annotation: ast.AST) -> list[str]:
    if not isinstance(annotation, ast.Subscript):
        return []
    if not isinstance(annotation.value, ast.Name) or annotation.value.id != "Literal":
        return []
    literal_slice = annotation.slice
    if isinstance(literal_slice, ast.Tuple):
        values = literal_slice.elts
    else:
        values = [literal_slice]
    return sorted(str(ast.literal_eval(value)) for value in values)


class ConfigMetadataDriftTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.config_sections = _config_sections()
        cls.schema_sections = _extract_schema_models()

    def test_expected_sections_exist(self) -> None:
        for section_name in SECTION_CLASS_NAMES:
            with self.subTest(section=section_name):
                self.assertIn(section_name, self.config_sections)
                self.assertIn(section_name, self.schema_sections)

    def test_schema_matches_config_yaml(self) -> None:
        for section_name, schema_fields in self.schema_sections.items():
            config_options = self.config_sections[section_name]
            self.assertEqual(set(schema_fields), set(config_options))
            for field_name, field_info in schema_fields.items():
                with self.subTest(section=section_name, field=field_name):
                    self.assertIn(field_name, config_options)
                    self.assertEqual(
                        _normalize_default(field_info["default"]),
                        _normalize_default(config_options[field_name]["default"]),
                    )

                    if field_info["values"]:
                        self.assertEqual(field_info["values"], config_options[field_name]["values"])

    def test_saddle_search_method_contract(self) -> None:
        self.assertEqual(
            self.schema_sections["Saddle Search"]["method"]["values"],
            ["artn", "basin_hopping", "bgsd", "dynamics", "min_mode"],
        )


if __name__ == "__main__":
    unittest.main()
