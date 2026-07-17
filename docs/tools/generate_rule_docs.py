"""Generate documentation for Snakemake rules."""

from __future__ import annotations

import argparse
import inspect
import re
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

RULE_PATTERN = re.compile(r"^\s*rule(?:\s+([A-Za-z0-9_]+))?\s*:")
USE_RULE_PATTERN = re.compile(
    r"^\s*use rule\s+([A-Za-z0-9_]+)\s+as\s+([A-Za-z0-9_]+)\s+with\s*:"
)
DOCSTRING_STARTS = ('"""', "'''")


@dataclass
class RuleDocumentation:
    name: str | None
    docstring: str | None
    dynamic_template: str | None
    source_path: Path
    line_number: int
    derived_from: str | None = None

    @property
    def display_name(self) -> str:
        if self.name:
            return self.name
        if self.dynamic_template:
            return self.dynamic_template
        return "<unnamed rule>"

    @property
    def is_dynamic(self) -> bool:
        return self.name is None and self.dynamic_template is not None


@dataclass
class FileDocumentation:
    path: Path
    docstring: str | None
    rules: list[RuleDocumentation]


def parse_docstring(lines: list[str], start_idx: int) -> tuple[str | None, int]:
    line = lines[start_idx]
    stripped = line.strip()
    quote = stripped[:3]
    _before, _, after = line.partition(quote)
    remainder = after
    content_lines: list[str] = []

    if remainder.endswith(quote) and len(remainder) > 3:
        inner = remainder[:-3]
        text = inner
        next_index = start_idx + 1
    elif remainder.strip() == "":
        idx = start_idx + 1
        while idx < len(lines):
            current_line = lines[idx]
            if quote in current_line:
                before_close, _, _ = current_line.partition(quote)
                content_lines.append(before_close)
                idx += 1
                break
            content_lines.append(current_line.rstrip("\n"))
            idx += 1
        text = "\n".join(content_lines)
        next_index = idx
    else:
        content_lines.append(remainder)
        idx = start_idx + 1
        while idx < len(lines):
            current_line = lines[idx]
            if quote in current_line:
                before_close, _, _ = current_line.partition(quote)
                content_lines.append(before_close)
                idx += 1
                break
            content_lines.append(current_line.rstrip("\n"))
            idx += 1
        text = "\n".join(content_lines)
        next_index = idx

    cleaned = inspect.cleandoc(text)
    return (cleaned or None), next_index


def normalize_dynamic_template(expr: str) -> str:
    expr = expr.strip()
    prefixes = ('f"', "f'", '"', "'")
    for prefix in prefixes:
        if expr.startswith(prefix) and expr.endswith(prefix[-1]):
            return expr[len(prefix) : -1]
    return expr


def extract_file_docstring(lines: list[str]) -> str | None:
    idx = 0
    while idx < len(lines):
        stripped = lines[idx].strip()
        if stripped == "" or stripped.startswith("#"):
            idx += 1
            continue
        if stripped.startswith(DOCSTRING_STARTS):
            docstring, _ = parse_docstring(lines, idx)
            return docstring
        break
    return None


def extract_rules_from_lines(lines: list[str], path: Path) -> list[RuleDocumentation]:
    rules: list[RuleDocumentation] = []
    idx = 0
    while idx < len(lines):
        line = lines[idx]

        use_match = USE_RULE_PATTERN.match(line)
        if use_match:
            # `use rule <base> as <alias> with:`: the derived rule inherits the
            # base's docstring at runtime, but the source has none here, so just
            # record it and point back to the base
            rules.append(
                RuleDocumentation(
                    name=use_match.group(2),
                    docstring=None,
                    dynamic_template=None,
                    source_path=path,
                    line_number=idx + 1,
                    derived_from=use_match.group(1),
                )
            )
            idx += 1
            continue

        match = RULE_PATTERN.match(line)
        if not match:
            idx += 1
            continue

        rule_name = match.group(1)
        rule_indent = len(line) - len(line.lstrip(" "))
        docstring: str | None = None
        search_idx = idx + 1
        while search_idx < len(lines):
            candidate = lines[search_idx]
            stripped = candidate.strip()
            if stripped == "" or stripped.startswith("#"):
                search_idx += 1
                continue
            indent = len(candidate) - len(candidate.lstrip(" "))
            if indent <= rule_indent:
                break
            if docstring is None and stripped.startswith(DOCSTRING_STARTS):
                docstring, search_idx = parse_docstring(lines, search_idx)
            break

        dynamic_template: str | None = None
        lookahead_idx = idx + 1
        while lookahead_idx < len(lines):
            lookahead_line = lines[lookahead_idx]
            if lookahead_idx != idx + 1 and RULE_PATTERN.match(lookahead_line):
                break
            stripped = lookahead_line.strip()
            if stripped.startswith("utils.set_last_rule_name"):
                dyn_match = re.match(
                    r"utils\.set_last_rule_name\(\s*workflow\s*,\s*(.+?)\s*\)(?:\s*#.*)?$",
                    stripped,
                )
                if dyn_match:
                    dynamic_template = normalize_dynamic_template(dyn_match.group(1))
                else:
                    dynamic_template = stripped
                break
            lookahead_idx += 1

        rules.append(
            RuleDocumentation(
                name=rule_name,
                docstring=docstring,
                dynamic_template=dynamic_template,
                source_path=path,
                line_number=idx + 1,
            )
        )
        idx += 1
    return rules


def extract_file_documentation(path: Path) -> FileDocumentation:
    lines = path.read_text().splitlines()
    file_docstring = extract_file_docstring(lines)
    rules = extract_rules_from_lines(lines, path)
    return FileDocumentation(path=path, docstring=file_docstring, rules=rules)


def generate_markdown(
    docs_by_file: dict[Path, FileDocumentation], repo_root: Path
) -> str:
    lines: list[str] = ["# Snakemake Rules Reference", ""]
    for path in sorted(docs_by_file):
        file_doc = docs_by_file[path]
        try:
            rel_path = path.relative_to(repo_root).as_posix()
        except ValueError:
            rel_path = path.as_posix()
        lines.append(f"## {rel_path} module")
        lines.append("")

        if file_doc.docstring:
            doc_lines = file_doc.docstring.splitlines()
            lines.extend(doc_lines)
            lines.append("")

        if not file_doc.rules:
            lines.append("_No rules found._")
            lines.append("")
            continue

        for rule in file_doc.rules:
            lines.append(f"**`{rule.display_name}`**")

            body_lines: list[str] = []
            if rule.is_dynamic and rule.dynamic_template:
                body_lines += [
                    ":::{note}",
                    "This rule is dynamically generated and expands in a series of rules depending on the simflow runtime configuration.",
                    ":::",
                ]

            if rule.derived_from:
                if body_lines:
                    body_lines.append("")
                body_lines += [
                    ":::{note}",
                    f"Derived from `{rule.derived_from}` via ``use rule ... with:``; "
                    "see that rule for the shared behavior.",
                    ":::",
                ]

            if rule.docstring:
                doc_lines = rule.docstring.splitlines()
                if doc_lines:
                    if body_lines:
                        body_lines.append("")
                    body_lines.extend(doc_lines)
            else:
                if body_lines:
                    body_lines.append("")
                body_lines.append("_No description provided._")

            if not body_lines:
                body_lines.append("_No description provided._")

            first_line, *remaining_lines = body_lines
            lines.append(f": {first_line}")
            for segment in remaining_lines:
                if segment == "":
                    lines.append("")
                    continue
                lines.append(f"  {segment}")

            lines.append("")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def collect_file_documentation(paths: Iterable[Path]) -> dict[Path, FileDocumentation]:
    docs_by_file: dict[Path, FileDocumentation] = {}
    for path in paths:
        docs_by_file[path] = extract_file_documentation(path)
    return docs_by_file


def build_argument_parser(
    default_rules_dir: Path, default_output_path: Path
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate Snakemake rule documentation."
    )
    parser.add_argument(
        "--rules-dir",
        type=Path,
        default=default_rules_dir,
        help="Directory containing .smk rule files (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=default_output_path,
        help="Output path for the generated Markdown (default: %(default)s)",
    )
    return parser


def main() -> None:
    script_path = Path(__file__).resolve()
    docs_dir = script_path.parents[1]
    repo_root = docs_dir.parent
    default_rules_dir = repo_root / "workflow" / "rules"
    default_output = docs_dir / "source" / "api" / "snakemake_rules.md"

    parser = build_argument_parser(default_rules_dir, default_output)
    args = parser.parse_args()

    rule_files = sorted(args.rules_dir.glob("*.smk"))
    docs_by_file = collect_file_documentation(rule_files)

    output_content = generate_markdown(docs_by_file, repo_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(output_content)


if __name__ == "__main__":
    main()
