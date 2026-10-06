"""The documentation build page starts from the docs-mta environment."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGE = ROOT / "docs" / "source" / "devdocs" / "docbuild.md"
BLOCK = "```{code-block} bash\n"


def first_bash_block(text: str) -> tuple[str, str]:
    """Return the text before the first bash block and that block's body."""
    start = text.find(BLOCK)
    if start < 0:
        raise AssertionError("docbuild.md has no bash command block")
    body = text[start + len(BLOCK) :]
    end = body.find("```")
    if end < 0:
        raise AssertionError("docbuild.md bash block is not closed")
    return text[:start], body[:end]


def test_docbuild_first_command_follows_eonclient_install():
    text = PAGE.read_text(encoding="utf-8")
    before, block = first_bash_block(text)
    commands = [
        line.strip()
        for line in block.splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]
    assert commands, "the first bash block has no command"
    assert commands[0].startswith("pixi run -e docs-mta")
    assert "eonclient" in before
    assert "install" in before.lower()
    assert "pdm" not in text.lower()
