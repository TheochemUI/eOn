import sys
from pathlib import Path

path = str(Path(__file__).resolve().parent.parent / "eon")
sys.path.insert(0, path)
