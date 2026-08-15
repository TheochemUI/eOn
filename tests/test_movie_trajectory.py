"""get_trajectory reads product-id from the dynamics.txt header."""

from __future__ import annotations

from pathlib import Path

from eon.movie import get_trajectory


def test_get_trajectory_uses_product_id_column(tmp_path: Path):
    path = tmp_path / "dynamics.txt"
    path.write_text(
        "step-number  reactant-id  process-id  product-id  step-time\n"
        "-----------------------------------------------------------\n"
        "           0             0            1           3  1.0\n"
        "           1             3            2           5  2.0\n"
    )
    assert get_trajectory(str(path)) == [0, 3, 5]


def test_get_trajectory_finds_product_id_when_columns_move(tmp_path: Path):
    path = tmp_path / "dynamics.txt"
    path.write_text(
        "step-number  product-id  reactant-id\n"
        "------------------------------------\n"
        "           0           4            0\n"
        "           1           7            4\n"
    )
    assert get_trajectory(str(path)) == [0, 4, 7]
