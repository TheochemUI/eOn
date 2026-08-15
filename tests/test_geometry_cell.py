"""CON cell angles are crystallographic alpha, beta, gamma."""

from __future__ import annotations

import numpy as np
import pytest

from eon.geometry.cell import box_to_length_angle, length_angle_to_box


def test_length_angle_roundtrip_triclinic():
    lengths = np.array([4.0, 5.0, 6.0])
    angles = np.array([80.0, 85.0, 95.0])  # alpha, beta, gamma
    box = length_angle_to_box(lengths, angles)
    back_len, back_ang = box_to_length_angle(box)
    np.testing.assert_allclose(back_len, lengths, atol=1e-10)
    np.testing.assert_allclose(back_ang, angles, atol=1e-10)


def test_box_vectors_match_crystallographic_definition():
    box = length_angle_to_box([10.0, 11.0, 12.0], [80.0, 85.0, 95.0])
    a, b, c = box

    def ang(u, v):
        return float(
            np.degrees(
                np.arccos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
            )
        )

    assert ang(b, c) == pytest.approx(80.0, abs=1e-8)
    assert ang(a, c) == pytest.approx(85.0, abs=1e-8)
    assert ang(a, b) == pytest.approx(95.0, abs=1e-8)
    # a along x, b in the xy plane at gamma from a
    assert a[1] == pytest.approx(0.0, abs=1e-12)
    assert a[2] == pytest.approx(0.0, abs=1e-12)
    assert b[2] == pytest.approx(0.0, abs=1e-12)
