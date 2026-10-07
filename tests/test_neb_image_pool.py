"""NEB image forces run on a bounded std::thread pool and rethrow."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_image_pool_is_std_thread_and_rethrows():
    pool = (ROOT / "client" / "ForEachImage.h").read_text(encoding="utf-8")
    assert "tbb" not in pool.lower()
    assert "std::thread" in pool
    assert "std::rethrow_exception" in pool
    case = (ROOT / "client" / "unit_tests" / "ForEachImageTest.cpp").read_text(
        encoding="utf-8"
    )
    assert 'REQUIRE_THROWS_WITH(eonc::forEachImage(n, work), "image 7 failed")' in case
    neb = (ROOT / "client" / "NudgedElasticBand.cpp").read_text(encoding="utf-8")
    assert "eonc::forEachImage" in neb
    meson = (ROOT / "client" / "meson.build").read_text(encoding="utf-8")
    assert "with_parallel_neb no longer needs TBB" in meson
