"""queued_search_count uses the communicator's own bundle size."""

from eon.communicator import Communicator


class _FakeComm:
    bundle_size = 4

    def get_queue_size(self):
        return 3

    def get_number_in_progress(self):
        return 2

    queued_search_count = Communicator.queued_search_count
    in_progress_search_count = Communicator.in_progress_search_count


def test_queued_search_count_uses_frozen_bundle_size():
    assert _FakeComm().queued_search_count() == 12


def test_in_progress_search_count_uses_frozen_bundle_size():
    assert _FakeComm().in_progress_search_count() == 8
