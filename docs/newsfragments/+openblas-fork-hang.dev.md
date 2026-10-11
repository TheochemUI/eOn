Pytest keeps one OpenBLAS thread in the test process so a fork of eonclient cannot leave server_lock held. The Ubuntu Python test job and the Python coverage job time out a hung test.
