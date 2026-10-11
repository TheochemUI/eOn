The Windows wheel import test puts torch, metatensor, and the flang runtime on PATH. Windows does not search LD_LIBRARY_PATH, so importing _core.pyd failed with a missing DLL.
