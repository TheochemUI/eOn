The Windows wheel import test registers torch and flang directories with the loader directly. Git bash does not hand a Windows PATH through, so importing _core.pyd still failed.
