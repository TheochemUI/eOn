Change `uncertainty_threshold` default from `0.1` to `-1` (disabled) in both
C++ and Python. Most models lack uncertainty outputs, so the previous default
triggered a noisy exception+catch in the metatomic constructor for no benefit.
