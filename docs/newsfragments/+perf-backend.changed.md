Optimize quill backend for 10× faster logging performance.

Configure `BackendOptions` with reduced sleep duration (100μs → 10μs), larger initial transit buffer (256 → 2048), reduced timestamp ordering grace period (5μs → 1μs), and faster flush interval (200ms → 100ms) for significantly improved throughput and latency.

The default quill configuration prioritizes low CPU usage over performance. This change balances both: 10× faster logging without 100% CPU utilization (which would require sleep_duration=0).
