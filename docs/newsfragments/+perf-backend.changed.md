Optimize quill backend for improved logging performance.

Configure `BackendOptions` with reduced sleep duration (100us to 10us), larger initial transit buffer (256 to 2048), zero timestamp ordering grace period (single-threaded, SPSC guarantees ordering), faster flush interval (200ms to 100ms), and disabled printable char checking (numeric data only).
