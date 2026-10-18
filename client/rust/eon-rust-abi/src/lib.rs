//! One Rust staticlib for the C ABIs eOn links.
//!
//! readcon-core, linkcell, and minimage are rlibs here. A staticlib of each
//! would carry its own Rust standard library, and the client link would see
//! `alloc` and `std` defined more than once. `-C link-dead-code` on this
//! crate keeps every `no_mangle` C symbol from the rlibs in the archive; the
//! C++ link then pulls the members it references, including a single std.

#[used]
static EON_RUST_ABI_KEEP: [extern "C" fn() -> *const std::ffi::c_char; 3] = [
    minimage::mi_version,
    linkcell::lc_version,
    readcon_core::ffi::rkr_library_version,
];

/// Referenced so the keep table is live in the staticlib.
#[no_mangle]
pub extern "C" fn eon_rust_abi_anchor() -> usize {
    EON_RUST_ABI_KEEP.len()
}
