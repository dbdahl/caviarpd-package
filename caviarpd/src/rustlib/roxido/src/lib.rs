#![allow(dead_code)]

// Help: https://docs.rs/libR-sys, https://github.com/hadley/r-internals

pub mod r;
pub use r::SEXPExt;
pub mod rinternals;
pub use rinternals::SEXP;

