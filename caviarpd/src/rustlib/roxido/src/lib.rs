#![allow(dead_code)]

// Help: https://docs.rs/libR-sys, https://github.com/hadley/r-internals
use libR_sys::*; 
use rand_pcg::Pcg64Mcg;

pub use libR_sys::SEXP;

pub use libR_sys;

pub struct SEXPMethods;

impl SEXPMethods {
    pub fn rng_seeded_from_r() -> Pcg64Mcg {
        unsafe {
            let name = b"sample.int\0";
            let func = Rf_install(name.as_ptr() as *const i8).protect();
            let call = Rf_lang3(
                func,
                SEXPMethods::integer(256).protect(),
                SEXPMethods::integer(16).protect(),
            );
            let result = Rf_eval(call, R_GlobalEnv).protect();
            let slice = result.as_integer_slice();
            let mut bytes: [u8; 16] = [0; 16];
            bytes
                .iter_mut()
                .zip(slice)
                .for_each(|(b, s)| *b = (*s - 1) as u8);
            Rf_unprotect(4);
            let seed = u128::from_le_bytes(bytes);
            Pcg64Mcg::new(seed)
        }
    }
    pub fn unprotect(i: i32) {
        unsafe { Rf_unprotect(i) }
    }
    pub fn integer(x: i32) -> SEXP {
        unsafe { Rf_ScalarInteger(x) }
    }
    pub fn double(x: f64) -> SEXP {
        unsafe { Rf_ScalarReal(x) }
    }
    pub fn logical(x: bool) -> SEXP {
        unsafe { Rf_ScalarLogical(x as i32) }
    }
    pub fn integer_vector(x: i32) -> SEXP {
        unsafe { Rf_ScalarInteger(x) }
    }
    pub fn double_vector(x: f64) -> SEXP {
        unsafe { Rf_ScalarReal(x) }
    }
    pub fn logical_vector(x: bool) -> SEXP {
        unsafe { Rf_ScalarLogical(x as i32) }
    }
    pub fn integer_matrix(nrow: i32, ncol: i32) -> SEXP {
        unsafe { Rf_allocMatrix(INTSXP, nrow, ncol) }
    }
    pub fn double_matrix(x: f64) -> SEXP {
        unsafe { Rf_ScalarReal(x) }
    }
    pub fn logical_matrix(x: bool) -> SEXP {
        unsafe { Rf_ScalarLogical(x as i32) }
    }
}

pub trait SEXPExt {
    fn protect(self) -> Self;
    fn as_integer(self) -> i32;
    fn as_double(self) -> f64;
    fn as_logical(self) -> i32;
    fn as_bool(self) -> bool;
    fn as_integer_slice_mut(self) -> &'static mut [i32];
    fn as_integer_slice(self) -> &'static [i32];
    fn as_double_slice_mut(self) -> &'static mut [f64];
    fn as_double_slice(self) -> &'static [f64];
    fn as_logical_slice_mut(self) -> &'static mut [i32];
    fn as_logical_slice(self) -> &'static [i32];
    fn as_raw_slice_mut(self) -> &'static mut [u8];
    fn as_raw_slice(self) -> &'static [u8];
    fn length(self) -> i32;
    fn nrow(self) -> i32;
    fn ncol(self) -> i32;
}

impl SEXPExt for SEXP {
    fn protect(self) -> Self {
        unsafe { Rf_protect(self) }
    }
    fn as_integer(self) -> i32 {
        unsafe { Rf_asInteger(self) }
    }
    fn as_double(self) -> f64 {
        unsafe { Rf_asReal(self) }
    }
    fn as_logical(self) -> i32 {
        unsafe { Rf_asLogical(self) }
    }
    fn as_bool(self) -> bool {
        unsafe { Rf_asLogical(self) != 0 }
    }
    fn as_integer_slice_mut(self) -> &'static mut [i32] {
        unsafe { std::slice::from_raw_parts_mut(INTEGER(self), Rf_xlength(self) as usize) }
    }
    fn as_integer_slice(self) -> &'static [i32] {
        unsafe { std::slice::from_raw_parts(INTEGER(self), Rf_xlength(self) as usize) }
    }
    fn as_double_slice_mut(self) -> &'static mut [f64] {
        unsafe { std::slice::from_raw_parts_mut(REAL(self), Rf_xlength(self) as usize) }
    }
    fn as_double_slice(self) -> &'static [f64] {
        unsafe { std::slice::from_raw_parts(REAL(self), Rf_xlength(self) as usize) }
    }
    fn as_logical_slice_mut(self) -> &'static mut [i32] {
        unsafe { std::slice::from_raw_parts_mut(LOGICAL(self), Rf_xlength(self) as usize) }
    }
    fn as_logical_slice(self) -> &'static [i32] {
        unsafe { std::slice::from_raw_parts(LOGICAL(self), Rf_xlength(self) as usize) }
    }
    fn as_raw_slice_mut(self) -> &'static mut [u8] {
        unsafe { std::slice::from_raw_parts_mut(RAW(self), Rf_xlength(self) as usize) }
    }
    fn as_raw_slice(self) -> &'static [u8] {
        unsafe { std::slice::from_raw_parts(RAW(self), Rf_xlength(self) as usize) }
    }
    fn length(self) -> i32 {
        unsafe { Rf_length(self) }
    }
    fn nrow(self) -> i32 {
        unsafe { Rf_nrows(self) }
    }
    fn ncol(self) -> i32 {
        unsafe { Rf_ncols(self) }
    }
}
