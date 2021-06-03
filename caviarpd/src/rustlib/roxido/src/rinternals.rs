#![allow(non_snake_case)]
#![allow(non_camel_case_types)]

// Help: https://github.com/hadley/r-internals

use std::os::raw::{c_char, c_double, c_int, c_uchar, c_uint, c_void};

pub type SEXP = *mut SEXPREC;
pub type Rboolean = u32; // R_ext/Boolean.h
pub type R_len_t = c_int;
pub type R_xlen_t = isize;
pub type Rbyte = c_uchar;
pub type DL_FUNC = Option<unsafe extern "C" fn() -> *mut c_void>;
pub type SEXPTYPE = c_uint;

pub const LGLSXP: c_uint = 10; // logical vectors
pub const INTSXP: c_uint = 13; // integer vectors
pub const REALSXP: c_uint = 14; // real variables
pub const RAWSXP: c_uint = 24; // raw bytes

#[repr(C)]
pub struct SEXPREC {
    _dummy: [u8; 0],
}

#[repr(C)]
pub struct DllInfo {
    _dummy: [u8; 0],
}

#[repr(C)]
pub struct R_CallMethodDef {
    pub name: *const c_char,
    pub fun: DL_FUNC,
    pub numArgs: c_int,
}

extern "C" {
    // Rinternals.h
    pub static mut R_GlobalEnv: SEXP;
    pub fn TYPEOF(x: SEXP) -> c_int;
    pub fn LOGICAL(x: SEXP) -> *mut c_int;
    pub fn INTEGER(x: SEXP) -> *mut c_int;
    pub fn REAL(x: SEXP) -> *mut c_double;
    pub fn RAW(x: SEXP) -> *mut Rbyte;
    pub fn Rf_isLogical(x: SEXP) -> Rboolean;
    pub fn Rf_isInteger(x: SEXP) -> Rboolean;
    pub fn Rf_isReal(x: SEXP) -> Rboolean;
    pub fn Rf_asLogical(x: SEXP) -> c_int;
    pub fn Rf_asInteger(x: SEXP) -> c_int;
    pub fn Rf_asReal(x: SEXP) -> c_double;
    pub fn Rf_ScalarLogical(x: c_int) -> SEXP;
    pub fn Rf_ScalarInteger(x: c_int) -> SEXP;
    pub fn Rf_ScalarReal(x: c_double) -> SEXP;
    pub fn Rf_allocVector(tipe: SEXPTYPE, len: R_xlen_t) -> SEXP;
    pub fn Rf_allocMatrix(tipe: SEXPTYPE, nrow: c_int, ncol: c_int) -> SEXP;
    pub fn Rf_allocArray(tipe: SEXPTYPE, dim: SEXP) -> SEXP;
    pub fn Rf_length(x: SEXP) -> R_len_t;
    pub fn Rf_xlength(x: SEXP) -> R_xlen_t;
    pub fn Rf_nrows(x: SEXP) -> c_int;
    pub fn Rf_ncols(x: SEXP) -> c_int;
    pub fn Rf_install(x: *const c_char) -> SEXP;
    pub fn Rf_eval(expression: SEXP, environment: SEXP) -> SEXP;
    pub fn R_tryEval(expression: SEXP, environment: SEXP, result: *mut c_int) -> SEXP;
    pub fn Rf_lang1(function: SEXP) -> SEXP;
    pub fn Rf_lang2(function: SEXP, x1: SEXP) -> SEXP;
    pub fn Rf_lang3(function: SEXP, x1: SEXP, x2: SEXP) -> SEXP;
    pub fn Rf_lang4(function: SEXP, x1: SEXP, x2: SEXP, x3: SEXP) -> SEXP;
    pub fn Rf_lang5(function: SEXP, x1: SEXP, x2: SEXP, x3: SEXP, x4: SEXP) -> SEXP;
    pub fn Rf_lang6(function: SEXP, x1: SEXP, x2: SEXP, x3: SEXP, x4: SEXP, x5: SEXP) -> SEXP;
    pub fn Rf_protect(x: SEXP) -> SEXP;
    pub fn Rf_unprotect(x: c_int);

    // R_ext/Print.h
    pub fn Rprintf(x: *const c_char, ...);

    // R_ext/Rdynload.h
    pub fn R_registerRoutines(
        info: *mut DllInfo,
        croutines: *const c_void,
        callRoutines: *const R_CallMethodDef,
        fortranRoutines: *const c_void,
        externalRoutines: *const c_void,
    ) -> c_int;
    pub fn R_useDynamicSymbols(info: *mut DllInfo, value: Rboolean) -> Rboolean;
    pub fn R_forceSymbols(info: *mut DllInfo, value: Rboolean) -> Rboolean;
}
