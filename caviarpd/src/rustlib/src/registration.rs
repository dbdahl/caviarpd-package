// Generated by cargo: do not edit by hand

// If usage of .Call()/.Kall() functions in the package's R code changes, update
// this file by rerunning "cargo::register_calls(DIR)", where DIR is the root
// directory of this package.

/*
// Below is commented-out skeleton code that you can copy to your
// 'src/rustlib/src/lib.rs' file. You can change the body and arguments
// names of the functions, but changing the function name necessitates
// a corresponding change in the R code.

mod registration;
use roxido::*;

#[roxido]
fn sample_epa(nSamples: SEXP, similarity: SEXP, unnamed1: SEXP, discount: SEXP, nCores: SEXP) -> SEXP {
    r::nil()
}

#[roxido]
fn caviarpd_n_clusters(nSamplesSearch: SEXP, similarity: SEXP, mass: SEXP, discount: SEXP, unnamed1: SEXP, unnamed2: SEXP, maxNClusters: SEXP, nCores: SEXP) -> SEXP {
    r::nil()
}
*/

use roxido::*;

#[no_mangle]
extern "C" fn R_init_caviarpd_librust(info: *mut rbindings::DllInfo) {
    let mut call_routines = Vec::with_capacity(2);
    let mut names = Vec::with_capacity(2);
    names.push(std::ffi::CString::new(".sample_epa").unwrap());
    call_routines.push(rbindings::R_CallMethodDef {
        name: names.last().unwrap().as_ptr(),
        fun: unsafe { std::mem::transmute(crate::sample_epa as *const u8) },
        numArgs: 5,
    });
    names.push(std::ffi::CString::new(".caviarpd_n_clusters").unwrap());
    call_routines.push(rbindings::R_CallMethodDef {
        name: names.last().unwrap().as_ptr(),
        fun: unsafe { std::mem::transmute(crate::caviarpd_n_clusters as *const u8) },
        numArgs: 8,
    });
    call_routines.push(rbindings::R_CallMethodDef {
        name: std::ptr::null(),
        fun: None,
        numArgs: 0,
    });
    unsafe {
        rbindings::R_registerRoutines(
            info,
            std::ptr::null(),
            call_routines.as_ptr(),
            std::ptr::null(),
            std::ptr::null(),
        );
        rbindings::R_useDynamicSymbols(info, 1);
        rbindings::R_forceSymbols(info, 1);
    }
}
