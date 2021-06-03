use roxido::*;

#[no_mangle]
extern "C" fn R_init_caviarpd_librust(info: *mut rbindings::DllInfo) {
    let mut call_routines = Vec::new();
    let mut names = Vec::new();
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

