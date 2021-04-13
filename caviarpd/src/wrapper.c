#define R_NO_REMAP
#define STRICT_R_HEADERS
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Import C headers for rust API
#include "rustlib/dahl-randompartition.h"

// Actual Wrappers
Rr_Sexp sexp_to_rr_sexp(SEXP sexp) {
  Rr_Sexp s;
  s.sexp_ptr = (void*) sexp;
  return s;
}

Rr_Sexp null_rr_sexp() {
  Rr_Sexp s;
  s.sexp_ptr = (void*) 0;
  return s;
}

Rr_Sexp_Vector_Intsxp rrAllocVectorINTSXP(int len) {
  Rr_Sexp_Vector_Intsxp s;
  s.sexp_ptr = (void*) PROTECT(Rf_allocVector(INTSXP, len));
  s.data_ptr = INTEGER((SEXP)s.sexp_ptr);
  s.len = len;
  return s;
}

SEXP new_EpaParameters(SEXP similarity_sexp, SEXP permutation_sexp, SEXP use_natural_permutation_sexp, SEXP mass_sexp, SEXP discount_sexp) {
  int n_items = Rf_nrows(similarity_sexp);
  similarity_sexp = PROTECT(Rf_coerceVector(similarity_sexp, REALSXP));
  double* similarity = REAL(similarity_sexp);
  permutation_sexp = PROTECT(Rf_coerceVector(permutation_sexp, INTSXP));
  int* permutation = INTEGER(permutation_sexp);
  int use_natural_permutation = Rf_asLogical(use_natural_permutation_sexp);
  double mass = Rf_asReal(mass_sexp);
  double discount = Rf_asReal(discount_sexp);
  EpaParameters *ptr = dahl_randompartition__epaparameters_new(n_items, similarity, permutation, use_natural_permutation, mass, discount);
  UNPROTECT(2);
  return R_MakeExternalPtr(ptr, R_UnboundValue, R_UnboundValue);
}

SEXP free_EpaParameters(SEXP ptr_sexp) {
  EpaParameters *ptr = R_ExternalPtrAddr(ptr_sexp);
  dahl_randompartition__epaparameters_free(ptr);
  R_ClearExternalPtr(ptr_sexp);
  return R_UnboundValue;
}

SEXP samplePartition(SEXP n_samples_sexp, SEXP n_items_sexp, SEXP seed_sexp, SEXP prior_id_sexp, SEXP prior_ptr_sexp, SEXP randomize_permutation_sexp) {
  int n_samples = Rf_asInteger(n_samples_sexp);
  int n_items = Rf_asInteger(n_items_sexp);
  SEXP partition_labels_sexp = PROTECT(Rf_allocMatrix(INTSXP, n_samples, n_items));
  int *partition_labels = INTEGER(partition_labels_sexp);
  int *seed = INTEGER(seed_sexp);
  int prior_id = Rf_asInteger(prior_id_sexp);
  void *prior_ptr = R_ExternalPtrAddr(prior_ptr_sexp);
  int randomize_permutation = Rf_asLogical(randomize_permutation_sexp);
  dahl_randompartition__sample_partition(n_samples, n_items, partition_labels, seed, prior_id, prior_ptr, randomize_permutation);
  UNPROTECT(1);
  return partition_labels_sexp;
}

// Standard R package stuff
static const R_CallMethodDef CallEntries[] = {
  {".new_EpaParameters", (DL_FUNC) &new_EpaParameters, 5},
  {".free_EpaParameters", (DL_FUNC) &free_EpaParameters, 1},
  {".samplePartition", (DL_FUNC) &samplePartition, 6},
  {NULL, NULL, 0}
};

void R_init_caviarpd(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
