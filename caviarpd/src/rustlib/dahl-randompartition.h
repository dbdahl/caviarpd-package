#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct EpaParameters EpaParameters;

typedef struct Rr_Sexp {
  const void *sexp_ptr;
} Rr_Sexp;

struct EpaParameters *dahl_randompartition__epaparameters_new(int32_t n_items,
                                                              double *similarity_ptr,
                                                              const int32_t *permutation_ptr,
                                                              int32_t use_natural_permutation,
                                                              double mass,
                                                              double discount);

void dahl_randompartition__epaparameters_free(struct EpaParameters *obj);

void dahl_randompartition__sample_partition(int32_t n_partitions,
                                            int32_t n_items,
                                            int32_t *partition_labels_ptr,
                                            const int32_t *seed_ptr,
                                            int32_t prior_id,
                                            const void *prior_ptr,
                                            bool randomize_permutation);

