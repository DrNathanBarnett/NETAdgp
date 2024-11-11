#include <R.h>
#include <Rinternals.h>

SEXP hello_world() {
  Rprintf("Hello, world!\n");
  return R_NilValue;
}
