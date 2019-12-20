TEST_LEVEL <- 1            # 0, 1, 2, 3, 4
TEST_EXTENSIVE <- FALSE    # FALSE or TRUE

cat(paste0("************ Starting tests at level ",TEST_LEVEL," (extensive=",TEST_EXTENSIVE,"). ************\n"))

requireLevel <- function(requiredLevel) {
  if ( TEST_LEVEL < requiredLevel ) skip(paste0("Test level is ",TEST_LEVEL," but ",requiredLevel," is required."))
}
