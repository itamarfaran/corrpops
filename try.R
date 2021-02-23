library(corrfuncs)

data('tga_arrays')

dim_subset <- 1:18

control_subset <- tga_arrays$control[dim_subset, dim_subset, ]
diagnosed_subset <- tga_arrays$diagnosed[dim_subset, dim_subset, ]

test_corr_mat(control_subset)
test_corr_mat(diagnosed_subset)


results <- estimate_model(
  control_arr = control_subset,
  diagnosed_arr = diagnosed_subset,
  bias_correction = TRUE)


compute_gee_variance(results, control_subset, diagnosed_subset)
