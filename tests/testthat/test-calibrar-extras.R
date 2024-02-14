
message("\nCommand line arguments  --------\n")

args = c("no-save", "no-restore", "no-site-file", "no-init-file", "no-environ", "interactive")

for(arg in args) {

  test_that("command arguments", {
    expect_no_error({.get_command_argument(commandArgs(), arg)})
  })
  
}

tmp = tempdir(check = TRUE)
cat(c("parameter1 = 2\nparameter2 = 3,4,5,6\n"), file=file.path(tmp, "test.osm"))

test_that("command arguments", {
  expect_no_error({.read_configuration(file=file.path(tmp, "test.osm"), 
                                       conf.key = "model.configuration")})
})


