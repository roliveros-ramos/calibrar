library(testthat)
library(calibrar)

message(" \n\nTests for calibrate() --------")

tfn = function(x) sum(x^2) + 10
tgr = function(x) 2*x
foreach::registerDoSEQ()

tmp = tempdir(check = TRUE)
dir.create(file.path(tmp, "master"), recursive = TRUE)
dir.create(file.path(tmp, "run"), recursive = TRUE)
cat("test\n", file=file.path(tmp, "master", "test.txt"))

tfn2 = function(x) {
  if(!file.exists("test.txt")) stop("File not found.")
  cat(x, "\n", file="test.txt", append = TRUE)
  return(sum(x^2) + 10)
}

# calibrate ---------------------------------------------------------------

# algorithms = c('LBFGSB3', 'Rvmmin', 'Rcgmin', 'spg', 'nlminb', 'L-BFGS-B', 
#                'nmkb', 'hjkb', 'mads', 'bobyqa', 'ahr-es', 'cma-es', 'SANN', 
#                'genSA', 'DE', 'soma','genoud', 'PSO', 'PSO2007', 'PSO2011', 
#                'hybridPSO')


algorithms = c("L-BFGS-B", "nlminb", "Rcgmin", "Rvmmin", "Nelder-Mead",
               "hjn", "spg", "LBFGSB3", "CMA-ES", "genSA", "DE", "soma", "genoud", 
               "PSO", "hybridPSO", "mads", "hjkb", "nmkb", "bobyqa", "AHR-ES")

# algorithms = algorithms[3]

for(alg in algorithms) {
  
  dir.create(file.path(tmp, "run", alg), recursive = TRUE)
  
  dir.create(file.path(tmp, "run", alg, "test1"), recursive = TRUE)
  message("\n",alg, " - ", date())
  test_that(sprintf("calibrate - algorithm test: '%s'.", alg), {
    expect_no_error(
      suppressMessages(calibrate(par=rep(0.5, 5), fn = tfn,
                                 lower=rep(-100, 5), upper=rep(100, 5),
                                 method=alg, phases = c(1,1,1,2,1),
                                 control=list(master=file.path(tmp, "master"),
                                              run=file.path(tmp, "run", alg, "test1")))))
  })
  
  test_that(sprintf("calibrate - algorithm test using gr: '%s'.", alg), {
    expect_no_error(
      suppressMessages(calibrate(par=rep(0.5, 5), fn = tfn, gr=tgr, 
                                 lower=rep(-100, 5), upper=rep(100, 5),
                                 method=alg, phases = c(1,1,1,2,1))))
  })
  
  # ommiting some test due to speed
  if(!(alg %in% c("hybridPSO", "DE"))) {
    dir.create(file.path(tmp, "run", alg, "test2"), recursive = TRUE)
    test_that(sprintf("calibrate - algorithm test (parallel): '%s'.", alg), {
      expect_no_error(
        suppressMessages(calibrate(par=rep(0.5, 5), fn = tfn2,
                                   lower=rep(-100, 5), upper=rep(100, 5),
                                   method=alg, phases = c(1,1,1,2,1), parallel=TRUE,
                                   control=list(master=file.path(tmp, "master"),
                                                run=file.path(tmp, "run", alg, "test2")))))
    })
  }
  
}

alg0 = "AHR-ES"
alg1 = "Rvmmin"
dir.create(file.path(tmp, "run", "mix", "test1"), recursive = TRUE)
message("\nmix", " - ", date())
test_that(sprintf("calibrate - algorithm test: %s + %s", alg0, alg1), {
  expect_no_error(
    suppressMessages(calibrate(par=rep(0.5, 5), fn = tfn2,
                               lower=rep(-100, 5), upper=rep(100, 5),
                               method=c(alg0, alg1), phases = c(1,1,1,2,1),
                               control=list(master=file.path(tmp, "master"),
                                            run=file.path(tmp, "run", "mix", "test1")))))
  })
  
test_that("high trace and methods", {
  expect_no_error(
    suppressMessages(opt <- calibrate(par=rep(0.5, 3), fn = sphereN,
                               lower=rep(-10, 3), upper=rep(10, 3),
                               method="AHR-ES", 
                               control=list(REPORT=1, trace=4))))
  expect_no_error(opt)
  expect_no_error(coef(opt))
  expect_no_error(summ <- summary(opt))
  expect_no_error(summ)
  expect_no_error(plot(opt))
  expect_no_error(dim(opt))
})

unlink(file.path(tmp, "run"), recursive = TRUE)
