context("Step 1 - checking input")
load("testthat_objects/sample_stan_data.RData")
load("testthat_objects/sample_stan_fit.RData")
data(compositions_null)
data(compositions_hard_null)
data(compositions_pos_spike)
data(compositions_neg_spike)

test_that("check_vector of non-numeric vector gives error", {
    test_check_vector <- function(x){
        check_vector(parm.name="test_vector", parm=x, p=length(x),
                     verbose=FALSE)
    }
    err_string <- "must be a vector"
    expect_error(test_check_vector("string"),             err_string)
    expect_error(test_check_vector(matrix(1:10, nrow=5)), err_string)
    expect_error(test_check_vector(list(12, 13)),         err_string)
})

test_that("check_vector gives warning with mismatched lengths",{
    test_check_vector <- function(x, p){
        check_vector(parm.name="test_vector", parm=x, p=p, verbose=FALSE)
    }
    expect_warning(test_check_vector(1:10, 5), "using first p elements")
    expect_warning(test_check_vector(1:10, 12), "recycling")
})

test_that("check_vector recycles with too-short lengths", {
    test_check_vector <- function(x, p){
        suppressWarnings(check_vector(parm.name="test_vector", parm=x, p=p,
                                      verbose=FALSE))
    }
    expect_equal(test_check_vector(1:10, 5), 1:5)
    expect_equal(test_check_vector(1:10, 12), c(1:10, 1:2))
})

test_that("check_n works with all lengths", {
    test_check_n <- function(x, p){
        suppressWarnings(check_n(n=x, p=p))
    }
    expect_equal(test_check_n(1:10, 5), 1:5)
    expect_equal(test_check_n(1:10, 10), 1:10)
    expect_equal(test_check_n(1:10, 12), c(1:10, 1:2))
})

test_that("check_L fails if L not numeric", {
    err_string <- "must be numeric"
    expect_error(check_L(L="string",   p=1), err_string)
    expect_error(check_L(L=list(1, 2), p=2), err_string)
})

test_that("check_L fails if L not pxp", {
    err_string <- "must be a square matrix with the same number of columns as C"
    expect_error(check_L(L=matrix(1:10, nrow=5), p=5), err_string)
    expect_error(check_L(L=matrix(1:16, nrow=4), p=5), err_string)
})

test_that("check_L fails if L not positive definite", {
    err_string <- "not positive definite"
    expect_error(check_L(L=matrix(-16:-1, nrow=4), p=4), err_string)
})

test_that("check_L fails if L not symmetric", {
    err_string <- "must be symmetric"
    L_test <- diag(3)
    L_test[1, 2] <- 0.1
    expect_error(check_L(L=L_test, p=3), err_string)
})

test_that("check_L warns if L is a vector with mismatched length", {
    expect_warning(check_L(L=1:10, p=5),  "only first p elements")
    expect_warning(check_L(L=1:10, p=12), "recycled")
})

test_that("check_L returns a matrix", {
    test_check_L <- function(x, p) {
        suppressWarnings(check_L(L=x, p=p))
    }
    l_vec <- 1:5
    expected_class <- is_a("matrix")
    expect_that(test_check_L(diag(l_vec), 5),  expected_class)
    expect_that(test_check_L(l_vec, 5),        expected_class)
    expect_that(test_check_L(l_vec, 3),        expected_class)
    expect_that(test_check_L(l_vec, 7),        expected_class)

})

test_that("check_L works overall", {
    test_check_L <- function(x, p){
        suppressWarnings(check_L(L=x, p=p))
    }
    l <- 1:5
    l_long <- c(l, l[1:2])
    l_short <- l[1:3]
    expect_equal(test_check_L(diag(l), length(l)), diag(l))
    expect_equal(test_check_L(l, length(l)),       diag(l))
    expect_equal(test_check_L(l, length(l_short)), diag(l_short))
    expect_equal(test_check_L(l, length(l_long)),  diag(l_long))
})

test_that("check_C accepts a data frame", {
    expect_is(check_C(Data$C), "matrix")
})

test_that("check_C gives error if row sums < 1", {
    zi <- Data$C
    zi[sample(nrow(zi), 20), 1] <- 0
    expect_error(check_C(zi), "Some of the subject totals are less than 1")
})

test_that("check_C handles zero-inflation", {
    zi <- Data$C
    zi[sample(nrow(zi), 20), 1] <- 0
    zi$extra <- ifelse(abs(rowSums(zi) - 1) > 1e-7, 1-rowSums(zi), 0)
    expect_warning(check_C(zi), "values of C are zero.")
    expect_is(suppressWarnings(check_C(zi)), "matrix")
})
test_that("check_C accepts a matrix", {
    expect_is(check_C(as.matrix(Data$C)), "matrix")
})

test_that("check_C works on all included compositional datasets", {
    expect_is(check_C(compositions_null), "matrix")
    expect_is(check_C(compositions_hard_null), "matrix")
    expect_is(check_C(compositions_pos_spike), "matrix")
    expect_is(check_C(compositions_neg_spike), "matrix")
})

test_that("check_conf_alpha fails if conf_alpha is NULL or non-numeric", {
    numeric_error <- "conf_alpha must be coercible to numeric type"
    expect_error(check_conf_alpha(NULL), "conf_alpha must not be NULL")
    expect_error(check_conf_alpha("ab"), numeric_error)
})


test_that("get_stanfit returns a given stanfit object", {
    expect_is(get_stanfit(Fit), "stanfit")
    expect_equal(get_stanfit(Fit), Fit)
})

test_that("get_stanfit returns stanfit object from list", {
    Fit_list <- list(b_fit=Fit, Data=Data)
    b_fit_list <- list(b_fit=Fit, Data=Data)
    expect_is(get_stanfit(Fit_list),   "stanfit")
    expect_is(get_stanfit(b_fit_list), "stanfit")
    
    expect_equal(get_stanfit(b_fit_list), Fit)
    expect_equal(get_stanfit(Fit_list),   Fit)
})

test_that("get_stanfit gives error if no stanfit object in list", {
    expect_error(get_stanfit(list(Data=Data)), "No 'stanfit' object found")
})

test_that("get_stanfit gives error if object not of class stanfit or list", {
    expect_error(get_stanfit(12), "class 'stanfit' or 'list'")
})

test_that("get_data returns a list if banoccfit has data element", {
    Data_list <- list(Fit=Fit, Data=Data)
    data_list <- list(Fit=Fit, data=Data)
    
    expect_is(get_data(Data_list), "list")
    expect_is(get_data(data_list), "list")

    expect_equal(get_data(Data_list), Data)
    expect_equal(get_data(data_list), Data)
})

test_that("get_data gives warning if list has no data element", {
    expect_warning(get_data(list(b_data=Data)), "no data element was found")
})

test_that("get_data returns NULL if banoccfit is 'stanfit' object", {
    expect_is(get_data(Fit), "NULL")
})

test_that("get_data returns NULL if list has no data element", {
    expect_is(suppressWarnings(get_data(list(b_data=Data))), "NULL")
})
