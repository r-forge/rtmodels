`EZ2.pe` <-
function (nu, z, a, s = 0.1) 
{
    .expr3 <- s * s
    .expr4 <- -2 * nu/.expr3
    .expr6 <- exp(.expr4 * a)
    .expr8 <- exp(.expr4 * z)
    .expr9 <- .expr6 - .expr8
    .expr10 <- .expr6 - 1
    .expr12 <- 2/.expr3
    .expr14 <- .expr6 * (.expr12 * a)
    .expr20 <- .expr10^2
    .expr27 <- .expr6 * .expr4
    .value <- .expr9/.expr10
    .grad <- array(0, c(length(.value), 4), list(NULL, c("nu", 
        "z", "a", "Ter")))
    .grad[, "nu"] <- -((.expr14 - .expr8 * (.expr12 * z))/.expr10 - 
        .expr9 * .expr14/.expr20)
    .grad[, "z"] <- -(.expr8 * .expr4/.expr10)
    .grad[, "a"] <- .expr27/.expr10 - .expr9 * .expr27/.expr20
    .grad[, "Ter"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}
