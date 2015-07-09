`EZ2.mrt` <-
function (nu, z, a, s = 0.1, Ter) 
{
    .expr3 <- s * s
    .expr4 <- -2 * nu/.expr3
    .expr6 <- exp(.expr4 * a)
    .expr7 <- .expr6 - 1
    .expr8 <- a/.expr7
    .expr9 <- .expr8/nu
    .expr11 <- exp(.expr4 * z)
    .expr13 <- 1/nu
    .expr16 <- .expr13 * a
    .expr20 <- 2/.expr3
    .expr22 <- .expr6 * (.expr20 * a)
    .expr24 <- .expr7^2
    .expr27 <- nu^2
    .expr35 <- 1/.expr27
    .expr48 <- .expr6 * .expr4
    .value <- .expr9 * .expr11 - .expr13 * z - .expr16/.expr7 + 
        Ter
    .grad <- array(0, c(length(.value), 4), list(NULL, c("nu", 
        "z", "a", "Ter")))
    .grad[, "nu"] <- (a * .expr22/.expr24/nu - .expr8/.expr27) * 
        .expr11 - .expr9 * (.expr11 * (.expr20 * z)) + .expr35 * 
        z + (.expr35 * a/.expr7 - .expr16 * .expr22/.expr24)
    .grad[, "z"] <- .expr9 * (.expr11 * .expr4) - .expr13
    .grad[, "a"] <- (1/.expr7 - a * .expr48/.expr24)/nu * .expr11 - 
        (.expr13/.expr7 - .expr16 * .expr48/.expr24)
    .grad[, "Ter"] <- 1
    attr(.value, "gradient") <- .grad
    .value
}
