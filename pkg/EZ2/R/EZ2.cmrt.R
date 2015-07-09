`EZ2.cmrt` <-
function (nu, z, a, s = 0.1) 
{
    .expr2 <- s^2
    .expr3 <- 4 * nu/.expr2
    .expr5 <- exp(.expr3 * a)
    .expr6 <- 2 * nu
    .expr7 <- z + a
    .expr10 <- exp(.expr6 * .expr7/.expr2)
    .expr12 <- .expr6/.expr2
    .expr14 <- exp(.expr12 * a)
    .expr18 <- exp(.expr6 * z/.expr2)
    .expr19 <- .expr5 + .expr10 - .expr14 - .expr18
    .expr23 <- 2 * .expr14 - 2 * .expr10
    .expr25 <- .expr19 * z + .expr23 * a
    .expr26 <- .expr25/nu
    .expr27 <- .expr14 - .expr18
    .expr28 <- .expr26/.expr27
    .expr30 <- -1 + .expr14
    .expr37 <- .expr10 * (2 * .expr7/.expr2)
    .expr41 <- .expr14 * (2/.expr2 * a)
    .expr45 <- .expr18 * (2 * z/.expr2)
    .expr60 <- .expr27^2
    .expr65 <- .expr30^2
    .expr68 <- .expr10 * .expr12
    .expr69 <- .expr18 * .expr12
    .expr73 <- 2 * .expr68
    .expr84 <- .expr14 * .expr12
    .value <- .expr28/.expr30
    .grad <- array(0, c(length(.value), 4), list(NULL, c("nu", 
        "z", "a", "Ter")))
    .grad[, "nu"] <- ((((.expr5 * (4/.expr2 * a) + .expr37 - 
        .expr41 - .expr45) * z + (2 * .expr41 - 2 * .expr37) * 
        a)/nu - .expr25/nu^2)/.expr27 - .expr26 * (.expr41 - 
        .expr45)/.expr60)/.expr30 - .expr28 * .expr41/.expr65
    .grad[, "z"] <- (((.expr68 - .expr69) * z + .expr19 - .expr73 * 
        a)/nu/.expr27 + .expr26 * .expr69/.expr60)/.expr30
    .grad[, "a"] <- (((.expr5 * .expr3 + .expr68 - .expr84) * 
        z + ((2 * .expr84 - .expr73) * a + .expr23))/nu/.expr27 - 
        .expr26 * .expr84/.expr60)/.expr30 - .expr28 * .expr84/.expr65
    .grad[, "Ter"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}
