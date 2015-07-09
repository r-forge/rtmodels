`EZ2.vrt` <-
function (nu, z, a, s = 0.1) 
{
    .expr4 <- s^2
    .expr5 <- -2 * nu/.expr4
    .expr7 <- exp(.expr5 * z)
    .expr8 <- .expr7 - 1
    .expr9 <- .expr8^2
    .expr10 <- -nu * .expr9
    .expr11 <- a^2
    .expr13 <- 4 * nu
    .expr14 <- .expr13 * .expr8
    .expr16 <- .expr10 * .expr11 - .expr14 * .expr11
    .expr18 <- exp(.expr5 * a)
    .expr19 <- .expr18 - 1
    .expr20 <- .expr19^2
    .expr23 <- -3 * nu
    .expr25 <- .expr13 * z
    .expr26 <- .expr25 * a
    .expr29 <- .expr23 * .expr11 + .expr26 + .expr4 * a
    .expr31 <- .expr29 * .expr8 + .expr26
    .expr35 <- .expr16/.expr20 + .expr31/.expr19 - .expr4 * z
    .expr36 <- nu^3
    .expr38 <- 2/.expr4
    .expr40 <- .expr7 * (.expr38 * z)
    .expr53 <- .expr18 * (.expr38 * a)
    .expr57 <- .expr20^2
    .expr61 <- 4 * z * a
    .expr80 <- .expr13 * a
    .expr82 <- .expr7 * .expr5
    .expr98 <- 2 * a
    .expr103 <- .expr18 * .expr5
    .value <- .expr35/.expr36
    .grad <- array(0, c(length(.value), 4), list(NULL, c("nu", 
        "z", "a", "Ter")))
    .grad[, "nu"] <- (((nu * (2 * (.expr40 * .expr8)) - .expr9) * 
        .expr11 - (4 * .expr8 - .expr13 * .expr40) * .expr11)/.expr20 + 
        .expr16 * (2 * (.expr53 * .expr19))/.expr57 + (((.expr61 - 
        3 * .expr11) * .expr8 - .expr29 * .expr40 + .expr61)/.expr19 + 
        .expr31 * .expr53/.expr20))/.expr36 - .expr35 * (3 * 
        nu^2)/.expr36^2
    .grad[, "z"] <- ((.expr80 * .expr8 + .expr29 * .expr82 + 
        .expr80)/.expr19 - (nu * (2 * (.expr82 * .expr8)) * .expr11 + 
        .expr13 * .expr82 * .expr11)/.expr20 - .expr4)/.expr36
    .grad[, "a"] <- ((.expr10 * .expr98 - .expr14 * .expr98)/.expr20 - 
        .expr16 * (2 * (.expr103 * .expr19))/.expr57 + (((.expr23 * 
        .expr98 + .expr25 + .expr4) * .expr8 + .expr25)/.expr19 - 
        .expr31 * .expr103/.expr20))/.expr36
    .grad[, "Ter"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}
