`EZ2.cvrt` <-
function (nu, z, a, s = 0.1) 
{
    .expr2 <- -4 * nu
    .expr4 <- s^2
    .expr5 <- 2 * nu/.expr4
    .expr7 <- exp(.expr5 * a)
    .expr8 <- .expr2 * .expr7
    .expr10 <- 2 * z
    .expr13 <- exp(.expr10 * nu/.expr4)
    .expr14 <- -1 + .expr13
    .expr15 <- .expr8 * .expr14
    .expr17 <- 4 * nu/.expr4
    .expr19 <- exp(.expr17 * a)
    .expr20 <- .expr19 - .expr13
    .expr21 <- .expr15 * .expr20
    .expr22 <- a^2
    .expr25 <- 2 * (z + a)
    .expr28 <- exp(.expr25 * nu/.expr4)
    .expr29 <- 4 * .expr28
    .expr30 <- .expr29 * nu
    .expr31 <- .expr7 - 1
    .expr32 <- .expr31^2
    .expr33 <- .expr30 * .expr32
    .expr34 <- z^2
    .expr37 <- 8 * .expr28
    .expr38 <- .expr37 * nu
    .expr39 <- .expr38 * .expr32
    .expr40 <- .expr39 * a
    .expr43 <- 2 * .expr4
    .expr44 <- .expr43 * .expr7
    .expr45 <- .expr44 * .expr14
    .expr46 <- .expr45 * .expr31
    .expr48 <- -.expr7 + .expr13
    .expr49 <- .expr46 * .expr48
    .expr52 <- .expr4 * .expr32
    .expr54 <- 4 * z
    .expr57 <- exp(.expr54 * nu/.expr4)
    .expr58 <- -.expr19 + .expr57
    .expr59 <- .expr52 * .expr58
    .expr61 <- .expr21 * .expr22 - .expr33 * .expr34 + .expr40 * 
        z + .expr49 * a - .expr59 * z
    .expr62 <- .expr61/.expr32
    .expr63 <- nu^3
    .expr64 <- .expr62/.expr63
    .expr65 <- .expr7 - .expr13
    .expr66 <- .expr65^2
    .expr70 <- .expr7 * (2/.expr4 * a)
    .expr76 <- .expr13 * (.expr10/.expr4)
    .expr82 <- .expr19 * (4/.expr4 * a)
    .expr88 <- .expr28 * (.expr25/.expr4)
    .expr94 <- 2 * (.expr70 * .expr31)
    .expr132 <- .expr32^2
    .expr147 <- .expr66^2
    .expr150 <- .expr13 * .expr5
    .expr156 <- .expr28 * .expr5
    .expr159 <- 4 * .expr156 * nu * .expr32
    .expr166 <- 8 * .expr156 * nu * .expr32
    .expr191 <- .expr7 * .expr5
    .expr195 <- .expr19 * .expr17
    .expr203 <- 2 * (.expr191 * .expr31)
    .value <- .expr64/.expr66
    .grad <- array(0, c(length(.value), 4), list(NULL, c("nu", 
        "z", "a", "Ter")))
    .grad[, "nu"] <- ((((((.expr2 * .expr70 - 4 * .expr7) * .expr14 + 
        .expr8 * .expr76) * .expr20 + .expr15 * (.expr82 - .expr76)) * 
        .expr22 - ((4 * .expr88 * nu + .expr29) * .expr32 + .expr30 * 
        .expr94) * .expr34 + ((8 * .expr88 * nu + .expr37) * 
        .expr32 + .expr38 * .expr94) * a * z + (((.expr43 * .expr70 * 
        .expr14 + .expr44 * .expr76) * .expr31 + .expr45 * .expr70) * 
        .expr48 + .expr46 * (.expr76 - .expr70)) * a - (.expr4 * 
        .expr94 * .expr58 + .expr52 * (.expr57 * (.expr54/.expr4) - 
        .expr82)) * z)/.expr32 - .expr61 * .expr94/.expr132)/.expr63 - 
        .expr62 * (3 * nu^2)/.expr63^2)/.expr66 - .expr64 * (2 * 
        ((.expr70 - .expr76) * .expr65))/.expr147
    .grad[, "z"] <- ((.expr8 * .expr150 * .expr20 - .expr15 * 
        .expr150) * .expr22 - (.expr159 * .expr34 + .expr33 * 
        .expr10) + (.expr166 * a * z + .expr40) + (.expr44 * 
        .expr150 * .expr31 * .expr48 + .expr46 * .expr150) * 
        a - (.expr52 * (.expr57 * .expr17) * z + .expr59))/.expr32/.expr63/.expr66 + 
        .expr64 * (2 * (.expr150 * .expr65))/.expr147
    .grad[, "a"] <- (((.expr2 * .expr191 * .expr14 * .expr20 + 
        .expr15 * .expr195) * .expr22 + .expr21 * (2 * a) - (.expr159 + 
        .expr30 * .expr203) * .expr34 + ((.expr166 + .expr38 * 
        .expr203) * a + .expr39) * z + (((.expr43 * .expr191 * 
        .expr14 * .expr31 + .expr45 * .expr191) * .expr48 - .expr46 * 
        .expr191) * a + .expr49) - (.expr4 * .expr203 * .expr58 - 
        .expr52 * .expr195) * z)/.expr32 - .expr61 * .expr203/.expr132)/.expr63/.expr66 - 
        .expr64 * (2 * (.expr191 * .expr65))/.expr147
    .grad[, "Ter"] <- 0
    attr(.value, "gradient") <- .grad
    .value
}
