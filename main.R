# Jan Augustyn - gr. 9
# FEM - 4.1 Równanie transportu ciepła

repeat {
    n <- as.integer(readline("Enter number of elements: ")) # nolint
    if (n > 0) break
}

# n <- 4

# Dziedzina x należy do [a,b]
a <- 0
b <- 2

# liczba przedziałów
num_bfn <- n + 1

# długość pojedynczego przedziału
interval_len <- (b - a) / n


# helper function
between <- function(x, a, b) {
    return(a <= x && x <= b)
}


# początek, środek i koniec przedziału funkcji bazowej ei
ex <- function(i) {
    xi <- (i - 1) * interval_len
    return(c(xi - interval_len, xi, xi + interval_len))
}

# współczynniki liniowe funkcji bazowych
a1 <- (1 / interval_len)
a2 <- (-a1)

# funkcje bazowe
e <- function(i, x) {
    xs <- ex(i)
    if (between(x, xs[1], xs[2]))
        return(a1 * x - i + 2)
    if (between(x, xs[2], xs[3]))
        return(a2 * x + i)
    return(0)
}

# pochodne funkcji bazowych
e_prim <- function(i, x) {
    xs <- ex(i)
    if (between(x, xs[1], xs[2])){
        return(a1)
    }
    if (between(x, xs[2], xs[3])){
        return(a2)
    }
    return(0)
}


# całkowanie kwadraturą Gaussa-Legendre'a dla 2 punktów Gaussa                                                              # nolint                                                            # nolint
def_integral <- function(F, a, b) {                                     # nolint
    A <- c(1, 1)                                                        # nolint
    c <- (b - a) / 2
    d <- (b + a) / 2
    u <- c(-1 / sqrt(3), 1 / sqrt(3))
    return(c * sum(sapply(c * u + d, F) * A))                           # nolint
}

# prawa strona sformulowania wariacyjnego L(v)
L <- function(i, x) {
    return(20 * e(i, x))
}

# iloczyn pochodnych funkcji bazowych
e_prims_product <- function(i, j) {
    return(
        function(x) {
            return(e_prim(i, x) * e_prim(j, x))
        }
    )
}

# lewa strona sformulowania wariacyjnego B(u,v)
B <- function(i, j) {
    lower <- max(0, ex(i)[1], ex(i)[1])
    upper <- min(ex(i)[3], ex(i)[3])
    return(e(i, 0) * e(j, 0) - def_integral(e_prims_product(i, j), lower, upper)) # nolint
}

# znalezienie rozwiązania równania
solution <- function() {
    M <- matrix(0, num_bfn, num_bfn)
    for (i in 1:num_bfn){
        for (j in 1:num_bfn){
            M[i, j] <- B(i, j)
        }
    }

    C <- matrix(0, num_bfn, 1)
    for (i in 1:num_bfn){
        C[i, 1] <- L(i, 0)
    }

    U <- solve(M, C)
    
    return(
        function(x) {
            return(sum(mapply(
                function(i, x) {U[i] * e(i, x)},
                seq(1, num_bfn, 1), x))
            )
        }
    )
}

# tworzenie wykresu rozwiązania równania
plot_result <- function() {
    u <- solution()
    x <- seq(a, b, 1 / (100 * n))
    y <- mapply(u, x)
    plot(x, y,
        main = 'Solution of the heat transfer equation',
        type='l',
        xlab='x',
        ylab='y = f(x)',
        cex.lab=1.5,
        cex.axis=1.5,
        cex.main=2
    )
}

# tworzenie wykresu funkcji bazowych
plot_basis <- function() {
    x <- seq(a, b, 1 / (100 * n))
    y <- mapply(e, 1, x)
    plot(x, y, 
        main = 'basic functions',
        xlab='x',
        ylab='y = f(x)',
        type='l',
        cex.axis=1.5,
        cex.main=2
    )
    for (i in 2: num_bfn){
        lines(x, mapply(e, i, x))
    }
}


# plot_basis()
plot_result()