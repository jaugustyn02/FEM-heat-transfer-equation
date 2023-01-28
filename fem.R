# Jan Augustyn 410105 - gr. 9
# FEM - 4.1 Równanie transportu ciepła

# Wczytanie liczby elementów od użytkownika
repeat {
    n <- as.integer(readline("Enter number of elements: "))
    if (n > 0) break
}

# Dziedzina x
a <- 0
b <- 2

# Liczba funkcji bazowych
num_e <- n + 1

# Długość pojedynczego przedziału
interval_len <- (b - a) / n

# Współczynniki liniowe funkcji bazowych
a1 <- 1 / interval_len
a2 <- -a1

# Wartości i-tej funkcji bazowej w punkcie x
e <- function(i, x) {
    t <- (x / interval_len) + 1
    if (i - 1 <= t && t <= i)
        return(a1 * x - i + 2)
    if (i < t && t <= i + 1)
        return(a2 * x + i)
    return(0)
}

# Wartość pochodnej i-tej funkcji bazowej w punkcie x
e_prim <- function(i, x) {
    t <- (x / interval_len) + 1
    if (i - 1 <= t && t <= i)
        return(a1)
    if (i < t && t <= i + 1)
        return(a2)
    return(0)
}

# Całkowanie kwadraturą Gaussa-Legendre'a dla 2 punktów Gaussa
def_integral <- function(f, a, b) {
    c <- (b - a) / 2
    d <- (b + a) / 2
    u <- c(-1 / sqrt(3), 1 / sqrt(3))
    return(c * sum(sapply(d + c * u, f)))
}

# Iloczyn pochodnych i-tej oraz j-tej funkcji bazowej
e_prims_product <- function(i, j) {
    return(
        function(x) {
            return(e_prim(i, x) * e_prim(j, x))
        }
    )
}

# Wartość lewej strony sformulowania wariacyjnego B(u,v)
# gdzie u(x), v(x) to odpowiednio i-ta, j-ta funkcja bazowa
B <- function(i, j) {
    lower <- max(0, i - 2, j - 2) * interval_len
    upper <- min(i, j) * interval_len
    return(
        e(i, 0) * e(j, 0) - def_integral(e_prims_product(i, j), lower, upper)
    )
}

# Wartość prawej strony sformulowania wariacyjnego - L(v)
# gdzie v(x) to i-ta funkcja bazowa
L <- function(i, x) {
    return(20 * e(i, x))
}

# Wyliczanie rozwiązania równania transportu ciepła
solution <- function() {
    # Główna macierz układu równań
    m <- matrix(0, num_e, num_e)
    for (i in 1:num_e)
        for (j in 1:num_e)
            m[i, j] <- B(i, j)

    # Macierz wyrazów wolnych
    c <- matrix(0, num_e, 1)
    for (i in 1:num_e)
        c[i, 1] <- L(i, 0)

    # Macierz wynikowa - rozwiązanie układu równań
    w <- solve(m, c)

    # Kombinacja liniowa funkcji bazowych
    result <- function(x) {
        sum <- 0
        for (i in 1:num_e)
            sum <- sum + w[i] * e(i, x)
        return(sum)
    }
    return(result)
}

# Rysowanie wykresu rozwiązania równania
plot_result <- function(u) {
    x <- seq(0, 2, 1 / (100 * n))
    y <- mapply(u, x)
    plot(x, y,
        main = "Solution of the heat transfer equation",
        type = "l",
        xlab = "x",
        ylab = "y = f(x)",
        cex.lab = 1.75,
        cex.axis = 1.5,
        cex.main = 2.25
    )
}

plot_result(solution())