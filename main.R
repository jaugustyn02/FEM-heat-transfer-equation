# Jan Augustyn - gr. 9
# FEM - 4.1 Równanie transportu ciepła

repeat {
    n <- as.integer(readline("Enter number of elements: ")) # nolint
    if (n > 0) break
}

# Dziedzina x należy do [a,b]
a <- 0
b <- 2

# liczba przedziałów
num_bfn <- n + 1
# Długość pojedynczego przedziału
interval_len <- (b - a) / n


# helper functions

between <- function(x, a, b) {
    return(a <= x && x <= b)
}

######################################

# początek, środek i koniec przedziału funkcji bazowej ei
ex <- function(i) {
    xi <- (i - 1) * interval_len
    return(c(xi - interval_len, xi, xi + interval_len))
}

# --0----2/3------4/3-----2--------->

a1 <- 1 / interval_len
a2 <- -a1

# funkcje bazowe
e <- function(i, x) {
    xs <- ex(i)
    if (between(x, xs[1], xs[2]))
        return(a1 * x - i)
    if (between(x, xs[2], xs[3]))
        return(a2 * x + i)
    return(0)
}

# pochodne funkcji bazowych
e_prim <- function(i, x) {
    xs <- ex(i)
    if (between(x, xs[1], xs[2]))
        return(a1)
    if (between(x, xs[2], xs[3]))
        return(a2)
    return(0)
}

# Całkowanie kwadraturą Gaussa-Legendre'a dla 2 punktów Gaussa                                                              # nolint                                                            # nolint
def_integral <- function(F, a, b) {                                     # nolint
    A <- c(1, 1)                                                        # nolint
    c <- (b - a) / 2
    d <- (b + a) / 2
    u <- c(-1 / sqrt(3), 1 / sqrt(3))
    return(c * sum(sapply(c * u + d, F) * A))                           # nolint
}

# Prawa strona sformulowania wariacyjnego L(v)
L <- function(i, x) {
    return(20 * e(i, x))
}

e_prims_product <- function(i, j) {
    return(
        function(x) {
            return(e_prim(i, x) * e_prim(j, x))
        }
    )
}

# Lewa strona sformulowania wariacyjnego B(u,v)
B <- function(i, j) {
    lower <- min(ex(i), ex(j), 0)
    upper <- max(ex(i), ex(j), 0)
    return(e(i, 0) * e(j, 0) - def_integral(e_prims_product(i, j), lower, upper)) # nolint
}

solution <- function() {
    M <- matrix(0, num_bfn, num_bfn)
    for (i in 1:num_bfn){
        for (j in 1:num_bfn){
            M[i, j] <- B(j, i)
        }
    }

    C <- matrix(0, num_bfn, 1)
    for (i in 1:num_bfn){
        C[i, 1] <- L(i, 0)
    }


    W <- solve(M, C)

    print(W)

    return(
        function(x){
            return(sum(W * rep(e(i, x), num_bfn)))
        }
    )
}

plot_result <- function(){
  u = solution()
  plot(seq(a, b, 1/(100*num_bfn)), 
       mapply(u, seq(a, b, 1/(100 * num_bfn))),
       main = 'Rozwiązanie równania transportu ciepła',
       xlab=' ',
       ylab=' ',
       type='l')
}

# plot_result()

# effect <- function(){
#       #tworzenie macierzy M
#       M <- matrix(0, nrow = n, ncol = n)

#       for (i in 1:n){
#           for (j in 1:n){
#               M[i,j] <- L(j-1, i-1)
#           }
#       }
#       # tworzenie macierzy V
#       V <- vector()
#       for (i in 1:n){
#           V = c(V, P(i-1))
#       }
#       #funkcja rozwiązuje układ macierzy M * u = V
#       u = solve(M, V)
    
#       # kombinacja liniowa funkcji bazowych
#       effect_function <- function(x){
#           result = 0
#           for (i in 1:n){
#           result = result + u[i] * e(x, i - 1)
#           }
#           return(result)
#       }
#       return(effect_function)
#   }

#   #funkcja rysuje wykres dla równania transportu ciepła
#   plot_effect <- function(){
#     plot(seq(0, 2, 1/(100*n)), 
#         mapply(effect(), seq(0, 2, 1/(100*n))),
#         main = 'solution of the heat transfer equation',
#         type='l',
#         xlab='values x',
#         ylab='values y=f(x)',
#         cex.lab=1.5,
#         cex.axis=1.5,
#         cex.main=2)
#   }


  #funkcja rysuje wykres funkcji bazowych
  plot_basis <- function(){
    plot(seq(0, 2, 1/(100*n)), mapply(e, 1, seq(0, 2, 1/(100*n))), 
        main = 'basic functions',
        xlab='',
        ylab='',
        type='l',
        cex.axis=1.5,
        cex.main=2)
    for (i in 1:num_bfn){
      lines(seq(0, 2, 1/(10*n)), mapply(e, i, seq(0, 2, 1/(10*n))))
    }
  }

#   # wywolanie funkcji/narysowania wykresu funkcji bazowych
  
plot_basis()

#   # wywolanie funkcji/narysowania wykresu rozwiazania rownania transportu ciepla

# plot_effect()