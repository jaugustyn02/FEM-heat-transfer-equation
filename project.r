# zadanie z Metody Elementow Skonczonych
# Jakub Kroczek
# 1.1 Przeplyw ciepla

#przypisanie liczby elementów oraz wybór wykresu
n_line <- readline('Enter the number of intervals: ');
writeLines('Enter the number:\n1) Draw the graph of the basis functions\n2) Draw a graph of the heat transfer equation');
l_line <- readline('Your choice: ');
n <- as.integer(n_line);

if (n>0) {


  # funkcja prawa strona podcalkowa

  f1 <- function(x,i){
      if(x<=1){
          return((100*x/(x+1))*e(x,i))
      }
      else if (x>1) {
        return(50*e(x,i))
      }
      return(0)
  }

  # funkcja zwraca funkcję bazową dla x oraz i
  e <- function(x, i){
    if (x > 2*(i-1)/n && x <= 2*i/n){
      return (n/2*x - i + 1)
    }
    else if ((x > 2*i/n) && x < 2*(i+1)/n){
      return (-n/2*x + i + 1)
    }
    else {
      return (0)
    }
  }

  #pochodna funkcji bazowej dla x oraz y
  e_p <- function(x, i){
    if (x > 2*(i-1)/n && x <= 2*i/n){
      return (n/2)
    }
    else if (x > (2*i)/n && x < (2*(i+1))/n){
      return (-n/2)
    }
    else {
      return (0)
    }
  }


  # liczenie calki metodą całkowania numerycznego kwadraturą Gaussa-Legendre'a z 3 punktami kwadratury

  calka <- function(f, a, b, i){
    if(i<=-1){
      return(
          (b-a)/2 * (
              (8/9)*f((a+b)/2)+
              (5/9)*f(((b-a)/2)*sqrt((3/5)) + (a+b)/2)+
              (5/9)*f(((b-a)/2)*(-sqrt(3/5)) + (a+b)/2)
          )
      )
    }
    return(
      (b-a)/2 * (
          (8/9)*f((a+b)/2,i)+
          (5/9)*f(((b-a)/2)*sqrt((3/5)) + (a+b)/2,i)+
          (5/9)*f(((b-a)/2)*(-sqrt(3/5)) + (a+b)/2,i)
      )
  )
  }

  # Lewa strona sformulowania wariacyjnego B(u,v)

  L <- function(i,j){
      return(
        calka (
          function(x){
            return(e_p(x, i)*e_p(x, j))
          },
          max(0, (2*(i-1))/n, (2*(j-1))/n), min((2*(i+1))/n, (2*(j+1))/n),-1) - e(0,i)*e(0,j)
        ) # nolint: infix_spaces_linter.
  }

  # Prawa strona sformulowania wariacyjnego L(v)
  P <- function(i){
      return(-20 * e(0, i) + calka(f1, max(0, 2 * (i - 1)/n), min(2, 2 * (i + 1) / n), i)) # nolint
  }

  #funkcja wyliczajaca rownanie transportu ciepla

  effect <- function(){
      #tworzenie macierzy M
      M <- matrix(0, nrow = n, ncol = n)

      for (i in 1:n){
          for (j in 1:n){
              M[i,j] <- L(j-1, i-1)
          }
      }
      # tworzenie macierzy V
      V <- rep(0, n)
      for (i in 1:n){
          V[i] = P(i-1)
      }
      # print(V)
      #funkcja rozwiązuje układ macierzy M * u = V
      u = solve(M, V)
      print(u)
      # return(function(x) {return(x)})
      # kombinacja liniowa funkcji bazowych
      effect_function <- function(x){
          result = 0
          for (i in 1:n){
          result = result + u[i] * e(x, i - 1)
          }
          return(result)
      }
      return(effect_function)
      
  }

  #funkcja rysuje wykres dla równania transportu ciepła
  plot_effect <- function(){
    plot(seq(0, 2, 1/(100*n)), 
        mapply(effect(), seq(0, 2, 1/(100*n))),
        main = 'solution of the heat transfer equation',
        type='l',
        xlab='values x',
        ylab='values y=f(x)',
        cex.lab=1.5,
        cex.axis=1.5,
        cex.main=2)
  }


  #funkcja rysuje wykres funkcji bazowych
  plot_basis <- function(){
    plot(seq(0, 2, 1/(100*n)), mapply(e, seq(0, 2, 1/(100*n)), 1), 
        main = 'basic functions',
        xlab='',
        ylab='',
        type='l',
        cex.axis=1.5,
        cex.main=2)
    for (i in 0:n){
      lines(seq(0, 2, 1/(10*n)), mapply(e, seq(0, 2, 1/(10*n)), i))
    }
  }

  # wywolanie funkcji/narysowania wykresu funkcji bazowych
  if(l_line=='1'){
    plot_basis()
  }
  # wywolanie funkcji/narysowania wykresu rozwiazania rownania transportu ciepla
  if(l_line=='2'){
    plot_effect()
  }
}else {
   cat("number of basis function must be greater than 0\n")
}
