# From Devroye:
# 1. If U ~ Unif(0,pi) and E ~ Exp(1) then
#    ( ( sin(alpha*U) / sin(U) )^{1/(1-alpha)} * sin((1-alpha)*U) / (E * sin(alpha*U) ) )^{(1-alpha)/alpha}
#  = ( sin(alpha*U) / sin(U) )^(1/alpha) * ( sin((1-alpha)*U) / (E * sin(alpha*U) ) )^((1-alpha)/alpha)
# is Stable(alpha, 1)
# 2. If S ~ Stable(alpha/2, 1) and N ~ N(0,1) then
#   N * sqrt(2*S) ~ Stable(alpha, 0)

rstable <- function(n, alpha) {
    U <- runif(n, 0, pi)
    E <- rexp(n)
    N <- rnorm(n)
    a <- alpha/2
    S <- ( sin(a*U) / sin(U) )^(1/a) * ( sin((1-a)*U) / (E * sin(a*U) ) )^((1-a)/a)
    return(N * sqrt(2*S))
}

rstable2 <- function (n, alpha) {
    a <- alpha/2
    U <- runif(n, 0, pi)
    return(
        rnorm(n) * sqrt(2 * 
            ( sin(a*U) / sin(U) )^(1/a) * ( sin((1-a)*U) / (rexp(n) * sin(a*U) ) )^((1-a)/a)
        )
    )
}
