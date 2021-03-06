### Carga de datos
datos <- read.csv("datos.csv")
n <- dim(datos)[1]
y <- datos$poblacion
ln_y <- log(datos$poblacion)
tt <- datos$anio
ln_tt <- log(datos$anio)
uno <- rep.int(1, n)

### Ajuste modelo lineal
# z_lin = b_lin + a_lin * tt

# Defino matriz asociada
A_lin <- matrix(c(uno, tt), nrow = n)

# Sea C = t(A) %*% A,  d = t(A) %*% y x los coeficientes lineales del modelo, tenemos
# C %*% x = d y podemos obtener x como C^-1 %*% d, que se implementa como solve(C,d)
C_lin <- t(A_lin) %*% A_lin
d_lin <- t(A_lin) %*% y
x_lin <- solve(C_lin, d_lin)

# Extraigo los coeficientes
b_lin <- x_lin[1]
a_lin <- x_lin[2]

# Calculo la estimacion lineal de y, z_lin, y el error asociado al modelo
z_lin <- A_lin %*% x_lin
err_lin <- (y-z_lin)^2
et_lin <- sqrt(sum(err_lin))

# Estimo la poblacion segun el modelo lineal para el 2030
est_2030_lin <- b_lin + a_lin * 2030

### Ajuste modelo polinomial
# z_pot = b_pot * t^a_pot
# Que se lineariza como
# ln(z_pot) = ln(b_pot) + a_pot * ln(t)

# Defino matriz asociada
A_pot <- matrix(c(uno, ln_tt), nrow = n)

# Sea C = t(A) %*% A,  d = t(A) %*% y x los coeficientes lineales del modelo, tenemos
# C %*% x = d y podemos obtener x como C^-1 %*% d, que se implementa como solve(C,d)
C_pot <- t(A_pot) %*% A_pot
d_pot <- t(A_pot) %*% ln_y
x_pot <- solve(C_pot, d_pot)

# Extraigo los coeficientes
b_pot <- x_pot[1]
a_pot <- x_pot[2]
# Recupero los coeficientes del modelo no linearizado
b_pot_orig <- exp(b_pot)


# Calculo la estimacion polinomial de y, z_pot, y el error asociado al modelo
z_pot <- b_pot_orig * tt^a_pot
err_pot <- (y-z_pot)^2
et_pot <- sqrt(sum(err_pot))

# Estimo la poblacion segun el modelo polinomial para el 2030
est_2030_pot <- b_pot_orig * 2030^a_pot

### Ajuste modelo exponencial
# z_exp = b_exp * e^(t*a_exp)
# Que se lineariza como
# ln(z_exp) = ln(b_exp) + a_exp * t

# Defino matriz asociada
A_exp <- matrix(c(uno, tt), nrow = n)

# Sea C = t(A) %*% A,  d = t(A) %*% y x los coeficientes lineales del modelo, tenemos
# C %*% x = d y podemos obtener x como C^-1 %*% d, que se implementa como solve(C,d)
C_exp <- t(A_exp) %*% A_exp
d_exp <- t(A_exp) %*% ln_y
x_exp <- solve(C_exp, d_exp)

# Extraigo los coeficientes
b_exp <- x_exp[1]
a_exp <- x_exp[2]
# Recupero los coeficientes del modelo no linearizado
b_exp_orig <- exp(b_exp)


# Calculo la estimacion lineal de y, z_lin, y el error asociado al modelo
z_exp <- b_exp_orig * exp(tt*a_exp)
err_exp <- (y-z_exp)^2
et_exp <- sqrt(sum(err_exp))

# Estimo la poblacion segun el modelo exponencial para el 2030
est_2030_exp <- b_exp_orig * exp(a_exp * 2030)


# Preparo grafico
# Funciones de estimacion
fun_lin <- function(x) { b_lin + a_lin * x}
fun_pot <- function(x) { b_pot_orig * x^a_pot}
fun_exp <- function(x) { b_exp_orig * exp(a_exp*x)}

grafico <- ggplot(datos, aes(x = anio, y = poblacion)) + geom_point()

grafico +
  stat_function(
    fun = fun_lin,
    color = "blue"
  ) +
  stat_function(
    fun = fun_pot,
    color = "red"
  ) +
  stat_function(
    fun = fun_exp,
    color = "green"
  )
