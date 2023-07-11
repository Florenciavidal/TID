
################# TID-2023, Florencia Vidal

rm(list=ls())
###### simulación de transacciones
nweeks=1200
nprod=4
nruns=2000
elast=matrix(c(-2.93,0.15,0.18,0.17,0.38,-1.84,0.03,-0.01,0.89,0.32,-3.72,0.36,1.24,0.92,-0.13,-3.82), nrow=nprod, ncol=nprod)
elast=t(elast)
runelast<-elast
elastsd=matrix(c(1,0.75,0.65,0.39,0.22,0.96,0.3,0.22,0.43,0.98,1.22,0.42,0.76,0.9,0.77,0.85), nrow=nprod, ncol=nprod)

Prices=matrix(0,nrow=nweeks, ncol=nprod)
for (i in 1:nprod ) {

  Prices[,i]=round(rnorm(nweeks, mean=(50-9*i), sd=3))
    }
alpha=c(1.6,0.96,0.929,0.869)
alpha=c(10,8,18,15)
Sales=round(exp(matrix(rnorm(nweeks*nprod,mean=alpha,sd=elastsd), nrow=nweeks, ncol=nprod,byrow = TRUE)+log(Prices)%*%elast))

lp <- as.matrix.data.frame(log(Prices))
colnames(lp)=c("Price1", "Price2", "Price3", "Price4")
colnames(Sales)=c("Sales1", "Sales2", "Sales3", "sales4")
#Estim_elast=solve(t(Prices)%*%Prices)%*%(t(Prices)%*%Sales)
#aux=solve(t(Prices)%*%Prices)
demanda <- cbind(Sales,lp)
demanda <- as.data.frame(demanda)

#View(demanda)

model <- list(1,2,3,4)
for (i in 1:nprod ) {
  aux2=paste("m",i,sep="_")
    assign(aux2,lm(log(demanda[,i])~Price1+Price2+Price3+Price4,data=demanda))
#  assign(model[i],aux2)
}

summary(m_1)


model[[1]] <- m_1
model[[2]] <- m_2
model[[3]] <- m_3
model[[4]] <- m_4

Precio <- Prices[1,]

f_venta <- function(Precio) {
  lp=as.data.frame(t(log(Precio)))
  colnames(lp)=c("Price1", "Price2", "Price3", "Price4")  
    venta <- 0
  for (i in 1:nprod ) {
    v<- Precio[i]*round(exp(predict(model[[i]],data.frame(lp))))
    venta <- venta + v
  }
  return(venta)
}

Precio
f_venta(Precio)

########OPTIMIZACION SIN RESTRICCIONES
Priceopt<-optim(Precio,f_venta,control=list(fnscale=-1))

Priceopt$par
Priceopt$value

############# Restricciones con respecto a precios iniciales
Precioini <- colMeans(Prices)
Precioini

Precio3 <-optim(Precioini,f_venta, lower=Precio*0.8, upper=Precio*1.2,method="L-BFGS-B",control=list(fnscale=-1))
Precio3$par
Precio3$value

Precio3$value/f_venta(Precio)
Precio3$par/Precio

############# Restricciones con respectos a precios mínimos y máximos
preciomin <- apply(Prices, 2, function(x) min(x, na.rm = TRUE))
preciomax <- apply(Prices, 2, function(x) max(x, na.rm = TRUE))

Precio4 <-optim(Precioini,f_venta, lower=preciomin, upper=preciomax,method="L-BFGS-B",control=list(fnscale=-1))
Precio4$par
Precio4$value

Precio4$value/f_venta(Precio)
Precio4$par/Precio

##### Restricciones a precio promedio

ui <- matrix(1,nrow=2,ncol=4)
ui[2,] <- c(-1,-1,-1,-1)
ci <- c(round(mean(rowMeans(Prices))*4)*0.9,-round(mean(rowMeans(Prices))*4)*1.1)

Precio2 <- rep(round(mean(colMeans(Prices))),4)
Priceopt<-constrOptim(Precio2,f_venta,NULL,ui=ui,ci=ci,control=list(fnscale=-1))

Priceopt$par
Priceopt$value

################################################## FLO

########## Aggregando todas las restricciones
u2 <- diag(1,4,4)
u3 <- -u2

uu <- rbind(ui,u2,u3)
cc <- rbind(as.matrix(ci),as.matrix(preciomin),-as.matrix(preciomax))
dim(uu)
dim(cc)

Precioini <- colMeans(Prices)

start_time1 <- Sys.time()
Priceopt1<-constrOptim(Precioini,f_venta,NULL,ui=uu,ci=cc,control=list(fnscale=-1))
end_time1 <- Sys.time()
tiempo1= end_time1 - start_time1

# cuanto se desmora solver 1
tiempo1

# Imprimo los resultados
Priceopt1$par
Priceopt1$value

Priceopt1$value/f_venta(Precio)
Priceopt1$par/Precio

# Simulaciones
f_optim <- function(ui,ci,niter) {
  message=FALSE
  warning=FALSE
  suppressMessages(library(limSolve))
  aux=diag(4)
  aux2=rbind(ui,aux)
  aux3=rbind(as.matrix(ci), as.matrix(rep(0,4)))
  dim(aux2)
  precio1=xsample(G=aux2, H=aux3)
  precio2=precio1$X[nrow(precio1$X),]
  
  precio=0
  pi=0
  iter=0
  pi_values = numeric(niter)
  for(i in 1:niter) {
    precio2=precio1$X[nrow(precio1$X)-i,]
    p_prod <- constrOptim(precio2, f_venta, NULL, ui=ui, ci=ci, control=list(fnscale=-1))
    print(c(i, pi,p_prod$value))
    pi_values[i] = p_prod$value
    if (p_prod$value> pi) {
      pi= p_prod$value
      iter=i
      precio=p_prod$par
      value_over_f_venta = pi / f_venta(Precio)
      par_over_precio = precio / Precio
    } 
  }
  pi_sd = sd(pi_values)
  print(c(iter,pi, pi_sd, precio, value_over_f_venta, par_over_precio))
  result=list(iter=iter, pi=pi, pi_sd=pi_sd, pi_values=pi_values)
  return(result)
}

# Simulaciones
start_time <- Sys.time()
constrOptim1 <- f_optim(uu,cc,100)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
constrOptim2 <- f_optim(uu,cc,400)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
constrOptim3 <- f_optim(uu,cc,700)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
constrOptim4<- f_optim(uu,cc,1000)
end_time <- Sys.time()
end_time - start_time


########## nloptr()
# instalo nloptr()
#install.packages("nloptr")
library(nloptr)

# fo negativa
f_venta_neg <- function(Precio) {
  lp=as.data.frame(t(log(Precio)))
  colnames(lp)=c("Price1", "Price2", "Price3", "Price4")  
  venta <- 0
  for (i in 1:nprod ) {
    v<- Precio[i]*round(exp(predict(model[[i]],data.frame(lp))))
    venta <- venta + v
  }
  return(-venta)
}

# funciones de restricción
eval_g_ineq <- function(Precio) {
  return(as.numeric(uu %*% Precio - cc))
}

# resuelve problema con nloptr
start_time2 <- Sys.time()
Priceopt2 <- nloptr(x0 = Precioini, eval_f = f_venta_neg, eval_g_ineq = eval_g_ineq, opts = list("algorithm" = "NLOPT_LN_COBYLA", "print_level" = 0))
end_time2 <- Sys.time()
tiempo2= end_time2 - start_time2

# cuanto se desmora solver 1
tiempo2

# Imprimo los resultados
Priceopt2$solution
abs(Priceopt2$objective)

Priceopt2$solution / Precio
Priceopt2$objective / f_venta_neg(Precio) 

# Compruebo que las restricciones sean de esta forma restricciones >= 0
restricciones <- uu %*% Priceopt2$solution - cc
print(restricciones >= 0)

# Simulaciones
f_optim <- function(ui,ci,niter) {
  message=FALSE
  warning=FALSE
  suppressMessages(library(limSolve))
  aux=diag(4)
  aux2=rbind(ui,aux)
  aux3=rbind(as.matrix(ci), as.matrix(rep(0,4)))
  dim(aux2)
  precio1=xsample(G=aux2, H=aux3)
  precio2=precio1$X[nrow(precio1$X),]
  
  precio=0
  pi=0
  iter=0
  pi_values = numeric(niter)
  for(i in 1:niter) {
    precio2=precio1$X[nrow(precio1$X)-i,]
    p_prod <- nloptr(
      x0 = precio2,
      eval_f = f_venta_neg,  
      eval_g_ineq = eval_g_ineq,
      opts = list("algorithm" = "NLOPT_LN_COBYLA", "print_level" = 0)
    )
    print(c(i, pi,p_prod$objective))
    pi_values[i] = abs(p_prod$objective)
    if (p_prod$objective < pi) {
      pi= p_prod$objective
      iter=i
      precio=p_prod$solution
      value_over_f_venta = pi / f_venta_neg(Precio)
      par_over_precio = precio / Precio
    } 
  }
  pi_sd = sd(abs(pi_values))
  print(c(iter,pi, pi_sd, precio, value_over_f_venta, par_over_precio))
  result=list(iter=iter, pi=pi, pi_sd=pi_sd, pi_values=pi_values)
  return(result)
}

# Simulaciones
start_time <- Sys.time()
nloptr1 <- f_optim(uu,cc,100)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
nloptr2 <- f_optim(uu,cc,400)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
nloptr3 <- f_optim(uu,cc,700)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
nloptr4<- f_optim(uu,cc,1000)
end_time <- Sys.time()
end_time - start_time




########## Rsolnp
# instalo Rsolnp
# install.packages("Rsolnp")
library(Rsolnp)

# Definir los límites inferiores para las restricciones de desigualdad,  cada restricción de desigualdad debe ser mayor o igual a 0, creo que tambien podría poner -inf enves de 0 porque ya los estableci en restricciones
ineq_lb <- rep(0, 10)
# Definir los límites superiores para las restricciones de desigualdad, puse inf porque los limites superiores ya estan definidos en las restricciones
ineq_ub <- rep(Inf, 10)

# resuelve problema con Rsolnp
start_time3 <- Sys.time()
Priceopt3 <- solnp(pars = Precioini, fun = f_venta_neg, ineqfun = eval_g_ineq, ineqLB = ineq_lb, ineqUB = ineq_ub)
end_time3 <- Sys.time()
tiempo3= end_time3 - start_time3

#cuanto se desmora solver 1
tiempo3

# Imprimo los resultados
Priceopt3$pars
f_venta(Priceopt3$pars) 

Priceopt3$pars / Precio
f_venta(Priceopt3$pars) / f_venta(Precio) 

# Compruebo que las restricciones sean de esta forma restricciones >= 0
restricciones <- uu %*% Priceopt3$pars - cc
print(restricciones >= 0)

# Simulaciones
f_optim <- function(ui,ci,niter) {
  message=FALSE
  warning=FALSE
  suppressMessages(library(limSolve))
  aux=diag(4)
  aux2=rbind(ui,aux)
  aux3=rbind(as.matrix(ci), as.matrix(rep(0,4)))
  dim(aux2)
  precio1=xsample(G=aux2, H=aux3)
  precio2=precio1$X[nrow(precio1$X),]
  
  precio=0
  pi=0
  iter=0
  pi_values = numeric(niter)
  for(i in 1:niter) {
    precio2=precio1$X[nrow(precio1$X)-i,]
    p_prod <- solnp(pars = precio2, fun = f_venta_neg, ineqfun = eval_g_ineq, ineqLB = ineq_lb, ineqUB = ineq_ub)
    print(c(i, pi,f_venta_neg(p_prod$par)))
    pi_values[i] = abs(f_venta_neg(p_prod$par))
    if (f_venta_neg(p_prod$par) < pi) {
      pi= f_venta_neg(p_prod$par)
      iter=i
      precio=p_prod$par
      value_over_f_venta = pi / f_venta_neg(Precio)
      par_over_precio = precio / Precio
    } 
  }
  pi_sd = sd(abs(pi_values))
  print(c(iter,pi, pi_sd, precio, value_over_f_venta, par_over_precio))
  result=list(iter=iter, pi=pi, pi_sd=pi_sd, pi_values=pi_values)
  return(result)
}

# Simulaciones
start_time <- Sys.time()
solnp1 <- f_optim(uu,cc,100)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
solnp2 <- f_optim(uu,cc,400)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
solnp3 <- f_optim(uu,cc,700)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
solnp4<- f_optim(uu,cc,1000)
end_time <- Sys.time()
end_time - start_time



########## optimx
#install.packages("optimx")
library(optimx)

# Define la función objetivo con restricciones
f_venta_res <- function(Precio) {
  # Evalúa las restricciones
  restricciones <- eval_g_ineq(Precio)
  
  # Si alguna restricción es violada, devuelve -Infinito
  if (any(restricciones < 0)) {
    return(-100)
  }
  
  # Si no se violan las restricciones, calcula y devuelve el valor de la función objetivo
  lp = as.data.frame(t(log(Precio)))
  colnames(lp) <-  c("Price1", "Price2", "Price3", "Price4")  
  venta = 0
  for (i in 1:nprod ) {
    v <- Precio[i] * round(exp(predict(model[[i]], data.frame(lp))))
    venta <- venta + v
  }
  return(venta)
}

# Utiliza optimx para resolver el problema de optimización L-BFGS-B
start_time4 <- Sys.time()
Priceopt4 <- optimx(Precioini, f_venta_res, gr = NULL, method = "BFGS", 
                 control = list(maxit = 1000,maximize = TRUE))
end_time4 <- Sys.time()
tiempo4= end_time4 - start_time4
#cuanto se desmora solver 1
tiempo4

#### convierto 
# Extrae los valores de las variables en el óptimo
precios_optimos <- coef(Priceopt4)

# Formatea los valores con el formato deseado
convertido1 <- format(precios_optimos, nsmall = 4)

# Convierte los valores formateados en un vector numérico
convertido2 <- as.numeric(convertido1)

# Imprimo los resultados
convertido2
f_venta_res(convertido2)

convertido2 / Precio
f_venta_res(convertido2) / f_venta(Precio) 

# Simulaciones
f_optim <- function(ui,ci,niter) {
  message=FALSE
  warning=FALSE
  suppressMessages(library(limSolve))
  aux=diag(4)
  aux2=rbind(ui,aux)
  aux3=rbind(as.matrix(ci), as.matrix(rep(0,4)))
  dim(aux2)
  precio1=xsample(G=aux2, H=aux3)
  precio2=precio1$X[nrow(precio1$X),]
  
  precio=0
  pi=0
  iter=0
  pi_values = numeric(niter)
  for(i in 1:niter) {
    precio2=precio1$X[nrow(precio1$X)-i,]
    p_prod <- optimx(precio2, f_venta_res, gr = NULL, method = "BFGS", 
                     control = list(maxit = 1000, maximize = TRUE))
    precios_optimos <- coef(p_prod)
    convertido1 <- format(precios_optimos, nsmall = 4)
    convertido2 <- as.numeric(convertido1)
    print(c(i, pi,f_venta_res(convertido2)))
    pi_values[i] = f_venta_res(convertido2)
    if (f_venta_res(convertido2) > pi) {
      pi= f_venta_res(convertido2)
      iter=i
      precio=convertido2
      value_over_f_venta = pi / f_venta_res(Precio)
      par_over_precio = precio / Precio
    } 
  }
  pi_sd = sd(pi_values)
  print(c(iter,pi, pi_sd, precio, value_over_f_venta, par_over_precio))
  result=list(iter=iter, pi=pi, pi_sd=pi_sd, pi_values=pi_values)
  return(result)
}

# Simulaciones
start_time <- Sys.time()
optimx1 <- f_optim(uu,cc,100)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx2 <- f_optim(uu,cc,400)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx3 <- f_optim(uu,cc,700)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx4<- f_optim(uu,cc,1000)
end_time <- Sys.time()
end_time - start_time

# Utiliza optimx con algoritmo Nelder-Mead para resolver el problema de optimización 
start_time4 <- Sys.time()
Priceopt4 <- optimx(Precioini, f_venta_res, gr = NULL, method = "Nelder-Mead", 
                    control = list(maxit = 1000,maximize = TRUE))
end_time4 <- Sys.time()
tiempo4= end_time4 - start_time4

#cuanto se desmora solver 1
tiempo4

#### convierto 
# Extrae los valores de las variables en el óptimo
precios_optimos <- coef(Priceopt4)

# Formatea los valores con el formato deseado
convertido1 <- format(precios_optimos, nsmall = 4)

# Convierte los valores formateados en un vector numérico
convertido2 <- as.numeric(convertido1)

# Imprimo los resultados
convertido2
f_venta_res(convertido2)

convertido2 / Precio
f_venta_res(convertido2) / f_venta(Precio) 

# Simulaciones
f_optim <- function(ui,ci,niter) {
  message=FALSE
  warning=FALSE
  suppressMessages(library(limSolve))
  aux=diag(4)
  aux2=rbind(ui,aux)
  aux3=rbind(as.matrix(ci), as.matrix(rep(0,4)))
  dim(aux2)
  precio1=xsample(G=aux2, H=aux3)
  precio2=precio1$X[nrow(precio1$X),]
  
  precio=0
  pi=0
  iter=0
  pi_values = numeric(niter)
  for(i in 1:niter) {
    precio2=precio1$X[nrow(precio1$X)-i,]
    p_prod <- optimx(precio2, f_venta_res, gr = NULL, method = "BFGS", 
                     control = list(maxit = 1000, maximize = TRUE))
    precios_optimos <- coef(p_prod)
    convertido1 <- format(precios_optimos, nsmall = 4)
    convertido2 <- as.numeric(convertido1)
    print(c(i, pi,f_venta_res(convertido2)))
    pi_values[i] = f_venta_res(convertido2)
    if (f_venta_res(convertido2) > pi) {
      pi= f_venta_res(convertido2)
      iter=i
      precio=convertido2
      value_over_f_venta = pi / f_venta_res(Precio)
      par_over_precio = precio / Precio
    } 
  }
  pi_sd = sd(pi_values)
  print(c(iter,pi, pi_sd, precio, value_over_f_venta, par_over_precio))
  result=list(iter=iter, pi=pi, pi_sd=pi_sd, pi_values=pi_values)
  return(result)
}

# Simulaciones
start_time <- Sys.time()
optimx1 <- f_optim(uu,cc,100)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx2 <- f_optim(uu,cc,400)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx3 <- f_optim(uu,cc,700)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx4<- f_optim(uu,cc,1000)
end_time <- Sys.time()
end_time - start_time



# Utiliza optimx con algoritmo L-BFGS-B  para resolver el problema de optimización
start_time4 <- Sys.time()
Priceopt4 <- optimx(Precioini, f_venta_res, gr = NULL, method = "Nelder-Mead", 
                    control = list(maxit = 1000, maxit = 1000, maximize = TRUE))
end_time4 <- Sys.time()
tiempo4= end_time4 - start_time4
Priceopt4
#cuanto se desmora solver 1
tiempo4

#### convierto 
# Extrae los valores de las variables en el óptimo
precios_optimos <- coef(Priceopt4)

# Formatea los valores con el formato deseado
convertido1 <- format(precios_optimos, nsmall = 4)

# Convierte los valores formateados en un vector numérico
convertido2 <- as.numeric(convertido1)

# Imprimo los resultados
convertido2
f_venta_res(convertido2)

convertido2 / Precio
f_venta_res(convertido2) / f_venta(Precio) 

# Simulaciones
f_optim <- function(ui,ci,niter) {
  message=FALSE
  warning=FALSE
  suppressMessages(library(limSolve))
  aux=diag(4)
  aux2=rbind(ui,aux)
  aux3=rbind(as.matrix(ci), as.matrix(rep(0,4)))
  dim(aux2)
  precio1=xsample(G=aux2, H=aux3)
  precio2=precio1$X[nrow(precio1$X),]
  
  precio=0
  pi=0
  iter=0
  pi_values = numeric(niter)
  for(i in 1:niter) {
    precio2=precio1$X[nrow(precio1$X)-i,]
    p_prod <- optimx(precio2, f_venta_res, gr = NULL, method = "BFGS", 
                     control = list(maxit = 1000, maximize = TRUE))
    precios_optimos <- coef(p_prod)
    convertido1 <- format(precios_optimos, nsmall = 4)
    convertido2 <- as.numeric(convertido1)
    print(c(i, pi,f_venta_res(convertido2)))
    pi_values[i] = f_venta_res(convertido2)
    if (f_venta_res(convertido2) > pi) {
      pi= f_venta_res(convertido2)
      iter=i
      precio=convertido2
      value_over_f_venta = pi / f_venta_res(Precio)
      par_over_precio = precio / Precio
    } 
  }
  pi_sd = sd(pi_values)
  print(c(iter,pi, pi_sd, precio, value_over_f_venta, par_over_precio))
  result=list(iter=iter, pi=pi, pi_sd=pi_sd, pi_values=pi_values)
  return(result)
}

# Simulaciones
start_time <- Sys.time()
optimx1 <- f_optim(uu,cc,100)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx2 <- f_optim(uu,cc,400)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx3 <- f_optim(uu,cc,700)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx4<- f_optim(uu,cc,1000)
end_time <- Sys.time()
end_time - start_time



# Utiliza optimx con algoritmo nlm  para resolver el problema de optimización
start_time4 <- Sys.time()
Priceopt4 <- optimx(Precioini, f_venta_res, gr = NULL, method = "nlm", 
                    control = list(maxit = 1000, maximize = TRUE))
end_time4 <- Sys.time()
tiempo4= end_time4 - start_time4
Priceopt4
#cuanto se desmora solver 1
tiempo4

#### convierto 
# Extrae los valores de las variables en el óptimo
precios_optimos <- coef(Priceopt4)

# Formatea los valores con el formato deseado
convertido1 <- format(precios_optimos, nsmall = 4)

# Convierte los valores formateados en un vector numérico
convertido2 <- as.numeric(convertido1)

# Imprimo los resultados
convertido2
f_venta_res(convertido2)

convertido2 / Precio
f_venta_res(convertido2) / f_venta(Precio) 

# Simulaciones
f_optim <- function(ui,ci,niter) {
  message=FALSE
  warning=FALSE
  suppressMessages(library(limSolve))
  aux=diag(4)
  aux2=rbind(ui,aux)
  aux3=rbind(as.matrix(ci), as.matrix(rep(0,4)))
  dim(aux2)
  precio1=xsample(G=aux2, H=aux3)
  precio2=precio1$X[nrow(precio1$X),]
  
  precio=0
  pi=0
  iter=0
  pi_values = numeric(niter)
  for(i in 1:niter) {
    precio2=precio1$X[nrow(precio1$X)-i,]
    p_prod <- optimx(precio2, f_venta_res, gr = NULL, method = "BFGS", 
                     control = list(maxit = 1000, maximize = TRUE))
    precios_optimos <- coef(p_prod)
    convertido1 <- format(precios_optimos, nsmall = 4)
    convertido2 <- as.numeric(convertido1)
    print(c(i, pi,f_venta_res(convertido2)))
    pi_values[i] = f_venta_res(convertido2)
    if (f_venta_res(convertido2) > pi) {
      pi= f_venta_res(convertido2)
      iter=i
      precio=convertido2
      value_over_f_venta = pi / f_venta_res(Precio)
      par_over_precio = precio / Precio
    } 
  }
  pi_sd = sd(pi_values)
  print(c(iter,pi, pi_sd, precio, value_over_f_venta, par_over_precio))
  result=list(iter=iter, pi=pi, pi_sd=pi_sd, pi_values=pi_values)
  return(result)
}

# Simulaciones
start_time <- Sys.time()
optimx1 <- f_optim(uu,cc,100)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx2 <- f_optim(uu,cc,400)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx3 <- f_optim(uu,cc,700)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx4<- f_optim(uu,cc,1000)
end_time <- Sys.time()
end_time - start_time


# Utiliza optimx con algoritmo nlminb  para resolver el problema de optimización
start_time4 <- Sys.time()
Priceopt4 <- optimx(Precioini, f_venta_res, gr = NULL, method = "nlminb", 
                    control = list(maxit = 1000, maximize = TRUE))
end_time4 <- Sys.time()
tiempo4= end_time4 - start_time4
Priceopt4
#cuanto se desmora solver 1
tiempo4

#### convierto 
# Extrae los valores de las variables en el óptimo
precios_optimos <- coef(Priceopt4)

# Formatea los valores con el formato deseado
convertido1 <- format(precios_optimos, nsmall = 4)

# Convierte los valores formateados en un vector numérico
convertido2 <- as.numeric(convertido1)

# Imprimo los resultados
convertido2
f_venta_res(convertido2)

convertido2 / Precio
f_venta_res(convertido2) / f_venta(Precio) 

# Simulaciones
f_optim <- function(ui,ci,niter) {
  message=FALSE
  warning=FALSE
  suppressMessages(library(limSolve))
  aux=diag(4)
  aux2=rbind(ui,aux)
  aux3=rbind(as.matrix(ci), as.matrix(rep(0,4)))
  dim(aux2)
  precio1=xsample(G=aux2, H=aux3)
  precio2=precio1$X[nrow(precio1$X),]
  
  precio=0
  pi=0
  iter=0
  pi_values = numeric(niter)
  for(i in 1:niter) {
    precio2=precio1$X[nrow(precio1$X)-i,]
    p_prod <- optimx(precio2, f_venta_res, gr = NULL, method = "BFGS", 
                     control = list(maxit = 1000, maximize = TRUE))
    precios_optimos <- coef(p_prod)
    convertido1 <- format(precios_optimos, nsmall = 4)
    convertido2 <- as.numeric(convertido1)
    print(c(i, pi,f_venta_res(convertido2)))
    pi_values[i] = f_venta_res(convertido2)
    if (f_venta_res(convertido2) > pi) {
      pi= f_venta_res(convertido2)
      iter=i
      precio=convertido2
      value_over_f_venta = pi / f_venta_res(Precio)
      par_over_precio = precio / Precio
    } 
  }
  pi_sd = sd(pi_values)
  print(c(iter,pi, pi_sd, precio, value_over_f_venta, par_over_precio))
  result=list(iter=iter, pi=pi, pi_sd=pi_sd, pi_values=pi_values)
  return(result)
}

# Simulaciones
start_time <- Sys.time()
optimx1 <- f_optim(uu,cc,100)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx2 <- f_optim(uu,cc,400)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx3 <- f_optim(uu,cc,700)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx4<- f_optim(uu,cc,1000)
end_time <- Sys.time()
end_time - start_time


# Utiliza optimx con algoritmo Rcgmin para resolver el problema de optimización
start_time4 <- Sys.time()
Priceopt4 <- optimx(Precioini, f_venta_res, gr = NULL, method = "Rcgmin", 
                    control = list(maxit = 1000, maximize = TRUE))
end_time4 <- Sys.time()
tiempo4= end_time4 - start_time4
Priceopt4
#cuanto se desmora solver 1
tiempo4

#### convierto 
# Extrae los valores de las variables en el óptimo
precios_optimos <- coef(Priceopt4)

# Formatea los valores con el formato deseado
convertido1 <- format(precios_optimos, nsmall = 4)

# Convierte los valores formateados en un vector numérico
convertido2 <- as.numeric(convertido1)

# Imprimo los resultados
convertido2
f_venta_res(convertido2)

convertido2 / Precio
f_venta_res(convertido2) / f_venta(Precio) 

# Simulaciones
f_optim <- function(ui,ci,niter) {
  message=FALSE
  warning=FALSE
  suppressMessages(library(limSolve))
  aux=diag(4)
  aux2=rbind(ui,aux)
  aux3=rbind(as.matrix(ci), as.matrix(rep(0,4)))
  dim(aux2)
  precio1=xsample(G=aux2, H=aux3)
  precio2=precio1$X[nrow(precio1$X),]
  
  
  precio=0
  pi=0
  iter=0
  pi_values = numeric(niter)
  for(i in 1:niter) {
    precio2=precio1$X[nrow(precio1$X)-i,]
    p_prod <- optimx(precio2, f_venta_res, gr = NULL, method = "BFGS", 
                     control = list(maxit = 1000, maximize = TRUE))
    precios_optimos <- coef(p_prod)
    convertido1 <- format(precios_optimos, nsmall = 4)
    convertido2 <- as.numeric(convertido1)
    print(c(i, pi,f_venta_res(convertido2)))
    pi_values[i] = f_venta_res(convertido2)
    if (f_venta_res(convertido2) > pi) {
      pi= f_venta_res(convertido2)
      iter=i
      precio=convertido2
      value_over_f_venta = pi / f_venta_res(Precio)
      par_over_precio = precio / Precio
    } 
  }
  pi_sd = sd(pi_values)
  print(c(iter,pi, pi_sd, precio, value_over_f_venta, par_over_precio))
  result=list(iter=iter, pi=pi, pi_sd=pi_sd, pi_values=pi_values)
  return(result)
}

# Simulaciones
start_time <- Sys.time()
optimx1 <- f_optim(uu,cc,100)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx2 <- f_optim(uu,cc,400)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx3 <- f_optim(uu,cc,700)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx4<- f_optim(uu,cc,1000)
end_time <- Sys.time()
end_time - start_time


# Utiliza optimx con algoritmo Rvmmin para resolver el problema de optimización
start_time4 <- Sys.time()
Priceopt4 <- optimx(Precioini, f_venta_res, gr = NULL, method = "Rvmmin", 
                    control = list(maxit = 1000, maximize = TRUE))
end_time4 <- Sys.time()
tiempo4= end_time4 - start_time4

#cuanto se desmora solver 1
tiempo4

#### convierto 
# Extrae los valores de las variables en el óptimo
precios_optimos <- coef(Priceopt4)

# Formatea los valores con el formato deseado
convertido1 <- format(precios_optimos, nsmall = 4)

# Convierte los valores formateados en un vector numérico
convertido2 <- as.numeric(convertido1)

# Imprimo los resultados
convertido2
f_venta_res(convertido2)

convertido2 / Precio
f_venta_res(convertido2) / f_venta(Precio) 

# Simulaciones
f_optim <- function(ui,ci,niter) {
  message=FALSE
  warning=FALSE
  suppressMessages(library(limSolve))
  aux=diag(4)
  aux2=rbind(ui,aux)
  aux3=rbind(as.matrix(ci), as.matrix(rep(0,4)))
  dim(aux2)
  precio1=xsample(G=aux2, H=aux3)
  precio2=precio1$X[nrow(precio1$X),]
  
  precio=0
  pi=0
  iter=0
  pi_values = numeric(niter)
  for(i in 1:niter) {
    precio2=precio1$X[nrow(precio1$X)-i,]
    p_prod <- optimx(precio2, f_venta_res, gr = NULL, method = "BFGS", 
                     control = list(maxit = 1000, maximize = TRUE))
    precios_optimos <- coef(p_prod)
    convertido1 <- format(precios_optimos, nsmall = 4)
    convertido2 <- as.numeric(convertido1)
    print(c(i, pi,f_venta_res(convertido2)))
    pi_values[i] = f_venta_res(convertido2)
    if (f_venta_res(convertido2) > pi) {
      pi= f_venta_res(convertido2)
      iter=i
      precio=convertido2
      value_over_f_venta = pi / f_venta_res(Precio)
      par_over_precio = precio / Precio
    } 
  }
  pi_sd = sd(pi_values)
  print(c(iter,pi, pi_sd, precio, value_over_f_venta, par_over_precio))
  result=list(iter=iter, pi=pi, pi_sd=pi_sd, pi_values=pi_values)
  return(result)
}

# Simulaciones
start_time <- Sys.time()
optimx1 <- f_optim(uu,cc,100)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx2 <- f_optim(uu,cc,400)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx3 <- f_optim(uu,cc,700)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
optimx4<- f_optim(uu,cc,1000)
end_time <- Sys.time()
end_time - start_time


########## nlminb
# Función objetivo negativa con restricciones
f_venta_neg_res <- function(Precio) {
  # Evalúa las restricciones
  restricciones <- eval_g_ineq(Precio)
  
  # Si alguna restricción es violada, devuelve Infinito
  if (any(restricciones < 0)) {
    return(100)
  }
  
  # Si no se violan las restricciones, calcula y devuelve el valor de la función objetivo
  lp = as.data.frame(t(log(Precio)))
  colnames(lp) <-  c("Price1", "Price2", "Price3", "Price4")  
  venta = 0
  for (i in 1:nprod ) {
    v <- Precio[i] * round(exp(predict(model[[i]], data.frame(lp))))
    venta <- venta + v
  }
  return(-venta)
}

# resuelve problema con nlminb
start_time5 <- Sys.time()
Priceopt5 <- nlminb(start = Precioini, objective = f_venta_neg_res, gradient = NULL)
end_time5 <- Sys.time()
tiempo5= end_time5 - start_time5

#cuanto se desmora solver 1
tiempo5

# Imprimo los resultados # Priceopt6$objective
Priceopt5$par
abs(f_venta_neg_res(Priceopt5$par))

Priceopt5$par / Precio
f_venta_neg_res(Priceopt5$par) / f_venta_neg(Precio) 

# Simulaciones
f_optim <- function(ui,ci,niter) {
  message=FALSE
  warning=FALSE
  suppressMessages(library(limSolve))
  aux=diag(4)
  aux2=rbind(ui,aux)
  aux3=rbind(as.matrix(ci), as.matrix(rep(0,4)))
  dim(aux2)
  precio1=xsample(G=aux2, H=aux3)
  precio2=precio1$X[nrow(precio1$X),]
  
  precio=0
  pi=0
  iter=0
  pi_values = numeric(niter)
  for(i in 1:niter) {
    precio2=precio1$X[nrow(precio1$X)-i,]
    p_prod <- nlminb(start = precio2 , objective = f_venta_neg_res, gradient = NULL)
    print(c(i, pi,f_venta_neg_res(p_prod$par)))
    pi_values[i] =abs(f_venta_neg_res(p_prod$par))
    if (f_venta_neg_res(p_prod$par) < pi) {
      pi= f_venta_neg_res(p_prod$par)
      iter=i
      precio=p_prod$par
      value_over_f_venta = pi / f_venta_neg(Precio)
      par_over_precio = precio / Precio
      
    } 
  }
  pi_sd = sd(abs(pi_values))
  print(c(iter,pi, pi_sd, precio, value_over_f_venta, par_over_precio))
  result=list(iter=iter, pi=pi, pi_sd=pi_sd, pi_values=pi_values)
  return(result)
}

# Simulaciones
start_time <- Sys.time()
nlminb1 <- f_optim(uu,cc,100)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
nlminb2 <- f_optim(uu,cc,400)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
nlminb3 <- f_optim(uu,cc,700)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
nlminb4<- f_optim(uu,cc,1000)
end_time <- Sys.time()
end_time - start_time


##################################### gráficos de barra

# install.packages("ggplot2")
# Cargar la biblioteca necesaria
library(ggplot2)

# 1 iteración grafico de barra
# Crear un marco de datos con los nombres de los solvers y las ganancias correspondientes
datos <- data.frame(
  Solver = c('constrOptim', 'nloptr: COBYLA', 'Rsolnp', 'optimx: BFGS', 'optimx: Nelder-Mead', 'optimx: L-BFGS-B', 'optimx: nlm', 'optimx: nlminb', 'optimx: Rcgmin', 'optimx: Rvmmin', 'nlminb'),
  Ganancia = c(142659, 71648, 75285, 111860, 142441, 142441, 66381, 66378, 66381, 66390, 66401)
)

# Convertir las ganancias a miles de pesos
datos$Ganancia <- datos$Ganancia / 1000

# Agregar una nueva columna para el color de las barras
# Los tres solvers con las mayores ganancias tendrán colores distintos. Los demás serán de color gris claro.
datos$Color <- ifelse(rank(-datos$Ganancia) <= 3, c("#FFD700", "#C0C0C0", "#CD7F32")[rank(-datos$Ganancia)], "skyblue")

# Agregar una nueva columna para el ranking
datos$Ranking <- rank(-datos$Ganancia)

# Definir una función para ajustar los tamaños de la fuente
adjust_font_size <- function(ranking) {
  return(7 - log(ranking))
}

# Reordenar los factores de la columna Solver según la Ganancia
datos$Solver <- factor(datos$Solver, levels = datos$Solver[order(datos$Ganancia, decreasing = TRUE)])

# Ajustar el tamaño de la fuente utilizando la función definida
datos$FontSize <- sapply(datos$Ranking, adjust_font_size)

# Generar el gráfico de barras
ggplot(datos, aes(x = Solver, y = Ganancia, fill = Color)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Ranking, "º"), size = FontSize), vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"), # números del eje y en negro
        plot.title = element_text(size = 20, hjust = 0.5, margin = margin(t = 20, b = 20), face = "bold"), # título del gráfico en negrita
        legend.position = "none",
        axis.title.x = element_text(face = "bold"), # título del eje x en negrita
        axis.title.y = element_text(face = "bold")) + # título del eje y en negrita
  ylab("Ganancia (en miles de pesos)") +
  xlab("Solver") +
  ggtitle("Ganancia máxima por solver para 1 iteración") +
  coord_cartesian(ylim = c(60, 145)) + # Ajusta los límites según lo necesites
  scale_fill_identity()

# 100 iteraciones grafico de barra
# Crear un marco de datos con los nombres de los solvers y las ganancias correspondientes
datos <- data.frame(
  Solver = c('constrOptim', 'nloptr: COBYLA', 'Rsolnp', 'optimx: BFGS', 'optimx: Nelder-Mead', 'optimx: L-BFGS-B', 'optimx: nlm', 'optimx: nlminb', 'optimx: Rcgmin', 'optimx: Rvmmin', 'nlminb'),
  Ganancia = c(1156006, 71691, 1035574, 1137108, 1140981, 1138870, 1146600, 1129985, 1124319, 1141542, 1153153)
)

# Convertir las ganancias a miles de pesos
datos$Ganancia <- datos$Ganancia / 1000

# Agregar una nueva columna para el color de las barras
# Los tres solvers con las mayores ganancias tendrán colores distintos. Los demás serán de color gris claro.
datos$Color <- ifelse(rank(-datos$Ganancia) <= 3, c("#FFD700", "#C0C0C0", "#CD7F32")[rank(-datos$Ganancia)], "skyblue")

# Agregar una nueva columna para el ranking
datos$Ranking <- rank(-datos$Ganancia)

# Definir una función para ajustar los tamaños de la fuente
adjust_font_size <- function(ranking) {
  return(7 - log(ranking))
}

# Reordenar los factores de la columna Solver según la Ganancia
datos$Solver <- factor(datos$Solver, levels = datos$Solver[order(datos$Ganancia, decreasing = TRUE)])

# Ajustar el tamaño de la fuente utilizando la función definida
datos$FontSize <- sapply(datos$Ranking, adjust_font_size)

# Generar el gráfico de barras
ggplot(datos, aes(x = Solver, y = Ganancia, fill = Color)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Ranking, "º"), size = FontSize), vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"), # números del eje y en negro
        plot.title = element_text(size = 20, hjust = 0.5, margin = margin(t = 20, b = 20), face = "bold"), # título del gráfico en negrita
        legend.position = "none",
        axis.title.x = element_text(face = "bold"), # título del eje x en negrita
        axis.title.y = element_text(face = "bold")) + # título del eje y en negrita
  ylab("Ganancia (en miles de pesos)") +
  xlab("Solver") +
  ggtitle("Ganancia máxima por solver para 100 iteraciones") +
  coord_cartesian(ylim = c(100, 1200)) + # Ajusta los límites según lo necesites
  scale_fill_identity()

# 400 iteraciones grafico de barra
# Crear un marco de datos con los nombres de los solvers y las ganancias correspondientes
datos <- data.frame(
  Solver = c('constrOptim', 'nloptr: COBYLA', 'Rsolnp', 'optimx: BFGS', 'optimx: Nelder-Mead', 'optimx: L-BFGS-B', 'optimx: nlm', 'optimx: nlminb', 'optimx: Rcgmin', 'optimx: Rvmmin', 'nlminb'),
  Ganancia = c(1156008, 71693, 1056482, 1149065, 1154298, 1150344, 1149840, 1142783, 1153088, 1138106, 1153153)
)

# Convertir las ganancias a miles de pesos
datos$Ganancia <- datos$Ganancia / 1000

# Agregar una nueva columna para el color de las barras
# Los tres solvers con las mayores ganancias tendrán colores distintos. Los demás serán de color gris claro.
datos$Color <- ifelse(rank(-datos$Ganancia) <= 3, c("#FFD700", "#C0C0C0", "#CD7F32")[rank(-datos$Ganancia)], "skyblue")

# Agregar una nueva columna para el ranking
datos$Ranking <- rank(-datos$Ganancia)

# Definir una función para ajustar los tamaños de la fuente
adjust_font_size <- function(ranking) {
  return(7 - log(ranking))
}

# Reordenar los factores de la columna Solver según la Ganancia
datos$Solver <- factor(datos$Solver, levels = datos$Solver[order(datos$Ganancia, decreasing = TRUE)])

# Ajustar el tamaño de la fuente utilizando la función definida
datos$FontSize <- sapply(datos$Ranking, adjust_font_size)

# Generar el gráfico de barras
ggplot(datos, aes(x = Solver, y = Ganancia, fill = Color)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Ranking, "º"), size = FontSize), vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"), # números del eje y en negro
        plot.title = element_text(size = 20, hjust = 0.5, margin = margin(t = 20, b = 20), face = "bold"), # título del gráfico en negrita
        legend.position = "none",
        axis.title.x = element_text(face = "bold"), # título del eje x en negrita
        axis.title.y = element_text(face = "bold")) + # título del eje y en negrita
  ylab("Ganancia (en miles de pesos)") +
  xlab("Solver") +
  ggtitle("Ganancia máxima por solver para 400 iteraciones") +
  coord_cartesian(ylim = c(100, 1200)) + # Ajusta los límites según lo necesites
  scale_fill_identity()

# 700 iteraciones grafico de barra
# Crear un marco de datos con los nombres de los solvers y las ganancias correspondientes
datos <- data.frame(
  Solver = c('constrOptim', 'nloptr: COBYLA', 'Rsolnp', 'optimx: BFGS', 'optimx: Nelder-Mead', 'optimx: L-BFGS-B', 'optimx: nlm', 'optimx: nlminb', 'optimx: Rcgmin', 'optimx: Rvmmin', 'nlminb'),
  Ganancia = c(1156009, 71693, 1076736, 1150898, 1154473, 1153610, 1153649, 1148501, 1147933, 1152093, 924990)
)

# Convertir las ganancias a miles de pesos
datos$Ganancia <- datos$Ganancia / 1000

# Agregar una nueva columna para el color de las barras
# Los tres solvers con las mayores ganancias tendrán colores distintos. Los demás serán de color gris claro.
datos$Color <- ifelse(rank(-datos$Ganancia) <= 3, c("#FFD700", "#C0C0C0", "#CD7F32")[rank(-datos$Ganancia)], "skyblue")

# Agregar una nueva columna para el ranking
datos$Ranking <- rank(-datos$Ganancia)

# Definir una función para ajustar los tamaños de la fuente
adjust_font_size <- function(ranking) {
  return(7 - log(ranking))
}

# Reordenar los factores de la columna Solver según la Ganancia
datos$Solver <- factor(datos$Solver, levels = datos$Solver[order(datos$Ganancia, decreasing = TRUE)])

# Ajustar el tamaño de la fuente utilizando la función definida
datos$FontSize <- sapply(datos$Ranking, adjust_font_size)

# Generar el gráfico de barras
ggplot(datos, aes(x = Solver, y = Ganancia, fill = Color)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Ranking, "º"), size = FontSize), vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"), # números del eje y en negro
        plot.title = element_text(size = 20, hjust = 0.5, margin = margin(t = 20, b = 20), face = "bold"), # título del gráfico en negrita
        legend.position = "none",
        axis.title.x = element_text(face = "bold"), # título del eje x en negrita
        axis.title.y = element_text(face = "bold")) + # título del eje y en negrita
  ylab("Ganancia (en miles de pesos)") +
  xlab("Solver") +
  ggtitle("Ganancia máxima por solver para 700 iteraciones") +
  coord_cartesian(ylim = c(100, 1200)) + # Ajusta los límites según lo necesites
  scale_fill_identity()

# 1000 iteraciones grafico de barra
# Crear un marco de datos con los nombres de los solvers y las ganancias correspondientes
datos <- data.frame(
  Solver = c('constrOptim', 'nloptr: COBYLA', 'Rsolnp', 'optimx: BFGS', 'optimx: Nelder-Mead', 'optimx: L-BFGS-B', 'optimx: nlm', 'optimx: nlminb', 'optimx: Rcgmin', 'optimx: Rvmmin', 'nlminb'),
  Ganancia = c(1156009, 71692, 1064694, 1153153, 1149248, 1148660, 1153649, 1153510, 1149831, 1151376, 1002168)
)

# Convertir las ganancias a miles de pesos
datos$Ganancia <- datos$Ganancia / 1000

# Agregar una nueva columna para el color de las barras
# Los tres solvers con las mayores ganancias tendrán colores distintos. Los demás serán de color gris claro.
datos$Color <- ifelse(rank(-datos$Ganancia) <= 3, c("#FFD700", "#C0C0C0", "#CD7F32")[rank(-datos$Ganancia)], "skyblue")

# Agregar una nueva columna para el ranking
datos$Ranking <- rank(-datos$Ganancia)

# Definir una función para ajustar los tamaños de la fuente
adjust_font_size <- function(ranking) {
  return(7 - log(ranking))
}

# Reordenar los factores de la columna Solver según la Ganancia
datos$Solver <- factor(datos$Solver, levels = datos$Solver[order(datos$Ganancia, decreasing = TRUE)])

# Ajustar el tamaño de la fuente utilizando la función definida
datos$FontSize <- sapply(datos$Ranking, adjust_font_size)

# Generar el gráfico de barras
ggplot(datos, aes(x = Solver, y = Ganancia, fill = Color)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Ranking, "º"), size = FontSize), vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"), # números del eje y en negro
        plot.title = element_text(size = 20, hjust = 0.5, margin = margin(t = 20, b = 20), face = "bold"), # título del gráfico en negrita
        legend.position = "none",
        axis.title.x = element_text(face = "bold"), # título del eje x en negrita
        axis.title.y = element_text(face = "bold")) + # título del eje y en negrita
  ylab("Ganancia (en miles de pesos)") +
  xlab("Solver") +
  ggtitle("Ganancia máxima por solver para 1000 iteraciones") +
  coord_cartesian(ylim = c(100, 1158)) + # Ajusta los límites según lo necesites
  scale_fill_identity()

# La ganancia más alta encontrada por cada solver
# Crear un marco de datos con los nombres de los solvers y las ganancias correspondientes
datos <- data.frame(
  Solver = c('constrOptim', 'nloptr: COBYLA', 'Rsolnp', 'optimx: BFGS', 'optimx: Nelder-Mead', 'optimx: L-BFGS-B', 'optimx: nlm', 'optimx: nlminb', 'optimx: Rcgmin', 'optimx: Rvmmin', 'nlminb'),
  Ganancia = c(1156009, 71692, 1076736, 1153153, 1154473, 1153610, 1153649, 1153510, 1153088, 1152093, 1153153)
)

# Convertir las ganancias a miles de pesos
datos$Ganancia <- datos$Ganancia / 1000

# Agregar una nueva columna para el color de las barras
# Los tres solvers con las mayores ganancias tendrán colores distintos. Los demás serán de color gris claro.
datos$Color <- ifelse(rank(-datos$Ganancia) <= 3, c("#FFD700", "#C0C0C0", "#CD7F32")[rank(-datos$Ganancia)], "skyblue")

# Agregar una nueva columna para el ranking
datos$Ranking <- rank(-datos$Ganancia)

# Definir una función para ajustar los tamaños de la fuente
adjust_font_size <- function(ranking) {
  return(7 - log(ranking))
}

# Reordenar los factores de la columna Solver según la Ganancia
datos$Solver <- factor(datos$Solver, levels = datos$Solver[order(datos$Ganancia, decreasing = TRUE)])

# Ajustar el tamaño de la fuente utilizando la función definida
datos$FontSize <- sapply(datos$Ranking, adjust_font_size)

# Generar el gráfico de barras
ggplot(datos, aes(x = Solver, y = Ganancia, fill = Color)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Ranking, "º"), size = FontSize), vjust = -0.3) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"), # números del eje y en negro
        plot.title = element_text(size = 20, hjust = 0.5, margin = margin(t = 20, b = 20), face = "bold"), # título del gráfico en negrita
        legend.position = "none",
        axis.title.x = element_text(face = "bold"), # título del eje x en negrita
        axis.title.y = element_text(face = "bold")) + # título del eje y en negrita
  ylab("Ganancia (en miles de pesos)") +
  xlab("Solver") +
  ggtitle("Ganancia máxima encontrada por solver") +
  coord_cartesian(ylim = c(100, 1250)) + # Ajusta los límites según lo necesites
  scale_fill_identity()

#######################

# librería
library(tidyverse)

# Datos
data <- data.frame(
  Solver = rep(c("constrOptim", "nloptr: COBYLA", "Rsolnp", "optimx: BFGS", 
                 "optimx: Nelder-Mead", "optimx: L-BFGS-B", "optimx: nlm", 
                 "optimx: Rcgmin", "optimx: Rvmmin", "optimx: nlminb", "nlminb"), each = 5),
  Iteraciones = rep(c(1, 100, 400, 700, 1000), times = 11),
  GananciaMax = c(142.6589, 1156.006, 1156.008, 1156.009, 1156.009,
                  71.64808, 71.69095, 71.69254, 71.69256, 71.69288,
                  75.28458, 1035.574, 1056.482, 1076.736, 1064.694,
                  111.8599, 1137.108, 1149.065, 1150.898, 1153.153,
                  142.4414, 1140.981, 1154.298, 1154.473, 1149.248,
                  142.4414, 1138.870, 1150.344, 1153.610, 1148.660,
                  66.37831, 1129.985, 1142.783, 1148.501, 1153.510,
                  66.38058, 1124.319, 1153.088, 1147.933, 1149.831,
                  66.38967, 1141.542, 1138.106, 1152.093, 1151.376,
                  66.40127, 1153.153, 1153.153, 924.9904, 1002.168,
                  66.37831, 1129.985, 1142.783, 1148.501, 1153.510) 
)

# Gráfico
ggplot(data, aes(x = Iteraciones, y = GananciaMax, group = Solver, color = Solver)) +
  geom_line() +
  geom_point(size = 4) +
  labs(x = "Número de Iteraciones", y = "Ganancia Máxima (en miles)", color = "Solver") +
  scale_y_continuous(limits = c(1040, 1158)) +
  theme_minimal()+
  ggtitle("Desempeño de Solver en función del Número de Iteraciones y Ganancia Máxima")


#
ggplot(data, aes(x = Iteraciones, y = GananciaMax, group = Solver, color = Solver)) +
  geom_line() +
  geom_point(size = 4) +
  facet_wrap(~ Solver, scales = "free") +  # Crear un gráfico para cada solver
  labs(x = "Número de Iteraciones", y = "Ganancia Máxima (en miles)", color = "Solver") +
  theme_minimal()


