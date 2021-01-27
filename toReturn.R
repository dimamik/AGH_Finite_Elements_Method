library(cubature)

e <- function(n,i){
  return (
    
    Vectorize(
      function(x){
        tmp_n <- n/2
        if (1 - abs(tmp_n * (x - i / tmp_n)) >= 0){
          return (1 - abs(tmp_n * (x - i / tmp_n)))
        }
        else{
          return (0)
        }
        
        
      }
    )
    
  )
}


e_prim <- function(n,i){
  return (
    
    Vectorize(
      function(x){
        tmp_n <- n/2
        if ((i - 1) / tmp_n <= x  && x< (i / tmp_n)){
          return (tmp_n)
        }
        else if (i / tmp_n <= x && x< (i + 1) / tmp_n){
          return (-tmp_n)
        }
        else{
          return (0)
        }
        
      }
    )
    
  )
}


E <- function(x){
  if (x<=1 && x>=0){
    return (3)
  }  
  else if (x>1 && x<=2){
    return (5)
  }
  else {
    return (0)
  }
}

B <- function(i,j,n){
  h<- 2/n
  return (
    cubintegrate(
      function(x){
        return (e_prim(n,i)(x) * e_prim(n,j)(x) * (E(x)))
      },
      max(0,i*h-h),min(2,i*h+h),method="hcubature"
      
    )$integral - (3 * e(n,i)(0) * e(n,j)(0))     
  )
}

L <- function(i,n){
  return (-10 * E(0) * e(n,i)(0))
}


Build_B_M <- function(n){

B_M <- matrix(nrow=n-1,ncol = n-1)

for (i in 1:(n-1))
  for (j in 1:(n-1)){
    if (abs(i-j)<=1){
      B_M[i,j] <- B(i-1,j-1,n)
    }
    else{
      B_M[i,j] = 0
    }
    
  }
  return (B_M)

}


Build_L_M <- function(n){
  L_M <- matrix(nrow=n-1,ncol = 1)
  for (i in 1:(n-1))
    L_M[i] <- L(i-1,n)
  return (L_M)
}

draw_results <- function(points,n,u_v){
x_v = c()
y_v = c()

for (curr in 1:points){
  
  x <- curr*2/points
  y <- 0
  for (v in 1:length(u_v)){
    y = y + u_v[v]*e(n,v-1)(x)
  }
  
  x_v = append(x_v,x)
  y_v = append(y_v,y)

  
}
 plot(x_v,y_v)
}

solve_equation <- function(n,points){
  B_M <- Build_B_M(n)
  L_M <- Build_L_M(n)
  u_v <- solve(B_M,L_M)
  draw_results(points,n,u_v)
  print(n)
}
solve_equation(100,100)


# Model Building

build_model <- function(n,points){
  B_M <- Build_B_M(n)
  L_M <- Build_L_M(n)
  u_v <- solve(B_M,L_M)
  x_v = c()
  y_v = c()
  for (curr in 1:points){
    
    x <- curr*2/points
    y <- 0
    
    for (v in 1:length(u_v)){
      y = y + u_v[v]*e(n,v-1)(x)
    }
    
    x_v = append(x_v,x)
    y_v = append(y_v,y)
  
  }
  frame <- data.frame(x_v,y_v)
 
  return (frame)
}
frame <- build_model(150,100)
model <- lm(frame$y_v ~ frame$x_v)
summary(model)
plot(frame)
abline(model)
