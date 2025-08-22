# Equilibrium Rates
lequil = function(p1, p2, L, m, M, t, R)
{
  K1=(M*(1-t) - L*(R-p1))^2 + 4*(M*(1+t)+L*(R-p2))*L*(p1-p2)
  K2=(M*(1-t) + L*(R-p2))^2 - 4*(M*(1+t)+L*(R-p1))*L*(p1-p2)
   if (p1 >= R- (R-p2)*t*M/(M+(R-p2)*L) & (p1<=R))
             {
                l1 <- 0 
                l2 <- (R-p2)*L*m/(M+(R-p2)*L)
             } else if  (
               (p1 <=  R - (R-p2)*t*M/(M-(R-p2)*L) & (p2 >= R - R*M/(t*M+R*L)) &
                  (p2 <= R-(1-t)*M/L)) ||
               ((p1<=p2) & (p2>=R-(1-t)*M/L))
             )
             {
                l1=(R-p1)*L*m/(t*M+(R-p1)*L)
                l2=0
             } else if ((p1>p2) & (p1 < R - (R-p2)*t*M/(M+(R-p2)*L)))
             {
                l2 = m*(L*(R-p1)- M *(1-t) + sqrt(K1))/(2*(M*(1+t) + L*(R-p2) ))
                l1 = ((R-p1)*m*L - t*M*l2)/(t*M+L*(R-p1))
             } else if (
               ((p1<=p2) & (p2 <= R - R*M/(t*M+R*L))) ||
               ( (p1>=  R - (R-p2)*t*M/(M-(R-p2)*L)) & (p1<=p2) & (p2 >=R - R*M/(t*M+R*L))  &
                  (p2<=R-(1-t)*M/L)  )
               )
             {
                l1 = m*(L*(R-p2)+ M *(1-t) + sqrt(K2))/(2*(M*(1+t) + L*(R-p1) ))
                l2 = ((R-p2)*m*L - M*l1)/(M+L*(R-p2))
             } else if ( (p1==p2) & (p2<R-(1-t)*M/L) )
             {
                l1 = m*( L*(R-p1) + M*(1-t))/(L*(R-p1)+ M*(1+t))
                l2 = m*( L*(R-p2) - M*(1-t))/(L*(R-p2)+ M*(1+t))
             } else 
             {
               l1=-1
               l2=-1
             }
  le=c(l1,l2)
  return(le)
}
#
# Equilibrium Prices
# 
pequil = function(l1, l2, L, m, M, t, R)
{
  if (l1==0)
  {
    p1=R-t*M*l2/(L*m)
    p2=R-M*l2/((m-l2)*L)
  } else if (l2==0)
  {
    p1=R - t*M*l1/((m-l1)*L)
    p2=R - M*l1/(L*m)
  } else if (l1 < m - t*(m-l2))
  {
    p1=R-M*t*(l1+l2)/((m-l1)*L)
    p2= p1 - M*l2*(1/(m-l2) - t/(m-l1))/L
  } else if (l1 > m - t*(m-l2) )
  {
    p2=R-M*(l1+l2)/((m-l2)*L)
    p1=p2 + M*l1*(1/(m-l2) - t/(m-l1))/L
  } else
  {
    p1=M*(l1+l1)/((m-l2)*L)
    p2=p1
  }
  pe=c(p1,p2)
  return(pe)
}
#
# p20 : Price so that l2e(p1, p20)=l2
#
p20 = function(p1, l2, L, m, M, t, R, eps)
{
  # Use bisection
  # If for p2=0 l2e<l2 then it is infeasible, p2=-1
  le0=lequil(p1, 0, L, m, M, t, R)
  print(le0)
  l20=le0[2]
  print(l20)
  if (l20 < l2) x=-1 else
  {
    x0=0
    x1=R
    f0=l20
    f1=0
    x=(x0+x1)/2
    lex=lequil(p1, x, L, m, M, t, R)
    f2=lex[2]
    while(abs(f2-l2)>eps)
    {
      if (f2>l2) x0=x else x1=x
      x=(x0+x1)/2
      lex=lequil(p1, x, L, m, M, t, R)
      f2=lex[2]
    }
  }
  l10=lex[1]
  p20=x
  sol=c(p20,l10)
  return(sol)
}

#
# Alternative check 
# Solve p1e(x,l2)=p1 and l2e(p1,y)=l2
#
