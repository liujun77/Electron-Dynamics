// *****************************************************************************
// *** integration.cpp                                                       ***
// *** Jun Liu                                                               ***
// *****************************************************************************





long numero integ(numero a, numero b, long numero f(numero x), int n)
{
 if(n<7)  n=7;
 int i;  numero  sum=0, h=(b-a)/n;
 if(a==b)  return 0;

 sum=(  17 * (f(a)     + f(b    ))  + 59 * (f(a+h)   + f(b-h  ))
      + 43 * (f(a+2*h) + f(b-2*h))  + 49 * (f(a+3*h) + f(b-3*h))) / 48;
 for(i=4; i<=n-4; i++)  sum+=f(a+i*h);

 return sum*h;
}
