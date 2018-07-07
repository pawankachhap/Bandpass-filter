#Chebyshev Bandpass filter response using
#Coupling Matrix method

global n fo BW RL Lar f
global i=sqrt(-1)

#Filter specifications
n = input("Filter order : ")
fo = input("Center frequency(GHz) : ")
BW = input("Bandwidth(MHz) : ")
RL = input("Return Loss(dB) : ")
Lar = -10*log10(1-10^(-RL/10))    #ripple level in passband

#parameters used in Coupling matrix method
go=1
B=log(coth(Lar/17.37))
y=sinh(B/(2*n))

a=zeros(1,n)
b=zeros(1,n)
for k=[1:n]
  a(k)= sin((((2*k)-1)*pi)/(2*n))
  b(k)= ((y.^2)+(sin((k*pi)/n))^2)
endfor
g=zeros(1,n+1)
g(1)=(2*a(1))/y
for k=[2:n]
  g(k)=(4*a(k-1)*a(k))/(b(k-1)*g(k-1))
endfor

if (mod(n,2)==0)
  g(n+1)=(coth(B/4))^2
else
  g(n+1)= 1
endif

#The termination impedance matrix R (n x n)
R=zeros(n,n)
R(1,1)= 1/(go*g(1))
R(n,n)= 1/(g(n)*g(n+1))

#Symmetric Coupling Matrix M (n x n)
M=zeros(n,n)
for j=[1:n-1]
  M(j,j+1)=1/sqrt(g(j)*g(j+1))
  M(j+1,j)=M(j,j+1)
endfor

#Defining Identity matrix I (nxn)
I=zeros(n,n)
for j=[1:n]
  I(j,j)=1
endfor

#Lambda value
l_lim=((fo*1000)-(2*BW))*0.001
u_lim=((fo*1000)+(2*BW))*0.001
ele=((u_lim-l_lim)*2000)+1
f=[l_lim:0.0005:u_lim]
lambda=zeros(1,ele)
for k=[1:ele]
  lambda(k)=((fo*1000/BW)*((f(k)./fo)-(fo./f(k))))
endfor

#S-parameters for the bandpass filter S(1,1) and S(2,1)
t11=zeros(1,ele)
t21=zeros(1,ele)
temp3=zeros(n,n)

for k=[1:ele]
  temp=zeros(n,n)
  temp= ((lambda(k).*I)- (i*R)+ M)
  temp3=inv(temp)
  t11(k)= 1+ (2*i*R(1,1))*(temp3(1,1))
  t21(k)= (-2*i*sqrt(R(1,1)*R(n,n)))*(temp3(n,1))
endfor

#conversion into dB
S11=zeros(1,ele)
S21=zeros(1,ele)
for k=[1:ele]
  S11(k)=20*log10(abs(t11(k)))
  S21(k)=20*log10(abs(t21(k)))
endfor

#Plots
disp("Couplin Matrix [M]")
disp(M)

figure(1)
subplot(2,1,1)
plot(f,S11)
xlabel ("Frequency GHz")
ylabel ("S11")
title ("Plot of S11 parameter")

subplot(2,1,2)
plot(f,S21)
xlabel ("Frequency GHz")
ylabel ("S21")
title ("Plot of S21 parameter")

figure(2)
plot(f,S11,f,S21)
xlabel ("Frequency GHz")
ylabel("S(1,1) & S(2,1) dB")
title("Chebyshev filter response")