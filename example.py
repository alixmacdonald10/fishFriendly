D  = outer diameter of leading edge

d = hub diameter of leading edge

Ri = d/2

Ro = D/2

Nblade = number of blades

 

N = shaft speed (rpm)

Omega = 2*pi*N/60

 

Lf = fish length

Bf = fish height

Compute Leff = Leff(Lf,Bf) (NEN8775 section 9.5)

 

%Initialise

n = number of intervals in Q (flow rate)

m = number of intervals in H (head)

nr = 100

Q = 0 : Qmax/n : Qmax

H = 0 : Hmax/m : Hmax

Pth(i,j) = 0 ; i=1 : n ; j=1 : m  

fMR(i,j) = 0 ; i=1 : n ; j=1 : m  

Pm(i,j) = 0 ; i=1 : n ; j=1 : m  

 

for i = 1:n

     for j  1:m

          Q = Q(i)

          H = H(j)

          A = 0

          dr = (Ro – Ri)/nr

          Pth = Pth(i,j)

          fMR = fMR(i,j)

          Pm = Pm(i,j)

          for ir = 1:nr

               r = Ri + (ir-1)*dr + dr/2

               dA = 2*pi*r*dr

               A = A + dA

               vm = vm(r)  (for a uniform flow:  vm = Q / (pi*(Ro^2-Ri^2) )

               dPth = dPth(Leff,vm,Omega,r,Nblade)  (NEN8775 eq. 4)

               dPth = max(0,min(1,dPth))

               beta = beta(r) (NEN8775 figure 15)

               delta = delta(r) (NEN8775 figure 16)

               d = d(r) (NEN8775 figure 17)

               vstrike = vstrike(vm, Omega,r,beta,delta) (NEN8775 eq. 17)

               dfMR = dfMR(Lf,d,vstrike)  (NEN8775 eq. 14 or 15)

               dfMR = max(0,min(1,dfMR))

               dPm = dPth * dfMR

               Pth = Pth + dPth*vm*dA

               fMR = fMR + dfMR* vm*dA

               Pm = Pm + dPm* vm*dA

          end

          Pth = Pth / Q

          fMR = fMR / Q

          Pm = Pm / Q

     end

end