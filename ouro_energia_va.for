        Program fibonacci_ouro
        Implicit none
        Integer ::  a, b, p, c, ka, kb
        Integer ::  n, np, j, f, nva, nvb
        integer , parameter :: dp = 8
        real (kind = dp) :: eta, pi, ea, eb, eg, va, vb, E_max, E_min 
        real (kind = dp) :: de, E, Gya, Gyb, Gyg, Dosa, Dosb
        real (kind = dp) :: Dosg, DosT, va_min, va_max, vb_min, vb_max 
        real (kind = dp) :: d_va, d_vb
        Complex (kind = dp) :: t, Ec, eac, ebc, egc, vac, vbc, d, vavb
        Complex (kind = dp) :: vb2, va2, neac, Ga, Gb, Gg, f1, f2
        
        open(unit=44,file='Energia_VA.dat')
        
        a = 1
        b = 0
        n = 15
        c = 0
        
        !write(*,*) b
        !write(*,*) a
        
        do p = 1, n 
          c = a+b
          b = a
          a = c
        end do
       
        write(*,*) c
        
        t = (1+sqrt(complex(5.0_dp,0.0_dp)))/2
        eta = 0.005 !Alterar em relação ao np.
        pi = 4.d0*atan(1.d0)
        ea = 0.0_dp
        eb = 0.0_dp
        eg = 0.0_dp
        nva = 201
        nvb = 101
        
        vb_max = 5.0_dp
        vb_min = -5.0_dp
        d_vb = (vb_max -vb_min)/(nvb -1)
          
          do kb = 1, nvb, 50
            vb = vb_min+(kb-1)*d_vb
            write(*,*) "VB = ", vb
            
            va_max = 10.0_dp
            va_min = -10.0_dp
            d_va = (va_max -va_min)/(nva -1)
            
            do ka = 1, nva
              va = va_min+(ka-1)*d_va
              
              np = 201
              E_max = 10.0_dp
              E_min = -10.0_dp
              de = (E_max - E_min)/(np-1)
              
              do j = 1, np
                E = E_min+(j-1)*de
                
                Ec = dcmplx(E,eta)
                eac = dcmplx(ea,0.0_dp)
                ebc = dcmplx(eb,0.0_dp)
                egc = dcmplx(eg,0.0_dp)
                vac = dcmplx(va,0.0_dp)
                vbc = dcmplx(vb,0.0_dp)
                
                do f = 1, (n-1)
                  
                  d = Ec-ebc
                  vavb = vac*vbc
                  vb2 = vbc*vbc
                  va2 = vac*vac
                  neac = eac
                  vbc = vac
                  vac = vavb/d
                  ebc = egc+(vb2/d) 
                  eac = egc+((vb2+va2)/d)
                  egc = neac+(va2/d)
               
                end do
                
                f1 = (2*t-3)
                f2 = (2-t)
               
                Ga = f1/(Ec-eac)
                Gb = f2/(Ec-ebc)
                Gg = f2/(Ec-egc)
               
                Gya = imag(Ga)
                Gyb = imag(Gb)
                Gyg = imag(Gg)
               
                Dosa =(-1/pi)*Gya
                Dosb =(-1/pi)*Gyb
                Dosg =(-1/pi)*Gyg
                DosT = Dosa+Dosb+Dosg
                
                If (DosT > 0.5) then
                 
                  write(44,*) E, vb, va
               
                end if
              end do
            end do
          end do
        close(44)

        End Program fibonacci_ouro
