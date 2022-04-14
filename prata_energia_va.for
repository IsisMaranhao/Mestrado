        Program fibonacci_prata
        Implicit none
        Integer ::  k
        Integer ::  n, np, j, f, nva
        integer , parameter :: dp = 8
        real (kind = dp) :: eta, pi, ea, eb, eg, va, vb, E_max, E_min 
        real (kind = dp) :: de, E, Gya, Gyb, Gyg, Dosa, Dosb
        real (kind = dp) :: Dosg, DosT, va_min, va_max, d_va
        Complex (kind = dp) :: t, Ec, eac, ebc, egc, vac, vbc, d
        Complex (kind = dp) :: vb2, va2, neac, Ga, Gb, Gg, f1, f2
        Complex (kind = dp) :: nebc, negc, nvac, nvbc
        real,dimension(100000000):: S
        open(unit=11,file='dos_prata.dat')
        open(unit=22,file='conj_cantor.dat')
        open(unit=44,file='Energia_VA.dat')
        n = 20
        t = (1+sqrt(complex(5.0_dp,0.0_dp)))/2
        eta = 5.d-4 !Alterar em relação ao np.
        pi = 4.d0*atan(1.d0)
        ea = 1.0_dp
        eb = -1.0_dp
        eg = -1.0_dp
        va_max = 4.0_dp
        va_min = -4.0_dp
        nva = 401
        d_va = (va_max -va_min)/(nva -1)
        do k = 1, nva
        va = va_min+(k-1)*d_va
        vb = 1.0_dp
        np = 2001 !Alterar número de pontos (Energia)
        E_max = 3.0_dp
        E_min = -7.0_dp
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
                va2 = vac*vac
                vb2 = vbc*vbc
                nvac = vac
                nvbc = vbc
                nebc = ebc
                neac = eac
                negc = egc
                d= ((Ec-neac)*(Ec-nebc))-va2
                vac = (va2*nvbc)/d
                vbc = nvac
                eac = negc+((va2*(Ec-nebc) + vb2*(Ec-neac))/d)
                ebc = negc+((vb2*(Ec-neac))/d)
                egc = neac+((va2*(Ec-nebc))/d)
            end do
            f1 = (1/t)
            f2 = ((t-1)/(2*t))
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
            If (DosT > 0.009) then
                S(j)= 1.0
                write(44,*) E, va
            end if
            write(11,*) E, DosT
            write(22,*) S(j), E
        end do
        end do
        close(11)
        close(22)
        End Program fibonacci_prata
