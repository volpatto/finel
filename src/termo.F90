        !> Module that contains thermodynamics function and parameters.
        !! @author Diego T. Volpatto
        module mtermo

            implicit none

            real*8 :: zcoef(10)

            contains

            !> Initializes interpolation coefficients to compute
            !! compressibility factor.
            !! @author Diego T. Volpatto
            subroutine init_zcoef

                implicit none

                zcoef(1) =  9.99921144125617722409060661448165774345397949218750d-01
                zcoef(2) = -1.19919829040115086206268166821135856547897446944262d-08
                zcoef(3) =  2.95290097079410864772724180594317997860333752145959d-16
                zcoef(4) =  9.42231835327529024372629441226386885009463920848079d-24
                zcoef(5) = -2.46929568577390678055712987160325876321312377285599d-31
                zcoef(6) = -7.40016953399249021667823040632784926032215181166411d-39
                zcoef(7) =  4.21756086831535086775256556143298098377555103824524d-46
                zcoef(8) = -7.90995787006734393072263251251053413138923617742665d-54
                zcoef(9) =  6.96584174374744927048653426883823314263335677413527d-62
                zcoef(10) = -2.42926517319393920651606665298459161434867679132808d-70

            endsubroutine

            !> Computes compressibility factor.
            !! @param p     [in] Pressure
            !! @param Z     [out] Compressibility factor
            !! @author Diego T. Volpatto
            subroutine Z_p(p, Z)

                implicit none

                real*8 :: p, Z

                 Z = zcoef(10)*(p**9)+zcoef(9)*(p**8)+zcoef(8)*(p**7) &
                 + zcoef(7)*(p**6)+ zcoef(6)*(p**5) + zcoef(5)*(p**4) &
                 + zcoef(4)*(p**3) + zcoef(3)*(p**2)+ zcoef(2)*p + zcoef(1)

            endsubroutine

        endmodule
