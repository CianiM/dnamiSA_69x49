

 real(wp) ::  d1_rhs_rho_dx_0_nxp4p0p0nyp4m1k,d1_rhs_rho_dx_0_nxp4p0m1nyp4m1k,d1_rhs_rho_dx_0_nxp4p0m2nyp4m1k &
            ,d1_rhs_rho_dx_0_nxp4p0nyp4m1k &
            ,d1_rhs_rho_dx_1_nxp4p0p0nyp4m1k,d1_rhs_rho_dx_1_nxp4p0m1nyp4m1k,d1_rhs_rho_dx_1_nxp4p0m2nyp4m1k &
            ,d1_rhs_rho_dx_1_nxp4p0nyp4m1k &
            ,d1_rhs_rho_dx_2_nxp4p0p0nyp4m1k,d1_rhs_rho_dx_2_nxp4p0m1nyp4m1k,d1_rhs_rho_dx_2_nxp4p0m2nyp4m1k &
            ,d1_rhs_rho_dx_2_nxp4p0nyp4m1k &
            ,d1_rhs_rho_dx_3_nxp4p0p0nyp4m1k,d1_rhs_rho_dx_3_nxp4p0m1nyp4m1k,d1_rhs_rho_dx_3_nxp4p0m2nyp4m1k &
            ,d1_rhs_rho_dx_3_nxp4p0nyp4m1k &
            ,d1_rhs_rho_dy_0_nxp4p0nyp4m1m1k,d1_rhs_rho_dy_0_nxp4p0nyp4m1p0k,d1_rhs_rho_dy_0_nxp4p0nyp4m1p1k &
            ,d1_rhs_rho_dy_0_nxp4p0nyp4m1k &
            ,d2_rhs_u_dxdx_6_0_nxp4p0p0nyp4m1k_nxp4p0p0p0nyp4m1k,d2_rhs_u_dxdx_6_0_nxp4p0p0nyp4m1k_nxp4p0p0m1nyp4m1k,d2_rhs_u_dxdx_6_0_nxp4p0p0nyp4m1k_nxp4p0p0m2nyp4m1k &
            ,d2_rhs_u_dxdx_6_0_nxp4p0p0nyp4m1k &
            ,d2_rhs_u_dxdx_6_0_nxp4p0m1nyp4m1k_nxp4p0m1m1nyp4m1k,d2_rhs_u_dxdx_6_0_nxp4p0m1nyp4m1k_nxp4p0m1p0nyp4m1k,d2_rhs_u_dxdx_6_0_nxp4p0m1nyp4m1k_nxp4p0m1p1nyp4m1k &
            ,d2_rhs_u_dxdx_6_0_nxp4p0m1nyp4m1k &
            ,d2_rhs_u_dxdx_6_0_nxp4p0m2nyp4m1k_nxp4p0m2m1nyp4m1k,d2_rhs_u_dxdx_6_0_nxp4p0m2nyp4m1k_nxp4p0m2p0nyp4m1k,d2_rhs_u_dxdx_6_0_nxp4p0m2nyp4m1k_nxp4p0m2p1nyp4m1k &
            ,d2_rhs_u_dxdx_6_0_nxp4p0m2nyp4m1k &
            ,d2_rhs_u_dxdy_6_0_nxp4p0p0nyp4m1k_nxp4p0p0nyp4m1m1k,d2_rhs_u_dxdy_6_0_nxp4p0p0nyp4m1k_nxp4p0p0nyp4m1p0k,d2_rhs_u_dxdy_6_0_nxp4p0p0nyp4m1k_nxp4p0p0nyp4m1p1k &
            ,d2_rhs_u_dxdy_6_0_nxp4p0p0nyp4m1k &
            ,d2_rhs_u_dxdy_6_0_nxp4p0m1nyp4m1k_nxp4p0m1nyp4m1m1k,d2_rhs_u_dxdy_6_0_nxp4p0m1nyp4m1k_nxp4p0m1nyp4m1p0k,d2_rhs_u_dxdy_6_0_nxp4p0m1nyp4m1k_nxp4p0m1nyp4m1p1k &
            ,d2_rhs_u_dxdy_6_0_nxp4p0m1nyp4m1k &
            ,d2_rhs_u_dxdy_6_0_nxp4p0m2nyp4m1k_nxp4p0m2nyp4m1m1k,d2_rhs_u_dxdy_6_0_nxp4p0m2nyp4m1k_nxp4p0m2nyp4m1p0k,d2_rhs_u_dxdy_6_0_nxp4p0m2nyp4m1k_nxp4p0m2nyp4m1p1k &
            ,d2_rhs_u_dxdy_6_0_nxp4p0m2nyp4m1k &
            ,d1_rhs_u_dx_6_nxp4p0p0nyp4m1k,d1_rhs_u_dx_6_nxp4p0m1nyp4m1k,d1_rhs_u_dx_6_nxp4p0m2nyp4m1k &
            ,d1_rhs_u_dx_6_nxp4p0nyp4m1k &
            ,d1_rhs_u_dy_0_nxp4p0nyp4m1m1k,d1_rhs_u_dy_0_nxp4p0nyp4m1p0k,d1_rhs_u_dy_0_nxp4p0nyp4m1p1k &
            ,d1_rhs_u_dy_0_nxp4p0nyp4m1k &
            ,d2_rhs_u_dydx_1_0_nxp4p0nyp4m1m1k_nxp4p0p0nyp4m1m1k,d2_rhs_u_dydx_1_0_nxp4p0nyp4m1m1k_nxp4p0m1nyp4m1m1k,d2_rhs_u_dydx_1_0_nxp4p0nyp4m1m1k_nxp4p0m2nyp4m1m1k &
            ,d2_rhs_u_dydx_1_0_nxp4p0nyp4m1m1k &
            ,d2_rhs_u_dydx_1_0_nxp4p0nyp4m1p1k_nxp4p0p0nyp4m1p1k,d2_rhs_u_dydx_1_0_nxp4p0nyp4m1p1k_nxp4p0m1nyp4m1p1k,d2_rhs_u_dydx_1_0_nxp4p0nyp4m1p1k_nxp4p0m2nyp4m1p1k &
            ,d2_rhs_u_dydx_1_0_nxp4p0nyp4m1p1k &
            ,d2_rhs_u_dydy_1_0_nxp4p0nyp4m1m1k_nxp4p0nyp4m1m1m1k,d2_rhs_u_dydy_1_0_nxp4p0nyp4m1m1k_nxp4p0nyp4m1m1p0k,d2_rhs_u_dydy_1_0_nxp4p0nyp4m1m1k_nxp4p0nyp4m1m1p1k &
            ,d2_rhs_u_dydy_1_0_nxp4p0nyp4m1m1k &
            ,d2_rhs_u_dydy_1_0_nxp4p0nyp4m1p1k_nxp4p0nyp4m1p1p0k,d2_rhs_u_dydy_1_0_nxp4p0nyp4m1p1k_nxp4p0nyp4m1p1m1k,d2_rhs_u_dydy_1_0_nxp4p0nyp4m1p1k_nxp4p0nyp4m1p1m2k &
            ,d2_rhs_u_dydy_1_0_nxp4p0nyp4m1p1k &
            ,d1_rhs_u_dy_1_nxp4p0nyp4m1m1k,d1_rhs_u_dy_1_nxp4p0nyp4m1p0k,d1_rhs_u_dy_1_nxp4p0nyp4m1p1k &
            ,d1_rhs_u_dy_1_nxp4p0nyp4m1k &
            ,d1_rhs_v_dx_4_nxp4p0p0nyp4m1k,d1_rhs_v_dx_4_nxp4p0m1nyp4m1k,d1_rhs_v_dx_4_nxp4p0m2nyp4m1k &
            ,d1_rhs_v_dx_4_nxp4p0nyp4m1k &
            ,d2_rhs_v_dxdx_5_0_nxp4p0p0nyp4m1k_nxp4p0p0p0nyp4m1k,d2_rhs_v_dxdx_5_0_nxp4p0p0nyp4m1k_nxp4p0p0m1nyp4m1k,d2_rhs_v_dxdx_5_0_nxp4p0p0nyp4m1k_nxp4p0p0m2nyp4m1k &
            ,d2_rhs_v_dxdx_5_0_nxp4p0p0nyp4m1k &
            ,d2_rhs_v_dxdx_5_0_nxp4p0m1nyp4m1k_nxp4p0m1m1nyp4m1k,d2_rhs_v_dxdx_5_0_nxp4p0m1nyp4m1k_nxp4p0m1p0nyp4m1k,d2_rhs_v_dxdx_5_0_nxp4p0m1nyp4m1k_nxp4p0m1p1nyp4m1k &
            ,d2_rhs_v_dxdx_5_0_nxp4p0m1nyp4m1k &
            ,d2_rhs_v_dxdx_5_0_nxp4p0m2nyp4m1k_nxp4p0m2m1nyp4m1k,d2_rhs_v_dxdx_5_0_nxp4p0m2nyp4m1k_nxp4p0m2p0nyp4m1k,d2_rhs_v_dxdx_5_0_nxp4p0m2nyp4m1k_nxp4p0m2p1nyp4m1k &
            ,d2_rhs_v_dxdx_5_0_nxp4p0m2nyp4m1k &
            ,d2_rhs_v_dxdy_5_0_nxp4p0p0nyp4m1k_nxp4p0p0nyp4m1m1k,d2_rhs_v_dxdy_5_0_nxp4p0p0nyp4m1k_nxp4p0p0nyp4m1p0k,d2_rhs_v_dxdy_5_0_nxp4p0p0nyp4m1k_nxp4p0p0nyp4m1p1k &
            ,d2_rhs_v_dxdy_5_0_nxp4p0p0nyp4m1k &
            ,d2_rhs_v_dxdy_5_0_nxp4p0m1nyp4m1k_nxp4p0m1nyp4m1m1k,d2_rhs_v_dxdy_5_0_nxp4p0m1nyp4m1k_nxp4p0m1nyp4m1p0k,d2_rhs_v_dxdy_5_0_nxp4p0m1nyp4m1k_nxp4p0m1nyp4m1p1k &
            ,d2_rhs_v_dxdy_5_0_nxp4p0m1nyp4m1k &
            ,d2_rhs_v_dxdy_5_0_nxp4p0m2nyp4m1k_nxp4p0m2nyp4m1m1k,d2_rhs_v_dxdy_5_0_nxp4p0m2nyp4m1k_nxp4p0m2nyp4m1p0k,d2_rhs_v_dxdy_5_0_nxp4p0m2nyp4m1k_nxp4p0m2nyp4m1p1k &
            ,d2_rhs_v_dxdy_5_0_nxp4p0m2nyp4m1k &
            ,d1_rhs_v_dx_5_nxp4p0p0nyp4m1k,d1_rhs_v_dx_5_nxp4p0m1nyp4m1k,d1_rhs_v_dx_5_nxp4p0m2nyp4m1k &
            ,d1_rhs_v_dx_5_nxp4p0nyp4m1k &
            ,d1_rhs_v_dy_0_nxp4p0nyp4m1m1k,d1_rhs_v_dy_0_nxp4p0nyp4m1p0k,d1_rhs_v_dy_0_nxp4p0nyp4m1p1k &
            ,d1_rhs_v_dy_0_nxp4p0nyp4m1k &
            ,d2_rhs_v_dydx_1_0_nxp4p0nyp4m1m1k_nxp4p0p0nyp4m1m1k,d2_rhs_v_dydx_1_0_nxp4p0nyp4m1m1k_nxp4p0m1nyp4m1m1k,d2_rhs_v_dydx_1_0_nxp4p0nyp4m1m1k_nxp4p0m2nyp4m1m1k &
            ,d2_rhs_v_dydx_1_0_nxp4p0nyp4m1m1k &
            ,d2_rhs_v_dydx_1_0_nxp4p0nyp4m1p1k_nxp4p0p0nyp4m1p1k,d2_rhs_v_dydx_1_0_nxp4p0nyp4m1p1k_nxp4p0m1nyp4m1p1k,d2_rhs_v_dydx_1_0_nxp4p0nyp4m1p1k_nxp4p0m2nyp4m1p1k &
            ,d2_rhs_v_dydx_1_0_nxp4p0nyp4m1p1k &
            ,d2_rhs_v_dydy_1_0_nxp4p0nyp4m1m1k_nxp4p0nyp4m1m1m1k,d2_rhs_v_dydy_1_0_nxp4p0nyp4m1m1k_nxp4p0nyp4m1m1p0k,d2_rhs_v_dydy_1_0_nxp4p0nyp4m1m1k_nxp4p0nyp4m1m1p1k &
            ,d2_rhs_v_dydy_1_0_nxp4p0nyp4m1m1k &
            ,d2_rhs_v_dydy_1_0_nxp4p0nyp4m1p1k_nxp4p0nyp4m1p1p0k,d2_rhs_v_dydy_1_0_nxp4p0nyp4m1p1k_nxp4p0nyp4m1p1m1k,d2_rhs_v_dydy_1_0_nxp4p0nyp4m1p1k_nxp4p0nyp4m1p1m2k &
            ,d2_rhs_v_dydy_1_0_nxp4p0nyp4m1p1k &
            ,d1_rhs_v_dy_1_nxp4p0nyp4m1m1k,d1_rhs_v_dy_1_nxp4p0nyp4m1p0k,d1_rhs_v_dy_1_nxp4p0nyp4m1p1k &
            ,d1_rhs_v_dy_1_nxp4p0nyp4m1k &
            ,d2_rhs_et_dxdx_9_0_nxp4p0p0nyp4m1k_nxp4p0p0p0nyp4m1k,d2_rhs_et_dxdx_9_0_nxp4p0p0nyp4m1k_nxp4p0p0m1nyp4m1k,d2_rhs_et_dxdx_9_0_nxp4p0p0nyp4m1k_nxp4p0p0m2nyp4m1k &
            ,d2_rhs_et_dxdx_9_0_nxp4p0p0nyp4m1k &
            ,d2_rhs_et_dxdx_9_0_nxp4p0m1nyp4m1k_nxp4p0m1m1nyp4m1k,d2_rhs_et_dxdx_9_0_nxp4p0m1nyp4m1k_nxp4p0m1p0nyp4m1k,d2_rhs_et_dxdx_9_0_nxp4p0m1nyp4m1k_nxp4p0m1p1nyp4m1k &
            ,d2_rhs_et_dxdx_9_0_nxp4p0m1nyp4m1k &
            ,d2_rhs_et_dxdx_9_0_nxp4p0m2nyp4m1k_nxp4p0m2m1nyp4m1k,d2_rhs_et_dxdx_9_0_nxp4p0m2nyp4m1k_nxp4p0m2p0nyp4m1k,d2_rhs_et_dxdx_9_0_nxp4p0m2nyp4m1k_nxp4p0m2p1nyp4m1k &
            ,d2_rhs_et_dxdx_9_0_nxp4p0m2nyp4m1k &
            ,d1_rhs_et_dx_9_nxp4p0p0nyp4m1k,d1_rhs_et_dx_9_nxp4p0m1nyp4m1k,d1_rhs_et_dx_9_nxp4p0m2nyp4m1k &
            ,d1_rhs_et_dx_9_nxp4p0nyp4m1k &
            ,d1_rhs_et_dy_0_nxp4p0nyp4m1m1k,d1_rhs_et_dy_0_nxp4p0nyp4m1p0k,d1_rhs_et_dy_0_nxp4p0nyp4m1p1k &
            ,d1_rhs_et_dy_0_nxp4p0nyp4m1k &
            ,d2_rhs_et_dydy_1_0_nxp4p0nyp4m1m1k_nxp4p0nyp4m1m1m1k,d2_rhs_et_dydy_1_0_nxp4p0nyp4m1m1k_nxp4p0nyp4m1m1p0k,d2_rhs_et_dydy_1_0_nxp4p0nyp4m1m1k_nxp4p0nyp4m1m1p1k &
            ,d2_rhs_et_dydy_1_0_nxp4p0nyp4m1m1k &
            ,d2_rhs_et_dydy_1_0_nxp4p0nyp4m1p1k_nxp4p0nyp4m1p1p0k,d2_rhs_et_dydy_1_0_nxp4p0nyp4m1p1k_nxp4p0nyp4m1p1m1k,d2_rhs_et_dydy_1_0_nxp4p0nyp4m1p1k_nxp4p0nyp4m1p1m2k &
            ,d2_rhs_et_dydy_1_0_nxp4p0nyp4m1p1k &
            ,d1_rhs_et_dy_1_nxp4p0nyp4m1m1k,d1_rhs_et_dy_1_nxp4p0nyp4m1p0k,d1_rhs_et_dy_1_nxp4p0nyp4m1p1k &
            ,d1_rhs_et_dy_1_nxp4p0nyp4m1k &
            ,d2_rhs_nut_dydy_0_0_nxp4p0nyp4m1m1k_nxp4p0nyp4m1m1m1k,d2_rhs_nut_dydy_0_0_nxp4p0nyp4m1m1k_nxp4p0nyp4m1m1p0k,d2_rhs_nut_dydy_0_0_nxp4p0nyp4m1m1k_nxp4p0nyp4m1m1p1k &
            ,d2_rhs_nut_dydy_0_0_nxp4p0nyp4m1m1k &
            ,d2_rhs_nut_dydy_0_0_nxp4p0nyp4m1p1k_nxp4p0nyp4m1p1p0k,d2_rhs_nut_dydy_0_0_nxp4p0nyp4m1p1k_nxp4p0nyp4m1p1m1k,d2_rhs_nut_dydy_0_0_nxp4p0nyp4m1p1k_nxp4p0nyp4m1p1m2k &
            ,d2_rhs_nut_dydy_0_0_nxp4p0nyp4m1p1k &
            ,d1_rhs_nut_dy_1_nxp4p0nyp4m1m1k,d1_rhs_nut_dy_1_nxp4p0nyp4m1p0k,d1_rhs_nut_dy_1_nxp4p0nyp4m1p1k &
            ,d1_rhs_nut_dy_1_nxp4p0nyp4m1k &
            ,d1_rhs_nut_dy_2_nxp4p0nyp4m1m1k,d1_rhs_nut_dy_2_nxp4p0nyp4m1p0k,d1_rhs_nut_dy_2_nxp4p0nyp4m1p1k &
            ,d1_rhs_nut_dy_2_nxp4p0nyp4m1k &
            ,d1_rhs_nut_dy_0_nxp4p0nyp4m1m1k,d1_rhs_nut_dy_0_nxp4p0nyp4m1p0k,d1_rhs_nut_dy_0_nxp4p0nyp4m1p1k &
            ,d1_rhs_nut_dy_0_nxp4p0nyp4m1k 