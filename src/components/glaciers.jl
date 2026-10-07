# Three-reservoir glacier component (Ladrillo v1.0).
#
# Reservoirs (RGI regions):
#   R19    Antarctic periphery                 {19}
#   SLOWG  slow-responding glacier regions     {03, 09, 07, 06}
#   FASTG  the remaining 13 regions
#
# Per reservoir b, driven by its own area-weighted surface temperature T_b (K rel. 1850-1900), lagged one year:
#   S_eq,b = max(a_b (1 - exp(-b_b (T_b - T_off,b))), 0)          equilibrium committed loss, floored at 0
#   T_eq   = T_off,b - log(1 - S_b / a_b) / b_b                    temperature at which the current stock is in equilibrium
#   d      = T_b - T_eq                                            SIGNED distance from equilibrium
#   mult   = min(kappa_b |d|^nu_b, 1)
#   S_b   <- S_b + mult (S_eq,b - S_b)
# Because mult <= 1 the step is a convex combination of S_b and S_eq,b, so with S_eq,b >= 0 the stock cannot regrow past
# its 1850 extent. Regrowth (d < 0) runs at the same rate law as melt.
#
# OUTPUT CONTRACT: `gsic_sea_level` (the sum of the three reservoirs, m SLE) is the slot variable MimiBRICK's
# global_sea_level and Antarctic sea-level feedback read. `gsic_hind` = SLOWG + FASTG is the hindcast scope (the
# observational glacier target excludes the Antarctic periphery before 2019).

@defcomp ladrillo_glaciers begin
    gic_a_R19       = Parameter()
    gic_b_R19       = Parameter()
    gic_T_off_R19   = Parameter()
    gic_kappa_R19   = Parameter()
    gic_nu_R19      = Parameter()
    gic_a_SLOWG     = Parameter()
    gic_b_SLOWG     = Parameter()
    gic_T_off_SLOWG = Parameter()
    gic_kappa_SLOWG = Parameter()
    gic_nu_SLOWG    = Parameter()
    gic_a_FASTG     = Parameter()
    gic_b_FASTG     = Parameter()
    gic_T_off_FASTG = Parameter()
    gic_kappa_FASTG = Parameter()
    gic_nu_FASTG    = Parameter()
    gic_sl0         = Parameter()          # initial cumulative loss per reservoir (0)
    glacier_surface_temperature_R19   = Parameter(index=[time])
    glacier_surface_temperature_SLOWG = Parameter(index=[time])
    glacier_surface_temperature_FASTG = Parameter(index=[time])

    gsic_r19       = Variable(index=[time])
    gsic_slowg     = Variable(index=[time])
    gsic_fastg     = Variable(index=[time])
    gsic_hind      = Variable(index=[time])    # SLOWG + FASTG
    gsic_sea_level = Variable(index=[time])    # all three reservoirs: the slot contract

    function run_timestep(p, v, d, t)
        if is_first(t)
            v.gsic_r19[t]   = p.gic_sl0
            v.gsic_slowg[t] = p.gic_sl0
            v.gsic_fastg[t] = p.gic_sl0
        else
            v.gsic_r19[t] = _glacier_step(v.gsic_r19[t-1], p.glacier_surface_temperature_R19[t-1],
                                          p.gic_a_R19, p.gic_b_R19, p.gic_T_off_R19, p.gic_kappa_R19, p.gic_nu_R19)
            v.gsic_slowg[t] = _glacier_step(v.gsic_slowg[t-1], p.glacier_surface_temperature_SLOWG[t-1],
                                            p.gic_a_SLOWG, p.gic_b_SLOWG, p.gic_T_off_SLOWG,
                                            p.gic_kappa_SLOWG, p.gic_nu_SLOWG)
            v.gsic_fastg[t] = _glacier_step(v.gsic_fastg[t-1], p.glacier_surface_temperature_FASTG[t-1],
                                            p.gic_a_FASTG, p.gic_b_FASTG, p.gic_T_off_FASTG,
                                            p.gic_kappa_FASTG, p.gic_nu_FASTG)
        end
        v.gsic_hind[t] = v.gsic_slowg[t] + v.gsic_fastg[t]
        v.gsic_sea_level[t] = v.gsic_hind[t] + v.gsic_r19[t]
    end
end

"""One reservoir, one annual step of the law in the file header."""
function _glacier_step(S, T, a, b, T_off, kappa, nu)
    S_eq = max(a * (1 - exp(-b * (T - T_off))), 0.0)
    frac_left = max(1 - S / a, 1e-12)
    T_eq = T_off - log(frac_left) / b
    d = T - T_eq
    mult = min(kappa * abs(d)^nu, 1.0)
    return S + mult * (S_eq - S)
end
