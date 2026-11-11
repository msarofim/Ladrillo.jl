# Two-basin Greenland component.
#
# Basins (Mouginot et al. 2019 sectors):
#   active = SW + CW + CE + SE + NW
#   high   = NO + NE
# Both basins read ONE regional driver T (`greenland_surface_temperature`, K rel. 1850-1900, the amplified regional
# series built in model.jl) and share the same A+B shape parameters. Each basin's commitment is scaled by its FIXED
# volume share k_b; the high basin's two channel rates are scaled by the sampled rate factor s_high (the active basin
# is the reference, s = 1).
#
# Per basin b and channel c in {fast (surface mass balance), slow (dynamic discharge)}, with T lagged one year:
#   eq_b  = min(max(k_b (c1 T + c0), 0), k_b v0)                    committed loss, clamped to the basin's own volume
#   r_c,b = min(max(s_b (alpha_c T + beta_c), 1e-9), 1)             relaxation rate per year
#   X_c,b <- X_c,b + (share_c eq_b - X_c,b) r_c,b,   share_fast = f, share_slow = 1 - f
#
# ABOVE-THRESHOLD DISCHARGE CHANNEL (the "tap"). A prior-specified reservoir on the whole sheet, opened by GLOBAL mean
# temperature (`gis_tap_gmt`, K rel. 1850-1900; NOT the amplified regional driver). A soft ramp
#   u = min(max((GMT - onset) / ramp_w, 0), 1)
# feeds a two-stage cascade with total mean delay tau (each stage relaxes at rate 2/tau), and V * stage-2 of extra
# loss is released, clamped to the whole sheet's remaining headroom. It is likelihood-inert over the calibration
# window (onset 4.69 K vs a calibration that tops out at 1.39 K) and is set from TAP_CELL, never fitted. V = 0 is OFF.
#
# OUTPUT CONTRACT: `greenland_sea_level` = gis_fast + gis_slow (m SLE) is the slot variable MimiBRICK reads; the tap
# rides in gis_slow (it is dynamic discharge, not surface mass balance).

"""The shipped above-threshold discharge cell (chosen 2026-08-23): onset in GLOBAL warming, released volume, total mean
delay, ramp width. Two cascade stages and the whole-sheet clamp are structural in `ladrillo_greenland`."""
const TAP_CELL = (onset_K = 4.69, V_m = 5.64, tau_yr = 800.0, ramp_w_K = 1.0)

"""Mouginot 2019 sector volumes (m SLE) and the FIXED two-basin volume shares derived from them. The NW sector (`mid`)
is merged into the active basin; the shares are fractions of the sector sum, applied to whatever v0 the model uses."""
const GIS_SECTOR_VOL_M = (south = 3.35, mid = 1.27, high = 2.73)
const GIS_K_ACTIVE = GIS_SECTOR_VOL_M.south / sum(GIS_SECTOR_VOL_M) + GIS_SECTOR_VOL_M.mid / sum(GIS_SECTOR_VOL_M)
const GIS_K_HIGH   = GIS_SECTOR_VOL_M.high / sum(GIS_SECTOR_VOL_M)

@defcomp ladrillo_greenland begin
    gis_c1      = Parameter()   # equilibrium sensitivity, m SLE per K of the regional driver
    gis_c0      = Parameter()   # committed loss at zero driver anomaly, m SLE
    gis_v0      = Parameter()   # total Greenland volume, m SLE (the commitment cap)
    gis_f       = Parameter()   # fast (surface-mass-balance) share of the commitment
    gis_alpha_f = Parameter()   # fast rate, per yr per K
    gis_beta_f  = Parameter()   # fast rate at zero anomaly, per yr
    gis_alpha_s = Parameter()   # slow (dynamic) rate, per yr per K
    gis_beta_s  = Parameter()   # slow rate at zero anomaly, per yr
    gis_g       = Parameter()   # fraction of the 1850 commitment already realised
    gis_k_active = Parameter()  # volume shares, FIXED (Mouginot geometry)
    gis_k_high   = Parameter()
    gis_s_high   = Parameter()  # high-basin rate scale (sampled, linear; the posterior carries log10)
    greenland_surface_temperature = Parameter(index=[time])   # REGIONAL driver, shared by both basins

    gis_tap_v      = Parameter()              # m SLE released by the tap; 0 = OFF
    gis_tap_onset  = Parameter()              # K, GLOBAL mean temperature rel. 1850-1900
    gis_tap_tau    = Parameter()              # yr, total mean delay of the two-stage cascade
    gis_tap_ramp_w = Parameter()              # K, width of the soft ramp
    gis_tap_gmt    = Parameter(index=[time])  # GLOBAL driver, K rel. 1850-1900

    gis_fast_active = Variable(index=[time])
    gis_fast_high   = Variable(index=[time])
    gis_slow_active = Variable(index=[time])
    gis_slow_high   = Variable(index=[time])
    gis_tap_s1      = Variable(index=[time])  # cascade stage 1, in [0, 1]
    gis_tap_s2      = Variable(index=[time])  # cascade stage 2, in [0, 1]
    gis_tap_wanted  = Variable(index=[time])  # V * stage 2, m
    gis_tap_applied = Variable(index=[time])  # after the whole-sheet headroom clamp, m
    gis_fast        = Variable(index=[time])
    gis_slow        = Variable(index=[time])  # includes the tap
    greenland_sea_level = Variable(index=[time])

    function run_timestep(p, v, d, t)
        if is_first(t)
            T = p.greenland_surface_temperature[t]
            eq_a = _gis_eq(T, p.gis_k_active, p.gis_c1, p.gis_c0, p.gis_v0)
            eq_h = _gis_eq(T, p.gis_k_high,   p.gis_c1, p.gis_c0, p.gis_v0)
            fa = p.gis_g * p.gis_f * eq_a
            fh = p.gis_g * p.gis_f * eq_h
            sa = p.gis_g * (1 - p.gis_f) * eq_a
            sh = p.gis_g * (1 - p.gis_f) * eq_h
            v.gis_tap_s1[t] = 0.0
            v.gis_tap_s2[t] = 0.0
        else
            Tm = p.greenland_surface_temperature[t-1]
            eq_a = _gis_eq(Tm, p.gis_k_active, p.gis_c1, p.gis_c0, p.gis_v0)
            eq_h = _gis_eq(Tm, p.gis_k_high,   p.gis_c1, p.gis_c0, p.gis_v0)
            fa = v.gis_fast_active[t-1] + (p.gis_f * eq_a - v.gis_fast_active[t-1]) *
                 _gis_rate(Tm, p.gis_alpha_f, p.gis_beta_f)
            fh = v.gis_fast_high[t-1] + (p.gis_f * eq_h - v.gis_fast_high[t-1]) *
                 _gis_rate(Tm, p.gis_s_high, p.gis_alpha_f, p.gis_beta_f)
            sa = v.gis_slow_active[t-1] + ((1 - p.gis_f) * eq_a - v.gis_slow_active[t-1]) *
                 _gis_rate(Tm, p.gis_alpha_s, p.gis_beta_s)
            sh = v.gis_slow_high[t-1] + ((1 - p.gis_f) * eq_h - v.gis_slow_high[t-1]) *
                 _gis_rate(Tm, p.gis_s_high, p.gis_alpha_s, p.gis_beta_s)
            u = min(max((p.gis_tap_gmt[t-1] - p.gis_tap_onset) / p.gis_tap_ramp_w, 0.0), 1.0)
            rate = 2.0 / p.gis_tap_tau
            v.gis_tap_s1[t] = v.gis_tap_s1[t-1] + (u - v.gis_tap_s1[t-1]) * rate
            v.gis_tap_s2[t] = v.gis_tap_s2[t-1] + (v.gis_tap_s1[t-1] - v.gis_tap_s2[t-1]) * rate
        end
        v.gis_fast_active[t], v.gis_fast_high[t] = fa, fh
        v.gis_slow_active[t], v.gis_slow_high[t] = sa, sh
        v.gis_tap_wanted[t] = p.gis_tap_v * v.gis_tap_s2[t]
        head = max(p.gis_v0 - (fa + fh + sa + sh), 0.0)
        v.gis_tap_applied[t] = min(v.gis_tap_wanted[t], head)
        v.gis_fast[t] = fa + fh
        v.gis_slow[t] = sa + sh + v.gis_tap_applied[t]
        v.greenland_sea_level[t] = v.gis_fast[t] + v.gis_slow[t]
    end
end

"""Basin committed loss: the whole-sheet law scaled by the volume share k, clamped to [0, k v0]."""
_gis_eq(T, k, c1, c0, v0) = min(max(k * (c1 * T + c0), 0.0), k * v0)
"""Relaxation rate per year, floored at 1e-9 and capped at 1 AFTER the rate scale is applied."""
_gis_rate(T, alpha, beta) = min(max(alpha * T + beta, 1e-9), 1.0)
_gis_rate(T, s, alpha, beta) = min(max(s * (alpha * T + beta), 1e-9), 1.0)
