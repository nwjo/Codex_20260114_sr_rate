/*
 * Based on Chatterjee et al., 2001
 * Rewritten by NWJO, KATECH (Refactored for efficiency)
 * Date: 2025-10-29 (Refactored Version)
 * coverage-dependent k^(s) applied
 * Update: Added Desorption Reactions (R8-R59) as continuous else-if chain
*/

#include "udf.h"
#include <math.h>
#include <string.h>

/* ---------- Print Gates ---------- */
/* Using an array for cleaner management */
#define NUM_REACTIONS 65
static int print_gate[NUM_REACTIONS];

DEFINE_ADJUST(r1_reset_print_gate, d)
{
    /* Reset all print gates */
    int i;
    for (i = 0; i < NUM_REACTIONS; i++) print_gate[i] = 0;
}

/* ---------- UDM slots ---------- */
enum {
  UDM_RCO_NET   = 0,        /* r_CO,net (without base, eta)         [kmol/m^2/s] */
  UDM_RCO_NETP  = 1,        /* r_CO,net with Cco+(1+eps)            [kmol/m^2/s] */
  UDM_DRDC      = 2,        /* dr/dCco (coverage fixed)             [m/s]        */
  UDM_KAPP      = 3,        /* k_app,CO = a_w * dr/dCco             [1/s]        */
  UDM_PHI       = 4,        /* Thiele modulus phi                   [-]         */
  UDM_ETA       = 5         /* eta = tanh(phi)/phi                  [-]         */
};

/* guards */
#ifndef MAX
#define MAX(a,b) (( (a) > (b) ) ? (a) : (b))
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef CLAMP01
#define CLAMP01(x) ((x) < 0.0 ? 0.0 : ((x) > 1.0 ? 1.0 : (x)))
#endif

#ifndef STREQ
#define STREQ(a,b) (strcmp((a),(b))==0)
#endif

#define EPS 1.0E-30

/*** species index mapping (Combined) ***/
#define IDX_O2 0
#define IDX_C3H6 1
#define IDX_H2 2
#define IDX_H2O 3
#define IDX_CO2 4
#define IDX_CO 5
#define IDX_NO 6
//#define IDX_NO2 7
//#define IDX_N2 8

#define IDX_Pt_Vac 9
#define IDX_O_Pt 10
#define IDX_C3H6_Pt 11
#define IDX_H_Pt 12
#define IDX_H2O_Pt 13
#define IDX_CO2_Pt 14
#define IDX_CO_Pt 15
#define IDX_NO_Pt 16
#define IDX_N_Pt 17
#define IDX_C3H5_Pt 18
// #define IDX_C2H5_Pt 19
// #define IDX_CH2_Pt 20
// #define IDX_CH3_Pt 21
#define IDX_OH_Pt 22
// #define IDX_CH_Pt 23
// #define IDX_C_Pt 24
// #define IDX_CC2H5 25
// #define IDX_CH3CO 26
#define IDX_Rh_Vac 27
#define IDX_O_Rh 28
#define IDX_CO_Rh 29
#define IDX_NO_Rh 30
#define IDX_N_Rh 31

/*** Constants ***/
#define SITE_DEN_TOT 2.72E-8                /* kmol/m2 */
#define Pt_Frac 0.75
#define SITE_DEN_Pt (SITE_DEN_TOT*Pt_Frac)
#define SITE_DEN_Rh (SITE_DEN_TOT*(1.0-Pt_Frac))

// Molecular Weight [kg/kmol]
#define MW_O2   31.9988
#define MW_C3H6 42.081
#define MW_H2   2.01594
#define MW_H2O  18.01534
#define MW_CO2  44.00995
#define MW_CO   28.01055
#define MW_NO   30.0061

// Sticking Coefficients S0
#define S0_O2       7.0E-2      // R-1
#define S0_C3H6     9.8E-1      // R-2
#define S0_C3H6_O   5.00E-2     // R-3
#define S0_H2       4.60E-2     // R-4
#define S0_H2O      7.50E-1     // R-5
#define S0_CO2      5.00E-3     // R-6
#define S0_CO       8.40E-1     // R-7
#define S0_NO       8.50E-1     // R-48
#define S1_O2       1.00E-2     // R-53
#define S1_CO       5.00E-1     // R-54
#define S1_NO       5.00E-1     // R-55

// Reaction Orders (Site requirement)
#define q_R1    2.0
#define q_R2    2.0
#define q_R3    2.0 // fixed: 1.0 -> 2.0
#define q_R4    2.0
#define q_R5    1.0
#define q_R6    1.0
#define q_R7    1.0
#define q_R48   1.0
#define q_R53   2.0
#define q_R54   1.0
#define q_R55   1.0

// Washcoat Factor
#define Wash_F 70.0

/* Thiele Constants */
#define D_EFF_FIXED (3.40776e-7)
#define DELTA_W     (2.5e-5)
#define A_V         (2600)

#define NU_GAS_R7   1.0

/* ===== Helpers ===== */

/* Calculate k_surface with coverage dependency */
static inline real k_surface_covdep(real A, real beta, real Ea_J_per_kmol, real T,
                                    const int *idx_site, const real *mu, const real *eps_J_per_kmol,
                                    int Ns, const real yi[])
{
    const real Tpos  = MAX(EPS, T);
    const real invRT = 1.0 / (UNIVERSAL_GAS_CONSTANT * Tpos);
    real ln_k = log(MAX(EPS, A)) + beta*log(Tpos) - Ea_J_per_kmol*invRT;

    for (int i = 0; i < Ns; ++i) {
        const int  si   = idx_site ? idx_site[i] : 0;
        const real th   = MAX(1.0e-20, CLAMP01(yi[si]));
        const real mui  = mu  ? mu[i]  : 0.0;
        const real epsi = eps_J_per_kmol ? eps_J_per_kmol[i] : 0.0;
        ln_k += mui * log(th) + (epsi * th) * invRT;
    }
    return exp(ln_k);
}

/* Calculate gas concentration [kmol/m3] */
static inline real gas_conc_cell(cell_t c0, Thread *t0, real yi_k, real MW_k)
{
    const real T = MAX(EPS, C_T(c0, t0));
    const real P_abs = C_P(c0, t0) + RP_Get_Real("operating-pressure");
    const real C_tot = P_abs / MAX(EPS, UNIVERSAL_GAS_CONSTANT * T);
    (void)MW_k;
    return C_tot * yi_k;
}

/* Thiele Modulus & Eta calculation */
static inline void thiele_eta_from_powerlaw(real rpp, real Cg, real nu_g, real *phi, real *eta)
{
    real kv_lin = 0.0;
    if (nu_g > 0.0 && Cg > EPS) kv_lin = A_V * rpp * (nu_g / Cg);
    const real ph = (kv_lin > 0.0) ? (DELTA_W * sqrt(kv_lin / D_EFF_FIXED)) : 0.0;
    const real et = (ph < 1.0e-6) ? (1.0 - ph*ph/3.0) : tanh(ph)/ph;
    *phi = ph;
    *eta = et;
}

/* Sticking Coefficient to Pre-exponential Factor A */
static inline real A_from_sticking(real S0, real M_kg_per_kmol, real Gamma_kmol_m2, real q_site)    // ok?
{
    const real root = sqrt(UNIVERSAL_GAS_CONSTANT / MAX(EPS, 2.0 * M_PI * M_kg_per_kmol));
    const real invG = pow(1.0 / MAX(EPS, Gamma_kmol_m2), q_site);
    return S0 * root * invG;
}

/* ======================================================================= */
/* Reaction Constants Definitions */
/* ======================================================================= */

/* Define reusable parameters for coverage dependency */
static const int  idx_null[] = {0};
static const real val_null[] = {0.0};

/* Adsorption Parameters */
static const int  idx_site_r3[] = { IDX_Pt_Vac };
static const real mu_r3[]       = { -0.9 };
static const int  idx_site_r4[] = { IDX_Pt_Vac };
static const real mu_r4[]       = { -1.0 };
static const int  idx_site_r53[] = { IDX_Rh_Vac };
static const real mu_r53[]       = { -1.0 };

#define A1_k    6.084394E+14
#define A2_k    7.427986E+15
#define A3_k    3.789789E+14
#define A4_k    1.592985E+15
#define A5_k    2.363175E+08
#define A6_k    1.007976E+06
#define A7_k    2.122632E+08
#define A48_k   2.075130E+08
#define A53_k   8.691992E+13
#define A54_k   1.263471E+08
#define A55_k   1.220665E+08

// /* --- DESORPTION CONSTANTS (Added from rev19) --- */
// /* R8: 2O(s) -> O2 */
// static const int  idx_site_r8[] = {IDX_O_Pt};
// static const real eps_r8[]      = {9.0E7};
// #define NS_R8 1             // number of coverage dependent species
// #define A8_k 3.7E20
// #define B8_beta 0.0
// #define Ea8_Jpm 2.322E8

// /* R9: C3H6(s) -> C3H6 */
// #define NS_R9 0
// #define A9_k 1E13
// #define B9_beta 0.0
// #define Ea9_Jpm 7.27E7

// /* R10: C3H5(s) -> ... */
// #define NS_R10 0
// #define A10_k 3.7E20
// #define B10_beta 0.0
// #define Ea10_Jpm 3.1E7

// /* R11: 2H(s) -> H2 */
// static const int  idx_site_r11[] = {IDX_H_Pt};
// static const real eps_r11[]      = {6.0E6};         // fixed 6.0E7 -> 6.0E6
// #define NS_R11 1
// #define A11_k 3.7E20
// #define B11_beta 0.0
// #define Ea11_Jpm 6.74E7

// /* R12: H2O(s) -> H2O */
// #define NS_R12 0
// #define A12_k 1.0E13
// #define B12_beta 0.0
// #define Ea12_Jpm 4.03E7

// /* R13: CO(s) -> CO */
// static const int  idx_site_r13[] = {IDX_CO_Pt};
// static const real eps_r13[]      = {3.3E7};
// #define NS_R13 1
// #define A13_k 1.0E13
// #define B13_beta 0.0
// #define Ea13_Jpm 1.364E8

// /* R14: CO2(s) -> CO2 */
// #define NS_R14 0
// #define A14_k 1.0E13
// #define B14_beta 0.0
// #define Ea14_Jpm 2.71E7

// /* R49: NO(s) -> NO (Pt) */
// #define A49_k 1.0E16
// #define Ea49_Jpm 1.4E8

// /* R50: 2N(s) -> N2 (Pt) */
// static const int  idx_site_r50[] = {IDX_CO_Pt};
// static const real eps_r50[]      = {7.5E7};
// #define NS_R50 1
// #define A50_k 3.7E20
// #define B50_beta 0.0
// #define Ea50_Jpm 1.139E8

// /* R56: 2O(Rh) -> O2 */
// #define A56_k 3.0E20        /*miss, fixed*/
// #define B56_beta 0.0
// #define Ea56_Jpm 2.933E8

// /* R57: CO(Rh) -> CO */
// static const int  idx_site_r57[] = {IDX_CO_Rh, IDX_N_Rh};
// static const real eps_r57[]      = {1.88E7, 4.19E7};
// #define NS_R57 2
// #define A57_k 1.0E14
// #define B57_beta 0.0
// #define Ea57_Jpm 1.323E8

// /* R58: NO(Rh) -> NO */
// #define A58_k 5.0E13
// #define B58_beta 0.0
// #define Ea58_Jpm 1.089E8

// /* R59: 2N(Rh) -> N2 */
// static const int  idx_site_r59[] = {IDX_N_Rh};
// static const real eps_r59[]      = {1.67E7};
// #define NS_R59 1
// #define A59_k 1.11E18
// #define B59_beta 0.0
// #define Ea59_Jpm 1.369E8

/* ======================================================================= */
/* Helper: Reaction 7 Base Rate & Eta Calculation (The Reference)          */
/* ======================================================================= */
static inline void reaction7_base_and_eta(cell_t c0, Thread *t0,
                                          real Tw, const real *yi,
                                          real Cv_R7,
                                          real *rate7_base,
                                          real *phi7, real *eta7)
{
    /* Based on Reaction-7 (CO adsorption) */
    const real k7  = k_surface_covdep(A7_k, 0.5, 0.0, Tw, idx_null, val_null, val_null, 0, yi);
    const real Cg7 = gas_conc_cell(c0, t0, yi[IDX_CO], MW_CO);
    const real r7  = k7 * Cg7 * Cv_R7;      /* base rate r7'' */

    real phi, eta;
    thiele_eta_from_powerlaw(r7, Cg7, NU_GAS_R7, &phi, &eta);

    if (rate7_base) *rate7_base = r7;
    if (phi7)       *phi7       = phi;
    if (eta7)       *eta7       = eta;
}

/* ======================================================================= */
/* Helper: Common Calculation Logic for ADSORPTION Reactions               */
/* ======================================================================= */
static inline void calc_apply_rate_common(
    /* Reaction ID/Name info */
    int r_id_num, const char* r_short_name,
    /* Solver Pointers */
    cell_t c0, Thread *t0, real Tw, const real *yi, real *rr,
    /* Kinetic Params */
    real A, real beta, real Ea,
    const int *idx_site, const real *mu, const real *eps, int Ns,
    /* Species Info */
    int idx_gas, real MW_gas, real Cv_site,
    /* Reference Params (for Eta) */
    real Cv_R7,
    /* Logic Flags */
    int is_ref_reaction /* 1 if this IS reaction-7 */
)
{
    /* 1. Calculate Rate Constant k */
    const real k = k_surface_covdep(A, beta, Ea, Tw, idx_site, mu, eps, Ns, yi);

    /* 2. Calculate Gas Concentration */
    const real Cg = gas_conc_cell(c0, t0, yi[idx_gas], MW_gas);

    /* 3. Calculate Base Rate [kmol/m2-s] */
    real rate_base = k * Cg * Cv_site;

    /* 4. Determine Effectiveness Factor (eta) */
    real phi, eta;

    if (is_ref_reaction) {
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7, &rate_base, &phi, &eta);

        if (Cg > EPS) {
            C_UDMI(c0,t0, UDM_RCO_NET)    = rate_base;
            C_UDMI(c0,t0, UDM_RCO_NETP)   = 0;
            C_UDMI(c0,t0, UDM_DRDC)       = rate_base/Cg;
            C_UDMI(c0,t0, UDM_KAPP)       = A_V * rate_base/Cg;
            C_UDMI(c0,t0, UDM_PHI)        = phi;
            C_UDMI(c0,t0, UDM_ETA)        = eta;
        } else {
            C_UDMI(c0,t0, UDM_RCO_NET)    = 0;
            C_UDMI(c0,t0, UDM_RCO_NETP)   = 0;
            C_UDMI(c0,t0, UDM_DRDC)       = 0;
            C_UDMI(c0,t0, UDM_KAPP)       = 0;
            C_UDMI(c0,t0, UDM_PHI)        = 0;
            C_UDMI(c0,t0, UDM_ETA)        = 0;
        }

    } else {
        real r7_dummy;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7, &r7_dummy, &phi, &eta);
    }

    /* 5. Apply Washcoat Factor and Eta */
    *rr = rate_base * Wash_F * eta;

    /* 6. Logging */
    if (r_id_num < NUM_REACTIONS && print_gate[r_id_num] == 0) {
        #if RP_HOST
            Message("\n[%s] rate_base = %.6e | rate_wash = %.6e | rate_eta = %.6e (eta=%.3e)\n",
                    r_short_name, rate_base, rate_base * Wash_F, rate_base * Wash_F * eta, eta);
        #endif
        #if RP_NODE
            if (myid == 0) {
                 Message0("\n[%s] rate_base = %.6e | rate_wash = %.6e | rate_eta = %.6e (eta=%.3e)\n",
                        r_short_name, rate_base, rate_base * Wash_F, rate_base * Wash_F * eta, eta);
            }
        #endif
        print_gate[r_id_num] = 1;
    }
}


/* ======================================================================= */
/* Main UDF                                                                */
/* ======================================================================= */
DEFINE_SR_RATE(chatterjee_pt_ads_des_flat, f, fthread, r, mw, yi, rr)
{
    cell_t  c0  = F_C0(f, fthread);
    Thread *t0  = THREAD_T0(fthread);
    const real Tw  = F_T(f, fthread);

    /* Vacant Site Concentrations */
    const real theta_pt_vac = CLAMP01(yi[IDX_Pt_Vac]);
    const real theta_rh_vac = CLAMP01(yi[IDX_Rh_Vac]);
    const real theta_O_Pt = CLAMP01(yi[IDX_O_Pt]);

    /* Pre-calculate site terms */
    const real term_pt = MAX(EPS, SITE_DEN_Pt * theta_pt_vac);
    const real term_rh = MAX(EPS, SITE_DEN_Rh * theta_rh_vac);
    const real term_O_Pt = MAX(EPS, SITE_DEN_Pt * theta_O_Pt);

    /* Cv for Reaction 7 (Ref) */
    const real Cv_R7_val = term_pt;

    /* --- REACTION DISPATCHER (Continuous else if chain) --- */

    if (STREQ(r->name, "reaction-1")) {
        calc_apply_rate_common(1, "r1", c0, t0, Tw, yi, rr,
            A1_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_O2, MW_O2, term_pt*term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-2")) {
        calc_apply_rate_common(2, "r2", c0, t0, Tw, yi, rr,
            A2_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_C3H6, MW_C3H6, term_pt*term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-3")) {
        calc_apply_rate_common(3, "r3", c0, t0, Tw, yi, rr,
            A3_k, 0.5, 0.0, idx_site_r3, mu_r3, val_null, 1,
            IDX_C3H6, MW_C3H6, term_pt*term_O_Pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-4")) {
        calc_apply_rate_common(4, "r4", c0, t0, Tw, yi, rr,
            A4_k, 0.5, 0.0, idx_site_r4, mu_r4, val_null, 1,
            IDX_H2, MW_H2, term_pt*term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-5")) {
        calc_apply_rate_common(5, "r5", c0, t0, Tw, yi, rr,
            A5_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_H2O, MW_H2O, term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-6")) {
        calc_apply_rate_common(6, "r6", c0, t0, Tw, yi, rr,
            A6_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_CO2, MW_CO2, term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-7")) {
        calc_apply_rate_common(7, "r7", c0, t0, Tw, yi, rr,
            A7_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_CO, MW_CO, term_pt, Cv_R7_val, 1);
    }
    else if (STREQ(r->name, "reaction-48")) {
        calc_apply_rate_common(48, "r48", c0, t0, Tw, yi, rr,
            A48_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_NO, MW_NO, term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-53")) {
        calc_apply_rate_common(53, "r53", c0, t0, Tw, yi, rr,
            A53_k, 0.5, 0.0, idx_site_r53, mu_r53, val_null, 1,
            IDX_O2, MW_O2, term_rh*term_rh, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-54")) {
        calc_apply_rate_common(54, "r54", c0, t0, Tw, yi, rr,
            A54_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_CO, MW_CO, term_rh, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-55")) {
        calc_apply_rate_common(55, "r55", c0, t0, Tw, yi, rr,
            A55_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_NO, MW_NO, term_rh, Cv_R7_val, 0);
    }
    /* --- DESORPTION REACTIONS (Added from rev19, Flat Structure) --- */
    // else if (STREQ(r->name, "reaction-8")) {
    //     /* R8: 2O(s) -> O2 (Pt) */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k8 = k_surface_covdep(A8_k, B8_beta, Ea8_Jpm, Tw, idx_site_r8, val_null, eps_r8, NS_R8, yi);
    //     const real theta_O = CLAMP01(yi[IDX_O_Pt]);
    //     const real rate_base = k8 * theta_O * theta_O * SITE_DEN_Pt * SITE_DEN_Pt;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[8]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r8] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r8] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[8]=1;
    //     }
    // }
    // else if (STREQ(r->name, "reaction-9")) {
    //     /* R9: C3H6(s) -> C3H6 (Pt) */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k9 = k_surface_covdep(A9_k, B9_beta, Ea9_Jpm, Tw, idx_null, val_null, val_null, NS_R9, yi);
    //     const real theta = CLAMP01(yi[IDX_C3H6_Pt]);
    //     const real rate_base = k9 * theta * SITE_DEN_Pt;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[9]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r9] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r9] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[9]=1;
    //     }
    // }
    // else if (STREQ(r->name, "reaction-10")) {
    //     /* R10 */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k10 = k_surface_covdep(A10_k, B10_beta, Ea10_Jpm, Tw, idx_null, val_null, val_null, NS_R10, yi);
    //     const real thC3H5 = CLAMP01(yi[IDX_C3H5_Pt]);
    //     const real thOH = CLAMP01(yi[IDX_OH_Pt]);
    //     const real rate_base = k10 * thC3H5 * SITE_DEN_Pt * thOH * SITE_DEN_Pt;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[10]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r10] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r10] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[10]=1;
    //     }
    // }
    // else if (STREQ(r->name, "reaction-11")) {
    //     /* R11: 2H(s) -> H2 (Pt) */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k11 = k_surface_covdep(A11_k, B11_beta, Ea11_Jpm, Tw, idx_site_r11, val_null, eps_r11, NS_R11, yi);
    //     const real thH = CLAMP01(yi[IDX_H_Pt]);
    //     const real rate_base = k11 * thH * thH * SITE_DEN_Pt * SITE_DEN_Pt;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[11]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r11] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r11] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[11]=1;
    //     }
    // }
    // else if (STREQ(r->name, "reaction-12")) {
    //     /* R12: H2O(s) -> H2O (Pt) */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k12 = k_surface_covdep(A12_k, B12_beta, Ea12_Jpm, Tw, idx_null, val_null, val_null, NS_R12, yi);
    //     const real thH2O = CLAMP01(yi[IDX_H2O_Pt]);
    //     const real rate_base = k12 * thH2O * SITE_DEN_Pt;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[12]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r12] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r12] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[12]=1;
    //     }
    // }
    // else if (STREQ(r->name, "reaction-13")) {
    //     /* R13: CO(s) -> CO (Pt) */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k13 = k_surface_covdep(A13_k, B13_beta, Ea13_Jpm, Tw, idx_site_r13, val_null, eps_r13, NS_R13, yi);
    //     const real thCO = CLAMP01(yi[IDX_CO_Pt]);
    //     const real rate_base = k13 * thCO * SITE_DEN_Pt;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[13]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r13] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r13] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[13]=1;
    //     }
    // }
    // else if (STREQ(r->name, "reaction-14")) {
    //     /* R14: CO2(s) -> CO2 (Pt) */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k14 = k_surface_covdep(A14_k, B14_beta, Ea14_Jpm, Tw, idx_null, val_null, val_null, NS_R14, yi);
    //     const real thCO2 = CLAMP01(yi[IDX_CO2_Pt]);
    //     const real rate_base = k14 * thCO2 * SITE_DEN_Pt;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[14]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r14] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r14] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[14]=1;
    //     }
    // }
    // else if (STREQ(r->name, "reaction-49")) {
    //     /* R49: NO(s) -> NO (Pt) */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k49 = k_surface_covdep(A49_k, 0.0, Ea49_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
    //     const real thNO = CLAMP01(yi[IDX_NO_Pt]);
    //     const real rate_base = k49 * thNO * SITE_DEN_Pt;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[49]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r49] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r49] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[49]=1;
    //     }
    // }
    // else if (STREQ(r->name, "reaction-50")) {
    //     /* R50: 2N(s) -> N2 (Pt) */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k50 = k_surface_covdep(A50_k, B50_beta, Ea50_Jpm, Tw, idx_site_r50, val_null, eps_r50, NS_R50, yi);
    //     const real thN = CLAMP01(yi[IDX_N_Pt]);
    //     const real rate_base = k50 * thN * thN * SITE_DEN_Pt * SITE_DEN_Pt;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[50]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r50] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r50] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[50]=1;
    //     }
    // }
    // else if (STREQ(r->name, "reaction-56")) {
    //     /* R56: 2O(Rh) -> O2 */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k56 = k_surface_covdep(A56_k, B56_beta, Ea56_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
    //     const real thO_rh = CLAMP01(yi[IDX_O_Rh]);
    //     const real rate_base = k56 * thO_rh * thO_rh * SITE_DEN_Rh * SITE_DEN_Rh;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[56]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r56] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r56] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[56]=1;
    //     }
    // }
    // else if (STREQ(r->name, "reaction-57")) {
    //     /* R57: CO(Rh) -> CO */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k57 = k_surface_covdep(A57_k, B57_beta, Ea57_Jpm, Tw, idx_site_r57, val_null, eps_r57, NS_R57, yi);
    //     const real thCO_rh = CLAMP01(yi[IDX_CO_Rh]);
    //     const real rate_base = k57 * thCO_rh * SITE_DEN_Rh;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[57]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r57] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r57] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[57]=1;
    //     }
    // }
    // else if (STREQ(r->name, "reaction-58")) {
    //     /* R58: NO(Rh) -> NO */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k58 = k_surface_covdep(A58_k, B58_beta, Ea58_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
    //     const real thNO_rh = CLAMP01(yi[IDX_NO_Rh]);
    //     const real rate_base = k58 * thNO_rh * SITE_DEN_Rh;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[58]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r58] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r58] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[58]=1;
    //     }
    // }
    // else if (STREQ(r->name, "reaction-59")) {
    //     /* R59: 2N(Rh) -> N2 */
    //     real r7_dum, phi, eta;
    //     reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
    //     const real k59 = k_surface_covdep(A59_k, B59_beta, Ea59_Jpm, Tw, idx_site_r59, val_null, eps_r59, NS_R59, yi);
    //     const real thN_rh = CLAMP01(yi[IDX_N_Rh]);
    //     const real rate_base = k59 * thN_rh * thN_rh * SITE_DEN_Rh * SITE_DEN_Rh;
    //     *rr = rate_base * Wash_F * eta;

    //     if(print_gate[59]==0) {
    //         #if RP_NODE
    //         if(myid==0) Message0("\n[r59] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         #if RP_HOST
    //         Message("\n[r59] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
    //         #endif
    //         print_gate[59]=1;
    //     }
    // }
}
