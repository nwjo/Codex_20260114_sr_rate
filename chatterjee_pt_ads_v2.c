/*
 * Based on Chatterjee et al., 2001
 * Rewritten by NWJO, KATECH (Refactored for efficiency)
 * Date: 2025-10-29 (Refactored Version)
 * coverage-dependent k^(s) applied
*/

#include "udf.h"
#include <math.h>
#include <string.h>

/* ---------- Print Gates ---------- */
/* Using an array for cleaner management */
#define NUM_REACTIONS 60
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

/*** species index mapping (Hardcoded as requested) ***/
#define IDX_O2 0
#define IDX_C3H6 1
#define IDX_H2 2
#define IDX_H2O 3
#define IDX_CO2 4
#define IDX_CO 5
#define IDX_NO 6

#define IDX_Pt_Vac 9
#define IDX_O_Pt 10
#define IDX_C3H6_Pt 11
#define IDX_C3H5_Pt 18
#define IDX_OH_Pt 22
#define IDX_Rh_Vac 27

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
#define S0_O2       7.0E-2
#define S0_C3H6     9.8E-1
#define S0_C3H6_O   5.00E-2
#define S0_H2       4.60E-2
#define S0_H2O      7.50E-1
#define S0_CO2      5.00E-3
#define S0_CO       8.40E-1
#define S0_NO       8.50E-1
#define S1_O2       1.00E-2
#define S1_CO       5.00E-1
#define S1_NO       5.00E-1

// Reaction Orders (Site requirement)
#define q_R1    2.0
#define q_R2    2.0
#define q_R3    2.0
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
    real ln_k = log(MAX(EPS, A)) + beta * log(Tpos) - Ea_J_per_kmol * invRT;

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
    const real rho = C_R(c0, t0);
    return (rho * yi_k) / MAX(EPS, MW_k);
}

/* Thiele Modulus & Eta calculation */
static inline void thiele_eta_from_powerlaw(real rpp, real Cg, real nu_g, real *phi, real *eta)
{
    real kv_lin = 0.0;
    if (nu_g > 0.0 && Cg > EPS) kv_lin = A_V * rpp * (nu_g / Cg);
    const real ph = (kv_lin > 0.0) ? (DELTA_W * sqrt(kv_lin / D_EFF_FIXED)) : 0.0;
    const real et = (ph < 1.0e-6) ? (1.0 - ph * ph / 3.0) : tanh(ph) / ph;
    *phi = ph;
    *eta = et;
}

/* Sticking Coefficient to Pre-exponential Factor A */
static inline real A_from_sticking(real S0, real M_kg_per_kmol, real Gamma_kmol_m2, real q_site)
{
    const real root = sqrt(UNIVERSAL_GAS_CONSTANT / MAX(EPS, 2.0 * M_PI * M_kg_per_kmol));
    const real invG = pow(1.0 / MAX(EPS, Gamma_kmol_m2), q_site);
    return S0 * root * invG;
}

/* ======================================================================= */
/* Reaction Constants Definitions */
/* ======================================================================= */

/* Define reusable parameters for coverage dependency (Hardcoded arrays) */
/* Empty/Null Defaults */
static const int  idx_null[] = {0};
static const real val_null[] = {0.0};

/* R3: C3H6 + O(s) + Pt(s) -> ... */
static const int  idx_site_r3[] = { IDX_Pt_Vac };
static const real mu_r3[]       = { -0.9 };

/* R4: H2 adsorption */
static const int  idx_site_r4[] = { IDX_Pt_Vac };
static const real mu_r4[]       = { -1.0 };

/* R53: O2 on Rh */
static const int  idx_site_r53[] = { IDX_Rh_Vac };
static const real mu_r53[]       = { -1.0 };

/* A-Factors (Hardcoded from provided table, SI units) */
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
/* Helper: Common Calculation Logic for All Reactions                      */
/* ======================================================================= */
static inline void calc_apply_rate_common(
    /* Reaction ID/Name info */
    int r_id_num, const char *r_short_name,
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
        /* If this IS reaction 7, calculate its own parameters */
        /* Note: rate_base calculated above is effectively r7 logic,
           but we use the helper to ensure consistency for phi/eta calc */
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7, &rate_base, &phi, &eta);

        /* Special handling for R7 UDM updates */
        if (Cg > EPS) {
            C_UDMI(c0,t0, UDM_RCO_NET)    = rate_base;
            C_UDMI(c0,t0, UDM_RCO_NETP)   = 0;
            C_UDMI(c0,t0, UDM_DRDC)       = rate_base / Cg;
            C_UDMI(c0,t0, UDM_KAPP)       = A_V * rate_base / Cg;
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
        /* For all other reactions, use Reaction 7 as reference for eta */
        real r7_dummy;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7, &r7_dummy, &phi, &eta);
    }

    /* 5. Apply Washcoat Factor and Eta */
    *rr = rate_base * Wash_F * eta;

    /* 6. Logging (Thread-safe logic) */
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
DEFINE_SR_RATE(chatterjee_pt_ads_v2, f, fthread, r, mw, yi, rr)
{
    cell_t  c0  = F_C0(f, fthread);
    Thread *t0  = THREAD_T0(fthread);
    const real Tw  = F_T(f, fthread);

    /* Vacant Site Concentrations */
    const real theta_pt_vac = CLAMP01(yi[IDX_Pt_Vac]);
    const real theta_rh_vac = CLAMP01(yi[IDX_Rh_Vac]);

    /* Pre-calculate site terms.
       Note: pow(x, 1.0) is optimized to x. pow(x, 2.0) is x*x.
       Using helper variables for clarity. */
    const real term_pt = MAX(EPS, SITE_DEN_Pt * theta_pt_vac);
    const real term_rh = MAX(EPS, SITE_DEN_Rh * theta_rh_vac);

    /* Calculate Cv specifically for Reaction 7 (needed for reference eta) */
    const real Cv_R7_val = term_pt; // q_R7 is 1.0

    /* --- Reaction Dispatcher --- */

    if (STREQ(r->name, "reaction-1")) {
        // O2 + 2Pt -> 2O-Pt (q=2)
        calc_apply_rate_common(1, "r1", c0, t0, Tw, yi, rr,
            A1_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_O2, MW_O2, term_pt * term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-2")) {
        // C3H6 + 2Pt -> ... (q=2)
        calc_apply_rate_common(2, "r2", c0, t0, Tw, yi, rr,
            A2_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_C3H6, MW_C3H6, term_pt * term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-3")) {
        // C3H6 + O(s) + Pt -> ... (q=2 for Pt)
        calc_apply_rate_common(3, "r3", c0, t0, Tw, yi, rr,
            A3_k, 0.5, 0.0, idx_site_r3, mu_r3, val_null, 1,
            IDX_C3H6, MW_C3H6, term_pt * term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-4")) {
        // H2 + 2Pt -> ... (q=2)
        calc_apply_rate_common(4, "r4", c0, t0, Tw, yi, rr,
            A4_k, 0.5, 0.0, idx_site_r4, mu_r4, val_null, 1,
            IDX_H2, MW_H2, term_pt * term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-5")) {
        // H2O + Pt -> ... (q=1)
        calc_apply_rate_common(5, "r5", c0, t0, Tw, yi, rr,
            A5_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_H2O, MW_H2O, term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-6")) {
        // CO2 + Pt -> ... (q=1)
        calc_apply_rate_common(6, "r6", c0, t0, Tw, yi, rr,
            A6_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_CO2, MW_CO2, term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-7")) {
        // CO + Pt -> ... (q=1) - REFERENCE REACTION
        calc_apply_rate_common(7, "r7", c0, t0, Tw, yi, rr,
            A7_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_CO, MW_CO, term_pt, Cv_R7_val, 1); /* Flag=1 for Ref */
    }
    else if (STREQ(r->name, "reaction-48")) {
        // NO + Pt -> ... (q=1)
        calc_apply_rate_common(48, "r48", c0, t0, Tw, yi, rr,
            A48_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_NO, MW_NO, term_pt, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-53")) {
        // O2 + 2Rh -> ... (q=2)
        calc_apply_rate_common(53, "r53", c0, t0, Tw, yi, rr,
            A53_k, 0.5, 0.0, idx_site_r53, mu_r53, val_null, 1,
            IDX_O2, MW_O2, term_rh * term_rh, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-54")) {
        // CO + Rh -> ... (q=1)
        calc_apply_rate_common(54, "r54", c0, t0, Tw, yi, rr,
            A54_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_CO, MW_CO, term_rh, Cv_R7_val, 0);
    }
    else if (STREQ(r->name, "reaction-55")) {
        // NO + Rh -> ... (q=1)
        calc_apply_rate_common(55, "r55", c0, t0, Tw, yi, rr,
            A55_k, 0.5, 0.0, idx_null, val_null, val_null, 0,
            IDX_NO, MW_NO, term_rh, Cv_R7_val, 0);
    }
    else {
        return;
    }
}
