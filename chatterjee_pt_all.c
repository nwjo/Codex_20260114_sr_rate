/*
 * Based on Chatterjee et al., 2001
 * Rewritten by NWJO, KATECH (Refactored for efficiency)
 * Date: 2025-10-29 (Refactored Version)
 * coverage-dependent k^(s) applied
 * Update: Added full surface reaction set (R1–R61) as continuous else-if chain
 * Notes:
 *  - gas_conc_cell uses species molecular weight for ideal-gas concentration.
 *  - surface intermediates indices expanded for hydrocarbon oxidation network.
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
#define IDX_C2H3_Pt 19
#define IDX_CH2_Pt 20
#define IDX_CH3_Pt 21
#define IDX_OH_Pt 22
#define IDX_CH_Pt 23
#define IDX_C_Pt 24
#define IDX_CC2H5_Pt 25
#define IDX_CH3CO_Pt 26
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
    const real R_specific = UNIVERSAL_GAS_CONSTANT / MW_k;
    const real rho = P_abs / MAX(EPS, R_specific * T);
    return (rho * yi_k) / MAX(EPS, MW_k);
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

/* --- DESORPTION CONSTANTS (Added from rev19) --- */
/* R8: 2O(s) -> O2 */
static const int  idx_site_r8[] = {IDX_O_Pt};
static const real eps_r8[]      = {9.0E7};
#define NS_R8 1             // number of coverage dependent species
#define A8_k 3.7E20
#define B8_beta 0.0
#define Ea8_Jpm 2.322E8

/* R9: C3H6(s) -> C3H6 */
#define NS_R9 0
#define A9_k 1E13
#define B9_beta 0.0
#define Ea9_Jpm 7.27E7

/* R10: C3H5(s) -> ... */
#define NS_R10 0
#define A10_k 3.7E20
#define B10_beta 0.0
#define Ea10_Jpm 3.1E7

/* R11: 2H(s) -> H2 */
static const int  idx_site_r11[] = {IDX_H_Pt};
static const real eps_r11[]      = {6.0E6};         // fixed 6.0E7 -> 6.0E6
#define NS_R11 1
#define A11_k 3.7E20
#define B11_beta 0.0
#define Ea11_Jpm 6.74E7

/* R12: H2O(s) -> H2O */
#define NS_R12 0
#define A12_k 1.0E13
#define B12_beta 0.0
#define Ea12_Jpm 4.03E7

/* R13: CO(s) -> CO */
static const int  idx_site_r13[] = {IDX_CO_Pt};
static const real eps_r13[]      = {3.3E7};
#define NS_R13 1
#define A13_k 1.0E13
#define B13_beta 0.0
#define Ea13_Jpm 1.364E8

/* R14: CO2(s) -> CO2 */
#define NS_R14 0
#define A14_k 1.0E13
#define B14_beta 0.0
#define Ea14_Jpm 2.71E7

/* R15: C3H5(s) + 5O(s) -> 5OH(s) + 3C(s) */
#define A15_k 3.7E16
#define Ea15_Jpm 9.5E7

/* R16: C3H6(s) + H(s) -> CC2H5(s) */
#define A16_k 1.0E12
#define Ea16_Jpm 7.54E7

/* R17: CC2H5(s) + H(s) -> C3H6(s) */
#define A17_k 3.7E20
#define Ea17_Jpm 4.88E7

/* R18: CC2H5(s) + Pt(s) -> C2H3(s) + CH2(s) */
#define A18_k 3.7E20
#define Ea18_Jpm 1.082E8

/* R19: C2H3(s) + CH2(s) + Pt(s) -> CC2H5(s) */
#define A19_k 3.7E20
#define Ea19_Jpm 3.2E6

/* R20: C2H3(s) + Pt(s) -> CH3(s) + C(s) */
#define A20_k 3.7E20
#define Ea20_Jpm 4.6E7

/* R21: CH3(s) + C(s) -> C2H3(s) + Pt(s) */
#define A21_k 3.7E20
#define Ea21_Jpm 4.69E7

/* R22: CH3(s) + Pt(s) -> CH2(s) + H(s) */
#define A22_k 1.26E21
#define Ea22_Jpm 7.04E7

/* R23: CH2(s) + H(s) -> CH3(s) + Pt(s) */
#define A23_k 3.09E21
#define Ea23_Jpm 0.0

/* R24: CH2(s) + Pt(s) -> CH(s) + H(s) */
#define A24_k 7.0E21
#define Ea24_Jpm 5.92E7

/* R25: CH(s) + H(s) -> CH2(s) + Pt(s) */
#define A25_k 3.09E21
#define Ea25_Jpm 0.0

/* R26: CH(s) + Pt(s) -> C(s) + H(s) */
#define A26_k 3.09E21
#define Ea26_Jpm 0.0

/* R27: C(s) + H(s) -> CH(s) + Pt(s) */
#define A27_k 1.25E21
#define Ea27_Jpm 1.38E8

/* R28: C2H3(s) + O(s) -> Pt(s) + CH3CO(s) */
#define A28_k 3.7E18
#define Ea28_Jpm 6.23E7

/* R29: CH3CO(s) + Pt(s) -> C2H3(s) + O(s) */
static const int  idx_site_r29[] = {IDX_O_Pt};
static const real eps_r29[]      = {-4.5E7};
#define NS_R29 1
#define A29_k 3.7E20
#define Ea29_Jpm 1.967E8

/* R30: CH3(s) + CO(s) -> Pt(s) + CH3CO(s) */
#define A30_k 3.7E20
#define Ea30_Jpm 8.29E7

/* R31: CH3CO(s) + Pt(s) -> CH3(s) + CO(s) */
#define A31_k 3.7E20
#define Ea31_Jpm 0.0

/* R32: CH3(s) + O(s) -> CH2(s) + OH(s) */
#define A32_k 3.7E20
#define Ea32_Jpm 3.66E7

/* R33: CH2(s) + OH(s) -> CH3(s) + O(s) */
#define A33_k 3.7E20
#define Ea33_Jpm 2.51E7

/* R34: CH2(s) + O(s) -> CH(s) + OH(s) */
#define A34_k 3.7E20
#define Ea34_Jpm 2.51E7

/* R35: CH(s) + OH(s) -> CH2(s) + O(s) */
#define A35_k 3.7E20
#define Ea35_Jpm 2.52E7

/* R36: CH(s) + O(s) -> C(s) + OH(s) */
#define A36_k 3.7E20
#define Ea36_Jpm 2.51E7

/* R37: C(s) + OH(s) -> CH(s) + O(s) */
#define A37_k 3.7E20
#define Ea37_Jpm 2.248E8

/* R38: O(s) + H(s) -> OH(s) + Pt(s) */
#define A38_k 3.7E20
#define Ea38_Jpm 1.15E7

/* R39: OH(s) + Pt(s) -> O(s) + H(s) */
#define A39_k 5.77E21
#define Ea39_Jpm 7.49E7

/* R40: H(s) + OH(s) -> H2O(s) + Pt(s) */
#define A40_k 3.7E20
#define Ea40_Jpm 1.74E7

/* R41: H2O(s) + Pt(s) -> H(s) + OH(s) */
#define A41_k 3.66E20
#define Ea41_Jpm 7.36E7

/* R42: OH(s) + OH(s) -> H2O(s) + O(s) */
#define A42_k 3.7E20
#define Ea42_Jpm 4.82E7

/* R43: H2O(s) + O(s) -> OH(s) + OH(s) */
#define A43_k 2.35E19
#define Ea43_Jpm 4.1E7

/* R44: CO(s) + O(s) -> CO2(s) + Pt(s) */
static const int  idx_site_r44[] = {IDX_CO_Pt, IDX_NO_Pt};
static const real eps_r44[]      = {3.3E7, -9.0E7};
#define NS_R44 2
#define A44_k 3.7E19
#define Ea44_Jpm 1.08E8

/* R45: CO2(s) + Pt(s) -> CO(s) + O(s) */
static const int  idx_site_r45[] = {IDX_O_Pt};
static const real eps_r45[]      = {-4.5E7};
#define NS_R45 1
#define A45_k 3.7E20
#define Ea45_Jpm 1.651E8

/* R46: C(s) + O(s) -> CO(s) + Pt(s) */
static const int  idx_site_r46[] = {IDX_CO_Pt};
static const real eps_r46[]      = {-3.3E7};
#define NS_R46 1
#define A46_k 3.7E20
#define Ea46_Jpm 0.0

/* R47: CO(s) + Pt(s) -> C(s) + O(s) */
static const int  idx_site_r47[] = {IDX_O_Pt};
static const real eps_r47[]      = {-4.5E7};
#define NS_R47 1
#define A47_k 3.7E20
#define Ea47_Jpm 2.185E8

/* R49: NO(s) -> NO (Pt) */
#define A49_k 1.0E16
#define Ea49_Jpm 1.4E8

/* R50: 2N(s) -> N2 (Pt) */
static const int  idx_site_r50[] = {IDX_CO_Pt};
static const real eps_r50[]      = {7.5E7};
#define NS_R50 1
#define A50_k 3.7E20
#define B50_beta 0.0
#define Ea50_Jpm 1.139E8

/* R51: NO(s) + Pt(s) -> N(s) + O(s) */
static const int  idx_site_r51[] = {IDX_CO_Pt};
static const real eps_r51[]      = {-3.0E6};
#define NS_R51 1
#define A51_k 5.0E19
#define Ea51_Jpm 1.078E8

/* R52: N(s) + O(s) -> NO(s) + Pt(s) */
static const int  idx_site_r52[] = {IDX_O_Pt};
static const real eps_r52[]      = {4.5E7};
#define NS_R52 1
#define A52_k 3.7E20
#define Ea52_Jpm 1.281E8

/* R56: 2O(Rh) -> O2 */
#define A56_k 3.0E20        /*miss, fixed*/
#define B56_beta 0.0
#define Ea56_Jpm 2.933E8

/* R57: CO(Rh) -> CO */
static const int  idx_site_r57[] = {IDX_CO_Rh, IDX_N_Rh};
static const real eps_r57[]      = {1.88E7, 4.19E7};
#define NS_R57 2
#define A57_k 1.0E14
#define B57_beta 0.0
#define Ea57_Jpm 1.323E8

/* R58: NO(Rh) -> NO */
#define A58_k 5.0E13
#define B58_beta 0.0
#define Ea58_Jpm 1.089E8

/* R59: 2N(Rh) -> N2 */
static const int  idx_site_r59[] = {IDX_N_Rh};
static const real eps_r59[]      = {1.67E7};
#define NS_R59 1
#define A59_k 1.11E18
#define B59_beta 0.0
#define Ea59_Jpm 1.369E8

/* R60: CO(Rh) + O(Rh) -> CO2 + Rh(s) + Rh(s) */
#define A60_k 3.7E19
#define Ea60_Jpm 5.99E7

/* R61: NO(Rh) + Rh(s) -> N(Rh) + O(Rh) */
#define A61_k 2.22E21
#define Ea61_Jpm 7.95E7
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
    else if (STREQ(r->name, "reaction-8")) {
        /* R8: 2O(s) -> O2 (Pt) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k8 = k_surface_covdep(A8_k, B8_beta, Ea8_Jpm, Tw, idx_site_r8, val_null, eps_r8, NS_R8, yi);
        const real theta_O = CLAMP01(yi[IDX_O_Pt]);
        const real rate_base = k8 * theta_O * theta_O * SITE_DEN_Pt * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[8]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r8] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r8] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[8]=1;
        }
    }
    else if (STREQ(r->name, "reaction-9")) {
        /* R9: C3H6(s) -> C3H6 (Pt) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k9 = k_surface_covdep(A9_k, B9_beta, Ea9_Jpm, Tw, idx_null, val_null, val_null, NS_R9, yi);
        const real theta = CLAMP01(yi[IDX_C3H6_Pt]);
        const real rate_base = k9 * theta * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[9]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r9] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r9] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[9]=1;
        }
    }
    else if (STREQ(r->name, "reaction-10")) {
        /* R10 */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k10 = k_surface_covdep(A10_k, B10_beta, Ea10_Jpm, Tw, idx_null, val_null, val_null, NS_R10, yi);
        const real thC3H5 = CLAMP01(yi[IDX_C3H5_Pt]);
        const real thOH = CLAMP01(yi[IDX_OH_Pt]);
        const real rate_base = k10 * thC3H5 * SITE_DEN_Pt * thOH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[10]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r10] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r10] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[10]=1;
        }
    }
    else if (STREQ(r->name, "reaction-11")) {
        /* R11: 2H(s) -> H2 (Pt) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k11 = k_surface_covdep(A11_k, B11_beta, Ea11_Jpm, Tw, idx_site_r11, val_null, eps_r11, NS_R11, yi);
        const real thH = CLAMP01(yi[IDX_H_Pt]);
        const real rate_base = k11 * thH * thH * SITE_DEN_Pt * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[11]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r11] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r11] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[11]=1;
        }
    }
    else if (STREQ(r->name, "reaction-12")) {
        /* R12: H2O(s) -> H2O (Pt) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k12 = k_surface_covdep(A12_k, B12_beta, Ea12_Jpm, Tw, idx_null, val_null, val_null, NS_R12, yi);
        const real thH2O = CLAMP01(yi[IDX_H2O_Pt]);
        const real rate_base = k12 * thH2O * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[12]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r12] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r12] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[12]=1;
        }
    }
    else if (STREQ(r->name, "reaction-13")) {
        /* R13: CO(s) -> CO (Pt) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k13 = k_surface_covdep(A13_k, B13_beta, Ea13_Jpm, Tw, idx_site_r13, val_null, eps_r13, NS_R13, yi);
        const real thCO = CLAMP01(yi[IDX_CO_Pt]);
        const real rate_base = k13 * thCO * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[13]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r13] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r13] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[13]=1;
        }
    }
    else if (STREQ(r->name, "reaction-14")) {
        /* R14: CO2(s) -> CO2 (Pt) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k14 = k_surface_covdep(A14_k, B14_beta, Ea14_Jpm, Tw, idx_null, val_null, val_null, NS_R14, yi);
        const real thCO2 = CLAMP01(yi[IDX_CO2_Pt]);
        const real rate_base = k14 * thCO2 * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[14]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r14] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r14] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[14]=1;
        }
    }

    else if (STREQ(r->name, "reaction-15")) {
        /* R15: C3H5(s) + 5O(s) -> 5OH(s) + 3C(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k15 = k_surface_covdep(A15_k, 0.0, Ea15_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thC3H5 = CLAMP01(yi[IDX_C3H5_Pt]);
        const real thO = CLAMP01(yi[IDX_O_Pt]);
        const real term_O = MAX(EPS, SITE_DEN_Pt * thO);
        const real rate_base = k15 * thC3H5 * SITE_DEN_Pt * pow(term_O, 5.0);
        *rr = rate_base * Wash_F * eta;

        if (print_gate[15] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r15] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r15] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[15] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-16")) {
        /* R16: C3H6(s) + H(s) -> CC2H5(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k16 = k_surface_covdep(A16_k, 0.0, Ea16_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thC3H6 = CLAMP01(yi[IDX_C3H6_Pt]);
        const real thH = CLAMP01(yi[IDX_H_Pt]);
        const real rate_base = k16 * thC3H6 * SITE_DEN_Pt * thH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[16] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r16] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r16] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[16] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-17")) {
        /* R17: CC2H5(s) + H(s) -> C3H6(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k17 = k_surface_covdep(A17_k, 0.0, Ea17_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCC2H5 = CLAMP01(yi[IDX_CC2H5_Pt]);
        const real thH = CLAMP01(yi[IDX_H_Pt]);
        const real rate_base = k17 * thCC2H5 * SITE_DEN_Pt * thH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[17] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r17] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r17] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[17] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-18")) {
        /* R18: CC2H5(s) + Pt(s) -> C2H3(s) + CH2(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k18 = k_surface_covdep(A18_k, 0.0, Ea18_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCC2H5 = CLAMP01(yi[IDX_CC2H5_Pt]);
        const real rate_base = k18 * thCC2H5 * SITE_DEN_Pt * term_pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[18] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r18] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r18] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[18] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-19")) {
        /* R19: C2H3(s) + CH2(s) + Pt(s) -> CC2H5(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k19 = k_surface_covdep(A19_k, 0.0, Ea19_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thC2H3 = CLAMP01(yi[IDX_C2H3_Pt]);
        const real thCH2 = CLAMP01(yi[IDX_CH2_Pt]);
        const real rate_base = k19 * thC2H3 * SITE_DEN_Pt * thCH2 * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[19] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r19] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r19] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[19] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-20")) {
        /* R20: C2H3(s) + Pt(s) -> CH3(s) + C(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k20 = k_surface_covdep(A20_k, 0.0, Ea20_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thC2H3 = CLAMP01(yi[IDX_C2H3_Pt]);
        const real rate_base = k20 * thC2H3 * SITE_DEN_Pt * term_pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[20] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r20] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r20] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[20] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-21")) {
        /* R21: CH3(s) + C(s) -> C2H3(s) + Pt(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k21 = k_surface_covdep(A21_k, 0.0, Ea21_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH3 = CLAMP01(yi[IDX_CH3_Pt]);
        const real thC = CLAMP01(yi[IDX_C_Pt]);
        const real rate_base = k21 * thCH3 * SITE_DEN_Pt * thC * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[21] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r21] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r21] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[21] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-22")) {
        /* R22: CH3(s) + Pt(s) -> CH2(s) + H(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k22 = k_surface_covdep(A22_k, 0.0, Ea22_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH3 = CLAMP01(yi[IDX_CH3_Pt]);
        const real rate_base = k22 * thCH3 * SITE_DEN_Pt * term_pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[22] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r22] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r22] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[22] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-23")) {
        /* R23: CH2(s) + H(s) -> CH3(s) + Pt(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k23 = k_surface_covdep(A23_k, 0.0, Ea23_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH2 = CLAMP01(yi[IDX_CH2_Pt]);
        const real thH = CLAMP01(yi[IDX_H_Pt]);
        const real rate_base = k23 * thCH2 * SITE_DEN_Pt * thH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[23] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r23] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r23] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[23] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-24")) {
        /* R24: CH2(s) + Pt(s) -> CH(s) + H(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k24 = k_surface_covdep(A24_k, 0.0, Ea24_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH2 = CLAMP01(yi[IDX_CH2_Pt]);
        const real rate_base = k24 * thCH2 * SITE_DEN_Pt * term_pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[24] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r24] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r24] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[24] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-25")) {
        /* R25: CH(s) + H(s) -> CH2(s) + Pt(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k25 = k_surface_covdep(A25_k, 0.0, Ea25_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH = CLAMP01(yi[IDX_CH_Pt]);
        const real thH = CLAMP01(yi[IDX_H_Pt]);
        const real rate_base = k25 * thCH * SITE_DEN_Pt * thH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[25] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r25] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r25] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[25] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-26")) {
        /* R26: CH(s) + Pt(s) -> C(s) + H(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k26 = k_surface_covdep(A26_k, 0.0, Ea26_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH = CLAMP01(yi[IDX_CH_Pt]);
        const real rate_base = k26 * thCH * SITE_DEN_Pt * term_pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[26] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r26] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r26] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[26] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-27")) {
        /* R27: C(s) + H(s) -> CH(s) + Pt(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k27 = k_surface_covdep(A27_k, 0.0, Ea27_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thC = CLAMP01(yi[IDX_C_Pt]);
        const real thH = CLAMP01(yi[IDX_H_Pt]);
        const real rate_base = k27 * thC * SITE_DEN_Pt * thH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[27] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r27] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r27] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[27] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-28")) {
        /* R28: C2H3(s) + O(s) -> Pt(s) + CH3CO(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k28 = k_surface_covdep(A28_k, 0.0, Ea28_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thC2H3 = CLAMP01(yi[IDX_C2H3_Pt]);
        const real thO = CLAMP01(yi[IDX_O_Pt]);
        const real rate_base = k28 * thC2H3 * SITE_DEN_Pt * thO * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[28] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r28] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r28] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[28] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-29")) {
        /* R29: CH3CO(s) + Pt(s) -> C2H3(s) + O(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k29 = k_surface_covdep(A29_k, 0.0, Ea29_Jpm, Tw, idx_site_r29, val_null, eps_r29, NS_R29, yi);
        const real thCH3CO = CLAMP01(yi[IDX_CH3CO_Pt]);
        const real rate_base = k29 * thCH3CO * SITE_DEN_Pt * term_pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[29] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r29] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r29] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[29] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-30")) {
        /* R30: CH3(s) + CO(s) -> Pt(s) + CH3CO(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k30 = k_surface_covdep(A30_k, 0.0, Ea30_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH3 = CLAMP01(yi[IDX_CH3_Pt]);
        const real thCO = CLAMP01(yi[IDX_CO_Pt]);
        const real rate_base = k30 * thCH3 * SITE_DEN_Pt * thCO * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[30] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r30] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r30] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[30] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-31")) {
        /* R31: CH3CO(s) + Pt(s) -> CH3(s) + CO(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k31 = k_surface_covdep(A31_k, 0.0, Ea31_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH3CO = CLAMP01(yi[IDX_CH3CO_Pt]);
        const real rate_base = k31 * thCH3CO * SITE_DEN_Pt * term_pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[31] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r31] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r31] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[31] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-32")) {
        /* R32: CH3(s) + O(s) -> CH2(s) + OH(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k32 = k_surface_covdep(A32_k, 0.0, Ea32_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH3 = CLAMP01(yi[IDX_CH3_Pt]);
        const real thO = CLAMP01(yi[IDX_O_Pt]);
        const real rate_base = k32 * thCH3 * SITE_DEN_Pt * thO * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[32] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r32] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r32] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[32] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-33")) {
        /* R33: CH2(s) + OH(s) -> CH3(s) + O(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k33 = k_surface_covdep(A33_k, 0.0, Ea33_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH2 = CLAMP01(yi[IDX_CH2_Pt]);
        const real thOH = CLAMP01(yi[IDX_OH_Pt]);
        const real rate_base = k33 * thCH2 * SITE_DEN_Pt * thOH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[33] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r33] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r33] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[33] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-34")) {
        /* R34: CH2(s) + O(s) -> CH(s) + OH(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k34 = k_surface_covdep(A34_k, 0.0, Ea34_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH2 = CLAMP01(yi[IDX_CH2_Pt]);
        const real thO = CLAMP01(yi[IDX_O_Pt]);
        const real rate_base = k34 * thCH2 * SITE_DEN_Pt * thO * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[34] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r34] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r34] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[34] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-35")) {
        /* R35: CH(s) + OH(s) -> CH2(s) + O(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k35 = k_surface_covdep(A35_k, 0.0, Ea35_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH = CLAMP01(yi[IDX_CH_Pt]);
        const real thOH = CLAMP01(yi[IDX_OH_Pt]);
        const real rate_base = k35 * thCH * SITE_DEN_Pt * thOH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[35] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r35] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r35] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[35] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-36")) {
        /* R36: CH(s) + O(s) -> C(s) + OH(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k36 = k_surface_covdep(A36_k, 0.0, Ea36_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCH = CLAMP01(yi[IDX_CH_Pt]);
        const real thO = CLAMP01(yi[IDX_O_Pt]);
        const real rate_base = k36 * thCH * SITE_DEN_Pt * thO * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[36] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r36] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r36] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[36] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-37")) {
        /* R37: C(s) + OH(s) -> CH(s) + O(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k37 = k_surface_covdep(A37_k, 0.0, Ea37_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thC = CLAMP01(yi[IDX_C_Pt]);
        const real thOH = CLAMP01(yi[IDX_OH_Pt]);
        const real rate_base = k37 * thC * SITE_DEN_Pt * thOH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[37] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r37] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r37] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[37] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-38")) {
        /* R38: O(s) + H(s) -> OH(s) + Pt(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k38 = k_surface_covdep(A38_k, 0.0, Ea38_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thO = CLAMP01(yi[IDX_O_Pt]);
        const real thH = CLAMP01(yi[IDX_H_Pt]);
        const real rate_base = k38 * thO * SITE_DEN_Pt * thH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[38] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r38] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r38] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[38] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-39")) {
        /* R39: OH(s) + Pt(s) -> O(s) + H(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k39 = k_surface_covdep(A39_k, 0.0, Ea39_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thOH = CLAMP01(yi[IDX_OH_Pt]);
        const real rate_base = k39 * thOH * SITE_DEN_Pt * term_pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[39] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r39] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r39] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[39] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-40")) {
        /* R40: H(s) + OH(s) -> H2O(s) + Pt(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k40 = k_surface_covdep(A40_k, 0.0, Ea40_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thH = CLAMP01(yi[IDX_H_Pt]);
        const real thOH = CLAMP01(yi[IDX_OH_Pt]);
        const real rate_base = k40 * thH * SITE_DEN_Pt * thOH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[40] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r40] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r40] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[40] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-41")) {
        /* R41: H2O(s) + Pt(s) -> H(s) + OH(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k41 = k_surface_covdep(A41_k, 0.0, Ea41_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thH2O = CLAMP01(yi[IDX_H2O_Pt]);
        const real rate_base = k41 * thH2O * SITE_DEN_Pt * term_pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[41] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r41] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r41] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[41] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-42")) {
        /* R42: OH(s) + OH(s) -> H2O(s) + O(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k42 = k_surface_covdep(A42_k, 0.0, Ea42_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thOH = CLAMP01(yi[IDX_OH_Pt]);
        const real rate_base = k42 * thOH * SITE_DEN_Pt * thOH * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[42] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r42] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r42] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[42] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-43")) {
        /* R43: H2O(s) + O(s) -> OH(s) + OH(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k43 = k_surface_covdep(A43_k, 0.0, Ea43_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thH2O = CLAMP01(yi[IDX_H2O_Pt]);
        const real thO = CLAMP01(yi[IDX_O_Pt]);
        const real rate_base = k43 * thH2O * SITE_DEN_Pt * thO * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[43] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r43] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r43] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[43] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-44")) {
        /* R44: CO(s) + O(s) -> CO2(s) + Pt(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k44 = k_surface_covdep(A44_k, 0.0, Ea44_Jpm, Tw, idx_site_r44, val_null, eps_r44, NS_R44, yi);
        const real thCO = CLAMP01(yi[IDX_CO_Pt]);
        const real thO = CLAMP01(yi[IDX_O_Pt]);
        const real rate_base = k44 * thCO * SITE_DEN_Pt * thO * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[44] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r44] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r44] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[44] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-45")) {
        /* R45: CO2(s) + Pt(s) -> CO(s) + O(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k45 = k_surface_covdep(A45_k, 0.0, Ea45_Jpm, Tw, idx_site_r45, val_null, eps_r45, NS_R45, yi);
        const real thCO2 = CLAMP01(yi[IDX_CO2_Pt]);
        const real rate_base = k45 * thCO2 * SITE_DEN_Pt * term_pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[45] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r45] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r45] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[45] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-46")) {
        /* R46: C(s) + O(s) -> CO(s) + Pt(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k46 = k_surface_covdep(A46_k, 0.0, Ea46_Jpm, Tw, idx_site_r46, val_null, eps_r46, NS_R46, yi);
        const real thC = CLAMP01(yi[IDX_C_Pt]);
        const real thO = CLAMP01(yi[IDX_O_Pt]);
        const real rate_base = k46 * thC * SITE_DEN_Pt * thO * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[46] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r46] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r46] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[46] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-47")) {
        /* R47: CO(s) + Pt(s) -> C(s) + O(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k47 = k_surface_covdep(A47_k, 0.0, Ea47_Jpm, Tw, idx_site_r47, val_null, eps_r47, NS_R47, yi);
        const real thCO = CLAMP01(yi[IDX_CO_Pt]);
        const real rate_base = k47 * thCO * SITE_DEN_Pt * term_pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[47] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r47] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r47] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[47] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-51")) {
        /* R51: NO(s) + Pt(s) -> N(s) + O(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k51 = k_surface_covdep(A51_k, 0.0, Ea51_Jpm, Tw, idx_site_r51, val_null, eps_r51, NS_R51, yi);
        const real thNO = CLAMP01(yi[IDX_NO_Pt]);
        const real rate_base = k51 * thNO * SITE_DEN_Pt * term_pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[51] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r51] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r51] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[51] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-52")) {
        /* R52: N(s) + O(s) -> NO(s) + Pt(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k52 = k_surface_covdep(A52_k, 0.0, Ea52_Jpm, Tw, idx_site_r52, val_null, eps_r52, NS_R52, yi);
        const real thN = CLAMP01(yi[IDX_N_Pt]);
        const real thO = CLAMP01(yi[IDX_O_Pt]);
        const real rate_base = k52 * thN * SITE_DEN_Pt * thO * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[52] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r52] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r52] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[52] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-60")) {
        /* R60: CO(Rh) + O(Rh) -> CO2 + Rh(s) + Rh(s) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k60 = k_surface_covdep(A60_k, 0.0, Ea60_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thCO_rh = CLAMP01(yi[IDX_CO_Rh]);
        const real thO_rh = CLAMP01(yi[IDX_O_Rh]);
        const real rate_base = k60 * thCO_rh * SITE_DEN_Rh * thO_rh * SITE_DEN_Rh;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[60] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r60] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r60] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[60] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-61")) {
        /* R61: NO(Rh) + Rh(s) -> N(Rh) + O(Rh) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k61 = k_surface_covdep(A61_k, 0.0, Ea61_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thNO_rh = CLAMP01(yi[IDX_NO_Rh]);
        const real rate_base = k61 * thNO_rh * SITE_DEN_Rh * term_rh;
        *rr = rate_base * Wash_F * eta;

        if (print_gate[61] == 0) {
            #if RP_NODE
            if (myid == 0) Message0("\n[r61] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r61] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[61] = 1;
        }
    }
    else if (STREQ(r->name, "reaction-49")) {
        /* R49: NO(s) -> NO (Pt) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k49 = k_surface_covdep(A49_k, 0.0, Ea49_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thNO = CLAMP01(yi[IDX_NO_Pt]);
        const real rate_base = k49 * thNO * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[49]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r49] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r49] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[49]=1;
        }
    }
    else if (STREQ(r->name, "reaction-50")) {
        /* R50: 2N(s) -> N2 (Pt) */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k50 = k_surface_covdep(A50_k, B50_beta, Ea50_Jpm, Tw, idx_site_r50, val_null, eps_r50, NS_R50, yi);
        const real thN = CLAMP01(yi[IDX_N_Pt]);
        const real rate_base = k50 * thN * thN * SITE_DEN_Pt * SITE_DEN_Pt;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[50]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r50] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r50] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[50]=1;
        }
    }
    else if (STREQ(r->name, "reaction-56")) {
        /* R56: 2O(Rh) -> O2 */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k56 = k_surface_covdep(A56_k, B56_beta, Ea56_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thO_rh = CLAMP01(yi[IDX_O_Rh]);
        const real rate_base = k56 * thO_rh * thO_rh * SITE_DEN_Rh * SITE_DEN_Rh;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[56]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r56] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r56] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[56]=1;
        }
    }
    else if (STREQ(r->name, "reaction-57")) {
        /* R57: CO(Rh) -> CO */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k57 = k_surface_covdep(A57_k, B57_beta, Ea57_Jpm, Tw, idx_site_r57, val_null, eps_r57, NS_R57, yi);
        const real thCO_rh = CLAMP01(yi[IDX_CO_Rh]);
        const real rate_base = k57 * thCO_rh * SITE_DEN_Rh;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[57]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r57] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r57] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[57]=1;
        }
    }
    else if (STREQ(r->name, "reaction-58")) {
        /* R58: NO(Rh) -> NO */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k58 = k_surface_covdep(A58_k, B58_beta, Ea58_Jpm, Tw, idx_null, val_null, val_null, 0, yi);
        const real thNO_rh = CLAMP01(yi[IDX_NO_Rh]);
        const real rate_base = k58 * thNO_rh * SITE_DEN_Rh;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[58]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r58] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r58] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[58]=1;
        }
    }
    else if (STREQ(r->name, "reaction-59")) {
        /* R59: 2N(Rh) -> N2 */
        real r7_dum, phi, eta;
        reaction7_base_and_eta(c0, t0, Tw, yi, Cv_R7_val, &r7_dum, &phi, &eta);
        const real k59 = k_surface_covdep(A59_k, B59_beta, Ea59_Jpm, Tw, idx_site_r59, val_null, eps_r59, NS_R59, yi);
        const real thN_rh = CLAMP01(yi[IDX_N_Rh]);
        const real rate_base = k59 * thN_rh * thN_rh * SITE_DEN_Rh * SITE_DEN_Rh;
        *rr = rate_base * Wash_F * eta;

        if(print_gate[59]==0) {
            #if RP_NODE
            if(myid==0) Message0("\n[r59] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            #if RP_HOST
            Message("\n[r59] rate_base=%.6e | rate_wash=%.6e | rate_eta=%.6e (eta=%.3e)\n", rate_base, rate_base*Wash_F, *rr, eta);
            #endif
            print_gate[59]=1;
        }
    }
}
