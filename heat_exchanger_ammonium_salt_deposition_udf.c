#include "udf.h"
#include <math.h>

/*
 * Heat exchanger ammonium-salt deposition UDF.
 *
 * UDM(0): cumulative salt thickness on fluid cells next to exchanger wall [m].
 *
 * Hook suggestions in Fluent:
 * 1) DEFINE_SOURCE(ammonium_salt_deposition_sink): species source term (negative sink)
 * 2) DEFINE_ADJUST(update_ammonium_salt_fouling): update UDM(0) each time step
 * 3) DEFINE_PROFILE(heat_exchanger_wall_fouling_resistance): wall fouling resistance profile
 */

#define NH4_SALT_SPECIES_INDEX 0

#define GAS_CONSTANT 8.31446261815324      /* J/(mol*K) */
#define PRE_EXP_FACTOR 2.50e2              /* 1/s, tuned parameter */
#define ACTIVATION_ENERGY 3.80e4           /* J/mol */
#define SATURATION_TEMPERATURE 420.0       /* K */
#define DRIVING_TEMPERATURE_EPS 0.0        /* K */
#define SPECIES_ORDER 1.0                  /* dimensionless */
#define SALT_DENSITY 1700.0                /* kg/m^3 */
#define SALT_CONDUCTIVITY 0.35             /* W/(m*K) */
#define DEPOSITION_YIELD 1.0               /* kg salt / kg precursor */
#define MAX_SALT_THICKNESS 2.0e-3          /* m */

static real non_negative(real value)
{
    return (value > 0.0) ? value : 0.0;
}

static real deposition_rate_kg_m3_s(cell_t c, Thread *t)
{
    const real temp = C_T(c, t);
    const real y_precursor = non_negative(C_YI(c, t, NH4_SALT_SPECIES_INDEX));
    const real temp_driving = non_negative(SATURATION_TEMPERATURE - temp - DRIVING_TEMPERATURE_EPS);

    if (temp <= 1.0 || y_precursor <= 0.0 || temp_driving <= 0.0)
    {
        return 0.0;
    }

    return PRE_EXP_FACTOR
         * exp(-ACTIVATION_ENERGY / (GAS_CONSTANT * temp))
         * pow(y_precursor, SPECIES_ORDER)
         * temp_driving
         * C_R(c, t);
}

DEFINE_SOURCE(ammonium_salt_deposition_sink, c, t, dS, eqn)
{
    const real rate = deposition_rate_kg_m3_s(c, t);
    dS[eqn] = 0.0;

    /* Negative source -> precursor is consumed to form deposited ammonium salt. */
    return -DEPOSITION_YIELD * rate;
}

DEFINE_ADJUST(update_ammonium_salt_fouling, domain)
{
    Thread *t;
    cell_t c;

    const real dt = RP_Get_Real("physical-time-step");
    if (dt <= 0.0)
    {
        return;
    }

    thread_loop_c(t, domain)
    {
        if (!FLUID_THREAD_P(t))
        {
            continue;
        }

        begin_c_loop(c, t)
        {
            real thickness = C_UDMI(c, t, 0);
            const real growth_rate = deposition_rate_kg_m3_s(c, t) / SALT_DENSITY;

            thickness += growth_rate * dt;
            if (thickness > MAX_SALT_THICKNESS)
            {
                thickness = MAX_SALT_THICKNESS;
            }

            C_UDMI(c, t, 0) = non_negative(thickness);
        }
        end_c_loop(c, t)
    }
}

DEFINE_PROFILE(heat_exchanger_wall_fouling_resistance, t, i)
{
    face_t f;
    Thread *t0 = THREAD_T0(t);

    begin_f_loop(f, t)
    {
        const cell_t c0 = F_C0(f, t);
        const real thickness = non_negative(C_UDMI(c0, t0, 0));
        const real resistance = thickness / SALT_CONDUCTIVITY;

        F_PROFILE(f, t, i) = resistance;
    }
    end_f_loop(f, t)
}
