import operator
import functools
import mpmath

ELECTRON_MASS = 0.511  # MeV/c^2

# https://stackoverflow.com/a/7948307
def prod(iterable):
    return functools.reduce(operator.mul, iterable, 1)


def product_particle_retention_terms(
        distance,
        electron_charge,
        positron_charge,
        vacuum_permativity,
        electron_number_density,
        ionisation_potential,
        initial_energy,
):
    for n in range(0, distance + 1):
        a = (electron_number_density * (positron_charge**2) * (electron_charge**4)) / \
            (8 * mpmath.pi * (vacuum_permativity**2))
        b = 4 / ionisation_potential

        constant_of_integration = mpmath.ei(2 * mpmath.ln(b * initial_energy)) / (a * (b**2))
        exp_integral_inv = exponential_integral_inverse(a * (b**2) * (constant_of_integration - n))
        denominator = mpmath.e ** (exp_integral_inv / 2)
        energy_proportionality_correction = mpmath.sqrt(ELECTRON_MASS / (2 * (electron_number_density**2)))

        yield (1 - (energy_proportionality_correction * mpmath.sqrt(b / denominator)))


def exponential_integral_inverse(x):
    raise RuntimeError("todo")  # TODO fill in

def particles_after_annihilation(
        initial_particles,
        distance,
        electron_charge,
        positron_charge,
        vacuum_permativity,
        electron_number_density,
        ionisation_potential,
):
    terms = product_particle_retention_terms(
        distance,
        electron_charge,
        positron_charge,
        vacuum_permativity,
        electron_number_density,
        ionisation_potential,
    )

    return initial_particles * prod(terms)
