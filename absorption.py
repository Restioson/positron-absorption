import operator
import functools
import mpmath
import csv
import halleys_method

ELECTRON_MASS = 0.511  # MeV/c^2


# https://stackoverflow.com/a/7948307
def prod(iterable):
    return functools.reduce(operator.mul, iterable, 1)


def product_term(
        n,
        electron_charge,
        positron_charge,
        vacuum_permativity,
        electron_number_density,
        ionisation_potential,
        initial_energy,
):
    a = (electron_number_density * (positron_charge ** 2) * (electron_charge ** 4)) / \
        (8 * mpmath.pi * (vacuum_permativity ** 2))
    b = 4 / ionisation_potential

    constant_of_integration = mpmath.ei(2 * mpmath.ln(b * initial_energy)) / (a * (b ** 2))
    exp_integral_inv = exponential_integral_inverse(a * (b ** 2) * (constant_of_integration - n))
    denominator = mpmath.e ** (exp_integral_inv / 2)
    energy_proportionality_correction = mpmath.sqrt(ELECTRON_MASS / (2 * (electron_number_density ** 2)))

    return 1 - (mpmath.sqrt(b / denominator)) / energy_proportionality_correction


# Calculate the terms from scratch
def product_particle_retention_terms(
        distance,
        electron_charge,
        positron_charge,
        vacuum_permativity,
        electron_number_density,
        ionisation_potential,
        initial_energy,
):
    print(f"Looping {distance + 1} times")
    for n in range(0, distance + 1):
        yield product_term(
            n,
            electron_charge,
            positron_charge,
            vacuum_permativity,
            electron_number_density,
            ionisation_potential,
            initial_energy
        )


def exponential_integral_inverse(p):
    initial = None

    if p > 3.5:
        initial = p * mpmath.ln(p)
    elif 0.75 < p <= 3.5:
        initial = p + 1
    elif -0.5 < p <= 0.75:
        initial = 1.45137 + (p * 0.37251)
    elif -43.8 < p <= -0.5:
        initial = 1 + mpmath.e**(p + mpmath.euler)
    else:
        return 1

    return mpmath.ln(halleys_method.solve_fx(initial, p))


# Returns the particles after annihilation from scratch
def particles_after_annihilation(
        initial_particles,
        dist,
        width,
        electron_charge,
        positron_charge,
        vacuum_permativity,
        electron_number_density,
        ionisation_potential,
        initial_energy,
):
    terms = product_particle_retention_terms(
        round(dist / mpmath.mpf(width)),
        mpmath.mpf(electron_charge),
        mpmath.mpf(positron_charge),
        mpmath.mpf(vacuum_permativity),
        mpmath.mpf(electron_number_density),
        mpmath.mpf(ionisation_potential),
        mpmath.mpf(initial_energy),
    )

    term_prod = prod(terms)
    print(f"Term prod: {term_prod}")

    return mpmath.mpf(initial_particles) * term_prod


def loop_inner(initial_particles, dist, prev_product):
    width = mpmath.mpf('1')
    electron_charge = mpmath.mpf('-1.602176634e-19')
    positron_charge = mpmath.mpf('+1.602176634e-19')
    vacuum_permativity = mpmath.mpf("8.85e-12")
    electron_number_density = mpmath.mpf("3.333e-23")
    ionisation_potential = 75
    initial_energy = 4e6

    term = product_term(
        round(dist / mpmath.mpf(width)),
        mpmath.mpf(electron_charge),
        mpmath.mpf(positron_charge),
        mpmath.mpf(vacuum_permativity),
        mpmath.mpf(electron_number_density),
        mpmath.mpf(ionisation_potential),
        mpmath.mpf(initial_energy),
    )

    product = term * prev_product
    return mpmath.mpf(initial_particles) * product, product


mpmath.mp.dps = 100
previous_product = 1

initial = 10**26

with open('out.csv', 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=';')

    for distance in range(0, 10**3 + 1):
        print(f"Particles after annihilation, distance = {distance}:")

        particles_left, previous_product = loop_inner(initial, distance, previous_product)

        mpmath.mp.dps = 26
        print(particles_left - initial)
        writer.writerow([distance, particles_left - initial])
        mpmath.mp.dps = 100
