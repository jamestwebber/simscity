#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import numpy as np
import scipy.special as ssp

import simscity.latent


def projection(n_latent: int, sparsity: float, scale: int | float) -> np.ndarray:
    """Generates a linear weighting from a latent space to feature space,
    potentially with some coefficients set to zero. Returns the weighting.

    :param n_latent: dimensionality of the latent space
    :param sparsity: probability that a program is important for the drug. An array
                     of size ``(n_latent,)`` sets the rate per program
    :param scale: scaling factor for program weighting. If an array
                  of size ``(n_latent,)`` is given, different values are used
                  for each of the programs
    :return: array of shape (n_latent,) with weighting
    """

    z_weights = simscity.latent.gen_weighting(
        n_latent, 1, sparsity=sparsity, scale=scale
    ).flatten()

    return z_weights


def doses(scale: int | float, n_conditions: int) -> np.ndarray:
    """
    Generates an array of uniformly-spaced values with a bit of random noise
    added in.

    :param scale: scale of the expected drug response data
    :param n_conditions: number of conditions (doses) desired
    :return: array shape (n_conditions,) with thresholds
    """

    dose_thresholds = np.linspace(
        -3 * scale, 3 * scale, n_conditions
    ) + np.random.normal(size=n_conditions, scale=1.0 / (n_conditions**2))

    return dose_thresholds


def response(
    latent_exp: np.ndarray, z_weights: np.ndarray, doses: np.ndarray
) -> np.ndarray:
    """Given an array of samples from a latent space, the weighting for a drug,
    and the dose thresholds, this calculates the expected outcome for each
    sample according to a logistic function

    :param latent_exp: array of samples with shape (n_samples, n_latent)
    :param z_weights: (n_latent,) weights for projection into the drug space
    :param doses: (n_conditions,) array of dose thresholds
    :return: (n_samples, n_conditions) array of outcomes
    """

    return ssp.expit(np.dot(latent_exp, z_weights)[..., None] + doses)
