"""Référence classique : LBM D1Q2 advection-diffusion, bords périodiques.

f[0] = f+ (vers la droite), f[1] = f- (vers la gauche).
Diffusivité D = tau - 1/2 (unités de réseau, c_s^2 = 1).
"""
import numpy as np


def equilibrium(phi, u):
    return np.array([0.5 * phi * (1 + u), 0.5 * phi * (1 - u)])


def collide(f, tau, u):
    phi = f.sum(axis=0)
    return f - (f - equilibrium(phi, u)) / tau


def stream(f):
    out = np.empty_like(f)
    out[0] = np.roll(f[0], +1)
    out[1] = np.roll(f[1], -1)
    return out


def step(f, tau, u):
    return stream(collide(f, tau, u))


def gaussian_init(N, center=None, width=None, u=0.0):
    x = np.arange(N)
    center = N / 4 if center is None else center
    width = N / 16 if width is None else width
    phi = np.exp(-((x - center) / width) ** 2)
    return equilibrium(phi, u)


def run(f0, tau, u, T):
    f = f0.copy()
    history = [f.copy()]
    for _ in range(T):
        f = step(f, tau, u)
        history.append(f.copy())
    return history
