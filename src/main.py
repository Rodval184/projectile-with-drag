# -*- coding: utf-8 -*-
"""
Hammer throw simulation with quadratic air drag.

Numerical integration (RK4) of projectile motion with drag.
The required initial velocity is obtained using the bisection method
to match a target horizontal range (world record).

Author: Benjamín Rodríguez Valdez
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------- Physical parameters ----------------
m = 7.26          # mass [kg]
R = 0.06          # radius [m]
rho = 1.2         # air density [kg/m^3]
A = np.pi * R**2  # cross-sectional area [m^2]
g = 9.81          # gravity [m/s^2]

theta = np.deg2rad(45.0)  # launch angle
y0 = 2.0                 # initial height [m]

TARGET_RANGE = 86.74     # world record [m]


def acceleration(vx, vy, Cd):
    v = np.hypot(vx, vy)
    if v == 0:
        return 0.0, -g
    Fd = 0.5 * rho * A * Cd * v**2
    ax = -(Fd / m) * (vx / v)
    ay = -g - (Fd / m) * (vy / v)
    return ax, ay


def simulate(v0, Cd, dt=1e-3, tmax=60.0):
    vx = v0 * np.cos(theta)
    vy = v0 * np.sin(theta)
    x, y = 0.0, y0

    xs, ys = [x], [y]
    t = 0.0

    while y > 0 and t < tmax:
        ax1, ay1 = acceleration(vx, vy, Cd)

        vx2 = vx + ax1 * dt / 2
        vy2 = vy + ay1 * dt / 2
        ax2, ay2 = acceleration(vx2, vy2, Cd)

        vx3 = vx + ax2 * dt / 2
        vy3 = vy + ay2 * dt / 2
        ax3, ay3 = acceleration(vx3, vy3, Cd)

        vx4 = vx + ax3 * dt
        vy4 = vy + ay3 * dt
        ax4, ay4 = acceleration(vx4, vy4, Cd)

        x += dt * (vx + 2*(vx + ax1*dt/2) + 2*(vx + ax2*dt/2) + (vx + ax3*dt)) / 6
        y += dt * (vy + 2*(vy + ay1*dt/2) + 2*(vy + ay2*dt/2) + (vy + ay3*dt)) / 6

        vx += dt * (ax1 + 2*ax2 + 2*ax3 + ax4) / 6
        vy += dt * (ay1 + 2*ay2 + 2*ay3 + ay4) / 6

        t += dt
        xs.append(x)
        ys.append(y)

    return np.array(xs), np.array(ys)


def find_v0(Cd, vmin=1.0, vmax=100.0, tol=1e-3):
    def f(v):
        x, y = simulate(v, Cd)
        return x[-1] - TARGET_RANGE

    a, b = vmin, vmax
    while f(a) * f(b) > 0:
        b += 50
        if b > 500:
            raise RuntimeError("Velocity not found.")

    for _ in range(100):
        c = 0.5 * (a + b)
        if abs(f(c)) < tol:
            return c
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
    return c


def main():
    Cd_values = [0.0, 0.5, 0.75]

    plt.figure(figsize=(10, 5))
    for Cd in Cd_values:
        v0 = find_v0(Cd)
        x, y = simulate(v0, Cd)
        plt.plot(x, y, label=f"Cd={Cd}, v0={v0:.2f} m/s")

    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.title("Hammer throw trajectories with air drag")
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()
