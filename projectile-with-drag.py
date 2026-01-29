# -*- coding: utf-8 -*-
"""
Created on Sat Nov  1 01:33:27 2025

@author: jamin
"""

import numpy as np
import matplotlib.pyplot as plt

# ---------------- Parámetros físicos ----------------
m = 7.26                  # masa del martillo [kg]
R = 0.06                  # radio del martillo [m]
rho = 1.2                 # densidad del aire [kg/m^3]
A = np.pi * R**2          # área transversal del martillo [m^2]
g = 9.81                  # gravedad [m/s^2]

theta = np.deg2rad(45.0)  # ángulo de lanzamiento (45°)
y0 = 2.0                  # altura inicial [m]

alcance_deseado = 86.74   # récord mundial (m)



def simul(v0, Cd, max_t=60.0, dt=0.001):

    vx = v0 * np.cos(theta)
    vy = v0 * np.sin(theta)
    x, y = 0.0, y0

    tray_x = [x]
    tray_y = [y]

    t = 0.0

    def f(x, y, vx, vy):

        v = np.hypot(vx, vy)
        if v == 0:
            return 0.0, -g
        Fd = 0.5 * rho * A * Cd * v**2
        ax = -(Fd/m)*(vx/v)
        ay = -g - (Fd/m)*(vy/v)
        return ax, ay

    while y > 0 and t < max_t:

        x_prev, y_prev, t_prev = x, y, t


        ax1, ay1 = f(x, y, vx, vy)
        vx1, vy1 = vx, vy

        ax2, ay2 = f(x + vx1*dt/2, y + vy1*dt/2, vx + ax1*dt/2, vy + ay1*dt/2)
        vx2, vy2 = vx + ax1*dt/2, vy + ay1*dt/2

        ax3, ay3 = f(x + vx2*dt/2, y + vy2*dt/2, vx + ax2*dt/2, vy + ay2*dt/2)
        vx3, vy3 = vx + ax2*dt/2, vy + ay2*dt/2

        ax4, ay4 = f(x + vx3*dt, y + vy3*dt, vx + ax3*dt, vy + ay3*dt)
        vx4, vy4 = vx + ax3*dt, vy + ay3*dt

        # actualizar x y y (RK4)
        x += dt*(vx1 + 2*vx2 + 2*vx3 + vx4)/6
        y += dt*(vy1 + 2*vy2 + 2*vy3 + vy4)/6

        # actualizar vx y vy (RK4)
        vx += dt*(ax1 + 2*ax2 + 2*ax3 + ax4)/6
        vy += dt*(ay1 + 2*ay2 + 2*ay3 + ay4)/6

        t += dt
        tray_x.append(x)
        tray_y.append(y)

        # Interpolación del impacto cuando cruza y = 0
        if y <= 0:
            alpha = y_prev / (y_prev - y)
            x_hit = x_prev + alpha*(x - x_prev)
            t_hit = t_prev + alpha*(t - t_prev)

            tray_x[-1] = x_hit
            tray_y[-1] = 0.0
            return tray_x, tray_y, t_hit

    return [0], [0], -1   # indica simulación inválida (v0 muy baja)



def biseccion(f,xmenos,xmas,Nmax,err): # x+,x−,Nmax,error
    for it in range(0,Nmax):
        x=(xmas+xmenos)/2.0   # Punto medio
        print('iteracion',it,'x=',x,'f(x)=',f(x))
        if (f(xmas)*f(x)>0.0):    # Raíz del otro lado
            xmas=x                  # Cambiamos de x+ a x
        else:
            xmenos=x                # Cambiamos de x− a x
        if abs(f(x))<err:        # Prueba de convergencia
            print('\n Raíz encontrada con precisión',err)
            break
        if it==Nmax-1:
            print('\n Raíz no encontrada después de Nmax iteracciones\n')
    return x

err=1e-3         # Precisión de la raíz
a=0.0;b=7.0     # Intervalo de búsqueda de raíz [a,b]
Nmax=100         # Número máximo de iteracciones

def v0(Cd, vmin=1.0, vmax=100.0):

    def f(v):
        x_list, _, t = simul(v, Cd)
        if t < 0:
            return -alcance_deseado
        return x_list[-1] - alcance_deseado

    a, b = vmin, vmax
    fa, fb = f(a), f(b)

    while fa * fb > 0:   # ampliar intervalo automáticamente
        b += 50
        fb = f(b)
        if b > 500:
            raise RuntimeError("No se encontró velocidad, revisa parámetros.")

    return biseccion(f, a, b,Nmax,err)


CDS = [0.0, 0.5, 0.75]
resultados = {}

for Cd in CDS:
    vreq = v0(Cd)
    x, y, tvuelo = simul(vreq, Cd)
    resultados[Cd] = (vreq, x, y, tvuelo)
    print(f"\nCd={Cd:.2f}  →  v0={vreq:.4f} m/s,  alcance={x[-1]:.3f} m, t={tvuelo:.3f} s")

plt.figure(figsize=(10,5))
for Cd, (vreq, x, y, tf) in resultados.items():
    plt.plot(x, y, label=f"Cd={Cd}, v0={vreq:.2f} m/s")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title("Trayectorias del martillo (45°)")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10,5))
for Cd, (vreq, x, y, tf) in resultados.items():
    tiempos = np.linspace(0, tf, len(y))
    plt.plot(tiempos, y, label=f"Cd={Cd}")
plt.xlabel("t [s]")
plt.ylabel("y [m]")
plt.title("Altura vs tiempo")
plt.legend()
plt.grid()
plt.show()

