# Hammer Throw Simulation with Air Drag

Numerical simulation of hammer throw projectile motion including quadratic air resistance.

## 📌 Description
This project models the motion of a hammer throw using Newtonian mechanics with quadratic air drag.
The equations of motion are integrated numerically using the Runge–Kutta 4th order method (RK4).

The initial launch velocity is obtained via the bisection method in order to match a target horizontal range corresponding to the current world record.

## 🛠️ Tools
- Python
- NumPy
- Matplotlib

## 📊 Methodology
- Quadratic drag force proportional to velocity squared
- Numerical integration using RK4
- Root-finding via bisection method
- Comparison of trajectories for different drag coefficients

## ▶️ How to run
```bash
pip install -r requirements.txt
python src/main.py

