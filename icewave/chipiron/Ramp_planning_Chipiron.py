#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np


# USER INPUT


ramp_type = input("Ramp type (linear/log): ").strip().lower()

omega_start = float(input("Starting omega [rad/s]: "))

omega_end = float(input("Ending omega [rad/s]: "))

n_points = int(input("Number of points: "))

time_per_point = float(input("Time per point [s]: "))



# CHECK INPUTS



if ramp_type not in ["linear", "log"]:
    raise ValueError("Ramp type must be 'linear' or 'log'.")

if n_points < 2:
    raise ValueError("Number of points must be at least 2.")

if time_per_point <= 0:
    raise ValueError("Time per point must be > 0.")

if ramp_type == "log":
    if omega_start <= 0 or omega_end <= 0:
        raise ValueError("Logarithmic ramps require omega > 0.")



# GENERATE OMEGA VALUES


if ramp_type == "linear":

    omega_values = np.linspace(
        omega_start,
        omega_end,
        n_points
    )

elif ramp_type == "log":

    omega_values = np.geomspace(
        omega_start,
        omega_end,
        n_points
    )



# CALCULATIONS


total_time_s = n_points * time_per_point
total_time_min = total_time_s / 60
total_time_h = total_time_s / 3600

delta_omega = np.diff(omega_values)

if ramp_type == "log":

    n_decades = abs(
        np.log10(omega_end) - np.log10(omega_start)
    )

    points_per_decade = (n_points - 1) / n_decades

    ratio = omega_values[1] / omega_values[0]

else:

    n_decades = None
    points_per_decade = None
    ratio = None


# OUTPUT


print("\n" + "=" * 60)
print("RAMP SUMMARY")
print("=" * 60)

print(f"Ramp type              : {ramp_type}")
print(f"Starting omega         : {omega_start:g} rad/s")
print(f"Ending omega           : {omega_end:g} rad/s")
print(f"Number of points       : {n_points}")
print(f"Time per point         : {time_per_point:g} s")

print("\n--- TIME ---")

print(f"Total measurement time : {total_time_s:g} s")
print(f"                        {total_time_min:.2f} min")
print(f"                        {total_time_h:.2f} h")

print("\n--- OMEGA ---")

if ramp_type == "linear":

    print(f"Delta omega            : {delta_omega[0]:.6g} rad/s")

else:

    print(f"Number of decades      : {n_decades:.3f}")
    print(f"Points per decade      : {points_per_decade:.2f}")
    print(f"Multiplication factor  : {ratio:.6g}")

print("\n--- OMEGA VALUES ---")

for i, omega in enumerate(omega_values, start=1):
    print(f"{i:3d} : {omega:.8g} rad/s")

print("=" * 60)