import os
import sys
import math

import numpy as np
import matplotlib
matplotlib.use("Agg")  # без GUI
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

import pycosmos


# ---------------------------------------------------------------- константы
MU_EARTH = 398600.4418          # км^3/с^2
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")


# ---------------------------------------------------------------- утилиты
def orbital_period(a, mu=MU_EARTH):
    return 2.0 * math.pi * math.sqrt(a ** 3 / mu)


def specific_energy(r, v, mu=MU_EARTH):
    return 0.5 * float(np.dot(v, v)) - mu / float(np.linalg.norm(r))


def angular_momentum(r, v):
    return np.cross(r, v)


def propagate(r0, v0, t_array, mu=MU_EARTH):
    """Распространяет орбиту на массив времён. Возвращает (rs, vs)."""
    n = len(t_array)
    rs = np.empty((n, 3))
    vs = np.empty((n, 3))
    for i, t in enumerate(t_array):
        r, v = pycosmos.predict(r0, v0, 0.0, float(t), mu)
        rs[i] = np.asarray(r)
        vs[i] = np.asarray(v)
    return rs, vs


def sample_orbit(kep: dict, mu=MU_EARTH, n_points=400):
    """Сэмплирует позиции по всей орбите, меняя истинную аномалию nu."""
    positions = np.empty((n_points, 3))
    for i, nu in enumerate(np.linspace(0, 2 * math.pi, n_points, endpoint=False)):
        k = pycosmos.Kepler(
            kep["a"], kep["e"], kep["inc"], kep["raan"], kep["argp"], nu
        )
        r, _ = pycosmos.kep_to_dec(k, mu)
        positions[i] = np.asarray(r)
    return positions


# ---------------------------------------------------------------- тесты-ассерты
def test_kepler_class():
    k = pycosmos.Kepler()
    assert isinstance(k, pycosmos.Kepler)

    k = pycosmos.Kepler(7000.0, 0.1, 0.2, 0.3, 0.4, 0.5)
    assert abs(k.a - 7000.0) < 1e-12
    k.a = 8000.0
    assert abs(k.a - 8000.0) < 1e-12

    s = repr(k)
    assert "Kepler" in s and "a=" in s


def test_simple_functions():
    # mean_motion
    a = 7000.0
    n = pycosmos.mean_motion(a, MU_EARTH)
    assert abs(n - math.sqrt(MU_EARTH / a ** 3)) < 1e-15

    # solve_kepler для e=0
    M = 0.7
    E = pycosmos.solve_kepler(M, 0.0)
    assert abs(E - M) < 1e-12

    # solve_kepler удовлетворяет уравнению M = E - e sin E
    M, e = 1.2, 0.4
    E = pycosmos.solve_kepler(M, e)
    assert abs(M - (E - e * math.sin(E))) < 1e-10

    # ist_v_ecrt / ecsrt_v_ist — обратные
    e = 0.3
    E = 0.8
    nu = pycosmos.ecsrt_v_ist(E, e)
    E2 = pycosmos.ist_v_ecrt(nu, e)
    assert abs(math.sin(E) - math.sin(E2)) < 1e-10
    assert abs(math.cos(E) - math.cos(E2)) < 1e-10


def test_kep_dec_roundtrip():
    cases = [
        dict(a=7000.0, e=0.0, inc=0.0, raan=0.0, argp=0.0, nu=0.0),
        dict(a=10000.0, e=0.3,
             inc=math.radians(45.0),
             raan=math.radians(60.0),
             argp=math.radians(30.0),
             nu=math.radians(10.0)),
    ]
    for case in cases:
        k = pycosmos.Kepler(**case)
        r, v = pycosmos.kep_to_dec(k, MU_EARTH)
        k2 = pycosmos.dec_to_kep(np.array(r), np.array(v), MU_EARTH)
        assert abs(k2.a - k.a) / k.a < 1e-8, f"a drift: {k2.a} vs {k.a}"
        assert abs(k2.e - k.e) < 1e-8, f"e drift: {k2.e} vs {k.e}"
        assert abs(k2.inc - k.inc) < 1e-8, f"inc drift"
        assert abs(math.cos(k2.raan) - math.cos(k.raan)) < 1e-6


def test_predict():
    kep = dict(a=6778.0, e=0.0005,
               inc=math.radians(51.6),
               raan=math.radians(30.0),
               argp=0.0, nu=0.0)
    k = pycosmos.Kepler(**kep)
    r0, v0 = pycosmos.kep_to_dec(k, MU_EARTH)
    r0 = np.asarray(r0)
    v0 = np.asarray(v0)
    T = orbital_period(k.a)

    # t=0 → исходное состояние
    r, v = pycosmos.predict(r0, v0, 0.0, 0.0, MU_EARTH)
    assert np.allclose(r, r0, atol=1e-6)
    assert np.allclose(v, v0, atol=1e-8)

    # t=T → почти исходное состояние
    r, v = pycosmos.predict(r0, v0, 0.0, T, MU_EARTH)
    assert np.allclose(r, r0, rtol=1e-4, atol=1e-2), \
        f"r после периода: {r}, ожидали {r0}"


def test_energy_moment_conservation(rs, vs):
    energies = np.array([specific_energy(r, v) for r, v in zip(rs, vs)])
    eps0 = energies[0]
    drift_e = np.abs((energies - eps0) / eps0).max()
    assert drift_e < 1e-6, f"Energy drift = {drift_e:.2e}"

    h = np.array([angular_momentum(r, v) for r, v in zip(rs, vs)])
    h_norm = np.linalg.norm(h, axis=1)
    h0 = h_norm[0]
    drift_h = np.abs(h_norm - h0).max() / h0
    assert drift_h < 1e-6, f"|h| drift = {drift_h:.2e}"

    cos_angle = h @ h[0] / (h_norm * h0)
    assert cos_angle.min() > 1 - 1e-8, "направление h уплыло"

    return energies, h_norm


def test_kepler_equal_areas(rs, t):
    dt = t[1] - t[0]
    dr = np.diff(rs, axis=0)
    r_mid = 0.5 * (rs[:-1] + rs[1:])
    area_rates = 0.5 * np.linalg.norm(np.cross(r_mid, dr), axis=1) / dt
    rel_std = area_rates.std() / area_rates.mean()
    assert rel_std < 1e-3, f"Area rate variation = {rel_std:.2e}"


def test_orbit_geometry(positions):
    # радиус между перицентром и апоцентром
    a, e = 10000.0, 0.3
    r = np.linalg.norm(positions, axis=1)
    assert abs(r.min() - a * (1 - e)) < 1e-2
    assert abs(r.max() - a * (1 + e)) < 1e-2

    # все точки в одной плоскости (нормали коллинеарны)
    h = np.cross(positions, np.roll(positions, -1, axis=0))
    h_norm = h / np.linalg.norm(h, axis=1, keepdims=True)
    dots = np.abs(h_norm @ h_norm[0])
    assert dots.min() > 1 - 1e-6, "орбита не плоская"


# ---------------------------------------------------------------- графики
def plot_energy_and_moment(t, energies, h_norm):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    axes[0].plot(t, energies - energies[0], color="C0")
    axes[0].set_ylabel("ΔE, км²/с²")
    axes[0].set_title("Сохранение удельной энергии")
    axes[0].grid(True)

    axes[1].plot(t, h_norm, color="C1")
    axes[1].set_xlabel("Время, с")
    axes[1].set_ylabel("|h|, км²/с")
    axes[1].set_title("Сохранение модуля орбитального момента")
    axes[1].grid(True)

    path = os.path.join(OUTPUT_DIR, "energy_moment.png")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"  [saved] {path}")


def plot_orbit_3d(positions):
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    x, y, z = positions.T

    fig = plt.figure(figsize=(14, 6))

    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax.plot(x, y, z, lw=1.5)
    ax.scatter([0], [0], [0], color="blue", marker="o", s=80, label="Земля")
    ax.set_xlabel("X, км")
    ax.set_ylabel("Y, км")
    ax.set_zlabel("Z, км")
    ax.set_title("3D орбита")
    ax.legend()
    ax.set_box_aspect([1, 1, 0.5])

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.plot(x, y, lw=1.5)
    ax2.scatter([0], [0], color="blue", marker="o", s=80)
    ax2.set_aspect("equal")
    ax2.set_xlabel("X, км")
    ax2.set_ylabel("Y, км")
    ax2.set_title("Проекция на XY")
    ax2.grid(True)

    path = os.path.join(OUTPUT_DIR, "orbit_3d.png")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"  [saved] {path}")


def plot_orbit_evolution():
    """МКС-подобная орбита за 1.5 периода, раскраска по времени."""
    kep = dict(a=6778.0, e=0.0005,
               inc=math.radians(51.6),
               raan=math.radians(30.0),
               argp=0.0, nu=0.0)
    k = pycosmos.Kepler(**kep)
    r0, v0 = pycosmos.kep_to_dec(k, MU_EARTH)
    r0 = np.asarray(r0)
    v0 = np.asarray(v0)
    T = orbital_period(k.a)

    t = np.linspace(0, 1.5 * T, 600)
    positions = np.empty((len(t), 3))
    for i, ti in enumerate(t):
        r, _ = pycosmos.predict(r0, v0, 0.0, float(ti), MU_EARTH)
        positions[i] = np.asarray(r)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")

    sc = ax.scatter(
        positions[:, 0], positions[:, 1], positions[:, 2],
        c=t / T, cmap="viridis", s=4,
    )
    ax.scatter([0], [0], [0], color="red", marker="o", s=100, label="Земля")
    ax.set_xlabel("X, км")
    ax.set_ylabel("Y, км")
    ax.set_zlabel("Z, км")
    ax.set_title("Эволюция орбиты за 1.5 периода")
    cbar = fig.colorbar(sc, ax=ax, shrink=0.6)
    cbar.set_label("t / T")
    ax.legend()

    path = os.path.join(OUTPUT_DIR, "orbit_evolution.png")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    print(f"  [saved] {path}")


# ---------------------------------------------------------------- main
def main():
    print("=" * 60)
    print("pycosmos self-test")
    print("=" * 60)

    print("\n[1] Класс Kepler и простые функции ...")
    test_kepler_class()
    test_simple_functions()
    print("  OK")

    print("\n[2] Round-trip kep <-> dec ...")
    test_kep_dec_roundtrip()
    print("  OK")

    print("\n[3] Предсказание орбиты (predict) ...")
    test_predict()
    print("  OK")

    print("\n[4] Распространение эллиптической орбиты на 3 периода ...")
    k = pycosmos.Kepler(
        10000.0, 0.3,
        math.radians(45.0), math.radians(60.0),
        math.radians(30.0), 0.0,
    )
    r0, v0 = pycosmos.kep_to_dec(k, MU_EARTH)
    r0 = np.asarray(r0)
    v0 = np.asarray(v0)
    T = orbital_period(k.a)
    t = np.linspace(0, 3 * T, 1500)
    rs, vs = propagate(r0, v0, t)
    print(f"  точек: {len(t)}, период: {T:.1f} с")

    print("\n[5] Сохранение энергии и момента ...")
    energies, h_norm = test_energy_moment_conservation(rs, vs)
    drift_e = np.abs((energies - energies[0]) / energies[0]).max()
    drift_h = np.abs(h_norm - h_norm[0]).max() / h_norm[0]
    print(f"  max |ΔE/E| = {drift_e:.2e}")
    print(f"  max |Δ|h||/|h| = {drift_h:.2e}")
    print("  OK")

    print("\n[6] Второй закон Кеплера (равные площади) ...")
    test_kepler_equal_areas(rs, t)
    print("  OK")

    print("\n[7] Геометрия орбиты (границы радиуса, плоскость) ...")
    kep_ell = dict(a=10000.0, e=0.3,
                   inc=math.radians(45.0),
                   raan=math.radians(60.0),
                   argp=math.radians(30.0),
                   nu=0.0)
    positions = sample_orbit(kep_ell)
    test_orbit_geometry(positions)
    print("  OK")

    print("\n[8] Построение графиков ...")
    plot_energy_and_moment(t, energies, h_norm)
    plot_orbit_3d(positions)
    plot_orbit_evolution()

    print("\n" + "=" * 60)
    print("Все проверки пройдены.")
    print(f"Графики: {OUTPUT_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        print(f"\n!!! ПРОВАЛ ТЕСТА: {e}", file=sys.stderr)
        sys.exit(1)