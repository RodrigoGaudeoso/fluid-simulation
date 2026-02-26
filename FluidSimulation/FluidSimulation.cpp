#include "FluidSimulation.h"
#include "raylib.h"
#include <vector>
#include <algorithm>
#include <cmath>

inline int FluidSimulation::IX(int x, int y, int N) {
    return x + y * N;
}

const std::vector<float>& FluidSimulation::GetDensity() const { return density; }

const int FluidSimulation::GetSize() const { return N; }

void FluidSimulation::AddDensity(int x, int y, float amount) {
    density[IX(x, y, N)] += amount;
    density[IX(x + 1, y, N)] += amount;
    density[IX(x - 1, y, N)] += amount;
    density[IX(x, y + 1, N)] += amount;
    density[IX(x, y - 1, N)] += amount;
}

void FluidSimulation::AddVelocity(int x, int y, float ax, float ay) {
    const int idx = IX(x, y, N);
    Vx[idx] += ax;
    Vy[idx] += ay;
}

void FluidSimulation::Step() {
    Diffuse(1, Vx0, Vx, visc, dt, iter, N);
    Diffuse(2, Vy0, Vy, visc, dt, iter, N);

    std::vector<float> p(N * N, 0.0f), div(N * N, 0.0f);
    Project(Vx0, Vy0, p, div, iter, N);

    Advect(1, Vx, Vx0, Vx0, Vy0, dt, N);
    Advect(2, Vy, Vy0, Vx0, Vy0, dt, N);

    p.assign(N * N, 0.0f);
    div.assign(N * N, 0.0f);
    Project(Vx, Vy, p, div, iter, N);

    Diffuse(0, s, density, diff, dt, iter, N);
    Advect(0, density, s, Vx, Vy, dt, N);
}

void FluidSimulation::Clear() {
    std::fill(density.begin(), density.end(), 0.0f);
    std::fill(Vx.begin(), Vx.end(), 0.0f);
    std::fill(Vy.begin(), Vy.end(), 0.0f);
}

void FluidSimulation::SetBnd(int b, std::vector<float>& x, int N) {
    for (int y = 1; y < N - 1; y++) {
        x[IX(0, y, N)] = (b == 1) ? -x[IX(1, y, N)] : x[IX(1, y, N)];
        x[IX(N - 1, y, N)] = (b == 1) ? -x[IX(N - 2, y, N)] : x[IX(N - 2, y, N)];
    }

    for (int x0 = 1; x0 < N - 1; x0++) {
        x[IX(x0, 0, N)] = (b == 2) ? -x[IX(x0, 1, N)] : x[IX(x0, 1, N)];
        x[IX(x0, N - 1, N)] = (b == 2) ? -x[IX(x0, N - 2, N)] : x[IX(x0, N - 2, N)];
    }

    x[IX(0, 0, N)] = 0.5f * (x[IX(1, 0, N)] + x[IX(0, 1, N)]);
    x[IX(0, N - 1, N)] = 0.5f * (x[IX(1, N - 1, N)] + x[IX(0, N - 2, N)]);
    x[IX(N - 1, 0, N)] = 0.5f * (x[IX(N - 2, 0, N)] + x[IX(N - 1, 1, N)]);
    x[IX(N - 1, N - 1, N)] = 0.5f * (x[IX(N - 2, N - 1, N)] + x[IX(N - 1, N - 2, N)]);
}

void FluidSimulation::LinSolve(int b,
    std::vector<float>& x,
    const std::vector<float>& x0,
    float a, float c,
    int iter, int N) {

    const float cRecip = 1.0f / c;

    for (int it = 0; it < iter; it++) {
        for (int y = 1; y < N - 1; y++) {
            for (int x1 = 1; x1 < N - 1; x1++) {
                x[IX(x1, y, N)] =
                    (x0[IX(x1, y, N)] +
                        a * (x[IX(x1 + 1, y, N)] +
                            x[IX(x1 - 1, y, N)] +
                            x[IX(x1, y + 1, N)] +
                            x[IX(x1, y - 1, N)])) * cRecip;
            }
        }
        SetBnd(b, x, N);
    }
}

void FluidSimulation::Diffuse(int b,
    std::vector<float>& x,
    const std::vector<float>& x0,
    float diff, float dt,
    int iter, int N) {

    const float a = dt * diff * (N - 2) * (N - 2);
    LinSolve(b, x, x0, a, 1.0f + 4.0f * a, iter, N);
}

void FluidSimulation::Project(std::vector<float>& velocX,
    std::vector<float>& velocY,
    std::vector<float>& p,
    std::vector<float>& div,
    int iter, int N) {

    for (int y = 1; y < N - 1; y++) {
        for (int x = 1; x < N - 1; x++) {
            div[IX(x, y, N)] = -0.5f * (
                velocX[IX(x + 1, y, N)] - velocX[IX(x - 1, y, N)] +
                velocY[IX(x, y + 1, N)] - velocY[IX(x, y - 1, N)]
                ) / static_cast<float>(N);

            p[IX(x, y, N)] = 0.0f;
        }
    }

    SetBnd(0, div, N);
    SetBnd(0, p, N);
    LinSolve(0, p, div, 1.0f, 4.0f, iter, N);

    for (int y = 1; y < N - 1; y++) {
        for (int x = 1; x < N - 1; x++) {
            velocX[IX(x, y, N)] -= 0.5f * (p[IX(x + 1, y, N)] - p[IX(x - 1, y, N)]) * N;
            velocY[IX(x, y, N)] -= 0.5f * (p[IX(x, y + 1, N)] - p[IX(x, y - 1, N)]) * N;
        }
    }

    SetBnd(1, velocX, N);
    SetBnd(2, velocY, N);
}

void FluidSimulation::Advect(int b,
    std::vector<float>& d,
    const std::vector<float>& d0,
    const std::vector<float>& velocX,
    const std::vector<float>& velocY,
    float dt, int N) {

    const float dt0 = dt * (N - 2);
    const float Nf = static_cast<float>(N);

    for (int y = 1; y < N - 1; y++) {
        for (int x = 1; x < N - 1; x++) {
            const int idx = IX(x, y, N);

            float xf = static_cast<float>(x - dt0 * velocX[idx]);
            float yf = static_cast<float>(y - dt0 * velocY[idx]);

            xf = std::clamp(xf, 0.5f, Nf - 1.5f);
            yf = std::clamp(yf, 0.5f, Nf - 1.5f);

            int x0 = static_cast<int>(std::floor(xf));
            int x1 = x0 + 1;
            int y0 = static_cast<int>(std::floor(yf));
            int y1 = y0 + 1;

            float s1 = xf - x0, s0 = 1.0f - s1;
            float t1 = yf - y0, t0 = 1.0f - t1;

            d[idx] =
                s0 * (t0 * d0[IX(x0, y0, N)] + t1 * d0[IX(x0, y1, N)]) +
                s1 * (t0 * d0[IX(x1, y0, N)] + t1 * d0[IX(x1, y1, N)]);
        }
    }

    SetBnd(b, d, N);
}