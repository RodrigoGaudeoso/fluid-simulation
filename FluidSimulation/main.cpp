#include "FluidSimulation.h"
#include "raylib.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

static unsigned char ClampToByte(float v) {
    v = std::clamp(v, 0.0f, 255.0f);
    return static_cast<unsigned char>(v);
}
static unsigned char GetColor(unsigned char c, unsigned char rgb) {
    return static_cast<unsigned char>(static_cast<int>(c) - static_cast<int>(rgb));
}

int main() {
    const int screenW = 900;
    const int screenH = 900;
    const int N = 150;
    const int cellPx = screenW / N;

    InitWindow(screenW, screenH, "Fluid Simulation");
    SetTargetFPS(60);
                              
    FluidSimulation fluidSim(N, 0.1f, 0.001f, 0.0001f, 4); //size //timestep //diffusion //viscosity //iterations

    while (!WindowShouldClose()) {
        Vector2 m = GetMousePosition();
        int gx = std::clamp(static_cast<int>(m.x / cellPx), 1, N - 2);
        int gy = std::clamp(static_cast<int>(m.y / cellPx), 1, N - 2);
    
        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON)) {
            fluidSim.AddDensity(gx, gy, 240.0f);

            Vector2 md = GetMouseDelta();
            fluidSim.AddVelocity(gx, gy, md.x * 0.03f, md.y * 0.03f);
        }

        if (IsKeyPressed(KEY_C)) {
            fluidSim.Clear();
        }

        fluidSim.Step();

        BeginDrawing();
        ClearBackground(RAYWHITE);

        const auto& dens = fluidSim.GetDensity();
        unsigned char color = static_cast<char>(GetRandomValue(0, 255));

        int midX = N / 2;
        int midY = N / 2;
        float midD = dens[fluidSim.IX(midX, midY, N)];
 
        for (int y = 0; y < N; y++) {
            for (int x = 0; x < N; x++) {
                float d = dens[fluidSim.IX(x, y, N)];
                unsigned char c = ClampToByte(d * 5.0f);

                DrawRectangle(x * cellPx, y * cellPx, cellPx, cellPx, Color{ c, c, c, 255 });
            }
        }

        DrawText("Fluid Simulation using Navier-Stokes equations | LMB: add dye | C: clear", 10, 10, 18, RAYWHITE);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
