#pragma once
#include <vector>

class FluidSimulation {
public:
    FluidSimulation(int size, float timestep, float diffusion, float viscosity, int iterations = 4)
        :   N(size), 
            dt(timestep), 
            diff(diffusion), 
            visc(viscosity), 
            iter(iterations),
            s(N* N, 0.0f),
            density(N* N, 0.0f),
            Vx(N* N, 0.0f),
            Vy(N* N, 0.0f),
            Vx0(N* N, 0.0f),
            Vy0(N* N, 0.0f)
    {};
    ~FluidSimulation() = default;
 
    inline int IX(int x, int y, int N);
    const std::vector<float>& GetDensity() const;
    const int GetSize() const;
    void AddDensity(int x, int y, float amount);
    void AddVelocity(int x, int y, float ax, float ay);
    void Step();
    void Clear();
   

    void SetBnd(int b, std::vector<float>& x, int N);
    void LinSolve(int b, std::vector<float>& x, const std::vector<float>& x0, float a, float c, int iter, int N);
    void Diffuse(int b, std::vector<float>& x, const std::vector<float>& x0, float diff, float dt, int iter, int N);
    void Project(std::vector<float>& velocX, std::vector<float>& velocY, std::vector<float>& p, std::vector<float>& div, int iter, int N);
    void Advect(int b, std::vector<float>& d, const std::vector<float>& d0, const std::vector<float>& velocX, const std::vector<float>& velocY, float dt, int N);

private:
    int N = 0;
    float dt = 0.1f;
    float diff = 0.0f;
    float visc = 0.0f;
    int iter = 4;

    std::vector<float> s, density;
    std::vector<float> Vx, Vy;
    std::vector<float> Vx0, Vy0;
};
