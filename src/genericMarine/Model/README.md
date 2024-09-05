# ModelZonalAdvection

This sample YAML configuration provides an example of the ZonalAdvection model
configuration for simulating zonal (east-west) advection of variables like sea
surface height anomaly. This setup includes parameters for advection and
diffusion, offering a reasonable default configuration to produce stable results.

```yaml
model:
  name: ZonalAdvection
  tstep: PT1H
  variables: [seaSurfaceHeightAnomaly]
  advection:
    speed:
      latitude: [ 60.0,  30.0,  10.0,  2.0, 1, -1, -2.0, -10.0, -30.0, -60.0]
      value:    [ -0.0, -0.05, -0.15, -0.8, 2.7, 2.7, -0.8, -0.15, -0.05,   0.0]
    coastal damping:
      distance: 100.0e3
      amount: 0.75
    boundary condition:
      a: 0.0
      b: 0.0
    asselin filter: 0.2
  diffusion:
    coefficient smoothing: 1
    #Kh: 1.0e2
    Kh_smag scale: 2.0e4
    Kh_smag max: 1.0e5
    Ah: 2.0e10
```

## Advection

The horizontal advection equation:

```math
\frac{\partial \phi}{\partial t} = -u \frac{\partial \phi}{\partial x} - v\frac{\partial\phi}{\partial y}
```

describes the rate of change of the scalar $\phi$ in terms of the velocity
components $u$ and $v$. The underlying `ModelAdvectionBase` class is capable of
full 2D advection, however, only zonal velocities are specified for this
`ModelZonalAdvection` class.

The equation has been discretized using a leapfrog time scheme and
centered-in-space diferences:

```math
\frac{\phi^{n+1}_{i,j} - \phi^{n-1}_{i,j}}{2 \Delta t} =  - \frac{u}{2 \Delta x}(\phi^{n}_{i+1,j}-\phi^{n}_{i-1,j}) - \frac{v}{2 \Delta y}(\phi^{n}_{i,j+1}-\phi^{n}_{i,j-1})
```

This leapfrog scheme is prone to oscillations, also known as the "computational
mode." To prevent this, an Asselin filter is used. After calculating the
advection tendencies, the current timestep is smoothed with the previous and
next timesteps to remove spurious oscillations. The smoothing operation is:

```math
\phi^{n} \rightarrow \phi^{n} + \epsilon ( \phi^{n+1} - 2\phi^{n} + \phi^{n-1})
```

where $\epsilon$ is the Asselin filter parameter, controlling the amount of
smoothing applied.

### Parameters

- `speed` : a series of latitude/value values that specifies the zonal advection
  velocity at various latitudes. The values between the given latitudes are
  linearly interpolated.
- `coastal damping` : the speed from above can be damped in areas near the
  coast. A small amount of damping is usually beneficial for the stability of
  the advection model.
  - `distance` : Distance from the coast over which the speed is linearly
    damped.
  - `amount` : The maximum amount to damp the speeds. 1.0 is fully damped at the
    coast (speed is 0), 0.0 is no impact.
- `boundary condition` : Sets the Dirichlet boundary condition for incoming flow. (The
  outgoing flow is a fixed Neumann condition to maintain a zero spatial
  derivative). The value at the boundary where there is inflow $= a
  \cdot \phi_i + b$, where $\phi_i$ is the value of the nearest valid neighbor.
  For example, the two extreme cases are:
  - To set inflow boundary condition to a constant value, set `a=0, b=<constant value>`.
  - To set inflow boundary condition to be equal to the nearest neighbor, set
    `a=1, b=0`
- `asselin filter` : Strength of the Asselin time filter to mitigate advection
  computational mode **(default: 0.2)**

## Diffusion
Diffusion in this model is used to smooth out small-scale features and stabilize
the simulation. Three types of diffusion are implemented, each using different
orders of the Laplacian operator.

### 1. Harmonic Diffusion
The harmonic diffusion equation is given by:

```math
\frac{\partial \phi}{\partial t} = K_h \nabla^2\phi
```

This equation describes second-order diffusion, where $K_h$ is the diffusion
coefficient and $\nabla^2\phi$ is the Laplacian of the scalar field. The
equation is discretized as:

```math
\frac{\phi^{n+1} -\phi^{n-1}}{2 \Delta t} = K_h (\frac{\phi^{n-1}_{i+1,j} -2\phi^{n-1}_{i,j}+\phi^{n-1}_{i-1,j}}{\Delta x^2} + \frac{\phi^{n-1}_{i,j+1} - 2\phi^{n-1}_{i,j}+\phi^{n-1}_{i,j-1}}{\Delta y^2})
```

To be compatible with the leapfrog time scheme of the advection, note that the
previous timestep is used for calculating the spatial derivatives for diffusion

To ensure numerical stability, the diffusion coefficient $K_h$ is automatically
adjusted to respect the CFL (Courant–Friedrichs–Lewy) condition:

```math
K_h \le 0.5 \frac{\Delta x ^2}{\Delta y}
```

Note, the basic harmonic diffusion is usually overly diffusive, you'll probably
want to use some combination of the other two types of diffusion.

### 2. Smagorinsky Diffusion
Smagorinsky diffusion is a nonlinear form of harmonic diffusion where the
diffusion coefficient $K_{h\_smag}$ depends on the local flow's strain rate:

```math
K_{h\_smag} = l_s^2 \sqrt{T^2 + S^2}
```

```math
T = \frac{\partial u}{\partial x} - \frac{\partial v}{\partial y}
```

```math
S = \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x}
```

Where $T$ is the horizontal tension strain and $S$ is the horizontal shearing
strain. The length scale $l_s$ controls the strength of the Smagorinsky
diffusion. This type of diffusion is less diffusive than basic harmonic
diffusion in regions with weak velocity gradients but stronger where velocity gradients are
larger. This will help the advection from producing unwanted artifacts in areas
of high shear.

The resulting Smagorinsky diffusion parameter $K_{h\_smag}$ is added to the base $K_h$.

### 3. Biharmonic Diffusion
Biharmonic diffusion is a higher-order diffusion process:

```math
\frac{\partial \phi} {\partial t} = A_h\nabla^4 \phi
```

Where $A_h$ is the biharmonic diffusion coefficient, and $\nabla^4\phi$ is the
fourth-order Laplacian operator (Laplacian of the Laplacian). This diffusion
type is designed to remove small-scale noise more efficiently than harmonic
diffusion while preserving larger-scale features.

Similar to harmonic diffusion, the biharmonic coefficient $A_h$ is constrained
by the CFL condition for numerical stability:

```math
A_h \le \frac{1}{16} \frac{\Delta x ^ 4} {\Delta t}
```

### Parameters

- `coefficient smoothing` : The number of iterations of smoothing for the $K_h$
  and $A_h$ parameters to ensure diffusion  coefficients transition smoothly
  across the domain **(default: 1)**
- `Kh` : The harmonic diffusion coefficient **(units: $m^2/s$, default: 0.0)**
- `Kh_smag scale`: The scaling factor for the Smagorinsky diffusion **(units:
  $m$, default: 0.0)**
- `Kh_smag max`: Maximum value for Smagorinsky diffusion $K_{h\_smag}$ to
  prevent excessively large values **(units: $m^4/s$, default: 0.0)**
- `Ah` : The biharmonic diffusion coefficient **(units: $m^4/s$, default: 0.0)**