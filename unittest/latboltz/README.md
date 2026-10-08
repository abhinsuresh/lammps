# Multicomponent lattice-Boltzmann unit tests

This document describes the implemented moment-reconstruction and
gradient/Laplacian stencil tests.

| Unit | Status | Function under test |
| --- | --- | --- |
| Local moment reconstruction | Implemented | `FixLbMulticomponent::calc_moments()` |
| Local spatial derivatives | Implemented | `FixLbMulticomponent::calc_gradient_laplacian()` |

Source: [multicomponent implementation](../../src/LATBOLTZ/fix_lb_multicomponent.cpp),
[test source](test_fix_lb_multicomponent.cpp), and [CTest registration](CMakeLists.txt).

## 1. Implemented: moment reconstruction

### Physical meaning and equations

Each D3Q19 site stores three sets of 19 populations: `f_lb`, `g_lb`, and `k_lb`.
D3Q19 has one rest direction, six axial directions, and twelve face-diagonal
directions. Let $\mathbf{e}_i$ denote direction $i$, for $i=0,\ldots,18$.

The local density and composition fields are the zeroth moments:

$$
\rho = \sum_{i=0}^{18} f_i, \qquad
\phi = \sum_{i=0}^{18} g_i, \qquad
\psi = \sum_{i=0}^{18} k_i.
$$

The momentum density is the first moment of the fluid populations. Velocity
is momentum density divided by mass density:

$$
\mathbf{j} = \sum_{i=0}^{18} f_i\mathbf{e}_i, \qquad
\mathbf{u} = \frac{\mathbf{j}}{\rho}.
$$

The corresponding component densities are

$$
\rho_1 = \frac{\rho+\phi-\psi}{2}, \qquad
\rho_2 = \frac{\rho-\phi-\psi}{2}, \qquad
\rho_3 = \psi.
$$

These are component densities; divide by total density to obtain component
fractions. `calc_moments()` also stores the bulk pressure
$p_0(\rho,\phi,\psi)$ (Eq. (43) of Semprebon et al.). With
$\kappa_1=\kappa_2=\kappa_3=0$ every free-energy term vanishes and

$$
p_0 = c_s^2\rho = \frac{\rho}{3}.
$$

The test asserts $\rho$, $\phi$, $\psi$, all three velocity components, and
$p_0$ directly.

### Controlled input and expected result

The fixture creates the fix with `kappa1 0 kappa2 0 kappa3 0`, so the pressure
reference is independent of the free-energy expression. At the tested site the
test assigns every population from a table keyed by lattice vector, so the
inputs do not depend on the index order of `e19`. The test asserts that each of
the 19 D3Q19 directions is assigned exactly once.

| $\mathbf{e}_i$ | $f_i$ | $g_i$ | $k_i$ |
| --- | ---: | ---: | ---: |
| $(0,0,0)$ | 1.30 | 0.10 | 0.20 |
| $(+1,0,0)$ | 0.40 | 0.05 | 0 |
| $(-1,0,0)$ | 0.10 | 0 | 0.03 |
| $(0,+1,0)$ | 0.20 | 0 | 0 |
| $(0,-1,0)$ | 0.50 | $-0.04$ | 0 |
| $(0,0,+1)$ | 0.60 | 0 | 0.05 |
| $(0,0,-1)$ | 0.20 | 0 | 0 |
| $(+1,+1,0)$ | 0.10 | 0.03 | 0 |
| $(-1,-1,0)$ | 0.05 | 0 | 0 |
| $(+1,-1,0)$ | 0.07 | 0 | 0 |
| $(-1,+1,0)$ | 0.03 | 0 | 0.04 |
| $(+1,0,+1)$ | 0.08 | 0 | 0 |
| $(-1,0,-1)$ | 0.02 | 0.14 | 0 |
| $(+1,0,-1)$ | 0.04 | 0 | 0.06 |
| $(-1,0,+1)$ | 0.06 | 0 | 0 |
| $(0,+1,+1)$ | 0.09 | 0 | 0.02 |
| $(0,-1,-1)$ | 0.01 | 0 | 0 |
| $(0,+1,-1)$ | 0.05 | 0 | 0 |
| $(0,-1,+1)$ | 0.10 | 0.02 | 0 |

The independent, hand-calculated answers follow by grouping rest, axial, and
diagonal contributions. For the density,

$$
\rho = \underbrace{1.3}_{\text{rest}}
     + \underbrace{(0.4+0.1+0.2+0.5+0.6+0.2)}_{\text{axial}=2.0}
     + \underbrace{(0.25+0.20+0.25)}_{\text{diagonal}=0.7} = 4,
$$

where the diagonal sums are taken over the $xy$, $xz$, and $yz$ planes. The
composition fields are

$$
\phi = 0.10+0.05-0.04+0.03+0.14+0.02 = 0.30, \qquad
\psi = 0.20+0.03+0.05+0.04+0.06+0.02 = 0.40.
$$

Each momentum component is the axial difference plus the diagonal
differences in the two planes containing that axis:

$$
\begin{aligned}
j_x &= (0.4-0.1) + \underbrace{(0.17-0.08)}_{xy} + \underbrace{(0.12-0.08)}_{xz} = 0.43,\\
j_y &= (0.2-0.5) + \underbrace{(0.13-0.12)}_{xy} + \underbrace{(0.14-0.11)}_{yz} = -0.26,\\
j_z &= (0.6-0.2) + \underbrace{(0.14-0.06)}_{xz} + \underbrace{(0.19-0.06)}_{yz} = 0.61,
\end{aligned}
$$

so that

$$
\mathbf{u} = \mathbf{j}/\rho = (0.1075,\,-0.065,\,0.1525), \qquad
p_0 = 4/3.
$$

Why these inputs:

- Every diagonal is nonzero, so a loop that skips diagonals (e.g. stops at
  $i<7$) gives $\rho=3.3$ and fails.
- Opposite populations are unequal on axial and diagonal links, so a sign
  error in any `e19` component changes $\mathbf{j}$.
- $g$ and $k$ are spread over rest, axial, and diagonal directions, so summing
  only part of the stencil changes $\phi$ or $\psi$. The negative $g$ value
  catches clamping or absolute-value errors.
- $\rho\neq1$ exposes a missing division by density; $\phi\neq\psi$ exposes
  swapped composition fields.

These prescribed populations are numerical test inputs and need not be an
equilibrium distribution.

### Implementation walkthrough

1. **Access declaration:** `FixLbMulticomponent` declares
   `friend class FixLbMulticomponentTestAccess`. No production method is made
   public. The helper is defined in the test source in `LAMMPS_NS`.
2. **Fixture setup:** `FixLbMulticomponentTest::SetUp()` constructs LAMMPS,
   creates a periodic, atom-free $8\times8\times8$ box with lattice spacing
   and timestep equal to one, creates the fix with all $\kappa=0$, and
   validates its pointer.
3. **Choose a site:** `interior_site()` returns
   `{halo_extent[0], halo_extent[1], halo_extent[2]}`, i.e. local index
   $(2,2,2)$: the first owned site, just inside the halo. Its neighbours at
   index 1 are allocated halo cells.
4. **Supply inputs:** the test maps the table above onto D3Q19 indices, then
   `set_populations()` overwrites all 19 entries of each population family at
   that site. This removes dependence on the random mixture initialization at
   the tested site.
5. **Execute production code:** `calculate_moments()` calls `calc_moments()`
   and returns copies of `density_lb`, `phi_lb`, `psi_lb`, `pressure_lb`, and
   `u_lb`. It does not recalculate the expected answer.
6. **Compare:** `ReconstructsDensityCompositionAndVelocity` uses seven
   `EXPECT_NEAR` checks with absolute tolerance $10^{-13}$.
7. **Clean up:** `TearDown()` resets the owning `std::unique_ptr<LAMMPS>`.
   LAMMPS destroys its fix and associated resources.

No timesteps or halo exchanges are required: `calc_moments()` reads only
populations at the selected site. The shared test main initializes MPI;
the test runs with one MPI rank.

### Coverage and limitations

This test checks local zeroth and first moments over all 19 directions,
momentum signs on axial and diagonal links, field assignment, division by
density, and the ideal-gas part of the pressure. It does not check the
$\kappa$-dependent pressure terms, collision, streaming, global conservation,
or MPI communication.

### Mutation check

To confirm the test can fail, each of the following temporary edits was made
in `calc_moments()`, rebuilt, and run; every one made the test fail. The source
was then restored and both registered tests passed again.

| Temporary edit | Assertions that failed |
| --- | --- |
| Loop bound `i < numvel` → `i < 7` (skip diagonals) | all seven ($\rho=3.3$) |
| `j[1] += fi * e19[i][1]` → `e19[i][0]` | $u_y$ (0.1075 instead of $-0.065$) |
| Swap the `phi_lb` / `psi_lb` assignments | $\phi$ and $\psi$ |
| Remove `/ rho` from $u_x$ | $u_x$ (0.43 instead of 0.1075) |

This documents that validation run, not a guarantee about subsequent source
edits.

## 2. Implemented: gradient and Laplacian

### Physical relevance

`calc_equilibrium()` applies `calc_gradient_laplacian()` to $\rho$, $\phi$,
and $\psi$. Their derivatives enter chemical potentials and interfacial
stress terms. A derivative error can change diffusion and interface dynamics
even when total mass remains conserved.

### Discrete equations used by the code

For a scalar field $q$ at lattice site $\mathbf{n}$, the implementation uses

$$
(\nabla_{\mathrm{LB}}q)_\alpha(\mathbf{n})
= 3\sum_{i=0}^{18} w_i\,q(\mathbf{n}+\mathbf{e}_i)e_{i\alpha},
\qquad \alpha\in\{x,y,z\},
$$

$$
\nabla_{\mathrm{LB}}^2q(\mathbf{n})
= 6\sum_{i=0}^{18}w_i
\left[q(\mathbf{n}+\mathbf{e}_i)-q(\mathbf{n})\right].
$$

The D3Q19 weights are

$$
w_i =
\begin{cases}
1/3 & \text{rest direction},\\
1/18 & \text{six axial directions},\\
1/36 & \text{twelve face-diagonal directions}.
\end{cases}
$$

These are derivatives in lattice coordinates: the function contains no
explicit factors of `dx_lb`. For unchanged field units and physical spacing
$h$, physical derivatives are approximated by

$$
\nabla_{\mathrm{phys}}q \approx \frac{1}{h}\nabla_{\mathrm{LB}}q,
\qquad
\nabla_{\mathrm{phys}}^2q \approx \frac{1}{h^2}\nabla_{\mathrm{LB}}^2q.
$$

The unit tests use $h=1$. Conversion of field units is a separate
concern.

### Why polynomial fields provide independent answers

The stencil satisfies

$$
\sum_i w_i\mathbf{e}_i=0, \qquad
\sum_i w_i e_{i\alpha}e_{i\beta}=\frac{1}{3}\delta_{\alpha\beta},
$$

where $\delta_{\alpha\beta}$ is one for equal indices and zero otherwise.
Opposite directions also cancel odd-order terms. Expanding a quadratic field
about the tested site therefore gives its exact analytic gradient and
Laplacian, up to floating-point rounding. General smooth fields have
discretization error and should not be required to match continuum derivatives
to roundoff.

### Test cases

The test uses local coordinates $(X,Y,Z)$ relative to the chosen site, so the tested
site is $(0,0,0)$, and evaluates each polynomial over its complete stencil neighborhood.

| Field $q(X,Y,Z)$ | Expected gradient at the site | Expected Laplacian | Purpose |
| --- | --- | ---: | --- |
| $5$ | $(0,0,0)$ | 0 | Constant cancellation |
| $X$, $Y$, $Z$, separately | Corresponding unit vector | 0 | Each derivative direction |
| $5+2X-3Y+4Z$ | $(2,-3,4)$ | 0 | Signs and axis assignment |
| $X^2$, $Y^2$, $Z^2$, separately | $(0,0,0)$ | 2 | Curvature along each axis |
| $X^2+2Y^2+3Z^2$ | $(0,0,0)$ | 12 | Unequal curvature contributions |
| $(X+1)(Y+2)$ | $(2,1,0)$ | 0 | Mixed term in the xy plane |
| $(Y+1)(Z+2)$ | $(0,2,1)$ | 0 | Mixed term in the yz plane |
| $(Z+1)(X+2)$ | $(1,0,2)$ | 0 | Mixed term in the zx plane |

For example,

$$
q=5+2X-3Y+4Z+X^2+2Y^2+3Z^2
$$

has, at the tested site,

$$
\nabla q=(2,-3,4), \qquad \nabla^2q=2+4+6=12.
$$

This combined case is useful in addition to the isolated cases, which make
failures easier to diagnose. Every case in the table, plus this combined case
(13 polynomials), runs for each of the `rho`, `phi`, and `psi` field/output
triples, giving 39 parameterized tests.

**Negative control.** For $q=X^3$ the continuum gradient at the origin is
zero, but the second-order stencil gives

$$
(\nabla_{\mathrm{LB}}q)_x = 3\sum_i w_i e_{ix}^4
= 3\left(2\cdot\tfrac{1}{18} + 8\cdot\tfrac{1}{36}\right) = 1,
\qquad \nabla_{\mathrm{LB}}^2 q = 6\sum_i w_i e_{ix}^3 = 0.
$$

`CubicShowsSecondOrderTruncation` asserts these discrete values. It pins the
known truncation error and documents why only polynomials up to second order
are compared with continuum derivatives.

### Implementation walkthrough

1. **Access:** the existing `FixLbMulticomponentTestAccess` friend helper
   gains `Field {Rho, Phi, Psi}` and maps it to the same field/output triples
   that `calc_equilibrium()` passes: (`density_lb`, `density_gradient`,
   `laplace_rho`), (`phi_lb`, `phi_gradient`, `laplace_phi`), and
   (`psi_lb`, `psi_gradient`, `laplace_psi`). `calc_gradient_laplacian()`
   stays private.
2. **Fixture:** `FixLbDerivativeTest` derives from `FixLbMulticomponentTest`
   (same one-rank $8\times8\times8$ box) so the derivative tests get their own
   suite name for the CTest filter.
3. **Supply inputs:** `set_scalar_field()` evaluates the polynomial on the
   whole $3\times3\times3$ block around `interior_site()` = $(2,2,2)$. D3Q19
   reads only the centre, axial neighbours, and face diagonals; filling the
   block keeps setup simple. Neighbours at index 1 are allocated halo cells.
4. **Sentinel outputs:** before the call, `fill_derivatives()` sets every
   gradient component and Laplacian on that block to $10^3$. This replaces the
   planned "call twice" check: if the function accumulated into its outputs
   instead of resetting them, the tested site would be off by $10^3$.
5. **Execute production code:** `calculate_derivatives()` calls
   `calc_gradient_laplacian()` with the selected triple and returns copies of
   the three gradient components and the Laplacian at the site.
6. **Compare:** `FixLbDerivativePolynomialTest.MatchesAnalyticDerivatives` is
   parameterized over (field, polynomial) and checks the four values against
   the analytic derivatives with absolute tolerance $10^{-13}$. It also checks
   that all 26 neighbours still hold the sentinel, i.e. the function writes
   only to the tested site.

The local polynomials need not be periodic over the whole box because this
test neither exchanges halos nor advances the fluid. Boundary and MPI halo
behaviour require separate tests. Because the helper passes the triple
explicitly, the test does not validate the wiring inside `calc_equilibrium()`;
that needs a separate integration test.

### Mutation check

Each of the following temporary edits to `calc_gradient_laplacian()` was
rebuilt and run against the 40 derivative tests; every one made tests fail.
The source was then restored and all three registered tests passed again.

| Temporary edit | Tests failed (of 40) |
| --- | ---: |
| Gradient factor `3.` → `2.` | 25 |
| Laplacian factor `6.` → `3.` | 15 |
| Gradient uses `e19[i][0]` for every direction | 25 |
| Laplacian output not reset to zero | 39 |
| Gradient outputs not reset to zero | 39 |
| Laplacian without subtracting the centre value | 18 |
| $x$ neighbour taken as `x - e19[i][0]` | 16 |

The unaffected cases pass as expected; e.g. the Laplacian factor does not
change any case whose Laplacian is zero. This documents that validation run,
not a guarantee about subsequent source edits.

## 3. Build and run the implemented tests

From a clean CMake-compatible source tree with an MPI toolchain:

```sh
cmake -S cmake -B build-latboltz -D BUILD_MPI=ON -D PKG_LATBOLTZ=ON -D ENABLE_TESTING=ON
cmake --build build-latboltz --target test_fix_lb_multicomponent -j 4
ctest --test-dir build-latboltz -R '^FixLBMulticomponent:moments$' --output-on-failure
ctest --test-dir build-latboltz -R '^FixLBMulticomponent:derivatives$' --output-on-failure
```

Run all three registered LATBOLTZ tests with:

```sh
ctest --test-dir build-latboltz -R '^FixLBMulticomponent:' --output-on-failure
```

`FixLBMulticomponent:moments` and `FixLBMulticomponent:derivatives` use one MPI
rank; the `FixLBMulticomponent:plain` conservation test uses two. GoogleTest
filters select the intended tests from their shared executable:

| CTest name | GoogleTest filter |
| --- | --- |
| `FixLBMulticomponent:plain` | `FixLBMulticomponent.plain` |
| `FixLBMulticomponent:moments` | `FixLbMulticomponentTest.*` |
| `FixLBMulticomponent:derivatives` | `FixLbDerivativeTest.*:*/FixLbDerivativePolynomialTest.*` |

A filter that matches no test still exits successfully in GoogleTest 1.12, so
after renaming a fixture, run with `-V` and confirm the test count is non-zero.
The executable still receives the existing YAML file because the shared
force-style test main expects it; the unit tests use their C++ tolerances, not
the YAML epsilon.
