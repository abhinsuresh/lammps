# Multicomponent lattice-Boltzmann unit tests

This document describes the implemented moment-reconstruction test and the
proposed first-priority test of the gradient and Laplacian stencil.

| Unit | Status | Function under test |
| --- | --- | --- |
| Local moment reconstruction | Implemented | `FixLbMulticomponent::calc_moments()` |
| Local spatial derivatives | Proposed; not implemented | `FixLbMulticomponent::calc_gradient_laplacian()` |

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
fractions. The current test asserts $\rho$, $\phi$, $\psi$, and all three
velocity components directly.

### Controlled input and expected result

All populations initially equal zero. The test then assigns:

| Direction | $f_i$ | $g_i$ | $k_i$ |
| --- | ---: | ---: | ---: |
| Rest | 2.0 | 0.3 | 0.4 |
| $+x$ | 0.4 | 0 | 0 |
| $-x$ | 0.1 | 0 | 0 |
| $+y$ | 0.2 | 0 | 0 |
| $-y$ | 0.5 | 0 | 0 |
| $+z$ | 0.6 | 0 | 0 |
| $-z$ | 0.2 | 0 | 0 |
| All diagonals | 0 | 0 | 0 |

Thus the independent, hand-calculated answers are

$$
\rho = 2+0.4+0.1+0.2+0.5+0.6+0.2 = 4,
\qquad \phi=0.3, \qquad \psi=0.4,
$$

$$
\mathbf{j} = (0.4-0.1,\;0.2-0.5,\;0.6-0.2)
             = (0.3,-0.3,0.4),
$$

$$
\mathbf{u} = (0.075,-0.075,0.100).
$$

Unequal opposite populations check momentum signs. Density differs from one
to expose missing velocity normalization. Distinct values of $\phi$ and $\psi$
help expose swapped composition fields. These prescribed populations are
numerical test inputs and need not be an equilibrium distribution.

### Implementation walkthrough

1. **Access declaration:** `FixLbMulticomponent` declares
   `friend class FixLbMulticomponentTestAccess`. No production method is made
   public. The helper is defined in the test source in `LAMMPS_NS`.
2. **Fixture setup:** `LBMomentReconstruction::SetUp()` constructs LAMMPS,
   creates a periodic, atom-free $8\times8\times8$ box with lattice spacing
   and timestep equal to one, creates the fix, and validates its pointer.
3. **Choose a site:** `owned_site()` returns the midpoint of the local
   allocated lattice. For this fixture it is inside the owned region,
   away from halo cells.
4. **Supply inputs:** `set_populations()` overwrites all 19 entries of each
   population family at that site. This removes dependence on the random
   mixture initialization at the tested site.
5. **Execute production code:** `reconstruct()` calls `calc_moments()` and
   returns copies of `density_lb`, `phi_lb`, `psi_lb`, and `u_lb`. It does
   not recalculate the expected answer.
6. **Compare:** `ReconstructsDensityCompositionAndVelocity` uses six
   `EXPECT_NEAR` checks with absolute tolerance $10^{-13}$.
7. **Clean up:** `TearDown()` resets the owning `std::unique_ptr<LAMMPS>`.
   LAMMPS destroys its fix and associated resources.

No timesteps or halo exchanges are required: `calc_moments()` reads only
populations at the selected site. The shared test main initializes MPI;
the test runs with one MPI rank.

### Coverage and limitations

This test checks local population sums, axial momentum signs, field assignment,
and division by density. It does not check nonzero diagonal contributions,
distributed `g`/`k` populations, pressure correctness, collision, streaming,
global conservation, or MPI communication. Although `calc_moments()` also
updates pressure, no pressure assertion is currently included.

During initial validation, removing division by density from the x-velocity
calculation in a temporary source copy caused this test to fail: the actual
value became 0.3 instead of 0.075. The temporary change was restored, and the
reconstruction and existing conservation tests passed. This documents that
validation run, not a guarantee about subsequent source edits.

## 2. Proposed: gradient and Laplacian

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

The proposed unit tests use $h=1$. Conversion of field units is a separate
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

### Proposed test cases

Use local coordinates $(X,Y,Z)$ relative to the chosen site, so the tested
site is $(0,0,0)$. Evaluate each polynomial over its complete stencil neighborhood.

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
failures easier to diagnose.

### Implementation plan

1. Reuse the small one-rank fixture and friend-access mechanism. Keep
   `calc_gradient_laplacian()` private.
2. Extend the test helper to fill a selected scalar field in the complete
   $3\times3\times3$ neighborhood around an owned interior site. D3Q19 reads
   only the center, axial neighbors, and face diagonals; initializing the
   whole neighborhood keeps setup simple.
3. Add a helper that calls the production derivative function and returns
   copies of its three gradient components and one Laplacian value.
4. Parameterize the analytic cases above. Expected values must come from
   polynomial derivatives, not a second implementation of the stencil sum.
5. Start with absolute tolerance $10^{-13}$ for these small local values.
   Reassess roundoff bounds if field magnitudes are increased.
6. Call the function twice at the same site with different fields, ending
   with a constant field. Confirm outputs are reset rather than accumulated.
7. Exercise the `rho`, `phi`, and `psi` field/output pairs. Checking the
   automatic wiring through `calc_equilibrium()` can be a separate integration
   test; direct stencil tests alone do not validate that wiring.
8. Register a separate one-rank CTest entry, proposed name
   `FixLBMulticomponent:derivatives`, with a matching GoogleTest filter.

The local polynomials need not be periodic over the whole box because this
test neither exchanges halos nor advances the fluid. It requires a valid
neighborhood around the selected site. Boundary and MPI halo behavior require
separate tests.

### Acceptance criteria

- Every analytic case passes for all three gradient components and the Laplacian.
- Repeated calls overwrite previous results.
- A deliberate sign or normalization error makes an appropriate case fail.
- Existing reconstruction and conservation tests still pass.
- No production numerical behavior is changed merely to accommodate the tests.

## 3. Build and run the implemented tests

From a clean CMake-compatible source tree with an MPI toolchain:

```sh
cmake -S cmake -B build-latboltz -D BUILD_MPI=ON -D PKG_LATBOLTZ=ON -D ENABLE_TESTING=ON
cmake --build build-latboltz --target test_fix_lb_multicomponent -j 4
ctest --test-dir build-latboltz -R '^FixLBMulticomponent:moments$' --output-on-failure
```

Run both currently registered LATBOLTZ tests with:

```sh
ctest --test-dir build-latboltz -R '^FixLBMulticomponent:' --output-on-failure
```

`FixLBMulticomponent:moments` uses one MPI rank; the existing
`FixLBMulticomponent:plain` conservation test uses two. GoogleTest filters
select the intended test from their shared executable. The executable still
receives the existing YAML file because the shared force-style test main
expects it; the moment test uses its C++ tolerance, not the YAML epsilon.

The derivative test name above is proposed and will not run until implemented
and registered in CMake.
