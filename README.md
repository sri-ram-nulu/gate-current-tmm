## Part 2: Calculating Quantum Transmission Coefficient T(E)
**Developer:** Sri Ram Nulu

The objective is to calculate the probability that an electron can tunnel from the channel/oxide interface to the oxide/gate interface using the **Transfer Matrix Method** to solve 1-D time-independent Schrodinger Equation, for a given electron energy E, effective masses (m*) of electron in 3 regions, and potential profile accross the channel, which acts as a barrier.
### 2.1 Mathematical Approach
#### Wave Number ($k_i$)
In a region $j$ with constant potential $V_i$ and effective mass $m_i^*$, the wave number is defined using:
```math
-\frac{\hbar^2}{2m_i^*} \frac{d^2\psi(x)}{dx^2} + V_i\psi(x) = E\psi(x) \implies \frac{d^2\psi}{dx^2} = -\frac{2m_i^*}{\hbar^2}(E - V_i)\psi \implies k_i = \frac{\sqrt{2m_i^*(E - V_i)}}{\hbar}
```
#### Wave Function
In region j, $\qquad \psi(x) = A e^{ik_ix} + B e^{-ik_ix}$
This can be expressed in a matrix vector form in terms of forward-propagating component ($\psi_{i}^+$; incident) and the backward-propagating component ($\psi_{i}^-$; reflected): 
```math
\mathbf{\Psi}(x) = \begin{pmatrix} \psi_{i}^+(x) \\ \psi_{i}^-(x) \end{pmatrix} = \begin{pmatrix} A_i e^{j k_i x} \\ B_i e^{-j k_i x} \end{pmatrix}
```

For $E>V_i$, the kinetic energy of electron is greater than the barrier potential energy, yielding an oscillatory propagating wave through the channel. For $E<V_i$, it leads to an exponential decay in the incident electron flux.
#### Space Propagation Matrix $M_0(k_i, \Delta x_i)$
Within a flat slice j of width $\Delta x_i$, both the ends of the slice can be related by Space Propagation Matrix as follows
```math
\begin{pmatrix} A_i e^{j k_i (x + \Delta x_i)} \\ B_i e^{-j k_i (x + \Delta x_i)} \end{pmatrix} = \begin{pmatrix} e^{j k_i \Delta x_i} & 0 \\ 0 & e^{-j k_i \Delta x_i} \end{pmatrix} \begin{pmatrix} A_i e^{j k_i x} \\ B_i e^{-j k_i x} \end{pmatrix} \implies M_0 = \begin{pmatrix} e^{j k_i \Delta x_i} & 0 \\ 0 & e^{-j k_i \Delta x_i} \end{pmatrix}
```
#### Interface Boundary Matrix $`M_s( \frac{k_i}{{m_i}^*} , \frac{k_{i+1}}{m_{i+1}^*} )`$
At the interface of slice i and slice i+1, the wavefunctions just before and after the boundary can be expressed from the 2 boundary conditions in Quantum Mechanics:
* Continuity $`\implies \psi_i(x) = \psi_{i+1}(x) \implies A_i e^{j k_i x} + B_i e^{-j k_i x} = A_{i+1} e^{j k_{i+1} x} + B_{i+1} e^{-j k_{i+1} x}`$ 
* Differentiability $`\implies \frac{1}{m_i^*}\frac{d\psi_i}{dx} = \frac{1}{m_{i+1}^*}\frac{d\psi_{i+1}}{dx} \implies \frac{j k_i}{m_i^*}\left( A_i e^{j k_i x} - B_i e^{-j k_i x} \right) = \frac{j k_{i+1}}{m_{i+1}^*}\left( A_{i+1} e^{j k_{i+1} x} - B_{i+1} e^{-j k_{i+1} x} \right)`$

The above 2 conditions can be integrated into a matrix after transformation as shown below:
```math
\begin{pmatrix} \psi_{i+1}^+(x) \\ \psi_{i+1}^-(x) \end{pmatrix} = \frac{1}{2} \begin{pmatrix} 1 + \mathcal{R} & 1 - \mathcal{R} \\ 1 - \mathcal{R} & 1 + \mathcal{R} \end{pmatrix} \begin{pmatrix} \psi_{i}^+(x) \\ \psi_{i}^-(x) \end{pmatrix} \implies M_s=\frac{1}{2} \begin{pmatrix} 1 + \mathcal{R} & 1 - \mathcal{R} \\ 1 - \mathcal{R} & 1 + \mathcal{R} \end{pmatrix}, \quad \text{where } \mathcal{R} = \frac{(\frac{k_i}{m_i^*})}{(\frac{k_{i+1}}{m_{i+1}^*})}
```


#### Cascading Matrices
The continuous channel barrier profile is chopped into $N$ tiny flat steps. The program multiplies new layer operations **on the left** of the running matrix variable (`M = mmul(New_Matrix, M)`), cascading from the left interface ($`0 \rightarrow 1`$) to the right ($`N \rightarrow N+1`$):
```math
M_{\text{total}} = M_{s(N \rightarrow N+1)} \cdot M_{0(N)} \cdots M_{0(1)} \cdot M_{s(0 \rightarrow 1)}
```

This maps the initial coefficients on the far-left directly onto the final coefficients on the far-right:
```math
\begin{pmatrix} A_{N+1} e^{j k_{N+1} x_{N}} \\ B_{N+1} e^{-j k_{N+1} x_{N}} \end{pmatrix} = \begin{pmatrix} M_{00} & M_{01} \\ M_{10} & M_{11} \end{pmatrix} \begin{pmatrix} A_{0} e^{j k_{0} x_{0}} \\ B_{0} e^{-j k_{0} x_{0}} \end{pmatrix}
```

### Extraction of the Transmission Coefficient

Assuming an electron is incident only from the channel, there is no backward-propagating wave reflecting from the gate ($`B_{N+1} = 0`$). Under this boundary constraint, the total cascaded matrix system reduces to:

```math
\begin{pmatrix} A_{N+1} e^{j k_{N+1} x_{N}} \\ 0 \end{pmatrix} = \begin{pmatrix} M_{00} & M_{01} \\ M_{10} & M_{11} \end{pmatrix} \begin{pmatrix} A_{0} e^{j k_{0} x_{0}} \\ B_{0} e^{-j k_{0} x_{0}} \end{pmatrix}
```

The above matrix can be unrolled into 2 equations:

1. $`A_{N+1} e^{j k_{N+1} x_{N}} = M_{00} \left(A_{0} e^{j k_{0} x_{0}}\right) + M_{01} \left(B_{0} e^{-j k_{0} x_{0}}\right)`$
2. $`0 = M_{10} \left(A_{0} e^{j k_{0} x_{0}}\right) + M_{11} \left(B_{0} e^{-j k_{0} x_{0}}\right)`$


**Reflection Coefficient ($r$) & Transmission Coefficient ($`t`$):**
```math 
r = \frac{B_{0} e^{-j k_{0} x_{0}}}{A_{0} e^{j k_{0} x_{0}}} = -\frac{M_{10}}{M_{11}} \quad t = \frac{A_{N+1} e^{j k_{N+1} x_{N}}}{A_{0} e^{j k_{0} x_{0}}} = \frac{\det(M)}{M_{11}}
```

#### Calculating Physical Transmission Probability (T(E))
Because the incident channel material differs from the receiving gate material, the electrons possesses different effective masses. Transmission Probability T can be obtained by finding the ratio of probabilistic current densities in the path.
```math 
T = \frac{J_{\text{transmitted}}}{J_{\text{incident}}} = |t|^2 \times \frac{\hbar k_{N+1} / m_{N+1}^*}{\hbar k_{0} / m_{0}^*}
```

### 2.2 Numerical Implementation & Software Architecture
The computation code is entirely written in C language with plotting for verification using python.
#### Custom Data Structures
`Complex` struct is created to store real and imaginary parts of a complex number in a single element/identifier. `cmat2` struct is used to store a 2x2 matrix of `Complex` elements.
#### Complex Arithmetic Handling 
Created `cadd`, `cmul`, `cdiv`, `cexpj`, `c_abs2`, `c_abs` functions to handle the operations on `Complex` elements. The custom data structures and complex functions are based on `double` data type only.
#### Program I/O and Data Processing
The program initializes by parsing `Eminmax.txt` and `input_parameters.txt`to fetch the required parameters: `barrier_width`, `m_left` (channel), `m_barrier` (oxide), `m_right` (gate), and `Emin`, `Emax`, `steps_for_E`, required for finding T(E) for a range of electron energies (Emin, Emax).
Following the parameter configuration, the program parses the `poisson_sio2.txt` containing the electrostatic potential profile accross the oxide, generated by the `Poisson Solver` (Part 1)
**Energy Sweep Execution Loop** computes the transmission coefficient T(E) for every discrete energy step E in the range (Emin, Emax) with step size as `steps_for_E`. The value E and the corresponding T(E) are then written into `transmission_data.txt`.

### 2.3 Verification and Validation

To verify the numerical stability and physical accuracy of the TMM solver, the program was benchmarked against three distinct base potential profiles across an energy sweep from $`0.1\text{ eV}$ to $4.0\text{ eV}`$. 
**The inputs and results (including plots) for each test case can be seen in the T_E_test_{i} folders.**

* **Test Case 1 (Symmetric Rectangular):** A uniform potential barrier ($`U = 1.25\text{ eV}`$) with identical boundary masses.
* **Test Case 2 (Asymmetric Step-Down):** A rough triangular approximation cascading downward ($`1.875 \rightarrow 1.25 \rightarrow 0.625\text{ eV}`$) with asymmetric contact masses ($`m_{\text{left}} \neq m_{\text{right}}`$).
* **Test Case 3 (Asymmetric Step-Up):** The exact spatial inversion of Test Case 2 ($`0.625 \rightarrow 1.25 \rightarrow 1.875\text{ eV}`$), where both the potentials and boundary masses are completely mirrored ($`m_{\text{left, Test 3}} = m_{\text{right, Test 2}}`$ and vice versa).

#### Key Physical Observations
1. **Transmission Probability T(E) is always $\le$ 1** 
2. **Regime Behaviours:** In test case 1, Transmission probability plot perfectly demonstrates evanescent tunneling regime behaviour for $`E < V`$ and propagation regime behaviour for $`E > V`$.
3. **Quantum Symmetry and Reciprocity:** In test case 2 & 3, despite the structural inversion of the system, the system satisfies quantum reciprocity ($`T_{\text{left}\rightarrow\text{right}} = T_{\text{right}\rightarrow\text{left}}`$). Therefore the TMM solver follows the "Reciprocity Theorem".
