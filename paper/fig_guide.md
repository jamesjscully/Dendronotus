# Figure Guide: Dendronotus CPG Paper

This document maps paper figures to their source MATLAB files.

**Directory Key:**
- `jack_figs/` = Primary working directory for figure generation
- `m_files/` = Older/reference MATLAB files
- `old figs/` = Archived versions

---

## Original 8-Figure Structure (manuscript_6.tex)

| Fig # | File | Label | Content |
|-------|------|-------|---------|
| 1 | fig1.jpg | fig:bodysize_collage | Summary/Circuit architecture |
| 2 | fig2.jpg | fig:antiphase_pattern | E/I mechanism |
| 3 | fig3.jpg | fig:fig1 | CPG Assembly (5 steps) |
| 4 | fig4.jpg | fig:perturbation_responses | Perturbations - Validation I (6 panels, S&K 2016 Fig 5) |
| 5 | fig5.jpg | fig:group_fig7 | Symmetric perturbations - Validation II (4 panels, S&K 2016 Fig 7) |
| 6 | fig6.jpg | fig:group_fig5 | Synaptic timescales |
| 7 | fig7.jpg | fig:si1_control | Si1 control (4 panels) |
| 8 | fig8.jpg | fig:neuromodulation_comparison | Neuromodulation comparison (2 panels) |

### Current 6-Figure Structure (figures.tex)
- Fig 5 (symmetric perturbations) → moved to supplementary
- Fig 7 + Fig 8 → combined into single Figure 5

---

## Figure 1: Summary / Graphical Abstract
**File:** `fig1.jpg`
**Label:** `fig:bodysize_collage`

| Panel | Description | MATLAB Source | Output |
|-------|-------------|---------------|--------|
| A | Circuit schematic with Si1 | *Manual illustration* | — |
| B1 | Fast burst frequency (~2s period) | `SIN41_32_assembly_line_als_reg_burst_2sec.m` | `norm_swim_2sec.jpg` |
| B2 | Medium burst frequency (~5s period) | `SIN41_32_assembly_line_als_reg_burst_5sec.m` | `norm_swim_5sec.jpg` |
| B3 | Slow burst frequency (~12s period) | `SIN41_32_assembly_line_als_reg_burst_12sec.m` | `norm_swim_1o_12sec.jpg` |
| C1-C3 | Body size conceptual illustration | *Manual illustration (nudibranch)* | — |

**Caption:** Circuit architecture and burst frequency adaptability in the *Dendronotus iris* swim CPG.

---

## Figure 2: E/I Module Mechanism
**File:** `fig2.jpg`
**Label:** `fig:antiphase_pattern`

| Panel | Description | MATLAB Source | Output |
|-------|-------------|---------------|--------|
| A1-A2 | Circuit schematics showing alternating active modules | *Manual illustration* | — |
| B | Voltage traces (2L, 3R, 2R, 3L) + synaptic release variables | `SIN41_32_coupled_Jack.m` | `normal_swim.jpg` |
| C1-C2 | Whole-body undulation sketch | *Manual illustration (nudibranch)* | — |

**Caption:** Mechanism of left--right alternation: two reciprocally inhibitory E/I sub-networks generate anti-phase bursting.

---

## Figure 3: CPG Assembly
**File:** `fig3.jpg`
**Label:** `fig:fig1` / `fig:assembly`

| Panel | Description | MATLAB Source | Output |
|-------|-------------|---------------|--------|
| A | Uncoupled neurons (Si3 tonic, Si2 quiescent) | `SIN41_32_coupled_Jack_assembly_line_ALS_step1.m` | `assembly_step1.jpg` |
| B | Si3→Si2 fast excitation added | `SIN41_32_coupled_Jack_assembly_line_ALS_step2.m` | `assembly_step2.jpg` |
| C | Si2→Si3 slow inhibition added | `SIN41_32_coupled_Jack_assembly_line_ALS_step3.m` | `assembly_step3.jpg` |
| D | Si3↔Si3 reciprocal inhibition added (rhythm emerges) | `SIN41_32_coupled_Jack_assembly_line_ALS_step4.m` | `assembly_step4.jpg` |
| E | Full network (Si2↔Si2 inhibition + electrical coupling) | `SIN41_32_coupled_Jack_assembly_line_ALS_step5.m` | `assembly_step5.jpg` |

**Caption:** CPG assembly in five sequential steps.

---

## Figure 4: Targeted Perturbations (Validation I)
**File:** `fig4.jpg`
**Label:** `fig:perturbation_responses`
**Reference:** Sakurai & Katz (2016) Fig 5A-C, E-G

| Panel | Description | S&K Panel | MATLAB Source | Output |
|-------|-------------|-----------|---------------|--------|
| A | Depolarize Si3L | 5A | `SIN41_32_coupled_Fig5A_als.m` | `5a_als.jpg` |
| B | Hyperpolarize Si2R | 5B | `SIN41_32_coupled_Fig5B_als.m` | `5b_als.jpg` |
| C | Bilateral Si2 hyperpolarization | 5C | `SIN41_32_coupled_Fig5C_als.m` | `5c_als.jpg` |
| D | Depolarize Si3R | 5G | `SIN41_32_coupled_Fig5G_als.m` | `5g_als.jpg` |
| E | Hyperpolarize Si3R | 5E | `SIN41_32_coupled_Fig5E_als.m` | `5e_als.jpg` |
| F | Bilateral Si3 hyperpolarization | 5F | `old figs/SIN41_32_coupled_Fig5F.m` | `5f_als1.jpg` |

**Note:** All perturbation files use similar base parameters with only the applied current (Iapp) and timing (t1, t2) differing.

**Caption:** Validation of model resilience through targeted perturbations.

---

## Figure 5: Symmetric Perturbations (Validation II)
**File:** `fig5.jpg`
**Label:** `fig:group_fig7`
**Reference:** Sakurai & Katz (2016) Fig 7Ai-Di

| Panel | Description | S&K Panel | MATLAB Source | Output |
|-------|-------------|-----------|---------------|--------|
| A | Bilateral Si2 depolarization (schematic + traces) | 7A | `SIN41_32_coupled_Fig7A_als.m` | `7a_als1.jpg` |
| B | Bilateral Si2 depolarization (continued traces) | 7B | `SIN41_32_coupled_Fig7B_als.m` | `7b_als.jpg` |
| C | Bilateral Si3 depolarization (schematic + traces) | 7C | `SIN41_32_coupled_Fig7C_als.m` | `7c_als.jpg` |
| D | Bilateral Si3 depolarization (continued traces) | 7D | `SIN41_32_coupled_Fig7D_als.m` | `7d_als1.jpg` |

**Caption:** Symmetric perturbation experiments validating the swim CPG model (Validation list II).

---

## Figure 6: Synaptic Timescales
**File:** `fig6.jpg`
**Label:** `fig:group_fig5` / `fig:synaptic_timescales`
**Reference:** Sakurai & Katz (2022) Fig 7

| Panel | Description | MATLAB Source | Output |
|-------|-------------|---------------|--------|
| A | Circuit schematic with current injection electrode to Si2R | *Manual illustration* | — |
| B | 2-cell response: Si2L, Si2R with multiple synaptic rate traces | `ZOO_synapses_Yosef_1.m` | `ZOO_synapses.jpg` |
| C | 4-cell response: 3 graded current pulses to Si2R | `ZOO_synapses_Yosef_1.m` (different pulse protocol) | — |
| D | 4-cell response: 2 current pulses showing dual-timescale dynamics | `ZOO_synapses_Yosef_1.m` (different pulse protocol) | — |

**Note:** `ZOO_synapses_Yosef_1.m` runs a 2-cell model with current pulses at t1=10s, t2=15s (small), t5=25s, t6=31s (medium), t3=42s, t4=50s (large). The figure panels appear to extract different time windows.

**Caption:** Validation of slow synaptic dynamics and electrical coupling mechanisms.

---

## Figure 7: Si1 Control
**File:** `fig7.jpg`
**Label:** `fig:si1_control`
**Reference:** Sakurai & Katz (2019) Figs 2A, 3A, 3B

| Panel | Description | S&K Panel | MATLAB Source | Output |
|-------|-------------|-----------|---------------|--------|
| A | Circuit schematic with Si1 neurons | Fig 2A | *Manual illustration* | — |
| B | Gradual Si1 ramp up → rhythm initiation → ramp down | Fig 3A | `S1_drives_CPG_Fig3A_als.m` | `s1_Fig3A_drives_cpg1.jpg` |
| C | Si1 hyperpolarization during ongoing rhythm | Fig 3B | `S1_drives_CPG_NOneuroM_fig3B.m` | `s1_drives_CPG_Fig3B_No_NM_als.jpg` |
| D | Sustained Si1 activity maintaining rhythm | — | `S1_drives_CPG_neuroM_fig2A_als.m` | `s1_Fig2A_drives_cpg_als.jpg` |

**Caption:** Role of Si1–Si3 interactions in rhythmogenesis, episode control, and network resilience.

---

## Figure 8: Neuromodulation Comparison
**File:** `fig8.jpg`
**Label:** `fig:neuromodulation_comparison`

| Panel | Description | MATLAB Source | Output |
|-------|-------------|---------------|--------|
| A | Presynaptic modulation (α parameter): smooth, graded frequency control | `si1_control_neuroMod_a.m` | — |
| B | Postsynaptic modulation (conductance g): threshold-like transitions | `si1_control_neuroMod_g.m` | — |

**Additional analysis:**
- `alpha_neuromod_scatter.m` → burst frequency vs α scatter analysis
- `compare_three_models.m` → systematic comparison of modulation strategies

**Caption:** Neuromodulation effects on burst frequency control.

---

## Parameter Reference

### Base 4-cell network parameters (typical values from perturbation files):
```matlab
% Electrical coupling
gelec = 0.002-0.005;

% Si3→Si2 excitation (fast α-synapse)
alpha4 = 0.05-0.06;  beta4 = 0.001-0.008;
g41 = 0.005-0.008;   g32 = g41;

% Si2→Si3 inhibition (slow logistic synapse)
alpha1 = 0.012-0.015;  beta1 = 0.001;
g14 = 0.01;            g23 = g14;

% Si3↔Si3 reciprocal inhibition
alpha34 = 0.005-0.01;  beta34 = 0.005;
g34 = 0.01-0.015;      g43 = g34;

% Si2↔Si2 reciprocal inhibition
alpha12 = 0.005;  beta12 = 0.005;
g21 = 0.009-0.01; g12 = g21;

% Calcium shifts (determines intrinsic excitability)
Ca_shift4 = -90 (Si3: high excitability)
Ca_shift1 = -20 to -34 (Si2: near threshold)
```

---

## File Naming Conventions

- `Fig5X_als.m` → Perturbation experiments (S&K 2016 Fig 5)
- `Fig7X_als.m` → Symmetric perturbations (S&K 2016 Fig 7)
- `S1_drives_CPG_*.m` → Si1 control experiments
- `si1_control_neuroMod_*.m` → Neuromodulation experiments
- `assembly_line_*.m` → Network assembly steps
- `*_reg_burst_Xsec.m` → Different burst period variants

---

## Remaining Questions

1. **Figure 6 panels C-D:** The `ZOO_synapses_Yosef_1.m` file runs a 2-cell simulation, but panels C-D show 4-cell responses. Is there a separate file for the 4-cell current injection experiments, or are these panels composited from different runs?

2. **Figure consolidation for Current Biology:** Should we maintain the 8-figure structure or use the consolidated 6-figure version?
