# Thesis: The *Dendronotus* CPG as a Controllable Oscillator

## Core Statement

**The *Dendronotus iris* swim CPG achieves what most small circuits cannot: smooth, single-knob frequency control that doesn't break under perturbation. The mechanism is an E/I module poised near oscillation, with synaptic kinetics matched to the rhythm's timescale. Presynaptic modulation of release rate—not postsynaptic conductance—provides the control lever.**

---

## Supporting Structure

### 1. The E/I module is the oscillatory engine

The contralateral Si3→Si2 loop with delayed inhibitory feedback sets the rhythm. The half-center oscillators (Si2–Si2, Si3–Si3) stabilize and structure alternation; they don't drive it.

### 2. The engine is poised, not locked

Each E/I module sits near the oscillatory threshold—close enough to be recruited by coupling, far enough to be tunable. This is the operating regime where control is possible without fragility.

### 3. Timescale compatibility is the design constraint

The module works because its slow synaptic dynamics evolve on the same timescale as the burst period. Fast synapses cannot do this job. The rhythm lives in the slow variables.

### 4. Neuromodulation enables control, not existence

The rhythm can run without modulation. What Si1 provides is *range*: the ability to shift frequency smoothly across behavioral demands. Gating and frequency control are separable functions.

### 5. Presynaptic kinetics beat postsynaptic conductance

Modulating release rate (α) produces graded, monotonic frequency changes. Modulating conductance (g) produces threshold-like jumps and narrower operating ranges. If you want smooth control, tune kinetics.

---

## Context (not thesis, but motivation)

Single-cell bursters are robust but hard to tune without side effects.  
Two-cell HCOs alternate but tend toward brittle control modes.  
This architecture—near-oscillatory E/I modules with slow synaptic variables—offers both robustness and controllability from a single biological input (Si1 firing rate).

That is the theoretical contribution: a plausible design principle for controllable CPGs.

---

## Open Questions

1. **Operational definition of "near-oscillatory"**: Damped oscillations? Subthreshold fluctuations? What simulation shows this cleanly for each isolated module?

2. **Timescale compatibility as falsifiable claim**: Can we show failure when synaptic kinetics are too fast or too slow? A parameter sweep would make this concrete.

3. **Quantifying "smoother"**: What metric—monotonicity, dynamic range, failure rate—best captures the presynaptic advantage? Even a limited comparison is valuable.

4. **Essential validations**: Which perturbation results must stay in main text to support "robust to perturbation"? Which can move to Supplement?

5. **Si1→frequency relationship**: What's the cleanest figure showing the model reproduces the biological Si1 firing rate → burst frequency curve?
