# Current Biology Reviewer Report

**Manuscript:** "Slow synaptic dynamics and neuromodulatory control generate robust half-center rhythms from non-bursting neurons"

---

## General Assessment

This is a well-written computational study of the *Dendronotus iris* swim CPG that makes interesting mechanistic claims about how slow synaptic dynamics enable rhythm generation in circuits lacking intrinsic pacemaker neurons. The central finding—that kinetic (presynaptic) modulation provides smoother frequency control than conductance (postsynaptic) modulation—is potentially significant for understanding CPG flexibility. However, I have concerns about the validation depth, the novelty claims, and aspects of the presentation that require attention before this work is suitable for Current Biology's broad readership.

---

## Major Concerns

### 1. Validation is qualitative, not quantitative

The authors acknowledge this limitation (Discussion line 45-46), but it significantly weakens the paper's claims. The perturbation experiments (Fig. 4) reproduce *qualitative* response patterns from Sakurai & Katz (2016), but no statistical comparison is provided. Key questions remain:

- What are the burst periods, duty cycles, and phase relationships in the model vs. experiment?
- How robust is the model to parameter variation? The authors claim "robustness" but provide no quantitative metric (coefficient of variation, bifurcation margins, parameter volumes).
- Without fitting or systematic comparison, how can readers assess whether the model is predictive vs. merely descriptive?

**Recommendation:** At minimum, provide a table comparing key quantitative measures (burst duration, period, duty cycle, phase lag) between model and experimental values where available.

### 2. The core mechanistic claim requires experimental validation

The central claim—that kinetic modulation provides smooth frequency control while conductance modulation is threshold-like—is derived entirely from model behavior. This is a testable prediction, but no experimental evidence is provided or even proposed in concrete terms.

- Has any attempt been made to determine which modulatory mechanism operates in *Dendronotus*?
- What specific experiments would distinguish these mechanisms?
- Can the authors point to any indirect evidence favoring one mechanism?

For Current Biology, computational predictions ideally would be accompanied by at least preliminary experimental support, or the predictions should be presented more cautiously.

### 3. Novelty relative to existing CPG modeling literature

The concept that slow synaptic dynamics contribute to rhythm generation is well-established (Manor et al. 1997, Nadim & Bucher 2014—both cited). The E/I oscillator framework dates to Wilson & Cowan (1972). The parallel to *Tritonia* facilitation mechanisms is extensively noted.

**What is genuinely new here?** The authors should more clearly articulate:
- How does "network hysteresis" differ from established concepts of network-level bistability in CPGs?
- What predictions does this model make that differ from existing CPG models?
- The claim that synaptic kinetics are "co-equal partners" with cellular properties (Discussion line 58)—is this distinct from the existing literature?

### 4. Missing analysis of parameter sensitivity

The model contains numerous parameters (synaptic conductances, time constants, modulation gains). How sensitive is the core result—smooth kinetic control vs. abrupt conductance control—to these choices?

- Is the qualitative difference maintained across parameter space?
- Could different parameter choices produce the opposite result?

Without this analysis, the generality of the findings is uncertain.

---

## Minor Concerns

### 5. Figure quality and labeling
- Fig. 1 (assembly) would benefit from consistent color coding across panels
- Several figures reference comparison with experimental data (Figs. 2, 4, 5), but side-by-side experimental panels are not shown—this makes validation difficult for readers unfamiliar with the original papers

### 6. Title accuracy
The title states "half-center rhythms" but the main mechanistic claim is that E/I modules, not half-centers, are the primary oscillatory engine (Results, line 17: "stabilizers—the E/I modules supply the oscillatory drive"). Consider revising for consistency.

### 7. Missing detail in Methods
- The STAR Methods mention time constant τ_x is "100 or 235ms"—when is each value used and why?
- Synaptic conductance values (g_syn) are not provided in the key resources table or methods
- The GitHub repository placeholder needs to be completed before publication

### 8. Si1 representation
The model uses a single Si1 neuron while the biological system has bilateral Si1L/R neurons (acknowledged in limitations). Could asymmetric Si1 activity produce interesting dynamics relevant to behavioral flexibility? This seems worth brief discussion.

### 9. Introduction could be tightened
The introduction effectively motivates the problem but is somewhat dense. The four bullet points (lines 34-40) read like a methods overview rather than highlighting the conceptual advances. Consider restructuring to emphasize the key biological insights.

### 10. "Network hysteresis" definition
This term is introduced (Introduction line 10, Results line 6) but never formally defined. What distinguishes network hysteresis from cellular hysteresis mechanistically? A clear definition early in the paper would help readers.

### 11. Broader implications understated
The Discussion mentions developmental scaling and evolutionary adaptation briefly but doesn't develop these ideas. Given Current Biology's broad readership, more attention to how these findings might generalize beyond nudibranch CPGs would strengthen the paper's impact.

---

## Additional Questions for Authors

1. The model predicts that blocking/slowing the Si2→Si3 logistic synapses should eliminate bursting while preserving tonic activity. Has this been tested experimentally?

2. What happens to the model if Si3 neurons are given weak intrinsic bursting capability? Does the smooth kinetic control break down?

3. The authors note convergence with *Tritonia* mechanisms. Could the model architecture be applied to *Tritonia* data to test its generality?

---

## Summary

This is a competent computational study that develops a plausible mechanistic account of rhythm generation in a well-characterized CPG. The writing is clear and the model is systematically developed. However, the lack of quantitative validation, absence of experimental support for the key kinetic vs. conductance prediction, and incremental novelty relative to existing CPG theory are significant limitations.

**Recommendation:** Major revision. The paper would be strengthened by (1) quantitative model-experiment comparisons where data permit, (2) parameter sensitivity analysis for the key finding, and (3) either preliminary experimental tests or clearer articulation of testable predictions with specific proposed experiments.
