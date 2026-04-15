# Novelty Assessment Report: Manuscript 7 - Dendronotus CPG Model

## Executive Summary

This report assesses the novelty of the major claims in manuscript_7 against the historical CPG literature (1970s-2000s and beyond). The assessment identifies which claims represent genuine advances versus restatements of established knowledge.

---

## Major Claims Analyzed

### 1. "Network Hysteresis" as a Novel Mechanism

**Manuscript Claim:** The paper introduces "network hysteresis" as a mechanism whereby slow synaptic dynamics substitute for cellular bistability to generate rhythm in non-bursting neurons.

**Literature Assessment:**

| Aspect | Status |
|--------|--------|
| Terminology "network hysteresis" | **GENUINELY NOVEL** |
| Underlying concept | **PARTIALLY KNOWN** |

**Evidence:**
- The term "network hysteresis" appears to be **newly coined by your lab** in the 2024 Frontiers paper (Scully et al., 2024). Web searches return only your group's publications using this specific terminology.
- However, the underlying phenomena are well-established:
  - **Synaptic depression mediating bistability:** Bhalla & Bhargava (2001) showed that "synaptic depression mediates bistability in neuronal networks with recurrent inhibitory connectivity" and that "hysteresis and bistability occurred...when both synapses were depressing."
  - **Wilson-Cowan (1972):** Already described "the existence of multiple stable states, and hysteresis, in the population response" for E/I networks.
  - **Marder & Calabrese (1996):** Reviewed extensively how "synaptic and cellular processes interact to play specific roles in shaping motor patterns."

**Recommendation:** The *terminology* is novel, but the manuscript currently implies the *mechanism* itself is new. Consider revising to: "We term this phenomenon *network hysteresis* to distinguish it from cellular hysteresis mechanisms..." This frames it as naming/clarifying an existing concept rather than discovering it.

---

### 2. Network Oscillators Without Intrinsic Bursting

**Manuscript Claim:** Circuits lacking intrinsic pacemakers can generate robust rhythms through network interactions alone.

**Literature Assessment:**

| Aspect | Status |
|--------|--------|
| Basic concept | **WELL-ESTABLISHED (1980s - Getting's work)** |
| Application to Dendronotus | **INCREMENTAL** |

**Evidence:**

**CRITICAL: Peter Getting's Tritonia work (1983, 1989) already established this:**
- **Getting (1989)** paper was literally titled **"A network oscillator underlying swimming in Tritonia"**
- Getting explicitly showed that **"none of the CPG neurons has been observed to have intrinsic bursting properties. Rather, the rhythmic bursting during the swim motor pattern seems to arise through conventional synaptic interactions, making the CPG a rare example of a 'Network Oscillator'"**
- This is the **same organism family** (nudibranchs) and **same behavior** (swimming) as Dendronotus

- **Getting (1983):** "Mechanisms of pattern generation underlying swimming in Tritonia. III. Intrinsic and synaptic mechanisms for delayed excitation" - described slow synaptic mechanisms
- **Skinner, Kopell & Marder (1994):** Showed "four different mechanisms that lead to oscillations in a network of two reciprocally inhibitory cells"
- **Wang & Rinzel (1992):** Analyzed "alternating and synchronous rhythms in reciprocally inhibitory model neurons"

**The Problem:** The manuscript implies this is a gap ("Classical CPG models depend on intrinsically bursting neurons...but many circuits lack these features") when Getting already characterized the Tritonia swim CPG as a network oscillator without intrinsic bursters **40 years ago**.

**Recommendation:** Significantly revise the framing. The introduction should acknowledge: "Network oscillators lacking intrinsic bursters have been characterized, most notably in the Tritonia swim CPG (Getting 1983, 1989). However, the mechanisms enabling *smooth frequency control* in such circuits remain less understood..."

---

### 3. E/I Modules as Rhythmogenic Building Blocks

**Manuscript Claim:** Excitatory-inhibitory modules operating in an "almost-oscillatory" regime constitute the core oscillatory engine.

**Literature Assessment:**

| Aspect | Status |
|--------|--------|
| E/I oscillations concept | **ESTABLISHED (Wilson-Cowan 1972)** |
| "Almost-oscillatory" framing | **NOVEL TERMINOLOGY** |
| Application to CPG timescales | **PARTIALLY NOVEL** |

**Evidence:**
- **Wilson & Cowan (1972):** Already showed that E/I interactions "presented a simple mechanism for oscillations in firing rate activity." Their model is the foundational work on E/I oscillations.
- **However:** Wilson-Cowan models traditionally operate on fast (gamma-range) timescales. The manuscript correctly notes (Introduction line 19): "standard E/I models employ fast synapses tuned for gamma rhythms, not CPG timescales."
- The term "almost-oscillatory" appears to be **novel terminology** for what bifurcation theorists call "subcritical" or "near Hopf" dynamics.

**Recommendation:** The novelty here is in *applying* E/I module concepts to slow CPG timescales and framing them as "almost-oscillatory." This is legitimate but should be framed more carefully. Consider: "While E/I oscillations are well-characterized at cortical timescales (Wilson & Cowan 1972), their role in slow CPG rhythms has been less explored. We show that..."

---

### 4. Presynaptic vs. Postsynaptic Modulation for Frequency Control

**Manuscript Claim:** Presynaptic modulation of release rate provides smooth frequency control, whereas postsynaptic conductance modulation yields threshold-like transitions.

**Literature Assessment:**

| Aspect | Status |
|--------|--------|
| Presynaptic/postsynaptic modulation distinction | **EXTENSIVELY STUDIED** |
| Smooth vs. threshold prediction | **PARTIALLY NOVEL** |
| Quantitative demonstration in this system | **NOVEL** |

**Evidence:**
- **Nadim & Bucher (2014):** Already reviewed that "neuromodulators frequently target presynaptic release machinery...in ways that can strongly reshape circuit dynamics." This is a paper you cite.
- **Katz, Getting & Frost (1994):** Demonstrated in Tritonia that "serotonergic dorsal swim interneurons heterosynaptically increase the amount of transmitter released" - showing presynaptic modulation controls CPG dynamics.
- **Sakurai & Katz (2003):** Showed "spike timing-dependent serotonergic neuromodulation of synaptic strength" - directly relevant presynaptic modulation work.
- **Stomatogastric literature (Marder lab):** Extensively documented both presynaptic and postsynaptic modulation mechanisms.

**However:** The specific prediction that presynaptic modulation produces *smooth* control while postsynaptic produces *threshold-like* transitions appears to be a **novel computational prediction** worth highlighting.

**Recommendation:** Frame as: "While both presynaptic and postsynaptic modulation are known to affect CPG dynamics (Nadim & Bucher 2014), the differential effects on frequency control curves have not been systematically compared. Our model predicts..."

---

### 5. Slow Synaptic Dynamics as Primary Frequency Regulators

**Manuscript Claim:** Slow synaptic kinetics (hundreds of milliseconds to seconds) can substitute for cellular bursting mechanisms and serve as primary regulators of rhythm.

**Literature Assessment:**

| Aspect | Status |
|--------|--------|
| Slow synaptic dynamics affect CPG rhythm | **WELL-ESTABLISHED** |
| Can substitute for cellular mechanisms | **ESTABLISHED (Getting 1983)** |
| Specific timescale analysis | **CONTRIBUTES TO LITERATURE** |

**Evidence:**

**CRITICAL: Getting already described this in Tritonia:**
- **Getting (1983):** "Intrinsic and synaptic mechanisms for **delayed excitation**" - described how slow synaptic mechanisms generate rhythm
- In Tritonia: "The conventional synaptic potentials evoked by a DSI onto C2 and onto the DFNs consist of a **fast excitatory postsynaptic potential (EPSP) followed by a prolonged slow EPSP**"
- C2 has a "**delayed excitatory effect**" on VSIs - this is the slow synaptic dynamics operating on behavioral timescales
- The DSIs use **multiaction synapses** with both fast and slow components

**Other supporting literature:**
- **Zucker & Regehr (2002):** "Repeated presynaptic firing can facilitate or depress transmission over tens to hundreds of milliseconds, timescales comparable to burst durations."
- **Manor et al. (1997):** Studied "temporal dynamics of convergent excitation in an oscillatory network" in the stomatogastric system.
- **Respiratory CPG models:** Studies have shown "a synaptic depression/facilitation mechanism is sufficient for neurons to generate network rhythms, without the need for intrinsically rhythmic neurons."

**Recommendation:** Getting's "delayed excitation" is conceptually similar to your "slow synaptic dynamics." The manuscript should explicitly connect to this lineage: "Building on Getting's characterization of delayed excitation in the Tritonia swim CPG (1983), we show that slow synaptic kinetics in the Dendronotus circuit..."

---

### 6. Distributed Robustness / No Single Essential Neuron

**Manuscript Claim:** The network tolerates perturbations because no single neuron is essential for rhythm generation.

**Literature Assessment:**

| Aspect | Status |
|--------|--------|
| CPG degeneracy concept | **WELL-ESTABLISHED (Marder)** |
| Distributed robustness concept | **ESTABLISHED** |
| Demonstration in Dendronotus | **NOVEL VALIDATION** |

**Evidence:**
- **Marder (2011):** "Variability, compensation, and modulation in neurons and circuits" - extensively reviewed degeneracy.
- **Marder & Calabrese (1996):** Established that "CPG circuits exhibit degeneracy: similar output despite varied synaptic and cellular properties."
- **Wagner (2005):** "Distributed robustness versus redundancy as causes of mutational robustness" - established the theoretical framework.

**Recommendation:** The concept is not novel, but demonstrating it computationally in Dendronotus and matching experimental perturbation data is valuable. Frame as: "Consistent with the principle of degeneracy in CPG circuits (Marder 2011), our model reproduces..."

---

### 7. Smooth Frequency Control from a Single Command Signal (Si1)

**Manuscript Claim:** Si1 neurons provide a unified control signal that enables smooth, graded frequency modulation of the CPG output.

**Literature Assessment:**

| Aspect | Status |
|--------|--------|
| Single command → graded frequency | **WELL-ESTABLISHED (1960s onward)** |
| Mechanism in network oscillators | **PARTIALLY NOVEL** |
| Dendronotus-specific demonstration | **NOVEL** |

**Evidence:**
- **MLR (Mesencephalic Locomotor Region) - foundational work:**
  - "Electrical stimulation of the MLR triggers forward locomotion in a graded fashion as a function of stimulation intensity" (Grillner et al.)
  - "Low-level stimulation leads to slow (low frequency) movements, and high-level stimulation to faster (higher frequency) movements"
  - This principle was established in the 1960s-1970s in cats, then confirmed across all vertebrates including lamprey.

- **Lamprey CPG:**
  - "The power of swimming...increased as the intensity or frequency of the stimulating current were increased"
  - "The frequency of oscillation can be modulated by the tonic input, with the frequency of oscillation increasing with the input level"

- **Invertebrate command neurons:**
  - Leech, Tritonia, Clione all have command-like neurons that activate swimming
  - Serotonergic modulation in Tritonia accelerates swim frequency
  - In Clione, "modulation of swimming speed involves changes at network and cellular levels"

- **Zebrafish:**
  - "The dorsal-ventral position of excitatory interneurons is correlated with the minimal swimming frequency at which the neuron is active, suggesting...smoothly graded shifts in locomotor speed"

**The Core Issue:** The concept that a single descending drive can smoothly grade CPG frequency is a **foundational principle of motor control**, not a novel discovery. This has been known since the MLR stimulation experiments of the 1960s.

**What IS potentially novel:**
1. The *mechanism* by which this works in a network oscillator lacking intrinsic bursters
2. The specific finding that *presynaptic* (not postsynaptic) modulation enables smooth control
3. The demonstration in this particular simple circuit

**Recommendation:** The manuscript should NOT present single-command frequency control as novel. Instead, frame the question as: "While graded frequency control from descending drive is well-established (Grillner), the *synaptic mechanisms* enabling this in network oscillators remain unclear. We show that presynaptic modulation of release kinetics provides the substrate for smooth frequency control in the Dendronotus CPG."

---

### 8. Progressive Network Assembly Analysis

**Manuscript Claim:** Sequentially activating synapses reveals which interactions are necessary/sufficient for rhythm generation.

**Literature Assessment:**

| Aspect | Status |
|--------|--------|
| Methodological approach | **STANDARD PRACTICE** |
| Results for Dendronotus | **NOVEL** |

**Evidence:**
- This is a standard computational neuroscience approach for circuit analysis.
- The specific results about which connections are necessary/sufficient in Dendronotus are novel.

**Recommendation:** Don't overclaim the method; emphasize the results.

---

## Summary Table

| Claim | Novelty Level | Recommendation |
|-------|---------------|----------------|
| "Network hysteresis" terminology | **HIGH** | Keep, but clarify it's naming a phenomenon |
| Network oscillators without bursters | **LOW** | Getting (1989) already characterized Tritonia this way |
| E/I modules at CPG timescales | **MEDIUM** | Frame as novel application of known principles |
| Presynaptic vs. postsynaptic frequency effects | **MEDIUM-HIGH** | Novel prediction, acknowledge prior modulation work |
| Slow synaptic dynamics as regulators | **LOW** | Getting's "delayed excitation" (1983) is similar |
| Distributed robustness | **LOW** | Acknowledge Marder's degeneracy concept |
| **Single command → graded frequency** | **LOW** | Well-known since 1960s MLR work |
| **Mechanism enabling smooth control** | **MEDIUM-HIGH** | Novel mechanistic insight |
| Dendronotus-specific findings | **MEDIUM** | Novel for this species, but parallels Tritonia |

---

## POTENTIALLY NOVEL: The Controllability Argument

### The Core Insight

The manuscript could make a stronger, more defensible novelty claim by emphasizing:

**"The absence of cellular hysteresis enables controllability"**

| Neuron Type | PIR Properties | Frequency Control |
|-------------|----------------|-------------------|
| Bistable/hysteretic | PIR has "set" properties tied to intrinsic currents | Frequency relatively locked by cellular properties |
| Non-hysteretic | PIR is frequency/duration sensitive | Frequency tunable via synaptic parameters |

### Building on Skinner, Kopell & Marder (1994)

They showed four oscillation mechanisms:
- **Intrinsic release/escape** → frequency INSENSITIVE to synaptic threshold
- **Synaptic release/escape** → frequency SENSITIVE to synaptic threshold

Your extension: Non-hysteretic neurons necessarily rely on synaptic mechanisms, which makes frequency inherently sensitive to synaptic parameters → **controllable**.

### The Blueprint Argument

The Dendronotus circuit demonstrates a **design principle** for building controllable network oscillators:

1. Use non-bursting neurons (no cellular hysteresis)
2. Their PIR becomes frequency/duration sensitive
3. Add slow, frequency-sensitive synapses
4. Combine into E/I and HCO building blocks (as in Getting)
5. Result: A system where a single command signal can smoothly tune frequency

### How to Frame This

**Instead of:** "We show that network oscillators can work without intrinsic bursters" (already known from Getting)

**Use:** "We demonstrate that the ABSENCE of cellular hysteresis is a design feature that enables controllability. Non-hysteretic neurons exhibit frequency/duration-sensitive PIR, which - when driven by slow frequency-sensitive synapses - creates E/I and HCO modules whose output frequency can be smoothly tuned by a single command signal. The Dendronotus circuit provides a minimal blueprint for this principle."

---

## CRITICAL: The Getting Precedent (1983-1989)

Peter Getting's work on the Tritonia swim CPG established many of the concepts your manuscript presents:

| Your Manuscript Claims | Getting Already Showed (Tritonia) |
|------------------------|----------------------------------|
| Network oscillator without intrinsic bursters | "None of the CPG neurons has intrinsic bursting properties...a rare example of a 'Network Oscillator'" (1989) |
| Slow synaptic dynamics generate rhythm | "Delayed excitation" and "prolonged slow EPSP" mechanisms (1983) |
| Building blocks concept | "Cellular, synaptic, and network building blocks" (1989) |
| E/I interactions drive rhythm | DSI excites C2, C2 has delayed excitation to VSI, VSI inhibits DSI/C2 |

**This matters because:**
1. Tritonia and Dendronotus are **both nudibranchs**
2. Both exhibit **swimming behavior**
3. Both use **network oscillators** without intrinsic bursting
4. Your co-author (Katz) has extensive publications on both systems

**The manuscript needs to explicitly acknowledge this lineage** rather than implying these are novel insights.

---

## Key Revisions Recommended

### Introduction Revisions

**Current (line 4-5):** "Classical CPG models depend on intrinsically bursting neurons or reciprocal inhibition between bistable cells, but many circuits lack these features."

**Suggested:** "While pacemaker-driven and bistable half-center oscillators have been extensively studied (Skinner et al. 1994; Wang & Rinzel 1992), the mechanisms by which circuits lacking these features achieve smooth frequency control remain less well characterized."

---

**Current (line 10):** "We introduce the concept of *network hysteresis* to explain how this circuit produces tunable rhythm without a cellular pacemaker."

**Suggested:** "We propose the term *network hysteresis* to describe how slow synaptic dynamics, rather than cellular bistability, create the history-dependence underlying rhythm generation in this circuit."

---

### Discussion Revisions

**Current (line 4-5):** "Previous work established that synaptic properties shape CPG output (Manor1997, NadimBucher2014), but these studies examined circuits containing intrinsically bursting neurons or focused on synaptic amplitude rather than temporal dynamics."

**This is good** - it appropriately acknowledges prior work while carving out the novel contribution.

---

### Summary Revisions

**Current:** "How do networks of non-bursting neurons generate robust, frequency-tunable rhythms?"

**Suggested:** "While network oscillators based on non-bursting neurons are well-established (Skinner et al. 1994), the mechanisms enabling smooth frequency control in such circuits remain poorly understood. How does the Dendronotus swim CPG achieve this?"

---

### Introduction Revision for Frequency Control

**Current (line 7):** "The central mechanistic question is how circuits vary rhythm frequency smoothly without degrading phase structure or stability."

**Issue:** This is fine as a question, but be careful not to imply that smooth frequency modulation from command signals is itself novel.

**Suggested addition after this sentence:** "While graded locomotor control from descending brainstem commands is well-established across vertebrates and invertebrates (Grillner; Katz & Frost 1994), the synaptic mechanisms enabling smooth frequency tuning in network oscillators lacking intrinsic pacemakers remain less understood."

---

## Genuinely Novel Contributions (Highlight These)

1. **The term "network hysteresis"** as a unifying concept for synaptic-based rhythm generation

2. **Controllability through lack of cellular hysteresis** (POTENTIALLY NOVEL FRAMING):
   - Skinner, Kopell & Marder (1994) showed that "intrinsic escape/release" mechanisms make frequency **insensitive** to synaptic parameters, while "synaptic escape/release" makes frequency **sensitive** to synaptic parameters
   - Your contribution: Neurons WITHOUT cellular hysteresis/bistability have frequency/duration-sensitive PIR (unlike bistable neurons with "set" intrinsic frequencies)
   - This sensitivity, combined with slow frequency-sensitive synapses, creates a system optimized for controllability
   - **This is a novel framing** - emphasizing that the ABSENCE of intrinsic bistability is a *design feature* enabling tunability, not just a description of the circuit

3. **Quantitative prediction** that presynaptic modulation yields smooth frequency control vs. postsynaptic yielding threshold-like transitions

4. **Specific circuit analysis** of Dendronotus showing E/I modules + Si3 reciprocal inhibition = minimal rhythm generator

5. **Model validation** against multiple experimental perturbations without parameter fitting

6. **Blueprint for controllability**: The Dendronotus circuit as a minimal example showing HOW to achieve tunable network oscillations - non-hysteretic neurons + slow frequency-sensitive synapses + E/I/HCO building blocks

7. **Evolutionary comparison** suggesting conserved facilitation mechanisms across Dendronotus and Tritonia despite independent evolution of swimming

---

## Sources Consulted

### Critical Getting References (1983-1989)
- **Getting PA (1983a)** Mechanisms of pattern generation underlying swimming in Tritonia. II. Network reconstruction. J Neurophysiol 49:1017-1035
- **Getting PA (1983b)** Mechanisms of pattern generation underlying swimming in Tritonia. III. Intrinsic and synaptic mechanisms for delayed excitation. J Neurophysiol 49:1036-1050
- **Getting PA (1989a)** Emerging principles governing the operation of neural networks. Ann Rev Neurosci 12:185-204
- **Getting PA (1989b)** A network oscillator underlying swimming in Tritonia. In: Neuronal and Cellular Oscillators (Jacklet JW, ed), pp 215-236. Marcel Dekker

### Other Key References
- Marder E, Calabrese RL (1996) Principles of rhythmic motor pattern generation. Physiol Rev
- Skinner FK, Kopell N, Marder E (1994) Mechanisms for oscillation and frequency control. J Comp Neurosci
- Wang XJ, Rinzel J (1992) Alternating and synchronous rhythms. Neural Computation
- Wilson HR, Cowan JD (1972) Excitatory and inhibitory interactions. Biophys J
- Katz PS, Getting PA, Frost WN (1994) Dynamic neuromodulation. Nature
- Nadim F, Bucher D (2014) Neuromodulation of synaptic transmission. Curr Opin Neurobiol
- Sakurai A, Katz PS (2016, 2019, 2022) Dendronotus swim CPG studies. J Neurophysiol
- Scully J et al. (2024) Pairing cellular and synaptic dynamics. Front Netw Physiol
- Grillner S et al. (multiple) Lamprey locomotor CPG and MLR control
- Shik ML, Severin FV, Orlovskii GN (1966) Control of walking and running. Biofizika
- Ijspeert AJ (2008) Central pattern generators for locomotion control. Neural Networks
