1) Placement/role of the “scaling” figure and how to frame it

Jack: Two ways to address the figure issue: (A) reframe it as not a visual abstract and keep it in the Introduction, or (B) change placement/framing.

Paul: Disagrees with keeping it as an intro “visual abstract.” The scaling result is better as an end-of-paper speculative prediction: the model can explore how behavior changes as the animal scales; larger animals plausibly move more slowly. Present it explicitly as a model prediction/speculation.

Jack: Notes Dr. Scholnikov prefers keeping the figure because he likes it; Jack anticipates needing to “navigate” that preference.

Paul: The figure may be “nice” but is not representative of the core paper.

Jack: If retained earlier, could try to justify it as aligned with the paper’s motivation: adaptability to changing conditions (ability to change with changing conditions).

Paul: This connects to title/positioning: current framing feels like “running before walking.” Title/abstract should reflect what was actually shown.

2) What the paper’s real centerpiece is (and how it should be framed)

Paul: The core achievement: showing you can account for a CPG / half-center oscillator behavior using measured properties—contradicting the claim (attributed to Andre) that a half-center oscillator can’t exist without intrinsic bursters (or wouldn’t be self-sustaining).

Jack: Pushback/clarification: that claim is not true; a half-center oscillator is (in a sense) defined by non-bursting elements producing rhythm via circuit dynamics.

Paul: Correct: the mistaken belief was about self-sustaining behavior; the empirical system contradicts it.

3) Critical missing result: what happens when you isolate one “half”

Paul: The paper currently speculates that SI2/SI3 “work as a group” but doesn’t clearly show the key finding in this manuscript: what happens when you isolate one half/module.

Paul: You should add results showing each half as a network oscillator (or near oscillating / on the verge), and how coupling those halves yields the full rhythm.

Paul: This is “the whole point” of one of your central figures: compare with electrical coupling vs without electrical coupling and complete the story by showing behavior of each half on its own.

Jack: Agrees this should go into the Results.

4) Clarifying the “2–3 mixture” / EI module and why it matters

Jack: Cites prior long paper section (modeling portion) showing a half-center circuit and also a 2–3 isolated case.

Paul: The important case is the 2–3 mixture (EI group). They shouldn’t burst on their own under some assumptions; in your model they can under certain coupling strengths.

Paul: Key difference in the example: no electrical coupling in that configuration. (Coupling language is confusing—see below.)

Jack: Claims oscillations can arise via EI networks without electrical coupling; can mimic experimentally (e.g., dynamic clamp) by adjusting synaptic strength.

Paul: This is important specifically for the dendritonus (context) story and must be made explicit: the paper doesn’t currently articulate how the CPG emerges from these non-bursting units.

5) Main conceptual message to “bring out”: dynamics + near-burstiness

Paul: The manuscript should emphasize:

how a half-center oscillator emerges from non-bursting cells;

each “half” is close to bursting (or near oscillation);

compatibility of timescales (slow synaptic dynamics) is central—“speeds of the currents need to be compatible.”

Paul: The paper currently doesn’t “hit the nail”: needs to foreground that slow synaptic dynamics promote resilience/regular bursting and that the dynamics (not just wiring) matter.

6) Abstract framing critique (overclaiming variability/perturbations)

Paul: First abstract sentence claims “despite variability in components and external perturbations…” but you didn’t actually do the breadth that sentence implies; it reads like the paper is going in another direction.

Jack: You did include external perturbations in the model—based on perturbations from the 2016 paper and some from ~2020.

Paul: Okay—but the abstract should still reflect what is actually demonstrated and emphasize the key insights above.

7) Figure overhaul: current central figure is hard to parse and conceptually inconsistent

Jack: The figure is large and doesn’t accentuate the key point; proposes redoing it to make the modeling pipeline clear: build two-cell networks → reduce coupling until non-oscillatory (or near) → assemble full circuit.

Paul: Agrees: current presentation is confusing and needs restructuring.

7a) Synaptic-variable traces are conceptually confusing

Paul: Confused by synaptic-variable traces shown when no synapse exists (coupling set to zero). If coupling is zero, “postsynaptic” is the wrong description.

Akira: Better to describe as “release of transmitter” rather than postsynaptic effect.

Paul: If there is no synapse, it cannot be postsynaptic.

Jack: Clarifies intent: the trace under each voltage is meant to be the output-associated synaptic variable driven by that neuron, not an input; but acknowledges wording and mapping are unclear.

Paul: Additional inconsistency: SI3 makes different-sign synapses onto different targets; SI2 makes both fast and slow—this is not represented clearly as “outputs” in the current visualization. The figure does not show which synapse each variable corresponds to.

7b) Visual encoding problems: shading, color, labeling

Paul: Shading difference (light/dark) is too subtle; hard to see and interpret; mid-panel “zoom” doesn’t fix clarity.

Paul: Avoid calling connections “diagonal” (nothing is literally diagonal). Use module terminology.

Paul: Suggestion: explicitly label two modules as alpha and beta, shade modules distinctly, and refer to them consistently.

7c) Axis/scalebar and visual clutter

Paul: Voltage axis labeling is inconsistent across panels; time is shown but voltage scale isn’t clear. Better: use scale bars (e.g., 1 s and some mV) on each panel; label the values once.

Paul: 5 s time window is not useful; 1 s is easier for interpretation.

Paul: Remove bounding boxes; they add clutter without information.

7d) Ordering and bilateral connections

Paul: Reorder neuron classes: “1s before 2s/3s” (put SI1 at top) for readability and consistency with prior model figures.

Paul: Some drawn connections are wrong because biologically they are bilateral, not unilateral. Either note bilaterality or redraw.

Jack: Proposes merging left/right where appropriate (e.g., one node labeled “1L/R”) since modeling used one representative cell anyway.

Paul + Akira: Agree; alpha/beta module depiction exists in prior work and is clearer.

8) Experimental-manipulation figures: don’t force readers to hunt prior papers

Paul: When showing perturbations/manipulations, you’re asking the reader to go back to prior papers to see the experiments—that’s too much. Better:

reproduce the relevant published figure(s) beside the model results, with proper permission/citation; or

include modified versions that mimic them; and/or

move less-central manipulations to Supplementary.

Jack: Notes it may become a large figure set; considers splitting into multiple figures (one per experiment), but recognizes this could bloat the paper.

Paul: Use Supplementary strategically: pick the most important manipulations for the main text; others can be supplemental.

9) Journal/venue strategy (visibility over convenience)

Jack: Initial thought: Frontiers Network Science (mentions “free paper”/fee reasons from advisor).

Paul: Strongly disagrees choosing a venue for fee convenience; prioritize where the audience will actually read it. Suggests Journal of Neurophysiology (where prior work was published) and argues the modeling would fit.

Paul: If the message is strong and focused, could even shop higher (e.g., Journal of Neuroscience / Current Biology), but then tighten the manuscript and push extensive figures to Supplementary.

Paul: Warning: don’t make papers undigestible by trying to include everything; focus and clarity matter.

10) “Fact checks” terminology

Paul: “Fact checks” reads negative and informal.

Jack: That wording is advisor-driven; open to change.

Paul: Use “model validation” / “validation through perturbation.”

11) Neuromodulation section: clarify mechanism and communicate insight

Jack: Two neuromodulation approaches in the model:

Conductance-based modulation: changing electrical properties by modifying a conductance term (GM) scaling a current (interpretable as changing receptor density/postsynaptic strength).

Kinetic modulation: changing the rate/alpha of transmitter release (presynaptic mechanism).

Paul: Wants the exact equation/definition; after review:

Conductance route corresponds to a postsynaptic mechanism (e.g., receptor density).

Kinetic route corresponds to a presynaptic mechanism (e.g., calcium sensitivity/release machinery).

Calling them presynaptic vs postsynaptic will be immediately intelligible to neurophysiologists.

Paul: Believes the kinetic/presynaptic story is likely closer to biology; this distinction should be explicit in the paper.

11a) Robustness / “smooth control” claim needs careful support

Jack: Subjective modeling observation: kinetic modulation provides smoother control; conductance modulation works but can be “clunkier,” sometimes breaks, and has a narrower operating range. Unsure how to visualize convincingly; results not super strong.

Paul: Even if evidence is limited, it’s still a valuable model-derived insight; can be stated as something learned from the model.

Paul: Another key insight: neuromodulation can provide frequency control (not just turning on rhythm), via dynamic (not static) changes in synaptic strength.

Jack: Suggests adding a plot (e.g., scatter) linking SI1 firing frequency to burst frequency; notes a similar relationship appears in the 2022 paper and could be reproduced with model curves/examples.

Paul: Yes—show model examples alongside the biological figure to amplify meaning.

11b) Higher-level conceptual claim: circuits are “modulated into existence”

Paul: Broader framing: circuits aren’t “turned on like a light switch”; they must be modulated into existence, with presynaptic control of synapses as a key locus. This is a major conceptual contribution, especially emphasizing dynamics.

Jack: Agrees; timescale compatibility is foundational—fast synapses can’t simply be strength-tuned to create slow behavior; some slow process must drive the rhythm.

Paul: This conceptual framing also ties back to the scaling prediction: a single control point could allow rhythmic period to scale with animal size.

12) Agreed “important points” list for abstract + paper focus (as drafted verbally)

The oscillator halves/modules are near-bursty / near oscillation.

The model shows how a half-center oscillator emerges from non-bursting cells.

Slow synaptic dynamics / timescale compatibility is central to robustness and rhythm generation.

Presynaptic modulation of synaptic release kinetics enables robust control of periodicity/frequency.

The model implies a single point of control that could support scaling with growth/size (present as prediction/speculation).

Validation via perturbations exists, but should be framed as model validation and potentially moved partly to Supplementary if it distracts.