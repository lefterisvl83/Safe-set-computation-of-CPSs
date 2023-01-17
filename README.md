# Safe set computation of Cyber Physical Control Systems subject to multi-attack actions.

### Repository

The repository contains matlab scripts for computing the maximal safe set of a cyber physical control system subject to persistent and/or intermittent attack patterns. 


### Cyber Physical Control Systems

Cyber-physical control systems represent a broad spectrum of safety-critical applications, ranging from power generation and distribution networks to autonomous mobility and industrial processes. Due to their extent and intrinsic link to society, secure operation of such schemes is vital. Vulnerability to cyber attacks typically depends on the degree of integrating unsafe communication channels between computation, sensing, and actuation modules that control the underlying physical process.

### Overall System description

We model the overall CPS under attack as a constrained switching system with the switching signal forming a regular language, generated by a nondeterministic directed graph. Each node of the graph is associated with a set of states that evolve with time according to the modes assigned to the corresponding outgoing edges. Each labelled edge describes either an attack-free operation or a specific malicious action carried out over a subset of unsafe channels. This approach to attack modelling permits rigorous description a large class of non-deterministic attack patterns. 

### Attack patterns

We consider attack policies that enable attackers to embed logic. Individual attack operations are made up of two main ingredients: the targeted channel(s) and the set of logic rules (e.g., dwell-time, attack channels). Logic rules are expressed via a regular language. The overall attack policy is described by a directed labelled graph. An edge indicates a set of attack operations, each acting on a specific system signal (e.g., measurements readings, actuation) over an unsafe channel. An edge also signifies the transition of the physical process in a single time step under the set of underlying attack actions, and, thus, is associated with a specific dynamic mode.

### Impact metrics

#### Safe set

By modelling the overall attack scheme as a constrained switching system, we characterise the set of all initial states that cannot be driven to an unsafe state under any allowable attack. We call this the _safe set_ of the attacked CPS. This is an infinite-reachability, dynamic programming problem: The maximal safe set can be retrieved by computing, in a recursive fashion, the fixed point of the sequence of sets $\{S_i\}_{i \in \{1,2,\ldots\}}$ with
$S_{i+1} = \textnormal{Pre}(S_i) \cap S_0,$ where $S_0 = X_0$ denotes the state-constraints set, and $\textnormal{Pre}(S_i)$ is the _preimage map_, that is the set of states $x$ for which, for all permissible attack patterns, the successor state $x^+\in S_i$. 

#### Minkowski and Lebesgue measure

### Reference

Detailed description of the overall system model and the attack patterns considered can be found in https://doi.org/10.48550/arXiv.2211.12196.

### Requirements

The computation of safe sets requires the MPT3 toolbox https://www.mpt3.org/ for the generation of Polyhedra objects.

### Acknowledgements

This work is the result of collaboration between Queen’s University Belfast, UK, University of Cork, Ireland, and Rochester Institute of Technology, USA.





