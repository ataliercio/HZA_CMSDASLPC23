# HZA_CMSDASLPC23
Long exercise for a search of two neutral Higgs bosons through the H -> Z A -> llττ  process

## Reccomended short exercise

## Introduction

The final state llττ can also be exploited in the search for a new resonance decaying into a lighter, new resonance and a Z boson, where the light resonance would decay into a pair of τ leptons. This scenario is mainly inspired by models with two Higgs doublets (2HDM) [23]. In a 2HDM model two doublets of Higgs fields are present, as opposed to the single one of the SM, giving rise to two neutral CP-even scalars (h, H), one neutral CP-odd pseudo-scalar (A), and two charged scalars (H+, H-). The light h scalar is generally identified with the 125 GeV scalar resonance discovered in 2012 and in order to be consistent with the SM-like properties of the latter, one must consider a limited part of the parameter space (the so-called “alignment limit”). The 2HDM models are attractive because they allow generation of the observed asymmetry in the Universe between matter and anti-matter [24]. Another important motivation is Supersymmetry [25], which belongs to the broad class of 2HDMs. Supersymmetric models can provide an elegant solution to the naturalness issue. Supersymmetry is also attractive because it generally includes a candidate for particle dark matter and can possibly lay down the bases for the unification of forces. Axion models [26], which would explain how the strong interaction does not violate the CP symmetry, would also give rise to an effective low-energy theory with two Higgs doublets. Finally, it has also been noted [27] that certain realizations of 2HDMs can accommodate the muon g-2 [28] anomaly without violating the present theoretical and experimental constraints.
According to Ref [29], the pseudo-scalar boson A would decay into HZ with the highest branching ratio. For this reason, A→ZH appears as the most promising exotic decay channel and, in addition, H→ ττ would represent a clean final state. This search has not yet been investigated with Run-2 data and would differ from the hZ channel described previously by the presence of an intermediate heavy resonance (A or H) that will boost the H (or A) and the Z bosons, and by the different and unknown mass of the (pseudo)scalar boson H (A).

## Event selection

The final state that we are investigatin has 4 leptons: 2 muons coming from the Z decay and 2 taus coming from the A decay.

The first step of the analysis is the reconstruction of the Z into mumu, in order to do so you have to selected:

- 2 muons ...

After the Z selection, we have to select the taus coming from the A, depending on the tau decay process and selecting the three most sensitive categories we can divide the events into:

- 1mu1tau region
- 1ele1tau region
- 2tau region

You have to implement both of this selections in the analysis macro.

## Fake rate estimation

The reducible background (Z + X) originates from processes that contain one or more non-prompt leptons. The main sources of non-prompt leptons are non-isolated electrons and muons coming from decays of heavy-flavour mesons, mis-reconstructed jets (usually originating from light-flavour quarks) and electrons from γ conversions. In the following discussion we will consider a fake muon any jet mis-reconstructed as a lepton and any lepton originating from a heavy meson decay. Similarly, any electron originating from a photon conversion will be considered a fake electron.
The rate of these background processes is estimated by measuring the fe and fμ probabilities for fake electrons and fake muons which do pass the loose selection criteria to also pass the final selection criteria (defined in Section 5.3). These probabilities, referred to as fake rates, are applied in dedicated control samples in order to extract the expected background yield in the signal region.

## ML application

After the FR estimation

## Extracting results

To extract the final limit we have to use combine + combine harvester:

```
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_2_13
cd CMSSW_10_2_13/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
git clone https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
scram b
```



