from bamboo.analysismodules import NanoAODHistoModule
from bamboo.treedecorators import NanoAODDescription
from bamboo.scalefactors import binningVariables_nano,lumiPerPeriod_default

from bamboo import treedecorators as td
from bamboo import treefunctions as op
from bamboo import scalefactors

from bamboo.plots import Plot, EquidistantBinning, CutFlowReport
from bamboo import treefunctions as op

from bamboo.treeoperations import Const

from itertools import chain
from functools import partial

import os
from bamboo.root import loadDependency 

if os.environ["VIRTUAL_ENV"] == "":
    print("$VIRTUAL_ENV is not set. Please, activate your virtual environment")
    exit(-1)

print(os.environ["VIRTUAL_ENV"])
print("{0}/../ClassicSVfit/interface".format(os.environ["VIRTUAL_ENV"]))

loadDependency( headers="SVfitBambooProducer.h",
                includePath="{0}/../ClassicSVfit/interface".format(os.environ["VIRTUAL_ENV"]),
                dynamicPath="{0}/../build-ClassicSVfit-FastMTT/src".format(os.environ["VIRTUAL_ENV"]), 
                libraries="ClassicSVfit")
                
svFitBambooProducer = op.define("SVfitBambooProducer", 'auto <<name>> = SVfitBambooProducer();')


class category:
    def __init__(self,nMuons=0, nElectrons=0, nTaus=0):
        #FIXME: add checks on total number of leptons

        self.mu  = nMuons
        self.ele = nElectrons
        self.tau = nTaus  
        self._name = self.__str__()

    def __str__(self):
        strmu  = "" if not self.mu else f"{self.mu}mu"
        strele = "" if not self.ele else f"{self.ele}ele"
        strtau = "" if not self.tau else f"{self.tau}tau"  

        return strmu+strele+strtau

    def nMuons(self):
        return self.mu
    def nElectrons(self):
        return self.ele
    def nTaus(self):
        return self.tau
    def name(self):
        return self._name


class ZhToLLTauTau(NanoAODHistoModule):
    """ Module for Zh lltautau analysis"""
    def __init__(self, args):
        super(ZhToLLTauTau, self).__init__(args)

    def isMC(self, sampleName):
        return sampleName.split("_")[0] not in ("SingleMuon", "DoubleMuon", "SingleEGamma")
    def prepareTree(self, tree, sample=None, sampleCfg=None, backend=None):
        from bamboo.treedecorators import nanoRochesterCalc
        from bamboo.analysisutils import configureJets
        from bamboo.analysisutils import configureRochesterCorrection
        tree,noSel,be,lumiArgs = NanoAODHistoModule.prepareTree(self, tree, sample=sample, sampleCfg=sampleCfg, description=NanoAODDescription.get('v7', year='2018', isMC=True, systVariations=[nanoRochesterCalc]), backend=backend)
        era = sampleCfg["era"]
        if era == "2018":
            configureRochesterCorrection(tree._Muon, "./RoccoR2018UL.txt", isMC=self.isMC(sample), backend=be)

        return tree,noSel,be,lumiArgs
 
    def definePlots(self,t, noSel, sample=None, sampleCfg=None):

        plots = []

        if self.isMC(sample):
            noSel = noSel.refine("mcWeight", weight=t.genWeight)

       
        yields = CutFlowReport("yields", printInLog=True)
        era = sampleCfg["era"]
        # Raccomanded HLT Path Muon POG: https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2018#Recommended_trigger_paths_for_20
        if era == "2018":
            singleMuonTriggers= [ t.HLT.IsoMu24 ] # Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8


        isoMuFilterMask = 0xA
        #dimuonTriggers = [t.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v]
        # TESTING
        hasTriggerFired = noSel.refine("passSingleMuonHLT", cut=op.OR(*(singleMuonTriggers))) #+ mueTriggers+ emuTriggers)))
        #noSel_dimuon = noSel.refine("passDiMuonHLT", cut=op.OR(*(dimuonTriggers)))#+ mueTriggers+ emuTriggers)))
        # noSel_diele  = noSel.refine("passDiElectronHLT", cut=op.OR(*(dielectronTriggers)))#+ mueTriggers+ emuTriggers)))

        plots = []

        triggerObj = op.select(t.TrigObj, lambda trgObj: op.AND( trgObj.id == 13,
                                                                 (trgObj.filterBits & isoMuFilterMask) )) # + 1mu     trgObj.filterBits == 8,   # Iso

        muons = op.sort(op.select(t.Muon, lambda mu : op.AND(
            mu.pt > 10.,
            mu.mediumId,
            op.abs(mu.eta) < 2.4,
            op.abs(mu.pfRelIso04_all) < 0.40,
            op.abs(mu.dxy) < 0.5,
            op.abs(mu.dz ) < 1.,
            op.sqrt(mu.dxy**2 + mu.dz**2)/op.sqrt(mu.dxyErr**2+mu.dzErr**2) < 4, ## SIP3D
            )), lambda mu : -mu.pt)
        electrons = op.sort(op.select(t.Electron, lambda el : op.AND(
            el.pt > 10.,
            el.mvaFall17V2noIso_WP90,
            op.abs(el.eta) < 2.5,
            op.abs(el.pfRelIso03_all) < 0.40,
            op.abs(el.dxy) < 0.5,
            op.abs(el.dz ) < 1.,
            op.sqrt(el.dxy**2 + el.dz**2)/op.sqrt(el.dxyErr**2+el.dzErr**2) < 4, ## SIP3D
            )), lambda el : -el.pt)

        taus = op.select(t.Tau, lambda tau : op.AND(tau.p4.Pt() > 20., 
            op.abs(tau.p4.Eta()) < 2.4 ,
            tau.idDeepTau2017v2p1VSe >= 1,
            tau.idDeepTau2017v2p1VSmu >= 1, 
            tau.idDeepTau2017v2p1VSjet >= 1,
            tau.decayMode, 
            tau.dz < 0.2))

        cleanedTaus = op.select(taus, lambda it : op.AND(
                                                  op.NOT(op.rng_any(electrons, lambda ie : op.deltaR(it.p4, ie.p4) < 0.3 )),
                                                  op.NOT(op.rng_any(muons, lambda im : op.deltaR(it.p4, im.p4) < 0.3 ))
                                ))
 
        #deepTau wp ?
        wp = {"Vloose": 4,
              "loose": 8,
              "medium": 16}

        tau_wps = {}

        taus_vsmu = op.select(cleanedTaus, lambda tau : tau.idDeepTau2017v2p1VSmu >= 8)
        
        for jetKey in wp:
            
            for eKey in wp:

                print("VSjet wp : ", jetKey)
                print("VSe   wp : ", eKey)

                name = "VSjet"+jetKey+"VSe"+eKey
                tau_wps[name] = op.select(taus_vsmu, lambda tau : op.AND(tau.idDeepTau2017v2p1VSjet >= wp[jetKey],
                                                                         tau.idDeepTau2017v2p1VSe >= wp[eKey]))                 
            

         
        mZ = 91.1876
        pairsMuMu = op.combine(muons, N=2, pred=lambda l1,l2 : op.AND(l1.charge != l2.charge,
                                                                      l1.p4.Pt()>24))
        bestZ = op.rng_min_element_by(pairsMuMu, lambda ll : op.abs(op.invariant_mass(ll[0].p4, ll[1].p4)-mZ))        

        hasZmm      = hasTriggerFired.refine("hasZmm",  cut=[op.rng_len(pairsMuMu) > 0] )
        hasZmmTight = hasZmm.refine("hasZmmTight",  cut=[ bestZ[0].pfRelIso04_all < 0.15,
                                                          bestZ[0].tightId,
                                                          bestZ[1].pfRelIso04_all < 0.15,
                                                          bestZ[0].mediumId] )
       
        hasZmmTight_mZcut = hasZmm.refine("hasZmmTight_mZcut",  cut=[ bestZ[0].pfRelIso04_all < 0.15,
                                                                    bestZ[0].tightId,
                                                                    bestZ[1].pfRelIso04_all < 0.15,
                                                                    bestZ[0].mediumId,
                                                                    op.invariant_mass(bestZ[0].p4, bestZ[1].p4) > 76.,
                                                                    op.invariant_mass(bestZ[0].p4, bestZ[1].p4) < 106.] )
 
        plots += self.plotPairs(bestZ, hasZmm, "osMuMu")
        plots += self.plotPairs(bestZ, hasZmmTight, "osMuMu_Tight")
        plots += self.plotPairs(bestZ, hasZmmTight_mZcut, "osMuMu_Tight_mZcut")

        categories = []
        categories.append(category(nMuons=1, nElectrons=0, nTaus=1))
        categories.append(category(nMuons=0, nElectrons=1, nTaus=1))
        categories.append(category(nMuons=0, nElectrons=0, nTaus=2))

        plots.append(yields)
        yields.add(noSel, title="NoSel")
        yields.add(hasTriggerFired, title="TriggerFired")
        yields.add(hasZmm, title="Zloose")
        yields.add(hasZmmTight, title="Ztight")
        yields.add(hasZmmTight_mZcut, title="Ztight_mZcut")


        for cat in categories:

            otherLeptons = op.select(muons, partial(lambda l,oz=None : op.AND(l.idx != oz[0].idx, l.idx != oz[1].idx), oz=bestZ))

            for key in tau_wps: 
                print(key)
                ttau = tau_wps[key]     
 
                if cat.nMuons() > 0:
                   oslep3lep4 = op.combine((otherLeptons, ttau), pred=lambda mu,tau : mu.charge != tau.charge)
                elif cat.nMuons() == 0 and cat.nElectrons() == 0:
                   oslep3lep4 = op.combine(ttau, N=2, pred=lambda l1, l2 : l1.charge != l2.charge)
                else:
                   oslep3lep4 = op.combine((electrons, ttau), pred=lambda ele,tau : ele.charge != tau.charge)
                

            
            #bestH = op.rng_max_element_by(oslep3lep4, lambda ll : ll[0].p4.Pt()+ll[1].p4.Pt())

                #hasMoreLeps  = hasZmmTight.refine(f"hasMoreLeps_{cat}",  cut=[op.rng_len(otherLeptons) >= cat.nMuons()] )
                hasSeconPair = hasZmmTight.refine(f"hasSeconPairCategory_{cat}_{key}", cut=[op.rng_len(oslep3lep4) > 0,
                                                                             op.rng_len(otherLeptons) == cat.nMuons(),
                                                                             op.rng_len(ttau) >= cat.nTaus() ] )

                hasSeconPair_mZcut = hasZmmTight_mZcut.refine(f"hasSeconPairCategory_{cat}_{key}_mZcut", cut=[op.rng_len(oslep3lep4) > 0,
                                                                             op.rng_len(otherLeptons) == cat.nMuons(),
                                                                             op.rng_len(ttau) >= cat.nTaus() ] )
 
                yields.add(hasSeconPair, title=f"With a higgs pair cadidate {cat}, wp {key}")
                yields.add(hasSeconPair_mZcut, title=f"With a higgs pair cadidate {cat}, wp {key} - Z mass cut ")


            #plots.append(Plot.make1D(f"h_{cat}_3lep_pT", otherLeptons[0].p4.Pt(), hasMoreLeps, EquidistantBinning(100, 10., 200.), title=" 3rd lepton pT", xTitle= "pT (GeV)"))
            #mass = svFitBambooProducer.computeSVfit(t.MET.covXX, t.MET.covXY, t.MET.covXY, t.MET.covYY,
            #                                        t.MET.pt*op.cos(t.MET.phi), t.MET.pt*op.cos(t.MET.phi),
            #                                        bestH[0].p4.Pt(), bestH[0].p4.Eta(), bestH[0].p4.Phi(), bestH[0].p4.M(),
            #                                        bestH[1].p4.Pt(), bestH[1].p4.Eta(), bestH[1].p4.Phi(), bestH[1].p4.M(),
            #                                        Const("std::string", f"\"{cat}\"") )
            #massFast = svFitBambooProducer.computeFastMTT(t.MET.covXX, t.MET.covXY, t.MET.covXY, t.MET.covYY,
            #                                              t.MET.pt*op.cos(t.MET.phi), t.MET.pt*op.cos(t.MET.phi),
            #                                              bestH[0].p4.Pt(), bestH[0].p4.Eta(), bestH[0].p4.Phi(), bestH[0].p4.M(),
            #                                              bestH[1].p4.Pt(), bestH[1].p4.Eta(), bestH[1].p4.Phi(), bestH[1].p4.M(),
            #                                              Const("std::string", f"\"{cat}\"") )


            #plots.append(Plot.make1D(f"h_SVfit_mll_{cat}", op.c_float(mass), hasSeconPair, EquidistantBinning(100, 20., 200.), title=" pair SVfit mass", xTitle= "mll (GeV)"))
            #plots.append(Plot.make1D(f"h_SVfit_Fast_mll_{cat}", op.c_float(massFast), hasSeconPair, EquidistantBinning(100, 20., 200.), title=" pair SVfit mass", xTitle= "mll (GeV)"))

        
            #yields.add(hasMoreLeps, title="With additional leptons")
            #yields.add(hasSeconPair, title=f"With a higgs pair cadidate {cat}")
            #yields.add(hasSeconPair_mZcut, title=f"With a higgs pair cadidate {cat} - Z mass cut ")

        return plots

    def plotPairs(self, pair, sel, category):

            plots = []
            plots.append(Plot.make1D(f"h_{category}_mll", op.invariant_mass(pair[0].p4, pair[1].p4), sel, EquidistantBinning(100, 20., 200.), title=" pair invariant mass", xTitle= "mll (GeV)"))
            plots.append(Plot.make1D(f"h_{category}_lep1_pT", pair[0].p4.Pt(), sel, EquidistantBinning(100, 20., 200.), title=" lepton1 pT", xTitle= "pT (GeV)"))
            plots.append(Plot.make1D(f"h_{category}_lep2_pT", pair[1].p4.Pt(), sel, EquidistantBinning(100, 10., 200.), title=" lepton2 pT", xTitle= "pT (GeV)"))
            plots.append(Plot.make1D(f"h_{category}_lep1_eta", pair[0].p4.Eta(),  sel, EquidistantBinning(100, -3, 3), title=" lepton1 Eta", xTitle= "eta"))
            plots.append(Plot.make1D(f"h_{category}_lep2_eta", pair[1].p4.Eta(),  sel, EquidistantBinning(100, -3, 3), title=" lepton2 Eta", xTitle= "eta"))
            plots.append(Plot.make1D(f"h_{category}_deltaR", op.deltaR(pair[0].p4,pair[1].p4),  sel, EquidistantBinning(20, 0, 10), title=" Delta R", xTitle= "deltaR") )
            plots.append(Plot.make1D(f"h_{category}_deltaPhi", op.deltaPhi(pair[0].p4,pair[1].p4),  sel, EquidistantBinning(50, -3.3, 3.3), title=" Delta R", xTitle= "deltaR") )

            return plots

