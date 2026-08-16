File {
  Grid      = "n2_msh.tdr"
  Plot      = "deep_off_des.tdr"
  Parameter = "sdevice.par"
  Current   = "deep_off.plt"
  Output    = "deep_off.log"
}

Electrode {
  { Name="source"    Voltage=0.0 }
  { Name="drain"     Voltage=0.0 }
  { Name="gate"      Voltage=-0.5 }
  { Name="substrate" Voltage=0.0 }
}

Physics {
  eQuantumPotential
  EffectiveIntrinsicDensity(OldSlotboom)
  Mobility(
    DopingDep
    eHighFieldsaturation(GradQuasiFermi)
    hHighFieldsaturation(GradQuasiFermi)
    Enormal
  )
  Recombination(
    SRH(DopingDep TempDependence)
  )
}

Plot {
  eDensity hDensity
  TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
  eMobility/Element hMobility/Element
  eQuasiFermi hQuasiFermi
  ElectricField/Vector Potential SpaceCharge
  Doping DonorConcentration AcceptorConcentration
  SRH
  BandGap BandGapNarrowing Affinity
  ConductionBand ValenceBand eQuantumPotential
}

Math {
  Extrapolate
  RelErrControl
  Digits=5
  ErrRef(electron)=1.e10
  ErrRef(hole)=1.e10
  Iterations=20
  Notdamped=100
}

Solve {
  NewCurrentPrefix="deep_off_init_"
  Coupled(Iterations=100) { Poisson eQuantumPotential }
  Coupled { Poisson Electron Hole eQuantumPotential }

  Quasistationary(
    InitialStep=0.01 MinStep=1e-5 MaxStep=0.1
    Goal { Name="drain" Voltage=0.1 }
  ) { Coupled { Poisson Electron Hole eQuantumPotential } }

  Plot(FilePrefix="deep_off")
}
