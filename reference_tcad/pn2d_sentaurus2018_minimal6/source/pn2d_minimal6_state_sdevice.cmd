File {
  Grid      = "pn2d_minimal6.tdr"
  Doping    = "pn2d_minimal6.tdr"
  Parameter = "models.par"
  Plot      = "pn2d_minimal6_state___BIAS_TAG__.tdr"
  Current   = "pn2d_minimal6_state___BIAS_TAG__.plt"
  Output    = "pn2d_minimal6_state___BIAS_TAG__.log"
}

Electrode {
  { Name="Anode"   Voltage=0.0 }
  { Name="Cathode" Voltage=0.0 }
}

Physics {
  Mobility(
    DopingDependence
    HighFieldSaturation
  )
  Recombination(
    SRH
    Avalanche(VanOverstraeten)
  )
  EffectiveIntrinsicDensity(
    OldSlotboom
  )
}

Plot {
  Potential
  eQuasiFermi
  hQuasiFermi
  eDensity
  hDensity
  ElectricField
  ElectricField/Vector
  eCurrent
  hCurrent
  TotalCurrent
  eCurrentDensity/Vector
  hCurrentDensity/Vector
  TotalCurrentDensity/Vector
  Doping
  DonorConcentration
  AcceptorConcentration
  SpaceCharge
  SRHRecombination
  eAlphaAvalanche
  hAlphaAvalanche
  AvalancheGeneration
  eMobility
  hMobility
  eVelocity
  hVelocity
  eIonIntegral
  hIonIntegral
  MeanIonIntegral
}

Math {
  Extrapolate
  RelErrControl
  Digits=5
  Iterations=80
  NotDamped=100
  Method=Blocked
}

Solve {
  Coupled(Iterations=100) { Poisson }
  Coupled(Iterations=100) { Poisson Electron Hole }
  Quasistationary(
    InitialStep=1e-4
    MinStep=1e-10
    MaxStep=0.05
    Increment=1.2
    Decrement=2.0
    Goal { Name="Anode" Voltage=__TARGET_BIAS_V__ }
  ) {
    Coupled { Poisson Electron Hole }
  }
}
