File {
  Grid      = "pn2d_minimal6.tdr"
  Doping    = "pn2d_minimal6.tdr"
  Parameter = "models.par"
  Plot      = "pn2d_minimal6_gate_des.tdr"
  Output    = "pn2d_minimal6_gate_des.log"
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
    Auger
    Avalanche(VanOverstraeten)
  )
  EffectiveIntrinsicDensity(
    OldSlotboom
  )
}

Plot {
  Doping
  DonorConcentration
  AcceptorConcentration
}

Math {
  RelErrControl
  Digits=5
  Iterations=100
  NotDamped=100
  Method=Blocked
}

Solve {
  Coupled(Iterations=100 LineSearchDamping=1e-4) {
    Poisson Electron Hole
  }
}
