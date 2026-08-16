File {
  Grid       = "n2_msh.tdr"
  Plot       = "eq231_theta_zero_probe.tdr"
  Parameter  = "sentaurus_singledevice_eq231_theta_zero.par"
  Current    = "eq231_theta_zero_probe.plt"
  Output     = "eq231_theta_zero_probe.log"
  NewtonPlot = "eq231_theta_zero_newton_%d.tdr"
}

Electrode {
  { Name="source"    Voltage=0.0 }
  { Name="drain"     Voltage=0.1 }
  { Name="gate"      Voltage=2.2 }
  { Name="substrate" Voltage=0.0 }
}

Physics {
  eQuantumPotential(ParameterSetName="theta_zero")
  EffectiveIntrinsicDensity(OldSlotboom)
}

Plot {
  ElectrostaticPotential
  ConductionBand
  eQuantumPotential
}

Math {
  RelErrControl
  Digits=12
  Iterations=1
  Notdamped=100
  CNormPrint
  NewtonPlot(Residual Error Update Plot)
}

Solve {
  Load(FilePrefix="lin_state_0020")
  Coupled(Iterations=1) { eQuantumPotential }
}
