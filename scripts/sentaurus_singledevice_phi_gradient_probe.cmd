File {
  Grid      = "n2_msh.tdr"
  Plot      = "eq231_phi_gradient_probe.tdr"
  Parameter = "sdevice.par"
  Current   = "eq231_phi_gradient_probe.plt"
  Output    = "eq231_phi_gradient_probe.log"
}

Electrode {
  { Name="source" Voltage=0.0 }
  { Name="drain" Voltage=0.1 }
  { Name="gate" Voltage=2.2 }
  { Name="substrate" Voltage=0.0 }
}

Physics {
  eQuantumPotential
  EffectiveIntrinsicDensity(OldSlotboom)
}

Plot {
  ElectrostaticPotential
  ElectricField/Vector
  eQuantumPotential
}

Math {
  RelErrControl
  Digits=5
}

Solve {
  Load(FilePrefix="lin_phi_as_potential")
  Plot(FilePrefix="eq231_phi_gradient_probe")
}
