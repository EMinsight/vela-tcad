File {
  Grid      = "n2_msh.tdr"
  Plot      = "eq231_box_measure_probe.tdr"
  Parameter = "sdevice.par"
  Current   = "eq231_box_measure_probe.plt"
  Output    = "eq231_box_measure_probe.log"
}

Electrode {
  { Name="source" Voltage=0.0 }
  { Name="drain" Voltage=0.0 }
  { Name="gate" Voltage=-0.5 }
  { Name="substrate" Voltage=0.0 }
}

Physics {
  EffectiveIntrinsicDensity(OldSlotboom)
}

Plot {
  BM_AngleElements
  BM_CoeffIntersectionNonDelaunayElements
  BM_ElementVolume
  BM_IntersectionNonDelaunayElements
  BM_VolumeIntersectionNonDelaunayElements
}

Math {
  AverageBoxMethod
  BoxMeasureFromFile(GrdNumbering)
}

Solve {
  Coupled(Iterations=1) { Poisson }
  Plot(FilePrefix="eq231_box_measure_probe")
}
