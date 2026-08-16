set dataset [load_file lin_state_0020_des.tdr -name Eq231Probe]
export_variables -dataset $dataset -filename eq231_probe_variables.csv -overwrite {
    X Y
    ElectrostaticPotential ElectricField-X ElectricField-Y
    eQuasiFermiPotential eGradQuasiFermi-X eGradQuasiFermi-Y
    ConductionBandEnergy ElectronAffinity BandGap BandgapNarrowing
    eQuantumPotential eDensity
}
exit
