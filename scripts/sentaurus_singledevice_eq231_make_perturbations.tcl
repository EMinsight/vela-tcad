# Generate one-state eQuantumPotential perturbations directly inside the
# Sentaurus VM.  Keeping this operation next to the licensed solver avoids
# transferring proprietary TDR state files into the VM.

set input_file "lin_state_0020_des.tdr"
set delta 1.0e-5

# Each entry is: output node {{region_index region_value_offset} ...}
set cases {
  {lin_eq231_perturbed_1848_des.tdr 1848 {{2 1}}}
  {lin_eq231_perturbed_1847_des.tdr 1847 {{2 0} {4 239}}}
  {lin_eq231_perturbed_3540_des.tdr 3540 {{2 95} {4 410}}}
  {lin_eq231_perturbed_1849_des.tdr 1849 {{2 2}}}
  {lin_eq231_perturbed_1852_des.tdr 1852 {{2 4}}}
  {lin_eq231_perturbed_1859_des.tdr 1859 {{2 9}}}
  {lin_eq231_perturbed_1861_des.tdr 1861 {{2 11}}}
  {lin_eq231_perturbed_1845_des.tdr 1845 {{4 237}}}
  {lin_eq231_perturbed_1846_des.tdr 1846 {{4 238}}}
  {lin_eq231_perturbed_1851_des.tdr 1851 {{2 3} {4 241}}}
  {lin_eq231_perturbed_1847r2_des.tdr 1847r2 {{2 0}}}
  {lin_eq231_perturbed_1847r4_des.tdr 1847r4 {{4 239}}}
  {lin_eq231_perturbed_1851r2_des.tdr 1851r2 {{2 3}}}
  {lin_eq231_perturbed_1851r4_des.tdr 1851r4 {{4 241}}}
  {lin_eq231_perturbed_2074_des.tdr 2074 {{0 240}}}
  {lin_eq231_perturbed_2071_des.tdr 2071 {{0 238}}}
  {lin_eq231_perturbed_2072_des.tdr 2072 {{0 239} {4 355}}}
  {lin_eq231_perturbed_2075_des.tdr 2075 {{0 241} {4 357}}}
  {lin_eq231_perturbed_2119_des.tdr 2119 {{0 250}}}
  {lin_eq231_perturbed_2077_des.tdr 2077 {{0 243}}}
  {lin_eq231_perturbed_2121_des.tdr 2121 {{0 252}}}
}

foreach perturbation $cases {
  lassign $perturbation output_file global_node changes
  TdrFileOpen $input_file
  set ngeo [TdrFileGetNumGeometry $input_file]
  set applied 0

  foreach change $changes {
    lassign $change target_region target_offset
    set found 0
    for {set igeo 0} {$igeo < $ngeo} {incr igeo} {
      set nstates [TdrGeometryGetNumState $input_file $igeo]
      set nregions [TdrGeometryGetNumRegion $input_file $igeo]
      if {$target_region >= $nregions} {
        continue
      }
      for {set istate 0} {$istate < $nstates} {incr istate} {
        set ndata [TdrRegionGetNumDataset $input_file $igeo $target_region $istate]
        for {set idata 0} {$idata < $ndata} {incr idata} {
          set dsname [TdrDatasetGetName $input_file $igeo $target_region $istate $idata]
          if {![string equal -nocase $dsname "eQuantumPotential"]} {
            continue
          }
          set nvalue [TdrDatasetGetNumValue $input_file $igeo $target_region $istate $idata]
          if {$target_offset >= $nvalue} {
            error "offset $target_offset is outside dataset $dsname (size $nvalue)"
          }
          set original [TdrDataGetComponent $input_file $igeo $target_region $istate $idata $target_offset 0 0]
          set modified [expr {$original + $delta}]
          TdrDataSetComponent $input_file $igeo $target_region $istate $idata $target_offset 0 0 $modified
          set regname [TdrRegionGetName $input_file $igeo $target_region]
          puts "node=$global_node region=$target_region name=$regname offset=$target_offset before=$original after=$modified"
          incr applied
          set found 1
        }
      }
    }
    if {!$found} {
      error "eQuantumPotential not found for node $global_node region $target_region"
    }
  }

  if {$applied != [llength $changes]} {
    error "node $global_node expected [llength $changes] changes, applied $applied"
  }
  TdrFileSave $input_file $output_file
  TdrFileClose $input_file
  puts "saved $output_file"
}
