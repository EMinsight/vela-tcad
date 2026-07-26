Title "Skewed Tri3 constrained diagnostic mesh"

Controls {
}

IOControls {
  EnableSections
}

Delaunizer {
  type = constrained
  maxAngle = 180
}

Definitions {
  Constant "P.Doping" {
    Species = "BoronActiveConcentration"
    Value = 1e+17
  }
  Constant "N.Doping" {
    Species = "PhosphorusActiveConcentration"
    Value = 1e+17
  }
  Refinement "Global.Mesh" {
    MaxElementSize = ( 2 0.5 )
    MinElementSize = ( 2 0.5 )
  }
}

Placements {
  Constant "P.Place" {
    Reference = "P.Doping"
    EvaluateWindow {
      Element = Rectangle [(0 0) (1 0.5)]
    }
  }
  Constant "N.Place" {
    Reference = "N.Doping"
    EvaluateWindow {
      Element = Rectangle [(1 0) (2 0.5)]
    }
  }
  Refinement "Global.Mesh.Place" {
    Reference = "Global.Mesh"
    RefineWindow = Rectangle [(0 0) (2 0.5)]
  }
}
