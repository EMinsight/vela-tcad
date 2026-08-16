; Two-dimensional n-type silicon Schottky diode translated from the Charon
; schottky_contacts/n_diode benchmark.  Coordinates are in micrometres.

(sde:clear)
(sdegeo:set-default-boolean "ABA")

(define L 1.0)
(define H 1.0)

(sdegeo:create-rectangle
  (position 0.0 0.0 0.0)
  (position L H 0.0)
  "Silicon"
  "R.Silicon")

(sdegeo:define-contact-set
  "anode" 4.0 (color:rgb 1.0 0.0 0.0) "##")
(sdegeo:define-contact-set
  "cathode" 4.0 (color:rgb 0.0 0.0 1.0) "##")

(sdegeo:set-current-contact-set "anode")
(sdegeo:define-2d-contact
  (find-edge-id (position 0.0 (/ H 2.0) 0.0))
  "anode")

(sdegeo:set-current-contact-set "cathode")
(sdegeo:define-2d-contact
  (find-edge-id (position L (/ H 2.0) 0.0))
  "cathode")

; Charon n_diode: uniform donor concentration 1e16 cm^-3.
(sdedr:define-constant-profile
  "Const.Donor" "PhosphorusActiveConcentration" 1.0e16)
(sdedr:define-constant-profile-region
  "Place.Donor" "Const.Donor" "R.Silicon")

; The transport direction is x.  Resolve the Schottky depletion region more
; finely than the neutral bulk while keeping the reference inexpensive.
(sdedr:define-refinement-size
  "Global.Mesh" 0.05 0.10 0.02 0.05)
(sdedr:define-refinement-region
  "Place.Global" "Global.Mesh" "R.Silicon")

(sdedr:define-refeval-window
  "Anode.Window" "Rectangle"
  (position 0.0 0.0 0.0)
  (position 0.25 H 0.0))
(sdedr:define-refinement-size
  "Anode.Mesh" 0.02 0.10 0.005 0.05)
(sdedr:define-refinement-placement
  "Place.Anode" "Anode.Mesh" "Anode.Window")

(sde:build-mesh "snmesh" "" "schottky_n")
