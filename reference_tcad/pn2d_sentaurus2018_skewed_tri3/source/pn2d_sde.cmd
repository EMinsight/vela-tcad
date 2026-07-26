;==========================================================
; 2D PN junction - skewed Tri3 diagnostic mesh
; Sentaurus Structure Editor
;
; Purpose:
;   Supply acute scalene and obtuse Tri3 cells for avalanche
;   current-support diagnostics. The internal region edges preserve
;   the intended eight-triangle topology. This deck is diagnostic-only.
;
; Unit:
;   Length = um
;   Doping = cm^-3
;==========================================================

(sde:clear)
(sdegeo:set-default-boolean "ABA")

(define L 2.0)
(define H 0.5)
(define XJ 1.0)

;----------------------------------------------------------
; Eight explicit triangles. The PN junction remains aligned
; with the shared x=1 region boundary.
;----------------------------------------------------------

(sdegeo:create-polygon
  (list
    (position 0.0 0.0 0.0)
    (position 1.0 0.0 0.0)
    (position 0.4 0.18 0.0))
  "Silicon" "R.P0")

(sdegeo:create-polygon
  (list
    (position 1.0 0.0 0.0)
    (position 1.0 0.5 0.0)
    (position 0.4 0.18 0.0))
  "Silicon" "R.P1")

(sdegeo:create-polygon
  (list
    (position 1.0 0.5 0.0)
    (position 0.0 0.5 0.0)
    (position 0.4 0.18 0.0))
  "Silicon" "R.P2")

(sdegeo:create-polygon
  (list
    (position 0.0 0.5 0.0)
    (position 0.0 0.0 0.0)
    (position 0.4 0.18 0.0))
  "Silicon" "R.P3")

(sdegeo:create-polygon
  (list
    (position 1.0 0.0 0.0)
    (position 2.0 0.0 0.0)
    (position 1.6 0.32 0.0))
  "Silicon" "R.N0")

(sdegeo:create-polygon
  (list
    (position 2.0 0.0 0.0)
    (position 2.0 0.5 0.0)
    (position 1.6 0.32 0.0))
  "Silicon" "R.N1")

(sdegeo:create-polygon
  (list
    (position 2.0 0.5 0.0)
    (position 1.0 0.5 0.0)
    (position 1.6 0.32 0.0))
  "Silicon" "R.N2")

(sdegeo:create-polygon
  (list
    (position 1.0 0.5 0.0)
    (position 1.0 0.0 0.0)
    (position 1.6 0.32 0.0))
  "Silicon" "R.N3")

;----------------------------------------------------------
; Contacts
;----------------------------------------------------------

(sdegeo:define-contact-set
  "Anode" 4.0 (color:rgb 1 0 0) "##")
(sdegeo:define-contact-set
  "Cathode" 4.0 (color:rgb 0 0 1) "##")

(sdegeo:set-current-contact-set "Anode")
(sdegeo:define-2d-contact
  (find-edge-id (position 0.0 (/ H 2.0) 0.0))
  "Anode")

(sdegeo:set-current-contact-set "Cathode")
(sdegeo:define-2d-contact
  (find-edge-id (position L (/ H 2.0) 0.0))
  "Cathode")

;----------------------------------------------------------
; Doping
;----------------------------------------------------------

(sdedr:define-refeval-window
  "P.Window" "Rectangle"
  (position 0.0 0.0 0.0)
  (position XJ H 0.0))
(sdedr:define-refeval-window
  "N.Window" "Rectangle"
  (position XJ 0.0 0.0)
  (position L H 0.0))

(sdedr:define-constant-profile
  "P.Doping" "BoronActiveConcentration" 1e17)
(sdedr:define-constant-profile
  "N.Doping" "PhosphorusActiveConcentration" 1e17)
(sdedr:define-constant-profile-placement
  "P.Place" "P.Doping" "P.Window")
(sdedr:define-constant-profile-placement
  "N.Place" "N.Doping" "N.Window")

;----------------------------------------------------------
; Keep the explicit region vertices; do not introduce a
; structured right-triangle lattice.
;----------------------------------------------------------

(sdedr:define-refinement-size
  "Global.Mesh"
  2.0 0.5
  2.0 0.5)

(sdedr:define-refinement-window
  "Global.Window" "Rectangle"
  (position 0.0 0.0 0.0)
  (position L H 0.0))

(sdedr:define-refinement-placement
  "Global.Mesh.Place"
  "Global.Mesh"
  "Global.Window")

(sde:build-mesh "snmesh" "" "pn2d")
