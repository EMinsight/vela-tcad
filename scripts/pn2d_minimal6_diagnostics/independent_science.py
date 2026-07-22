"""Independent raw-observation scientific replay for the Minimal6 verifier.

This module intentionally does not import the audited field, transport,
avalanche, replacement, or report evaluators. It is verifier-only.
"""
from __future__ import annotations
import math
import statistics
from .inverse_contracts import AcceptanceThresholds, Identifiability, SampleStatus, SupportKind

Q = 1.602176634e-19
DEPENDENCIES = (
    "gradient_recovery", "mobility", "current_semantics",
    "impact_driving_field", "alpha_law", "geometric_integration",
    "source_to_node_mapping",
)

def _key(value):
    return type(value).__name__, str(value)

def _first_invalid(statuses):
    return next((status for status in statuses if status is not SampleStatus.VALID), None)

def _usable(row, unit):
    if row is None:
        return None, SampleStatus.MISSING_FIELD
    if row.status is not SampleStatus.VALID:
        return None, row.status
    if row.value_si is None:
        return None, SampleStatus.MISSING_FIELD
    if row.unit_si != unit:
        return None, SampleStatus.INVALID_UNIT
    try:
        value = float(row.value_si)
    except (TypeError, ValueError):
        return None, SampleStatus.NONFINITE
    return (value, SampleStatus.VALID) if math.isfinite(value) else (None, SampleStatus.NONFINITE)

def _gradient(points, values):
    (x0, y0), (x1, y1), (x2, y2) = points
    f0, f1, f2 = values
    det = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)
    if not math.isfinite(det) or abs(det) <= 1e-300:
        raise ValueError("independent replay found degenerate triangle")
    return (((f1-f0)*(y2-y0) - (f2-f0)*(y1-y0))/det,
            ((x1-x0)*(f2-f0) - (x2-x0)*(f1-f0))/det)

def _area(points, signed=False):
    (x0,y0),(x1,y1),(x2,y2) = points
    twice = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)
    if not math.isfinite(twice) or abs(twice) <= 1e-300:
        raise ValueError("independent replay found degenerate triangle")
    if signed and twice <= 0:
        raise ValueError("independent generation replay requires CCW triangles")
    return .5 * (twice if signed else abs(twice))

def _topology(mesh, topology):
    selected = mesh[topology] if topology in mesh else mesh
    cells = {str(cell): tuple(str(node) for node in nodes)
             for cell, nodes in selected["triangles"].items()}
    if "edges" in selected:
        edges = {str(edge): tuple(str(node) for node in nodes)
                 for edge, nodes in selected["edges"].items()}
    else:
        pairs = {tuple(sorted((nodes[i], nodes[(i+1)%3]), key=_key))
                 for nodes in cells.values() for i in range(3)}
        edges = {str(i): pair for i, pair in enumerate(sorted(
            pairs, key=lambda pair: (_key(pair[0]), _key(pair[1]))))}
    return cells, edges

def _weighted(cell_vectors, cells, coordinates, support_nodes):
    adjacent = [cell for cell in sorted(cells, key=_key)
                if all(node in cells[cell] for node in support_nodes)]
    areas = {cell: _area(tuple(coordinates[node] for node in cells[cell]))
             for cell in adjacent}
    total = sum(areas.values())
    if total <= 0:
        raise ValueError("independent replay support has no adjacent area")
    return (sum(areas[c]*cell_vectors[c][0] for c in adjacent)/total,
            sum(areas[c]*cell_vectors[c][1] for c in adjacent)/total)

def _tangent(start, end):
    dx, dy = end[0]-start[0], end[1]-start[1]
    length = math.hypot(dx,dy)
    if length == 0 or not math.isfinite(length):
        raise ValueError("independent replay found invalid edge")
    return dx/length, dy/length, length

def _project(vector, start, end):
    tx,ty,_ = _tangent(start,end)
    return vector[0]*tx + vector[1]*ty

def _edge_vector(component, start, end):
    tx,ty,_ = _tangent(start,end)
    return component*tx, component*ty

def _vector_error(candidate, reference, relative):
    if candidate is None or reference is None:
        return SampleStatus.MISSING_FIELD,None,SampleStatus.MISSING_FIELD,None,None
    values = tuple(float(value) for value in (*candidate,*reference))
    if not all(math.isfinite(value) for value in values):
        return SampleStatus.NONFINITE,None,SampleStatus.NONFINITE,None,None
    cmag,rmag = math.hypot(*candidate),math.hypot(*reference)
    if rmag == 0:
        return SampleStatus.GEOMETRIC_ZERO,None,SampleStatus.DIRECTION_UNDEFINED,None,None
    if not relative and cmag <= 0:
        return SampleStatus.BELOW_FLOOR,None,SampleStatus.DIRECTION_UNDEFINED,None,None
    error = abs(cmag-rmag)/rmag if relative else abs(math.log10(cmag/rmag))
    prediction = math.log10(cmag/rmag) if cmag > 0 else None
    if cmag == 0:
        return SampleStatus.VALID,error,SampleStatus.DIRECTION_UNDEFINED,None,prediction
    cosine = (candidate[0]*reference[0]+candidate[1]*reference[1])/(cmag*rmag)
    angle = math.degrees(math.acos(max(-1,min(1,cosine))))
    return SampleStatus.VALID,error,SampleStatus.VALID,angle,prediction

def _row(candidate,state,support,support_id,carrier,metric,factor,status,error,prediction,split,missing=()):
    result = {"candidate":candidate,"solver":state[0],"topology":state[1],
              "bias_V":state[2],"support_kind":support,"support_id":support_id,
              "carrier":carrier,"metric":metric,"factor":factor,"status":status,
              "error":error,"prediction_dex":prediction,"split":split}
    if missing:
        result["missing_independent_factors"] = tuple(missing)
    return result

def _vector_rows(candidate,state,support,support_id,carrier,value,reference,split,
                 field=False,invalid=None,reference_invalid=None):
    if invalid is not None:
        ms=ds=invalid; me=angle=prediction=None
    elif reference_invalid is not None:
        ms=reference_invalid
        ds=(SampleStatus.DIRECTION_UNDEFINED if reference_invalid in {
            SampleStatus.BELOW_FLOOR,SampleStatus.GEOMETRIC_ZERO} else reference_invalid)
        me=angle=prediction=None
    else:
        ms,me,ds,angle,prediction = _vector_error(value,reference,field)
    factor = "gradient_recovery" if field else "current_semantics"
    mm = "field_magnitude_relative" if field else "transport_abs_dex"
    dm = "field_direction_deg" if field else "transport_direction_deg"
    return [_row(candidate,state,support,support_id,carrier,mm,factor,ms,me,prediction,split),
            _row(candidate,state,support,support_id,carrier,dm,factor,ds,angle,prediction,split)]

def _groups(rows, quantities):
    groups={}
    for row in rows:
        if row.support_kind is SupportKind.NODE and row.quantity in quantities:
            groups.setdefault((row.solver,row.topology,float(row.bias_V)),[]).append(row)
    return groups

def _index(rows):
    return {(str(row.support_id),row.quantity,row.component):row for row in rows}

def _coordinates(index):
    nodes=sorted({key[0] for key in index if key[1]=="coordinate"},key=_key)
    coordinates={}
    for node in nodes:
        x,xs=_usable(index.get((node,"coordinate","x")),"m")
        y,ys=_usable(index.get((node,"coordinate","y")),"m")
        if _first_invalid((xs,ys)):
            raise ValueError("independent replay requires canonical coordinates")
        coordinates[node]=(x,y)
    return nodes,coordinates

def _field_evidence(rows,mesh,discovery):
    quantities={"coordinate","ElectrostaticPotential","ElectricField"}
    evidence=[]
    for state,state_rows in sorted(_groups(rows,quantities).items()):
        index=_index(state_rows)
        nodes,coordinates=_coordinates(index)
        cells,edges=_topology(mesh,state[1])
        potential,pstatus,fields,fstatus={},{},{},{}
        for node in nodes:
            potential[node],pstatus[node]=_usable(index.get((node,"ElectrostaticPotential","component0")),"V")
            ex,sx=_usable(index.get((node,"ElectricField","component0")),"V/m")
            ey,sy=_usable(index.get((node,"ElectricField","component1")),"V/m")
            invalid=_first_invalid((sx,sy))
            fields[node],fstatus[node]=((None,invalid) if invalid else ((ex,ey),SampleStatus.VALID))
        gradients,gstatus={},{}
        for cell,cell_nodes in cells.items():
            invalid=_first_invalid(pstatus[node] for node in cell_nodes)
            if invalid:
                gradients[cell],gstatus[cell]=None,invalid
            else:
                gx,gy=_gradient(tuple(coordinates[node] for node in cell_nodes),
                                tuple(potential[node] for node in cell_nodes))
                gradients[cell],gstatus[cell]=(-gx,-gy),SampleStatus.VALID
        split="discovery" if (state[1],state[2]) in discovery else "holdout"
        for cell in sorted(cells,key=_key):
            cell_nodes=cells[cell]
            invalid_ref=_first_invalid(fstatus[node] for node in cell_nodes)
            reference=(None if invalid_ref else tuple(
                sum(fields[node][i] for node in cell_nodes)/3 for i in (0,1)))
            invalid=(gstatus[cell] if gstatus[cell] is not SampleStatus.VALID else invalid_ref)
            evidence.extend(_vector_rows("triangle_minus_grad_psi",state,"cell",cell,None,
                                         gradients[cell],reference,split,field=True,invalid=invalid))
        for node in nodes:
            adjacent=[cell for cell in cells if node in cells[cell]]
            invalid=_first_invalid(gstatus[cell] for cell in adjacent)
            value=None if invalid else _weighted(gradients,cells,coordinates,(node,))
            if invalid is None and fstatus[node] is not SampleStatus.VALID:
                invalid=fstatus[node]
            evidence.extend(_vector_rows("node_area_weighted_minus_grad_psi",state,"node",node,None,
                                         value,fields[node],split,field=True,invalid=invalid))
        for edge in sorted(edges,key=_key):
            start_node,end_node=edges[edge]
            invalid=_first_invalid(gstatus[cell] for cell in cells
                                   if start_node in cells[cell] and end_node in cells[cell])
            value=(None if invalid else _weighted(
                gradients,cells,coordinates,(start_node,end_node)))
            invalid_ref=_first_invalid(fstatus[node] for node in (start_node,end_node))
            reference=(None if invalid_ref else tuple(
                .5*(fields[start_node][i]+fields[end_node][i]) for i in (0,1)))
            evidence.extend(_vector_rows("edge_area_weighted_minus_grad_psi",state,"edge",edge,None,
                                         value,reference,split,field=True,
                                         invalid=invalid or invalid_ref))
            invalid_scalar=_first_invalid(pstatus[node] for node in (start_node,end_node)) or invalid_ref
            candidate=projected=None
            if invalid_scalar is None:
                start,end=coordinates[start_node],coordinates[end_node]
                tx,ty,length=_tangent(start,end)
                signed=-(potential[end_node]-potential[start_node])/length
                scalar=reference[0]*tx+reference[1]*ty
                candidate=(signed*tx,signed*ty); projected=(scalar*tx,scalar*ty)
            evidence.extend(_vector_rows("signed_edge_minus_delta_psi_over_h",state,"edge",edge,None,
                                         candidate,projected,split,field=True,invalid=invalid_scalar))
    return evidence

def _bernoulli(value):
    if abs(value)<1e-8:
        return 1-.5*value+value*value/12
    if value>50:
        return value*math.exp(-value)
    if value < -50:
        return -value
    return value/math.expm1(value)

def _transport_maps(index,nodes,prefix):
    qf,qs,density,ds,mobility,ms,current,js={},{},{},{},{},{},{},{}
    for node in nodes:
        qf[node],qs[node]=_usable(index.get((node,f"{prefix}QuasiFermiPotential","component0")),"V")
        density[node],ds[node]=_usable(index.get((node,f"{prefix}Density","component0")),"m^-3")
        mobility[node],ms[node]=_usable(index.get((node,f"{prefix}Mobility","component0")),"m^2*V^-1*s^-1")
        jx,sx=_usable(index.get((node,f"{prefix}CurrentDensity","component0")),"A/m^2")
        jy,sy=_usable(index.get((node,f"{prefix}CurrentDensity","component1")),"A/m^2")
        invalid=_first_invalid((sx,sy))
        current[node],js[node]=((None,invalid) if invalid else ((jx,jy),SampleStatus.VALID))
        if ds[node] is SampleStatus.VALID and density[node]<=0: ds[node]=SampleStatus.BELOW_FLOOR
        if ms[node] is SampleStatus.VALID and mobility[node]<=0: ms[node]=SampleStatus.BELOW_FLOOR
        if js[node] is SampleStatus.VALID and math.hypot(*current[node])<=0: js[node]=SampleStatus.BELOW_FLOOR
    return qf,qs,density,ds,mobility,ms,current,js

def _transport_evidence(rows,mesh,thermal_voltage,discovery):
    quantities={"coordinate","ElectrostaticPotential","eQuasiFermiPotential",
                "hQuasiFermiPotential","eDensity","hDensity","eMobility","hMobility",
                "eCurrentDensity","hCurrentDensity"}
    evidence=[]
    for state,state_rows in sorted(_groups(rows,quantities).items()):
        index=_index(state_rows); nodes,coordinates=_coordinates(index)
        cells,edges=_topology(mesh,state[1])
        split="discovery" if (state[1],state[2]) in discovery else "holdout"
        psi,psi_status={},{}
        for node in nodes:
            psi[node],psi_status[node]=_usable(index.get((node,"ElectrostaticPotential","component0")),"V")
        for carrier,prefix,sign in (("electron","e",-1.),("hole","h",1.)):
            qf,qs,density,ds,mobility,ms,current,js=_transport_maps(index,nodes,prefix)
            cell_gradients,gs={},{}
            for cell in sorted(cells,key=_key):
                invalid=_first_invalid(qs[node] for node in cells[cell])
                if invalid:
                    cell_gradients[cell],gs[cell]=None,invalid
                else:
                    cell_gradients[cell]=_gradient(
                        tuple(coordinates[node] for node in cells[cell]),
                        tuple(qf[node] for node in cells[cell]))
                    gs[cell]=SampleStatus.VALID
            complete=all(gs[cell] is SampleStatus.VALID for cell in cells)
            node_gradients=({node:_weighted(cell_gradients,cells,coordinates,(node,))
                             for node in nodes} if complete else {})
            edge_gradients=({edge:_weighted(cell_gradients,cells,coordinates,edges[edge])
                             for edge in edges} if complete else {})
            specs=(("triangle_qf_gradient_current","cell"),
                   ("node_area_weighted_qf_gradient_current","node"),
                   ("edge_area_weighted_qf_gradient_current","edge"),
                   ("signed_edge_qf_difference_current","edge"))
            for candidate,support in specs:
                ids=sorted(cells if support=="cell" else nodes if support=="node" else edges,key=_key)
                for support_id in ids:
                    snodes=(cells[support_id] if support=="cell" else
                            (support_id,) if support=="node" else edges[support_id])
                    if support=="cell": gradient=cell_gradients[support_id]
                    elif support=="node": gradient=node_gradients.get(support_id)
                    elif candidate=="signed_edge_qf_difference_current":
                        start_node,end_node=snodes
                        gradient=(None if _first_invalid(qs[node] for node in snodes) else
                                  _edge_vector((qf[end_node]-qf[start_node])/
                                  _tangent(coordinates[start_node],coordinates[end_node])[2],
                                  coordinates[start_node],coordinates[end_node]))
                    else: gradient=edge_gradients.get(support_id)
                    invalid=_first_invalid(tuple(qs[node] for node in snodes)+
                                           tuple(ds[node] for node in snodes)+
                                           tuple(ms[node] for node in snodes))
                    reference_invalid=_first_invalid(js[node] for node in snodes)
                    dvalue=(None if _first_invalid(ds[node] for node in snodes) else
                            statistics.mean(density[node] for node in snodes))
                    mvalue=(None if _first_invalid(ms[node] for node in snodes) else
                            statistics.mean(mobility[node] for node in snodes))
                    reference=(None if any(current[node] is None for node in snodes) else
                               tuple(statistics.mean(current[node][i] for node in snodes) for i in (0,1)))
                    if candidate=="signed_edge_qf_difference_current" and reference is not None:
                        reference=_edge_vector(_project(reference,coordinates[snodes[0]],coordinates[snodes[1]]),
                                               coordinates[snodes[0]],coordinates[snodes[1]])
                    value=(None if invalid is not None or gradient is None else
                           (sign*Q*mvalue*dvalue*gradient[0],sign*Q*mvalue*dvalue*gradient[1]))
                    evidence.extend(_vector_rows(candidate,state,support,support_id,carrier,
                                                 value,reference,split,invalid=invalid,
                                                 reference_invalid=reference_invalid))
                    if invalid is not None:
                        missing=("mobility",) if _first_invalid(ms[node] for node in snodes) else ()
                        evidence.append(_row(candidate,state,support,support_id,carrier,
                                             "transport_abs_dex","current_semantics",
                                             invalid,None,None,split,missing))
            candidate="current_inverted_qf_gradient"
            for node in nodes:
                invalid=_first_invalid((ds[node],ms[node],js[node]))
                reference=node_gradients.get(node)
                value=(None if invalid is not None or reference is None else
                       (current[node][0]/(sign*Q*mobility[node]*density[node]),
                        current[node][1]/(sign*Q*mobility[node]*density[node])))
                evidence.extend(_vector_rows(candidate,state,"node",node,carrier,
                                             value,reference,split,invalid=invalid))
                if invalid is not None:
                    missing=("mobility",) if ms[node] is not SampleStatus.VALID else ()
                    evidence.append(_row(candidate,state,"node",node,carrier,
                                         "transport_abs_dex","current_semantics",
                                         invalid,None,None,split,missing))
            for candidate in ("signed_edge_sg_density_current","signed_edge_drift_diffusion_current"):
                for edge in sorted(edges,key=_key):
                    start_node,end_node=edges[edge]; pair=(start_node,end_node)
                    invalid=_first_invalid(tuple(psi_status[node] for node in pair)+
                                           tuple(ds[node] for node in pair)+
                                           tuple(ms[node] for node in pair))
                    reference_invalid=_first_invalid(js[node] for node in pair)
                    reference=(None if any(current[node] is None for node in pair) else
                               tuple(statistics.mean(current[node][i] for node in pair) for i in (0,1)))
                    start,end=coordinates[start_node],coordinates[end_node]
                    if reference is not None:
                        reference=_edge_vector(_project(reference,start,end),start,end)
                    value=None
                    if invalid is None:
                        length=_tangent(start,end)[2]; mu=statistics.mean(mobility[node] for node in pair)
                        if candidate=="signed_edge_sg_density_current":
                            u=(psi[end_node]-psi[start_node])/thermal_voltage
                            difference=(_bernoulli(u)*density[end_node]-_bernoulli(-u)*density[start_node]
                                        if carrier=="electron" else
                                        _bernoulli(-u)*density[start_node]-_bernoulli(u)*density[end_node])
                            scalar=Q*mu*thermal_voltage*difference/length
                        else:
                            gp=(psi[end_node]-psi[start_node])/length
                            gn=(density[end_node]-density[start_node])/length
                            navg=statistics.mean(density[node] for node in pair)
                            scalar=(Q*mu*(-navg*gp+thermal_voltage*gn) if carrier=="electron"
                                    else Q*mu*(navg*gp-thermal_voltage*gn))
                        value=_edge_vector(scalar,start,end)
                    evidence.extend(_vector_rows(candidate,state,"edge",edge,carrier,value,reference,
                                                 split,invalid=invalid,reference_invalid=reference_invalid))
                    if invalid is not None:
                        missing=("mobility",) if _first_invalid(ms[node] for node in pair) else ()
                        evidence.append(_row(candidate,state,"edge",edge,carrier,
                                             "transport_abs_dex","current_semantics",
                                             invalid,None,None,split,missing))
    return evidence

def _invert_alpha(alpha,parameters,carrier):
    if alpha is None or not math.isfinite(alpha) or alpha<=0:
        return None,SampleStatus.BELOW_FLOOR
    gamma,switch=parameters["gamma"],parameters["switch_field_V_m"]
    found=[]
    for branch in ("low","high"):
        prefactor,critical=parameters[carrier][branch]
        if alpha>=gamma*prefactor: continue
        field=-gamma*critical/math.log(alpha/(gamma*prefactor))
        if ((branch=="low" and field<switch) or (branch=="high" and field>=switch)):
            found.append(field)
    return ((found[0],SampleStatus.VALID) if len(found)==1 else
            (None,SampleStatus.BRANCH_AMBIGUOUS))

def _forward_alpha(driver,parameters,carrier):
    if driver is None: return None,SampleStatus.MISSING_FIELD
    if driver<=0: return None,SampleStatus.BELOW_FLOOR
    gamma,switch=parameters["gamma"],parameters["switch_field_V_m"]
    branch="low" if driver<switch else "high"
    prefactor,critical=parameters[carrier][branch]
    alpha=gamma*prefactor*math.exp(-gamma*critical/driver)
    return ((alpha,SampleStatus.VALID) if alpha>0 and math.isfinite(alpha)
            else (None,SampleStatus.EXPONENTIAL_UNDERFLOW))

def _avalanche_evidence(rows,mesh,parameters,discovery):
    quantities={"coordinate","ElectricField","eQuasiFermiPotential","hQuasiFermiPotential",
                "eCurrentDensity","hCurrentDensity","eDensity","hDensity",
                "eAlphaAvalanche","hAlphaAvalanche","ImpactIonization"}
    evidence=[]
    for state,state_rows in sorted(_groups(rows,quantities).items()):
        index=_index(state_rows); nodes,coordinates=_coordinates(index)
        cells,_=_topology(mesh,state[1])
        split="discovery" if (state[1],state[2]) in discovery else "holdout"
        fields,fs,native,ns={},{},{},{}
        for node in nodes:
            ex,sx=_usable(index.get((node,"ElectricField","component0")),"V/m")
            ey,sy=_usable(index.get((node,"ElectricField","component1")),"V/m")
            invalid=_first_invalid((sx,sy))
            fields[node],fs[node]=((None,invalid) if invalid else ((ex,ey),SampleStatus.VALID))
            native[node],ns[node]=_usable(index.get((node,"ImpactIonization","component0")),"m^-3*s^-1")
        carrier_data={}
        for carrier,prefix in (("electron","e"),("hole","h")):
            qf,qs,current,js,density,ds,reference,rs={},{},{},{},{},{},{},{}
            for node in nodes:
                qf[node],qs[node]=_usable(index.get((node,f"{prefix}QuasiFermiPotential","component0")),"V")
                jx,sx=_usable(index.get((node,f"{prefix}CurrentDensity","component0")),"A/m^2")
                jy,sy=_usable(index.get((node,f"{prefix}CurrentDensity","component1")),"A/m^2")
                invalid=_first_invalid((sx,sy))
                current[node],js[node]=((None,invalid) if invalid else ((jx,jy),SampleStatus.VALID))
                if js[node] is SampleStatus.VALID and math.hypot(*current[node])<=0:
                    current[node],js[node]=None,SampleStatus.BELOW_FLOOR
                density[node],ds[node]=_usable(index.get((node,f"{prefix}Density","component0")),"m^-3")
                alpha,astatus=_usable(index.get((node,f"{prefix}AlphaAvalanche","component0")),"m^-1")
                reference[node],rs[node]=(_invert_alpha(alpha,parameters,carrier)
                                         if astatus is SampleStatus.VALID else (None,astatus))
            complete=all(qs[node] is SampleStatus.VALID for node in nodes)
            if complete:
                cg={cell:_gradient(tuple(coordinates[node] for node in cells[cell]),
                                   tuple(qf[node] for node in cells[cell])) for cell in cells}
                gradients={node:_weighted(cg,cells,coordinates,(node,)) for node in nodes}
                gs={node:SampleStatus.VALID for node in nodes}
            else:
                missing=_first_invalid(qs[node] for node in nodes) or SampleStatus.MISSING_FIELD
                gradients={node:None for node in nodes}; gs={node:missing for node in nodes}
            carrier_data[carrier]={"current":current,"js":js,"density":density,"ds":ds,
                                   "gradient":gradients,"gs":gs,"reference":reference,"rs":rs}
        for candidate in ("electric_field_magnitude","qf_gradient_magnitude",
                          "electric_field_current_aligned","qf_gradient_current_aligned"):
            candidate_generation={}
            for node in nodes:
                alphas={}; invalids=[]
                for carrier in ("electron","hole"):
                    data=carrier_data[carrier]
                    vector,status=((fields[node],fs[node]) if candidate.startswith("electric_field")
                                   else (data["gradient"][node],data["gs"][node]))
                    if candidate.endswith("current_aligned"):
                        if status is SampleStatus.VALID and data["js"][node] is SampleStatus.VALID:
                            current=data["current"][node]; magnitude=math.hypot(*current)
                            driver=abs((vector[0]*current[0]+vector[1]*current[1])/magnitude)
                            dstatus=SampleStatus.VALID
                        else:
                            driver=None; dstatus=status if status is not SampleStatus.VALID else data["js"][node]
                    else:
                        driver=math.hypot(*vector) if status is SampleStatus.VALID else None
                        dstatus=status
                    alphas[carrier],forward=_forward_alpha(driver,parameters,carrier)
                    invalid=_first_invalid((dstatus,forward,data["js"][node],data["rs"][node]))
                    if invalid:
                        invalids.append(invalid); alphas[carrier]=None
                invalid=_first_invalid(invalids); value=None
                if invalid is None:
                    value=(alphas["electron"]*math.hypot(*carrier_data["electron"]["current"][node])+
                           alphas["hole"]*math.hypot(*carrier_data["hole"]["current"][node]))/Q
                    candidate_generation[node]=value
                reference=native[node] if ns[node] is SampleStatus.VALID else None
                if invalid is not None: status,error,prediction=invalid,None,None
                elif value is None or reference is None:
                    status,error,prediction=SampleStatus.MISSING_FIELD,None,None
                elif value<=0 or reference<=0:
                    status,error,prediction=SampleStatus.BELOW_FLOOR,None,None
                else:
                    status=SampleStatus.VALID; prediction=math.log10(value/reference); error=abs(prediction)
                evidence.append(_row(candidate,state,"node",node,None,"local_generation_abs_dex",
                                     "impact_driving_field",status,error,prediction,split))
            status,error,prediction=SampleStatus.MISSING_FIELD,None,None
            if len(candidate_generation)==len(nodes) and all(ns[node] is SampleStatus.VALID for node in nodes):
                candidate_total=native_total=0.
                for cell in sorted(cells,key=_key):
                    cnodes=cells[cell]; area=_area(tuple(coordinates[node] for node in cnodes),signed=True)
                    candidate_total+=area*statistics.mean(candidate_generation[node] for node in cnodes)
                    native_total+=area*statistics.mean(native[node] for node in cnodes)
                if candidate_total>0 and native_total>0:
                    status=SampleStatus.VALID; prediction=math.log10(candidate_total/native_total); error=abs(prediction)
                else: status=SampleStatus.BELOW_FLOOR
            evidence.append(_row(candidate,state,"integrated","integrated",None,
                                 "integrated_generation_abs_dex","impact_driving_field",
                                 status,error,prediction,split))
    return evidence

def _percentile(values,fraction):
    ordered=sorted(values)
    if len(ordered)==1: return ordered[0]
    position=fraction*(len(ordered)-1); lower=math.floor(position); upper=math.ceil(position)
    weight=position-lower
    return ordered[lower]+weight*(ordered[upper]-ordered[lower])

def _summary(rows):
    values=sorted(abs(float(row["error"])) for row in rows
                  if row["status"] is SampleStatus.VALID and row["error"] is not None)
    return {"valid_count":len(values),
            "median_abs_error":statistics.median(values) if values else None,
            "p95_abs_error":_percentile(values,.95) if values else None}

def _limit(metric,thresholds):
    if metric=="field_magnitude_relative": return thresholds.field_median_relative,math.inf
    if metric=="field_direction_deg": return thresholds.field_median_angle_deg,math.inf
    if metric=="transport_abs_dex": return thresholds.gradient_median_abs_dex,thresholds.gradient_p95_abs_dex
    if metric=="transport_direction_deg": return thresholds.gradient_median_angle_deg,math.inf
    if metric=="integrated_generation_abs_dex":
        return thresholds.integrated_generation_abs_dex,thresholds.integrated_generation_abs_dex
    if metric=="local_generation_abs_dex":
        return thresholds.local_generation_abs_dex,thresholds.local_generation_abs_dex
    raise ValueError(f"independent verifier encountered unknown metric {metric}")

def _base(candidate,evidence,thresholds):
    rows=[row for row in evidence if row["candidate"]==candidate]
    if any(row.get("missing_independent_factors") for row in rows):
        return Identifiability.CONFOUNDED
    valid=[row for row in rows if row["status"] is SampleStatus.VALID]
    if (not any(row["split"]=="discovery" for row in valid) or
            not any(row["split"]=="holdout" for row in valid)):
        return Identifiability.INSUFFICIENT_DATA
    metrics=sorted({row["metric"] for row in valid})
    if not metrics: return Identifiability.INSUFFICIENT_DATA
    for metric in metrics:
        combined=[row for row in rows if row["metric"]==metric]
        holdout=[row for row in combined if row["split"]=="holdout"]
        for summary in (_summary(combined),_summary(holdout)):
            if not summary["valid_count"]: return Identifiability.INSUFFICIENT_DATA
            median_limit,p95_limit=_limit(metric,thresholds)
            if (summary["median_abs_error"]>median_limit or
                    summary["p95_abs_error"]>p95_limit):
                return Identifiability.REJECTED
    return Identifiability.IDENTIFIED

def _identity(row):
    return (row["factor"],row["metric"],row["solver"],row["split"],row["topology"],
            row["bias_V"],row["support_kind"],str(row["support_id"]),
            str(row.get("carrier")),"")

def _indistinguishable(first,second,evidence):
    maps={}
    for candidate in (first,second):
        mapping={}
        for row in (row for row in evidence if row["candidate"]==candidate
                    and row["status"] is SampleStatus.VALID):
            if row["prediction_dex"] is None: return False
            identity=_identity(row)
            if identity in mapping:
                raise ValueError("independent replay found duplicate candidate identity")
            mapping[identity]=float(row["prediction_dex"])
        maps[candidate]=mapping
    return (bool(maps[first]) and set(maps[first])==set(maps[second]) and
            all(abs(maps[first][key]-maps[second][key])<=1e-10 for key in maps[first]))

def _rank(evidence,thresholds):
    candidates=sorted({row["candidate"] for row in evidence})
    bases={candidate:_base(candidate,evidence,thresholds) for candidate in candidates}
    classifications={}
    for candidate in candidates:
        classification=bases[candidate]
        if classification is Identifiability.IDENTIFIED and any(
            bases[peer] is Identifiability.IDENTIFIED and
            _indistinguishable(candidate,peer,evidence)
            for peer in candidates if peer!=candidate):
            classification=Identifiability.CONSISTENT_NONUNIQUE
        classifications[candidate]=classification
    ranked=[]
    for candidate in candidates:
        rows=[row for row in evidence if row["candidate"]==candidate]
        valid=[row for row in rows if row["split"]=="discovery"
               and row["status"] is SampleStatus.VALID]
        score=(statistics.mean(abs(float(row["error"])) for row in valid) if valid else math.inf)
        classification=classifications[candidate]
        if classification is Identifiability.IDENTIFIED:
            reason="all combined and holdout gates pass"
        elif classification is Identifiability.CONSISTENT_NONUNIQUE:
            reason="passing prediction is indistinguishable within 1e-10 dex"
        elif classification is Identifiability.CONFOUNDED:
            missing=sorted({factor for row in rows for factor in row.get("missing_independent_factors",())})
            reason=f"independent factors are absent: {', '.join(missing)}"
        elif classification is Identifiability.INSUFFICIENT_DATA:
            reason="discovery or holdout lacks valid typed support"
        else: reason="combined or holdout acceptance gate failed"
        ranked.append((score,candidate,classification,reason,rows))
    return sorted(ranked,key=lambda item:(item[0],item[1]))

def recompute_science(rows,mesh,parameters,thermal_voltage,discovery_keys,thresholds=None):
    """Recompute full authoritative metrics, classifications, and replacement."""
    limits=thresholds or AcceptanceThresholds(); discovery=set(discovery_keys)
    evidence=[]
    evidence.extend(_field_evidence(rows,mesh,discovery))
    evidence.extend(_transport_evidence(rows,mesh,thermal_voltage,discovery))
    evidence.extend(_avalanche_evidence(rows,mesh,parameters,discovery))
    metrics=[]; classifications=[]
    for _,candidate,classification,reason,candidate_rows in _rank(evidence,limits):
        for split in ("discovery","holdout","combined"):
            split_rows=(candidate_rows if split=="combined" else
                        [row for row in candidate_rows if row["split"]==split])
            for metric in sorted({row["metric"] for row in split_rows}):
                summary=_summary([row for row in split_rows if row["metric"]==metric])
                direction=metric.endswith("_direction_deg")
                metrics.append({"candidate":candidate,"quantity":metric,"carrier":"both",
                    "split":split,"topology":"all","bias_V":None,"support_kind":"mixed",
                    "valid_count":summary["valid_count"],
                    "median_abs_error":None if direction else summary["median_abs_error"],
                    "p95_abs_error":None if direction else summary["p95_abs_error"],
                    "median_angle_deg":summary["median_abs_error"] if direction else None,
                    "classification":classification.value})
        classifications.append({"candidate":candidate,"classification":classification.value,
                                "claim_type":"identifiability","reason":reason})
    candidates=sorted({row["candidate"] for row in evidence}); observed={}
    for factor in DEPENDENCIES:
        values=[float(row["prediction_dex"]) for row in evidence
                if row["factor"]==factor and row["status"] is SampleStatus.VALID
                and row["prediction_dex"] is not None and math.isfinite(float(row["prediction_dex"]))]
        if values: observed[factor]=statistics.median(values)
    unavailable=next((factor for factor in DEPENDENCIES if factor not in observed),None)
    if unavailable is None:
        raise ValueError("independent replay unexpectedly has a complete replacement target")
    replacement={"status":SampleStatus.MISSING_FIELD.value,"unavailable_factor":unavailable,
        "dependency_order":list(DEPENDENCIES),"baseline":None,"one_factor":[],"forward":[],
        "reverse":[],"full_replacement":None,"adjacent_interactions":[],"closure":None,
        "evidence_source":"typed_candidate_evidence","evidence_candidates":candidates,
        "observed_prediction_dex":observed}
    return metrics,classifications,replacement
