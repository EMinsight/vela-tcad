#!/usr/bin/env python3
"""Independent PN2D minimal6 fixed-state operator audit."""
from __future__ import annotations
import argparse,csv,hashlib,json,math,platform
from pathlib import Path
from types import SimpleNamespace
import numpy as np
SCHEMA="vela.pn2d_minimal6_fixed_state_audit.v1"; BIASES=(0.0,-12.0,-19.0); FORMULA_LIMIT=5e-12; STATE_LIMIT=1e-12; VT=.025852
NODES={1:(0.,.5),2:(1.,.5),3:(2.,.5),4:(2.,0.),5:(0.,0.),6:(1.,0.)}
TRIS={"sketch":((1,5,2),(5,6,2),(2,6,4),(2,4,3)),"mirror":((1,5,6),(1,6,2),(2,6,3),(6,4,3))}; CONTACTS={"Anode":(1,5),"Cathode":(3,4)}
FIELDS={"ElectrostaticPotential":(1,"V"),"eQuasiFermiPotential":(1,"V"),"hQuasiFermiPotential":(1,"V"),"eDensity":(1,"m^-3"),"hDensity":(1,"m^-3"),"DonorConcentration":(1,"m^-3"),"AcceptorConcentration":(1,"m^-3"),"eMobility":(1,"m^2/(V*s)"),"hMobility":(1,"m^2/(V*s)"),"ElectricField":(2,"V/m"),"eCurrentDensity":(2,"A/m^2"),"hCurrentDensity":(2,"A/m^2"),"eAlphaAvalanche":(1,"m^-1"),"hAlphaAvalanche":(1,"m^-1")}
STATE=("node_id","x_um","y_um","psi_V","phin_V","phip_V","n_m3","p_m3","donors_m3","acceptors_m3","electron_mobility_m2_per_V_s","hole_mobility_m2_per_V_s","electric_field_x_V_per_m","electric_field_y_V_per_m","electron_current_x_A_per_m2","electron_current_y_A_per_m2","hole_current_x_A_per_m2","hole_current_y_A_per_m2","electron_alpha_per_m","hole_alpha_per_m")
VN=("node_id","psi_V","phin_V","phip_V","n_m3","p_m3")
VE=("edge_id","node0","node1","length_m","electron_raw_signed_flux_per_m2_s","hole_raw_signed_flux_per_m2_s","electron_midpoint_density_m3","hole_midpoint_density_m3","electron_impact_field_V_per_m","hole_impact_field_V_per_m","electron_alpha_per_m","hole_alpha_per_m","edge_area_m2")
BT=("cell_id","node0","node1","node2","signed_double_area_m2","grad_psi_x_V_per_m","grad_psi_y_V_per_m","grad_phin_x_V_per_m","grad_phin_y_V_per_m","grad_phip_x_V_per_m","grad_phip_y_V_per_m")
LS=("edge_id","node0","node1","truncated_partial_volume_m2","electron_cell_qf_field_V_per_m","hole_cell_qf_field_V_per_m","electron_midpoint_density_m3","hole_midpoint_density_m3","electron_mobility_m2_per_V_s","hole_mobility_m2_per_V_s","electron_alpha_per_m","hole_alpha_per_m","electron_flux_proxy_per_m2_s","hole_flux_proxy_per_m2_s","electron_source_integral_per_m_s","hole_source_integral_per_m_s")
VTRI=BT+tuple(f"local_edge{i}_{s}" for i in range(3) for s in LS)
class ContractError(RuntimeError): pass
def hybrid_error(a,e,abs_floor=1e-300): return abs(a-e)/max(abs(a),abs(e),abs_floor)
def classify_orientation_pair(s,m):
 d={"mirror_over_sketch":None,"signed_difference":m-s,"absolute_log10_ratio":None}
 if s==0 and m==0: d["zero_classification"]="both_zero"
 elif s==0: d["zero_classification"]="sketch_zero"
 elif m==0: d.update(mirror_over_sketch=0.,zero_classification="mirror_zero")
 else:
  r=m/s; d.update(mirror_over_sketch=r,absolute_log10_ratio=math.log10(abs(r)),zero_classification="neither_zero")
 return d
def bernoulli(x):
 if not math.isfinite(x): raise ContractError("non-finite Bernoulli argument")
 if abs(x)<1e-8: return 1-x/2+x*x/12-x**4/720
 if x>50:
  e=math.exp(-x) if x<745 else 0.; return x*e/(1-e) if e else 0.
 if x<-50:
  e=math.exp(x); return -x/(1-e)
 return x/math.expm1(x)
def triangle_gradient(points,values):
 a=np.array([[1.,float(x),float(y)] for x,y in points]);
 if abs(np.linalg.det(a))<=1e-300: raise ContractError("degenerate triangle")
 c=np.linalg.inv(a)@np.array(values); return float(c[1]),float(c[2])
def area2(p): return (p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[2][0]-p[0][0])*(p[1][1]-p[0][1])
def sg_electron_flux(n0,n1,dpsi,vt,mu,h):
 u=dpsi/vt; return -mu*vt/h*(bernoulli(-u)*n0-bernoulli(u)*n1)
def sg_hole_flux(p0,p1,dpsi,vt,mu,h):
 u=dpsi/vt; return -mu*vt/h*(bernoulli(u)*p0-bernoulli(-u)*p1)
def aux2(x):
 if x>=0:
  e=math.exp(-x) if x<745 else 0.; return e/(1+e)
 e=math.exp(x); return 1/(1+e)
def gss_logistic_midpoint(d0,d1,v0,v1,vt,carrier):
 if carrier=="electron": x=(v1-v0)/(2*vt)
 elif carrier=="hole": x=(v0-v1)/(2*vt)
 else: raise ContractError("carrier must be electron or hole")
 return d0*aux2(x)+d1*aux2(-x)
def canonical_projection(v,p0,p1):
 dx=p1[0]-p0[0];dy=p1[1]-p0[1];h=math.hypot(dx,dy)
 if h<=1e-300: raise ContractError("zero-length projection edge")
 return (v[0]*dx+v[1]*dy)/h
def genius_truncated_partial_volume(points,local_edge):
 p=[np.array(x,float) for x in points]; d=area2(p)
 if abs(d)<=1e-300:return 0.
 s=[float(x@x) for x in p]; c=np.array([(s[0]*(p[1][1]-p[2][1])+s[1]*(p[2][1]-p[0][1])+s[2]*(p[0][1]-p[1][1]))/(2*d),(s[0]*(p[2][0]-p[1][0])+s[1]*(p[0][0]-p[2][0])+s[2]*(p[1][0]-p[0][0]))/(2*d)])
 sides=((0,1),(1,2),(2,0)); lengths=[];dt=[];obt=-1
 for k,(i,j) in enumerate(sides):
  z=(2+k)%3; lengths.append(float(np.linalg.norm(p[i]-p[j]))); dist=float(np.linalg.norm((p[i]+p[j])/2-c))
  if float((p[i]-p[z])@(p[j]-p[z]))<0: dt.append(-dist);obt=k
  else:dt.append(dist)
 if obt>=0:
  i,j=sides[obt];z=(2+obt)%3;a,b,q=p[i],p[j],p[z]
  def ang(x,y):
   den=np.linalg.norm(x)*np.linalg.norm(y);return 0 if den<=1e-300 else math.acos(max(-1,min(1,float(x@y)/den)))
  c1=math.cos(ang(b-a,q-a));c2=math.cos(ang(a-b,q-b));dt[obt]=0
  if abs(c1)>1e-300:
   mid=(a+q)/2;m=a+(b-a)/np.linalg.norm(b-a)*(np.linalg.norm(mid-a)/c1);dt[(obt+2)%3]=float(np.linalg.norm(mid-m))
  if abs(c2)>1e-300:
   mid=(b+q)/2;m=b+(a-b)/np.linalg.norm(a-b)*(np.linalg.norm(mid-b)/c2);dt[(obt+1)%3]=float(np.linalg.norm(mid-m))
 return .5*lengths[local_edge]*max(0.,dt[local_edge])
def require_unique_keys(rows,keys):
 seen=set()
 for r in rows:
  k=tuple(r[x] for x in keys)
  if k in seen:raise ContractError(f"duplicate key {k}")
  seen.add(k)
def finite(v,label):
 try:x=float(v)
 except:raise ContractError(f"invalid numeric value for {label}")
 if not math.isfinite(x):raise ContractError(f"non-finite value for {label}")
 return x
def readcsv(path,header):
 with Path(path).open(encoding="utf-8",newline="") as f:
  rd=csv.DictReader(f)
  if tuple(rd.fieldnames or ())!=tuple(header):raise ContractError(f"wrong CSV schema in {Path(path).name}")
  rows=list(rd)
 for r in rows:
  for k,v in r.items():finite(v,f"{Path(path).name}:{k}")
 return rows
def edges(tris):return sorted({tuple(sorted((t[i],t[(i+1)%3]))) for t in tris for i in range(3)})
def gate(a,e,limit,label):
 x=hybrid_error(a,e)
 if x>=limit:raise ContractError(f"{label} {x:.17g} >= {limit:.17g}")
 return x
def validate_topology(data,tid):
 if data.get("topology_id")!=tid or tid not in TRIS:raise ContractError("wrong topology id")
 mapped={}
 for r in data.get("nodes",[]):
  xy=(finite(r.get("x_um"),"x_um"),finite(r.get("y_um"),"y_um"));m=[n for n,p in NODES.items() if abs(xy[0]-p[0])<1e-12 and abs(xy[1]-p[1])<1e-12]
  if len(m)!=1 or m[0] in mapped:raise ContractError("topology coordinate mapping is missing or duplicate")
  mapped[m[0]]=r
 tris=tuple(tuple(map(int,t)) for t in data.get("triangles",[]))
 if len(mapped)!=6 or tris!=TRIS[tid] or len(edges(tris))!=9:raise ContractError("wrong topology connectivity or counts")
 if data.get("contacts")!={k:list(v) for k,v in CONTACTS.items()}:raise ContractError("wrong topology contact ownership")
 for n in (2,6):
  if finite(mapped[n].get("donors_cm3"),"donors")!=1e17 or finite(mapped[n].get("acceptors_cm3"),"acceptors")!=1e17:raise ContractError("wrong topology compensated semantics")
 return tris
def validate_fields(path):
 d=json.loads(Path(path).read_text(encoding="utf-8"));fs=d.get("fields",{})
 for name,(comp,unit) in FIELDS.items():
  if name not in fs:raise ContractError(f"missing required field {name}")
  x=fs[name]
  if x.get("components")!=comp or x.get("unit")!=unit:raise ContractError(f"wrong unit or component count for {name}")
  if x.get("region")!="Silicon" or x.get("mapping_status")!="mapped":raise ContractError(f"incomplete field mapping for {name}")
 if d.get("global_node_mapping")!="exact_coordinate_1e-12_um":raise ContractError("incomplete field manifest node mapping")
def mapstate(rows):
 out={}
 for r in rows:
  xy=(finite(r["x_um"],"x"),finite(r["y_um"],"y"));m=[n for n,p in NODES.items() if abs(xy[0]-p[0])<1e-12 and abs(xy[1]-p[1])<1e-12]
  if len(m)!=1 or m[0] in out:raise ContractError("state coordinate mapping is missing or duplicate")
  out[m[0]]=r
 if len(out)!=6:raise ContractError("partial state matrix")
 return out
def annotate_reconstructed(rows):
 for tid,ts in TRIS.items():
  for bias in BIASES:
   group=[r for r in rows if r["topology_id"]==tid and r["bias_V"]==bias];lookup={(r["node0"],r["node1"]):r for r in group}
   for carrier in ("electron","hole"):
    vectors=[]
    for t in ts:
     matrix=[];values=[]
     for i in range(3):
      key=tuple(sorted((t[i],t[(i+1)%3])));e=lookup[key];matrix.append([e["dx_m"]/e["length_m"],e["dy_m"]/e["length_m"]]);values.append(e[f"vela_{carrier}_flux_per_m2_s"])
     vectors.append((t,np.linalg.lstsq(np.array(matrix),np.array(values),rcond=None)[0]))
    for edge,r in lookup.items():
     adjacent=[v for t,v in vectors if edge[0] in t and edge[1] in t];v=np.mean(adjacent,axis=0);r[f"vela_{carrier}_reconstructed_projection_per_m2_s"]=canonical_projection(v,(0.,0.),(r["dx_m"],r["dy_m"]));r[f"vela_{carrier}_reconstructed_magnitude_per_m2_s"]=float(np.linalg.norm(v))
def build_report(fixture):
 root=Path(fixture);man=json.loads((root/"manifest.json").read_text(encoding="utf-8"))
 if man.get("fixture_schema")!="vela.pn2d_minimal6_synthetic_fixture.v1" or len(man.get("topologies",[]))!=2:raise ContractError("wrong fixture topology matrix")
 nr=[];er=[];tr=[];seen=set()
 for top in man["topologies"]:
  tid=top.get("topology_id")
  if tid in seen:raise ContractError("duplicate topology id")
  seen.add(tid);tris=validate_topology(json.loads((root/top["topology_file"]).read_text(encoding="utf-8")),tid);states=top.get("states",[])
  if len(states)!=3 or tuple(x.get("bias_V") for x in states)!=BIASES:raise ContractError("states must contain exact biases 0, -12, -19 V")
  for ss in states:
   bias=float(ss["bias_V"]);validate_fields(root/ss["field_manifest"]);sent=mapstate(readcsv(root/ss["sentaurus_state_csv"],STATE));vn=readcsv(root/ss["vela_node_csv"],VN)
   if len(vn)!=6:raise ContractError("partial state matrix")
   vm={int(x["node_id"])+1:x for x in vn}
   if set(vm)!=set(NODES):raise ContractError("missing or duplicate Vela node key")
   for n in NODES:
    s=sent[n];v=vm[n];r={"topology_id":tid,"bias_V":bias,"node_id":n,"x_um":NODES[n][0],"y_um":NODES[n][1]}
    for name in ("psi_V","phin_V","phip_V","n_m3","p_m3"):
     a=finite(v[name],name);e=finite(s[name],name);r[f"sentaurus_{name}"]=e;r[f"vela_{name}"]=a;r[f"state_error_{name}"]=gate(a,e,STATE_LIMIT,"state parity error")
    for name in STATE[8:]:r[f"sentaurus_{name}"]=finite(s[name],name)
    nr.append(r)
   vr=readcsv(root/ss["vela_edge_csv"],VE)

   ix={}
   for x in vr:
    k=tuple(sorted((int(x["node0"])+1,int(x["node1"])+1)))
    if k in ix:raise ContractError("duplicate Vela edge key")
    ix[k]=x
   if len(vr)!=9:raise ContractError("edge row count mismatch")
   if set(ix)!=set(edges(tris)):raise ContractError("wrong topology edge matrix")
   for eid,(a,b) in enumerate(edges(tris)):
    s0,s1=sent[a],sent[b];p0=np.array(NODES[a])*1e-6;p1=np.array(NODES[b])*1e-6;dv=p1-p0;h=float(np.linalg.norm(dv));raw=ix[(a,b)];dpsi=float(s1["psi_V"])-float(s0["psi_V"])
    n0,n1=float(s0["n_m3"]),float(s1["n_m3"]);pden0,pden1=float(s0["p_m3"]),float(s1["p_m3"]);mun=(float(s0["electron_mobility_m2_per_V_s"])+float(s1["electron_mobility_m2_per_V_s"]))/2;mup=(float(s0["hole_mobility_m2_per_V_s"])+float(s1["hole_mobility_m2_per_V_s"]))/2
    pe=sg_electron_flux(n0,n1,dpsi,VT,mun,h);ph=sg_hole_flux(pden0,pden1,dpsi,VT,mup,h);me=gss_logistic_midpoint(n0,n1,float(s0["psi_V"]),float(s1["psi_V"]),VT,"electron");mh=gss_logistic_midpoint(pden0,pden1,float(s0["psi_V"]),float(s1["psi_V"]),VT,"hole")
    ve=float(raw["electron_raw_signed_flux_per_m2_s"]);vh=float(raw["hole_raw_signed_flux_per_m2_s"]);errs=[gate(float(raw["length_m"]),h,FORMULA_LIMIT,"formula error"),gate(ve,pe,FORMULA_LIMIT,"formula error"),gate(vh,ph,FORMULA_LIMIT,"formula error"),gate(float(raw["electron_midpoint_density_m3"]),me,FORMULA_LIMIT,"formula error"),gate(float(raw["hole_midpoint_density_m3"]),mh,FORMULA_LIMIT,"formula error")]
    je=[(float(s0["electron_current_x_A_per_m2"])+float(s1["electron_current_x_A_per_m2"]))/2,(float(s0["electron_current_y_A_per_m2"])+float(s1["electron_current_y_A_per_m2"]))/2];jh=[(float(s0["hole_current_x_A_per_m2"])+float(s1["hole_current_x_A_per_m2"]))/2,(float(s0["hole_current_y_A_per_m2"])+float(s1["hole_current_y_A_per_m2"]))/2];adj=sum(a in t and b in t for t in tris);cls="contact" if (a,b) in CONTACTS.values() else ("interior" if adj==2 else "boundary")
    er.append({"topology_id":tid,"bias_V":bias,"edge_id":eid,"node0":a,"node1":b,"edge_class":cls,"dx_m":float(dv[0]),"dy_m":float(dv[1]),"length_m":h,"delta_phin_over_h_V_per_m":(float(s1["phin_V"])-float(s0["phin_V"]))/h,"delta_phip_over_h_V_per_m":(float(s1["phip_V"])-float(s0["phip_V"]))/h,"bernoulli_argument":dpsi/VT,"bernoulli_plus":bernoulli(dpsi/VT),"bernoulli_minus":bernoulli(-dpsi/VT),"python_electron_flux_per_m2_s":pe,"python_hole_flux_per_m2_s":ph,"vela_electron_flux_per_m2_s":ve,"vela_hole_flux_per_m2_s":vh,"sentaurus_electron_projection_A_per_m2":canonical_projection(je,p0,p1),"sentaurus_hole_projection_A_per_m2":canonical_projection(jh,p0,p1),"sentaurus_electron_magnitude_A_per_m2":math.hypot(*je),"sentaurus_hole_magnitude_A_per_m2":math.hypot(*jh),"electron_midpoint_density_m3":me,"hole_midpoint_density_m3":mh,"electron_alpha_per_m":float(raw["electron_alpha_per_m"]),"hole_alpha_per_m":float(raw["hole_alpha_per_m"]),"edge_area_m2":float(raw["edge_area_m2"]),"electron_edge_source_integral_per_s":float(raw["electron_alpha_per_m"])*abs(ve)*float(raw["edge_area_m2"]),"hole_edge_source_integral_per_s":float(raw["hole_alpha_per_m"])*abs(vh)*float(raw["edge_area_m2"]),"sentaurus_vs_vela_electron_current_diagnostic_error":hybrid_error(canonical_projection(je,p0,p1),1.602176634e-19*ve),"sentaurus_vs_vela_hole_current_diagnostic_error":hybrid_error(canonical_projection(jh,p0,p1),1.602176634e-19*vh),"max_formula_error":max(errs)})
   vt=readcsv(root/ss["vela_triangle_csv"],VTRI)
   if len(vt)!=4:raise ContractError("triangle row count mismatch")
   ti={tuple(int(x[f"node{i}"])+1 for i in range(3)):x for x in vt}
   if len(ti)!=4 or set(ti)!=set(tris):raise ContractError("duplicate or wrong topology triangle key")
   for cid,t in enumerate(tris):
    raw=ti[t];pts=[tuple(np.array(NODES[n])*1e-6) for n in t];ar=area2(pts)
    if ar<=0:raise ContractError("reversed triangle orientation")
    gs={name:triangle_gradient(pts,[float(sent[n][name]) for n in t]) for name in ("psi_V","phin_V","phip_V")};es=[gate(float(raw["signed_double_area_m2"]),ar,FORMULA_LIMIT,"formula error")]
    for name,prefix in (("psi_V","grad_psi"),("phin_V","grad_phin"),("phip_V","grad_phip")):
     es+=[gate(float(raw[f"{prefix}_x_V_per_m"]),gs[name][0],FORMULA_LIMIT,"formula error"),gate(float(raw[f"{prefix}_y_V_per_m"]),gs[name][1],FORMULA_LIMIT,"formula error")]
    se=sh=0.;pv=[];en={n:0. for n in NODES};hn={n:0. for n in NODES}
    for local in range(3):
     x=genius_truncated_partial_volume(pts,local);pre=f"local_edge{local}_";a,b=t[local],t[(local+1)%3];s0,s1=sent[a],sent[b];h=math.dist(pts[local],pts[(local+1)%3]);eef=abs(float(s1["phin_V"])-float(s0["phin_V"]))/h;hef=abs(float(s1["phip_V"])-float(s0["phip_V"]))/h;em=gss_logistic_midpoint(float(s0["n_m3"]),float(s1["n_m3"]),float(s0["psi_V"]),float(s1["psi_V"]),VT,"electron");hm=gss_logistic_midpoint(float(s0["p_m3"]),float(s1["p_m3"]),float(s0["psi_V"]),float(s1["psi_V"]),VT,"hole");mun=(float(s0["electron_mobility_m2_per_V_s"])+float(s1["electron_mobility_m2_per_V_s"]))/2;mup=(float(s0["hole_mobility_m2_per_V_s"])+float(s1["hole_mobility_m2_per_V_s"]))/2;ep=mun*em*eef;hp=mup*hm*hef;esi=float(raw[pre+"electron_alpha_per_m"])*ep*x;hsi=float(raw[pre+"hole_alpha_per_m"])*hp*x
     es += [gate(float(raw[pre+"truncated_partial_volume_m2"]),x,FORMULA_LIMIT,"formula error"),gate(float(raw[pre+"electron_cell_qf_field_V_per_m"]),math.hypot(*gs["phin_V"]),FORMULA_LIMIT,"formula error"),gate(float(raw[pre+"hole_cell_qf_field_V_per_m"]),math.hypot(*gs["phip_V"]),FORMULA_LIMIT,"formula error"),gate(float(raw[pre+"electron_midpoint_density_m3"]),em,FORMULA_LIMIT,"formula error"),gate(float(raw[pre+"hole_midpoint_density_m3"]),hm,FORMULA_LIMIT,"formula error"),gate(float(raw[pre+"electron_flux_proxy_per_m2_s"]),ep,FORMULA_LIMIT,"formula error"),gate(float(raw[pre+"hole_flux_proxy_per_m2_s"]),hp,FORMULA_LIMIT,"formula error"),gate(float(raw[pre+"electron_source_integral_per_m_s"]),esi,FORMULA_LIMIT,"formula error"),gate(float(raw[pre+"hole_source_integral_per_m_s"]),hsi,FORMULA_LIMIT,"formula error")]
     pv.append(x);se+=esi;sh+=hsi;en[a]+=esi/2;en[b]+=esi/2;hn[a]+=hsi/2;hn[b]+=hsi/2
    rr={"topology_id":tid,"bias_V":bias,"cell_id":cid,"node0":t[0],"node1":t[1],"node2":t[2],"signed_double_area_m2":ar,"grad_psi_x_V_per_m":gs["psi_V"][0],"grad_psi_y_V_per_m":gs["psi_V"][1],"grad_phin_x_V_per_m":gs["phin_V"][0],"grad_phin_y_V_per_m":gs["phin_V"][1],"grad_phip_x_V_per_m":gs["phip_V"][0],"grad_phip_y_V_per_m":gs["phip_V"][1],"electron_source_integral_per_m_s":se,"hole_source_integral_per_m_s":sh,"total_source_integral_per_m_s":se+sh,"partial_volumes_m2":json.dumps(pv),"max_formula_error":max(es)};rr.update({f"vela_{k}":finite(v,k) for k,v in raw.items()});rr.update({f"electron_node{n}_source_partition_per_m_s":en[n] for n in NODES});rr.update({f"hole_node{n}_source_partition_per_m_s":hn[n] for n in NODES});tr.append(rr)
 annotate_reconstructed(er)
 require_unique_keys(nr,("topology_id","bias_V","node_id"));require_unique_keys(er,("topology_id","bias_V","node0","node1"));require_unique_keys(tr,("topology_id","bias_V","node0","node1","node2"))
 if (len(nr),len(er),len(tr))!=(36,54,24):raise ContractError("output row-count mismatch")
 orient=[]
 for b in BIASES:
  for car in ("electron","hole"):
   key=f"vela_{car}_flux_per_m2_s"
   for edge in sorted({(r["node0"],r["node1"]) for r in er}):
    d={r["topology_id"]:r[key] for r in er if r["bias_V"]==b and (r["node0"],r["node1"])==edge}
    if set(d)==set(TRIS):orient.append({"bias_V":b,"quantity":f"edge_{edge[0]}_{edge[1]}_{car}_flux",**classify_orientation_pair(d["sketch"],d["mirror"])})
  for car in ("electron","hole","total"):
   key=f"{car}_source_integral_per_m_s";d={tid:sum(r[key] for r in tr if r["bias_V"]==b and r["topology_id"]==tid) for tid in TRIS};orient.append({"bias_V":b,"quantity":f"integrated_{car}_avalanche_source",**classify_orientation_pair(d["sketch"],d["mirror"])})
 summary={"schema":SCHEMA,"status":"PASS","scope":"fixed-state operator audit, not a BV curve","row_counts":{"node_state":36,"edge_audit":54,"triangle_audit":24},"gates":{"passed":True,"state_parity_limit":STATE_LIMIT,"formula_limit":FORMULA_LIMIT,"sentaurus_vs_vela_current_source_threshold":None},"orientation_sensitivity":orient,"chart_map":["minimal6-topologies","minimal6-edge-current-audit-{0v,minus12v,minus19v}","minimal6-triangle-source-audit-{0v,minus12v,minus19v}"],"qa_notes":["Three biases are separate diagnostic samples.","Zero pairs use both_zero without log ratio.","Blue/orange and marker/style redundancy."]}
 return SimpleNamespace(node_rows=nr,edge_rows=er,triangle_rows=tr,summary=summary,manifest=man,fixture_root=root)
def writecsv(path,rows):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 with path.open("w",encoding="utf-8",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def plot_reports(rep,out):
 import matplotlib;matplotlib.use("Agg")
 import matplotlib.pyplot as plt
 fd=Path(out)/"figures";fd.mkdir(parents=True,exist_ok=True);blue="#2166ac";orange="#e08214"
 fig,axs=plt.subplots(1,2,figsize=(13,4.8),constrained_layout=True)
 for ax,(tid,ts) in zip(axs,TRIS.items()):
  for t in ts:
   p=np.array([NODES[n] for n in t+(t[0],)]);ax.plot(p[:,0],p[:,1],color="#555555");c=p[:-1].mean(0);ax.text(c[0],c[1],"T"+"-".join(map(str,t)),fontsize=8,ha="center",bbox={"facecolor":"white","alpha":.8,"edgecolor":"none"})
  for a,b in edges(ts):
   m=(np.array(NODES[a])+np.array(NODES[b]))/2;ax.text(m[0],m[1]+.025,f"E{a}-{b}",fontsize=7,ha="center")
  for n,(x,y) in NODES.items():ax.scatter(x,y,s=65,color=blue,edgecolor="black",zorder=3);ax.text(x,y-.055,f"N{n}",ha="center",fontsize=9,fontweight="bold")
  for n in (2,6):x,y=NODES[n];ax.scatter(x,y,s=155,facecolor="none",edgecolor=orange,lw=2,zorder=4)
  ax.plot([0,0],[0,.5],color=orange,lw=4,label="Anode contact");ax.plot([2,2],[0,.5],color=blue,lw=4,ls="--",label="Cathode contact")
  ax.set(title=f"{tid}: all nodes, edges, triangles, contacts\norange rings: compensated nodes 2 and 6",xlabel="x (um)",ylabel="y (um)",xlim=(-.12,2.12),ylim=(-.12,.64),aspect="equal");ax.legend(fontsize=8,ncol=2)
 fig.suptitle("PN2D minimal6 canonical topology audit")
 for ext in ("png","pdf"):fig.savefig(fd/f"minimal6-topologies.{ext}",dpi=180)
 plt.close(fig)
 for bias,slug in ((0.,"0v"),(-12.,"minus12v"),(-19.,"minus19v")):
  rows=[r for r in rep.edge_rows if r["bias_V"]==bias];x=np.arange(len(rows));labs=[f"{r['topology_id'][0].upper()} {r['node0']}-{r['node1']}" for r in rows]
  fig,ax=plt.subplots(figsize=(14,6),constrained_layout=True);ax.plot(x,[r["vela_electron_flux_per_m2_s"] for r in rows],color=blue,marker="o",label="Vela electron particle flux");ax.plot(x,[r["vela_hole_flux_per_m2_s"] for r in rows],color=orange,marker="s",ls="--",label="Vela hole particle flux");ax.axhline(0,color="#555",lw=.8);ax.set_yscale("symlog",linthresh=1e16);ax.set_xticks(x,labs,rotation=55,ha="right",fontsize=8);ax.set(xlabel="topology and canonical edge (discrete sample)",ylabel="signed particle flux (m$^{-2}$ s$^{-1}$)",title=f"PN2D minimal6 edge current at {bias:g} V\nfixed-state operator audit, not a BV curve");ax.legend(ncol=2);ax.grid(alpha=.2)
  for ext in ("png","pdf"):fig.savefig(fd/f"minimal6-edge-current-audit-{slug}.{ext}",dpi=180)
  plt.close(fig);rows=[r for r in rep.triangle_rows if r["bias_V"]==bias];x=np.arange(len(rows));labs=[f"{r['topology_id'][0].upper()} T{r['node0']}-{r['node1']}-{r['node2']}" for r in rows]
  fig,ax=plt.subplots(figsize=(12,5.5),constrained_layout=True);ax.plot(x,[r["electron_source_integral_per_m_s"] for r in rows],color=blue,marker="o",label="electron contribution");ax.plot(x,[r["hole_source_integral_per_m_s"] for r in rows],color=orange,marker="s",ls="--",label="hole contribution");ax.axhline(0,color="#555",lw=.8);ax.ticklabel_format(axis="y",style="sci",scilimits=(0,0));ax.set_xticks(x,labs,rotation=35,ha="right");ax.set(xlabel="topology and canonical triangle (discrete sample)",ylabel="source integral (m$^{-1}$ s$^{-1}$)",title=f"PN2D minimal6 triangle source at {bias:g} V\nfixed-state operator audit, not a BV curve");ax.legend();ax.grid(alpha=.2)
  for ext in ("png","pdf"):fig.savefig(fd/f"minimal6-triangle-source-audit-{slug}.{ext}",dpi=180)
  plt.close(fig)
def write_report(rep,out):
 out=Path(out);out.mkdir(parents=True,exist_ok=True);writecsv(out/"node_state.csv",rep.node_rows);writecsv(out/"edge_audit.csv",rep.edge_rows);writecsv(out/"triangle_audit.csv",rep.triangle_rows);(out/"summary.json").write_text(json.dumps(rep.summary,indent=2)+"\n",encoding="utf-8");files=sorted(p for p in rep.fixture_root.rglob("*") if p.is_file());hashes={p.relative_to(rep.fixture_root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files};output_manifest={"schema":SCHEMA,"topology_definitions":rep.manifest["topologies"],"bias_states_V":list(BIASES),"model_configuration":{"workflow":"immutable_fixed_state","thermal_voltage_V":VT,"solver_runs":False},"tool_versions":{"python":platform.python_version(),"numpy":np.__version__},"command_status":{"fixture_validation":"PASS","report_generation":"PASS"},"input_sha256":hashes,"row_counts":rep.summary["row_counts"],"gate_status":"PASS"};(out/"manifest.json").write_text(json.dumps(output_manifest,indent=2)+"\n",encoding="utf-8")
 md='''# PN2D Minimal6 Fixed-State Operator Audit

**Answer first:** all strict topology, state-parity, independent-formula, completeness, finiteness, and uniqueness gates pass for this synthetic six-state matrix. This is a fixed-state operator audit, not a BV curve, and supports no physical breakdown-voltage conclusion.

## Scope and definitions

Two exact six-node/four-triangle topologies are evaluated independently at the discrete diagnostic samples 0 V, -12 V, and -19 V. Canonical edges are sorted node pairs and canonical triangles are the approved CCW tuples. `both_zero` is reported without a logarithmic ratio.

## Method and gates

Sentaurus-like state/field exports are joined by exact coordinates to immutable Task4-schema Vela CSV outputs. Python independently recomputes inverse-matrix triangle gradients, stable Bernoulli factors, electron/hole SG fluxes, GSS logistic midpoint densities, canonical vector projections, and Genius-truncated partial volumes. State parity must be below 1e-12 and formula error below 5e-12; equality fails closed. Exact row counts are 36/54/24. Sentaurus-versus-Vela current and source differences are diagnostic and have no improvement threshold.

## Limitations

The fixture is synthetic and does not replace a real Sentaurus VM export. The biases are separate fixed-state samples, not a sweep or trend. No Newton, Gummel, continuation, or independent DC solve is performed, and no physical BV claim is permitted.

## Next step

Replay the same contract against the six validated real Sentaurus exports and Task4 production CSVs, retaining all fail-closed gates and reviewing orientation sensitivity before interpreting cross-product differences.
''';(out/"summary.md").write_text(md,encoding="utf-8");plot_reports(rep,out)
def stateval(n,bias):
 x,y=NODES[n];return {"node_id":n-1,"x_um":x,"y_um":y,"psi_V":bias*x/2+.02*y,"phin_V":.004*x+.002*y+.001*bias,"phip_V":-.003*x+.001*y+.0005*bias,"n_m3":1e21*(1+.1*x+.05*y+.01*abs(bias)),"p_m3":8e20*(1+.08*(2-x)+.03*y+.008*abs(bias)),"donors_m3":1e23 if n in (2,3,4,6) else 0.,"acceptors_m3":1e23 if n in (1,2,5,6) else 0.,"electron_mobility_m2_per_V_s":.135,"hole_mobility_m2_per_V_s":.048,"electric_field_x_V_per_m":-bias*.5e6,"electric_field_y_V_per_m":-.02e6,"electron_current_x_A_per_m2":(n+.1*abs(bias))*2e4,"electron_current_y_A_per_m2":n*1e3,"hole_current_x_A_per_m2":-(n+.05*abs(bias))*8e3,"hole_current_y_A_per_m2":n*5e2,"electron_alpha_per_m":1e3+abs(bias)*50,"hole_alpha_per_m":800+abs(bias)*40}
def make_synthetic_fixture(root):
 root=Path(root);root.mkdir(parents=True,exist_ok=True);tops=[]
 for tid,ts in TRIS.items():
  td=root/tid;td.mkdir(exist_ok=True);nodes=[{"node_id":n,"x_um":x,"y_um":y,"donors_cm3":1e17 if n in (2,3,4,6) else 0.,"acceptors_cm3":1e17 if n in (1,2,5,6) else 0.} for n,(x,y) in NODES.items()];(td/"topology.json").write_text(json.dumps({"topology_id":tid,"nodes":nodes,"triangles":[list(x) for x in ts],"contacts":{k:list(v) for k,v in CONTACTS.items()},"region":"Silicon"},indent=2)+"\n",encoding="utf-8");states=[]
  for bias,slug in ((0.,"0v"),(-12.,"minus12v"),(-19.,"minus19v")):
   sd=td/slug;sd.mkdir(exist_ok=True);fs={n:{"region":"Silicon","components":c,"unit":u,"mapping_status":"mapped"} for n,(c,u) in FIELDS.items()};(sd/"field_manifest.json").write_text(json.dumps({"bias_V":bias,"global_node_mapping":"exact_coordinate_1e-12_um","fields":fs},indent=2)+"\n",encoding="utf-8");st=[stateval(n,bias) for n in NODES];writecsv(sd/"sentaurus_state.csv",st);writecsv(sd/"vela_node_state.csv",[{k:r[k] for k in VN} for r in st]);ers=[]
   for eid,(a,b) in enumerate(edges(ts)):
    x,y=st[a-1],st[b-1];p0=np.array(NODES[a])*1e-6;p1=np.array(NODES[b])*1e-6;h=float(np.linalg.norm(p1-p0));dv=y["psi_V"]-x["psi_V"]
    ers.append({"edge_id":eid,"node0":a-1,"node1":b-1,"length_m":h,"electron_raw_signed_flux_per_m2_s":sg_electron_flux(x["n_m3"],y["n_m3"],dv,VT,.135,h),"hole_raw_signed_flux_per_m2_s":sg_hole_flux(x["p_m3"],y["p_m3"],dv,VT,.048,h),"electron_midpoint_density_m3":gss_logistic_midpoint(x["n_m3"],y["n_m3"],x["psi_V"],y["psi_V"],VT,"electron"),"hole_midpoint_density_m3":gss_logistic_midpoint(x["p_m3"],y["p_m3"],x["psi_V"],y["psi_V"],VT,"hole"),"electron_impact_field_V_per_m":abs(y["phin_V"]-x["phin_V"])/h,"hole_impact_field_V_per_m":abs(y["phip_V"]-x["phip_V"])/h,"electron_alpha_per_m":x["electron_alpha_per_m"],"hole_alpha_per_m":x["hole_alpha_per_m"],"edge_area_m2":h*1e-6})
   writecsv(sd/"vela_edge_audit.csv",ers);eids={tuple(sorted((r["node0"]+1,r["node1"]+1))):r["edge_id"] for r in ers};trs=[]
   for cid,t in enumerate(ts):
    pts=[tuple(np.array(NODES[n])*1e-6) for n in t];gs={name:triangle_gradient(pts,[st[n-1][name] for n in t]) for name in ("psi_V","phin_V","phip_V")};r={"cell_id":cid,"node0":t[0]-1,"node1":t[1]-1,"node2":t[2]-1,"signed_double_area_m2":area2(pts),"grad_psi_x_V_per_m":gs["psi_V"][0],"grad_psi_y_V_per_m":gs["psi_V"][1],"grad_phin_x_V_per_m":gs["phin_V"][0],"grad_phin_y_V_per_m":gs["phin_V"][1],"grad_phip_x_V_per_m":gs["phip_V"][0],"grad_phip_y_V_per_m":gs["phip_V"][1]};ef=math.hypot(*gs["phin_V"]);hf=math.hypot(*gs["phip_V"]);ea=1e3+ef*1e-4;ha=800+hf*1e-4
    for i in range(3):
     a,b=t[i],t[(i+1)%3];x,y=st[a-1],st[b-1];h=math.dist(pts[i],pts[(i+1)%3]);eef=abs(y["phin_V"]-x["phin_V"])/h;hef=abs(y["phip_V"]-x["phip_V"])/h;em=gss_logistic_midpoint(x["n_m3"],y["n_m3"],x["psi_V"],y["psi_V"],VT,"electron");hm=gss_logistic_midpoint(x["p_m3"],y["p_m3"],x["psi_V"],y["psi_V"],VT,"hole");pv=genius_truncated_partial_volume(pts,i);ep=.135*em*eef;hp=.048*hm*hef;vals=(eids[tuple(sorted((a,b)))],a-1,b-1,pv,ef,hf,em,hm,.135,.048,ea,ha,ep,hp,ea*ep*pv,ha*hp*pv)
     for s,v in zip(LS,vals):r[f"local_edge{i}_{s}"]=v
    trs.append(r)
   writecsv(sd/"vela_triangle_audit.csv",trs);rel=lambda p:p.relative_to(root).as_posix();states.append({"bias_V":bias,"field_manifest":rel(sd/"field_manifest.json"),"sentaurus_state_csv":rel(sd/"sentaurus_state.csv"),"vela_node_csv":rel(sd/"vela_node_state.csv"),"vela_edge_csv":rel(sd/"vela_edge_audit.csv"),"vela_triangle_csv":rel(sd/"vela_triangle_audit.csv")})
  tops.append({"topology_id":tid,"topology_file":f"{tid}/topology.json","states":states})
 (root/"manifest.json").write_text(json.dumps({"fixture_schema":"vela.pn2d_minimal6_synthetic_fixture.v1","description":"Synthetic Sentaurus-like fields and static Task4-schema Vela CSV outputs","topologies":tops},indent=2)+"\n",encoding="utf-8")
def main():
 p=argparse.ArgumentParser();p.add_argument("--fixture",type=Path);p.add_argument("--out-dir",type=Path);p.add_argument("--make-synthetic-fixture",type=Path);a=p.parse_args()
 if a.make_synthetic_fixture:make_synthetic_fixture(a.make_synthetic_fixture);return
 if not a.fixture or not a.out_dir:p.error("--fixture and --out-dir are required")
 r=build_report(a.fixture);write_report(r,a.out_dir);print(f"PASS {SCHEMA}: node={len(r.node_rows)} edge={len(r.edge_rows)} triangle={len(r.triangle_rows)} figures=14")
if __name__=="__main__":main()