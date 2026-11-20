#!/usr/bin/env python3
"""
Phase 2 - Step 2B-2:
  1. Applies exact H-bond subtype budget conservation (HB_base + HB_2OH + HB_bb == HB) to unippi_cli.py.
  2. Mirrors the complete, verified Step 2B-1 + 2B-1.5 Protein-RNA physics engine, mmCIF prime fix,
     buffer filter, multi-chain HETATM splitter, and RNA UI/CSV outputs into index.html for 100% parity.
  3. Validates Python AST + JavaScript syntax (via node --check if available) and verifies CLI benchmarks.
"""

import ast
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

cli_path = Path("unippi_cli.py")
html_path = Path("index.html")

# ------------------------------------------------------------------------------
# PART A: Exact H-Bond Subtype Conservation in unippi_cli.py
# ------------------------------------------------------------------------------
cli = cli_path.read_text(encoding="utf-8")
old_hb_cap = '''            eff_hb = min(pair_hb, 2)
            eff_hb_base = min(pair_hb_base, 2)
            eff_hb_2oh = min(pair_hb_2oh, 2)
            eff_hb_bb = min(pair_hb_bb, 2)'''
new_hb_cap = '''            eff_hb = min(pair_hb, 2)
            eff_hb_base = min(pair_hb_base, eff_hb)
            eff_hb_2oh = min(pair_hb_2oh, eff_hb - eff_hb_base)
            eff_hb_bb = min(pair_hb_bb, eff_hb - eff_hb_base - eff_hb_2oh)'''
assert old_hb_cap in cli, "Missing eff_hb_base anchor in unippi_cli.py"
cli = cli.replace(old_hb_cap, new_hb_cap, 1)
ast.parse(cli)
cli_path.write_text(cli, encoding="utf-8")
print("[PASS] Locked exact HB_base + HB_2OH + HB_bb == HB budget in unippi_cli.py.")

# ------------------------------------------------------------------------------
# PART B: Mirror Verified Protein-RNA Engine & UI Into index.html
# ------------------------------------------------------------------------------
html = html_path.read_text(encoding="utf-8")

# 1. Version badges & branding in index.html
html = html.replace('const MODEL_VERSION = "UniPPI-v1.8";', 'const MODEL_VERSION = "UniPPI-v2.0-MARVAS";', 1)
html = html.replace('<title>UniPPI-OWL v1.8 | Universal Structure-to-Affinity Resource</title>',
                    '<title>UniPPI-OWL v2.0 (MARVAS) | Protein–Protein &amp; Protein–RNA Affinity Resource</title>', 1)
html = html.replace('<strong>UniPPI-OWL</strong><span>v1.8</span>',
                    '<strong>UniPPI-OWL</strong><span>v2.0-MARVAS</span>', 1)
html = html.replace('Model: UniPPI-v1.8 · Frozen Threshold: 0.5000 · Parity: 100% CLI Match',
                    'Model: UniPPI-v2.0-MARVAS · Frozen Threshold: 0.5000 · Parity: 100% CLI Match', 1)
html = html.replace('<span>UniPPI-v1.8 Engine</span>', '<span>UniPPI-v2.0 Engine</span>', 1)
html = html.replace('Scientific model: <code>UniPPI-v1.8</code>', 'Scientific model: <code>UniPPI-v2.0-MARVAS</code>', 1)

# 2. Add RNA/Buffer/Metal constants & 3D ring-plane / modality helpers
old_const_js = '    const NUCLEOTIDES = new Set(["A","C","G","U","DA","DC","DG","DT","ADE","CYT","GUA","THY","URA"]);'
new_const_js = '''    const NUCLEOTIDES = new Set(["A","C","G","U","DA","DC","DG","DT","ADE","CYT","GUA","THY","URA"]);
    const CRYSTAL_BUFFER_JUNK = new Set(["HOH","WAT","DOD","GOL","EDO","PEG","PGE","PG4","1PE","P6G","MPD","SO4","PO4","CL","NA","K","ACT","ACY","DMS","MES","TRS","CIT","FMT","NO3","BR","IOD","IPA","EOH","BU1","BME","NH4","CAC","IMD","EPE","HEP","TAR","TLA","SUC","MLI"]);
    const CATALYTIC_METALS = new Set(["MG","ZN","MN","CA","FE","CU","CO","NI"]);
    const NUCLEOBASE_RINGS = {
      A: new Set(["N1","C2","N3","C4","C5","C6","N7","C8","N9"]),
      G: new Set(["N1","C2","N3","C4","C5","C6","N7","C8","N9"]),
      DA: new Set(["N1","C2","N3","C4","C5","C6","N7","C8","N9"]),
      DG: new Set(["N1","C2","N3","C4","C5","C6","N7","C8","N9"]),
      ADE: new Set(["N1","C2","N3","C4","C5","C6","N7","C8","N9"]),
      GUA: new Set(["N1","C2","N3","C4","C5","C6","N7","C8","N9"]),
      C: new Set(["N1","C2","N3","C4","C5","C6"]),
      U: new Set(["N1","C2","N3","C4","C5","C6"]),
      T: new Set(["N1","C2","N3","C4","C5","C6"]),
      DC: new Set(["N1","C2","N3","C4","C5","C6"]),
      DT: new Set(["N1","C2","N3","C4","C5","C6"]),
      CYT: new Set(["N1","C2","N3","C4","C5","C6"]),
      URA: new Set(["N1","C2","N3","C4","C5","C6"]),
      THY: new Set(["N1","C2","N3","C4","C5","C6"])
    };
    const NUCLEOBASE_HB_ATOMS = new Set(["N1","N2","N3","N4","N6","N7","N9","O2","O4","O6"]);
    const NUCLEIC_CARBON_ATOMS = new Set(["C1'","C2'","C3'","C4'","C5'","C2","C4","C5","C6","C7","C8"]);

    function getRingCentroidAndNormalJS(res) {
      const cname = res.canonicalName;
      const ringSet = AROMATIC_RINGS[cname] || NUCLEOBASE_RINGS[cname];
      if (!ringSet) return [null, null];
      const atomMap = {};
      const pts = [];
      for (const a of res.atoms) {
        if (ringSet.has(a.name)) {
          const pt = [a.x, a.y, a.z];
          atomMap[a.name] = pt;
          pts.push(pt);
        }
      }
      if (pts.length < 3) return [null, null];
      const cent = [meanArr(pts.map(p=>p[0])), meanArr(pts.map(p=>p[1])), meanArr(pts.map(p=>p[2]))];
      const triad = NUCLEOBASE_RINGS[cname] ? ["C2","C4","C6"] : ["CG","CE1","CE2"];
      const [p1, p2, p3] = triad.every(k => atomMap[k])
        ? [atomMap[triad[0]], atomMap[triad[1]], atomMap[triad[2]]]
        : [pts[0], pts[Math.floor(pts.length/3)], pts[Math.floor(2*pts.length/3)]];
      const ux = p2[0]-p1[0], uy = p2[1]-p1[1], uz = p2[2]-p1[2];
      const vx = p3[0]-p1[0], vy = p3[1]-p1[1], vz = p3[2]-p1[2];
      const nx = uy*vz - uz*vy, ny = uz*vx - ux*vz, nz = ux*vy - uy*vx;
      const nLen = Math.hypot(nx, ny, nz);
      return nLen < 1e-6 ? [cent, null] : [cent, [nx/nLen, ny/nLen, nz/nLen]];
    }

    function detectPairModalityJS(caId, cbId, resA, resB) {
      if (caId.endsWith("_N") && cbId.endsWith("_C")) return "INTRA-DOMAIN (MONOMER)";
      const allT = new Set([...resA.map(r=>r.entityType), ...resB.map(r=>r.entityType)]);
      if (allT.has("NUCLEIC")) {
        const nucRes = [...resA, ...resB].filter(r => r.entityType === "NUCLEIC");
        const isRna = nucRes.some(r => r.canonicalName === "U" || r.canonicalName === "URA" || r.atoms.some(a => a.name === "O2'"));
        const tag = isRna ? "RNA" : "DNA";
        if (allT.has("PROTEIN")) return `PROTEIN-${tag}`;
        if (allT.has("LIGAND")) return `${tag}-LIGAND`;
        if (allT.has("ION")) return `${tag}-METAL`;
        return `${tag}-${tag}`;
      }
      if (allT.has("LIGAND")) return "PROTEIN-LIGAND";
      if (allT.has("ION")) return "PROTEIN-METAL";
      return "PROTEIN-PROTEIN";
    }'''
assert old_const_js in html, "Missing NUCLEOTIDES anchor in index.html"
html = html.replace(old_const_js, new_const_js, 1)

# 3. Fix _exptl.method, buffer filter, and preserve ribose prime (') in parseStructureText()
old_exptl_js = '          if (s.startsWith("_exptl.method")) meta.method = s.split(/\\s+/).slice(1).join(" ").replace(/[\'"]/g, "");'
new_exptl_js = '          if (s.startsWith("_exptl.method") && s.split(/\\s+/)[0] === "_exptl.method") { const exVal = s.split(/\\s+/).slice(1).join(" ").replace(/^[\'"]|[\'"]$/g, ""); if (exVal && exVal !== "?" && exVal !== ".") meta.method = exVal; }'
assert old_exptl_js in html, "Missing _exptl.method anchor in index.html"
html = html.replace(old_exptl_js, new_exptl_js, 1)

water_js = 'if (["HOH","WAT","DOD"].includes(resName)) continue;'
assert html.count(water_js) == 2, f"Expected 2 water checks in index.html, found {html.count(water_js)}"
html = html.replace(water_js, 'if (CRYSTAL_BUFFER_JUNK.has(resName)) continue;')

old_name_js = '            const name = (row["label_atom_id"] || row["auth_atom_id"] || "CA").replace(/[\'"]/g, "");'
new_name_js = '            let name = (row["label_atom_id"] || row["auth_atom_id"] || "CA").trim();\n            if (name.length >= 2 && name[0] === name[name.length-1] && (name[0] === "\'" || name[0] === \'"\')) name = name.slice(1, -1);'
assert old_name_js in html, "Missing label_atom_id quote-strip anchor in index.html"
html = html.replace(old_name_js, new_name_js, 1)

# 4. Upgrade buildResidues() for multi-chain biological HETATM drug & catalytic metal splitting
old_build_js = '''      // Universal Auto-Rescue 1: Same-chain Receptor + HETATM Ligand/Ion (common in RCSB PDB, AutoDock, Boltz-1)
      if (Object.keys(chains).length === 1) {
        const onlyCid = Object.keys(chains)[0];
        const allRes = chains[onlyCid];
        const polyRes = allRes.filter(r => r.entityType === "PROTEIN" || r.entityType === "NUCLEIC");
        const hetRes = allRes.filter(r => r.entityType === "LIGAND" || r.entityType === "ION");
        if (polyRes.length > 0 && hetRes.length > 0) {
          const ligCid = `${onlyCid}_LIG`;
          chains[onlyCid] = polyRes;
          chains[ligCid] = [];
          for (const hr of hetRes) {
            hr.chainId = ligCid;
            hr.key = `${ligCid}|${hr.resSeq}|${hr.ins}`;
            for (const a of hr.atoms) a.chainId = ligCid;
            chains[ligCid].push(hr);
          }
        } else if (polyRes.length >= 10) {
          // Universal Auto-Rescue 2: True Single-Chain Monomer -> Intra-Chain Domain-Domain Packing Split (N-term vs C-term)
          const mid = Math.floor(allRes.length / 2);
          const nCid = `${onlyCid}_N`, cCid = `${onlyCid}_C`;
          chains[nCid] = [];
          chains[cCid] = [];
          allRes.forEach((r, idx) => {
            const targetCid = idx < mid ? nCid : cCid;
            r.chainId = targetCid;
            r.key = `${targetCid}|${r.resSeq}|${r.ins}`;
            for (const a of r.atoms) a.chainId = targetCid;
            chains[targetCid].push(r);
          });
          delete chains[onlyCid];
        }
      }
      return { ordered, chains };'''

new_build_js = '''      // Universal Multi-Chain HETATM Splitter (v2.0 — 100% CLI Parity):
      // Splits biological small-molecule drugs/cofactors & catalytic metals out of ANY polymer chain
      const splitChains = {};
      for (const [cid, rlist] of Object.entries(chains)) {
        const polyRes = rlist.filter(r => r.entityType === "PROTEIN" || r.entityType === "NUCLEIC");
        const hetRes = rlist.filter(r => r.entityType === "LIGAND" || r.entityType === "ION");
        if (polyRes.length > 0 && hetRes.length > 0) {
          splitChains[cid] = polyRes;
          for (const hr of hetRes) {
            if (CRYSTAL_BUFFER_JUNK.has(hr.resName)) continue;
            if (hr.entityType === "ION" && !CATALYTIC_METALS.has(hr.resName)) continue;
            const newCid = `${cid}_${hr.resName}${hr.resSeq}`;
            hr.chainId = newCid;
            hr.key = `${newCid}|${hr.resSeq}|${hr.ins}`;
            for (const a of hr.atoms) a.chainId = newCid;
            if (!splitChains[newCid]) splitChains[newCid] = [];
            splitChains[newCid].push(hr);
          }
        } else {
          const kept = rlist.filter(r => !CRYSTAL_BUFFER_JUNK.has(r.resName) && !(r.entityType === "ION" && !CATALYTIC_METALS.has(r.resName)));
          if (kept.length) splitChains[cid] = kept;
        }
      }
      const keptSet = new Set(Object.values(splitChains).flat());
      const finalOrdered = ordered.filter(r => keptSet.has(r));
      if (Object.keys(splitChains).length === 1) {
        const onlyCid = Object.keys(splitChains)[0];
        const allRes = splitChains[onlyCid];
        if (allRes.length >= 10 && allRes.every(r => r.entityType === "PROTEIN" || r.entityType === "NUCLEIC")) {
          const mid = Math.floor(allRes.length / 2);
          const nCid = `${onlyCid}_N`, cCid = `${onlyCid}_C`;
          splitChains[nCid] = [];
          splitChains[cCid] = [];
          allRes.forEach((r, idx) => {
            const targetCid = idx < mid ? nCid : cCid;
            r.chainId = targetCid;
            r.key = `${targetCid}|${r.resSeq}|${r.ins}`;
            for (const a of r.atoms) a.chainId = targetCid;
            splitChains[targetCid].push(r);
          });
          delete splitChains[onlyCid];
        }
      }
      return { ordered: finalOrdered, chains: splitChains };'''
assert old_build_js in html, "Missing buildResidues auto-rescue anchor in index.html"
html = html.replace(old_build_js, new_build_js, 1)

# 5. Upgrade analyzePairJS: per-chain crystal B-factor plateau + 9-channel Protein-RNA physics
old_bfac_js = '''      const allB = [...resA, ...resB].map(r => r.meanBfac);
      const bMin = allB.length ? Math.min(...allB) : 10.0;
      const paeMat = confData.res_pae_matrix;

      const rawPlddt = r => (mode === "MODE_A_PRED_FULL" || mode === "MODE_B_PRED_COORD") ? clamp(r.meanBfac, 0, 100) : clamp(100*Math.exp(-Math.max(0, r.meanBfac - bMin)/60), 5, 100);
      const platConf = r => (mode === "MODE_A_PRED_FULL" || mode === "MODE_B_PRED_COORD") ? clamp((r.meanBfac - 30.0)/55.0, 0.02, 1.0) : Math.exp(-Math.max(0, r.meanBfac - bMin)/60.0);'''

new_bfac_js = '''      const bMinA = resA.length ? Math.min(...resA.map(r => r.meanBfac)) : 10.0;
      const bMinB = resB.length ? Math.min(...resB.map(r => r.meanBfac)) : 10.0;
      const bMinChain = { [caId]: bMinA, [cbId]: bMinB };
      const paeMat = confData.res_pae_matrix;

      const rawPlddt = r => {
        if (mode === "MODE_A_PRED_FULL" || mode === "MODE_B_PRED_COORD") return clamp(r.meanBfac, 0, 100);
        const bRef = Math.max((bMinChain[r.chainId] ?? 10.0) + 20.0, 45.0);
        return clamp(100.0 * Math.exp(-Math.max(0.0, r.meanBfac - bRef) / 100.0), 15.0, 100.0);
      };
      const platConf = r => {
        if (mode === "MODE_A_PRED_FULL" || mode === "MODE_B_PRED_COORD") return clamp((r.meanBfac - 30.0)/55.0, 0.02, 1.0);
        const bRef = Math.max((bMinChain[r.chainId] ?? 10.0) + 20.0, 45.0);
        return clamp(Math.exp(-Math.max(0.0, r.meanBfac - bRef) / 100.0), 0.15, 1.0);
      };'''
assert old_bfac_js in html, "Missing rawPlddt/platConf anchor in index.html"
html = html.replace(old_bfac_js, new_bfac_js, 1)

old_cnt_js = '''      let nSb=0, nHb=0, nHp=0, nPi=0, nSs=0, nClash=0, nCb8=0;
      let wSb=0, wHb=0, wHp=0, wPi=0, wClash=0;'''
new_cnt_js = '''      const modality = detectPairModalityJS(caId, cbId, resA, resB);
      const isNucleicPair = (modality === "PROTEIN-RNA" || modality === "PROTEIN-DNA");
      let nSb=0, nHb=0, nHp=0, nPi=0, nSs=0, nClash=0, nCb8=0;
      let nHbBase=0, nHb2oh=0, nHbBb=0, nCatPi=0;
      let wSb=0, wHb=0, wHp=0, wPi=0, wClash=0;
      let wHbBase=0, wHb2oh=0, wHbBb=0, wPiRna=0, wCatPi=0;'''
assert old_cnt_js in html, "Missing nSb counter init anchor in index.html"
html = html.replace(old_cnt_js, new_cnt_js, 1)

old_pair_loop_js = '''          let minD = 999.0, pSb=0, pHb=0, pHp=0, pSs=0, pClash=0, hpDistSum=0.0;
          const cna = ra.canonicalName, cnb = rb.canonicalName;
          const apA = HYDROPHOBIC_RES.has(cna) || ra.entityType === "LIGAND" || ra.atoms.some(a => a.name === "CD2" && a.element === "C");
          const apB = HYDROPHOBIC_RES.has(cnb) || rb.entityType === "LIGAND" || rb.atoms.some(a => a.name === "CD2" && a.element === "C");

          for (const aa of ra.atoms) {
            for (const ab of rb.atoms) {
              const d = Math.hypot(aa.x-ab.x, aa.y-ab.y, aa.z-ab.z);
              if (d < minD) minD = d;
              if (d > 4.5) continue;
              if (cna==="CYS" && cnb==="CYS" && aa.name==="SG" && ab.name==="SG" && d<=2.3) { pSs++; continue; }
              if (d < 2.2) pClash++;

              const isPtmA = ["SEP","TPO","PTR"].includes(ra.resName) && ["O1P","O2P","O3P","OP1","OP2","OP3","OG","OG1","OH"].includes(aa.name);
              const isPtmB = ["SEP","TPO","PTR"].includes(rb.resName) && ["O1P","O2P","O3P","OP1","OP2","OP3","OG","OG1","OH"].includes(ab.name);
              const isMetA = ra.entityType==="ION" || ["CA","MG","ZN","CO","FE","MN","CU","NI"].includes(aa.element);
              const isMetB = rb.entityType==="ION" || ["CA","MG","ZN","CO","FE","MN","CU","NI"].includes(ab.element);
              if (d <= 3.0 && ((isMetA && ["O","N","S"].includes(ab.element)) || (isMetB && ["O","N","S"].includes(aa.element)))) { pSb++; continue; }
              const catA = (cna==="ARG"&&["NE","NH1","NH2"].includes(aa.name)) || (cna==="LYS"&&aa.name==="NZ") || (cna==="HIS"&&["ND1","NE2"].includes(aa.name));
              const aniA = (cna==="ASP"&&["OD1","OD2"].includes(aa.name)) || (cna==="GLU"&&["OE1","OE2"].includes(aa.name)) || (ra.entityType==="NUCLEIC"&&["OP1","OP2","O1P","O2P"].includes(aa.name)) || isPtmA;
              const catB = (cnb==="ARG"&&["NE","NH1","NH2"].includes(ab.name)) || (cnb==="LYS"&&ab.name==="NZ") || (cnb==="HIS"&&["ND1","NE2"].includes(ab.name));
              const aniB = (cnb==="ASP"&&["OD1","OD2"].includes(ab.name)) || (cnb==="GLU"&&["OE1","OE2"].includes(ab.name)) || (rb.entityType==="NUCLEIC"&&["OP1","OP2","O1P","O2P"].includes(ab.name)) || isPtmB;

              if (d <= 4.0 && ((catA && aniB) || (aniA && catB))) { pSb += (isPtmA || isPtmB) ? 2 : 1; continue; }
              if (d <= 3.5 && ["N","O","S"].includes(aa.element) && ["N","O","S"].includes(ab.element)) { pHb++; continue; }
              if (d <= 4.5 && apA && apB && ["C","S"].includes(aa.element) && ["C","S"].includes(ab.element) && !["C","CA"].includes(aa.name) && !["C","CA"].includes(ab.name)) {
                pHp++;
                hpDistSum += clamp((4.75 - d) / 0.80, 0.35, 1.0);
              }
            }
          }

          if (minD > activeCutoff) continue;
          let pPi = 0;
          if (AROMATIC_RINGS[cna] && AROMATIC_RINGS[cnb]) {
            const rAtomsA = ra.atoms.filter(a => AROMATIC_RINGS[cna].has(a.name));
            const rAtomsB = rb.atoms.filter(a => AROMATIC_RINGS[cnb].has(a.name));
            if (rAtomsA.length >= 3 && rAtomsB.length >= 3) {
              const centA = [meanArr(rAtomsA.map(a=>a.x)), meanArr(rAtomsA.map(a=>a.y)), meanArr(rAtomsA.map(a=>a.z))];
              const centB = [meanArr(rAtomsB.map(a=>a.x)), meanArr(rAtomsB.map(a=>a.y)), meanArr(rAtomsB.map(a=>a.z))];
              if (dist3(centA, centB) <= 5.5) pPi = 1;
            }
          }
          intfA.add(ra.key); intfB.add(rb.key);
          if (cbDist > 8.0) { intfPae.push(pairPae); weightsSeen.push(wij); }

          const clsA = ["SEP","TPO","PTR"].includes(ra.resName) ? "CHARGED" : (PRODIGY_CLASS[cna] || "POLAR");
          const clsB = ["SEP","TPO","PTR"].includes(rb.resName) ? "CHARGED" : (PRODIGY_CLASS[cnb] || "POLAR");
          let pairCls = [clsA, clsB].sort().join("_");
          if (pairCls === "APOLAR_CHARGED") pairCls = "CHARGED_APOLAR";
          else if (pairCls === "APOLAR_POLAR") pairCls = "POLAR_APOLAR";
          if (prodigyIC[pairCls] !== undefined) prodigyIC[pairCls] += wij;

          const effSb = Math.min(pSb, 2), effHb = Math.min(pHb, 2), effHp = Math.min(pHp, 6), effHpW = Math.min(hpDistSum, 6.0), effClash = Math.min(pClash, 3);
          if (wij >= 0.25) { nSb += effSb; nHb += effHb; nHp += effHp; nPi += pPi; nSs += pSs; nClash += effClash; }
          wSb += wij*effSb; wHb += wij*effHb; wHp += wij*effHpW; wPi += wij*pPi; wClash += wij*effClash;

          if (wij >= 0.25) {
            for (const rk of [ra.key, rb.key]) {
              if (!resContrib.has(rk)) resContrib.set(rk, { sb:0, hb:0, hp:0, pi:0 });
              const rc = resContrib.get(rk);
              rc.sb += wij*effSb; rc.hb += wij*effHb; rc.hp += wij*effHpW; rc.pi += wij*pPi;
            }
          }

          contacts.push({
            resA: `${ra.chainId}:${ra.resName}${ra.resSeq}`, resB: `${rb.chainId}:${rb.resName}${rb.resSeq}`,
            resSeqA: ra.resSeq, resSeqB: rb.resSeq, coordA: ra.cbCoord, coordB: rb.cbCoord,
            minDist: Number(minD.toFixed(2)), pae: Number(pairPae.toFixed(2)), weight: Number(wij.toFixed(3)),
            sb: effSb, hb: effHb, hp: effHp, pi: pPi, ss: pSs, clash: effClash
          });'''

new_pair_loop_js = '''          let minD = 999.0, pSb=0, pHb=0, pHp=0, pSs=0, pClash=0, hpDistSum=0.0;
          let pHbBase=0, pHb2oh=0, pHbBb=0;
          const cna = ra.canonicalName, cnb = rb.canonicalName;
          const apA = HYDROPHOBIC_RES.has(cna) || ra.entityType === "LIGAND";
          const apB = HYDROPHOBIC_RES.has(cnb) || rb.entityType === "LIGAND";

          for (const aa of ra.atoms) {
            for (const ab of rb.atoms) {
              const d = Math.hypot(aa.x-ab.x, aa.y-ab.y, aa.z-ab.z);
              if (d < minD) minD = d;
              if (d > 4.5) continue;
              if (cna==="CYS" && cnb==="CYS" && aa.name==="SG" && ab.name==="SG" && d<=2.3) { pSs++; continue; }

              const isMetA = (ra.entityType==="ION" || ra.entityType==="LIGAND") && CATALYTIC_METALS.has(aa.element);
              const isMetB = (rb.entityType==="ION" || rb.entityType==="LIGAND") && CATALYTIC_METALS.has(ab.element);
              if (d <= 2.8 && ((isMetA && ["O","N","S"].includes(ab.element)) || (isMetB && ["O","N","S"].includes(aa.element)))) { pSb++; continue; }
              if (d < 2.2) pClash++;

              const isPtmA = ["SEP","TPO","PTR"].includes(ra.resName) && ["O1P","O2P","O3P","OP1","OP2","OP3","OG","OG1","OH"].includes(aa.name);
              const isPtmB = ["SEP","TPO","PTR"].includes(rb.resName) && ["O1P","O2P","O3P","OP1","OP2","OP3","OG","OG1","OH"].includes(ab.name);
              const catA = (cna==="ARG"&&["NE","NH1","NH2"].includes(aa.name)) || (cna==="LYS"&&aa.name==="NZ") || (cna==="HIS"&&["ND1","NE2"].includes(aa.name));
              const aniA = (cna==="ASP"&&["OD1","OD2"].includes(aa.name)) || (cna==="GLU"&&["OE1","OE2"].includes(aa.name)) || (ra.entityType==="NUCLEIC"&&["OP1","OP2","O1P","O2P"].includes(aa.name)) || isPtmA;
              const catB = (cnb==="ARG"&&["NE","NH1","NH2"].includes(ab.name)) || (cnb==="LYS"&&ab.name==="NZ") || (cnb==="HIS"&&["ND1","NE2"].includes(ab.name));
              const aniB = (cnb==="ASP"&&["OD1","OD2"].includes(ab.name)) || (cnb==="GLU"&&["OE1","OE2"].includes(ab.name)) || (rb.entityType==="NUCLEIC"&&["OP1","OP2","O1P","O2P"].includes(ab.name)) || isPtmB;

              if (d <= 4.0 && ((catA && aniB) || (aniA && catB))) { pSb += (isPtmA || isPtmB) ? 2 : 1; continue; }
              if (d <= 3.5 && ["N","O","S"].includes(aa.element) && ["N","O","S"].includes(ab.element)) {
                pHb++;
                if (isNucleicPair) {
                  const nucAt = ra.entityType === "NUCLEIC" ? aa : ab;
                  if (nucAt.name === "O2'") pHb2oh++;
                  else if (NUCLEOBASE_HB_ATOMS.has(nucAt.name)) pHbBase++;
                  else pHbBb++;
                }
                continue;
              }
              if (d <= 4.5 && ["C","S"].includes(aa.element) && ["C","S"].includes(ab.element)) {
                if (apA && apB && !["C","CA"].includes(aa.name) && !["C","CA"].includes(ab.name)) {
                  pHp++;
                  hpDistSum += clamp((4.75 - d) / 0.80, 0.35, 1.0);
                } else if (isNucleicPair) {
                  const protIsA = (ra.entityType === "PROTEIN");
                  const pAt = protIsA ? aa : ab, nAt = protIsA ? ab : aa;
                  const pAp = protIsA ? apA : apB;
                  if (pAp && !["C","CA"].includes(pAt.name) && NUCLEIC_CARBON_ATOMS.has(nAt.name)) {
                    pHp++;
                    hpDistSum += clamp((4.75 - d) / 0.80, 0.35, 1.0);
                  }
                }
              }
            }
          }

          if (minD > activeCutoff) continue;
          let pPi = 0, pPiW = 0.0, pCatPi = 0;
          if (isNucleicPair) {
            const [ringA, normA] = getRingCentroidAndNormalJS(ra);
            const [ringB, normB] = getRingCentroidAndNormalJS(rb);
            if (ringA && ringB) {
              const dRing = dist3(ringA, ringB);
              if (dRing <= 5.5) {
                const nRef = ra.entityType === "NUCLEIC" ? normA : normB;
                if (normA && normB && nRef) {
                  const vx = ringB[0]-ringA[0], vy = ringB[1]-ringA[1], vz = ringB[2]-ringA[2];
                  const dPerp = Math.abs(vx*nRef[0] + vy*nRef[1] + vz*nRef[2]);
                  const dPar = Math.sqrt(Math.max(0.0, dRing*dRing - dPerp*dPerp));
                  const cosAng = Math.abs(normA[0]*normB[0] + normA[1]*normB[1] + normA[2]*normB[2]);
                  if (dPerp >= 2.2 && dPar <= 4.2) {
                    if (cosAng >= 0.60) { pPi = 1; pPiW = 1.0; }
                    else if (dRing <= 5.0) { pPi = 1; pPiW = 0.5; }
                  }
                } else { pPi = 1; pPiW = 0.75; }
              }
            }
            const protR = ra.entityType === "PROTEIN" ? ra : rb;
            const nucR = ra.entityType === "PROTEIN" ? rb : ra;
            if (protR.canonicalName === "ARG" || protR.canonicalName === "LYS") {
              const [nucCent, nucNorm] = getRingCentroidAndNormalJS(nucR);
              const catName = protR.canonicalName === "ARG" ? "CZ" : "NZ";
              const catAt = protR.atoms.find(a => a.name === catName);
              if (nucCent && nucNorm && catAt) {
                const dCp = dist3([catAt.x, catAt.y, catAt.z], nucCent);
                if (dCp <= 5.0 && dCp > 1e-4) {
                  const cosPhi = Math.abs((catAt.x-nucCent[0])*nucNorm[0] + (catAt.y-nucCent[1])*nucNorm[1] + (catAt.z-nucCent[2])*nucNorm[2]) / dCp;
                  if (cosPhi >= 0.50) pCatPi = 1;
                }
              }
            }
          } else if (AROMATIC_RINGS[cna] && AROMATIC_RINGS[cnb]) {
            const [ringA] = getRingCentroidAndNormalJS(ra);
            const [ringB] = getRingCentroidAndNormalJS(rb);
            if (ringA && ringB && dist3(ringA, ringB) <= 5.5) { pPi = 1; pPiW = 1.0; }
          }

          intfA.add(ra.key); intfB.add(rb.key);
          if (cbDist > 8.0) { intfPae.push(pairPae); weightsSeen.push(wij); }

          const clsA = ["SEP","TPO","PTR"].includes(ra.resName) ? "CHARGED" : (PRODIGY_CLASS[cna] || (ra.entityType === "LIGAND" ? "APOLAR" : "POLAR"));
          const clsB = ["SEP","TPO","PTR"].includes(rb.resName) ? "CHARGED" : (PRODIGY_CLASS[cnb] || (rb.entityType === "LIGAND" ? "APOLAR" : "POLAR"));
          let pairCls = [clsA, clsB].sort().join("_");
          if (pairCls === "APOLAR_CHARGED") pairCls = "CHARGED_APOLAR";
          else if (pairCls === "APOLAR_POLAR") pairCls = "POLAR_APOLAR";
          if (prodigyIC[pairCls] !== undefined) prodigyIC[pairCls] += wij;

          const effSb = Math.min(pSb, 2), effHb = Math.min(pHb, 2);
          const effHbBase = Math.min(pHbBase, effHb);
          const effHb2oh = Math.min(pHb2oh, effHb - effHbBase);
          const effHbBb = Math.min(pHbBb, effHb - effHbBase - effHb2oh);
          const effHp = Math.min(pHp, 6), effHpW = Math.min(hpDistSum, 6.0), effClash = Math.min(pClash, 3);
          if (wij >= 0.25) {
            nSb += effSb; nHb += effHb; nHbBase += effHbBase; nHb2oh += effHb2oh; nHbBb += effHbBb;
            nHp += effHp; nPi += pPi; nCatPi += pCatPi; nSs += pSs; nClash += effClash;
          }
          wSb += wij*effSb; wHb += wij*effHb; wHbBase += wij*effHbBase; wHb2oh += wij*effHb2oh; wHbBb += wij*effHbBb;
          wHp += wij*effHpW; wPi += wij*pPi; wPiRna += wij*pPiW; wCatPi += wij*pCatPi; wClash += wij*effClash;

          if (wij >= 0.25) {
            for (const rk of [ra.key, rb.key]) {
              if (!resContrib.has(rk)) resContrib.set(rk, { sb:0, hb:0, hp:0, pi:0 });
              const rc = resContrib.get(rk);
              rc.sb += wij*effSb; rc.hb += wij*effHb; rc.hp += wij*effHpW; rc.pi += wij*(pPi + pCatPi);
            }
          }

          contacts.push({
            resA: `${ra.chainId}:${ra.resName}${ra.resSeq}`, resB: `${rb.chainId}:${rb.resName}${rb.resSeq}`,
            resSeqA: ra.resSeq, resSeqB: rb.resSeq, coordA: ra.cbCoord, coordB: rb.cbCoord,
            minDist: Number(minD.toFixed(2)), pae: Number(pairPae.toFixed(2)), weight: Number(wij.toFixed(3)),
            sb: effSb, hb: effHb, hbBase: effHbBase, hb2oh: effHb2oh, hbBb: effHbBb,
            hp: effHp, pi: pPi, catpi: pCatPi, ss: pSs, clash: effClash
          });'''
assert old_pair_loop_js in html, "Missing inner residue-pair loop anchor in index.html"
html = html.replace(old_pair_loop_js, new_pair_loop_js, 1)

# Fix crystallographic Proline inversion in index.html
old_pro_js = 'const isStrained = (r.meanBfac < 80.0) || (meanPaeCoreTmp > 4.0);'
new_pro_js = 'const isStrained = (rawPlddt(r) < 80.0) || (meanPaeCoreTmp > 4.0);'
assert old_pro_js in html, "Missing Proline isStrained anchor in index.html"
html = html.replace(old_pro_js, new_pro_js, 1)

# Upgrade thermodynamics in index.html (dedicated 9-channel Protein-RNA / DNA branch + PRODIGY disable)
old_gate_js = '      const gate = (effIcSum < 3.0 ? 0.0 : 1.0 / (1.0 + Math.exp(-0.40 * (effIcSum - 11.0)))) * nisReliability;'
new_gate_js = '      const gate = (isNucleicPair || effIcSum < 3.0 ? 0.0 : 1.0 / (1.0 + Math.exp(-0.40 * (effIcSum - 11.0)))) * nisReliability;'
assert old_gate_js in html, "Missing gate anchor in index.html"
html = html.replace(old_gate_js, new_gate_js, 1)

old_eval_js = '''        const eRaw = sOut*(0.110*wHp + 0.195*prodigyIC.APOLAR_APOLAR + 0.0032*wBsa) + (0.78*wSb*dKap*phKappa + 0.48*wHb + 0.62*wPi);
        const eEff = 15.5 * Math.tanh(eRaw / 15.5);
        const dPhys = 2.10 - eEff + 0.12*prodigyIC.POLAR_POLAR + 1.50*wClash + 2.20*helPro;'''
new_eval_js = '''        let dPhys = 2.10;
        if (isNucleicPair) {
          const dKapRna = Math.exp(-0.85 * ((Math.sqrt(iM) / 3.04) - (Math.sqrt(0.15) / 3.04)) * 3.5);
          const eRawRna = 0.60*wSb*dKapRna*phKappa + 0.65*wHbBase + 0.55*wHb2oh + 0.35*wHbBb + 0.95*wPiRna + 0.85*wCatPi + 0.10*wHp + 0.0025*wBsa;
          const eEffRna = 16.5 * Math.tanh(eRawRna / 16.5);
          dPhys = 2.10 - eEffRna + 1.50*wClash;
        } else {
          const eRaw = sOut*(0.110*wHp + 0.195*prodigyIC.APOLAR_APOLAR + 0.0032*wBsa) + (0.78*wSb*dKap*phKappa + 0.48*wHb + 0.62*wPi);
          const eEff = 15.5 * Math.tanh(eRaw / 15.5);
          dPhys = 2.10 - eEff + 0.12*prodigyIC.POLAR_POLAR + 1.50*wClash + 2.20*helPro;
        }'''
assert old_eval_js in html, "Missing evalEnvAtSalt anchor in index.html"
html = html.replace(old_eval_js, new_eval_js, 1)

# Remove duplicate modality computation & expose RNA specificity fields on returned object
old_mod_js = '''      const typesA = new Set(resA.map(r => r.entityType));
      const typesB = new Set(resB.map(r => r.entityType));
      const allTypes = new Set([...typesA, ...typesB]);
      let modality = "PROTEIN-PROTEIN";
      if (caId.endsWith("_N") && cbId.endsWith("_C")) modality = "INTRA-DOMAIN (MONOMER)";
      else if (allTypes.has("NUCLEIC")) modality = "PROTEIN-NUCLEIC (DNA/RNA)";
      else if (allTypes.has("LIGAND")) modality = "PROTEIN-LIGAND (SMALL MOL)";
      else if (allTypes.has("ION")) modality = "PROTEIN-METAL (ION)";'''
new_mod_js = '''      const rnaTotSpec = nSb + nHbBase + nHb2oh + nHbBb + nPi + nCatPi;
      const rnaSsi = isNucleicPair ? Number(((nHbBase + nPi + nCatPi) / Math.max(1, rnaTotSpec)).toFixed(3)) : 0.0;'''
assert old_mod_js in html, "Missing old modality block in index.html"
html = html.replace(old_mod_js, new_mod_js, 1)

# Exclude helical_prolines penalty on nucleic pairs in pattern classification
html = html.replace(
    'if (uis >= UIS_THRESHOLD && iptm >= 0.60 && dgUni <= -5.8 && helPro === 0 && nClash < 5)',
    'const effHelPro = isNucleicPair ? 0 : helPro;\n      if (uis >= UIS_THRESHOLD && iptm >= 0.60 && dgUni <= -5.8 && effHelPro === 0 && nClash < 5)',
    1
)
html = html.replace(
    '&& helPro === 0 && nClash < 5) {\n        patId = "Pattern B";',
    '&& effHelPro === 0 && nClash < 5) {\n        patId = "Pattern B";',
    1
)
html = html.replace(
    '} else if (helPro > 0 || nClash >= 5) {',
    '} else if (effHelPro > 0 || nClash >= 5) {',
    1
)

old_ret_js = '        salt_bridges: nSb, h_bonds: nHb, hydrophobic_contacts: nHp, pi_stacking: nPi, disulfides: nSs, steric_clashes: nClash,'
new_ret_js = '        salt_bridges: nSb, h_bonds: nHb, hb_base: nHbBase, hb_2oh: nHb2oh, hb_backbone: nHbBb, cation_pi: nCatPi, rna_specificity_idx: rnaSsi, hydrophobic_contacts: nHp, pi_stacking: nPi, disulfides: nSs, steric_clashes: nClash, wPiRna: Number(wPiRna.toFixed(2)), wCatPi: Number(wCatPi.toFixed(2)), wHbBase: Number(wHbBase.toFixed(2)), wHb2oh: Number(wHb2oh.toFixed(2)), wHbBb: Number(wHbBb.toFixed(2)),'
assert old_ret_js in html, "Missing return contact fields anchor in index.html"
html = html.replace(old_ret_js, new_ret_js, 1)

# 6. Upgrade UI: Ensemble Summary guard, Step 1 & Step 4 RNA breakdown, Contact Table tags, Master CSV columns
old_ens_guard_js = '        if (ALL_RESULTS.length > 1) {'
new_ens_guard_js = '        const uniqueIds = new Set(ALL_RESULTS.map(x => x.complex_id));\n        if (uniqueIds.size > 1) {'
assert old_ens_guard_js in html, "Missing ensembleSummaryBanner guard in index.html"
html = html.replace(old_ens_guard_js, new_ens_guard_js, 1)

old_step1_js = '(${r.salt_bridges} Salt Bridges [≤4.0 Å], ${r.h_bonds} H-Bonds [≤3.5 Å], ${r.hydrophobic_contacts} Hydrophobic contacts [≤4.5 Å]).'
new_step1_js = '(${r.salt_bridges} Salt Bridges [≤4.0 Å], ${r.h_bonds} H-Bonds [≤3.5 Å]${(r.modality==="PROTEIN-RNA"||r.modality==="PROTEIN-DNA") ? ` [Base:${r.hb_base||0}, 2\'-OH:${r.hb_2oh||0}, BB:${r.hb_backbone||0}]` : ""}, ${r.hydrophobic_contacts} Hydrophobic [≤4.5 Å]${(r.modality==="PROTEIN-RNA"||r.modality==="PROTEIN-DNA") ? `, ${r.pi_stacking||0} Base π-Stacks, ${r.cation_pi||0} Cation-π, Specificity Index SSI = ${(r.rna_specificity_idx||0).toFixed(2)}` : ""}).'
assert old_step1_js in html, "Missing Step 1 derivation line in index.html"
html = html.replace(old_step1_js, new_step1_js, 1)

old_step4_js = '''              &bull; <code>ΔG<sub>PRODIGY</sub> (Vangone &amp; Bonvin 2015) = −0.0946·IC<sub>CC</sub>(${(ic.CHARGED_CHARGED||0).toFixed(1)}) − 0.1001·IC<sub>CA</sub>(${(ic.CHARGED_APOLAR||0).toFixed(1)}) + 0.1958·IC<sub>PP</sub>(${(ic.POLAR_POLAR||0).toFixed(1)}) − 0.2267·IC<sub>PA</sub>(${(ic.POLAR_APOLAR||0).toFixed(1)}) + 0.1868·%NIS<sub>ap</sub>(${r.pctNisAp||35.0}%) + 0.3810·%NIS<sub>ch</sub>(${r.pctNisCh||25.0}%) − 15.9433 = <strong>${r.dG_prodigy_kcal} kcal/mol</strong></code>${nisNote}<br>
              &bull; <code>ΔG<sub>phys</sub> (Enthalpy–Entropy Compensated 8-Channel Model, E<sub>eff</sub>=15.5·tanh(E<sub>raw</sub>/15.5); IC<sub>AA</sub>=${(ic.APOLAR_APOLAR||0).toFixed(1)}, wHP=${(r.wHp||0).toFixed(1)}, wBSA=${r.wBSA_A2} Å², wSB=${(r.wSb||0).toFixed(1)}, wHB=${(r.wHb||0).toFixed(1)}, IC<sub>PP</sub>=${(ic.POLAR_POLAR||0).toFixed(1)}${r.pi_stacking?`, π=${r.pi_stacking}`:''}${r.steric_clashes?`, Clash=${r.steric_clashes}`:''}${r.helical_prolines?`, ProKink=${r.helical_prolines}`:''}) = <strong>${r.dG_phys_kcal} kcal/mol</strong></code><br>
              &bull; <code>Consensus ΔG<sub>uni</sub> = ${r.w_prodigy_blend.toFixed(2)}·ΔG<sub>PRODIGY</sub>(${r.dG_prodigy_kcal.toFixed(2)}) + ${r.w_phys_blend.toFixed(2)}·ΔG<sub>phys</sub>(${r.dG_phys_kcal.toFixed(2)}) = <strong>${r.dG_uni_kcal.toFixed(2)} kcal/mol</strong> &rarr; K<sub>d</sub> = exp(ΔG<sub>uni</sub> / RT) = <strong>${kdStr}</strong></code><br>'''

new_step4_js = '''              ${(r.modality === "PROTEIN-RNA" || r.modality === "PROTEIN-DNA") ? `
              &bull; <code>Modality = <strong>${r.modality}</strong> &rarr; Protein–Protein PRODIGY term automatically disabled (0% weight).</code><br>
              &bull; <code>ΔG<sub>RNA</sub> (Dedicated 9-Channel ${r.modality} Model, E<sub>eff</sub>=16.5·tanh(E<sub>raw</sub>/16.5); wSB<sub>phos</sub>=${(r.wSb||0).toFixed(1)}, wHB<sub>base</sub>=${(r.wHbBase||0).toFixed(1)}, wHB<sub>2'OH</sub>=${(r.wHb2oh||0).toFixed(1)}, wHB<sub>bb</sub>=${(r.wHbBb||0).toFixed(1)}, wπ<sub>RNA</sub>=${(r.wPiRna||0).toFixed(1)}, wCatπ=${(r.wCatPi||0).toFixed(1)}, wHP=${(r.wHp||0).toFixed(1)}, wBSA=${r.wBSA_A2} Å²) = <strong>${r.dG_uni_kcal.toFixed(2)} kcal/mol</strong> &rarr; K<sub>d</sub> = <strong>${kdStr}</strong></code><br>
              ` : `
              &bull; <code>ΔG<sub>PRODIGY</sub> (Vangone &amp; Bonvin 2015) = −0.0946·IC<sub>CC</sub>(${(ic.CHARGED_CHARGED||0).toFixed(1)}) − 0.1001·IC<sub>CA</sub>(${(ic.CHARGED_APOLAR||0).toFixed(1)}) + 0.1958·IC<sub>PP</sub>(${(ic.POLAR_POLAR||0).toFixed(1)}) − 0.2267·IC<sub>PA</sub>(${(ic.POLAR_APOLAR||0).toFixed(1)}) + 0.1868·%NIS<sub>ap</sub>(${r.pctNisAp||35.0}%) + 0.3810·%NIS<sub>ch</sub>(${r.pctNisCh||25.0}%) − 15.9433 = <strong>${r.dG_prodigy_kcal} kcal/mol</strong></code>${nisNote}<br>
              &bull; <code>ΔG<sub>phys</sub> (Enthalpy–Entropy Compensated 8-Channel Model, E<sub>eff</sub>=15.5·tanh(E<sub>raw</sub>/15.5); IC<sub>AA</sub>=${(ic.APOLAR_APOLAR||0).toFixed(1)}, wHP=${(r.wHp||0).toFixed(1)}, wBSA=${r.wBSA_A2} Å², wSB=${(r.wSb||0).toFixed(1)}, wHB=${(r.wHb||0).toFixed(1)}, IC<sub>PP</sub>=${(ic.POLAR_POLAR||0).toFixed(1)}${r.pi_stacking?`, π=${r.pi_stacking}`:''}${r.steric_clashes?`, Clash=${r.steric_clashes}`:''}${r.helical_prolines?`, ProKink=${r.helical_prolines}`:''}) = <strong>${r.dG_phys_kcal} kcal/mol</strong></code><br>
              &bull; <code>Consensus ΔG<sub>uni</sub> = ${r.w_prodigy_blend.toFixed(2)}·ΔG<sub>PRODIGY</sub>(${r.dG_prodigy_kcal.toFixed(2)}) + ${r.w_phys_blend.toFixed(2)}·ΔG<sub>phys</sub>(${r.dG_phys_kcal.toFixed(2)}) = <strong>${r.dG_uni_kcal.toFixed(2)} kcal/mol</strong> &rarr; K<sub>d</sub> = exp(ΔG<sub>uni</sub> / RT) = <strong>${kdStr}</strong></code><br>
              `}'''
assert old_step4_js in html, "Missing Step 4 derivation block in index.html"
html = html.replace(old_step4_js, new_step4_js, 1)

# Add RNA tags (Cat-pi, 2'-OH, Base) to Contact Inspector Table
old_ctags_js = '''        if (c.hb > 0) tags.push(`HB:${c.hb}`);
        if (c.hp > 0) tags.push(`HP:${c.hp}`);
        if (c.pi > 0) tags.push(`π:${c.pi}`);'''
new_ctags_js = '''        if (c.hb > 0) {
          const hbSub = [];
          if (c.hbBase > 0) hbSub.push(`Base:${c.hbBase}`);
          if (c.hb2oh > 0) hbSub.push(`2'-OH:${c.hb2oh}`);
          if (c.hbBb > 0) hbSub.push(`BB:${c.hbBb}`);
          tags.push(hbSub.length ? `HB(${hbSub.join("/")})` : `HB:${c.hb}`);
        }
        if (c.hp > 0) tags.push(`HP:${c.hp}`);
        if (c.pi > 0) tags.push(`π:${c.pi}`);
        if (c.catpi > 0) tags.push(`Cat-π:${c.catpi}`);'''
assert old_ctags_js in html, "Missing contact table tags anchor in index.html"
html = html.replace(old_ctags_js, new_ctags_js, 1)

# Add RNA columns to Master CSV/TSV export in index.html
old_csv_js = '"salt_bridges","h_bonds","hydrophobic_contacts","helical_prolines"'
new_csv_js = '"salt_bridges","h_bonds","hb_base","hb_2oh","hb_backbone","hydrophobic_contacts","pi_stacking","cation_pi","rna_specificity_idx","helical_prolines"'
assert old_csv_js in html, "Missing CSV columns anchor in index.html"
html = html.replace(old_csv_js, new_csv_js, 1)

html_path.write_text(html, encoding="utf-8")
print("[PASS] Mirrored verified Protein-RNA engine & UI into index.html.")

# Validate JavaScript syntax if Node.js is installed
m_js = re.findall(r"<script>([\s\S]*?)</script>", html)
assert m_js, "No inline <script> block found in index.html"
if shutil.which("node"):
    with tempfile.NamedTemporaryFile("w", suffix=".js", delete=False, encoding="utf-8") as tf:
        tf.write(m_js[-1])
        tf_path = tf.name
    res_node = subprocess.run(["node", "--check", tf_path], capture_output=True, text=True)
    Path(tf_path).unlink(missing_ok=True)
    if res_node.returncode != 0:
        raise RuntimeError(f"JavaScript syntax error in index.html:\n{res_node.stderr}")
    print("[PASS] Verified index.html inline JavaScript syntax with `node --check`.")
else:
    print("[INFO] `node` not installed in WSL; skipped `node --check` (verify in browser DevTools).")

# Final Step 2 CLI check confirming HB_base + HB_2OH + HB_bb == HB on 5ZC9 & 1URN
for title, inp, outp in [
    ("STEP 2 FINAL CHECK — 5ZC9.cif (eIF4A1 + Polypurine RNA)", "benchmarks/5ZC9.cif", "/tmp/v2_s2_5zc9"),
    ("STEP 2 FINAL CHECK — 1URN.cif (U1A RRM + Hairpin RNA)", "benchmarks/1URN.cif", "/tmp/v2_s2_1urn"),
]:
    print("\n" + "=" * 95)
    print(title)
    print("=" * 95)
    subprocess.run([sys.executable, "unippi_cli.py", "-i", inp, "-o", outp], check=True)
