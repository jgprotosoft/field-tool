// Solveur éléments finis 2D plan (triangles linéaires)
//  - magnétostatique : ∇·(ν∇A) = −J − rot(ν·Br), Newton-Raphson pour les matériaux saturables
//  - électrostatique : ∇·(ε∇V) = −ρ
(function (G) {
  "use strict";
  const FT = (G.FT = G.FT || {});
  const MU0 = FT.MU0, EPS0 = FT.EPS0;

  // Géométrie élémentaire en mètres : gradients des fonctions de forme (b,c)/(2Δ)
  function elemGeom(mesh) {
    const nt = mesh.tri.length / 3, { x, y, tri } = mesh;
    const gb = new Float64Array(3 * nt), gc = new Float64Array(3 * nt), area = new Float64Array(nt);
    for (let e = 0; e < nt; e++) {
      const i = tri[3 * e], j = tri[3 * e + 1], k = tri[3 * e + 2];
      const xi = x[i] * 1e-3, yi = y[i] * 1e-3, xj = x[j] * 1e-3, yj = y[j] * 1e-3, xk = x[k] * 1e-3, yk = y[k] * 1e-3;
      gb[3 * e] = yj - yk; gb[3 * e + 1] = yk - yi; gb[3 * e + 2] = yi - yj;
      gc[3 * e] = xk - xj; gc[3 * e + 1] = xi - xk; gc[3 * e + 2] = xj - xi;
      area[e] = 0.5 * ((xj - xi) * (yk - yi) - (xk - xi) * (yj - yi));
    }
    return { gb, gc, area };
  }
  FT.elemGeom = elemGeom;

  // Structure creuse CSR pour les degrés de liberté libres
  function buildStructure(tri, dof, nd) {
    const nt = tri.length / 3;
    const rows = Array.from({ length: nd }, () => []);
    for (let e = 0; e < nt; e++)
      for (let a = 0; a < 3; a++) {
        const i = dof[tri[3 * e + a]];
        if (i < 0) continue;
        const r = rows[i];
        for (let b = 0; b < 3; b++) {
          const j = dof[tri[3 * e + b]];
          if (j >= 0 && r.indexOf(j) < 0) r.push(j);
        }
      }
    const rp = new Int32Array(nd + 1);
    for (let i = 0; i < nd; i++) rp[i + 1] = rp[i] + rows[i].length;
    const ci = new Int32Array(rp[nd]), dg = new Int32Array(nd);
    for (let i = 0; i < nd; i++) {
      const r = rows[i].sort((p, q) => p - q);
      for (let k = 0; k < r.length; k++) { ci[rp[i] + k] = r[k]; if (r[k] === i) dg[i] = rp[i] + k; }
    }
    const pos = new Int32Array(9 * nt).fill(-1);
    for (let e = 0; e < nt; e++)
      for (let a = 0; a < 3; a++) {
        const i = dof[tri[3 * e + a]];
        if (i < 0) continue;
        for (let b = 0; b < 3; b++) {
          const j = dof[tri[3 * e + b]];
          if (j < 0) continue;
          let lo = rp[i], hi = rp[i + 1] - 1;
          while (lo < hi) { const m = (lo + hi) >> 1; if (ci[m] < j) lo = m + 1; else hi = m; }
          pos[9 * e + 3 * a + b] = lo;
        }
      }
    return { rp, ci, dg, pos, va: new Float64Array(rp[nd]), nd };
  }

  // Gradient conjugué préconditionné Gauss-Seidel symétrique
  function pcg(M, b, x, tol, maxit) {
    const { rp, ci, va, dg, nd: n } = M;
    const r = new Float64Array(n), z = new Float64Array(n), p = new Float64Array(n), q = new Float64Array(n);
    const mul = (v, out) => {
      for (let i = 0; i < n; i++) { let s = 0; for (let k = rp[i]; k < rp[i + 1]; k++) s += va[k] * v[ci[k]]; out[i] = s; }
    };
    const prec = () => {
      for (let i = 0; i < n; i++) {
        let s = r[i];
        for (let k = rp[i]; k < dg[i]; k++) s -= va[k] * z[ci[k]];
        z[i] = s / va[dg[i]];
      }
      for (let i = n - 1; i >= 0; i--) {
        let s = 0;
        for (let k = dg[i] + 1; k < rp[i + 1]; k++) s += va[k] * z[ci[k]];
        z[i] -= s / va[dg[i]];
      }
    };
    mul(x, q);
    let bn = 0, rn = 0;
    for (let i = 0; i < n; i++) { r[i] = b[i] - q[i]; bn += b[i] * b[i]; rn += r[i] * r[i]; }
    if (bn === 0) { x.fill(0); return { it: 0, res: 0 }; }
    const target = tol * tol * bn;
    if (rn <= target) return { it: 0, res: Math.sqrt(rn / bn) };
    prec();
    p.set(z);
    let rz = 0;
    for (let i = 0; i < n; i++) rz += r[i] * z[i];
    let it = 0;
    for (; it < maxit; it++) {
      mul(p, q);
      let pq = 0;
      for (let i = 0; i < n; i++) pq += p[i] * q[i];
      const al = rz / pq;
      rn = 0;
      for (let i = 0; i < n; i++) { x[i] += al * p[i]; r[i] -= al * q[i]; rn += r[i] * r[i]; }
      if (rn <= target) { it++; break; }
      prec();
      let rz2 = 0;
      for (let i = 0; i < n; i++) rz2 += r[i] * z[i];
      const be = rz2 / rz;
      rz = rz2;
      for (let i = 0; i < n; i++) p[i] = z[i] + be * p[i];
    }
    return { it, res: Math.sqrt(rn / bn) };
  }

  // props[r] = {mur, Js, J (A/m²), brx, bry (T)} ; frontière extérieure A = 0
  FT.solveMagnetic = function (mesh, props) {
    const T0 = Date.now();
    const n = mesh.x.length, tri = mesh.tri, nt = tri.length / 3, reg = mesh.reg;
    const { gb, gc, area } = elemGeom(mesh);
    const dof = new Int32Array(n);
    let nd = 0;
    for (let i = 0; i < n; i++) dof[i] = mesh.outer[i] ? -1 : nd++;
    const M = buildStructure(tri, dof, nd);
    const nuE = new Float64Array(nt), dnuE = new Float64Array(nt);
    const f = new Float64Array(nd);
    let nonlinear = false;
    for (let e = 0; e < nt; e++) {
      const p = props[reg[e]];
      nuE[e] = 1 / (MU0 * p.mur);
      if (p.Js > 0) nonlinear = true;
      for (let a = 0; a < 3; a++) {
        const i = dof[tri[3 * e + a]];
        if (i < 0) continue;
        f[i] += (p.J * area[e]) / 3 + (nuE[e] * (p.brx * gc[3 * e + a] - p.bry * gb[3 * e + a])) / 2;
      }
    }
    const A = new Float64Array(n);
    const gxE = new Float64Array(nt), gyE = new Float64Array(nt);
    const update = (Av) => {
      for (let e = 0; e < nt; e++) {
        const a0 = Av[tri[3 * e]], a1 = Av[tri[3 * e + 1]], a2 = Av[tri[3 * e + 2]], d = 2 * area[e];
        const gx = (a0 * gb[3 * e] + a1 * gb[3 * e + 1] + a2 * gb[3 * e + 2]) / d;
        const gy = (a0 * gc[3 * e] + a1 * gc[3 * e + 1] + a2 * gc[3 * e + 2]) / d;
        gxE[e] = gx; gyE[e] = gy;
        const p = props[reg[e]];
        if (p.Js > 0) { const r = FT.nuOfB(p.mur, p.Js, Math.hypot(gx, gy)); nuE[e] = r[0]; dnuE[e] = r[1]; }
      }
    };
    const residual = (out) => {
      out.set(f);
      for (let e = 0; e < nt; e++) {
        const nu = nuE[e], gx = gxE[e], gy = gyE[e];
        for (let a = 0; a < 3; a++) {
          const i = dof[tri[3 * e + a]];
          if (i >= 0) out[i] -= (nu * (gx * gb[3 * e + a] + gy * gc[3 * e + a])) / 2;
        }
      }
      let s = 0;
      for (let i = 0; i < nd; i++) s += out[i] * out[i];
      return Math.sqrt(s);
    };
    const assemble = () => {
      const va = M.va, pos = M.pos;
      va.fill(0);
      for (let e = 0; e < nt; e++) {
        const nu = nuE[e], dn = dnuE[e], ar = area[e], gx = gxE[e], gy = gyE[e];
        for (let a = 0; a < 3; a++) {
          const qa = (gx * gb[3 * e + a] + gy * gc[3 * e + a]) / 2;
          for (let b = 0; b < 3; b++) {
            const k = pos[9 * e + 3 * a + b];
            if (k < 0) continue;
            let v = (nu * (gb[3 * e + a] * gb[3 * e + b] + gc[3 * e + a] * gc[3 * e + b])) / (4 * ar);
            if (dn > 0) v += ((2 * dn) / ar) * qa * ((gx * gb[3 * e + b] + gy * gc[3 * e + b]) / 2);
            va[k] += v;
          }
        }
      }
    };
    let fn = 0;
    for (let i = 0; i < nd; i++) fn += f[i] * f[i];
    fn = Math.sqrt(fn);
    const r = new Float64Array(nd), dx = new Float64Array(nd), At = new Float64Array(n);
    let iters = 0, cgIt = 0, converged = true, relRes = 0;
    if (fn > 0) {
      update(A);
      let rn = residual(r);
      const maxNewton = nonlinear ? 40 : 1;
      converged = !nonlinear;
      for (let it = 0; it < maxNewton; it++) {
        assemble();
        dx.fill(0);
        cgIt += pcg(M, r, dx, nonlinear ? 1e-6 : 1e-10, 8000).it;
        iters++;
        let al = 1, rn2 = 0;
        for (let ls = 0; ls < 8; ls++) {
          for (let i = 0; i < n; i++) At[i] = dof[i] >= 0 ? A[i] + al * dx[dof[i]] : 0;
          update(At);
          rn2 = residual(r);
          if (!nonlinear || rn2 < rn * (1 - 1e-4 * al)) break;
          al *= 0.5;
        }
        A.set(At);
        rn = rn2;
        relRes = rn / fn;
        if (nonlinear && relRes < 1e-6) { converged = true; break; }
      }
    }
    // post-traitement
    const Bx = new Float64Array(nt), By = new Float64Array(nt), Bm = new Float64Array(nt), Hm = new Float64Array(nt), wE = new Float64Array(nt);
    let W = 0, lam = 0, absI = 0;
    for (let e = 0; e < nt; e++) {
      const p = props[reg[e]];
      Bx[e] = gyE[e]; By[e] = -gxE[e];
      Bm[e] = Math.hypot(Bx[e], By[e]);
      const hx = nuE[e] * (Bx[e] - p.brx), hy = nuE[e] * (By[e] - p.bry);
      Hm[e] = Math.hypot(hx, hy);
      wE[e] = p.brx || p.bry ? 0.5 * (hx * (Bx[e] - p.brx) + hy * (By[e] - p.bry)) : FT.energyDensity(p.mur, p.Js, Bm[e], nuE[e]);
      W += wE[e] * area[e];
      if (p.J) {
        const am = (A[tri[3 * e]] + A[tri[3 * e + 1]] + A[tri[3 * e + 2]]) / 3;
        lam += p.J * am * area[e];
        absI += Math.abs(p.J) * area[e];
      }
    }
    return {
      type: "mag", A, Bx, By, Bm, Hm, nu: nuE, w: wE, energy: W, linkage: lam, absI, area,
      iters, cgIt, converged, relRes, nonlinear, ms: Date.now() - T0,
    };
  };

  // props[r] = {eps (F/m), rho (C/m³), V (nombre ou null)} ; bc : 'V0' ou 'none' sur la frontière extérieure
  FT.solveElectric = function (mesh, props, bc) {
    const T0 = Date.now();
    const n = mesh.x.length, tri = mesh.tri, nt = tri.length / 3, reg = mesh.reg;
    const { gb, gc, area } = elemGeom(mesh);
    const fixed = new Uint8Array(n), V = new Float64Array(n);
    let nfix = 0;
    if (bc === "V0") for (let i = 0; i < n; i++) if (mesh.outer[i]) { fixed[i] = 1; nfix++; }
    for (let e = 0; e < nt; e++) {
      const p = props[reg[e]];
      if (p.V != null) for (let a = 0; a < 3; a++) { const i = tri[3 * e + a]; fixed[i] = 1; V[i] = p.V; nfix++; }
    }
    let pinned = false;
    if (!nfix) { fixed[0] = 1; pinned = true; }
    const dof = new Int32Array(n);
    let nd = 0;
    for (let i = 0; i < n; i++) dof[i] = fixed[i] ? -1 : nd++;
    const M = buildStructure(tri, dof, nd);
    const f = new Float64Array(nd), va = M.va, pos = M.pos;
    for (let e = 0; e < nt; e++) {
      const p = props[reg[e]], ar = area[e];
      for (let a = 0; a < 3; a++) {
        const i = dof[tri[3 * e + a]];
        if (i < 0) continue;
        f[i] += (p.rho * ar) / 3;
        for (let b = 0; b < 3; b++) {
          const v = (p.eps * (gb[3 * e + a] * gb[3 * e + b] + gc[3 * e + a] * gc[3 * e + b])) / (4 * ar);
          const k = pos[9 * e + 3 * a + b];
          if (k >= 0) va[k] += v;
          else f[i] -= v * V[tri[3 * e + b]];
        }
      }
    }
    const x = new Float64Array(nd);
    const cg = pcg(M, f, x, 1e-10, 10000);
    for (let i = 0; i < n; i++) if (dof[i] >= 0) V[i] = x[dof[i]];
    const Ex = new Float64Array(nt), Ey = new Float64Array(nt), Em = new Float64Array(nt), Dm = new Float64Array(nt), wE = new Float64Array(nt);
    let W = 0;
    for (let e = 0; e < nt; e++) {
      const a0 = V[tri[3 * e]], a1 = V[tri[3 * e + 1]], a2 = V[tri[3 * e + 2]], d = 2 * area[e];
      Ex[e] = -(a0 * gb[3 * e] + a1 * gb[3 * e + 1] + a2 * gb[3 * e + 2]) / d;
      Ey[e] = -(a0 * gc[3 * e] + a1 * gc[3 * e + 1] + a2 * gc[3 * e + 2]) / d;
      Em[e] = Math.hypot(Ex[e], Ey[e]);
      const eps = props[reg[e]].eps;
      Dm[e] = eps * Em[e];
      wE[e] = 0.5 * eps * Em[e] * Em[e];
      W += wE[e] * area[e];
    }
    return { type: "elec", V, Ex, Ey, Em, Dm, w: wE, energy: W, area, cgIt: cg.it, converged: cg.res < 1e-6, pinned, ms: Date.now() - T0 };
  };

  // Valeur lissée aux sommets, par région (pas de lissage à travers les interfaces). Retourne 3 valeurs par triangle.
  FT.smoothPerRegion = function (mesh, val, area) {
    const nt = mesh.tri.length / 3, tri = mesh.tri, reg = mesh.reg, n = mesh.x.length;
    const out = new Float32Array(3 * nt);
    const byReg = Array.from({ length: mesh.nReg }, () => []);
    for (let e = 0; e < nt; e++) byReg[reg[e]].push(e);
    const sum = new Float64Array(n), wt = new Float64Array(n);
    for (const list of byReg) {
      for (const e of list) for (let a = 0; a < 3; a++) { const i = tri[3 * e + a]; sum[i] += val[e] * area[e]; wt[i] += area[e]; }
      for (const e of list) for (let a = 0; a < 3; a++) { const i = tri[3 * e + a]; out[3 * e + a] = sum[i] / wt[i]; }
      for (const e of list) for (let a = 0; a < 3; a++) { const i = tri[3 * e + a]; sum[i] = 0; wt[i] = 0; }
    }
    return out;
  };

  // Lignes de niveau d'un champ nodal : retourne Float32Array [x1,y1,x2,y2,...] (mm)
  FT.contours = function (mesh, u, nLev) {
    let mn = Infinity, mxv = -Infinity;
    for (let i = 0; i < u.length; i++) { if (u[i] < mn) mn = u[i]; if (u[i] > mxv) mxv = u[i]; }
    const out = [];
    if (!(mxv > mn)) return new Float32Array(0);
    const { x, y, tri } = mesh, nt = tri.length / 3, st = (mxv - mn) / nLev;
    for (let e = 0; e < nt; e++) {
      const i = tri[3 * e], j = tri[3 * e + 1], k = tri[3 * e + 2];
      const lo = Math.min(u[i], u[j], u[k]), hi = Math.max(u[i], u[j], u[k]);
      for (let l = Math.ceil((lo - mn) / st - 0.5); l <= Math.floor((hi - mn) / st - 0.5); l++) {
        const c = mn + (l + 0.5) * st;
        const pts = [];
        const cut = (a, b) => {
          if ((u[a] - c) * (u[b] - c) < 0) { const t = (c - u[a]) / (u[b] - u[a]); pts.push(x[a] + t * (x[b] - x[a]), y[a] + t * (y[b] - y[a])); }
        };
        cut(i, j); cut(j, k); cut(k, i);
        if (pts.length === 4) out.push(pts[0], pts[1], pts[2], pts[3]);
      }
    }
    return Float32Array.from(out);
  };

  // Statistiques d'une région : aire, moyenne/max du champ, énergie, courant, force (tenseur de Maxwell, méthode "coquille")
  FT.regionStats = function (mesh, res, r, depth) {
    const tri = mesh.tri, nt = tri.length / 3, reg = mesh.reg, area = res.area;
    const fld = res.type === "mag" ? res.Bm : res.Em;
    const inBody = new Uint8Array(mesh.x.length);
    let S = 0, sum = 0, mx = 0, W = 0;
    for (let e = 0; e < nt; e++) {
      if (reg[e] !== r) continue;
      S += area[e]; sum += fld[e] * area[e]; W += res.w[e] * area[e];
      if (fld[e] > mx) mx = fld[e];
      inBody[tri[3 * e]] = inBody[tri[3 * e + 1]] = inBody[tri[3 * e + 2]] = 1;
    }
    const { gb, gc } = elemGeom(mesh);
    let fx = 0, fy = 0, contact = false;
    for (let e = 0; e < nt; e++) {
      if (reg[e] === r) continue;
      let ggx = 0, ggy = 0, any = false;
      for (let a = 0; a < 3; a++) if (inBody[tri[3 * e + a]]) { ggx += gb[3 * e + a]; ggy += gc[3 * e + a]; any = true; }
      if (!any) continue;
      let k, vx, vy;
      if (res.type === "mag") {
        k = res.nu[e]; vx = res.Bx[e]; vy = res.By[e];
        if (Math.abs(k * FT.MU0 - 1) > 0.02) contact = true;
      } else {
        k = res.Dm[e] > 0 ? res.Dm[e] / res.Em[e] : FT.EPS0; vx = res.Ex[e]; vy = res.Ey[e];
        if (Math.abs(k / FT.EPS0 - 1) > 0.02) contact = true;
      }
      const h2 = (vx * vx + vy * vy) / 2;
      // F = −∫ T·∇g dΩ, ∇g = (ggx, ggy)/(2Δ)
      fx -= (k * ((vx * vx - h2) * ggx + vx * vy * ggy)) / 2;
      fy -= (k * (vx * vy * ggx + (vy * vy - h2) * ggy)) / 2;
    }
    return { area: S * 1e6, mean: S ? sum / S : 0, max: mx, energy: W * depth, fx: fx * depth, fy: fy * depth, contact };
  };

  // ---------- Liaison modèle -> problème ----------
  // S : {edges, labels, customMats, airMargin, meshN, ...}
  FT.buildMesh = function (S) {
    if (!S.edges.length) return null;
    const bb = FT.modelBBox(S.edges);
    const Dm = Math.max(bb.x1 - bb.x0, bb.y1 - bb.y0);
    if (!(Dm > 0)) return null;
    const m = (S.airMargin || 1) * Dm;
    const box = { x0: bb.x0 - m, y0: bb.y0 - m, x1: bb.x1 + m, y1: bb.y1 + m };
    const pslg = FT.buildPSLG(S.edges, box);
    const mesh = FT.generateMesh(pslg, { hmax: Dm / (S.meshN || 40), refine: 32, grad: 0.25 });
    mesh.locator = FT.buildLocator(mesh);
    const g = elemGeom(mesh);
    mesh.regArea = new Float64Array(mesh.nReg); // m²
    for (let e = 0; e < mesh.reg.length; e++) mesh.regArea[mesh.reg[e]] += g.area[e];
    FT.assignLabels(S, mesh);
    return mesh;
  };

  // Associe chaque étiquette de matériau à la région qui la contient
  FT.assignLabels = function (S, mesh) {
    mesh.regLabel = new Array(mesh.nReg).fill(null);
    for (const l of S.labels) {
      const t = mesh.locator.find(l.x, l.y);
      if (t >= 0) mesh.regLabel[mesh.reg[t]] = l;
    }
  };

  FT.regionProps = function (S, mesh, physics) {
    const out = [];
    for (let r = 0; r < mesh.nReg; r++) {
      const l = mesh.regLabel[r];
      const m = FT.getMaterial(S, l ? l.mat : "air");
      if (physics === "mag") {
        let J = 0, brx = 0, bry = 0;
        if (l && l.srcType === "J") J = (+l.srcVal || 0) * 1e6;
        if (l && l.srcType === "NI") J = (+l.srcVal || 0) / mesh.regArea[r];
        if (m.Br > 0) { const a = (((l && +l.magAngle) || 0) * Math.PI) / 180; brx = m.Br * Math.cos(a); bry = m.Br * Math.sin(a); }
        out.push({ mur: m.mur, Js: m.Js, J, brx, bry });
      } else {
        const fixedV = l && l.elecType === "V";
        out.push({
          eps: EPS0 * (m.cond && !fixedV ? 1e5 : m.epsr),
          rho: l && l.elecType === "rho" ? +l.elecVal || 0 : 0,
          V: fixedV ? +l.elecVal || 0 : null,
        });
      }
    }
    return out;
  };

  FT.solve = function (S, mesh) {
    FT.assignLabels(S, mesh);
    const props = FT.regionProps(S, mesh, S.physics);
    return S.physics === "mag" ? FT.solveMagnetic(mesh, props) : FT.solveElectric(mesh, props, S.bcElec || "V0");
  };
})(typeof window !== "undefined" ? window : globalThis);
